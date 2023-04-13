# =========================================================
# DistributeTraffic.py
# @author Nina Thomsen
# @date 02.02.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Distribute traffic between TAZ and create OD Matrices;
#        Assign inner cell traffic to network links
# =========================================================

import itertools

import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd


def pt_shares_int(a, b):
    """
    Determine shares of international transport for personal transport for a country based on area and share of land border

    :param a: area of country
    :param b: share of land border of country
    :return: dict with international transport shares
    :type a: float
    :type b: float
    :rtype: dict
    """
    # transform area
    a = (np.log((a-(-15000))/6.59)/0.437)-22.69

    x1 = (0.3 / (1 + ((0.3 / 0.095 - 1) * np.exp(a * 0.3 * -(-1.05))))) * b
    x2 = max(0.02, (0.35 + -0.04 * a) * b**0.5)
    x3 = min(0.95, max(0.1, (0.3 + 0.05 * a)))

    return {'pt_share_foreign': x1, 'pt_share_transit': x2, 'pt_ratio_int': x3}


def get_trip_matrix(mx_grav, taz, goal_col, taz_id='id'):
    """
    Get trip matrix based on total trips per taz and gravity model matrix

    :param mx_grav: matrix from gravity model
    :param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ; needs id column to identify coordinate in mx_grav
    :param goal_col: name of column in taz with total trips per taz
    :return: matrix with trips between taz
    :type mx_grav: np.array
    :type taz: gpd.GeoDataFrame or pd.DataFrame
    :type goal_col: str
    :rtype: np.array
    """
    # choice probabilities and trips
    mx_trips = np.zeros(mx_grav.shape)
    for t in taz[taz_id]:
        goal = float(taz.loc[taz[taz_id] == t, goal_col].iloc[0])
        # get probabilities for trips starting at t
        prob_a = mx_grav[t, :] / mx_grav[t, :].sum()
        prob_b = mx_grav[:, t] / mx_grav[:, t].sum()
        mx_trips[t, :] += prob_a * (goal / 2)
        mx_trips[:, t] += prob_b * (goal / 2)
    # fill NA
    mx_trips[np.isnan(mx_trips)] = 0
    return mx_trips


def assignment_single(net_g, net, o, d, weight, trips, from_='from_node', to_='to_node'):
    """
    Performs a shortest path network assignment between two nodes on a MultiDiGraph based on a specified weight and returns a GeoDataFrame with the resulting network load
    :param net_g: MultiDiGraph with network
    :param net: gpd.GeoDataFrame with network
    :param o: ID of origin node
    :param d: ID of destination node
    :param weight: name of weight attribute; is an attribute of links in net and net_g, e.g. travel time or length
    :param trips: number of vehicle trips to distribute between o and d
    :param from_: name of from node attribute in net, defaults to 'from_node'
    :param to_: name of to node attribute in net, defaults to 'to_node'
    :return: GeoDataFrame of net with trips on the links, on the path between o and d

    :type net_g: MultiDiGraph
    :type net: gpd.GeoDataFrame
    :type o: int
    :type d: int
    :type weight: str
    :type trips: float
    :type from_: str
    :type to_: str
    :rtype: gpd.GeoDataFrame
    """
    path_nodes = ox.distance.shortest_path(net_g, o, d, weight=weight)  # tt or length
    nodes_tuples = [(path_nodes[x], path_nodes[x + 1]) for x in range(len(path_nodes) - 1)]
    nodes_tuples = pd.DataFrame(nodes_tuples, columns=[from_, to_])
    nodes_tuples['__tmp'] = trips
    # add trips to network links between nodes tuples
    net_tmp = net.merge(nodes_tuples, on=[from_, to_], how='left')
    net_tmp.loc[net_tmp['__tmp'].isna(), '__tmp'] = 0
    return net_tmp


def assignment_multiple(net_g, net, o, d, k, weight, trips, from_='from_node', to_='to_node', id_='link_id'):
    """
    Performs a shortest path assignment for the k shortest paths between two nodes on a MultiDiGraph based on a specified weight and returns a GeoDataFrame with the resulting network load
    Paths are weighted based on the total weight per path
    :param net_g: MultiDiGraph with network
    :param net: gpd.GeoDataFrame with network
    :param o: ID of origin node
    :param d: ID of destination node
    :param k: int, number of paths to find
    :param weight: name of weight attribute; is an attribute of links in net and net_g, e.g. travel time or length
    :param trips: number of vehicle trips to distribute between o and d
    :param from_: name of from node attribute in net, defaults to 'from_node'
    :param to_: name of to node attribute in net, defaults to 'to_node'
    :param id_: name of id attribute in net, defaults to 'link_id'
    :return: GeoDataFrame of net with trips on the links, on the path between o and d

    :type net_g: MultiDiGraph
    :type net: gpd.GeoDataFrame
    :type o: int
    :type d: int
    :type k: int
    :type weight: str
    :type trips: float
    :type from_: str
    :type to_: str
    :type id_: str
    :rtype: gpd.GeoDataFrame
    """
    # transform to digraph for nx.shortest_simple_path function
    dg = nx.DiGraph(net_g)
    k_paths = nx.shortest_simple_paths(dg, o, d, weight)
    net_tmp = net[[from_, to_, id_, weight]].copy()
    net_tmp['__weight_all'] = 0
    weights = []
    for path_nodes in itertools.islice(k_paths, 0, k):
        nodes_tuples = [(path_nodes[x], path_nodes[x + 1]) for x in range(len(path_nodes) - 1)]
        nodes_tuples = pd.DataFrame(nodes_tuples, columns=[from_, to_])
        # merge net_tmp, aggregate weight (tt or length) for current path
        nodes_tuples = nodes_tuples.merge(net[[from_, to_, weight]], on=[from_, to_], how='left')
        weight_ = 1 / nodes_tuples[weight].sum()
        nodes_tuples['__weight'] = weight_
        weights.append(weight_)
        # merge weight to net_tmp and add to weight_all
        net_tmp = net_tmp.merge(nodes_tuples[[from_, to_, '__weight']], on=[from_, to_], how='left')
        net_tmp.loc[net_tmp['__weight'].isna(), '__weight'] = 0
        net_tmp['__weight_all'] += net_tmp['__weight']
        net_tmp.drop(columns=['__weight'], inplace=True)
    k_paths.close()
    net_tmp['__weight_all'] /= sum(weights)
    net_tmp['__tmp'] = net_tmp['__weight_all'] * trips
    return net_tmp


def create_network_graph(net_, nodes_, node_id='node_id', node_geo='geometry', net_from='from_node', net_to='to_node', net_id='link_id'):
    """
    Transform a GeoDataFrame of a network and the corresponding nodes to a MultiDiGraph
    Needs columns for node id in nodes, from node, to node and link id in network
    :param net_: GeoDataFrame with network
    :param nodes_: GeoDataFrame with nodes
    :param node_id: name of id column in nodes_, defaults to 'node_id'
    :param node_geo: name of geometry column in nodes_, defaults to 'geometry'
    :param net_from: name of from node column in net_, defaults to 'from_node'
    :param net_to: name of to node column in net_, defaults to 'to_node'
    :param net_id: name of id column in net_, defaults to 'link_id'
    :return: MultiDiGraph of net_

    :type net_: gpd.GeoDataFrame
    :type nodes_: gpd.GeoDataFrame
    :type node_id: str
    :type node_geo: str
    :type net_from: str
    :type net_to: str
    :type net_id: str
    :rtype: MultiDiGraph
    """
    nodes_.set_index(node_id, inplace=True)
    nodes_['x'] = nodes_[node_geo].x
    nodes_['y'] = nodes_[node_geo].y
    net_.set_index([net_from, net_to, net_id], inplace=True)

    net_g = ox.utils_graph.graph_from_gdfs(nodes_, net_)
    net_g = ox.utils_graph.get_largest_component(net_g)

    return net_g


class TargetValues:
    """
    Determine shares of international and national transport for a country and calculate corresponding total trips and vehicle kilometers travelled (VKT)
    """

    def __init__(self, country_layer, cn_col='cntr_code'):
        """

        :param country_layer: countries with target values for raod transport volumes
        :param cn_col: name of column with country code in country layer, default 'cntr_code'
        :type country_layer: pd.DataFrame or gpd.GeoDataFrame
        :type cn_col: str
        """
        self.country_layer = country_layer
        self.cn_col = cn_col

    def targets_personal_transport(self, cn, target_col='car_pkm', unit=1.e+6):
        """
        Calculate VKT for the segments short and long distance (national) and international personal transport for a country

        :param cn: country code
        :param target_col: name of column in self.country_layer with the total transport volume in person kilometers travelled, default 'car_pkm'
        :param unit: scaling factor for total transport volume, default 1.e+6
        :return: vkt per segment

        :type cn: str
        :type target_col: str
        :type unit: float
        :rtype: dict
        """
        target = float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, target_col].iloc[0])
        shares_int = pt_shares_int(a=float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, 'area'].iloc[0]),
                                   b=float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, 'border_share'].iloc[0]))

        f_in_c = target * shares_int['pt_share_foreign']
        int1 = target * shares_int['pt_share_foreign'] * shares_int['pt_ratio_int']
        inner = target - f_in_c - int1

        short_ = 0.65 * (inner + int1)
        long_ = 0.35 * (inner + int1) - int1

        trans_ = shares_int['pt_share_transit'] * f_in_c
        int_ = (f_in_c - trans_) + int1

        return {"short": unit * short_, "long": unit * long_, "int": unit * int_, "trans": unit * trans_}

    def targets_freight_transport(self, cn, target_col='freight_tkm', segments=None, shares_tkm=None, loads_t=None, seg_split=None, unit=1.e+6):
        """
        Calculate VKT and trips for the segments short and long distance (national) and international freight transport for a country

        !! Long-distance and international freight transport only cover 1 segment (heavy weight transport)
        !! Short-distance freight transport can cover multiple segments, inlcuding the long-distance segment (heavy weight transport)
        !! The distribution of tkm among segments is given with the shares_tkm parameter, short-distance tkm are then further spread among these segments

        :param seg_split: freight transport segment where there is traffic in all categories (short, long, international)
        :param cn: country code
        :param target_col: name of column in self.country_layer with the total transport volume in tkm
        :param segments: names of freight transport segments; default None refers to ['lcv', 'mft', 'hft']
        :param shares_tkm: shares of different freight transport segments of total tkm; default None refers to shares {'lcv': 0.08, 'mft': 0.09, 'hft': 0.83}
        :param loads_t: average load in t per freight transport segment; default None refers to load {'lcv': 0.5, 'mft': 3, 'hft': 10}
        :param unit: scaling factor for total transport volume, default 1.e+6
        :return: vkt and trips per segment

        :type seg_split: str
        :type cn: str
        :type target_col: str
        :type segments: list
        :type shares_tkm: dict
        :type loads_t: dict
        :type unit: float
        :rtype: dict
        """

        if shares_tkm is None:
            shares_tkm = {'lcv': 0.08, 'mft': 0.09, 'hft': 0.83}
        if loads_t is None:
            loads_t = {'lcv': 0.5, 'mft': 3, 'hft': 10}
        if seg_split is None:
            seg_split = 'hft'
        if segments is None:
            segments = shares_tkm.keys()

        # check if segments, loads and shares match
        if (sorted(segments) == sorted(shares_tkm.keys())) & (sorted(segments) == sorted(loads_t.keys())) & (seg_split in segments):

            target = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, target_col].iloc[0])

            cat_dict = {0: {'int_imex': 0.03, 'int_transit': 0.00, 'nat': 0.97},
                        1: {'int_imex': 0.15, 'int_transit': 0.02, 'nat': 0.83},
                        2: {'int_imex': 0.15, 'int_transit': 0.05, 'nat': 0.80},
                        3: {'int_imex': 0.20, 'int_transit': 0.10, 'nat': 0.70},
                        4: {'int_imex': 0.00, 'int_transit': 0.00, 'nat': 1.00}}

            border_share = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'border_share'].iloc[0])
            neighbors = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'neighbors'].iloc[0])
            border_crossings = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'border_crossings_count'].iloc[0])

            # find category for country
            cat_country = 2
            if (border_share < 0.1) & (neighbors <= 1):
                cat_country = 0
            elif neighbors <= 5:
                if (border_crossings < 4) | (border_share < 0.5):
                    cat_country = 1
                else:
                    cat_country = 2
            elif neighbors >= 5:
                if border_share < 0.75:
                    cat_country = 2
                else:
                    cat_country = 3
            elif border_share == 0:
                cat_country = 4

            # distances
            area = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'area'].iloc[0])
            dis_long = 0.365 * area**(1/2) + 61.488
            dis_short = 1.183 * area**(1/4) + 22.361
            dis_transit = 0.928 * area**(1/2) + 38.998

            # transform input tkm to t:
            share_national = cat_dict[cat_country]['nat']
            share_international = cat_dict[cat_country]['int_imex']
            share_transit = cat_dict[cat_country]['int_transit']

            t_tot = target / ((dis_long * 0.2 + dis_short * 0.8) * share_national + dis_long * share_international + dis_transit * share_transit)
            t_nat = t_tot * share_national
            t_int = t_tot * share_international
            t_tra = t_tot * share_transit

            # tkm
            tkm_long = t_nat * 0.2 * dis_long
            tkm_short = t_nat * 0.8 * dis_short
            tkm_imex = t_int * dis_long
            tkm_transit = t_tra * dis_transit
            tkm = sum([tkm_imex, tkm_short, tkm_long, tkm_transit])

            tkm_dict = {'short': tkm_short, 'long': tkm_long}

            # transform tkm to vkm by using average tkm shares and average load
            # short distance transport: split into vehicle segments
            vkm_short_segments = {}
            for s in segments:
                # total
                tkm_s = shares_tkm[s] * tkm
                # short distance only segments
                if s != seg_split:
                    vkm_ts = tkm_s / loads_t[s]
                else:
                    # determine short distance tkm
                    tkm_s_short = tkm_s - (tkm - tkm_short)
                    vkm_ts = tkm_s_short / loads_t[s]
                vkm_short_segments.update({s: vkm_ts * unit})
            vkm_short = sum([vkm_short_segments[s] for s in segments])

            # vkm long-distance and international
            vkm_long = tkm_long / loads_t[seg_split] * unit
            vkm_imex = tkm_imex / loads_t[seg_split] * unit
            vkm_transit = tkm_transit / loads_t[seg_split] * unit

            # result dictionary
            result = {'short': vkm_short, 'trips_short': vkm_short / dis_short,
                      'long': vkm_long, 'trips_long': vkm_long / dis_long,
                      'int': vkm_imex, 'trips_int': vkm_imex / dis_long,
                      'transit': vkm_transit}

            result.update({'short_segments_vkm': vkm_short_segments})

            return result
        else:
            raise KeyError('Segments and shares / loads do not match!')


class GravityModel:
    """
    Create OD trip matrices for personal and freight transport using a gravity model
    """

    def __init__(self, cost_matrix, taz, taz_cn="cntr_code"):
        """

        :param matrix: np.array with shape [len(taz),len(taz),2], cost matrix with travel times [:,:,0] and distances [:,:,1]
        :param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ
        :param taz_cn: name of column with country code in TAZ

        :type matrix: np.array
        :type taz: gpd.GeoDataFrame or pd.DataFrame
        :type taz_cn: str
        """
        self.matrix = cost_matrix
        self.taz = taz
        self.taz_cn = taz_cn

    def get_country_model(self, cn):
        """
        Get TAZ and cost matrix for country
        :param cn: country code
        :return: TAZ for country, cost matrix for country
        :type cn: str
        :rtype: pd.DataFrame or gpd.GeoDataFrame, np.array
        """
        taz_c = self.taz[self.taz[self.taz_cn] == cn].copy()
        taz_c_ids = taz_c['id']
        # matrix from country zones
        mtx_cn = self.matrix[taz_c_ids, :, :]
        # matrix to country zones
        mtx_cn = mtx_cn[:, taz_c_ids, :]
        # set new taz id for taz_c
        new_id = np.c_[taz_c_ids, range(len(taz_c_ids))]
        taz_c['id'] = new_id[:, 1]
        return taz_c, mtx_cn

    def matrix_international(self):
        """
        Adjust cost matrix for international assignment (no inner country trips) by setting costs for inner country relations to 9.e+12

        :return: international cost matrix with travel times and distances
        :rtype: np.array
        """
        countries = self.taz[self.taz_cn].unique()
        mx_int = self.matrix.copy()
        mx_fil = np.zeros(self.matrix.shape)
        for cn in countries:
            # find all IDs for cn and set values between them to 9.e+12
            taz_c_ids = self.taz.loc[self.taz[self.taz_cn]==cn, 'id'].tolist()
            mx_fil[taz_c_ids, :, :] += 1
            mx_fil[:, taz_c_ids, :] += 1
            mx_fil[mx_fil == 1] = 0
        mx_int[mx_fil > 1] = 9.e+12
        return mx_int

    def trip_distribution_pt(self, target, cn=None, taz_pop='population', alpha=0.5, gamma=-2.75, unit_dis=1000, unit_tt = 60, mob_rate=36., occ_rate=1.3):
        """
        Distribute personal transport using a gravity model and an input for total VKT
        Gravity model parameters are given as defaults and were estimated using the German NHTS MiD 2017

        :param target: total VKT
        :param cn: country code; default None means all countries in self.taz are included
        :param taz_pop: name of column with population in taz, default 'population'
        :param alpha: alpha parameter in gravity model, default 0.5
        :param gamma: gamma parameter in gravity model, default -2.75
        :param unit_dis: factor to transform unit in distance matrix to km, default 1000 (suggesting distance is in m)
        :param unit_tt: factor to transform unit in travel time matrix to min, default 60 (suggesting distance is in s)
        :param mob_rate: mobility rate of inhabitants, default 36
        :param occ_rate: average vehicle occupancy
        :return: np.array with OD trip matrix

        :type target: float
        :type cn: str
        :type taz_pop: str
        :type alpha: float
        :type gamma: float
        :type unit_dis: float
        :type unit_tt: float
        :type mob_rate: float
        :type occ_rate: float
        :rtype: np.array
        """
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
        else:
            taz = self.taz
            mtx = self.matrix_international()
            #mtx[mtx==9.e+12] = 0

        # transform units in mtx to min, km
        mtx[:, :, 0] /= unit_tt
        mtx[:, :, 1] /= unit_dis

        # trip generation
        taz['pt_goal'] = taz[taz_pop] * mob_rate

        # trip distribution
        # gravity model: gravity values
        mx_grav = np.zeros((len(taz), len(taz)))
        for o in taz['id']:
            pop_o = float(taz.loc[taz['id'] == o, taz_pop].iloc[0])
            for d in taz['id']:
                pop_d = float(taz.loc[taz['id'] == d, taz_pop].iloc[0])
                mx_grav[o, d] = (pop_o*pop_d)**alpha * mtx[o, d, 0]**gamma
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0
        # choice probabilities and trips
        mx_trips = get_trip_matrix(mx_grav, taz, 'pt_goal', taz_id='id')
        # scaling to match target
        vkt = (mx_trips*mtx[:, :, 1]).sum()
        scale_fac = target / vkt
        print('Scaling factor personal transport for {}: {}'.format(cn, scale_fac))
        mx_trips *= scale_fac

        # transform person trips to vehicle trips
        mx_trips /= occ_rate

        return mx_trips

    def trip_distribution_ft(self, target_trips, target_vkt=None, cn=None, beta=0.00421, unit_tt=60, unit_dis=1000, trips_key=''):
        """
        Distribute freight transport using a gravity model and an input for total trips and VKT
        Gravity model parameters are given as defaults and were estimated using microscopic truck data for Germany
        :param target_trips: float or dict, total trips or if cn is None dict with total trips per country and segment in the form of {cn: {segment: target_trips}}
        :param target_vkt: optional: total VKT used for comparing VKT from OD matrix and goal, default None
        :param cn: country code; default None means all countries in self.taz are included
        :param beta: gamma parameter in gravity model, default 0.00421
        :param unit_dis: factor to transform unit in distance matrix to km, default 1000 (suggesting distance is in m)
        :param unit_tt: factor to transform unit in travel time matrix to min, default 60 (suggesting distance is in s)
        :param trips_key: key for segment if cn is None target_trips is dict
        :return: np.array with OD trip matrix

        :type target_trips: float or dict
        :type target_vkt: float
        :type cn: None or str
        :type beta: float
        :type unit_dis: float
        :type unit_tt: float
        :type trips_key: str
        :rtype: np.array
        """
        t_targ = type(target_trips)
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
            index_col = 'index_nat'
            if (t_targ != float) & (t_targ != int):
                raise ValueError('target_trips has to be float / int if cn is not None, is {}'.format(t_targ))
        else:
            taz = self.taz
            mtx = self.matrix_international()
            index_col = 'index_int'
            if t_targ != dict:
                raise ValueError('target_trips has to be dict if cn is None, is {}'.format(t_targ))

        # transform units in mtx to min, km
        mtx[:, :, 0] /= unit_tt
        mtx[:, :, 1] /= unit_dis

        # trip generation
        if cn is not None:
            taz['ft_goal'] = taz[index_col]/taz[index_col].sum() * target_trips
        else:
            # trip generation per country using target values
            taz['ft_goal'] = 0
            for c in taz[self.taz_cn].unique():
                if c in target_trips.keys():
                    taz.loc[taz[self.taz_cn]==c, 'ft_goal'] = taz.loc[taz[self.taz_cn]==c, index_col]/taz.loc[taz[self.taz_cn]==c, index_col].sum() * target_trips[c][trips_key]

        # trip distribution
        # gravity model: gravity values
        mx_grav = np.zeros((len(taz), len(taz)))
        for o in taz['id']:
            ind_o = float(taz.loc[taz['id'] == o, index_col].iloc[0])
            for d in taz['id']:
                ind_d = float(taz.loc[taz['id'] == d, index_col].iloc[0])
                mx_grav[o, d] = ind_o * ind_d * np.exp(-beta*mtx[o, d, 1])
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0

        # get trips
        mx_trips = get_trip_matrix(mx_grav, taz, 'ft_goal', taz_id='id')

        if target_vkt is not None:
            # determine relation to target
            vkt = (mx_trips * mtx[:, :, 1]).sum()
            scale_fac = target_vkt / vkt
            print('Freight transport gravity model result for {}: Relation to target vkm {}'.format(cn, scale_fac))

        return mx_trips


class IntraZonal:
    """
    Distribute intrazonal trips on road network
    """

    def __init__(self, taz, net, from_node='from_node', to_node='to_node', link_id='ultimo_id'):
        """

        :param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ
        :param net: gpd.GeoDataFrame or pd.DataFrame with road network
        :param from_node: column name for from node in net, default 'from_node'
        :param to_node: column name for to node in net, default 'to_node'
        :param link_id: column name for unique id of the elements in net, default 'ultimo_id'
        :type taz: gpd.GeoDataFrame or pd.DataFrame
        :type net: gpd.GeoDataFrame or pd.DataFrame
        :type from_node: str
        :type to_node: str
        :type link_id: str
        """
        self.taz = taz
        self.net = net
        self.from_ = from_node
        self.to_ = to_node
        self.link_id = link_id

    def road_type_weighted_single(self, target, weights=None, veh_types=None, taz_id='nuts_id',
                                  index='population', sub_len='length_sub', net_type='type', occ_rate=1.):
        """
        Distribute total VKT per TAZ and assign loads to roads within this TAZ, weighted by road type

        :param target: target values for total VKT per vehicle type in the form of {veh_type: target_vkt}
        :param weights: weights of different road types; default None leads to pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0.75, 1.5, 3.5, 3.5]})
        :param veh_types: names of vehicle types to be considered; default None leads to ['car']
        :param taz_id: name of column with taz id in self.taz and self.net, defaults to 'nuts_id'
        :param index: name of column with attraction index to be used for trip generation in self.taz, defaults to 'population'
        :param sub_len: name of column with length of subordinate network in self.taz, defaults to 'length'
        :param net_type: name of column with road type in self.net, defaults to 'type'
        :param occ_rate: average vehicle occupancy (applied for personal transport, i.e. if veh_type == 'car')
        :return: gpd.GeoDataFrame or pd.DataFrame for network with traffic loads, taz with subordinate network VKT
        :type target: dict
        :type weights: None or pd.DataFrame
        :type veh_types: list
        :type taz_id: str
        :type index: str
        :type sub_len: str
        :type net_type: str
        :type occ_rate: float
        :rtype: gpd.GeoDataFrame or pd.DataFrame
        """
        if weights is None:
            weights = pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0.75, 1.5, 3.5, 3.5]})
        if veh_types is None:
            veh_types = ['car']

        # subordinate network length from TAZ
        sub = self.taz[[taz_id, sub_len]].copy()
        sub.rename(columns={sub_len: 'length'}, inplace=True)
        sub[net_type] = 0
        sub['weighted_length'] = sub['length']*float(weights.loc[weights[net_type]==0, 'weight'].iloc[0])
        # concat with network length per taz from net
        len_per_type = self.net.groupby([taz_id, net_type])['length'].sum().reset_index()
        len_per_type = pd.concat([len_per_type, sub[[taz_id, net_type, 'length']]])

        # weighted length per net type and taz
        net = self.net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length']*net['weight']

        w_len_per_type = net.groupby([taz_id, net_type])['weighted_length'].sum().reset_index()
        w_len_per_type = pd.concat([w_len_per_type, sub[[taz_id, net_type, 'weighted_length']]])

        net.drop(columns=['weight', 'weighted_length'], inplace=True)

        # weighted length per taz
        w_len_per_taz = w_len_per_type.groupby(taz_id)['weighted_length'].sum().reset_index()
        w_len_per_taz.rename(columns={'weighted_length':'sum_weighted_length'}, inplace=True)

        # merge len and w_len
        len_per_type = len_per_type.merge(w_len_per_type, on=[taz_id, net_type], how='left')

        # merge weighted length sum to taz
        taz_2 = self.taz.merge(w_len_per_taz, on=taz_id, how='left')

        del w_len_per_taz, w_len_per_type

        taz_result = self.taz.copy()
        for veh_type in veh_types:

            vkt_p_p = target[veh_type] / self.taz[index].sum()

            if veh_type == 'car':
                # apply occupancy rate
                vkt_p_p /= occ_rate

            # calculate vkt per taz
            taz_2['goal_taz_{}'.format(veh_type)] = taz_2[index] * vkt_p_p
            taz_2['goal_norm_{}'.format(veh_type)] = taz_2['goal_taz_{}'.format(veh_type)] / taz_2['sum_weighted_length']

            # merge to len_per_type
            len_per_type = len_per_type.merge(taz_2[[taz_id, 'goal_taz_{}'.format(veh_type), 'goal_norm_{}'.format(veh_type)]], on=taz_id, how='left')
            len_per_type['vkt'] = len_per_type['goal_norm_{}'.format(veh_type)] * len_per_type['weighted_length']
            len_per_type['trips'] = len_per_type['vkt'] / len_per_type['length']

            # add to network links and taz
            net = net.merge(len_per_type[[taz_id, net_type, 'trips']], on=[taz_id, net_type], how='left')
            net.loc[net['trips'].isna(), 'trips'] = 0
            net.rename(columns={'trips': '{}_short'.format(veh_type)}, inplace=True)

            sec = len_per_type[len_per_type[net_type] == 0]
            taz_result = taz_result.merge(sec[[taz_id, 'vkt']], on=taz_id, how='left')
            taz_result.rename(columns={'vkt': '{}_sub'.format(veh_type)}, inplace=True)

        return net, taz_result

    def road_type_weighted_multiple(self, target, matrix_dis, veh_types=None, weights=None, taz_id='nuts_id', taz_mx_id='id',
                                    index='index_nat', sub_len='length_sub', net_type='type', distance=55, cell_size=500,
                                    fac_cell=1.5, occ_rate=1.):
        """
        Distribute total VKT per TAZ and assign loads to roads within this TAZ and surrounding TAZ, weighted by road type

        :param target: target values for total VKT per vehicle type in the form of {veh_type: target_vkt}
        :param matrix_dis: distance matrix between TAZ (in km)
        :param veh_types: names of vehicle types to be considered; default None leads to ['lcv', 'mft', 'hft']
        :param weights: weights of different road types; default None leads to pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0, 5, 1.5, 0.5]})
        :param taz_id: name of column with taz id in self.taz and self.net, defaults to 'nuts_id'
        :param taz_mx_id: name of column with taz id relating to position in matrix_dis in self.taz, defaults to 'id'
        :param index: name of column with attraction index to be used for trip generation in self.taz, defaults to 'index_nat'
        :param sub_len: name of column with length of subordinate network in self.taz, defaults to 'length'
        :param net_type: name of column with road type in self.net, defaults to 'type'
        :param distance: max distance between TAZ to be included in km, defaults to 55km
        :param cell_size: max cell size of TAZ to force inclusion of surrounding TAZ, defaults to 500km²
        :param fac_cell: factor to apply to main TAZ during distribution, defaults to 1.5
        :param occ_rate: average vehicle occupancy (applied for personal transport, i.e. if veh_type == 'car')
        :return: gpd.GeoDataFrame or pd.DataFrame for network with traffic loads, taz with subordinate network VKT
        :type target: dict
        :type matrix_dis: np.array
        :type veh_types: list
        :type weights: pd.DataFrame
        :type taz_id: str
        :type taz_mx_id: str
        :type index: str
        :type sub_len: str
        :type net_type: str
        :type distance: float
        :type cell_size: float
        :type fac_cell: float
        :type occ_rate: float
        :rtype: gpd.GeoDataFrame or pd.DataFrame
        """

        if weights is None:
            weights = pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0, 5, 1.5, 0.5]})
        if veh_types is None:
            veh_types = ['lcv', 'mft', 'hft']

        net = self.net.copy()
        # target per taz
        taz_2 = self.taz[[taz_id, taz_mx_id, index, 'area', sub_len]].copy()
        for veh_type in veh_types:
            taz_2[veh_type] = taz_2[index] / taz_2[index].sum() * target[veh_type]
            if veh_type == 'car':
                # apply occupancy rate
                taz_2[veh_type] /= occ_rate
            net['{}_short'.format(veh_type)] = 0
            taz_2['{}_sub'.format(veh_type)] = 0

        # weighted length per net type and taz
        net = net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length'] * net['weight']
        taz_2['weighted_sub'] = taz_2[sub_len]*float(weights.loc[weights[net_type]==0, 'weight'].iloc[0])

        # iterate over taz
        for t in taz_2[taz_id]:
            # adjust link weight for roads within taz
            net['weighted_length2'] = net['weighted_length']
            net['weight2'] = net['weight']
            taz_2['weighted_sub2'] = taz_2['weighted_sub']
            taz_2['weight2'] = float(weights.loc[weights[net_type]==0, 'weight'].iloc[0])
            net.loc[net[taz_id] == t, 'weighted_length2'] = net.loc[net[taz_id] == t, 'weighted_length2'] * fac_cell
            net.loc[net[taz_id] == t, 'weight2'] = net.loc[net[taz_id] == t, 'weight2'] *fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] = taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] * fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weight2'] = taz_2.loc[taz_2[taz_id] == t, 'weight2'] * fac_cell

            # list of taz within distance
            t_id = int(taz_2.loc[taz_2[taz_id] == t, taz_mx_id].iloc[0])
            matrix_t = matrix_dis[t_id, :]
            sur_ids = np.where(matrix_t <= distance)[0].tolist()

            # for small taz: add surrounding TAZ within higher distances if no sur found for distance
            if (len(sur_ids) == 0) & (float(taz_2.loc[taz_2[taz_id] == t, 'area'].iloc[0]) < cell_size):
                sur_ids = np.where(matrix_t == matrix_t.min())[0].tolist()

            # add id of t to sur ids (all relevant cells) and get taz_ids
            sur_ids.append(t_id)
            sur_ids_taz = taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), taz_id].tolist()

            # aggregate weighted network length in sur_ids
            agg_w_len = net.loc[net[taz_id].isin(sur_ids_taz), 'weighted_length2'].sum()
            agg_w_len += taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), 'weighted_sub2'].sum()

            # calculate trips per km for veh_types
            for veh_type in veh_types:
                taz_vkt = float(taz_2.loc[taz_2[taz_id] == t, veh_type].iloc[0])
                taz_veh = taz_vkt / agg_w_len
                net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] = net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] + net.loc[net[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh
                taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), '{}_sub'.format(veh_type)] += taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh * taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), sub_len]

        net.drop(columns=['weight', 'weight2', 'weighted_length', 'weighted_length2'], inplace=True)

        cols_keep = ['{}_sub'.format(veh_type) for veh_type in veh_types]
        cols_keep.extend([taz_id])
        taz_2 = taz_2[cols_keep]

        return net, taz_2
