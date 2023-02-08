# =========================================================
# DistributeTraffic.py
# @author Nina Thomsen
# @date 02.02.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Distribute traffic between TAZ and create OD Matrices;
#        Assign inner cell traffic to network links
# =========================================================

import pandas as pd
import numpy as np


def pt_shares_int(a, b):
    """
    Determine shares of international transport for personal transport for a country based on area and share of land border

    @param a: area of country
    @param b: share of land border of country
    @return: dict with international transport shares
    """
    # transform area
    a = (np.log((a-(-15000))/6.59)/0.437)-22.69

    x1 = (0.3 / (1 + ((0.3 / 0.095 - 1) * np.exp(a * 0.3 * -(-1.05))))) * b
    x2 = max(0.02, (0.35 + -0.04 * a) * b**0.5)
    x3 = min(0.95, max(0.1, (0.3 + 0.05 * a)))

    return {'pt_share_foreign': x1, 'pt_share_transit': x2, 'pt_ratio_int': x3}


def get_trip_matrix(mx_grav, taz, goal_col):
    """
    Get trip matrix based on total trips per taz and gravity model matrix

    @param mx_grav: np.array, matrix from gravity model
    @param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ; needs id column to identify coordinate in mx_grav
    @param goal_col: str, name of column in taz with total trips per taz
    @return: np.array, matrix with trips between taz
    """
    # choice probabilities and trips
    mx_trips = np.zeros(mx_grav.shape)
    for t in taz['id']:
        goal = float(taz.loc[taz['id'] == t, goal_col])
        # get probabilities for trips starting at t
        prob_a = mx_grav[t, :] / mx_grav[t, :].sum()
        prob_b = mx_grav[:, t] / mx_grav[:, t].sum()
        mx_trips[t, :] += prob_a * (goal / 2)
        mx_trips[:, t] += prob_b * (goal / 2)
    # fill NA
    mx_trips[np.isnan(mx_trips)] = 0
    return mx_trips


class TargetValues:
    """
    Determine shares of international and national transport for a country and calculate corresponding total trips and vehicle kilometers travelled (VKT)
    """

    def __init__(self, country_layer, cn_col='cntr_code'):
        """

        @param country_layer: pd.DataFrame or gpd.GeoDataFrame, countries with target values for raod transport volumes
        @param cn_col: str, name of column with country code in country layer, default 'cntr_code'
        """
        self.country_layer = country_layer
        self.cn_col = cn_col

    def targets_personal_transport(self, cn, target_col='car_pkm', unit=1.e+6):
        """
        Calculate VKT for the segments short and long distance (national) and international personal transport for a country

        @param cn: str, country code
        @param target_col: str, name of column in self.country_layer with the total transport volume in person kilometers travelled, default 'car_pkm'
        @param unit: float or int, scaling factor for total transport volume, default 1.e+6
        @return: dict, vkt per segment
        """
        target = float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, target_col])
        shares_int = pt_shares_int(a=float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, 'area']),
                                   b=float(self.country_layer.loc[self.country_layer[self.cn_col]==cn, 'border_share']))

        f_in_c = target * shares_int['pt_share_foreign']
        int1 = target * shares_int['pt_share_foreign'] * shares_int['pt_ratio_int']
        inner = target - f_in_c - int1

        short_ = 0.65 * (inner + int1)
        long_ = 0.35 * (inner + int1) - int1

        trans_ = shares_int['pt_share_transit'] * f_in_c
        int_ = (f_in_c - trans_) + int1

        return {"short": unit * short_, "long": unit * long_, "int": unit * int_, "trans": unit * trans_}

    def targets_freight_transport(self, cn, target_col='car_pkm', segments=None, shares_tkm=None, loads_t=None, unit=1.e+6):
        """
        Calculate VKT and trips for the segments short and long distance (national) and international freight transport for a country

        @param cn: str, country code
        @param target_col: str, name of column in self.country_layer with the total transport volume in tkm
        @param segments: None or list, names of freight transport segments; default None refers to ['lcv', 'mft', 'hft']
        @param shares_tkm: None or dict, shares of different freight transport segments of total tkm; default None refers to shares {'lcv': 0.08, 'mft': 0.09, 'hft': 0.83}
        @param loads_t: None or dict, average load in t per freight transport segment; default None refers to load {'lcv': 0.5, 'mft': 3, 'hft': 10}
        @param unit: float or int, scaling factor for total transport volume, default 1.e+6
        @return: dict, vkt and trips per segment
        """

        if shares_tkm is None:
            shares_tkm = {'lcv': 0.08, 'mft': 0.09, 'hft': 0.83}
        if loads_t is None:
            loads_t = {'lcv': 0.5, 'mft': 3, 'hft': 10}

        # check if segments, loads and shares match
        if (sorted(segments) == sorted(shares_tkm.keys())) & (sorted(segments) == sorted(loads_t.keys())):

            target = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, target_col])

            cat_dict = {0: {'int_imex': 0.03, 'int_transit': 0.00, 'nat': 0.97},
                        1: {'int_imex': 0.15, 'int_transit': 0.02, 'nat': 0.83},
                        2: {'int_imex': 0.15, 'int_transit': 0.05, 'nat': 0.80},
                        3: {'int_imex': 0.20, 'int_transit': 0.10, 'nat': 0.70},
                        4: {'int_imex': 0.00, 'int_transit': 0.00, 'nat': 1.00}}

            border_share = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'border_share'])
            neighbors = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'neighbors'])
            border_crossings = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'border_crossings_count'])

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
            area = float(self.country_layer.loc[self.country_layer[self.cn_col] == cn, 'area'])
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

            # transform tkm t vkm by using average tkm shares and average load
            lcv = shares_tkm['lcv'] * tkm / tkm_short
            mft = shares_tkm['mft'] * tkm / tkm_short
            hft = (shares_tkm['hft'] * tkm - (tkm - tkm_short)) / tkm_short

            # vkm
            vkm_lcv = (tkm_short * lcv) / loads_t['lcv']
            vkm_mft = (tkm_short * mft) / loads_t['mft']
            vkm_hft = (tkm_short * hft) / loads_t['hft']
            vkm_short = vkm_lcv + vkm_mft + vkm_hft
            vkm_long = tkm_long / 10
            vkm_imex = tkm_imex / 10

            return {"lcv": unit * vkm_lcv,"mft": unit * vkm_mft,"hft": unit * vkm_hft,"trips_short": vkm_short / dis_short,
                    "long": unit * vkm_long,"trips_long": unit * (vkm_long/dis_long),"int": unit * vkm_imex,
                    "trips_int": unit * (vkm_imex/dis_long),"trans": unit * (tkm_transit/loads_t['hft']), 'trips_trans': unit*((tkm_transit/loads_t['hft'])/ dis_transit)}
        else:
            raise KeyError('Segments and shares / loads do not match!')


class GravityModel:
    """
    Create OD trip matrices for personal and freight transport using a gravity model
    """

    def __init__(self, cost_matrix, taz, taz_cn="cntr_code"):
        """

        @param matrix: np.array with shape [len(taz),len(taz),2], cost matrix with travel times [:,:,0] and distances [:,:,1]
        @param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ
        @param taz_cn: str, name of column with country code in TAZ
        """
        self.matrix = cost_matrix
        self.taz = taz
        self.taz_cn = taz_cn

    def get_country_model(self, cn):
        """
        Get TAZ and cost matrix for country
        @param cn: str, country code
        @return: TAZ for country (pd.DataFrame or gpd.GeoDataFrame), cost matrix for country (np.array)
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

        @return: international cost matrix with travel times and distances
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

    def trip_distribution_pt(self, target, cn=None, taz_pop='population', alpha=0.5, gamma=-2.75, unit_dis=1000, mob_rate=36.):
        """
        Distribute personal transport using a gravity model and an input for total VKT
        Gravity model parameters are given as defaults and were estimated using the German NHTS MiD 2017

        @param target: float, total VKT
        @param cn: None or str, country code; default None means all countries in self.taz are included
        @param taz_pop: str, name of column with population in taz, default 'population'
        @param alpha: float, alpha parameter in gravity model, default 0.5
        @param gamma: float, gamma parameter in gravity model, default -2.75
        @param unit_dis: int or float, factor to transform unit in distance matrix to km, default 1000 (suggesting distance is in m)
        @param mob_rate: float, mobility rate of inhabitants, default 36
        @return: np.array with OD trip matrix
        """
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
        else:
            taz = self.taz
            mtx = self.matrix_international()
            mtx[mtx==9.e+12] = 0
        # trip generation
        taz['pt_goal'] = taz[taz_pop] * mob_rate

        # trip distribution
        # gravity model: gravity values
        mx_grav = np.zeros((len(taz), len(taz)))
        for o in taz['id']:
            pop_o = float(taz.loc[taz['id'] == o, taz_pop])
            for d in taz['id']:
                pop_d = float(taz.loc[taz['id'] == d, taz_pop])
                mx_grav[o, d] = (pop_o*pop_d)**alpha * mtx[o, d, 0]**gamma
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0
        # choice probabilities and trips
        mx_trips = get_trip_matrix(mx_grav, taz, 'pt_goal')
        # scaling to match target
        vkt = (mx_trips*mtx[:, :, 1]).sum()/unit_dis
        scale_fac = target / vkt
        print('Scaling factor for {}: {}'.format(cn, scale_fac))
        mx_trips *= scale_fac

        return mx_trips

    def trip_distribution_ft(self, target_trips, target_vkt=None, cn=None, beta=0.00421, unit_dis=1000, trips_key=''):
        """
        Distribute freight transport using a gravity model and an input for total trips and VKT
        Gravity model parameters are given as defaults and were estimated using microscopic truck data for Germany
        @param target_trips: float or dict, total trips or if cn is None dict with total trips per country and segment in the form of {cn: {segment: target_trips}}
        @param target_vkt: None or float, optional: total VKT used for comparing VKT from OD matrix and goal, default None
        @param cn: None or str, country code; default None means all countries in self.taz are included
        @param beta: float, gamma parameter in gravity model, default 0.00421
        @param unit_dis: int or float, factor to transform unit in distance matrix to km, default 1000 (suggesting distance is in m)
        @param trips_key: str, key for segment if cn is None target_trips is dict
        @return: np.array with OD trip matrix
        """
        # TODO: add check for target_trips (dict in correct form or value, depending on cn)
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
            index_col = 'index_nat'
        else:
            taz = self.taz
            mtx = self.matrix_international()
            index_col = 'index_int'
        if unit_dis == 1000:
            # assure distance matrix is in km
            mtx[:, :, 1] /= 1000

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
            ind_o = float(taz.loc[taz['id'] == o, index_col])
            for d in taz['id']:
                ind_d = float(taz.loc[taz['id'] == d, index_col])
                mx_grav[o, d] = ind_o * ind_d * np.exp(-beta*mtx[o, d, 1])
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0

        # get trips
        mx_trips = get_trip_matrix(mx_grav, taz, 'ft_goal')

        if target_vkt is not None:
            # determine relation to target
            vkt = (mx_trips * mtx[:, :, 1]).sum() / unit_dis
            scale_fac = target_vkt / vkt
            print('Relation to target vkm for {}: {}'.format(cn, scale_fac))

        return mx_trips


class IntraZonal:
    """
    Distribute intrazonal trips on road network
    """

    def __init__(self, taz, net):
        """

        @param taz: gpd.GeoDataFrame or pd.DataFrame with TAZ
        @param net: gpd.GeoDataFrame or pd.DataFrame with road network
        """
        self.taz = taz
        self.net = net

    def road_type_weighted_single(self, target, weights=None, veh_types=None, taz_id='nuts_id',
                                  index='population', sub_len='length', net_type='type', urban=False):
        """
        Distribute total VKT per TAZ and assign loads to roads within this TAZ, weighted by road type
        @param target: dict, target values for total VKT per vehicle type in the form of {veh_type: target_vkt}
        @param weights: None or pd.DataFrame, weights of different road types; default None leads to pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0.75, 1.5, 3.5, 3.5]})
        @param veh_types: None or list, names of vehicle types to be considered; default None leads to ['car']
        @param taz_id: str, name of column with taz id in self.taz and self.net, defaults to 'nuts_id'
        @param index: str, name of column with attraction index to be used for trip generation in self.taz, defaults to 'population'
        @param sub_len: str, name of column with length of subordinate network in self.taz, defaults to 'length'
        @param net_type: str, name of column with road type in self.net, defaults to 'type'
        @param urban: not applicable yet (bool, if True, additional weight is placed on roads based on their location within or outside of urban areas, default False)
        @return: gpd.GeoDataFrame or pd.DataFrame for network with traffic loads, taz with subordinate network VKT
        """
        if weights is None:
            weights = pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0.75, 1.5, 3.5, 3.5]})
        if veh_types is None:
            veh_types = ['car']

        # subordinate network length from TAZ
        sub = self.taz[[taz_id, sub_len]].copy()
        sub.rename(columns={sub_len: 'length'}, inplace=True)
        sub[net_type] = 0
        sub['weighted_length'] = sub['length']*float(weights.loc[weights[net_type]==0, 'weight'])
        # concat with network length per taz from net
        len_per_type = self.net.groupby([taz_id, net_type])['length'].sum().reset_index()
        len_per_type = pd.concat([len_per_type, sub[[taz_id, net_type, 'length']]])

        # weighted length per net type and taz
        net = self.net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length']*net['weight']

        # todo: adjust urban weights
        if urban:
            net.loc[net['urban']==1, 'weighted_length'] = net.loc[net['urban']==1, 'weighted_length']*0.5

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
                                    index='index_nat', sub_len='length', net_type='type', urban=None, distance=55, cell_size=500,
                                    fac_cell=1.5):
        """
        Distribute total VKT per TAZ and assign loads to roads within this TAZ and surrounding TAZ, weighted by road type
        @param target: dict, target values for total VKT per vehicle type in the form of {veh_type: target_vkt}
        @param matrix_dis: np.array, distance matrix between TAZ (in km)
        @param veh_types: None or list, names of vehicle types to be considered; default None leads to ['lcv', 'mft', 'hft']
        @param weights: None or pd.DataFrame, weights of different road types; default None leads to pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0, 5, 1.5, 0.5]})
        @param taz_id: str, name of column with taz id in self.taz and self.net, defaults to 'nuts_id'
        @param taz_mx_id: str, name of column with taz id relating to position in matrix_dis in self.taz, defaults to 'id'
        @param index: str, name of column with attraction index to be used for trip generation in self.taz, defaults to 'index_nat'
        @param sub_len: str, name of column with length of subordinate network in self.taz, defaults to 'length'
        @param net_type: str, name of column with road type in self.net, defaults to 'type'
        @param urban: not applicable yet (bool, if True, additional weight is placed on roads based on their location within or outside of urban areas, default False)
        @param distance: float, max distance between TAZ to be included in km, defaults to 55km
        @param cell_size: float, max cell size of TAZ to force inclusion of surrounding TAZ, defaults to 500km²
        @param fac_cell: float, factor to apply to main TAZ during distribution, defaults to 1.5
        @return: gpd.GeoDataFrame or pd.DataFrame for network with traffic loads, taz with subordinate network VKT
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
            net['{}_short'.format(veh_type)] = 0
            taz_2['{}_sub'.format(veh_type)] = 0

        # weighted length per net type and taz
        net = net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length'] * net['weight']
        taz_2['weighted_sub'] = taz_2[sub_len]*float(weights.loc[weights[net_type]==0, 'weight'])

        # todo: adjust urban weights
        if urban:
            net.loc[net['urban']==1, 'weighted_length'] = net.loc[net['urban']==1, 'weighted_length']*0.5

        # iterate over taz
        for t in taz_2[taz_id]:
            # adjust link weight for roads within taz
            net['weighted_length2'] = net['weighted_length']
            net['weight2'] = net['weight']
            taz_2['weighted_sub2'] = taz_2['weighted_sub']
            taz_2['weight2'] = float(weights.loc[weights[net_type]==0, 'weight'])
            net.loc[net[taz_id] == t, 'weighted_length2'] = net.loc[net[taz_id] == t, 'weighted_length2'] * fac_cell
            net.loc[net[taz_id] == t, 'weight2'] = net.loc[net[taz_id] == t, 'weight2'] *fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] = taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] * fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weight2'] = taz_2.loc[taz_2[taz_id] == t, 'weight2'] * fac_cell

            # list of taz within distance
            t_id = int(taz_2.loc[taz_2[taz_id] == t, taz_mx_id])
            matrix_t = matrix_dis[t_id, :]
            sur_ids = np.where(matrix_t <= distance)[0].tolist()

            # for small taz: add surrounding TAZ within higher distances if no sur found for distance
            if (len(sur_ids)==0) & (float(taz_2.loc[taz_2[taz_id]==t, 'area'])<cell_size):
                sur_ids = np.where(matrix_t == matrix_t.min())[0].tolist()

            # add id of t to sur ids (all relevant cells) and get taz_ids
            sur_ids.append(t_id)
            sur_ids_taz = taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), taz_id].tolist()

            # aggregate weighted network length in sur_ids
            agg_w_len = net.loc[net[taz_id].isin(sur_ids_taz), 'weighted_length2'].sum()
            agg_w_len += taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), 'weighted_sub2'].sum()

            # calculate trips per km for veh_types
            for veh_type in veh_types:
                taz_vkt = float(taz_2.loc[taz_2[taz_id] == t, veh_type])
                taz_veh = taz_vkt / agg_w_len
                net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] = net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] + net.loc[net[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh
                taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), '{}_sub'.format(veh_type)] += taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh * taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), sub_len]

        net.drop(columns=['weight', 'weight2', 'weighted_length', 'weighted_length2'], inplace=True)
        taz_2.drop(columns=['weight2', 'weighted_sub', 'weighted_sub2'], inplace=True)

        return net, taz_2







