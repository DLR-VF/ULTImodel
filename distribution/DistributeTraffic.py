# =========================================================
# DistributeTraffic.py
# @author Nina Thomsen
# @date 02.02.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Distribute traffic between TAZ and create OD Matrix
# =========================================================

from datetime import datetime

import geopandas as gpd
import pandas as pd
import numpy as np

from tqdm import tqdm


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
    x2 = max(0.02, (0.35 + -0.04 * a) * b ^ 0.5)
    x3 = min(0.95, max(0.1, (0.3 + 0.05 * a)))

    return {'pt_share_foreign': x1, 'pt_share_transit': x2, 'pt_ratio_int': x3}


def get_trip_matrix(mx_grav, taz, goal_col):
    # choice probabilities and trips
    mx_trips = np.zeros(mx_grav.shape)
    for t in taz['id']:
        goal = taz.loc[taz['id'] == t, goal_col]
        # get probabilities for trips starting at t
        prob_a = mx_grav[t, :] / mx_grav[t, :].sum()
        prob_b = mx_grav[:, t] / mx_grav[:, t].sum()
        mx_trips[t, :] += prob_a * (goal / 2)
        mx_trips[:, t] += prob_b * (goal / 2)
    return mx_trips


class TargetValues:

    def __init__(self, country_layer, cn_col='cntr_code'):
        self.country_layer = country_layer
        self.cn_col = cn_col

    def targets_personal_transport(self, cn, target_col='car_pkm', unit=1000000.):
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

    def targets_freight_transport(self, cn, target_col='car_pkm', shares_tkm=None, loads_t=None, unit=1000000):

        if shares_tkm is None:
            shares_tkm = {'lcv': 0.08, 'mft': 0.09, 'hft': 0.83}
        if loads_t is None:
            loads_t = {'lcv': 0.5, 'mft': 3, 'hft': 10}

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

        return {"lcv": unit * vkm_lcv,"mft": unit * vkm_mft,"hft": unit * vkm_hft,"trips3": vkm_short / dis_short,
                "long": unit * vkm_long,"trips1": unit * (vkm_long/dis_long),"int": unit * vkm_imex,
                "trips2": unit * (vkm_imex/dis_long),"trans": unit * (tkm_transit/10)}


class GravityModel:

    def __init__(self, matrix, taz, taz_cn):
        self.matrix = matrix
        self.taz = taz
        self.taz_cn = taz_cn

    def get_country_model(self, cn):
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

    def trip_distribution_pt(self, target, cn=None, taz_pop='population', alpha=0.5, gamma=-2.75):
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
        else:
            taz, mtx = self.taz, self.matrix
        # trip generation
        trips_p_person = 36
        taz['pt_goal'] = taz[taz_pop] * trips_p_person

        # trip distribution
        # gravity model: gravity values
        mx_grav = np.zeros((len(taz), len(taz)))
        for o in taz['id']:
            pop_o = taz.loc[taz['id'] == o, taz_pop]
            for d in taz['id']:
                pop_d = taz.loc[taz['id'] == d, taz_pop]
                mx_grav[o, d] = (pop_o*pop_d)**alpha * mtx[o, d, 0]**gamma
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0
        # choice probabilities and trips
        mx_trips = get_trip_matrix(mx_grav, taz, 'pt_goal')
        '''for t in taz['id']:
            goal = taz.loc[taz['id'] == t, 'pt_goal']
            # get probabilities for trips starting at t
            prob_a = mx_grav[t, :] / mx_grav[t, :].sum()
            prob_b = mx_grav[:, t] / mx_grav[:, t].sum()
            mx_trips[t, :] += prob_a * (goal / 2)
            mx_trips[:, t] += prob_b * (goal / 2)'''

        # scaling to match target
        # TODO units km or m?
        vkt = mx_trips*mtx[:, :, 1]
        scale_fac = target / vkt
        print('Scaling factor for {}: {}'.format(cn, scale_fac))
        mx_trips *= scale_fac

        # determine cold starts per cell
        cold_starts = mx_trips.sum()

        return mx_trips, cold_starts

    def trip_distribution_ft(self, target_vkm, target_trips, cn=None, beta=0.00421):
        if cn is not None:
            taz, mtx = self.get_country_model(cn=cn)
            index_col = 'index_nat'
        else:
            taz, mtx = self.taz, self.matrix
            index_col = 'index_int'

        # trip generation
        taz['ft_goal'] = taz[index_col]/taz[index_col].sum() * target_trips

        # trip distribution
        # gravity model: gravity values
        mx_grav = np.zeros((len(taz), len(taz)))
        for o in taz['id']:
            ind_o = taz.loc[taz['id'] == o, index_col]
            for d in taz['id']:
                ind_d = taz.loc[taz['id'] == d, index_col]
                mx_grav[o, d] = ind_o * ind_d * np.exp(-beta*mtx[o, d, 1])
        # fill diagonal (no trips) and NaN
        np.fill_diagonal(mx_grav, 0)
        mx_grav[np.isnan(mx_grav)] = 0

        return get_trip_matrix(mx_grav, taz, 'ft_goal')


class IntraZonal:

    def __init__(self, taz, net):
        self.taz = taz
        self.net = net

    def road_type_weighted_single(self, target, segment_suffix='car', weights=None, taz_id='nuts_id', taz_pop='population', sub_len='length', net_type='type', urban=None):
        if weights is None:
            weights = pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0.75, 1.5, 3.5, 3.5]})

        vkt_p_p = target/self.taz[taz_pop].sum()

        # subordinate network length from TAZ
        sub = self.taz[[taz_id, sub_len]]
        sub.rename(columns={sub_len: 'length'}, inplace=True)
        sub[net_type] = 0
        sub['weighted_length'] = sub['length']*weights.loc[weights[net_type]==0, 'weight']
        # concat with network length per taz from net
        len_per_type = self.net.groupby([taz_id, net_type])['length'].sum().reset_index()
        len_per_type = pd.concat([len_per_type, sec[[taz_id, net_type, 'length']]])

        # weighted length per net type and taz
        net = self.net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length']*net['weight']

        # todo: adjust urban weights
        if urban:
            net.loc[net['urban']==1, 'weighted_length'] = net.loc[net['urban']==1, 'weighted_length']*0.5

        w_len_per_type = net.groupby([taz_id, net_type])['weighted_length'].sum().reset_index()
        w_len_per_type = pd.concat([w_len_per_type, sec[[taz_id, net_type, 'weighted_length']]])

        net.drop(columns=['weight', 'weighted_length'], inplace=True)

        # weighted length per taz
        w_len_per_taz = w_len_per_type.groupby(taz_id)['weighted_length'].sum().reset_index()
        w_len_per_taz.rename(columns={'weighted_length':'sum_weighted_length'}, inplace=True)

        # merge len and w_len
        len_per_type = len_per_type.merge(w_len_per_type, on=[taz_id, net_type], how='left')

        # merge weighted length sum to taz
        taz_2 = self.taz.merge(w_len_per_taz, on=taz_id, how='left')

        del w_len_per_taz, w_len_per_type

        # calculate vkt per taz
        taz_2['goal_taz'] = taz_2[taz_pop]*vkt_p_p
        taz_2['goal_norm'] = taz_2['goal_taz'] / taz_2['sum_weighted_length']

        # merge to len_per_type
        len_per_type = len_per_type.merge(taz_2[[taz_id, 'goal_taz', 'goal_norm']], on=taz_id, how='left')
        len_per_type['vkt'] = len_per_type['goal_norm'] * len_per_type['weighted_length']
        len_per_type['trips'] = len_per_type['goal_vkt'] / len_per_type['length']

        # add to network links and taz
        net = net.merge(len_per_type[[taz_id, net_type, 'trips']], on=[taz_id, type], how='left')
        net.loc[net['trips'].isna(), 'trips'] = 0
        net.rename(columns={'trips': '{}_short'.format(segment_suffix)}, inplace=True)

        sec = len_per_type[len_per_type[net_type] == 0]
        taz_result = self.taz.merge(sec[[taz_id, 'vkt']], on=taz_id, how='left')
        taz_result.rename(columns={'vkt': 'vkt_subordinate'}, inplace=True)

        return net, taz_result

    def road_type_weighted_multiple(self, target, matrix_dis, veh_types=None, weights=None, taz_id='nuts_id', taz_mx_id='id',
                                    index='index_nat', sub_len='length', net_type='type', urban=None, distance=55, cell_size=500,
                                    fac_cell=1.5):

        if weights is None:
            weights = pd.DataFrame({net_type: [0, 1, 2, 3], 'weight': [0, 5, 1.5, 0.5]})
        if veh_types is None:
            veh_types = ['lcv', 'mft', 'hft']

        net = self.net.copy()
        # target per taz
        taz_2 = self.taz[[taz_id, index, 'area', sub_len]]
        for veh_type in veh_types:
            taz_2[veh_type] = taz_2[index] / taz_2[index].sum() * target[veh_type]
            net['{}_short'.format(veh_type)] = 0
            taz_2['{}_sub'.format(veh_type)] = 0

        # weighted length per net type and taz
        net = self.net.merge(weights, on=net_type, how='left')
        net.loc[net['weight'].isna(), 'weight'] = 0
        net['weighted_length'] = net['length'] * net['weight']
        taz_2['weighted_sub'] = taz_2[sub_len]*weights.loc[weights[net_type]==0, 'weight']

        # todo: adjust urban weights
        if urban:
            net.loc[net['urban']==1, 'weighted_length'] = net.loc[net['urban']==1, 'weighted_length']*0.5

        # iterate over taz
        for t in taz_2[taz_id]:
            # adjust link weight for roads within taz
            net['weighted_length2'] = net['weighted_length']
            net['weight2'] = net['weight']
            taz_2['weighted_sub2'] = taz_2['weighted_sub']
            taz_2['weight2'] = weights.loc[weights[net_type]==0, 'weight']
            net.loc[net[taz_id] == t, 'weighted_length2'] = net.loc[net[taz_id] == t, 'weighted_length2'] * fac_cell
            net.loc[net[taz_id] == t, 'weight2'] = net.loc[net[taz_id] == t, 'weight2'] *fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] = taz_2.loc[taz_2[taz_id] == t, 'weighted_sub2'] * fac_cell
            taz_2.loc[taz_2[taz_id] == t, 'weight2'] = taz_2.loc[taz_2[taz_id] == t, 'weight2'] * fac_cell

            # list of taz within distance
            t_id = taz_2.loc[taz_2[taz_id] == t, taz_mx_id]
            matrix_t = matrix_dis[t_id, :]
            sur_ids = np.where(matrix_t <= distance)[0].tolist()

            # for small taz: add surrounding TAZ within higher distances if no sur found for distance
            if (len(sur_ids)==0) & (taz_2.loc[taz_2[taz_id]==t, 'area']<cell_size):
                sur_ids = np.where(matrix_t == matrix_t.min())[0].tolist()

            # add id of t to sur ids (all relevant cells) and get taz_ids
            sur_ids.append(t_id)
            sur_ids_taz = taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), taz_id].tolist()

            # aggregate weighted network length in sur_ids
            agg_w_len = net.loc[net[taz_id].isin(sur_ids_taz), 'weighted_length2'].sum()
            agg_w_len += taz_2.loc[taz_2[taz_mx_id].isin(sur_ids), 'weighted_sub2'].sum()

            # calculate trips per km for veh_types
            for veh_type in veh_types:
                taz_vkt = taz_2.loc[taz_2[taz_id] == t, veh_type]
                taz_veh = taz_vkt / agg_w_len
                net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] = net.loc[net[taz_id].isin(sur_ids_taz), '{}_short'.format(veh_type)] + net.loc[net[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh
                taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), '{}_sub'.format(veh_type)] = taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), '{}_sub'.format(veh_type)] + taz_2.loc[taz_2[taz_id].isin(sur_ids_taz), 'weight2'] * taz_veh

        net.drop(columns=['weight', 'weight2', 'weighted_length', 'weighted_length2'], inplace=True)
        taz_2.drop(columns=['weight2', 'weighted_sub', 'weighted_sub2'], inplace=True)

        return net, taz_2







