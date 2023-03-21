# =========================================================
# Matrices.py
# @author Nina Thomsen
# @date 11.01.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief generate travel time and distance matrices between points
# =========================================================

import time
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm


class Matrix:
    """
    Create cost matrices from OSRM
    """

    def __init__(self, conn, zone_col='zone', id_col='con_id', weight_col='weight', x_='x', y_='y', conn_geo='geometry'):
        """

        :param conn: DF / GDF of connector nodes; if DF, include x and y coordinates
        :param zone_col: name of column with zone ID in conn
        :param id_col: name of connector ID column in conn
        :param weight_col: name of weight column in conn
        :param x_: name of x coordinate column in conn (if conn does not have a geometry column)
        :param y_: name of y coordinate column in conn (if conn does not have a geometry column)
        :param conn_geo: name of geometry column in conn (if conn is a GeoDataFrame)
        :type conn: pd.DataFrame or gpd.GeoDataFrame
        :type zone_col: str
        :type id_col: str
        :type weight_col: str
        :type x_: str
        :type y_: str
        :type conn_geo: str
        """
        self.conn = conn
        if type(self.conn) == gpd.GeoDataFrame:
            self.conn.to_crs(epsg=4326, inplace=True)
            # transform to data frame with x and y coordinates
            self.conn[x_] = self.conn[conn_geo].x
            self.conn[y_] = self.conn[conn_geo].y
            self.conn = self.conn[[id_col, zone_col, weight_col, x_, y_]]

        self.zone_col = zone_col
        self.id_col = id_col
        self.weight_col = weight_col
        self.x_ = x_
        self.y_ = y_

        self.all_mtx = np.zeros((len(conn),len(conn),2))

    def osrm_request_nodes(self, save_np=None):
        """
        OSRM Request: extract travel durations and distances from OSM routing
        
        --> result: self.all_mtx as updated npy array with distances and durations between connector nodes

        :param save_np: default None; filename for saving result matrux, if given, a npy-file of self.all_mtx is saved, updated with the results per iteration
        :type save_np: str
        :return: None
        """
        # create coordinates
        coord = np.c_[self.conn[self.x_], self.conn[self.y_]]
        print("{} connector point coordinates".format(len(coord)))
        str_coords = []
        for i in range(len(coord)):
            coo = np.round(coord[i, :], 6)
            str_ = ','.join(coo.astype(str))
            str_coords.append(str_)

        # OSRM request
        url_str = 'http://router.project-osrm.org/table/v1/driving/{}?sources={}&destinations={}&annotations=duration,distance'

        # base mtx
        self.all_mtx = np.zeros((len(str_coords), len(str_coords), 2))  # from, to, dur/dis

        # 100x100 per request; first 100 rows through rest
        it = int(np.ceil(len(str_coords) / 100))
        st = 0
        time1 = time.time()
        print('start requests:', datetime.now())

        for i in tqdm(range(it)):
            st2 = 0
            en = min(st + 100, len(str_coords))
            starts = ';'.join(np.arange(0, en - st).astype(str))
            s_coords = ';'.join(str_coords[st:en])
            for ii in range(it):
                if i == ii:
                    en2 = en
                    dests = starts
                    coords_url = s_coords
                else:
                    lenstarts = en - st
                    en2 = min(st2 + 100, len(str_coords))
                    dests = ';'.join(np.arange(lenstarts, (en2 - st2) + lenstarts).astype(str))
                    d_coords = ';'.join(str_coords[st2:en2])
                    coords_url = '{};{}'.format(s_coords, d_coords)
                # request
                url_sliced = url_str.format(coords_url, starts, dests)
                r_sl = requests.get(url_sliced)
                # fill matrix
                self.all_mtx[st:en, st2:en2, 0] += np.array(r_sl.json().get('durations')).astype("float64")  # duration
                self.all_mtx[st:en, st2:en2, 1] += np.array(r_sl.json().get('distances')).astype("float64")  # distance
                st2 += 100
            st += 100
            if save_np is not None:
                np.save(save_np, self.all_mtx)
        print('total time requests:', time.time() - time1, 's',
              '({} min)'.format(str(round(float(time.time() - time1) / 60, 2))))
        print(datetime.now())

    def transform_to_taz(self):
        """
        Use connectors weights to generate the mean travel time and distances per taz, based on matrices between connectors
        :return: array with weighted mean travel times and distances and DataFrame with zone IDs and index numbers (index in matrix)
        :rtype: np.array, pd.DataFrame
        """
        # assign conn ID to zones
        self.conn[self.id_col] = np.arange(len(self.conn))
        zones = self.conn[self.zone_col].unique()
        zones_df = pd.DataFrame({"zone": zones}).sort_values("zone")
        zones_df["id"] = np.arange(len(zones_df))

        # weight matrix
        w = np.zeros((len(self.conn), len(self.conn)))

        for a in tqdm(range(len(self.conn))):
            for b in range(len(self.conn)):
                # calc od weight
                w[a, b] = float(self.conn.loc[self.conn[self.id_col] == a, self.weight_col]) * float(self.conn.loc[self.conn[self.id_col] == b, self.weight_col])

        # container for final values
        mx_z = np.zeros((len(zones), len(zones), 2))

        # fill matrix container with weighted mean for zone relations
        time1 = time.time()
        print('start aggregating zone matrix:', datetime.now())
        for o in tqdm(zones):
            id_o = int(zones_df.loc[zones_df['zone'] == o, 'id'])
            id_conn_o = list(self.conn.loc[self.conn[self.zone_col] == o, self.id_col])
            mtx_o = self.all_mtx[id_conn_o, :]
            w_o = w[id_conn_o, :]
            for d in zones:
                id_d = int(zones_df.loc[zones_df['zone'] == d, 'id'])
                id_conn_d = list(self.conn.loc[self.conn[self.zone_col] == d, self.id_col])
                mtx_od = mtx_o[:, id_conn_d]
                w_od = w_o[:, id_conn_d]
                # weighted mean
                mtx_od[:, :, 0] *= w_od
                mtx_od[:, :, 1] *= w_od
                mx_z[id_o, id_d, :] = np.nansum(mtx_od, axis=0).sum(axis=0)  # should be shape (2,)
        print('total time zones:', time.time() - time1, 's',
              '({} min)'.format(str(round(float(time.time() - time1) / 60, 2))))
        print(datetime.now())

        return mx_z, zones_df
