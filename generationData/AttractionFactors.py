# =========================================================
# AttractionFactors.py
# @author Nina Thomsen
# @date 11.01.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief generate attraction factors per TAZ
# =========================================================

import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
from datetime import datetime


class AttractionIndex:
    """
    Calculate population and industrial areas per TAZ, transform to attraction index
    """

    def __init__(self, taz, taz_geo="geom", taz_id="nuts_id"):
        """

        @param taz: GeoDataFrame with TAZ
        @param taz_geo: str; column name of geometry of taz
        @param taz_id: str; column name of id column in taz
        """
        self.taz = taz
        self.taz_geo = taz_geo
        self.taz_id = taz_id

    def population(self, pop_nodes, pop_val="VALUE"):
        """
        Aggregate population per TAZ based on point layer with population density
        @param pop_nodes: GeoDataFrame with population density (points)
        @param pop_val: str; name of column with population density in pop_nodes
        @return: self.taz as GeoDataFrame is updated with total population per cell
        """
        # overlay taz and population
        taz_pop = gpd.overlay(pop_nodes, self.taz)

        # aggregate pop per taz
        taz_pop_g = taz_pop.groupby(self.taz_id)[pop_val].sum().reset_index()
        taz_pop_g.rename(columns={pop_val: 'population'}, inplace=True)
        self.taz = self.taz.merge(taz_pop_g, how='left', on=self.taz_id)

    def industry(self):
        """
        Aggregate total area and count of industrial sites per TAZ based on OSM data with landuse=industrial
        @return: self.taz as GeoDataFrame is updated with industrial area data per cell
        """
        # get industrial sites from OSM using TAZ as polygon
        print(". . . Start extracting OSM industrial sites {}".format(datetime.now()))
        industry = ox.geometries.geometries_from_polygon(self.taz[self.taz_geo].unary_union, tags={"landuse": "industrial"})
        print(". . . Finished extracting OSM industrial sites {}".format(datetime.now()))
        # extract only Polygons / Multipolygons
        if len(industry.geom_type.unique()) > 1:
            industry = industry[industry.geom_type.isin(["MultiPolygon", "Polygon"])].copy()
        industry['id'] = list(range(len(industry)))
        industry = industry.reset_index
        # get area in m2 per polygon
        industry.to_crs(epsg=3035, inplace=True)
        industry['area'] = industry.area
        industry.to_crs(epsg=4326, inplace=True)
        industry = industry[['id', 'area']]
        # aggregate industrial area per TAZ, number of industrial areas per TAZ
        taz_industry = gpd.overlay(industry, self.taz)
        taz_ind_area = taz_industry.groupby(self.taz_id)['area'].sum().reset_index()
        taz_ind_area.rename(columns={'area': 'ind_area_sum'})
        taz_ind_count = taz_industry.groupby(self.taz_id)['id'].aggregate('count').reset_index()
        taz_ind_count.rename(columns={'id': 'ind_area_count'})
        taz_industry = pd.merge(taz_ind_area, taz_ind_count, on=self.taz_id)
        self.taz = self.taz.merge(taz_industry, how='left', on=self.taz_id)

    def attraction_index(self, scope=None, taz_cn='cntr_code', alpha=1.):
        """
        Create attraction index with Cobb-Douglas transformation
        @param scope: None or str; either look at all TAZ or single country (str of ISO-code)
        @param taz_cn: str; column name for country identifier in self.taz
        @param alpha: float; alpha parameter for Cobb Douglas formula
        @return: GeoDataFrame with taz in scope and attraction index for scope
        """
        # set scope
        if scope is None:
            taz_scope = self.taz.copy
            index_type = "index_int"
        else:
            taz_scope = self.taz[self.taz[taz_cn] == scope]
            index_type = "index_nat"

        # columns of final geo data frame
        cols_taz = self.taz.columns
        cols_taz.extend(index_type)

        # normalize pop, industrial area to all taz
        taz_scope['pop_n'] = taz_scope['population'] / np.nanmean(taz_scope['population'])
        taz_scope['ind_a_n'] = taz_scope['ind_area_sum'] / np.nanmean(taz_scope['ind_area_sum'])
        taz_scope['ind_c_n'] = taz_scope['ind_area_count'] / np.nanmean(taz_scope['ind_area_count'])

        taz_scope[index_type] = alpha * (taz_scope['pop_n']**0.5 * taz_scope['ind_a_n']**0.25 * taz_scope['ind_c_n']**0.25)

        return taz_scope[cols_taz]
