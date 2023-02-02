# =========================================================
# AttractionFactors.py
# @author Nina Thomsen
# @date 11.01.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief generate attraction factors per TAZ and country attributes
# =========================================================

import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
from datetime import datetime
from tqdm import tqdm


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
        self.industry_ = gpd.GeoDataFrame()

    def population_from_point(self, pop_nodes, pop_val="VALUE"):
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

    def get_industry(self, taz_cn="cntr_code", inplace=None):
        """
        Extract industrial sites (polygons) per TAZ from OSM with landuse=industrial

        @param taz_cn: str, column name with country in self.taz
        @param inplace: return GeoDataFrame; if True, GeoDataFrame is not returned
        @return: GeoDataFrame with industry area polygons
        """
        # get industrial sites from OSM using country TAZ as polygon
        countries = self.taz[taz_cn].unique()
        print(". . . Start extracting OSM industrial sites {}".format(datetime.now()))
        self.industry_ = gpd.GeoDataFrame(columns=['id', 'landuse', 'geometry'])
        id_st = 0
        for c in tqdm(countries):
            taz_c = self.taz[self.taz[taz_cn]==c]
            industry_c = ox.geometries.geometries_from_polygon(taz_c[self.taz_geo].unary_union, tags={"landuse": "industrial"})
            if len(industry_c.geom_type.unique()) > 1:
                industry_c = industry_c[industry_c.geom_type.isin(["MultiPolygon", "Polygon"])].copy()
            industry_c['id'] = list(range(len(industry_c)))
            industry_c['id'] = industry_c['id'] + id_st
            industry_c = industry_c[['id', 'landuse', 'geometry']]
            id_st += len(industry_c) + 1

            self.industry_ = pd.concat([self.industry_, industry_c])
            ind_crs = industry_c.crs
            del industry_c

        self.industry_.set_geometry('geometry', crs=ind_crs, inplace=True)

        if not inplace:
            return self.industry_

        print(". . . Finished extracting OSM industrial sites {}".format(datetime.now()))

    def industry_attributes(self, industry_gdf=None):
        """
        Aggregate total area and count of industrial sites per TAZ based on GDF with industrial areas

        @param industry_gdf: GeoDataFrame with industry areas (Polygons); if None, self.industry_ is used
        @return: self.taz as GeoDataFrame is updated with industrial area data per cell
        """
        if industry_gdf is None:
            industry = self.industry_.reset_index()
        elif type(industry_gdf) == gpd.GeoDataFrame:
            industry = industry_gdf
            industry['id'] = list(range(len(industry)))
        else:
            print('Wrong input type for industry_gdf {}'.format(str(type(industry_gdf))))
            industry = None

        # get area in km2 per polygon
        industry.to_crs(epsg=3035, inplace=True)
        industry['area_ind'] = industry.area / 1000**2
        industry.to_crs(epsg=4326, inplace=True)
        industry = industry[['id', 'area_ind', 'geometry']]
        industry.rename(columns={'id': 'id_ind'}, inplace=True)
        # aggregate industrial area per TAZ, number of industrial areas per TAZ
        taz_industry = gpd.overlay(industry, self.taz)
        taz_ind_area = taz_industry.groupby(self.taz_id)['area_ind'].sum().reset_index()
        taz_ind_area.rename(columns={'area_ind': 'ind_area_sum'}, inplace=True)
        taz_ind_count = taz_industry.groupby(self.taz_id)['id_ind'].aggregate('count').reset_index()
        taz_ind_count.rename(columns={'id_ind': 'ind_area_count'}, inplace=True)
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
            taz_scope = self.taz.copy()
            index_type = "index_int"
        else:
            taz_scope = self.taz[self.taz[taz_cn] == scope].copy()
            index_type = "index_nat"

        # columns of final geo data frame
        cols_taz = list(self.taz.columns)
        cols_taz.append(index_type)

        # normalize pop, industrial area to all taz
        taz_scope['pop_n'] = taz_scope['population'] / np.nanmean(taz_scope['population'])
        taz_scope['ind_a_n'] = taz_scope['ind_area_sum'] / np.nanmean(taz_scope['ind_area_sum'])
        taz_scope['ind_c_n'] = taz_scope['ind_area_count'] / np.nanmean(taz_scope['ind_area_count'])

        taz_scope[index_type] = alpha * (taz_scope['pop_n']**0.5 * taz_scope['ind_a_n']**0.25 * taz_scope['ind_c_n']**0.25)

        return taz_scope[cols_taz]


class BorderCrossingAtts:
    """
    Calculate country attributes like number of border crossing streets, neighboring countries that define the character of border-crossing road traffic
    """

    def __init__(self, taz, taz_cn="cntr_code", taz_geo="geom"):
        """

        @param taz: GeoDataFrame with TAZ
        @param taz_cn: str, column name of country in TAZ
        @param taz_geo: str, column name of geometry in TAZ
        """
        self.taz = taz
        self.countries = taz[taz_cn].unique()
        self.taz_cn = taz_cn
        self.taz_geo = taz_geo

        # get countries and init borders
        self.country_layer = taz.dissolve(by=taz_cn, aggfunc="sum").reset_index()
        self.country_layer = self.country_layer[[taz_cn, taz_geo]]

        self.border_layer = gpd.GeoDataFrame()

    def get_borderbuffer(self, buffer=5000):
        """
        Create GeoDataFrame with borders, using a defined buffer around these borders (i.e. with a 5000m buffer, there will be Polygons along borders with a width of 5000m)

        @param buffer: float, total buffer width in m, default 5000m
        @return: GeoDataFrame with buffer polygons around borders
        """
        border = self.country_layer.copy()
        # set buffer around polygon borders (width of border polygons)
        border.to_crs(epsg=3035, inplace=True)
        border[self.taz_geo] = border[self.taz_geo].buffer(buffer/2)
        # find borders for each country
        gdf_borderbuffer = gpd.GeoDataFrame()

        for country in self.countries:

            country_ = border[border[self.taz_cn] == country]
            country_ = country_[self.taz_geo].buffer(1)

            for index, row in border.iterrows():

                if row[self.taz_cn] != country:
                    intersec = country_.intersection(row[self.taz_geo].buffer(1))
                    if not intersec.values.is_empty[0]:
                        gdf_borderbuffer = gdf_borderbuffer.append(
                            {'country1': country, 'country2': row[self.taz_cn], 'geometry': intersec.values[0]},
                            ignore_index=True)

        # merge borders create and set geometry
        borders = pd.merge(gdf_borderbuffer, gdf_borderbuffer, left_on=['country1', 'country2'],
                           right_on=['country2', 'country1'])
        borders['geometry'] = [r['geometry_x'].union(r['geometry_y']) for i, r in borders.iterrows()]
        borders['border'] = borders['country1_x'] + borders['country1_y']
        borders = borders.set_geometry('geometry')
        borders.crs = 3035
        # remove duplicates
        borders['abc'] = ["".join(sorted(row['border'])) for i, row in borders.iterrows()]
        borders.drop_duplicates(subset="abc", inplace=True)
        # finalize
        borders = borders[['border', 'country1_x', 'country1_y', 'geometry']]
        borders.rename(columns={'country1_x': 'country1', 'country1_y': 'country2'}, inplace=True)
        borders.to_crs(epsg=4326, inplace=True)
        return borders

    def shared_borders(self, inplace=None):
        """
        Determine the share of land borders of total country border. Islands would have a share of 0, while countries without a cost would land at 1.

        @param inplace: return DataFrame; if True, DataFrame is not returned
        @return: pd.DataFrame with border length and share per country, merged to self.country_layer
        """
        self.border_layer = self.get_borderbuffer(buffer=5000)
        border_shares = pd.DataFrame()
        crs_taz = self.country_layer.crs

        for country in self.countries:
            ctr_poly = self.country_layer[self.country_layer[self.taz_cn] == country]
            ctr_border_buffer = self.border_layer[(self.border_layer['country1'] == country) | (self.border_layer['country2'] == country)]

            ctr_border_line = gpd.GeoDataFrame([1], geometry=[ctr_poly[self.taz_geo].values[0].boundary])
            ctr_border_line.crs = crs_taz

            ctr_border_clipped = gpd.clip(ctr_border_line, ctr_border_buffer)
            ctr_border_clipped.crs = crs_taz

            # calculate length of total and shared border
            ctr_border_line.to_crs(epsg=3035, inplace=True)
            ctr_border_clipped.to_crs(epsg=3035, inplace=True)
            border_length = ctr_border_line.length.values[0] / 1000
            if len(ctr_border_clipped) > 0:
                border_length_shared = ctr_border_clipped.length.values[0] / 1000
            else:
                border_length_shared = 0
            ctr_border_share = (border_length_shared / border_length)

            border_shares = border_shares.append(
                {'country': country, 'border_share': ctr_border_share, 'border_length': border_length,
                 'border_length_shared': border_length_shared}, ignore_index=True)

        # merge to countries
        self.country_layer = pd.merge(self.country_layer, border_shares, how='left', left_on=self.taz_cn, right_on='country')
        self.country_layer.drop(['country'], axis=1, inplace=True)

        if not inplace:
            return border_shares

    def border_streets(self, net, net_type="type", type_filter=None, inplace=None):
        """
        Count the number of border crossing streets per country

        @param net: GeoDataFrame with international road network
        @param net_type: str, name of column specifiying road type in net
        @param type_filter: list of road types to filter for border crossing roads; default None, meaning types [1,2]
        @param inplace: return DataFrame; if True, DataFrame is not returned
        @return: pd.DataFrame with number of border crossing streets per country
        """

        if type_filter is None:
            type_filter = [1, 2]
        # filter road types to count for border crossings
        net = net[net[net_type].isin(type_filter)]
        # get streets within border layer
        border_2m = self.get_borderbuffer(buffer=2)
        borderbuffer_streets = gpd.overlay(net, border_2m, how="union")
        # aggregate streets per country and border
        country1_grp = borderbuffer_streets.groupby('country1')['border'].aggregate('count').reset_index()
        country2_grp = borderbuffer_streets.groupby('country2')['border'].aggregate('count').reset_index()
        country1_grp = country1_grp.rename(columns={'country1': 'country'})
        country2_grp = country2_grp.rename(columns={'country2': 'country'})
        country_grp = pd.concat([country1_grp, country2_grp])
        country_grp = country_grp.groupby('country')['border'].sum().reset_index()
        country_grp.rename(columns={'border': 'border_crossings_count'}, inplace=True)

        # merge to countries
        self.country_layer = pd.merge(self.country_layer, country_grp, how='left', left_on=self.taz_cn, right_on='country')
        self.country_layer.drop(['country'], axis=1, inplace=True)

        if not inplace:
            return country_grp

    def count_neighbors(self, inplace=None):
        """
        Determine the number of direct neighbor countries (shared border)

        @param inplace: return DataFrame; if True, DataFrame is not returned
        @return: pd.DataFrame with number of neighbor countries per country
        """
        dict_ = {}
        for c in self.countries:
            tab = self.border_layer[(self.border_layer['country1'] == c) | (self.border_layer['country2'] == c)]
            adj = len(tab['border'].unique())
            dict_.update({c: {'country': c, 'neighbors': adj}})
        neighbors = pd.DataFrame.from_dict(dict_, orient='index')
        # merge to countries
        self.country_layer = pd.merge(self.country_layer, neighbors, how='left', left_on=self.taz_cn, right_on='country')
        self.country_layer.drop(['country'], axis=1, inplace=True)
        if not inplace:
            return neighbors

    def pop_area(self, pop_values=None, pop_cn="country", taz_pop="population", inplace=None):
        """
        Determine total population and area im km² per country. Population can be determined by aggregating population per taz or with and external input DataFrame

        @param pop_values: pd.DataFrame with population per country; if None, population will be aggregated from taz attributes, default None
        @param pop_cn: str, name of column with country name in pop_values
        @param taz_pop: str, name of column with population in self.taz
        @param inplace: return GeoDataFrame with border shares; if True, DataFrame is not returned
        @return: GeoDataFrame with countries and their population and area
        """
        if pop_values is None:
            # aggregate population per country from taz
            pop_agg = self.taz.groupby(self.taz_cn)[taz_pop].sum().reset_index()
            self.country_layer = pd.merge(self.country_layer, pop_agg, how='left', on=self.taz_cn)
        else:
            self.country_layer = pd.merge(self.country_layer, pop_values, how='left', left_on=self.taz_cn, right_on=pop_cn)

        # calculate area in km2
        self.country_layer.to_crs(epsg=3035, inplace=True)
        self.country_layer['area'] = self.country_layer[self.taz_geo].area / 1000**2
        self.country_layer.to_crs(epsg=4326, inplace=True)

        if not inplace:
            return self.country_layer
