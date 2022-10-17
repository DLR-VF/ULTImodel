# =========================================================
# CreateNetwork.py
# @author Nina Thomsen
# @date 29.09.2022
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Creates nodes between edges
# =========================================================

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import numpy as np
import osmnx as ox
from shapely.geometry import MultiPolygon
from shapely.ops import unary_union


class Edges:
    """Create edges from OSM"""

    def __init__(self, country, taz, taz_cn="country", taz_geo="geometry"):
        """

        :@param country: str; code or name of country
        :@param taz: GeoDataFrame including the TAZ
        :@param taz_geo: str; Name of geometry column in TAZ layer, default "geometry
        :@param taz_cn: str; Name of column in TAZ layer that defines the country of the TAZ
        """
        self.country = country
        self.edges = gpd.GeoDataFrame()
        self.taz = taz
        self.taz_cn = taz_cn
        self.taz_geo = taz_geo

    def get_edges(self, filter_highway):
        """
        Get network for selected highway types from OSM for specified country

        :param filter_highway: List of OSM highway types to include; only main types (e.g. motorway), respective
                               link-types (e.g. motorway_link) are included automatically
        :type filter_highway: list
        """
        filter_highway.extend(["{}_link".format(h) for h in filter_highway])
        fil_str = "|".join(filter_highway)

        filter_nw = '["area"!~"yes"]["highway"~"{}"]["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]["service"!~"parking|parking_aisle|driveway|private|emergency_access"]'.format(
            fil_str)

        # get polygon for country
        taz_cn = self.taz[self.taz[self.taz_cn] == self.country]
        if len(taz_cn) > 1:
            poly = MultiPolygon(unary_union(taz_cn.geometry))
        else:
            poly = taz_cn.iloc[0][self.taz_geo]
        G = ox.graph_from_polygon(poly, simplify=False, custom_filter=filter_nw)
        G = ox.simplify_graph(G)
        # option: separate secondary from here
        G = ox.simplification.consolidate_intersections(G, tolerance=0.002)
        G = ox.speed.add_edge_speeds(G, fallback=50, precision=0)
        self.edges = gpd.GeoDataFrame([x[2] for x in G.edges.data()])[["highway", "speed_kph", "geometry"]]
        # correct speeds
        self.edges.loc[self.edges["speed_kph"] > 150, "speed_kph"] = 50
        self.edges.crs = "epsg:4326"

    def set_attributes(self, taz_id, start_id=0):
        """
        Overlay self.edges with TAZ layer to split links at borders and assign TAZ-ID
        Set length, travel time, highway type and ID of self.edges GeoDataFrame

        :param taz_id: str; Name of ID-column in TAZ layer
        :param start_id: int; start of IDs for each edge, default 0
        :return:
        """

        self.edges = gpd.overlay(self.edges, self.taz[[taz_id, self.taz_geo]], how='union')
        self.edges = self.edges.explode()
        # calculate length in km
        self.edges = self.edges.to_crs(epsg=3035)
        self.edges["length"] = self.edges.length/1000
        self.edges = self.edges.to_crs(epsg=4326)
        # calculate travel time in s
        self.edges["tt"] = np.round((self.edges["length"] / self.edges["speed_kph"])*60**2)
        # set ID, type as numeric
        self.edges = self.edges.to_crs(epsg=4326)
        self.edges["ultimo_id"] = np.arange(len(self.edges)) + start_id
        self.edges["type"] = 1
        self.edges.loc[self.edges["highway"] == "trunk", "type"] = 2
        self.edges.loc[self.edges["highway"] == "primary", "type"] = 3
        self.edges.loc[self.edges["highway"] == "motorway_link", "type"] = 11
        self.edges.loc[self.edges["highway"] == "trunk_link", "type"] = 21
        self.edges.loc[self.edges["highway"] == "primary_link", "type"] = 31

    def set_nodes(self, nodes, order_col, nodeid_col="node_id", linkid_col="ultimo_id"):
        """
        Set start and end node for each edge in self.edges based on node GDF
        :param nodes: GeoDataFrame with start and end node per edge
        :param order_col: str; Name of column with order of nodes per edge in nodes
        :param nodeid_col: str; Name of column with node ID in nodes, defaul "node_id"
        :param linkid_col: str; Name of column with link ID in self.edges, default "ultimo_id"
        :return:
        """
        # merge start nodes
        st_nodes = nodes[nodes[order_col] == 0]
        self.edges = pd.merge(self.edges, st_nodes[[linkid_col, nodeid_col]], left_on="ultimo_id", right_on=linkid_col, how="left")
        self.edges.rename(columns={nodeid_col: "from_node"}, inplace=True)
        # merge start nodes
        en_nodes = nodes[nodes[order_col] == 1]
        self.edges = pd.merge(self.edges, en_nodes[[linkid_col, nodeid_col]], left_on="ultimo_id", right_on=linkid_col,
                              how="left")
        self.edges.rename(columns={nodeid_col: "to_node"}, inplace=True)
        self.edges = self.edges[["ultimo_id", "from_node", "to_node", "type", "nuts_id", "length", "speed_kph", "tt", "geometry"]]

    def subordinate_road_length(self, taz_id, taz_geo="geometry", sub_type="secondary"):
        """
        Calculate the subordinate network length per TAZ for a selected category

        :param taz_id: str; Name of ID-column in TAZ layer
        :param taz_geo: str; Name of geometry column in TAZ layer, default "geometry
        :param sub_type: str; highway type for subordinate network, default "secondary"
        :return: DataFrame with aggregated network length per TAZ
        """
        # get subordinate network from OSM
        taz_cn = self.taz[self.taz[self.taz_cn] == self.country]
        poly = MultiPolygon(unary_union(taz_cn.geometry))
        sub_edges = ox.geometries.geometries_from_polygon(poly, tags={"highway": sub_type})
        # remove geometries other than line
        if len(sub_edges.geom_type.unique()) > 1:
            sub_edges = sub_edges[sub_edges.geom_type.isin(["MultiLineString", "LineString"])].copy()
        # overlay with TAZ: assign taz id to edges
        sub_edges.crs = 'epsg:4326'
        sub_edges = gpd.overlay(sub_edges['geometry'], self.taz[[taz_id, taz_geo]], how='union')
        # calculate length per edge in km and aggregate length per taz
        sub_edges = sub_edges.to_crs(epsg=3035)
        sub_edges["length"] = sub_edges.length / 1000
        # return DataFrame with nuts_id as index
        return pd.DataFrame(sub_edges.groupby(taz_id)['length'].apply(sum))


class Nodes:
    """Create nodes at the ends of edges"""

    def __init__(self, edges):
        """
        :param edges: GeoDataFrame with LineString-Objects
        """
        self.edges = edges
        self.nodes = gpd.GeoDataFrame()
        self.nodes_unique = gpd.GeoDataFrame()

    def init_nodes(self):
        """init node GeoDataFrame"""
        self.nodes = gpd.GeoDataFrame(columns=['LinkID', 'order', 'geom'])

    def create_nodes(self, id_col, geom_col="geometry"):
        """
        Create nodes at the ends of each LineString in edges
        !! only works with LineStrings, not with MultiLineString
        !! use gpd.GeoDataFrame.explode() to transform MultiLineString to LineString

        :param id_col: str; name of column with unique identifier in edges GDF
        :param geom_col: str; name of geometry column in edges GDF, default "geometry"
        :return GeoDataFrame self.nodes
        """
        self.init_nodes()
        LinkId_list = self.edges[id_col].to_list()
        coords_list = [line.coords for i, line in enumerate(self.edges[geom_col])]
        for LinkId, coords in zip(LinkId_list, coords_list):
            index_num = self.nodes.shape[0]
            # 1st point
            self.nodes.loc[index_num, 'geom'] = Point(coords[0])
            self.nodes.loc[index_num, 'LinkID'] = LinkId
            self.nodes.loc[index_num, 'order'] = 0
            # last point
            self.nodes.loc[index_num + 1, 'geom'] = Point(coords[-1])
            self.nodes.loc[index_num + 1, 'LinkID'] = LinkId
            self.nodes.loc[index_num + 1, 'order'] = 1

    def remove_duplicates(self):
        """
        Remove duplicate nodes (i.e. at the same coordinates)
        :return: GeoDataFrame self.nodes.unique without duplicates
        """
        self.nodes['xy'] = [str(list(p.coords)) for p in self.nodes['geom']]
        # collect LinkIDs per point
        dup = pd.DataFrame(self.nodes.groupby(['xy'])['LinkID'].apply(list).reset_index())
        # transform to point data
        dup['_xy'] = [x[2:-2] for x in dup['xy']]
        dup[["x", "y"]] = [x.split(", ") for x in dup['_xy']]
        dup['geom'] = [Point(float(x), float(y)) for x, y in zip(dup['x'], dup['y'])]
        self.nodes_unique = gpd.GeoDataFrame(dup[["LinkID", "xy", "geom"]], geometry="geom")

    def set_node_id(self, id_col="node_id", start=0):
        """
        Set ID per unique node and assign this ID to all nodes
        :param id_col: str; name of ID column for nodes, default "node_id"
        :param start: int; start of IDs for each node, default 0
        :return: GeoDataFrame of self.nodes, self.nodes_unique with ID column
        """
        self.nodes_unique[id_col] = np.arange(len(self.nodes_unique)) + start
        # join unique ID to all nodes
        self.nodes = pd.merge(self.nodes, self.nodes_unique[["xy", id_col]], how="left", on="xy")
        self.nodes = self.nodes[[id_col, "LinkID", "order", "geom"]]
