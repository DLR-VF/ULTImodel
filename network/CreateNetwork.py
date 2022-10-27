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
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans


def ckdnearest(gdA, gdB, bcol, which):
    nA = np.array(list(zip(gdA.geometry.y, gdA.geometry.x)) )
    nB = np.array(list(zip(gdB.geometry.y, gdB.geometry.x)) )
    btree = cKDTree(nB)
    dist, idx = btree.query(nA,k=[which])
    dist = np.squeeze(dist)
    idx = np.squeeze(idx)
    #print (idx)
    df = pd.DataFrame.from_dict({'distance': dist,#.astype(int),
                                 bcol : gdB.loc[idx, bcol].values})
    return df


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

    def subordinate_road_length(self, taz_id, sub_type="secondary"):
        """
        Calculate the subordinate network length per TAZ for a selected category

        :param taz_id: str; Name of ID-column in TAZ layer
        :param taz_geo: str; Name of geometry column in TAZ layer, default "geometry
        :param sub_type: str; highway type for subordinate network, default "secondary"
        :return: DataFrame with aggregated network length per TAZ
        """
        # get subordinate network from OSM
        taz_cn = self.taz[self.taz[self.taz_cn] == self.country]
        if len(taz_cn) > 1:
            poly = MultiPolygon(unary_union(taz_cn.geometry))
        else:
            poly = taz_cn.iloc[0][self.taz_geo]
        sub_edges = ox.geometries.geometries_from_polygon(poly, tags={"highway": sub_type})
        # remove geometries other than line
        if len(sub_edges.geom_type.unique()) > 1:
            sub_edges = sub_edges[sub_edges.geom_type.isin(["MultiLineString", "LineString"])].copy()
        # overlay with TAZ: assign taz id to edges
        sub_edges.crs = 'epsg:4326'
        sub_edges = gpd.overlay(sub_edges[['geometry']], self.taz[[taz_id, self.taz_geo]], how='union')
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


class Connectors:
    """
    Connect the TAZ to the road network: find possible connector locations and identif corresponding nodes
    """

    def __init__(self):
        """


        """


    def find_connector_locations(self, taz, pop, taz_id="nuts_id", taz_geom="geom", epsg_pop=4326):
        """
        Identify suitable connector locations per TAZ depending on population density

        @param epsg_pop: int; CRS of population dataframe, default 4326
        @param taz: GeoDataFrame with the TAZ
        @param taz_geom: str; name of geometry column in TAZ GDF
        @param taz_id: str; name of ID column for TAZ IDs
        @param pop: GeoDataFrame of points with population density

        @return GeoDataFrame with connector locations; Point layer
        """
        cols = ["nuts_id","c_n", "weight", "geometry"]
        return_con = gpd.GeoDataFrame(columns=cols)

        taz = taz.to_crs(epsg=3035)
        taz['area'] = taz.area/1000**2
        taz["numcon"] = 3
        taz.loc[taz["area"] > 1000, "numcon"] = 4
        taz.loc[taz["area"] > 7500, "numcon"] = 5
        taz.loc[taz["area"] > 15000, "numcon"] = 6

        # overlay pop and taz
        taz = taz.to_crs(epsg=epsg_pop)
        pop_taz = gpd.overlay(pop, taz[[taz_id, taz_geom]])

        for t in taz[taz_id].unique():
            pop_t = pop_taz[pop_taz[taz_id]==t]
            if len(pop_t) > 0:
                n_c = int(taz.loc[taz[taz_id]==t, "numcon"])
                kmeans = KMeans(init="random", n_clusters=n_c, n_init=10, max_iter=300, random_state=42)
                kmeans.fit(list(zip(pop_t.geometry.x, pop_t.geometry.y)))
                # assign pop to cluster
                pop_t = pd.concat([pop_t.reset_index(), pd.Series(kmeans.labels_, name="c_n")], axis=1)
                # create GeoDataFrame with weights and coordinates of centers (connectors)
                conn = gpd.GeoDataFrame(pop_t.groupby("c_n")['VALUE'].apply(sum).reset_index(),
                                 geometry=[Point(xy) for xy in kmeans.cluster_centers_])
                conn['weight'] = conn['VALUE'] / conn['VALUE'].sum()
                conn['nuts_id'] = t
                conn = conn[cols]
                return_con = return_con.append(conn)

        return return_con

    def identify_connector_nodes(self, nodes, con, node_no="node_id", geom="geom", con_no="c_n", zone="nuts_id", weight="weight"):
        """
        Move connector locations to network nodes

        @param nodes: GeoDataFrame Points; network nodes
        @param con: GeoDataFrame Points; connector locations
        @param node_no: str; column name with node identifyer

        @return GeoDataFrame with connector nodes
        """
        # set CRS
        nodes = nodes.to_crs(epsg=3035)
        con = con.to_crs(epsg=3035)

        # find one node per connector location
        distances = ckdnearest(con, nodes, node_no, 1)
        output = con[[con_no, zone, weight]].merge(distances, left_index=True, right_index=True)

        # check for duplicates
        con_nodes = output[[node_no, zone, con_no]].groupby([zone, node_no], as_index=False).count()

        if len(con_nodes) != len(output):
            print('Duplicates detected!')
            output["distprod"] = output['distance'] * 1/output[weight]
            dup_no = con_nodes.loc[con_nodes[con_no] > 1][node_no].tolist()
            dup = output[output[node_no].isin(dup_no)].copy()
            min_dis = dup.groupby([zone, node_no])['distprod'].transform('min').unique()
            delete = dup.loc[~dup['distprod'].isin(min_dis), con_no].tolist()
            keep = dup.loc[dup['distprod'].isin(min_dis), con_no].tolist()
            print(len(delete), sum(output.loc[output[con_no].isin(delete), weight]))
            output = output[~output[con_no].isin(delete)].copy()
            # correct weights: add deleted weight to remaining connector (duplicates only removed for same TAZ)
            delete_weights = dup.groupby([zone, node_no])[weight].apply(sum).reset_index()
            delete_weights.columns = [zone, node_no, "newweight"]
            output = output.merge(delete_weights, on=[zone, node_no], how="left")
            output.loc[~output['newweight'].isna(), weight] = output['newweight']
            print(sum(output[weight]))
        else:
            print('No duplicates detected!')

        # create GDF with correct connector nodes
        connodes = nodes[[node_no, geom]].merge(output[[node_no, zone, con_no, weight]], on=node_no, how="right")

        return connodes
