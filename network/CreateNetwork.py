# =========================================================
# CreateNetwork.py
# @author Nina Thomsen
# @date 21.12.2022
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Creates nodes between edges
# =========================================================
from datetime import datetime

import geopandas as gpd
import pandas as pd
from sklearn.cluster import KMeans
from shapely.geometry import Point, LineString, MultiPolygon
import numpy as np
import osmnx as ox
from shapely.ops import unary_union
from scipy.spatial import cKDTree


def ckdnearest(gdA, gdB, bcol, n):
    nA = np.array(list(zip(gdA.geometry.y, gdA.geometry.x)))
    nB = np.array(list(zip(gdB.geometry.y, gdB.geometry.x)))
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=n)
    df = pd.DataFrame(columns=[bcol, 'distance', '{}_near'.format(bcol)])
    if n > 1:
        for i in range(idx.shape[1]):
            dfi = pd.DataFrame.from_dict({bcol: gdA[bcol], 'distance': dist[:, i],  # .astype(int),
                                          '{}_near'.format(bcol): gdB.reset_index().loc[idx[:, i], bcol].values})
            df = pd.concat([df, dfi], axis=0)
    else:
        dfi = pd.DataFrame.from_dict({bcol: gdA[bcol], 'distance': dist,  # .astype(int),
                                      '{}_near'.format(bcol): gdB.reset_index().loc[idx, bcol].values})
        df = pd.concat([df, dfi], axis=0)
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

    def get_polygon(self, taz_id, proj_crs=3035, buffer=2500):
        """
        Define polygon(s) for network extraction; split islands and exclaves and apply buffer to include all streets near the border

        @param taz_id: str; Name of ID-column in TAZ layer
        @param proj_crs: int; EPSG number of projected crs for defining buffer in meter, default 3035 (LAEA europe)
        @param buffer: int, float; buffer in meter
        @return: GeoDataFrame with polygons for extraction; multiple polygons if there are islands or exclaves
        """
        # start: filter country
        taz_cn = self.taz[self.taz[self.taz_cn] == self.country].copy()
        taz_crs = taz_cn.crs
        # get all separate polygons (islands, exclaves)
        t_e = taz_cn.explode(index_parts=True)

        # find islands belonging to the same TAZ id and create hull
        if len(t_e[taz_id].unique()) > 1:
            isl = t_e.groupby(taz_id)[self.taz_geo].agg('count').reset_index()
            isl = isl[isl[self.taz_geo] > 1]
            # create hull for same TAZ id islands
            for i in isl[taz_id]:
                hull = taz_cn.loc[taz_cn[taz_id] == i].unary_union.convex_hull
                taz_cn.loc[taz_cn[taz_id] == i, self.taz_geo] = hull

        taz_cn_e = taz_cn.dissolve().explode(index_parts=True)
        # apply buffer (include all roads in coastal regions, islands that are not far away and their bridges
        taz_cn_e.to_crs(epsg=proj_crs, inplace=True)
        taz_cn_e[self.taz_geo] = taz_cn_e[self.taz_geo].buffer(buffer)
        taz_cn_e.to_crs(taz_crs, inplace=True)

        # return polygons
        return taz_cn_e

    def get_edges(self, filter_highway, polygons):  # add poly attribute, either none or input
        """
        Get network for selected highway types from OSM for specified country

        :@param filter_highway: list; List of OSM highway types to include; only main types (e.g. motorway), respective
                               link-types (e.g. motorway_link) are included automatically
        :@param polygons: GeoDataFrame; polygons of country / region to extract network from, several polygons possible
        :@return self.edges includes routable network
        """
        filter_highway.extend(["{}_link".format(h) for h in filter_highway])
        fil_str = "|".join(filter_highway)

        filter_nw = '["area"!~"yes"]["highway"~"{}"]["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]["service"!~"parking|parking_aisle|driveway|private|emergency_access"]'.format(
            fil_str)

        # get polygon for country
        polygons = polygons
        n_poly = 0  # number of polygons with road network
        for i, poly in enumerate(polygons[self.taz_geo][:]):
            try:
                G = ox.graph_from_polygon(poly, simplify=False, custom_filter=filter_nw)
                G = ox.simplify_graph(G)
                G = ox.simplification.consolidate_intersections(G, tolerance=0.002)
                G = ox.speed.add_edge_speeds(G, fallback=50, precision=0)
                edges = gpd.GeoDataFrame([x[2] for x in G.edges.data()])[["highway", "speed_kph", "geometry"]]
                if i == 0:
                    self.edges = edges.copy()
                else:
                    self.edges = pd.concat([self.edges, edges.reset_index()])
                n_poly += 1
            except:
                print('No edges for polygon {}'.format(i))
        # correct speeds
        self.edges.loc[self.edges["speed_kph"] > 150, "speed_kph"] = 50
        self.edges.crs = "epsg:4326"

        return n_poly

    def set_attributes(self, taz_id, start_id=0):
        """
        Overlay self.edges with TAZ layer to split links at borders and assign TAZ-ID
        Set length, travel time, highway type and ID of self.edges GeoDataFrame

        :@param taz_id: str; Name of ID-column in TAZ layer
        :@param start_id: int; start of IDs for each edge, default 0
        :@return: self.edges is updated with necessary attributes
        """

        self.edges = gpd.overlay(self.edges, self.taz[[taz_id, self.taz_geo]], how='union')
        self.edges = self.edges.explode()
        # calculate length in km
        self.edges = self.edges.to_crs(epsg=3035)
        self.edges["length"] = self.edges.length / 1000
        self.edges = self.edges.to_crs(epsg=4326)
        # calculate travel time in s
        self.edges["tt"] = np.round((self.edges["length"] / self.edges["speed_kph"])*(60**2))
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
        :@param nodes: GeoDataFrame with start and end node per edge
        :@param order_col: str; Name of column with order of nodes per edge in nodes
        :@param nodeid_col: str; Name of column with node ID in nodes, default "node_id"
        :@param linkid_col: str; Name of column with link ID in nodes, default "ultimo_id"
        :@return: self.edges includes from and to nodes
        """
        # merge start nodes
        st_nodes = nodes[nodes[order_col] == 0]
        self.edges = pd.merge(self.edges, st_nodes[[linkid_col, nodeid_col]], left_on="ultimo_id", right_on=linkid_col,
                              how="left")
        self.edges.rename(columns={nodeid_col: "from_node"}, inplace=True)
        # merge start nodes
        en_nodes = nodes[nodes[order_col] == 1]
        self.edges = pd.merge(self.edges, en_nodes[[linkid_col, nodeid_col]], left_on="ultimo_id", right_on=linkid_col,
                              how="left")
        self.edges.rename(columns={nodeid_col: "to_node"}, inplace=True)
        self.edges = self.edges[
            ["ultimo_id", "from_node", "to_node", "type", "nuts_id", "length", "speed_kph", "tt", "geometry"]]

    def subordinate_road_length(self, taz_id, sub_type="secondary"):
        """
        Calculate the subordinate network length per TAZ for a selected category

        :@param taz_id: str; Name of ID-column in TAZ layer
        :@param taz_geo: str; Name of geometry column in TAZ layer, default "geometry
        :@param sub_type: str; highway type for subordinate network, default "secondary"
        :@return: DataFrame with aggregated network length per TAZ
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
        :@param edges: GeoDataFrame with LineString-Objects
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

        :@param id_col: str; name of column with unique identifier in edges GDF
        :@param geom_col: str; name of geometry column in edges GDF, default "geometry"
        :@return GeoDataFrame self.nodes
        """
        self.init_nodes()
        LinkId_list = self.edges[id_col].to_list()
        coords_list = [line.coords for i, line in enumerate(self.edges[geom_col])]
        # TODO convert to numpy array with arr = construct_1d_object_array_from_listlike(values)
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
        self.nodes.set_geometry('geom', inplace=True)
        self.nodes.crs = 4326

    def remove_duplicates(self):
        """
        Remove duplicate nodes (i.e. at the same coordinates)
        :@return: GeoDataFrame self.nodes.unique without duplicates
        """
        self.nodes['xy'] = [str(list(p.coords)) for p in self.nodes['geom']]
        # collect LinkIDs per point
        dup = pd.DataFrame(self.nodes.groupby(['xy'])['LinkID'].apply(list).reset_index())
        # transform to point data
        dup['_xy'] = [x[2:-2] for x in dup['xy']]
        dup[["x", "y"]] = [x.split(", ") for x in dup['_xy']]
        dup['geom'] = [Point(float(x), float(y)) for x, y in zip(dup['x'], dup['y'])]
        self.nodes_unique = gpd.GeoDataFrame(dup[["LinkID", "xy", "geom"]], geometry="geom")
        self.nodes_unique.crs = self.nodes.crs

    def set_node_id(self, id_col="node_id", start=0):
        """
        Set ID per unique node and assign this ID to all nodes
        :@param id_col: str; name of ID column for nodes, default "node_id"
        :@param start: int; start of IDs for each node, default 0
        :@return: GeoDataFrame of self.nodes, self.nodes_unique with ID column
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
        self.connectors = gpd.GeoDataFrame()

    def find_connector_locations(self, taz, pop, taz_id="nuts_id", taz_geom="geom", epsg_pop=4326):
        """
        Identify suitable connector locations per TAZ depending on population density

        @param epsg_pop: int; CRS of population dataframe, default 4326
        @param taz: GeoDataFrame with the TAZ
        @param taz_geom: str; name of geometry column in TAZ GDF
        @param taz_id: str; name of ID column for TAZ IDs
        @param pop: GeoDataFrame of points with population density

        @return self.connectors as GeoDataFrame with connector locations; Point layer
        """
        cols = ["nuts_id", "c_n", "weight", "geometry"]
        return_con = gpd.GeoDataFrame(columns=cols)

        taz = taz.to_crs(epsg=3035)
        taz['area'] = taz.area / 1000 ** 2
        taz["numcon"] = 3
        taz.loc[taz["area"] > 1000, "numcon"] = 4
        taz.loc[taz["area"] > 7500, "numcon"] = 5
        taz.loc[taz["area"] > 15000, "numcon"] = 6

        # overlay pop and taz
        taz = taz.to_crs(epsg=epsg_pop)
        pop_taz = gpd.overlay(pop, taz[[taz_id, taz_geom]])

        for t in taz[taz_id].unique():
            pop_t = pop_taz[pop_taz[taz_id] == t]
            if len(pop_t) > 0:
                n_c = int(taz.loc[taz[taz_id] == t, "numcon"])
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

        self.connectors = return_con
        self.connectors.crs = epsg_pop
        self.connectors['c_n'] = self.connectors.reset_index().index


    def identify_connector_nodes(self, nodes, node_no="node_id", node_geom="geom", zone="nuts_id",
                                 weight="weight"):
        """
        Move connector locations to network nodes

        @return GeoDataFrame with connector nodes
        @param nodes: GeoDataFrame Points; network nodes
        @param node_no: str; column name with node identifyer
        @param node_geom: str; name of geometry column in nodes GDF
        @param zone: str; column name with zone identifyer
        @param weight: str; column name with weight of connector
        @return: GeoDataFrame with connector nodes
        """
        # set CRS
        nodes = nodes.to_crs(epsg=3035)
        con = self.connectors.copy()
        con = con.to_crs(epsg=3035)

        # find one node per connector location
        distances = ckdnearest(con, nodes, node_no, 1)
        output = con[['c_n', zone, weight]].merge(distances, left_index=True, right_index=True)

        # check for duplicates
        con_nodes = output[[node_no, zone, 'c_n']].groupby([zone, node_no], as_index=False).count()

        if len(con_nodes) != len(output):
            print('Duplicates detected!')
            output["distprod"] = output['distance'] * 1 / output[weight]
            dup_no = con_nodes.loc[con_nodes['c_n'] > 1][node_no].tolist()
            dup = output[output[node_no].isin(dup_no)].copy()
            min_dis = dup.groupby([zone, node_no])['distprod'].transform('min').unique()
            delete = dup.loc[~dup['distprod'].isin(min_dis), 'c_n'].tolist()
            print(len(delete), sum(output.loc[output['c_n'].isin(delete), weight]))
            output = output[~output['c_n'].isin(delete)].copy()
            # correct weights: add deleted weight to remaining connector (duplicates only removed for same TAZ)
            delete_weights = dup.groupby([zone, node_no])[weight].apply(sum).reset_index()
            delete_weights.columns = [zone, node_no, "newweight"]
            output = output.merge(delete_weights, on=[zone, node_no], how="left")
            output.loc[~output['newweight'].isna(), weight] = output['newweight']
            print(sum(output[weight]))
        else:
            print('No duplicates detected!')

        # create GDF with correct connector nodes
        connodes = nodes[[node_no, node_geom]].merge(output[[node_no, zone, 'c_n', weight]], on=node_no, how="right")

        return connodes


class Ferries:
    """
    Create connections over water between islands and main land, using ferry routes and bridges as reference
    """

    def __init__(self, taz, scope=None, taz_cn='cntr_code', taz_geo='geom'):
        """

        @param taz: GeoDataFrame with TAZ cells
        @param scope: str or list; region of taz to be considered; default None includes all taz, str with country code reduces taz to country, list with country codes reduces taz to countries in list
        @param taz_cn: str; name of column with country code in taz
        @param taz_geo: str; name of geometry column in taz
        """
        self.taz = taz
        self.scope = scope
        self.taz_cn = taz_cn
        self.taz_geo = taz_geo
        # outputs
        if scope is None:
            self.region = taz
        elif type(scope) is str:
            self.region = taz[taz[taz_cn]==scope]
        else:
            self.region = taz[taz[taz_cn].isin(scope)]
        self.ferry = gpd.GeoDataFrame()
        self.nodes = gpd.GeoDataFrame()

    def find_ferry_routes(self, buffer_water=0.05):
        """
        Extract ferry routes from OSM for water body around region (defined by scope in __init__)
        @param buffer_water: float; buffer to use around land mass to extract water polygon, in degrees
        @return: self.ferry and self.nodes include ferry routes and their end nodes
        """
        # create polygon of water body for region
        hull = self.region.unary_union.convex_hull
        water = hull.difference(self.region[self.taz_geo].buffer(buffer_water).unary_union)
        # get ferry routes and bridges with roads for water body
        print(". . . extracting bridges and ferries from OSM")
        self.ferry = ox.geometries.geometries_from_polygon(water, tags={"route": "ferry", "highway": ["motorway", "trunk"]})
        print(". . . bridges and ferries extracted!")
        # remove non-LineString geometries
        if len(self.ferry.geom_type.unique()) > 1:
            self.ferry = self.ferry[self.ferry.geom_type.isin(["MultiLineString", "LineString"])].copy()
        # filter only ferries with motor vehicle or bridges with highway attribute
        self.ferry = self.ferry[(~self.ferry['highway'].isna()) | (self.ferry['motor_vehicle'] == "yes")].copy()
        # set attributes and reduce to relevant columns
        self.ferry['id'] = list(range(len(self.ferry)))
        self.ferry = self.ferry.reset_index()
        self.ferry['type'] = 9  # ferry
        self.ferry.loc[~self.ferry['highway'].isna(), 'type'] = 8  # bridge
        self.ferry = self.ferry[['id', 'geometry', 'type']]
        # create ferry end nodes
        nodes = Nodes(self.ferry)
        nodes.create_nodes('id')
        nodes.remove_duplicates()
        nodes.set_node_id()
        self.nodes = nodes.nodes.set_geometry('geom')
        self.nodes.crs = 4326

    def ferry_national(self, region_id='id', ferry_buffer=0.01):
        """
        Create GeoDataFrames of ferry routes and nodes with a country; filter other relations out of self.ferry and self.nodes

        @param region_id: str; name of column with identifier
        @param ferry_buffer: float; buffer to use to look for end nodes in coast region
        @return: GeoDataFrame ferry_routes with national ferry routes and GeoDataFrame ferry_nodes with respective nodes
        """
        # find nodes in border area
        self.region.to_crs(self.nodes.crs, inplace=True)
        reg_buf = self.region.copy()
        reg_buf[self.taz_geo] = reg_buf[self.taz_geo].buffer(ferry_buffer)
        border_ferry = gpd.overlay(self.nodes, reg_buf[[region_id, self.taz_geo]])
        # find end nodes within two different cells
        n_g = border_ferry.groupby('LinkID')['order', region_id].nunique().reset_index()
        filter_ferry = list(n_g.loc[(n_g['order'] > 1) & (n_g[region_id] > 1), 'LinkID'])
        # final result
        ferry_routes = self.ferry[self.ferry['id'].isin(filter_ferry)].copy()
        ferry_nodes = self.nodes[self.nodes['LinkID'].isin(filter_ferry)].copy()

        return ferry_routes, ferry_nodes

    def ferry_international(self, ferry_buffer=0.01):
        """
        Create GeoDataFrame of ferry routes and nodes between different countries, filter national relations out of self.ferry and self.nodes

        @param ferry_buffer: float; buffer to use to look for end nodes in coast region
        @return: GeoDataFrame ferry_routes with international ferry routes and GeoDataFrame ferry_nodes with respective nodes
        """
        # find nodes in border area
        self.region.to_crs(self.nodes.crs, inplace=True)
        reg_buf = self.region.copy()
        reg_buf[self.taz_geo] = reg_buf[self.taz_geo].buffer(0.01)
        border_ferry = gpd.overlay(self.nodes, reg_buf[[self.taz_cn, self.taz_geo]])
        # find end nodes within two different countries
        n_g = border_ferry.groupby('LinkID')['order', self.taz_cn].nunique().reset_index()
        filter_ferry = list(n_g.loc[(n_g['order'] > 1) & (n_g[self.taz_cn] > 1), 'LinkID'])
        # final result
        ferry_routes = self.ferry[self.ferry['id'].isin(filter_ferry)].copy()
        ferry_nodes = self.nodes[self.nodes['LinkID'].isin(filter_ferry)].copy()

        return ferry_routes, ferry_nodes

    def ferry_to_road(self, ferry_routes, ferry_nodes, network_roads, network_nodes, ferry_buffer_m=2500, speed_ferry=30,
                      roads_id="ultimo_id", roads_fr="from_node", roads_to="to_node", node_id="node_id"):
        """
        Connect ferry end nodes to road network and add as links (direct lines)

        @param ferry_routes: GeoDataFrame with ferry routes (lines)
        @param ferry_nodes: GeoDataFrame with ferry nodes (points)
        @param network_roads: GeoDataFrame with road network (lines)
        @param network_nodes: GeoDataFrame with unique road network nodes (points)
        @param ferry_buffer_m: float; buffer size around ferry nodes for finding connecting roads in meter
        @param speed_ferry: float; speed on ferry link in kph
        @param roads_id: str; name of column with road identifier in network_roads
        @param roads_fr: str; name of column with from node in network_roads
        @param roads_to: str; name of column with to node in network_roads
        @param node_id: str; name of column with node identifier in network_nodes
        @return: GeoDataFrame of road network, connected along ferry routes
        """
        # find roads within buffer of ferry nodes
        ferry_buffer = ferry_nodes.copy()
        ferry_buffer.to_crs(epsg=3035, inplace=True)
        ferry_buffer['geom'] = ferry_buffer['geom'].buffer(ferry_buffer_m)
        ferry_buffer.to_crs(epsg=4326, inplace=True)
        border_roads = gpd.overlay(network_roads, ferry_buffer)

        # iterate over ferry routes, connect with closest road node in both directions
        lines = gpd.GeoDataFrame(columns=network_roads.columns)
        id_st = int(network_roads[roads_id].max())

        for ferry in ferry_routes['id'].unique():
            roads_ferry = border_roads[border_roads['LinkID'] == ferry]
            nodes_ferry = ferry_nodes[ferry_nodes['LinkID'] == ferry]
            if len(roads_ferry) == 0:
                print('No road connection for ferry route {}'.format(ferry))
            else:
                # start node of ferry
                roads_ferry_0 = roads_ferry[roads_ferry['order'] == 0]
                node_ferry_0 = nodes_ferry[nodes_ferry['order'] == 0]
                # end node of ferry
                roads_ferry_1 = roads_ferry[roads_ferry['order'] == 1]
                node_ferry_1 = nodes_ferry[nodes_ferry['order'] == 1]
                if (len(roads_ferry_0) > 0) & (len(roads_ferry_1) > 1):
                    # length ferry route
                    len_ferry = ferry_routes[ferry_routes['id'] == ferry].copy()
                    len_ferry.to_crs(epsg=3035, inplace=True)
                    len_ferry['len'] = len_ferry.length / 1000
                    len_ferry = float(len_ferry['len'])
                    # type (ferry=9, bridge=8)
                    type_ferry = int(ferry_routes.loc[ferry_routes['id'] == ferry, 'type'])
                    if type_ferry == 8:
                        speed = 60
                    else:
                        speed = speed_ferry

                    # connect roads at start and end
                    # start from nodes 1-0
                    nodes_0f = roads_ferry_0[[roads_id, roads_fr]].rename(columns={roads_fr: node_id})
                    # find closest node to ferry node
                    nodes_0f = nodes_0f.merge(network_nodes[[node_id, 'geom']], how='left',
                                              left_on=node_id, right_on=node_id)
                    nodes_0f = nodes_0f.set_geometry('geom')
                    node_0f = ckdnearest(node_ferry_0, nodes_0f, node_id, 1)
                    node_0f = nodes_0f[nodes_0f[node_id] == int(node_0f['{}_near'.format(node_id)])]
                    # end to nodes 1-0
                    nodes_1t = roads_ferry_1[[roads_id, roads_to]].rename(columns={roads_to: node_id})
                    # find closest node to ferry node
                    nodes_1t = nodes_1t.merge(network_nodes[[node_id, 'geom']], how='left',
                                              left_on=node_id, right_on=node_id)
                    nodes_1t = nodes_1t.set_geometry('geom')
                    node_1t = ckdnearest(node_ferry_1, nodes_1t, node_id, 1)
                    node_1t = nodes_1t[nodes_1t[node_id] == int(node_1t['{}_near'.format(node_id)])]
                    node_pair_10 = [node_1t.iloc[0], node_0f.iloc[0]]
                    # start to nodes 0-1
                    nodes_0t = roads_ferry_0[[roads_id, roads_to]].rename(columns={roads_to: node_id})
                    # find closest node to ferry node
                    nodes_0t = nodes_0t.merge(network_nodes[[node_id, 'geom']], how='left',
                                              left_on=node_id, right_on=node_id)
                    nodes_0t = nodes_0t.set_geometry('geom')
                    node_0t = ckdnearest(node_ferry_0, nodes_0t, node_id, 1)
                    node_0t = nodes_0t[nodes_0t[node_id] == int(node_0t['{}_near'.format(node_id)])]
                    # end from nodes 0-1
                    nodes_1f = roads_ferry_1[[roads_id, roads_fr]].rename(columns={roads_fr: node_id})
                    # find clostest node to ferry node
                    nodes_1f = nodes_1f.merge(network_nodes[[node_id, 'geom']], how='left',
                                              left_on=node_id, right_on=node_id)
                    nodes_1f = nodes_1f.set_geometry('geom')
                    node_1f = ckdnearest(node_ferry_1, nodes_1f, node_id, 1)
                    node_1f = nodes_1f[nodes_1f[node_id] == int(node_1f['{}_near'.format(node_id)])]
                    node_pair_01 = [node_0t.iloc[0], node_1f.iloc[0]]
                    pairs = [node_pair_01, node_pair_10]
                    # create line
                    for l in pairs:
                        # get coordinates and create LineString
                        node_from = l[0]['geom']
                        node_to = l[1]['geom']
                        line = LineString([[node_from.x, node_from.y], [node_to.x, node_to.y]])
                        # add to Lines GDF
                        lines.loc[len(lines) + 1] = [id_st, l[0][node_id], l[1][node_id], type_ferry, 'FERRY',
                                                     len_ferry, speed, (len_ferry / speed) * 60 ** 2, line]
                        id_st += 1
                else:
                    print('no start / end connection for ferry route {}'.format(ferry))

        network_roads = pd.concat([network_roads, lines])

        return network_roads


class CombineNetworks:
    """
    Merge road networks of multiple countries / regions to create a routable international network
    """

    def __init__(self, taz, roads, taz_cn='cntr_code'):
        """

        @param taz: GDF; TAZ for regions
        @param roads: GDF; road network covering all regions
        @param taz_cn: str; column name for country in taz
        """
        self.taz = taz
        self.taz_cn = taz_cn
        self.network = roads
        self.network_int = gpd.GeoDataFrame()
        self.border_roads = gpd.GeoDataFrame()
        self.borders = gpd.GeoDataFrame()
        self.countries = taz[taz_cn].unique()
        self.dict_borders = {}

    def find_borders(self, buffer=2500, taz_geo='geometry'):
        """
        Create polygons with a defined buffer around the borders

        @param buffer: int, float; buffer width in m
        @param taz_geo: str; name of geometry column in TAZ GDF
        @return: GDF with border polygons
        """
        # dissolve zones per country
        taz_dis = self.taz.dissolve(by=self.taz_cn, aggfunc="sum").reset_index()
        # set buffer around polygon borders (width of border polygons)
        taz_dis.to_crs(epsg=3035, inplace=True)
        taz_dis[taz_geo] = taz_dis[taz_geo].buffer(buffer / 2)
        # find borders for each country
        gdf_borderbuffer = gpd.GeoDataFrame()

        for country in self.countries:

            country_ = taz_dis[taz_dis[self.taz_cn] == country]
            country_ = country_[taz_geo].buffer(1)

            for index, row in taz_dis.iterrows():

                if row[self.taz_cn] != country:
                    intersec = country_.intersection(row[taz_geo].buffer(1))
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
        borders.to_crs(epsg=4326, inplace=True)

        self.borders = borders

    def get_border_roads(self, roads_cn="cn"):
        """
        Find all roads within the border buffer

        @param roads_cn: str; column name with country attribute in roads GDF
        @return: self.border_roads; GDF with all roads within border buffer
        """
        if self.borders.crs == self.network.crs:
            self.border_roads = gpd.overlay(self.network, self.borders, how='union')
            # remove all streets that are not within the border
            self.border_roads = self.border_roads[~self.border_roads['border'].isna()].copy()
            # remove streets that are assigned to a wrong border
            self.border_roads = self.border_roads[(self.border_roads[roads_cn] == self.border_roads['country1_x']) |
                                                  (self.border_roads[roads_cn] == self.border_roads['country1_y'])]
        else:
            print("CRS not matching!\n    Borders {}\n    Roads {}".format(self.borders.crs, self.border_roads.crs))

    def connect_border_roads(self, nodes_int, id_st=900000000, roads_cn='cn', roads_type='type', filter_types=None,
                             roads_id="ultimo_id", roads_fr="from_node", roads_to="to_node", node_id="node_id",
                             node_geo="geom"):
        """
        Create new directed lines between the nearest end nodes of each country at the border

        @param nodes_int: GDF; contains nodes of network
        @param id_st: int; start ID for new lines
        @param roads_cn: str; column name for country in self.network
        @param roads_type: str; column name for road type in self.network
        @param filter_types: list; types of roads to be included for international connections
        @param roads_id: str; column name for id in self.network
        @param roads_fr: str; column name for the start node in self.network
        @param roads_to: str; column name for the destination node in self.network
        @param node_id: str; column name for id in nodes
        @param node_geo: str; column name for geometry in nodes
        @return: GDF with connected network, dictionary with number of connected roads
        """
        if filter_types is None:
            filter_types = [1, 2]
        lines = gpd.GeoDataFrame(columns=self.network.columns)

        id_st = id_st
        # iterate over borders
        print("Start connecting border roads at {}".format(datetime.now()))
        for b in self.borders['border'].unique():
            # find all fitting edges for both countries
            edges_b1 = self.border_roads[(self.border_roads[roads_cn] == b[:2]) & (self.border_roads['border'] == b) &
                                         (self.border_roads[roads_type].isin(filter_types))]
            edges_b2 = self.border_roads[(self.border_roads[roads_cn] == b[2:]) & (self.border_roads['border'] == b) &
                                         (self.border_roads[roads_type].isin(filter_types))]
            # check if there are edges in both countries
            if (len(edges_b1) > 0) & (len(edges_b2) > 0):
                # get end nodes in both countries
                nodes_b11 = edges_b1[[roads_id, roads_fr]].rename(columns={roads_fr: node_id})
                nodes_b11['dir'] = 'from'
                nodes_b12 = edges_b1[[roads_id, roads_to]].rename(columns={roads_to: node_id})
                nodes_b12['dir'] = 'to'
                nodes_b1 = pd.concat([nodes_b11, nodes_b12]).groupby('node_id')[roads_id].agg('count').reset_index()
                if len(nodes_b1[nodes_b1[roads_id] == nodes_b1[roads_id].min()]) == 1:
                    nodes_b1 = nodes_b1[nodes_b1[roads_id] == nodes_b1[roads_id].min()+1]
                else:
                    nodes_b1 = nodes_b1[nodes_b1[roads_id] == nodes_b1[roads_id].min()]
                nodes_b21 = edges_b2[[roads_id, roads_fr]].rename(columns={roads_fr: node_id})
                nodes_b21['dir'] = 'from'
                nodes_b22 = edges_b2[[roads_id, roads_to]].rename(columns={roads_to: node_id})
                nodes_b22['dir'] = 'to'
                nodes_b2 = pd.concat([nodes_b21, nodes_b22]).groupby('node_id')[roads_id].agg('count').reset_index()
                if len(nodes_b2[nodes_b2[roads_id] == nodes_b2[roads_id].min()]) == 1:
                    nodes_b2 = nodes_b2[nodes_b2[roads_id] == nodes_b2[roads_id].min()+1]
                else:
                    nodes_b2 = nodes_b2[nodes_b2[roads_id] == nodes_b2[roads_id].min()]
                nodes_dir = pd.concat([nodes_b11, nodes_b12, nodes_b21, nodes_b22])
                # get coordinates
                nodes_b1 = nodes_b1.merge(nodes_int[[node_id, node_geo]], how='left', on=node_id)
                nodes_b1 = nodes_b1.set_geometry(node_geo)
                nodes_b2 = nodes_b2.merge(nodes_int[[node_id, node_geo]], how='left', on=node_id)
                nodes_b2 = nodes_b2.set_geometry(node_geo)
                # get the two nearest nodes for each node
                nearest_b1 = ckdnearest(nodes_b1, nodes_b2, node_id, 2)
                nearest_b2 = ckdnearest(nodes_b2, nodes_b1, node_id, 2)
                # find nearest node pairs (nearest nodes match) and remove duplicate connections
                pairs_b = nearest_b1.merge(nearest_b2, left_on=[node_id, '{}_near'.format(node_id)],
                                           right_on=['{}_near'.format(node_id), node_id],
                                           how='inner')[['{}_x'.format(node_id), '{}_y'.format(node_id), 'distance_x']]
                # find correct direction
                pairs_b = pairs_b.merge(nodes_dir[[node_id, 'dir']], left_on='{}_x'.format(node_id), right_on=node_id,
                                        how='left')
                pairs_b = pairs_b.merge(nodes_dir[[node_id, 'dir']], left_on='{}_y'.format(node_id), right_on=node_id,
                                        how='left')
                pairs_b = pairs_b.loc[:, ~pairs_b.columns.duplicated()]
                pairs_b = pairs_b[pairs_b['dir_x'] != pairs_b['dir_y']]
                pairs_b['node_from'] = [
                    row['{}_x'.format(node_id)] if row['dir_x'] == 'to' else row['{}_y'.format(node_id)] for i, row in
                    pairs_b.iterrows()]
                pairs_b['node_to'] = [
                    row['{}_x'.format(node_id)] if row['dir_x'] == 'from' else row['{}_y'.format(node_id)] for i, row in
                    pairs_b.iterrows()]
                # remove duplicate connections (if a node is paired with multiple other nodes)
                pairs_b = pairs_b.sort_values(by='distance_x')
                pairs_b = pairs_b.drop_duplicates(subset='node_from')
                pairs_b = pairs_b.drop_duplicates(subset='node_to')
                pairs_b = [[row['node_from'], row['node_to']] for i, row in pairs_b.iterrows()]
                # create lines between node pairs
                for l in pairs_b:
                    # get coordinates and create LineString
                    node_from = nodes_int.loc[nodes_int[node_id] == l[0], 'geom']
                    node_to = nodes_int.loc[nodes_int[node_id] == l[1], 'geom']
                    line = LineString([[node_from.x, node_from.y], [node_to.x, node_to.y]])
                    # add to Lines GDF
                    lines.loc[len(lines) + 1] = [id_st, l[0], l[1], 999, b + '0', 0, 50, 0, line, np.nan]
                    id_st += 1
                self.dict_borders.update({b: len(pairs_b)})
            else:
                print('No border connection for {}'.format(b))
                self.dict_borders.update({b: 0})

        # concat border crossings and rest of the network
        self.network_int = pd.concat([self.network, lines])
        print("Finished connecting border roads at {}".format(datetime.now()))
