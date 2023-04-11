# Cost matrix

The cost of travel in the form of travel times and distances is a key element of ULTImodel's 
gravity model. The cost matrix is created using `generationData.Matrices`.

Preparation: read connector points `GeoDataFrame`.
```python
import generationData.Matrices as Matrices
import geopandas as gpd

connectors = gpd.GeoDataFrame.from_file('connectors.gpkg')
connectors.head()

>>> output

	node_id 	geometry 	                nuts_id 	c_n 	weight
0 	1095 	    POINT (12.74977 55.60094) 	DK011 	    0 	    0.000048
1 	1014 	    POINT (12.51936 55.70192) 	DK011 	    1 	    0.785301
2 	1078 	    POINT (12.60314 55.62947) 	DK011 	    2 	    0.214651
3 	890 	    POINT (12.39536 55.67467) 	DK012 	    3 	    0.493351
4 	777 	    POINT (12.23742 55.65711) 	DK012 	    4 	    0.103410  
```

Initialize the class `Matrix` using the connectors

```python
mx = Matrices.Matrix(conn=connectors, zone_col='nuts_id', conn_geo='geometry', id_col='c_n')
```

Cost matrix between all connector points using OSRM requests

```python
mx.osrm_request_nodes()

>>> output
43 connector point coordinates
start requests: 2023-03-17 17:50:51.312626

100%|████████████████████████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  2.82it/s]

total time requests: 0.36452651023864746 s (0.01 min)
2023-03-17 17:50:51.677153

mx.all_mtx.shape

>>> (43, 43, 2)
```

Aggregate to zone matrix using connector weights

```python
zone_matrix, zone_ids = mx.transform_to_taz()

>>> output

100%|██████████████████████████████████████████████████████████████████████████████████| 43/43 [00:01<00:00, 36.68it/s]

start aggregating zone matrix: 2023-03-17 17:53:32.124560

100%|█████████████████████████████████████████████████████████████████████████████████| 11/11 [00:00<00:00, 120.63it/s]

total time zones: 0.09325337409973145 s (0.0 min)
2023-03-17 17:53:32.217814

zone_matrix.shape

>>> (11, 11, 2)
```

`zone_ids` contains the index of each TAZ in the cost matrix.

```python
zone_ids.head()

>>> output

	zone 	id
0 	DK011 	0
1 	DK012 	1
2 	DK013 	2
3 	DK014 	3
4 	DK021 	4
```