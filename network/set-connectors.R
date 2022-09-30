# ULTIMO Preparation: DEFINE CONNECTORS TO NETWORK FOR EACH CELL

require (rgdal)
require (rgeos)
require (dplyr)
require (raster)
require (reshape2)
require (maptools)
require (sp)
library(RPostgreSQL)
library(config)
library(sf)

# connect to db model_ultimo
cred = get("ultimo")
con = dbConnect(cred$driver, user= cred$uid, password=cred$pwd, host=cred$server, port=cred$port, dbname=cred$database)

## Read population grid and zones ----

pop.250<-readOGR("V:\\Projekte\\ELK\\network\\4326pop-points1000eu.gpkg")#, p4s="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
pop.250@data$id <- c(1:nrow(pop.250@data))

cells.eu.base <- st_read(con, query = 'SELECT * FROM links_sp."cells_eu_cen"')
cells.eu.base <- as(cells.eu.base, "Spatial")

## Determine number of connectors ----
cells.eu.base$centers <- ifelse(cells.eu.base$flaeche < 100, 3,
                                ifelse(cells.eu.base$flaeche < 750, 4, 
                                       ifelse(cells.eu.base$flaeche < 1500, 5, 6)))

pop.250.geo.cells <- cbind(pop.250@data,over(pop.250,cells.eu.base))
key.250 <- subset(pop.250.geo.cells, is.na(pop.250.geo.cells$pop)==FALSE)
key.250 <- key.250 %>% group_by(eu_id) %>% mutate(poprank = order(order(pop, decreasing=TRUE)))
key.250 <- sp::merge(pop.250, key.250[c("id", "eu_id", "poprank")], by="id") #TODO check id column for pop.250

t<-data.frame(coords.x1=double(),coords.x2=double(), zone=factor(), clusterkm=factor(), pop=double(), linkweight=double())

it<-0
for (zoneid in (unique(key.250$eu_id))){
  it<-it+1 
  c <- cells.eu.base@data$centers[cells.eu.base@data$eu_id == zoneid]
  test.250 <- subset(key.250,key.250@data$eu_id==zoneid)
  rt<-try(length(test.250 [,1]))
  if(class(rt) %in% 'try-error') {print("e");next} 
    else{test.250 <- subset(test.250,test.250@data$poprank <=rt)}
  km<-try(kmeans(test.250@coords, centers = c),silent = TRUE)
  if(class(km) %in% 'try-error') {#next
    print("E");next} else {
      u<-cbind(km$centers, zoneid)
      erg<-as.data.frame(cbind(test.250@coords,km$cluster, test.250@data$pop,test.250@data$poprank))
      colnames(erg)[c(3:5)] <- c("clusterkm","pop","poprank")
      erg$is_link <-0
      topskm<- array(aggregate(poprank~clusterkm,erg,min) [,2])
      erg$is_link [erg$poprank %in% topskm] <-1 
      erg$zone <- zoneid
      u<-subset(erg,erg$is_link==1) [c("coords.x1","coords.x2", "zone", "clusterkm")]
      u<-base::merge(u,aggregate(pop~clusterkm,erg,sum),by="clusterkm")
      u$linkweight <- u[,5]/sum(u[,5])
      
      t<- rbind(t,u)
      print(zoneid)
    }
  
}

# t: write csv, import into Qgis with x and y coordinates
dbWriteTable(con, c("input","connectors_eu"), value = t, overwrite=T, row.names = FALSE)

# convert to spatial
t2 <- t
coordinates(t2) <- c("coords.x1", "coords.x2")
st_write(st_as_sf(t2), con, c('links_sp', 'connectors_eu_sp'))

#********************************************TBC**************************************

#tbc ----
# move to nearest nodes

# convert to shape ----

t$X <- t$X + 10000000
# t2 <- t
# coordinates(t2) <- c("coords.x1", "coords.x2")
# t2@proj4string<-CRS("+init=epsg:3035")

# writeOGR(t2, dsn = "_output", layer = "01-max-points", driver = "ESRI Shapefile", overwrite_layer = T, check_exists = T)

# SWITCH TO JUPYTER FOR MOVING TO NODES----
# read result

node.points <- readOGR("_output/points-to-nodes.gpkg")
node.points <- spTransform(node.points, CRS("+init=epsg:3035"))
#node.points@proj4string<-CRS("+init=epsg:3035")
node.points <- merge(node.points, t[,c('X','clusterkm')], by = 'X', all.x = TRUE)

# add connector column and missing connector points
t2 <- node.points
t2 <- t2[,c('NO', 'zone', 'linkweight', 'clusterkm', 'distance','node')]
t2$conn <- 0

# remove points that are too close to each other ----
t2$cut <- 0
q<- 0

control.points <- t2@data %>% count(zone)
control.points[control.points$n <= 1,] -> missing.points

for (ix in unique(t2$zone)){
  # get subset
  subset.t2 <- subset(t2, zone == ix)
  subset.t2$count <- 0
  subset.t2$mean.dis <- 0
  if (nrow(subset.t2) <= 2){
    q <- q+1
    dist.close <- F
  } else {
    dist.close <- T
  }
  
  while(dist.close) {
    distances <- data.frame(pointDistance(subset.t2, lonlat = FALSE))
    all.p <- sort((subset.t2@data$clusterkm))
    names(distances) <- all.p
    distances$index <- all.p
    all.p2 <- all.p[all.p > 0]
    
    # find points which are too close to the other points
    for (i in all.p2){
      #print(i)
      i <- as.character(i)
      mean.dis <- 9999
      count <- length(which(distances[i] <= 500 & distances[i] > 0))
      if (count > 0){
        # CHECK AGAIN!
        # k <- nrow(distances) - which(distances[i] < 500)
        # mean.dis <- sum(distances[i, k])/length(k)
        # rm
        mean.dis <- sum(distances[i])/(length(all.p)-1)
      }
      subset.t2@data$count[subset.t2@data$clusterkm == i] <- count
      subset.t2@data$mean.dis[subset.t2@data$clusterkm == i] <- mean.dis
      rm(count, mean.dis)
    }
    
    #TODO make sure there are at least two or three points left (city zones seem to have too little, use aggregate to check)
    # find one point which is too close to the rest of the points and not a connector
    m <- max(subset.t2@data$count[subset.t2$conn == 0])
    
    if (m == 0){
      dist.close <- F
    } else {
      pp <- subset.t2@data[subset.t2@data$count == m & subset.t2$conn == 0, ]
      if (nrow(pp) == 0){
        dist.close <- F
      } else {
        pp$weight.dis <- pp$linkweight * pp$mean.dis
        p.x <- pp$clusterkm[pp$weight.dis == min(pp$weight.dis)]
        weight.px <- sum(pp$linkweight[pp$clusterkm %in% p.x])
        t2@data$cut[t2@data$clusterkm %in% p.x & t2@data$zone == ix] <-1
        subset.t2@data$cut[subset.t2@data$clusterkm %in% p.x] <-1
        # add weight of removed point to nearest neighbor(s)
        nearest.neighbor.dis <- min(distances[distances$index == p.x, as.character(all.p[all.p != p.x])])
        p.y <- names(distances)[which(distances[distances$index == p.x,] == nearest.neighbor.dis)]
        t2@data$linkweight[t2@data$clusterkm %in% p.y & t2@data$zone == ix] <- t2@data$linkweight[t2@data$clusterkm %in% p.y & t2@data$zone == ix]+weight.px/length(p.y)
        subset.t2@data$linkweight[subset.t2@data$clusterkm %in% p.y] <- subset.t2@data$linkweight[subset.t2@data$clusterkm %in% p.y] + (weight.px/length(p.y))
        rm(pp, p.x, p.y, nearest.neighbor.dis, weight.px)
      }
      
    }
    rm(m, distances, all.p, all.p2)
    subset.t2 <- subset(subset.t2, cut == 0)
    if (nrow(subset.t2) <= 3){
      dist.close <- F
    }
  }
}
print(q)
rm(subset.t2)
rm(i, ix, dist.close, q)

t3 <- t2[t2@data$cut == 0,]
print(sum(t2$conn) == sum(t3$conn))
control.points <- t3@data %>% count(zone)
control.points[control.points$n <= 1,] -> missing.points
print(nrow(missing.points))
rm(control.points, missing.points)

# export
writeOGR(t3, dsn = "_output", layer = "pick-up-DE", driver = "ESRI Shapefile", overwrite_layer = T)

rm(t, t2, t3, t.base)


# check double connectors ----
t4$alt.conn <- 0

control.conn <- t4@data[t4@data$conn == 1,] %>% count(NO)
duplicates <- control.conn[control.conn$n > 1,]
print(nrow(duplicates))

check.dup <- T

while (check.dup){
  dub.conn <- unique(duplicates$NO)
  
  for (conn in dub.conn){
    df.conn <- t4@data[t4@data$NO == conn,]
    # keep zone if a) no other node point or b) shortest distance
    zones.keep <- as.character(df.conn$zone[df.conn$node == 0])
    if (length(zones.keep) == 0){
      df.conn$disweight <- df.conn$distance*df.conn$linkweight
      zones.keep <- as.character(df.conn$zone[df.conn$disweight == min(df.conn$disweight)])
    }
    zones.move <- as.character(df.conn$zone[!df.conn$zone %in% zones.keep])
    # find alternative connector points based on linkweight
    for (zone in zones.move){
      df.alternative <- t4@data[t4@data$zone == zone & t4@data$conn == 0,]
      conn.alternative <- df.alternative$NO[df.alternative$linkweight == max(df.alternative$linkweight)]
      # settings for new connector
      t4@data$conn[t4@data$NO == conn.alternative & t4@data$zone == zone] <- 1
      t4@data$alt.conn[t4@data$NO == conn.alternative & t4@data$zone == zone] <- 1
    }
    t4@data$conn[t4@data$NO == conn & t4@data$zone %in% zones.move] <- 0
  }
  prev.dup <- nrow(duplicates)
  control.conn <- t4@data[t4@data$conn == 1,] %>% count(NO)
  duplicates <- control.conn[control.conn$n > 1,]
  print(nrow(duplicates))
  if (nrow(duplicates) == 0){
    check.dup <- F
  } else {
    check.dup <- nrow(duplicates) != prev.dup
  }
}
t4$alt.conn[t4$alt.conn == 1 & t4$conn == 0] <- 0

print(paste(sum(t4$alt.conn), "connectors were moved"))
#TODO why do some zones lose connectors? Check with aggregate and merge
print(sum(t4$conn)==sum(t3$conn))

rm(all.zones, check.dup, conn, conn.alternative, dub.conn, prev.dup, zone, zones.keep, zones.move)
rm(duplicates, control.conn, df.alternative, df.conn)


# Export POI as shp and connectors as csv ----
t5 <- t4@data
t5 <- t5[t5$conn == 1, c("NO", "zone", "linkweight","distance", "area_km", "buffer", "alt.conn")]
write.csv(t5, "_output/connectors-rs-DE.csv")

t4@data <- t4@data[,c("NO", "zone", "linkweight", "area_km", "buffer")]
# export final with area
writeOGR(t4, dsn = "_output", layer = "pick-up-final-DE", driver = "ESRI Shapefile", overwrite_layer = T, check_exists = T)

rm(t3, t4, t5, node.points)
