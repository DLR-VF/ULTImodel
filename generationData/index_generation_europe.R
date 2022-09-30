# Attraction index for all European NUTS3-zones

df_zones <- read.csv("_data/europe_zones.csv")
df_zones <- df_zones[,c("FID","CNTR_CODE","flaeche","pop","X_mean","ind_area_sum","ind_area_count")]

# add IDs for Visum

ids <- read.csv("_data/list_zones_europe.att", skip = 12, sep = ";")[,c("X.ZONE.NO","NAME")]
df_zones <- merge(ids, df_zones, by.y = "FID", by.x = "NAME")
names(df_zones) <- c("NUTS_ID","id","CNTR_CODE","area","pop","pop_dens","ind_area","ind_count")

rm(ids)

# index values for each country

df_zones$index_pop <- 0
df_zones$index_area <- 0
df_zones$index_count <- 0

for (c in unique(df_zones$CNTR_CODE)){
  mean_pop   <- mean(df_zones$pop[df_zones$CNTR_CODE==c], na.rm = TRUE)
  mean_area  <- mean(df_zones$ind_area[df_zones$CNTR_CODE==c], na.rm = TRUE)
  mean_count <- mean(df_zones$ind_count[df_zones$CNTR_CODE==c], na.rm = TRUE)
  
  df_zones$index_pop[df_zones$CNTR_CODE == c] <- df_zones$pop[df_zones$CNTR_CODE == c]/mean_pop
  df_zones$index_area[df_zones$CNTR_CODE == c] <- df_zones$ind_area[df_zones$CNTR_CODE == c]/mean_area
  df_zones$index_count[df_zones$CNTR_CODE == c] <- df_zones$ind_count[df_zones$CNTR_CODE == c]/mean_count
}

df_zones$index_pop[is.na(df_zones$index_pop)] <- 0
df_zones$index_area[is.na(df_zones$index_area)] <- 0
df_zones$index_count[is.na(df_zones$index_count)] <- 0

df_zones$index_country <- df_zones$index_pop**0.5*df_zones$index_area**0.25*df_zones$index_count**0.25

# index values on european level

mean_pop <- mean(df_zones$pop, na.rm = TRUE)
mean_area <- mean(df_zones$ind_area, na.rm = TRUE)
mean_count <- mean(df_zones$ind_count, na.rm = TRUE)

df_zones$index_pop_eu <- df_zones$pop / mean_pop
df_zones$index_area_eu <- df_zones$ind_area / mean_area
df_zones$index_count_eu <- df_zones$ind_count / mean_count

df_zones$index_pop_eu[is.na(df_zones$index_pop_eu)] <- 0
df_zones$index_area_eu[is.na(df_zones$index_area_eu)] <- 0
df_zones$index_count_eu[is.na(df_zones$index_count_eu)] <- 0

df_zones$index_eu <- df_zones$index_pop_eu**0.5*df_zones$index_area_eu**0.25*df_zones$index_count_eu**0.25

# Export as csv

write.table(df_zones,"_data/europe_cells_index.csv", sep = ",", col.names = TRUE, row.names = FALSE)
rm(df_zones)
rm(mean_area, mean_count, mean_pop)
