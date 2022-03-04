################################################################################
# Calculating wind fetch.
# Heavily based on https://github.com/blasee/fetchR
# however blasee/fetchR failed to correctly produce the vectors for site 
# in Japan.
# 
# 2022 March 4
# Greg Nishihara

# Load packages. ###############################################################
library(tidyverse)
library(ggpubr)
library(sf)
library(magick)
library(ggrepel)
library(showtext)

# やっぱり Noto Sans がきれい。
if(!any(str_detect(font_families(), "notosans"))) {
  font_add_google("Noto Sans","notosans")
}
# 図のフォントがからだったので、ここで修正した
# １）theme_set() をつかってデフォルトのフォントをかえる
# ２）ggplot() の theme() からとんとの指定をはずす。
theme_pubr(base_size = 6, base_family = "notosans") |> theme_set()
showtext_auto()
Sys.setlocale("LC_TIME", "en_US.UTF-8") # This is to set the server time locate to en_US.UTF-8

################################################################################
# Define some functions to determine the intersection of the vectors used to 
# determine fetch. The original version in fetchR did not work for this map
# because it returned the worn intersections.
################################################################################
## Calculate intersections and identify closest shoreline for those 
## vectors that hit land

calc_intersection = function(origin, map_layer, fetch_limits) {
  # This function creates a line from the origin to the fetch limit.
  # Then it will determine the points where the line and the map polygons
  # intersects. There can be more than one intersection if the line 
  # has sufficient length. The trick is to determine the point closest to the
  # origin so that it can be used as the end point for the fetch. If
  # there is no intersection, then the end point for the fetch is the 
  # fetch limit.
  fl = fetch_limits[1:2]
  or = st_coordinates(origin)
  X = matrix(c(or, fl), ncol = 2, byrow = T)
  lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
  int = st_intersection(lst, map_layer) 
  
  xy  = st_coordinates(int) |> as.matrix()
  distance = function(x1, x2) {
    sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2 )
  }
  # If there are not intersections, then the xy matrix has dimensions [0,2]
  j = dim(xy)
  if(j[1] > 0) {
    dis = apply(xy, 1, distance, x2 = fl)
    max_dis = max(dis)
    n = which(near(dis, max_dis)) 
    xy = xy[n,]
    X = matrix(c(or, xy[1:2]), ncol = 2, byrow =T)
    lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
    len = st_length(lst) 
    z = lst |> st_as_sf() |> mutate(length = len)
  } else {
    X = matrix(c(or, fl), ncol = 2, byrow =T)
    lst = st_linestring(X) |> st_sfc(crs = st_crs(map_layer))
    len = st_length(lst) 
    z = lst |> st_as_sf() |> mutate(length = len) 
  }
  z
}

calc_circle = function(map_layer, max_dist=30, n_vectors = 9) {
  # Calculate the fetch limits.
  delta_theta =  360 / (n_vectors * 4)
  theta = seq(0, 360, by = delta_theta) |> head(n = -1)
  max_dist = max_dist * 1000 # Convert from km to m.
  max_dist = units::set_units(max_dist, "m")
  d_bff = st_buffer(map_layer, dist = max_dist, nQuadSegs = n_vectors) 
  fetch_ends = st_coordinates(d_bff) |> head(n = -1) 
  list(d_bff = d_bff, fetch_limits = fetch_ends[order(theta), ] )
}
################################################################################
# Prepare data set #############################################################
# GPS coordinates to determine wind fetcj ######################################
matsushimagps = matrix(rev(c(38.34549669653925, 141.0807915733725)), ncol = 2)
hirotagps     = matrix(rev(c(39.02402594131857, 141.78725806724896)), ncol = 2)
bisegps       = matrix(rev(c(26.704302654710496, 127.85974269102186)), ncol = 2)
arikawagps    = matrix(rev(c(32.98827976845565, 129.11838896005543)), ncol = 2)
# arikawagps    = matrix(rev(c(32.99602966223495, 129.11335509246598)), ncol = 2) # For Fetch calcs
tainouragps   = matrix(rev(c(32.95134175383013, 129.1096027426365)), ncol = 2)
omuragps      = matrix(rev(c(32+52/60+11.9/60/60, 129+58/60+24.5/60/60)), ncol = 2)

gps_info = rbind(matsushimagps, hirotagps, bisegps, arikawagps, tainouragps, omuragps) |> 
  as_tibble() |> 
  mutate(name = c("matsushimagps", "hirotagps", "bisegps", "arikawagps", "tainouragps", "omuragps")) |> 
  mutate(label = str_to_sentence(str_remove(name, pattern = "(gps)"))) |> 
  rename(long = V1, lat = V2)

gps_info = gps_info |> 
  mutate(label2 = str_to_sentence(label)) |> 
  mutate(label2 = str_glue("{label2} {ifelse(str_detect(label2, 'Bise'), 'Point', 'Bay')}"))

# Prepare the Coordinate Reference System to be EPSG:4326 (Which is WGS 84)
# See st_crs(4326) for details
gps_info = gps_info |> select(long, lat, name) |> st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant")

# Load the map shape files #####################################################
# The map uses the ITRF94 system (st_crs(map_poly))
# gsi_low = read_sf("~/Lab_Data/Japan_map_data/GSI/coastl_jpn.shp")
# gsi_low = read_sf("~/Lab_Data/Japan_map_data/GADM_old/JPN_adm1.shp")
map_poly = read_sf("~/Lab_Data/Japan_map_data/GSI/polbnda_jpn.shp")
map_poly = map_poly |> select(nam, geometry)

# Convert the CRS to EPSG:2450 #################################################
map_poly = st_transform(map_poly, st_crs(2450))
gps_info  = st_transform(gps_info, st_crs(2450))

################################################################################
# Do the analysis one location at a time. ######################################
ptsize = 1
max_dist = 10 # In km
n_vectors = 3*9 # The number of vectors in every quadrant.

location = "Hirota Bay"
polygon_layer = subset(map_poly, str_detect(nam, "Iwate")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "hiro"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()
max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p1 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

location = "Matsushima Bay"
polygon_layer = subset(map_poly, str_detect(nam, "Miyag")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "matsu"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()
max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p2 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

location = "Bise Point"
polygon_layer = subset(map_poly, str_detect(nam, "Okinawa")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "bise"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()
max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p3 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

location = "Omura Bay"
polygon_layer = subset(map_poly, str_detect(nam, "Nagasaki")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "omura"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()
max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p4 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

location = "Arikawa Bay"
polygon_layer = subset(map_poly, str_detect(nam, "Nagasaki")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "arik"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()

max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p5 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

location = "Tainoura Bay"
polygon_layer = subset(map_poly, str_detect(nam, "Nagasaki")) |> st_union() 
site_layer    = subset(gps_info, str_detect(name, "tain"))
fetch_limits = calc_circle(site_layer, max_dist = max_dist, n_vectors = n_vectors)
fout = fetch_limits$fetch_limits |> as_tibble() |> 
  mutate(fe  = map2(X,Y,function(x,y) cbind(x,y))) |> 
  mutate(geometry = map(fe, calc_intersection, origin = site_layer, map_layer = polygon_layer))
fout = fout |> select(geometry) |>  unnest(geometry) |> st_as_sf()
temp_layer = st_crop(polygon_layer, st_bbox(fetch_limits$d_bff))
mean_fetch = fout |> pull(length) |> mean() |> as.numeric()
sd_fetch = fout |> pull(length) |> sd() |> as.numeric()

max_fetch = fout |> pull(length) |> as.numeric()
man_n = sum(near(max_fetch, max_dist * 1000))
tot_n = length(max_fetch)

p6 = ggplot() + 
  geom_sf(data = temp_layer, color = NA) +
  geom_sf(data = fout) +
  geom_sf(data = site_layer, color = "red", size = ptsize) +
  labs(title = str_glue("The mean ± sd fetch for {location} is {format(mean_fetch, digits = 4)} ± {format(sd_fetch, digits = 4)} m."),
       subtitle = str_glue("{man_n} out of {tot_n} vectors were at the upper limit."))

library(patchwork)
(p1 + p2 + p3) / (p4 + p5 + p6)

pdfname = "Determine_fetch.pdf"
pngname = str_replace(pdfname, "pdf", "png")
ggsave(pdfname, width = 4*80, height = 3*80, units = "mm")
img = image_read(pdfname, density = 300)
img |> image_write(pngname, format = "png")























