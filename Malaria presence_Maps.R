# load required packages:
library(terra)  # spatial (raster + vector) data analysis and mapping
library(geodata)  # spatial data download
library(rnaturalearth)  # more spatial data download
library(car)  # 'pointLabel' function

## Load global country data
countries <- terra::vect("../inputs/World_Countries/TM_WORLD_BORDERS-0.3.shp")


## plot global map (countries)
plot(countries) # map may take a bit to appear (it's a lot of vertices)

# Subset for South Asia
plot(countries, "region")

papua <- subset(countries, countries$name =="Papua New Guinea")
papua
plot(papua)

southAsia <- subset(countries, countries$region == "142")
southAsia
plot(southAsia)


SouthAsiaR <- combineGeoms(southAsia, papua)
SouthAsiaR
plot(SouthAsiaR)

## Plot presence/absence P.falciparum malaria for Asia

plot(SouthAsiaR, col = "gray88", xlim =c(60, 150), background = "lightblue1")  # awkward-looking map because of islands on the other side of the 180-degrees meridian


pfalcAsia <- read.csv("../inputs/Pfalciparum_Asia.csv")
pfalcAsia_points <- terra::vect(pfalcAsia, geom = c("longitude", "latitude"), keepgeom = TRUE)
pfalcAsia_points

plot(pfalcAsia_points, y = "pf_pa",  pch = c(20,20), col = c("magenta","blue"), cex = 1, add=TRUE)

sbar(type = "bar", divs=4, ticks = TRUE, below = "km")
text(145, 51, "C", cex=3)



## Plot presence/absence P. vivax malaria for Asia


pvivAsia <- read.csv("../inputs/Pvivax_Asia.csv")
pvivAsia_points <- terra::vect(pvivAsia, geom = c("longitude", "latitude"), keepgeom = TRUE)
pvivAsia_points

plot(SouthAsiaR, col = "gray88", xlim =c(60, 150), background="lightblue1")

plot(pvivAsia_points, y = "pv_pa", pch = c(20,20), col = c("magenta","blue"), cex = 1, add=TRUE)

sbar(type = "bar", divs=4, ticks = TRUE, below = "km")
text(145, 51, "D", cex=3)



## Plot presence/absence P. falciparum malaria for Africa

head(countries)

Africa <- subset(countries, countries$region == "2")
Africa
plot(Africa, col = "gray88", background = "lightblue1")


pfalcAfrica <- read.csv("../inputs/Pfalciparum_Africa.csv")
pfalcAfrica_points <- terra::vect(pfalcAfrica, geom = c("longitude", "latitude"), keepgeom = TRUE)
pfalcAfrica_points

plot(pfalcAfrica_points, y = "pf_pa", pch = c(20,20), col = c("magenta","blue"), cex = 1, add=TRUE) ##pch=21  and alpha = 0-1(0.5)

sbar(type = "bar", divs=4, ticks = TRUE, below = "km")
text(60, 32, "A", cex=3)

legend(-21, -20, legend = c("Malaria absence","Malaria presence"),pch = c(20,20), xpd=NA, bg="white", 
       col=c("magenta","blue"), cex = 1.4)

north(xy=c(55, -37), type = 2, cex = 2.4)


## P. Plot presence/absence vivax Africa

head(countries)

Africa <- subset(countries, countries$region == "2")
Africa
plot(Africa, col = "gray88", background = "lightblue1")


pvivAfrica <- read.csv("../inputs/Pvivax_Africa.csv")
pvivAfrica_points <- terra::vect(pvivAfrica, geom = c("longitude", "latitude"), keepgeom = TRUE)
pvivAfrica_points

plot(pvivAfrica_points, y = "pv_pa", pch = c(20,20), col = c("magenta", "blue"), cex = 1, add=TRUE)

sbar(type = "bar", divs=4, ticks = TRUE, below = "km")
text(60, 32, "B", cex=3)

text(-12, -20, "Malaria Absence", font=2, col="magenta", cex = 1.4, halo = TRUE)
text(-12, -25, "Malaria Presence", font=2, col="cyan", cex= 1.4, halo = TRUE)

north( "bottomright", type = 15, cex = 2.2)







