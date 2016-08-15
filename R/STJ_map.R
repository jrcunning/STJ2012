# St. John Map

library(rgdal)
stjmap <- readOGR("data/coastline/zipfolder", layer="Coastal_Coastline_line", verbose=F)  # units are in meters
stjmap <- spTransform(stjmap, CRS("+proj=longlat +datum=WGS84"))  # transform to lat/long

sitecoords <- data.frame(
  name=c("Booby Rock", "Tektite", "Kiddel", "Groot Pan", "Yawzi Point", "Anna's Point", "Haulover Bay"),
  lat=c(18.30327, 18.31057, 18.30770, 18.31067, 18.31515, 18.36710, 18.35008),
  lon=c(-64.71057, -64.72203, -64.71398, -64.71813, -64.72513, -64.73337, -64.67933)
)

viers <- c(-64.722697, 18.322277)

#pdf(file="figures/map.pdf", height=3, width=4)
par(mar=c(0,0,0,0))
plot(stjmap, lwd=1, col="gray40",
     xlim=c(-64.75, -64.68), ylim=c(18.303, 18.370))
box()
text(-64.74, 18.34, expression(italic("       St. John\nU.S. Virgin Islands")), col="gray40")
axis(side=3, at=seq(-64.75, -64.65, 0.01), line=0, tck=0.02, labels=FALSE)
axis(side=3, at=seq(-64.75, -64.65, 0.04), line=-2, lwd=0, cex.axis=0.6)
axis(side=2, at=seq(18.3, 18.4, 0.01), line=0, tck=0.02, labels=FALSE)
axis(side=2, at=seq(18.3, 18.4, 0.02), line=-2, lwd=0, cex.axis=0.6)
with(sitecoords, {
  points(lon, lat, pch=21, bg="black")
  text(lon, lat, name, pos=c(4,2,4,4,2,1,3))
})
points(viers[1], viers[2], pch=2, col="gray40", cex=0.8)
text(viers[1], viers[2], expression(italic("VIERS")), col="gray40", pos=4, cex=0.7)
#dev.off()



