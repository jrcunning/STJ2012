# St. John Map

library(rgdal)
stjmap <- readOGR("data/coastline/zipfolder", layer="Coastal_Coastline_line")  # units are in meters
stjmap <- spTransform(stjmap, CRS("+proj=longlat +datum=WGS84"))  # transform to lat/long

sitecoords <- data.frame(
  name=c("Booby Rock", "Tektite", "Kiddel", "Groot Pan", "Yawzi Point", "Anna's Point", "Haulover Bay"),
  lat=c(18.30327, 18.31057, 18.30770, 18.31067, 18.31515, 18.36710, 18.35008),
  lon=c(-64.71057, -64.72203, -64.71398, -64.71813, -64.72513, -64.73337, -64.67933)
)

viers <- c(-64.722697, 18.322277)

pdf(file="figures/map.pdf", height=3, width=4)
par(mar=c(0,0,0,0))
plot(stjmap, lwd=1, col="gray40",
     xlim=c(-64.75, -64.68), ylim=c(18.303, 18.370))
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
dev.off()




text(x=seq(-157.88, -157.73, 0.04), y=21.403, labels=seq(-157.88, -157.73, 0.04), cex=0.5)
axis(side=3, at=seq(-157.88, -157.80, 0.04), line=0, tck=0.01, labels=FALSE)
text(x=seq(-157.88, -157.80, 0.04), y=21.527, labels=seq(-157.88, -157.80, 0.04), cex=0.5)
axis(side=2, at=seq(21.41, 21.6, 0.04), tck=0.01)
text(x=-157.892, y=seq(21.41, 21.5, 0.04), labels=seq(21.41, 21.5, 0.04), cex=0.5)
axis(side=4, at=seq(21.41, 21.6, 0.04), tck=0.01)
text(x=-157.718, y=seq(21.41, 21.6, 0.04), labels=seq(21.41, 21.6, 0.04), cex=0.5)
box()
points(reefcoords[,c(2,1)], pch=21, cex=1.5, col="black", bg=reefcols)
text(reefcoords[,c(2,1)], labels=c("Reef 44", "Reef 25", "HIMB"), pos=4)
par(new=T, mar=c(0,0,9,12.8))
plot(HI, xlim=c(-158.3, -157.6), ylim=c(21.35, 21.6), lwd=0.4, col="gray", bg="white")  # Oahu
rect(-157.9, 21.39, -157.7, 21.53)
box()
par(new=T, mar=c(9,12.8,0,0))
plot(img2)
box()
dev.off()
