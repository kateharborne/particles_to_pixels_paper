# Kate Harborne 24/10/23
# ============================================================================ #
# Making images: 
# Making some pretty plots of the galaxies and their kinematic maps. Also, some
# visualisations of the observations following voronoi binning.
# Figure 3 in the paper. 
# ============================================================================ #

library(SimSpin)
library(magicaxis)

# PATHS ========================================================================
ss_loc = "simspin_files/"
sa_loc = "sim_analysis/"
plot_loc = "plots/"
data_loc = "data/"
cube_loc = paste0(data_loc, "/cubes/flux_obs/")

source("functions/plot_images.R")
source("functions/plot_array.R")

# Loading data =================================================================

bulge = readRDS(paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"))
disk = readRDS(paste0(ss_loc, "disk_650_control_simspin_file.Rdata"))

bulge_photometry = readRDS(paste0(data_loc, "bulge_sm650_control_photo_values_500.Rdata"))
disk_photometry = readRDS(paste0(data_loc, "disk_sm650_control_photo_values_500.Rdata"))

bulge_analysis = readRDS(paste0(sa_loc,"control_bulge.Rdata"))
disk_analysis = readRDS(paste0(sa_loc,"control_disk.Rdata"))

bulge_ids = which(abs(bulge$star_part$x) <= 50 & 
                    abs(bulge$star_part$y) <= 50 & 
                    abs(bulge$star_part$z) <= 50)[sample(5128024, 1e5)]
disk_ids = which(abs(disk$star_part$x) <= 50 & 
                   abs(disk$star_part$y) <= 50 & 
                   abs(disk$star_part$z) <= 50)[sample(6499784, 1e5)]

bulge_ifu = readRDS(paste0(cube_loc, "bulge_sm650_control_ifu_obs.Rdata"))
disk_ifu = readRDS(paste0(cube_loc, "disk_sm650_control_ifu_obs.Rdata"))

# Initialising plot regions ====================================================
# Centres 
crow1 = c(0,0.6); crow2 = c(0.4,1);
ccol3 = c(0.2,0.8);

# Outers - col overlap 0.35 with centre, col overlap 0.2 with each other
h = plot_array(0.4, n_plots = 4)
w = plot_array(0.2, n_plots = 6)

ocol1 = w[1,]#c(0,0.23); 
ocol2 = w[2,]#c(0.13,0.43)
ocol3 = w[5,]#c(0.57,0.82); 
ocol4 = w[6,]#c(0.7, 1)

orow1 = h[2,]#c(0.11,0.56); 
orow2 = h[3,]#c(0.33,0.78); 
orow3 = h[4,]#c(0.6,1)

orow4 = h[1,]#c(0,0.45); 
orow5 = h[2,]#c(0.22,0.67); 
orow6 = h[3,]#c(0.44,0.89)

bulge_ifu_region = c(0,0.4,0.2,1)

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")

cex_val = 1.2

# Plot =========================================================================
pdf.options(fonts="Arial")
png(file=paste0(plot_loc, "mock_observation_example.png"), 
    width = 12, height = 9, units = "in", res=200)

# Central images of galaxies ---------------------------------------------------
# BULGE
par(fig=c(ccol3,crow2), mar=c(5.1,4.1,4.1,2.1), family="sans", pty="s")

magplot(0, type = "n", xlim = c(-50,50), ylim = c(-50,50), 
        grid = F, xaxt="n", yaxt = "n")
magaxis(side=c(1,2,3,4), labels = F)
points(bulge$star_part$x[bulge_ids], bulge$star_part$y[bulge_ids],
       pch = ".", col = rgb(51/255, 34/255, 136/255, alpha = 0.2))
lines(c(-40, -30), c(-43,-43), lwd=3, col="black")
text(-35, -36, "10 kpc", font=2, cex = cex_val)
text(24, 42, "Bulge Model", font=2, cex = cex_val)
text(24, 35, expression(italic("t")*" = 10 Gyr"), font=1, cex = cex_val)

plotrix::draw.ellipse(0, 0, 
                      bulge_photometry$major_axis_a_kpc, 
                      bulge_photometry$minor_axis_b_kpc, 
                      bulge_photometry$primary_axis+90, 
                      border = "#DDCC77", density = NULL, lwd=2)

plotrix::draw.ellipse(0, 0, 
                      bulge_analysis$HalfMassProperties$Radius_a, 
                      bulge_analysis$HalfMassProperties$Radius_b,  
                      bulge_photometry$primary_axis+90, 
                      border = "#882255", density = NULL, lwd=2)

plotrix::draw.ellipse(0, 0,                       
                      bulge_ifu$observation$aperture_size/2, 
                      bulge_ifu$observation$aperture_size/2,  
                      90,  
                      border = "black", density = NULL, lwd=3)



# DISK
par(fig=c(ccol3,crow1), mar=c(5.1,4.1,4.1,2.1), family="sans", new = T, pty="s")
magplot(0, type = "n", xlim = c(-50,50), ylim = c(-50,50), 
        grid = F, xaxt="n", yaxt = "n")
magaxis(side=c(1,2,3,4), labels = F)
points(disk$star_part$x[disk_ids], disk$star_part$z[disk_ids],
       pch = ".", col = rgb(51/255, 34/255, 136/255, alpha = 0.4))
lines(c(-40, -30), c(-43,-43), lwd=3, col="black")
text(-35, -36, "10 kpc", font=2, cex = cex_val)
text(24, 42, "Disk Model", font=2, cex = cex_val)
text(24, 35, expression(italic("t")*" = 10 Gyr"), font=1, cex = cex_val)

plotrix::draw.ellipse(0, 0, 
                      disk_photometry$major_axis_a_kpc, 
                      disk_photometry$minor_axis_b_kpc, 
                      disk_photometry$primary_axis+90, 
                      border = "#DDCC77", density = NULL, lwd=2)

plotrix::draw.ellipse(0, 0, 
                      disk_analysis$HalfMassProperties$Radius_a, 
                      disk_analysis$HalfMassProperties$Radius_c,  
                      disk_photometry$primary_axis+90, 
                      border = "#882255", density = NULL, lwd=2)

plotrix::draw.ellipse(0, 0, 
                      disk_ifu$observation$aperture_size/2, 
                      disk_ifu$observation$aperture_size/2, 
                      90,  
                      border = "black", density = NULL, lwd=3)


legend("bottomright", inset = c(0.02, 0.01),
       legend=c(expression("R"["50"]), expression("R"["e"]), "Field-of-view"),
       col=c("#882255","#DDCC77", "black"), lty=1, lwd=3, bty="n", cex = cex_val)


# Bulge mocks ------------------------------------------------------------------
par(new = T, mfrow = c(1, 1), xpd=T)
xl=-106.5; xr=-38
yb=65; yt=-40

# rect(-35, yb,
#      -1, yt,
#      63.5, yb,
#      col=rgb(0.01,0.01,0.01, alpha = 0.05), border=NULL)

lines(c(xl,xl), c(yb,yt), col="grey", lwd=3)
lines(c(xr,xr), c(yb,yt), col="grey", lwd=3)
lines(c(xl,xr), c(yb,yb), col="grey", lwd=3)
lines(c(xl,xr), c(yt,yt), col="grey", lwd=3)
lines(c(xr, -35), c(yt, -1), col="grey", lwd=3)
lines(c(xr, -35), c(yb, 64), col="grey", lwd=3)

radii = list("a" = bulge_photometry$major_axis_a_kpc/bulge_ifu$observation$sbin_size,
             "b" = bulge_photometry$minor_axis_b_kpc/bulge_ifu$observation$sbin_size,
             "ang" = bulge_photometry$primary_axis+90)

tenpc = 10/bulge_ifu$observation$sbin_size

mask = array(data = bulge_ifu$observation$aperture_region, 
             dim = c(bulge_ifu$observation$sbin,bulge_ifu$observation$sbin))
mask[mask == 0] = NA

plot_particles(bulge_ifu$raw_images$particle_image*mask, 
               new = T, 
               fig=c(ocol1,orow3), 
               labN = 1, titleshift = -7,
               units = "Particles", 
               radii=radii, radii_col = "#DDCC77", cex = cex_val)
lines(c((15-tenpc), 15), c(5,5), lwd=3, col="black")
text((15-(tenpc/2)), 8, "10 kpc", font=2, cex = cex_val)

plot_flux(bulge_ifu$observed_images$flux_image*mask, 
          new = T, 
          fig=c(ocol2,orow3), 
          labN = 1, titleshift = -7.5, 
          radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_velocity(bulge_ifu$observed_images$velocity_image*mask, 
              new = T, 
              fig=c(ocol1,orow2), 
              labN = 2, titleshift = -7.5,
              units = expression("v"["LOS"]*", km/s"), 
              radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_dispersion(bulge_ifu$observed_images$dispersion_image*mask, 
                new = T, 
                fig=c(ocol2,orow2), 
                labN = 2, titleshift = -7.5,
                units = expression('\u03C3'["LOS"]*", km/s"), 
                radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_h3(bulge_ifu$observed_images$h3_image*mask, 
        new = T, 
        fig=c(ocol1,orow1), 
        zlim = c(-0.1,0.1),
        labN = 2, titleshift = -7,
        units = expression("h"["3"]), 
        radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_h4(bulge_ifu$observed_images$h4_image*mask, 
        new = T, 
        fig=c(ocol2,orow1),
        zlim = c(-0.1,0.1),
        labN = 2, titleshift = -7,
        units = expression("h"["4"]), 
        radii=radii, radii_col = "#DDCC77", cex = cex_val)

# Disk mocks ------------------------------------------------------------------
par(new = T, mfrow = c(1, 1), xpd=T)
xl=21; xr=40.5
yb=-27; yt=33

# rect(-35, yb,
#      -1, yt,
#      63.5, yb,
#      col=rgb(0.01,0.01,0.01, alpha = 0.05), border=NULL)

lines(c(xl,xl), c(yb,yt), col="grey", lwd=3)
lines(c(xr,xr), c(yb,yt), col="grey", lwd=3)
lines(c(xl,xr), c(yb,yb), col="grey", lwd=3)
lines(c(xl,xr), c(yt,yt), col="grey", lwd=3)
lines(c(xl, 20.1), c(yt, 13), col="grey", lwd=3)
lines(c(xl, 20.1), c(yb, -23.5), col="grey", lwd=3)

radii = list("a" = disk_photometry$major_axis_a_kpc/disk_ifu$observation$sbin_size,
             "b" = disk_photometry$minor_axis_b_kpc/disk_ifu$observation$sbin_size,
             "ang" = disk_photometry$primary_axis+90)

mask = array(data = disk_ifu$observation$aperture_region, 
             dim = c(disk_ifu$observation$sbin,disk_ifu$observation$sbin))

mask[mask == 0] = NA

tenpc = 10/disk_ifu$observation$sbin_size

plot_particles(disk_ifu$raw_images$particle_image*mask, 
               new = T, 
               fig=c(ocol3,orow6), 
               labN = 1, titleshift = -7,
               units = "Particles", 
               radii=radii, radii_col = "#DDCC77", cex = cex_val)
lines(c(5, (5+tenpc)), c(5,5), lwd=3, col="black")
text((5+(tenpc/2)), 8, "10 kpc", font=2, cex = cex_val)

plot_flux(disk_ifu$observed_images$flux_image*mask, 
          new = T, 
          fig=c(ocol4,orow6), 
          labN = 1, titleshift = -7.5, 
          radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_velocity(disk_ifu$observed_images$velocity_image*mask, 
              new = T, 
              fig=c(ocol3,orow5), 
              labN = 2, titleshift = -7.5,
              units = expression("v"["LOS"]*", km/s"), 
              radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_dispersion(disk_ifu$observed_images$dispersion_image*mask, 
                new = T, 
                fig=c(ocol4,orow5), 
                labN = 2, titleshift = -7.5,
                units = expression('\u03C3'["LOS"]*", km/s"), 
                radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_h3(disk_ifu$observed_images$h3_image*mask, 
        new = T, 
        fig=c(ocol3,orow4), 
        zlim = c(-0.15,0.15),
        labN = 2, titleshift = -7,
        units = expression("h"["3"]), 
        radii=radii, radii_col = "#DDCC77", cex = cex_val)

plot_h4(disk_ifu$observed_images$h4_image*mask, 
        new = T, 
        fig=c(ocol4,orow4),
        zlim = c(-0.15,0.15),
        labN = 2, titleshift = -7,
        units = expression("h"["4"]), 
        radii=radii, radii_col = "#DDCC77", cex = cex_val)

dev.off()




