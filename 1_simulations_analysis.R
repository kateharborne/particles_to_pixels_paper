# Kate Harborne 17/10/23
# ============================================================================ #
# Analysis of simulations: 
# First, making the necessary SimSpin files and running sim_analysis to check 
# that the cumulative mass distribution between the models is similar. 
# Generating plots of the cumulative mass wrt radius for each flavour of model.
# This file constructs Figures 1 & 2 for the manuscript.
# ============================================================================ #

library(SimSpin)
library(magicaxis)
source("functions/plot_array.R")

# PATHS ========================================================================
loc = "simulations/hdf5/" # stored using Git LFS
ss_loc = "simspin_files/" 
sa_loc = "sim_analysis/"
plot_loc = "plots/"

# Making some simspin files (DISK) =============================================

disk_control_galaxy = make_simspin_file(filename = paste0(loc, "res_disk_sm650_control.hdf5"),
                                        write_to_file = T, overwrite = T,
                                        cores = 3, centre = c(0,0,0),
                                        output = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"))

disk_650_p0_galaxy =  make_simspin_file(filename = paste0(loc, "res_disk_sm650_p0.hdf5"), 
                                        write_to_file = T, 
                                        cores = 3, centre = c(0,0,0),
                                        output = paste0(ss_loc, "disk_650_p0_simspin_file.Rdata"))

disk_650_p16_galaxy =  make_simspin_file(filename = paste0(loc, "res_disk_sm650_p16.hdf5"), 
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "disk_650_p16_simspin_file.Rdata"))

disk_650_p50_galaxy =  make_simspin_file(filename = paste0(loc, "res_disk_sm650_p50.hdf5"), 
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "disk_650_p50_simspin_file.Rdata"))

disk_650_p84_galaxy =  make_simspin_file(filename = paste0(loc, "res_disk_sm650_p84.hdf5"), 
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "disk_650_p84_simspin_file.Rdata"))

disk_650_p100_galaxy = make_simspin_file(filename = paste0(loc, "res_disk_sm650_p100.hdf5"), 
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "disk_650_p100_simspin_file.Rdata"))

# Making some simspin files (BULGE) ============================================

bulge_control_galaxy = make_simspin_file(filename = paste0(loc, "res_bulge_sm650_control.hdf5"),
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"))

bulge_650_p0_galaxy =  make_simspin_file(filename = paste0(loc, "res_bulge_sm650_p0.hdf5"),
                                         write_to_file = T, 
                                         cores = 3, centre = c(0,0,0),
                                         output = paste0(ss_loc, "bulge_650_p0_simspin_file.Rdata"))

bulge_650_p16_galaxy =  make_simspin_file(filename = paste0(loc, "res_bulge_sm650_p16.hdf5"), 
                                          write_to_file = T, 
                                          cores = 3, centre = c(0,0,0),
                                          output = paste0(ss_loc, "bulge_650_p16_simspin_file.Rdata"))

bulge_650_p50_galaxy =  make_simspin_file(filename =  paste0(loc, "res_bulge_sm650_p50.hdf5"), 
                                          write_to_file = T, 
                                          cores = 3, centre = c(0,0,0),
                                          output = paste0(ss_loc, "bulge_650_p50_simspin_file.Rdata"))

bulge_650_p84_galaxy =  make_simspin_file(filename =  paste0(loc, "res_bulge_sm650_p84.hdf5"), 
                                          write_to_file = T, 
                                          cores = 3, centre = c(0,0,0),
                                          output = paste0(ss_loc, "bulge_650_p84_simspin_file.Rdata"))

bulge_650_p100_galaxy = make_simspin_file(filename =  paste0(loc, "res_bulge_sm650_p100.hdf5"), 
                                          write_to_file = T, 
                                          cores = 3, centre = c(0,0,0),
                                          output = paste0(ss_loc, "bulge_650_p100_simspin_file.Rdata"))

# Running sim_analysis (DISK) ==================================================

disk_bins = c(seq(0,40, by=0.5))

p0_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_p0_simspin_file.Rdata"),
                       type="stars", bin_breaks = disk_bins, verbose = F, shape = F)
p16_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_p16_simspin_file.Rdata"), 
                        type="stars", bin_breaks = disk_bins, verbose = F, shape = F)
p50_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_p50_simspin_file.Rdata"), 
                        type="stars", bin_breaks = disk_bins, verbose = F, shape = F)
p84_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_p84_simspin_file.Rdata"), 
                        type="stars", bin_breaks = disk_bins, verbose = F, shape = F)
p100_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_p100_simspin_file.Rdata"), 
                         type="stars", bin_breaks = disk_bins, verbose = F, shape = F)
control_disk = sim_analysis(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                            type="stars", bin_breaks = disk_bins, verbose = F, shape = F)

saveRDS(p0_disk, paste0(sa_loc,"p0_disk.Rdata"))
saveRDS(p16_disk, paste0(sa_loc,"p16_disk.Rdata"))
saveRDS(p50_disk, paste0(sa_loc,"p50_disk.Rdata"))
saveRDS(p84_disk, paste0(sa_loc,"p84_disk.Rdata"))
saveRDS(p100_disk, paste0(sa_loc,"p100_disk.Rdata"))
saveRDS(control_disk, paste0(sa_loc,"control_disk.Rdata"))

# Running sim_analysis (BULGE) =================================================
bulge_bins = c(seq(0,200, by=2))

p0_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_p0_simspin_file.Rdata"), 
                        type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)
p16_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_p16_simspin_file.Rdata"), 
                         type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)
p50_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_p50_simspin_file.Rdata"), 
                         type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)
p84_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_p84_simspin_file.Rdata"), 
                         type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)
p100_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_p100_simspin_file.Rdata"), 
                          type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)
control_bulge = sim_analysis(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"), 
                             type="stars", bin_breaks = bulge_bins, verbose = F, shape = F)

saveRDS(p0_bulge, paste0(sa_loc,"p0_bulge.Rdata"))
saveRDS(p16_bulge, paste0(sa_loc,"p16_bulge.Rdata"))
saveRDS(p50_bulge, paste0(sa_loc,"p50_bulge.Rdata"))
saveRDS(p84_bulge, paste0(sa_loc,"p84_bulge.Rdata"))
saveRDS(p100_bulge, paste0(sa_loc,"p100_bulge.Rdata"))
saveRDS(control_bulge, paste0(sa_loc,"control_bulge.Rdata"))

# Loading data =================================================================
p0_disk = readRDS(paste0(sa_loc,"p0_disk.Rdata"))
p16_disk = readRDS(paste0(sa_loc,"p16_disk.Rdata"))
p50_disk = readRDS(paste0(sa_loc,"p50_disk.Rdata"))
p84_disk = readRDS(paste0(sa_loc,"p84_disk.Rdata"))
p100_disk = readRDS(paste0(sa_loc,"p100_disk.Rdata"))
control_disk = readRDS(paste0(sa_loc,"control_disk.Rdata"))

p0_bulge = readRDS(paste0(sa_loc,"p0_bulge.Rdata"))
p16_bulge = readRDS(paste0(sa_loc,"p16_bulge.Rdata"))
p50_bulge = readRDS(paste0(sa_loc,"p50_bulge.Rdata"))
p84_bulge = readRDS(paste0(sa_loc,"p84_bulge.Rdata"))
p100_bulge = readRDS(paste0(sa_loc,"p100_bulge.Rdata"))
control_bulge = readRDS(paste0(sa_loc,"control_bulge.Rdata"))

# Plotting cumulative mass distributions =======================================

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")

w = plot_array(0.38, 4)

# Bulge plot -------------------------------------------------------------------
pdf.options(fonts="Arial") # Figure 1
cairo_pdf(file=paste0(plot_loc, "bulge_cumulative_mass_updated.pdf"), 
          width = 5.5, height = 10)

# Cumulative mass section ------------------------------------------------------
par(fig=c(0,1,w[4,1],w[4,2]), mar=c(5.1,4.1,4.1,2.1), family="sans")
magplot(control_bulge$RadialTrends_Cylindrical$Radius, 
        control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
        type="l", lwd = 2, col=tol_pal[1], 
        xlim = c(0, 100), xaxt="n", ylim = c(0,6.2),
        ylab = bquote("Cumulative Stellar Mass," ~ 10^10 ~ M[.('\u0298')]), )
lines(p100_bulge$RadialTrends_Cylindrical$Radius,
      p100_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[2])
lines(p84_bulge$RadialTrends_Cylindrical$Radius,
      p84_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[3])
lines(p50_bulge$RadialTrends_Cylindrical$Radius,
      p50_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[4])
lines(p16_bulge$RadialTrends_Cylindrical$Radius,
      p16_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[5])
lines(p0_bulge$RadialTrends_Cylindrical$Radius,
      p0_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[6])
abline(v=control_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)
text(85, 5, "Bulge Models", font = 2, cex=1.1)

par(fig=c(0,1,w[3,1],w[3,2]), mar=c(5.1,4.1,4.1,2.1), family="sans", new=T)
magplot(control_bulge$RadialTrends_Cylindrical$Radius,
        (control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
           control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
        type="l", lwd = 2, col=tol_pal[1], xlim = c(0, 100),
        ylim = c(-1.8,1.8), ylab="Residual, %", crunch=T, #xlab="Radius, kpc")
        xaxt="n")
lines(p100_bulge$RadialTrends_Cylindrical$Radius,
      (p100_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[2])
lines(p84_bulge$RadialTrends_Cylindrical$Radius,
      (p84_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[3])
lines(p50_bulge$RadialTrends_Cylindrical$Radius,
      (p50_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[4])
lines(p16_bulge$RadialTrends_Cylindrical$Radius,
      (p16_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[5])
lines(p0_bulge$RadialTrends_Cylindrical$Radius,
      (p0_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_bulge$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[6])
abline(v=control_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

par(fig=c(0,1,w[2,1],w[2,2]), mar=c(5.1,4.1,4.1,2.1), family="sans", new=TRUE)
magplot(control_bulge$RadialTrends_Cylindrical$Radius, 
        control_bulge$RadialTrends_Cylindrical$RotationalVelocity,
        type="l", lwd = 2, col=tol_pal[1], 
        xlim = c(0, 100), ylim = c(-25,25),
        xaxt="n",
        ylab = expression("V"['\u03C6']*", km/s"))
lines(p0_bulge$RadialTrends_Cylindrical$Radius,
      p0_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[6])
lines(p16_bulge$RadialTrends_Cylindrical$Radius,
      p16_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[5])
lines(p50_bulge$RadialTrends_Cylindrical$Radius,
      p50_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[4])
lines(p84_bulge$RadialTrends_Cylindrical$Radius,
      p84_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[3])
lines(p100_bulge$RadialTrends_Cylindrical$Radius,
      p100_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[2])
lines(control_bulge$RadialTrends_Cylindrical$Radius,
      control_bulge$RadialTrends_Cylindrical$RotationalVelocity,
      lwd=2, col=tol_pal[1])
abline(v=control_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

# Sigma_rot section ------------------------------------------------------------

par(fig = c(0,1,w[1,1],w[1,2]), new=TRUE, family="sans")
magplot(control_bulge$RadialTrends_Cylindrical$Radius, 
        control_bulge$RadialTrends_Cylindrical$RotationalDispersion,
        type="l", lwd = 2, col=tol_pal[1],
        xlim = c(0,100), xlab = "Radius, kpc",
        ylim = c(95,195), ylab = expression('\u03C3'['\u03C6']*", km/s"))
lines(p0_bulge$RadialTrends_Cylindrical$Radius,
      p0_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[6])
lines(p16_bulge$RadialTrends_Cylindrical$Radius,
      p16_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[5])
lines(p50_bulge$RadialTrends_Cylindrical$Radius,
      p50_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[4])
lines(p84_bulge$RadialTrends_Cylindrical$Radius,
      p84_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[3])
lines(p100_bulge$RadialTrends_Cylindrical$Radius,
      p100_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[2])
lines(control_bulge$RadialTrends_Cylindrical$Radius,
      control_bulge$RadialTrends_Cylindrical$RotationalDispersion,
      lwd=2, col=tol_pal[1])
abline(v=control_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_bulge$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

legend("topright", inset = c(0.05, 0),
       legend=c("control", "p100", "p84", "p50", "p16", "p0"),
       col=c(tol_pal[c(1,2,3,4,5,6)]), pch=15, bty="n",
       cex = 1)

# off ==========================================================================
dev.off()


# Disk plot --------------------------------------------------------------------
pdf.options(fonts="Arial") # Figure 2
cairo_pdf(file=paste0(plot_loc, "disk_cumulative_mass_updated.pdf"), 
          width = 5.5, height = 10)

# Cumulative Mass Section ------------------------------------------------------
par(fig=c(0,1,w[4,1], w[4,2]), mar=c(5.1,4.1,4.1,2.1), family="sans")
magplot(control_disk$RadialTrends_Cylindrical$Radius, 
        control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
        type="l", lwd = 2, col=tol_pal[1], xlim = c(0, 40),
        xaxt="n", 
        ylab = expression("Cumulative Stellar Mass, 10"^{10}*"M"['\u0298']))
lines(p100_disk$RadialTrends_Cylindrical$Radius,
      p100_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[2])
lines(p84_disk$RadialTrends_Cylindrical$Radius,
      p84_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[3])
lines(p50_disk$RadialTrends_Cylindrical$Radius,
      p50_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[4])
lines(p16_disk$RadialTrends_Cylindrical$Radius,
      p16_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[5])
lines(p0_disk$RadialTrends_Cylindrical$Radius,
      p0_disk$RadialTrends_Cylindrical$CumulativeMass/1e10,
      lwd=2, col=tol_pal[6])
abline(v=control_disk$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_disk$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_disk$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_disk$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_disk$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_disk$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)
text(35, 5.7, "Disc Models", font = 2, cex=1.1)

par(fig=c(0,1,w[3,1], w[3,2]), mar=c(5.1,4.1,4.1,2.1), family="sans", new=T)
magplot(control_disk$RadialTrends_Cylindrical$Radius,
        (control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
           control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
        type="l", lwd = 2, col=tol_pal[1], xlim = c(0, 40),
        ylim = c(-1.7,1.7), ylab="Residual, %", #xlab="Radius, kpc", 
        xaxt="n", crunch=T)
lines(p100_disk$RadialTrends_Cylindrical$Radius,
      (p100_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[2])
lines(p84_disk$RadialTrends_Cylindrical$Radius,
      (p84_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[3])
lines(p50_disk$RadialTrends_Cylindrical$Radius,
      (p50_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[4])
lines(p16_disk$RadialTrends_Cylindrical$Radius,
      (p16_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[5])
lines(p0_disk$RadialTrends_Cylindrical$Radius,
      (p0_disk$RadialTrends_Cylindrical$CumulativeMass/1e10-
         control_disk$RadialTrends_Cylindrical$CumulativeMass/1e10)*100,
      lwd=2, col=tol_pal[6])
abline(v=control_disk$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_disk$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_disk$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_disk$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_disk$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_disk$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

# Vrot section -----------------------------------------------------------------

par(fig=c(0,1,w[2,1], w[2,2]), mar=c(5.1,4.1,4.1,2.1), family="sans", new=TRUE)
magplot(control_disk$RadialTrends_Cylindrical$Radius[control_disk$RadialTrends_Cylindrical$NumberOfParticles>0], 
        control_disk$RadialTrends_Cylindrical$RotationalVelocity[control_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
        type="l", lwd = 2, col=tol_pal[1], 
        xlim = c(0, 40), ylim = c(35,290),
        xaxt="n",
        ylab = expression("V"['\u03C6']*", km/s"))
lines(p100_disk$RadialTrends_Cylindrical$Radius[p100_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p100_disk$RadialTrends_Cylindrical$RotationalVelocity[p100_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[2])
lines(p84_disk$RadialTrends_Cylindrical$Radius[p84_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p84_disk$RadialTrends_Cylindrical$RotationalVelocity[p84_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[3])
lines(p50_disk$RadialTrends_Cylindrical$Radius[p50_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p50_disk$RadialTrends_Cylindrical$RotationalVelocity[p50_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[4])
lines(p16_disk$RadialTrends_Cylindrical$Radius[p16_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p16_disk$RadialTrends_Cylindrical$RotationalVelocity[p16_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[5])
lines(p0_disk$RadialTrends_Cylindrical$Radius[p0_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p0_disk$RadialTrends_Cylindrical$RotationalVelocity[p0_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[6])
abline(v=control_disk$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_disk$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_disk$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_disk$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_disk$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_disk$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

# Sigma_rot section ------------------------------------------------------------

par(fig = c(0,1,w[1,1], w[1,2]), new=TRUE, family="sans")
magplot(control_disk$RadialTrends_Cylindrical$Radius[control_disk$RadialTrends_Cylindrical$NumberOfParticles>0], 
        control_disk$RadialTrends_Cylindrical$RotationalDispersion[control_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
        type="l", lwd = 2, col=tol_pal[1],
        xlim = c(0,40), xlab = "Radius, kpc",
        ylim = c(-5,100), ylab = expression('\u03C3'['\u03C6']*", km/s"))
lines(p100_disk$RadialTrends_Cylindrical$Radius[p100_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p100_disk$RadialTrends_Cylindrical$RotationalDispersion[p100_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[2])
lines(p84_disk$RadialTrends_Cylindrical$Radius[p84_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p84_disk$RadialTrends_Cylindrical$RotationalDispersion[p84_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[3])
lines(p50_disk$RadialTrends_Cylindrical$Radius[p50_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p50_disk$RadialTrends_Cylindrical$RotationalDispersion[p50_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[4])
lines(p16_disk$RadialTrends_Cylindrical$Radius[p16_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p16_disk$RadialTrends_Cylindrical$RotationalDispersion[p16_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[5])
lines(p0_disk$RadialTrends_Cylindrical$Radius[p0_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      p0_disk$RadialTrends_Cylindrical$RotationalDispersion[p0_disk$RadialTrends_Cylindrical$NumberOfParticles>0],
      lwd=2, col=tol_pal[6])
abline(v=control_disk$HalfMassProperties$RadiusCircular, col = tol_pal[1], lty=2)
abline(v=p100_disk$HalfMassProperties$RadiusCircular, col = tol_pal[2], lty=2)
abline(v=p84_disk$HalfMassProperties$RadiusCircular, col = tol_pal[3], lty=2)
abline(v=p50_disk$HalfMassProperties$RadiusCircular, col = tol_pal[4], lty=2)
abline(v=p16_disk$HalfMassProperties$RadiusCircular, col = tol_pal[5], lty=2)
abline(v=p0_disk$HalfMassProperties$RadiusCircular, col = tol_pal[6], lty=2)

legend("topright", inset = c(0.05, 0.02),
       legend=c("control", "p100", "p84", "p50", "p16", "p0"),
       col=c(tol_pal[c(1,2,3,4,5,6)]), pch=15, bty="n",
       cex = 1)

# off ==========================================================================
dev.off()


