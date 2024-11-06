# Kate Harborne 17/10/23
# ============================================================================ #
# Running analysis of different resolution models: 
# Using the range of different resolution models for both the disk and bulge, 
# by how much do the kinematics in a given pixel change?
# ============================================================================ #

library(SimSpin)
library(magicaxis)

# PATHS ========================================================================
loc = "simulations/"
ss_loc = "simspin_files/"
data_loc = "data/"
plot_loc = "plots/"
sa_loc = "sim_analysis/"
cube_loc = paste0(data_loc, "/cubes/")

ifu = telescope(type="SAMI")
sky_00 = observing_strategy(dist_kpc_per_arcsec = 0.75, 
                            inc_deg = 0, blur = F,
                            pointing_kpc = c(0,0))
sky_90 = observing_strategy(dist_kpc_per_arcsec = 0.5, 
                            inc_deg = 90, blur = F, 
                            pointing_kpc = c(0,0))

# Mock observations of BULGE ===================================================
bulge_control = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_control, paste0(cube_loc, "bulge_control_mom4.Rdata"))

plot_particles(bulge_control$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_control$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_control$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

bulge_control_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(bulge_control_2, paste0(cube_loc, "bulge_control_mom2.Rdata"))

plot_particles(bulge_control_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_control_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_control_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

bulge_0 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p0_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_0, paste0(cube_loc, "bulge_650_p0_mom4.Rdata"))

plot_particles(bulge_0$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_0$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_0$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

bulge_0_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p0_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(bulge_0_2, paste0(cube_loc, "bulge_650_p0_mom2.Rdata"))

plot_particles(bulge_0_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_0_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_0_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

bulge_16 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p16_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_16, paste0(cube_loc, "bulge_650_p16_mom4.Rdata"))

plot_particles(bulge_16$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_16$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_16$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

bulge_16_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p16_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(bulge_16_2, paste0(cube_loc, "bulge_650_p16_mom2.Rdata"))

plot_particles(bulge_16_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_16_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_16_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

bulge_50 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p50_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_50, paste0(cube_loc, "bulge_650_p50_mom4.Rdata"))

plot_particles(bulge_50$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_50$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_50$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)


bulge_50_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p50_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments=2)

saveRDS(bulge_50_2, paste0(cube_loc, "bulge_650_p50_mom2.Rdata"))

plot_particles(bulge_50_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_50_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_50_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

bulge_84 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p84_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_84, paste0(cube_loc, "bulge_650_p84_mom4.Rdata"))

plot_particles(bulge_84$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_84$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_84$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

bulge_84_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p84_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments=2)

saveRDS(bulge_84_2, paste0(cube_loc, "bulge_650_p84_mom2.Rdata"))

plot_particles(bulge_84_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_84_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_84_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

# ---

bulge_100 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p100_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T)

saveRDS(bulge_100, paste0(cube_loc, "bulge_650_p100_mom4.Rdata"))

plot_particles(bulge_100$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_100$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_100$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

bulge_100_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p100_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(bulge_100_2, paste0(cube_loc, "bulge_650_p100_mom2.Rdata"))

plot_particles(bulge_100_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(bulge_100_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(bulge_100_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)



# Mock observations of DISK ====================================================
disk_control = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_control, paste0(cube_loc, "disk_control_mom4.Rdata"))

plot_particles(disk_control$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_control$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_control$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_control_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(disk_control_2, paste0(cube_loc, "disk_control_mom2.Rdata"))

plot_particles(disk_control_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_control_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_control_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

# ---

disk_0 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p0_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_0, paste0(cube_loc, "disk_650_p0_mom4.Rdata"))

plot_particles(disk_0$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_0$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_0$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_0_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p0_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments=2)

saveRDS(disk_0_2, paste0(cube_loc, "disk_650_p0_mom2.Rdata"))

plot_particles(disk_0_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_0_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_0_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

# ---

disk_16 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p16_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_16, paste0(cube_loc, "disk_650_p16_mom4.Rdata"))

plot_particles(disk_16$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_16$raw_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_16$raw_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_16_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p16_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(disk_16_2, paste0(cube_loc, "disk_650_p16_mom2.Rdata"))

plot_particles(disk_16_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_16_2$raw_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_16_2$raw_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

disk_50 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p50_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_50, paste0(cube_loc, "disk_650_p50_mom4.Rdata"))

plot_particles(disk_50$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_50$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_50$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_50_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p50_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(disk_50_2, paste0(cube_loc, "disk_650_p50_mom2.Rdata"))

plot_particles(disk_50_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_50_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_50_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

disk_84 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p84_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_84, paste0(cube_loc, "disk_650_p84_mom4.Rdata"))

plot_particles(disk_84$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_84$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_84$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_84_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p84_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(disk_84_2, paste0(cube_loc, "disk_650_p84_mom2.Rdata"))

plot_particles(disk_84_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_84_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_84_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

#---

disk_100 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p100_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T)

saveRDS(disk_100, paste0(cube_loc, "disk_650_p100_mom4.Rdata"))

plot_particles(disk_100$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_100$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_100$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)

disk_100_2 = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p100_simspin_file.Rdata"),
                 telescope = ifu, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, moments = 2)

saveRDS(disk_100_2, paste0(cube_loc, "disk_650_p100_mom2.Rdata"))

plot_particles(disk_100_2$raw_images$particle_image, fig = c(0,0.33,0,1))
plot_velocity(disk_100_2$observed_images$velocity_image, fig = c(0.33,0.67,0,1), new=T)
plot_dispersion(disk_100_2$observed_images$dispersion_image, fig = c(0.67,1,0,1), new = T)


# Loading saved data - moments 4 ===============================================
bulge_control = readRDS(paste0(cube_loc, "bulge_control_mom4.Rdata"))
bulge_0 = readRDS(paste0(cube_loc, "bulge_650_p0_mom4.Rdata"))
bulge_16 = readRDS(paste0(cube_loc, "bulge_650_p16_mom4.Rdata"))
bulge_50 = readRDS(paste0(cube_loc, "bulge_650_p50_mom4.Rdata"))
bulge_84 = readRDS(paste0(cube_loc, "bulge_650_p84_mom4.Rdata"))
bulge_100 = readRDS(paste0(cube_loc, "bulge_650_p100_mom4.Rdata"))

disk_control = readRDS(paste0(cube_loc, "disk_control_mom4.Rdata"))
disk_0 = readRDS(paste0(cube_loc, "disk_650_p0_mom4.Rdata"))
disk_16 = readRDS(paste0(cube_loc, "disk_650_p16_mom4.Rdata"))
disk_50 = readRDS(paste0(cube_loc, "disk_650_p50_mom4.Rdata"))
disk_84 = readRDS(paste0(cube_loc, "disk_650_p84_mom4.Rdata"))
disk_100 = readRDS(paste0(cube_loc, "disk_650_p100_mom4.Rdata"))

# Masking regions in BULGE =====================================================

inner_b_mask = bulge_0$raw_images$particle_image
inner_b_mask[inner_b_mask < 100] = 0
inner_b_mask[inner_b_mask!=0] = 1 
outer_b_mask = bulge_0$raw_images$particle_image
outer_b_mask[outer_b_mask < 50 | outer_b_mask >= 100] = 0
outer_b_mask[outer_b_mask!=0] = 1 

plot_particles(bulge_100$raw_images$particle_image*inner_b_mask)
plot_particles(bulge_100$raw_images$particle_image*outer_b_mask)

# Compute difference between full image and sampled images - BULGE =============

bulge_res_data = 
  list("InnerN" = 
         list(
           "p0" = quantile((bulge_0$raw_images$particle_image)[inner_b_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((bulge_16$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((bulge_50$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((bulge_84$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((bulge_100$raw_images$particle_image)[inner_b_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((bulge_control$raw_images$particle_image)[inner_b_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "InnerVelocity" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$velocity_image -
                               bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$velocity_image-
                                 bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerDispersion" =
         list(
           "p0" = quantile(((bulge_0$observed_images$dispersion_image -
                               bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$dispersion_image-
                                 bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerH3" =
         list(
           "p0" = quantile(((bulge_0$observed_images$h3_image -
                               bulge_control$observed_images$h3_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$h3_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$h3_image-
                                 bulge_control$observed_images$h3_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$h3_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerH4" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$h4_image -
                               bulge_control$observed_images$h4_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$h4_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$h4_image-
                                 bulge_control$observed_images$h4_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$h4_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       "OuterN" = 
         list(
           "p0" = quantile((bulge_0$raw_images$particle_image)[outer_b_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((bulge_16$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((bulge_50$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((bulge_84$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((bulge_100$raw_images$particle_image)[outer_b_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((bulge_control$raw_images$particle_image)[outer_b_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "OuterVelocity" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$velocity_image -
                               bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$velocity_image-
                                 bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterDispersion" =
         list(
           "p0" = quantile(((bulge_0$observed_images$dispersion_image -
                               bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$dispersion_image-
                                 bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterH3" =
         list(
           "p0" = quantile(((bulge_0$observed_images$h3_image -
                               bulge_control$observed_images$h3_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$h3_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$h3_image-
                                bulge_control$observed_images$h3_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h3_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$h3_image-
                                 bulge_control$observed_images$h3_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$h3_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterH4" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$h4_image -
                               bulge_control$observed_images$h4_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$h4_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$h4_image-
                                bulge_control$observed_images$h4_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$h4_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$h4_image-
                                 bulge_control$observed_images$h4_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$h4_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         )
  )


saveRDS(bulge_res_data, paste0(data_loc, "bulge_650_res_data_mom4.Rdata"))

# Masking regions in the DISK ==================================================

centre_d_mask = disk_0$raw_images$particle_image
centre_d_mask[11:20,] = 0; centre_d_mask[,14:17] = 0; centre_d_mask[centre_d_mask!=0] = 1
inner_d_mask = disk_0$raw_images$particle_image*centre_d_mask
inner_d_mask[inner_d_mask < 100] = 0; inner_d_mask[inner_d_mask!=0] = 1 
outer_d_mask = disk_0$raw_images$particle_image*centre_d_mask
outer_d_mask[outer_d_mask < 50 | outer_d_mask >= 100] = 0; outer_d_mask[outer_d_mask!=0] = 1 

plot_particles(disk_100$raw_images$particle_image*inner_d_mask)
plot_particles(disk_100$raw_images$particle_image*outer_d_mask)

# Compute difference between full image and sampled images - DISK ==============

disk_res_data = 
  list("InnerN" = 
         list(
           "p0" = quantile((disk_0$raw_images$particle_image)[inner_d_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((disk_16$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((disk_50$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((disk_84$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((disk_100$raw_images$particle_image)[inner_d_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((disk_control$raw_images$particle_image)[inner_d_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "InnerVelocity" = 
         list(
           "p0" = quantile(((disk_0$observed_images$velocity_image -
                               disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                              (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$velocity_image-
                                 disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                                (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerDispersion" =
         list(
           "p0" = quantile(((disk_0$observed_images$dispersion_image -
                               disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                              (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$dispersion_image-
                                 disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                                (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerH3" =
         list(
           "p0" = quantile(((disk_0$observed_images$h3_image -
                               disk_control$observed_images$h3_image)[inner_d_mask>0]/
                              (disk_control$observed_images$h3_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h3_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h3_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h3_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$h3_image-
                                 disk_control$observed_images$h3_image)[inner_d_mask>0]/
                                (disk_control$observed_images$h3_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerH4" = 
         list(
           "p0" = quantile(((disk_0$observed_images$h4_image -
                               disk_control$observed_images$h4_image)[inner_d_mask>0]/
                              (disk_control$observed_images$h4_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h4_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h4_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[inner_d_mask>0]/
                               (disk_control$observed_images$h4_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$h4_image-
                                 disk_control$observed_images$h4_image)[inner_d_mask>0]/
                                (disk_control$observed_images$h4_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       "OuterN" = 
         list(
           "p0" = quantile((disk_0$raw_images$particle_image)[outer_d_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((disk_16$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((disk_50$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((disk_84$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((disk_100$raw_images$particle_image)[outer_d_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((disk_control$raw_images$particle_image)[outer_d_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "OuterVelocity" = 
         list(
           "p0" = quantile(((disk_0$observed_images$velocity_image -
                               disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                              (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$velocity_image-
                                 disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                                (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterDispersion" =
         list(
           "p0" = quantile(((disk_0$observed_images$dispersion_image -
                               disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                              (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$dispersion_image-
                                 disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                                (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterH3" =
         list(
           "p0" = quantile(((disk_0$observed_images$h3_image -
                               disk_control$observed_images$h3_image)[outer_d_mask>0]/
                              (disk_control$observed_images$h3_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h3_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h3_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$h3_image-
                                disk_control$observed_images$h3_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h3_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$h3_image-
                                 disk_control$observed_images$h3_image)[outer_d_mask>0]/
                                (disk_control$observed_images$h3_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterH4" = 
         list(
           "p0" = quantile(((disk_0$observed_images$h4_image -
                               disk_control$observed_images$h4_image)[outer_d_mask>0]/
                              (disk_control$observed_images$h4_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h4_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h4_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$h4_image-
                                disk_control$observed_images$h4_image)[outer_d_mask>0]/
                               (disk_control$observed_images$h4_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$h4_image-
                                 disk_control$observed_images$h4_image)[outer_d_mask>0]/
                                (disk_control$observed_images$h4_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         )
  )

saveRDS(disk_res_data, paste0(data_loc, "disk_650_res_data_mom4.Rdata"))

# Loading saved data - moments 2 ===============================================
bulge_control = readRDS(paste0(cube_loc, "bulge_control_mom2.Rdata"))
bulge_0 = readRDS(paste0(cube_loc, "bulge_650_p0_mom2.Rdata"))
bulge_16 = readRDS(paste0(cube_loc, "bulge_650_p16_mom2.Rdata"))
bulge_50 = readRDS(paste0(cube_loc, "bulge_650_p50_mom2.Rdata"))
bulge_84 = readRDS(paste0(cube_loc, "bulge_650_p84_mom2.Rdata"))
bulge_100 = readRDS(paste0(cube_loc, "bulge_650_p100_mom2.Rdata"))

disk_control = readRDS(paste0(cube_loc, "disk_control_mom2.Rdata"))
disk_0 = readRDS(paste0(cube_loc, "disk_650_p0_mom2.Rdata"))
disk_16 = readRDS(paste0(cube_loc, "disk_650_p16_mom2.Rdata"))
disk_50 = readRDS(paste0(cube_loc, "disk_650_p50_mom2.Rdata"))
disk_84 = readRDS(paste0(cube_loc, "disk_650_p84_mom2.Rdata"))
disk_100 = readRDS(paste0(cube_loc, "disk_650_p100_mom2.Rdata"))

# Masking regions in BULGE =====================================================

inner_b_mask = bulge_0$raw_images$particle_image
inner_b_mask[inner_b_mask < 100] = 0
inner_b_mask[inner_b_mask!=0] = 1 
outer_b_mask = bulge_0$raw_images$particle_image
outer_b_mask[outer_b_mask < 50 | outer_b_mask >= 100] = 0
outer_b_mask[outer_b_mask!=0] = 1 

plot_particles(bulge_100$raw_images$particle_image*inner_b_mask)
plot_particles(bulge_100$raw_images$particle_image*outer_b_mask)

# Compute difference between full image and sampled images - BULGE =============

bulge_res_data = 
  list("InnerN" = 
         list(
           "p0" = quantile((bulge_0$raw_images$particle_image)[inner_b_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((bulge_16$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((bulge_50$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((bulge_84$raw_images$particle_image)[inner_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((bulge_100$raw_images$particle_image)[inner_b_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((bulge_control$raw_images$particle_image)[inner_b_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "InnerVelocity" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$velocity_image -
                               bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$velocity_image-
                                 bulge_control$observed_images$velocity_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$velocity_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerDispersion" =
         list(
           "p0" = quantile(((bulge_0$observed_images$dispersion_image -
                               bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                              (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$dispersion_image-
                                 bulge_control$observed_images$dispersion_image)[inner_b_mask>0]/
                                (bulge_control$observed_images$dispersion_image)[inner_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterN" = 
         list(
           "p0" = quantile((bulge_0$raw_images$particle_image)[outer_b_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((bulge_16$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((bulge_50$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((bulge_84$raw_images$particle_image)[outer_b_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((bulge_100$raw_images$particle_image)[outer_b_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((bulge_control$raw_images$particle_image)[outer_b_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "OuterVelocity" = 
         list(
           "p0" = quantile(((bulge_0$observed_images$velocity_image -
                               bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$velocity_image-
                                bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$velocity_image-
                                 bulge_control$observed_images$velocity_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$velocity_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterDispersion" =
         list(
           "p0" = quantile(((bulge_0$observed_images$dispersion_image -
                               bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                              (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((bulge_16$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((bulge_50$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((bulge_84$observed_images$dispersion_image-
                                bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                               (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((bulge_100$observed_images$dispersion_image-
                                 bulge_control$observed_images$dispersion_image)[outer_b_mask>0]/
                                (bulge_control$observed_images$dispersion_image)[outer_b_mask>0]),
                             c(0.25,0.5,0.75))
         )
       
  )


saveRDS(bulge_res_data, paste0(data_loc, "bulge_650_res_data_mom2.Rdata"))

# Masking regions in the DISK ==================================================

centre_d_mask = disk_0$raw_images$particle_image
centre_d_mask[11:20,] = 0; centre_d_mask[,14:17] = 0; centre_d_mask[centre_d_mask!=0] = 1
inner_d_mask = disk_0$raw_images$particle_image*centre_d_mask
inner_d_mask[inner_d_mask < 100] = 0; inner_d_mask[inner_d_mask!=0] = 1 
outer_d_mask = disk_0$raw_images$particle_image*centre_d_mask
outer_d_mask[outer_d_mask < 50 | outer_d_mask >= 100] = 0; outer_d_mask[outer_d_mask!=0] = 1 

plot_particles(disk_100$raw_images$particle_image*inner_d_mask)
plot_particles(disk_100$raw_images$particle_image*outer_d_mask)

# Compute difference between full image and sampled images - DISK ==============

disk_res_data = 
  list("InnerN" = 
         list(
           "p0" = quantile((disk_0$raw_images$particle_image)[inner_d_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((disk_16$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((disk_50$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((disk_84$raw_images$particle_image)[inner_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((disk_100$raw_images$particle_image)[inner_d_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((disk_control$raw_images$particle_image)[inner_d_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "InnerVelocity" = 
         list(
           "p0" = quantile(((disk_0$observed_images$velocity_image -
                               disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                              (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$velocity_image-
                                 disk_control$observed_images$velocity_image)[inner_d_mask>0]/
                                (disk_control$observed_images$velocity_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "InnerDispersion" =
         list(
           "p0" = quantile(((disk_0$observed_images$dispersion_image -
                               disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                              (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$dispersion_image-
                                 disk_control$observed_images$dispersion_image)[inner_d_mask>0]/
                                (disk_control$observed_images$dispersion_image)[inner_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterN" = 
         list(
           "p0" = quantile((disk_0$raw_images$particle_image)[outer_d_mask > 0],
                           c(0.25, 0.5, 0.75)),
           "p16" = quantile((disk_16$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p50" = quantile((disk_50$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p84" = quantile((disk_84$raw_images$particle_image)[outer_d_mask > 0],
                            c(0.25, 0.5, 0.75)),
           "p100" = quantile((disk_100$raw_images$particle_image)[outer_d_mask > 0],
                             c(0.25, 0.5, 0.75)),
           "control" = quantile((disk_control$raw_images$particle_image)[outer_d_mask >0],
                                c(0.25, 0.5, 0.75)
           )
         ), 
       
       "OuterVelocity" = 
         list(
           "p0" = quantile(((disk_0$observed_images$velocity_image -
                               disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                              (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$velocity_image-
                                disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                               (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$velocity_image-
                                 disk_control$observed_images$velocity_image)[outer_d_mask>0]/
                                (disk_control$observed_images$velocity_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         ),
       
       "OuterDispersion" =
         list(
           "p0" = quantile(((disk_0$observed_images$dispersion_image -
                               disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                              (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                           c(0.25,0.5,0.75)),
           "p16" = quantile(((disk_16$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p50" = quantile(((disk_50$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p84" = quantile(((disk_84$observed_images$dispersion_image-
                                disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                               (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                            c(0.25,0.5,0.75)),
           "p100" = quantile(((disk_100$observed_images$dispersion_image-
                                 disk_control$observed_images$dispersion_image)[outer_d_mask>0]/
                                (disk_control$observed_images$dispersion_image)[outer_d_mask>0]),
                             c(0.25,0.5,0.75))
         )
  )

saveRDS(disk_res_data, paste0(data_loc, "disk_650_res_data_mom2.Rdata"))