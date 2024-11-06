# Kate Harborne 23/10/23
# ============================================================================ #
# Given the shot noise, how does the spin parameter and other measured 
# kinematics vary due to the effects of shot noise? And how do these then
# change if we using voronoi bins?
# Making figures 7 and 8
# ============================================================================ #

library(SimSpin)
library(magicaxis)

# PATHS ========================================================================
ss_loc = "simspin_files/"
plot_loc = "plots/"
data_loc = "data/"
cube_loc = paste0(data_loc, "/cubes/flux_obs/")

# FUNCTIONS ====================================================================
source("functions/get_photometric_properties.R")
source("functions/get_kinematic_properties.R")
source("functions/plot_array.R")
source("functions/shadowtext.R")

# INPUTS FOR SIMSPIN ===========================================================
photo = telescope(type="IFU", 
                  spatial_res = 0.2, 
                  fov = 500, 
                  aperture_shape = "square", 
                  wave_res = 200, 
                  filter = "r") # KIDS-like imaging, r-band
ifu = telescope(type="SAMI")

sky_00 = observing_strategy(dist_kpc_per_arcsec = 1.75, 
                            inc_deg = 0, blur = F,
                            pointing_kpc = c(0,0))
sky_90 = observing_strategy(dist_kpc_per_arcsec = 1, 
                            inc_deg = 90, blur = F, 
                            pointing_kpc = c(0,0))

# Mock photometery observations of BULGE =======================================
bulge_control_mass =
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"),
                 telescope = photo, observing_strategy = sky_00,
                 method = "velocity", mass_flag = T, cores = 3)
saveRDS(bulge_control_mass, file=paste0(data_loc,"bulge_sm650_control_photometery_mass_obs_500.Rdata"))

bulge_control_photo =
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"),
                 telescope = photo, observing_strategy = sky_00,
                 method = "velocity", mass_flag = F, cores = 3)
saveRDS(bulge_control_photo, file=paste0(data_loc,"bulge_sm650_control_photometery_obs_500.Rdata"))

plot_flux(bulge_control_photo$observed_images$flux_image)

# Measuring the half-light radius of this object 
bulge_control_mass = readRDS(paste0(data_loc,"bulge_sm650_control_photometery_mass_obs_500.Rdata"))
bulge_control_photo = readRDS(paste0(data_loc,"bulge_sm650_control_photometery_obs_500.Rdata"))

mass_frac = sum(bulge_control_mass$raw_images$mass_image)/(6.5e10)

bulge_photometry = get_photometry_properties(
  photo_flux_image = bulge_control_photo$observed_images$flux_image,
  obs_summary = bulge_control_photo$observation, 
  plot = T, 
  mass_frac = mass_frac)

saveRDS(bulge_photometry, paste0(data_loc, "bulge_sm650_control_photo_values_500.Rdata"))

# Considering the true kinematics within the measured radius -------------------
bulge_control_ifu =
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"),
                 telescope = ifu,
                 observing_strategy = sky_00,
                 method = "velocity",
                 mass_flag = F, moments = 2)
saveRDS(bulge_control_ifu, paste0(cube_loc, "bulge_sm650_control_ifu_obs.Rdata"))

bulge_control_ifu = readRDS(paste0(cube_loc, "bulge_sm650_control_ifu_obs.Rdata"))
bulge_mask = bulge_control_ifu$raw_images$particle_image
bulge_mask[bulge_mask < 200] = NA; bulge_mask[!is.na(bulge_mask)] = 1

plot_particles(bulge_control_ifu$raw_images$particle_image*bulge_mask)

bulge_control_kin = get_kinematic_properties(
  ifu_output = bulge_control_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=T
)

saveRDS(bulge_control_kin, paste0(data_loc, "/kinematics/bulge_sm650_control_kin_values.Rdata"))
ap_mask = array(data = bulge_control_ifu$observation$aperture_region, dim = c(30,30))
ap_mask[ap_mask==0]=NA
plot_flux(bulge_control_ifu$observed_images$flux_image*ap_mask, 
          radii = list(
            "a" = bulge_photometry$major_axis_a_kpc/
              bulge_control_ifu$observation$sbin_size,
            "b" = bulge_photometry$minor_axis_b_kpc/
              bulge_control_ifu$observation$sbin_size,
            "ang" = bulge_photometry$primary_axis))

# Make similar observations of the other resolution models ---------------------
bulge_p100_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p100_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_00,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(bulge_p100_ifu, paste0(cube_loc, "bulge_sm650_p100_ifu_obs.Rdata"))

bulge_p100_kin = get_kinematic_properties(
  ifu_output = bulge_p100_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(bulge_p100_kin, paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values.Rdata"))

bulge_p84_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p84_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_00,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(bulge_p84_ifu, paste0(cube_loc, "bulge_sm650_p84_ifu_obs.Rdata"))

bulge_p84_kin = get_kinematic_properties(
  ifu_output = bulge_p84_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(bulge_p84_kin, paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values.Rdata"))


bulge_p50_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p50_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_00,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(bulge_p50_ifu, paste0(cube_loc, "bulge_sm650_p50_ifu_obs.Rdata"))

bulge_p50_kin = get_kinematic_properties(
  ifu_output = bulge_p50_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(bulge_p50_kin, paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values.Rdata"))

bulge_p16_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p16_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_00,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(bulge_p16_ifu, paste0(cube_loc, "bulge_sm650_p16_ifu_obs.Rdata"))

bulge_p16_kin = get_kinematic_properties(
  ifu_output = bulge_p16_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(bulge_p16_kin, paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values.Rdata"))

bulge_p0_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p0_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_00,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(bulge_p0_ifu, paste0(cube_loc, "bulge_sm650_p0_ifu_obs.Rdata"))

bulge_p0_kin = get_kinematic_properties(
  ifu_output = bulge_p0_ifu,
  photometry_properties = bulge_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(bulge_p0_kin, paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values.Rdata"))

# Mock photometery observations of DISK ========================================
disk_control_mass =
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                 telescope = photo, observing_strategy = sky_90,
                 method = "velocity", mass_flag = T, cores = 3)
saveRDS(disk_control_mass, file=paste0(data_loc,"disk_sm650_control_photometery_mass_obs_500.Rdata"))

disk_control_photo =
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                 telescope = photo, observing_strategy = sky_90,
                 method = "velocity", mass_flag = F, cores=3)
saveRDS(disk_control_photo, file=paste0(data_loc,"disk_sm650_control_photometery_obs_500.Rdata"))

plot_flux(disk_control_photo$observed_images$flux_image)

# Measuring the half-light radius of this object 
disk_control_mass = readRDS(paste0(data_loc,"disk_sm650_control_photometery_mass_obs_500.Rdata"))
disk_control_photo = readRDS(paste0(data_loc,"disk_sm650_control_photometery_obs_500.Rdata"))

mass_frac = sum(disk_control_mass$raw_images$mass_image)/(6.5e10)

disk_photometry = get_photometry_properties(
  photo_flux_image = disk_control_photo$observed_images$flux_image,
  obs_summary = disk_control_photo$observation, 
  plot = T, 
  mass_frac = mass_frac)

saveRDS(disk_photometry, paste0(data_loc, "disk_sm650_control_photo_values_500.Rdata"))

# Considering the true kinematics within the measured radius -------------------
disk_control_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, cores = 3, moments = 2)
saveRDS(disk_control_ifu, paste0(cube_loc, "disk_sm650_control_ifu_obs.Rdata"))

disk_photometry = readRDS(paste0(data_loc, "disk_sm650_control_photo_values.Rdata"))
disk_mask = disk_control_ifu$raw_images$particle_image
disk_mask[disk_mask < 200] = NA; disk_mask[!is.na(disk_mask)] = 1

plot_particles(disk_control_ifu$raw_images$particle_image*disk_mask)

disk_control_kin = get_kinematic_properties(
  ifu_output = disk_control_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=T
)

saveRDS(disk_control_kin, paste0(data_loc, "/kinematics/disk_sm650_control_kin_values.Rdata"))

plot_flux(disk_control_ifu$observed_images$flux_image, 
          radii = list(
            "a" = disk_photometry$major_axis_a_kpc/
              disk_control_ifu$observation$sbin_size,
            "b" = disk_photometry$minor_axis_b_kpc/
              disk_control_ifu$observation$sbin_size,
            "ang" = disk_photometry$primary_axis+90))

# Make similar observations of the other resolution models ---------------------
disk_p100_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p100_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(disk_p100_ifu, paste0(cube_loc, "disk_sm650_p100_ifu_obs.Rdata"))

disk_p100_kin = get_kinematic_properties(
  ifu_output = disk_p100_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(disk_p100_kin, paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values.Rdata"))

disk_p84_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p84_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(disk_p84_ifu, paste0(cube_loc, "disk_sm650_p84_ifu_obs.Rdata"))

disk_p84_kin = get_kinematic_properties(
  ifu_output = disk_p84_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(disk_p84_kin, paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values.Rdata"))


disk_p50_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p50_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(disk_p50_ifu, paste0(cube_loc, "disk_sm650_p50_ifu_obs.Rdata"))

disk_p50_kin = get_kinematic_properties(
  ifu_output = disk_p50_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(disk_p50_kin, paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values.Rdata"))

disk_p16_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p16_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(disk_p16_ifu, paste0(cube_loc, "disk_sm650_p16_ifu_obs.Rdata"))

disk_p16_kin = get_kinematic_properties(
  ifu_output = disk_p16_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(disk_p16_kin, paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values.Rdata"))

disk_p0_ifu = 
  build_datacube(simspin_file = paste0(ss_loc, "disk_650_p0_simspin_file.Rdata"),
                 telescope = ifu, 
                 observing_strategy = sky_90,
                 method = "velocity", 
                 mass_flag = F, moments = 2)
saveRDS(disk_p0_ifu, paste0(cube_loc, "disk_sm650_p0_ifu_obs.Rdata"))

disk_p0_kin = get_kinematic_properties(
  ifu_output = disk_p0_ifu,
  photometry_properties = disk_photometry,
  particle_limit = 0, 
  plot=F
)
saveRDS(disk_p0_kin, paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values.Rdata"))

# Voronoi binning to meet a minimum number -------------------------------------
disk_photometry = readRDS(paste0(data_loc, "disk_sm650_control_photo_values_500.Rdata"))
bulge_photometry = readRDS(paste0(data_loc, "bulge_sm650_control_photo_values_500.Rdata"))
min_particle_number = seq(50,500,by=50)
part_res = c(0,16,50,84,100)

for(p in part_res){
  for (each in min_particle_number){
    
    assign(paste0("bulge_p", p,"_ifu_v", each),
           build_datacube(simspin_file = paste0(ss_loc, "bulge_650_p", p,"_simspin_file.Rdata"),
                          telescope = ifu,
                          observing_strategy = sky_00,
                          method = "velocity",
                          mass_flag = F, voronoi_bin = T,
                          vorbin_limit = as.numeric(each), moments = 2)
    )
    
    saveRDS(get(paste0("bulge_p", p,"_ifu_v", each)), 
            paste0(cube_loc, "bulge_sm650_p", p,"_ifu_obs_vorbin_", each,".Rdata"))
    
    assign(paste0("bulge_p", p,"_kin_v", each),
           get_kinematic_properties(
             ifu_output = get(paste0("bulge_p", p,"_ifu_v", each)),
             photometry_properties = bulge_photometry,
             particle_limit = 0, 
             plot=F, voronoi = T)
    )
    saveRDS(get(paste0("bulge_p", p,"_kin_v", each)), 
            paste0(data_loc, "/kinematics/bulge_sm650_p", p,"_kin_values_vorbin_", each,".Rdata"))
    
    cat("Done with bulge model p", p, " and vorbin ", each, ".\n")
    
    assign(paste0("disk_p", p,"_ifu_v", each),
           build_datacube(simspin_file = paste0(ss_loc, "disk_650_p", p,"_simspin_file.Rdata"),
                          telescope = ifu,
                          observing_strategy = sky_90,
                          method = "velocity",
                          mass_flag = F, voronoi_bin = T,
                          vorbin_limit = as.numeric(each), moments = 2)
    )
    
    saveRDS(get(paste0("disk_p", p,"_ifu_v", each)),
            paste0(cube_loc, "disk_sm650_p", p,"_ifu_obs_vorbin_", each,".Rdata"))
    
    assign(paste0("disk_p", p,"_kin_v", each),
           get_kinematic_properties(
             ifu_output = get(paste0("disk_p", p,"_ifu_v", each)),
             photometry_properties = disk_photometry,
             particle_limit = 0,
             plot=F, voronoi = T)
    )
    saveRDS(get(paste0("disk_p", p,"_kin_v", each)),
            paste0(data_loc, "/kinematics/disk_sm650_p", p,"_kin_values_vorbin_", each,".Rdata"))
    
    cat("Done with disk model p", p, " and vorbin ", each, ".\n")
    
  }
  
}

# Fix for vorbin unique pixel count --------------------------------------------

for(p in part_res){
  for (each in min_particle_number){
    
    ifu_output = get(paste0("disk_p", p, "_ifu_v", each))
    kin_output = get(paste0("disk_p", p,"_kin_v", each))
    
    kin_output$unique_pixels = length(unique(as.numeric(ifu_output$raw_images$voronoi_bins)))
    
    assign(paste0("disk_p", p,"_kin_v", each), kin_output)
    saveRDS(get(paste0("disk_p", p,"_kin_v", each)),
            paste0(data_loc, "/kinematics/disk_sm650_p", p,"_kin_values_vorbin_", each,".Rdata"))
    
    ifu_output = get(paste0("bulge_p", p, "_ifu_v", each))
    kin_output = get(paste0("bulge_p", p,"_kin_v", each))
    
    kin_output$unique_pixels = length(unique(as.numeric(ifu_output$raw_images$voronoi_bins)))
    
    assign(paste0("bulge_p", p,"_kin_v", each), kin_output)
    saveRDS(get(paste0("bulge_p", p,"_kin_v", each)),
            paste0(data_loc, "/kinematics/bulge_sm650_p", p,"_kin_values_vorbin_", each,".Rdata"))
    
  }
}

# Loading data for plotting ====================================================
bulge_control_kin = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_control_kin_values.Rdata"))
bulge_p100_kin =  readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values.Rdata"))
bulge_p84_kin =  readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values.Rdata"))
bulge_p50_kin =  readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values.Rdata"))
bulge_p16_kin =  readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values.Rdata"))
bulge_p0_kin =  readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values.Rdata"))

bulge_p0_kin_v50 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_50.Rdata"))
bulge_p0_kin_v100 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_100.Rdata"))
bulge_p0_kin_v150 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_150.Rdata"))
bulge_p0_kin_v200 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_200.Rdata"))
bulge_p0_kin_v250 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_250.Rdata"))
bulge_p0_kin_v300 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_300.Rdata"))
bulge_p0_kin_v350 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_350.Rdata"))
bulge_p0_kin_v400 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_400.Rdata"))
bulge_p0_kin_v450 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_450.Rdata"))
bulge_p0_kin_v500 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p0_kin_values_vorbin_500.Rdata"))

bulge_p16_kin_v50 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_50.Rdata"))
bulge_p16_kin_v100 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_100.Rdata"))
bulge_p16_kin_v150 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_150.Rdata"))
bulge_p16_kin_v200 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_200.Rdata"))
bulge_p16_kin_v250 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_250.Rdata"))
bulge_p16_kin_v300 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_300.Rdata"))
bulge_p16_kin_v350 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_350.Rdata"))
bulge_p16_kin_v400 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_400.Rdata"))
bulge_p16_kin_v450 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_450.Rdata"))
bulge_p16_kin_v500 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p16_kin_values_vorbin_500.Rdata"))

bulge_p50_kin_v50 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_50.Rdata"))
bulge_p50_kin_v100 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_100.Rdata"))
bulge_p50_kin_v150 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_150.Rdata"))
bulge_p50_kin_v200 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_200.Rdata"))
bulge_p50_kin_v250 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_250.Rdata"))
bulge_p50_kin_v300 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_300.Rdata"))
bulge_p50_kin_v350 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_350.Rdata"))
bulge_p50_kin_v400 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_400.Rdata"))
bulge_p50_kin_v450 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_450.Rdata"))
bulge_p50_kin_v500 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p50_kin_values_vorbin_500.Rdata"))

bulge_p84_kin_v50 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_50.Rdata"))
bulge_p84_kin_v100 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_100.Rdata"))
bulge_p84_kin_v150 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_150.Rdata"))
bulge_p84_kin_v200 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_200.Rdata"))
bulge_p84_kin_v250 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_250.Rdata"))
bulge_p84_kin_v300 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_300.Rdata"))
bulge_p84_kin_v350 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_350.Rdata"))
bulge_p84_kin_v400 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_400.Rdata"))
bulge_p84_kin_v450 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_450.Rdata"))
bulge_p84_kin_v500 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p84_kin_values_vorbin_500.Rdata"))

bulge_p100_kin_v50 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_50.Rdata"))
bulge_p100_kin_v100 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_100.Rdata"))
bulge_p100_kin_v150 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_150.Rdata"))
bulge_p100_kin_v200 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_200.Rdata"))
bulge_p100_kin_v250 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_250.Rdata"))
bulge_p100_kin_v300 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_300.Rdata"))
bulge_p100_kin_v350 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_350.Rdata"))
bulge_p100_kin_v400 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_400.Rdata"))
bulge_p100_kin_v450 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_450.Rdata"))
bulge_p100_kin_v500 = readRDS(paste0(data_loc, "/kinematics/bulge_sm650_p100_kin_values_vorbin_500.Rdata"))


disk_control_kin = readRDS(paste0(data_loc, "/kinematics/disk_sm650_control_kin_values.Rdata"))
disk_p100_kin =  readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values.Rdata"))
disk_p84_kin =  readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values.Rdata"))
disk_p50_kin =  readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values.Rdata"))
disk_p16_kin =  readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values.Rdata"))
disk_p0_kin =  readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values.Rdata"))

disk_p0_kin_v50 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_50.Rdata"))
disk_p0_kin_v100 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_100.Rdata"))
disk_p0_kin_v150 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_150.Rdata"))
disk_p0_kin_v200 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_200.Rdata"))
disk_p0_kin_v250 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_250.Rdata"))
disk_p0_kin_v300 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_300.Rdata"))
disk_p0_kin_v350 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_350.Rdata"))
disk_p0_kin_v400 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_400.Rdata"))
disk_p0_kin_v450 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_450.Rdata"))
disk_p0_kin_v500 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p0_kin_values_vorbin_500.Rdata"))

disk_p16_kin_v50 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_50.Rdata"))
disk_p16_kin_v100 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_100.Rdata"))
disk_p16_kin_v150 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_150.Rdata"))
disk_p16_kin_v200 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_200.Rdata"))
disk_p16_kin_v250 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_250.Rdata"))
disk_p16_kin_v300 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_300.Rdata"))
disk_p16_kin_v350 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_350.Rdata"))
disk_p16_kin_v400 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_400.Rdata"))
disk_p16_kin_v450 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_450.Rdata"))
disk_p16_kin_v500 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p16_kin_values_vorbin_500.Rdata"))

disk_p50_kin_v50 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_50.Rdata"))
disk_p50_kin_v100 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_100.Rdata"))
disk_p50_kin_v150 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_150.Rdata"))
disk_p50_kin_v200 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_200.Rdata"))
disk_p50_kin_v250 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_250.Rdata"))
disk_p50_kin_v300 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_300.Rdata"))
disk_p50_kin_v350 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_350.Rdata"))
disk_p50_kin_v400 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_400.Rdata"))
disk_p50_kin_v450 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_450.Rdata"))
disk_p50_kin_v500 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p50_kin_values_vorbin_500.Rdata"))

disk_p84_kin_v50 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_50.Rdata"))
disk_p84_kin_v100 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_100.Rdata"))
disk_p84_kin_v150 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_150.Rdata"))
disk_p84_kin_v200 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_200.Rdata"))
disk_p84_kin_v250 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_250.Rdata"))
disk_p84_kin_v300 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_300.Rdata"))
disk_p84_kin_v350 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_350.Rdata"))
disk_p84_kin_v400 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_400.Rdata"))
disk_p84_kin_v450 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_450.Rdata"))
disk_p84_kin_v500 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p84_kin_values_vorbin_500.Rdata"))

disk_p100_kin_v50 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_50.Rdata"))
disk_p100_kin_v100 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_100.Rdata"))
disk_p100_kin_v150 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_150.Rdata"))
disk_p100_kin_v200 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_200.Rdata"))
disk_p100_kin_v250 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_250.Rdata"))
disk_p100_kin_v300 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_300.Rdata"))
disk_p100_kin_v350 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_350.Rdata"))
disk_p100_kin_v400 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_400.Rdata"))
disk_p100_kin_v450 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_450.Rdata"))
disk_p100_kin_v500 = readRDS(paste0(data_loc, "/kinematics/disk_sm650_p100_kin_values_vorbin_500.Rdata"))

# Plots ------------------------------------------------------------------------

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")

part_seq = c(10, seq(50,500,by=50))
source("functions/plot_images.R")

# Visualising the binning process ----------------------------------------------
w = plot_array(0.36, 4)

bulge_000 = readRDS(paste0(cube_loc, "bulge_sm650_p0_ifu_obs.Rdata"))
bulge_100 = readRDS(paste0(cube_loc, "bulge_sm650_p0_ifu_obs_vorbin_100.Rdata"))
bulge_200 = readRDS(paste0(cube_loc, "bulge_sm650_p0_ifu_obs_vorbin_200.Rdata"))
bulge_300 = readRDS(paste0(cube_loc, "bulge_sm650_p0_ifu_obs_vorbin_300.Rdata"))

bulge_photometry = readRDS(paste0(data_loc, "bulge_sm650_control_photo_values_500.Rdata"))
ap_mask = matrix(data = bulge_000$observation$aperture_region, ncol = bulge_000$observation$sbin, nrow = bulge_000$observation$sbin)
ap_mask[ap_mask == 0] = NA

tenpc = 10/bulge_000$observation$sbin_size

pdf.options(fonts="Arial")
png(file=paste0(plot_loc, "voronoi_binning_visual.png"), 
    width = 11, height = 7, units = "in", res = 200)

par(fig = c(0,0.5,0.46,1),
    mar = c(3.5,3.5,0.5,0.5), pty= "m")
magplot(0, type="n",
        ylab = expression('\u0394'*'\u03BB'['R']),
        xlim = c(0,500), xaxt="n",
        ylim = c(-0.07, 0.2))
abline(h=0, lwd=3, col = tol_pal[1])
text(435, 0.175, "Bulge Models", font = 2, cex=1.1)
lines(part_seq,
      c((bulge_p0_kin$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v50$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v100$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v150$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v200$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v250$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v300$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v350$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v400$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v450$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p0_kin_v500$lambda_re - bulge_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[6])
lines(part_seq,
      c((bulge_p16_kin$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v50$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v100$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v150$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v200$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v250$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v300$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v350$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v400$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v450$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p16_kin_v500$lambda_re - bulge_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[5])
lines(part_seq,
      c((bulge_p50_kin$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v50$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v100$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v150$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v200$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v250$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v300$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v350$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v400$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v450$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p50_kin_v500$lambda_re - bulge_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[4])
lines(part_seq,
      c((bulge_p84_kin$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v50$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v100$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v150$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v200$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v250$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v300$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v350$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v400$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p84_kin_v450$lambda_re - bulge_control_kin$lambda_re), 
        (bulge_p84_kin_v500$lambda_re - bulge_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[3])
lines(part_seq,
      c((bulge_p100_kin$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v50$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v100$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v150$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v200$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v250$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v300$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v350$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v400$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v450$lambda_re - bulge_control_kin$lambda_re),
        (bulge_p100_kin_v500$lambda_re - bulge_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[2])

par(fig = c(0,0.5,0,0.54),
    mar = c(3.5,3.5,0.5,0.5), new = T, pty="m")
magplot(0, type="n",
        ylab = expression('\u0394'*'\u03BB'['R']),
        xlab = "Particles per voronoi bin",
        xlim = c(0,500), 
        ylim = c(-0.07, 0.2))
abline(h=0, lwd=3, col = tol_pal[1])
text(440, 0.175, "Disc Models", font = 2, cex=1.1)
lines(part_seq,
      c((disk_p0_kin$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v50$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v100$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v150$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v200$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v250$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v300$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v350$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v400$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v450$lambda_re - disk_control_kin$lambda_re),
        (disk_p0_kin_v500$lambda_re - disk_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[6])
lines(part_seq,
      c((disk_p16_kin$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v50$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v100$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v150$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v200$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v250$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v300$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v350$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v400$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v450$lambda_re - disk_control_kin$lambda_re),
        (disk_p16_kin_v500$lambda_re - disk_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[5])
lines(part_seq,
      c((disk_p50_kin$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v50$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v100$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v150$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v200$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v250$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v300$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v350$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v400$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v450$lambda_re - disk_control_kin$lambda_re),
        (disk_p50_kin_v500$lambda_re - disk_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[4])
lines(part_seq,
      c((disk_p84_kin$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v50$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v100$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v150$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v200$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v250$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v300$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v350$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v400$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v450$lambda_re - disk_control_kin$lambda_re),
        (disk_p84_kin_v500$lambda_re - disk_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[3])
lines(part_seq,
      c((disk_p100_kin$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v50$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v100$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v150$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v200$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v250$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v300$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v350$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v400$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v450$lambda_re - disk_control_kin$lambda_re),
        (disk_p100_kin_v500$lambda_re - disk_control_kin$lambda_re)),
      lwd = 3, col = tol_pal[2])

legend("topleft", inset = c(0.05, 0.05),
       legend=c("control", "p100", "p84", "p50", "p16", "p0"),
       col=c(tol_pal[c(1,2,3,4,5,6)]), pch=15, bty="n",
       cex = 1)

par(mar = c(3.5,3.5,0.5,0.5), new = T)
plot_velocity(bulge_000$observed_images$velocity_image*ap_mask,
              fig = c(w[3,1], w[3,2], 0.64, 1), new = T,
              radii = list("a" = bulge_photometry$major_axis_a_kpc/bulge_000$observation$sbin_size,
                           "b" = bulge_photometry$minor_axis_b_kpc/bulge_000$observation$sbin_size,
                           "ang" = bulge_photometry$primary_axis), labN = 3, zlim = c(-50,50))
shadowtext(x = 21, y = 27, labels = "No minimum")
lines(c((15-tenpc), 15), c(4,4), lwd=4, col="white")
lines(c((15-tenpc), 15), c(4,4), lwd=2, col="black")
shadowtext((15-(tenpc/2)), 7, labels="10 kpc")

par(mar = c(3.5,3.5,0.5,0.5))
plot_velocity(bulge_200$observed_images$velocity_image*ap_mask,
              fig = c(w[4,1], w[4,2], 0.64, 1), new=T,
              radii = list("a" = bulge_photometry$major_axis_a_kpc/bulge_000$observation$sbin_size,
                           "b" = bulge_photometry$minor_axis_b_kpc/bulge_000$observation$sbin_size,
                           "ang" = bulge_photometry$primary_axis), labN = 3, zlim = c(-50,50))
shadowtext(x = 21, y = 27, labels = expression("N"["target"]*" = 200"))

par(mar = c(3.5,3.5,0.5,0.5))
plot_dispersion(bulge_000$observed_images$dispersion_image*ap_mask,
                fig = c(w[3,1], w[3,2], 0.32, 0.68),  new=T,
                radii = list("a" = bulge_photometry$major_axis_a_kpc/bulge_000$observation$sbin_size,
                             "b" = bulge_photometry$minor_axis_b_kpc/bulge_000$observation$sbin_size,
                             "ang" = bulge_photometry$primary_axis), labN = 3, zlim = c(50,200))

par(mar = c(3.5,3.5,0.5,0.5))
plot_dispersion(bulge_200$observed_images$dispersion_image*ap_mask,
                fig = c(w[4,1], w[4,2], 0.32, 0.68), new=T,
                radii = list("a" = bulge_photometry$major_axis_a_kpc/bulge_000$observation$sbin_size,
                             "b" = bulge_photometry$minor_axis_b_kpc/bulge_000$observation$sbin_size,
                             "ang" = bulge_photometry$primary_axis), labN = 3, zlim = c(50,200))



par(fig = c(w[3,1], w[3,2], 0, 0.36), mar = c(3.5,3.5,0.5,0.5), new = T, pty="s")
magplot(x = seq(1, length(as.numeric(bulge_000$raw_images$particle_image[bulge_000$raw_images$particle_image > 0]))),
        sort(as.numeric(bulge_000$raw_images$particle_image[bulge_000$raw_images$particle_image > 0]), decreasing = T),
        type = "p", pch = 16, col = "black",
        ylim = c(0,1500),
        xlab = "Pixel ID",
        ylab = "Particles per voronoi bin")

b200 = 
  data.table::data.table("part" = as.numeric(bulge_200$raw_images$particle_image),
                         "bin" = as.numeric(bulge_200$raw_images$voronoi_bins))

b200_vbin = b200[, list(partN = sum(part)), by="bin"]

par(fig = c(w[4,1], w[4,2], 0, 0.36), mar = c(3.5,3.5,0.5,0.5), new = T, pty="s")
magplot(x = seq(1, length(b200_vbin$partN[b200_vbin$partN>0])),
        sort(b200_vbin$partN[b200_vbin$partN>0], decreasing = T),
        type = "p", pch = 16, col = "black",
        ylim = c(0,1500), 
        xlab = "Pixel ID", minorn=5, majorn = 3)
abline(h = 200, col = "grey48", lwd=2)

dev.off()

# Examining the LOSVD of pixels in increasing vorbin pixels --------------------
unbinned_output = readRDS("data/cubes/flux_obs/bulge_sm650_control_ifu_obs.Rdata")
unbinned_output$velocity_cube = unbinned_output$velocity_cube[,,11:120]
unbinned_output$observation$vbin_seq = unbinned_output$observation$vbin_seq[11:120]

vbinned_output_200 = readRDS("data/cubes/flux_obs/bulge_sm650_p0_ifu_obs_vorbin_200.Rdata")
vbinned_output_300 = readRDS("data/cubes/flux_obs/bulge_sm650_p0_ifu_obs_vorbin_300.Rdata")
vbinned_output_400 = readRDS("data/cubes/flux_obs/bulge_sm650_p0_ifu_obs_vorbin_400.Rdata")
vbinned_output_500 = readRDS("data/cubes/flux_obs/bulge_sm650_p0_ifu_obs_vorbin_500.Rdata")

unbinned_output_p0 = readRDS("data/cubes/flux_obs/bulge_sm650_p0_ifu_obs.Rdata")

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")


.losvd_out = function(x, vel, sig){
  k=1
  w = (x - vel)/sig
  measured_vlos = (k * exp(-0.5*(w^2))) 
  return(measured_vlos)
}

.losvd_fit = function(par, x, losvd){
  vel = par[1]
  sig = par[2]
  k  = 1
  w = (x - vel)/sig
  measured_vlos = (k * exp(-0.5*(w^2))) 
  return=sum((measured_vlos-losvd)^2)
}

draw_bin = function(heights,widths,col){
  
  binwidth = diff(widths)[1]/2
  x1 = widths - binwidth
  x2 = widths + binwidth
  col = as.double(col2rgb(col))/256
  
  for (each in 1:length(heights)){
    polygon(c(x1[each], x1[each], x2[each], x2[each], x1[each]), 
            c(0, heights[each], heights[each], 0, 0),
            col=rgb(col[1], col[2], col[3], alpha = 0.5), border = NA)
  }
}


xpix = 9
ypix = 8

unbinned_output$raw_images$particle_image[xpix,ypix]
unbinned_output_p0$raw_images$particle_image[xpix,ypix]

unbinned_output$raw_images$dispersion_image[xpix,ypix]
unbinned_output_p0$raw_images$dispersion_image[xpix,ypix]
vbinned_output_200$raw_images$dispersion_image[xpix,ypix]

c(sum(vbinned_output_500$raw_images$particle_image[vbinned_output_500$raw_images$voronoi_bins == vbinned_output_500$raw_images$voronoi_bins[xpix,ypix]]), 
  sum(vbinned_output_400$raw_images$particle_image[vbinned_output_400$raw_images$voronoi_bins == vbinned_output_400$raw_images$voronoi_bins[xpix,ypix]]),
  sum(vbinned_output_300$raw_images$particle_image[vbinned_output_300$raw_images$voronoi_bins == vbinned_output_300$raw_images$voronoi_bins[xpix,ypix]]), 
  sum(vbinned_output_200$raw_images$particle_image[vbinned_output_200$raw_images$voronoi_bins == vbinned_output_200$raw_images$voronoi_bins[xpix,ypix]]))

# Vorbin comparisons
velocity_spectrum_control = unbinned_output$velocity_cube[xpix,ypix,]#colSums(unbinned_output$velocity_cube, dim = 2, na.rm=T)/sum(mask, na.rm=T)
velocity_spectrum_control = velocity_spectrum_control / max(velocity_spectrum_control)

vel_ini = unbinned_output$raw_images$velocity_image[xpix,ypix]
sd_ini  = unbinned_output$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_control_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                               fn    = .losvd_fit,
                                               x     = unbinned_output$observation$vbin_seq,
                                               losvd = velocity_spectrum_control,
                                               method="BFGS", control=list(reltol=1e-9))$par

velocity_spectrum_control_fit = .losvd_out(x = seq(-650,650, by=5),
                                           vel = velocity_spectrum_control_optim[1],
                                           sig = velocity_spectrum_control_optim[2])

# unbinned p0
velocity_spectrum_p0 = unbinned_output_p0$velocity_cube[xpix,ypix,]
velocity_spectrum_p0 = velocity_spectrum_p0 / max(velocity_spectrum_p0)

vel_ini = unbinned_output_p0$raw_images$velocity_image[xpix,ypix]
sd_ini  = unbinned_output_p0$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_p0_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                          fn    = .losvd_fit,
                                          x     = unbinned_output_p0$observation$vbin_seq,
                                          losvd = velocity_spectrum_p0,
                                          method="BFGS", control=list(reltol=1e-9))$par

velocity_spectrum_p0_fit = .losvd_out(x = seq(-650,650, by=5),
                                      vel = velocity_spectrum_p0_optim[1],
                                      sig = velocity_spectrum_p0_optim[2])

#vorbin 200
velocity_spectrum_v200 = vbinned_output_200$velocity_cube[xpix,ypix,]
velocity_spectrum_v200 = velocity_spectrum_v200 / max(velocity_spectrum_v200)
vel_ini = vbinned_output_200$raw_images$velocity_image[xpix,ypix]
sd_ini = vbinned_output_200$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_v200_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                            fn    = .losvd_fit,
                                            x     = vbinned_output_200$observation$vbin_seq,
                                            losvd = velocity_spectrum_v200,
                                            method="BFGS", control=list(reltol=1e-9))$par
velocity_spectrum_v200_fit = .losvd_out(x = seq(-650,650, by=5),
                                        vel = velocity_spectrum_v200_optim[1],
                                        sig = velocity_spectrum_v200_optim[2])

#vorbin 300
velocity_spectrum_v300 = vbinned_output_300$velocity_cube[xpix,ypix,]
velocity_spectrum_v300 = velocity_spectrum_v300 / max(velocity_spectrum_v300)
vel_ini = vbinned_output_300$raw_images$velocity_image[xpix,ypix]
sd_ini = vbinned_output_300$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_v300_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                            fn    = .losvd_fit,
                                            x     = vbinned_output_300$observation$vbin_seq,
                                            losvd = velocity_spectrum_v300,
                                            method="BFGS", control=list(reltol=1e-9))$par
velocity_spectrum_v300_fit = .losvd_out(x = seq(-650,650, by=5),
                                        vel = velocity_spectrum_v300_optim[1],
                                        sig = velocity_spectrum_v300_optim[2])


#vorbin 400
velocity_spectrum_v400 = vbinned_output_400$velocity_cube[xpix,ypix,]
velocity_spectrum_v400 = velocity_spectrum_v400 / max(velocity_spectrum_v400)
vel_ini = vbinned_output_400$raw_images$velocity_image[xpix,ypix]
sd_ini = vbinned_output_400$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_v400_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                            fn    = .losvd_fit,
                                            x     = vbinned_output_400$observation$vbin_seq,
                                            losvd = velocity_spectrum_v400,
                                            method="BFGS", control=list(reltol=1e-9))$par
velocity_spectrum_v400_fit = .losvd_out(x = seq(-650,650, by=5),
                                        vel = velocity_spectrum_v400_optim[1],
                                        sig = velocity_spectrum_v400_optim[2])

#vorbin 500
velocity_spectrum_v500 = vbinned_output_500$velocity_cube[xpix,ypix,]
velocity_spectrum_v500 = velocity_spectrum_v500 / max(velocity_spectrum_v500)
vel_ini = vbinned_output_500$raw_images$velocity_image[xpix,ypix]
sd_ini = vbinned_output_500$raw_images$dispersion_image[xpix,ypix]
velocity_spectrum_v500_optim = stats::optim(par   = c(vel_ini,sd_ini),
                                            fn    = .losvd_fit,
                                            x     = vbinned_output_500$observation$vbin_seq,
                                            losvd = velocity_spectrum_v500,
                                            method="BFGS", control=list(reltol=1e-9))$par
velocity_spectrum_v500_fit = .losvd_out(x = seq(-650,650, by=5),
                                        vel = velocity_spectrum_v500_optim[1],
                                        sig = velocity_spectrum_v500_optim[2])

# Anderson-Darling Normality test ----------------------------------------------
unbinned_distribution = rep(unbinned_output$observation$vbin_seq, ((velocity_spectrum_control/max(velocity_spectrum_control_fit))*100))
unbinned_distribution_p0  = rep(unbinned_output_p0$observation$vbin_seq, (velocity_spectrum_p0/max(velocity_spectrum_p0_fit))*100)
vbinned_distribution_v200 = rep(vbinned_output_200$observation$vbin_seq,((velocity_spectrum_v200/max(velocity_spectrum_v200_fit))*100))

library(tseries)

t=maghist(unbinned_distribution, freq = T, plot=F, breaks = (unbinned_output$observation$vbin_seq + diff(unbinned_output$observation$vbin_seq)[1]/2))
magplot(unbinned_output$observation$vbin_seq,
        (velocity_spectrum_control/max(velocity_spectrum_control_fit))*max(t$counts), pch = 16,
        xlim = c(-600,600), ylim = c(0, 120), col = tol_pal[1],
        ylab = "Normalised LOSVD",
        xaxt="n", cex.lab=1.2)
maghist(unbinned_distribution, freq = T, add=T, breaks = (unbinned_output$observation$vbin_seq + diff(unbinned_output$observation$vbin_seq)[1]/2))
jarque.bera.test(unbinned_distribution)

maghist(unbinned_distribution_p0)
jarque.bera.test(unbinned_distribution_p0)

maghist(vbinned_distribution_v200)
jarque.bera.test(vbinned_distribution_v200)

# Figure showing LOSVD =========================================================
pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "voronoi_binning_effect_LOSVD_new.pdf"), 
          width = 14, height = 10)

par(fig = c(0,0.55,0.24,1), cex=1.2, cex.lab=1.2)
magplot(unbinned_output$observation$vbin_seq,
        velocity_spectrum_control/max(velocity_spectrum_control_fit), pch = 16,
        xlim = c(-600,600), ylim = c(0, 1.2), col = tol_pal[1],
        ylab = "Normalised LOSVD",
        xaxt="n", cex.lab=1.2)
draw_bin(widths = unbinned_output$observation$vbin_seq,
         heights = velocity_spectrum_control/max(velocity_spectrum_control_fit), 
         col = tol_pal[1])
lines(seq(-650,650, by=5),
      velocity_spectrum_control_fit/max(velocity_spectrum_control_fit),
      col = "white", lwd = 4)
lines(seq(-650,650, by=5),
      velocity_spectrum_control_fit/max(velocity_spectrum_control_fit),
      col = tol_pal[1], lwd = 2)
text(390,1.15, "Unbinned observations", font =2)


points(unbinned_output_p0$observation$vbin_seq,
       velocity_spectrum_p0/max(velocity_spectrum_p0_fit), pch = 16,
       xlim = c(-600,600), col = tol_pal[6])
draw_bin(widths = unbinned_output_p0$observation$vbin_seq,
         heights = velocity_spectrum_p0/max(velocity_spectrum_p0_fit), 
         col = tol_pal[6])
lines(seq(-650,650, by=5),
      velocity_spectrum_p0_fit/max(velocity_spectrum_p0_fit),
      col = "white", lwd = 4)
lines(seq(-650,650, by=5),
      velocity_spectrum_p0_fit/max(velocity_spectrum_p0_fit),
      col = tol_pal[6], lwd = 2)

points(unbinned_output$observation$vbin_seq,
       velocity_spectrum_control/max(velocity_spectrum_control), pch = 16,
       xlim = c(-600,600), col = tol_pal[1])

legend("topleft", inset = c(0.02, 0.02),
       legend=c("control pixel - 2951 particles", "p0 pixel - 89 particles"),
       col=c(tol_pal[c(1,6)]), pch=15, bty="n",
       cex = 1)

text(-450, 1, expression(paste(sigma)["sd"]*" = 148 km/s"), col =  tol_pal[1],  font = 2)
text(-420, 0.925, expression(paste(sigma)["2-mom"]*" = 147 km/s"), col =  tol_pal[1], font = 2)
text(-420, 0.85, expression("JB "*italic("p")*"-value = 0.90"), col =  tol_pal[1], font = 2)

text(450, 1, expression(paste(sigma)["sd"]*" = 140 km/s"), col =  tol_pal[6],  font = 2)
text(420, 0.925, expression(paste(sigma)["2-mom"]*" = 71 km/s"), col =  tol_pal[6],  font = 2)
text(420, 0.85, expression("JB "*italic("p")*"-value = 0.99"), col =  tol_pal[6], font = 2)


# ------------------------------------------------------------------------------

par(fig = c(0,0.55,0,0.45), new=T, cex=1.2, cex.lab=1.2)
magplot(unbinned_output$observation$vbin_seq,
        .losvd_out(x = unbinned_output$observation$vbin_seq,
                   vel = velocity_spectrum_p0_optim[1],
                   sig = velocity_spectrum_p0_optim[2]) - velocity_spectrum_control,
        pch = 16, col = tol_pal[6], xlab=expression("v"["LOS"]*", km/s"),
        ylab = expression('\u0394'["X - control"]),
        xlim = c(-600,600), ylim = c(-0.55, 0.3), cex.lab=1.2)

abline(h=0, col = "black", lwd=2)

abline(h=sd(.losvd_out(x = unbinned_output$observation$vbin_seq,
                       vel = velocity_spectrum_p0_optim[1],
                       sig = velocity_spectrum_p0_optim[2]) - velocity_spectrum_control), 
       col = tol_pal[6], lwd=2)


# ------------------------------------------------------------------------------
vpal = c(rgb(126, 41, 84, maxColorValue = 256),
         rgb(159, 74,150, maxColorValue = 256),
         rgb(194,106,119, maxColorValue = 256),
         rgb(220,205,125, maxColorValue = 256))


par(fig = c(0.45,1,0.24,1), new=T, cex=1.2, cex.lab=1.2)
magplot(unbinned_output$observation$vbin_seq,
        velocity_spectrum_control/max(velocity_spectrum_control_fit), pch = 16,
        xlim = c(-650,650), ylim = c(0, 1.2), col = tol_pal[1],
        xaxt="n", yaxt="n")
draw_bin(widths = unbinned_output$observation$vbin_seq,
         heights = velocity_spectrum_control/max(velocity_spectrum_control_fit), 
         col = tol_pal[1])
lines(seq(-650,650, by=5),
      velocity_spectrum_control_fit/max(velocity_spectrum_control_fit),
      col = "white", lwd = 4)
lines(seq(-650,650, by=5),
      velocity_spectrum_control_fit/max(velocity_spectrum_control_fit),
      col = tol_pal[1], lwd = 2)
text(360,1.15, "Voronoi-binned observations", font =2)

points(vbinned_output_200$observation$vbin_seq,
       velocity_spectrum_v200/max(velocity_spectrum_v200_fit), pch = 16,
       xlim = c(-600,600), col = vpal[4])
draw_bin(widths = vbinned_output_200$observation$vbin_seq,
         heights = velocity_spectrum_v200/max(velocity_spectrum_v200_fit), 
         col = vpal[4])
lines(seq(-650,650, by=5),
      velocity_spectrum_v200_fit/max(velocity_spectrum_v200_fit),
      col = "white", lwd = 4)
lines(seq(-650,650, by=5),
      velocity_spectrum_v200_fit/max(velocity_spectrum_v200_fit),
      col = vpal[4], lwd = 2)

legend("topleft", inset = c(0.02, 0.02),
       legend=c("control pixel - 2951 particles", "p0 pixel - 203 particles"),
       col=c(tol_pal[1],vpal[4]), pch=15, bty="n",
       cex = 1)

text(-480, 1, expression(paste(sigma)["sd"]*" = 146 km/s"), col =  vpal[4], font = 2)
text(-470, 0.925, expression(paste(sigma)["2-mom"]*" = 126 km/s"), col =  vpal[4], font = 2)
text(-470, 0.85, expression("JB "*italic("p")*"-value = 0.19"), col =  vpal[4], font = 2)

# ------------------------------------------------------------------------------

par(fig = c(0.45,1,0,0.45), new=T, cex=1.2, cex.lab=1.2)
magplot(unbinned_output$observation$vbin_seq,
        .losvd_out(x = unbinned_output$observation$vbin_seq,
                   vel = velocity_spectrum_v200_optim[1],
                   sig = velocity_spectrum_v200_optim[2]) - velocity_spectrum_control,
        pch = 15, col = vpal[4], xlab=expression("v"["LOS"]*", km/s"),
        xlim = c(-600,600), ylim = c(-0.55, 0.3), yaxt="n", cex.lab=1.2)

points(unbinned_output$observation$vbin_seq,
       .losvd_out(x = unbinned_output$observation$vbin_seq,
                  vel = velocity_spectrum_v300_optim[1],
                  sig = velocity_spectrum_v300_optim[2]) - velocity_spectrum_control,
       pch = 15, col = vpal[3])

points(unbinned_output$observation$vbin_seq,
       .losvd_out(x = unbinned_output$observation$vbin_seq,
                  vel = velocity_spectrum_v400_optim[1],
                  sig = velocity_spectrum_v400_optim[2]) - velocity_spectrum_control,
       pch = 15, col = vpal[2])

points(unbinned_output$observation$vbin_seq,
       .losvd_out(x = unbinned_output$observation$vbin_seq,
                  vel = velocity_spectrum_v500_optim[1],
                  sig = velocity_spectrum_v500_optim[2]) - velocity_spectrum_control,
       pch = 15, col = vpal[1])

points(unbinned_output$observation$vbin_seq,
       .losvd_out(x = unbinned_output$observation$vbin_seq,
                  vel = velocity_spectrum_v200_optim[1],
                  sig = velocity_spectrum_v200_optim[2]) - velocity_spectrum_control,
       pch = 15, col = vpal[4])

abline(h=0, col = "black", lwd=2)

abline(h=sd(.losvd_out(x = unbinned_output$observation$vbin_seq,
                       vel = velocity_spectrum_v300_optim[1],
                       sig = velocity_spectrum_v300_optim[2])- velocity_spectrum_control) , 
       col = vpal[3], lwd=2)

abline(h=sd(.losvd_out(x = unbinned_output$observation$vbin_seq,
                       vel = velocity_spectrum_v400_optim[1],
                       sig = velocity_spectrum_v400_optim[2])-velocity_spectrum_control), 
       col = vpal[2], lwd=2)

abline(h=sd(.losvd_out(x = unbinned_output$observation$vbin_seq,
                       vel = velocity_spectrum_v500_optim[1],
                       sig = velocity_spectrum_v500_optim[2])-velocity_spectrum_control), 
       col = vpal[1], lwd=2)

abline(h=sd(.losvd_out(x = unbinned_output$observation$vbin_seq,
                       vel = velocity_spectrum_v200_optim[1],
                       sig = velocity_spectrum_v200_optim[2])-velocity_spectrum_control), 
       col = vpal[4], lwd=2)

legend("bottomright", inset = c(0.05, 0.05),
       legend=c("200", "300", "400", "500"),
       col=c(vpal[c(4,3,2,1)]), pch=15, bty="n",
       cex = 1)

dev.off()

# Quantitative assessment on what portion of lambda R change is dictated by sigma ------------------

(1/c((bulge_p0_kin_v50$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v100$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v150$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v200$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v250$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v300$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v350$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v400$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v450$sigma_m - bulge_control_kin$sigma_m),
     (bulge_p0_kin_v500$sigma_m - bulge_control_kin$sigma_m)) ) /
  c((bulge_p0_kin_v50$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v100$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v150$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v200$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v250$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v300$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v350$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v400$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v450$lambda_re - bulge_control_kin$lambda_re),
    (bulge_p0_kin_v500$lambda_re - bulge_control_kin$lambda_re)) 

(1/c((disk_p0_kin_v50$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v100$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v150$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v200$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v250$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v300$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v350$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v400$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v450$sigma_m - disk_control_kin$sigma_m),
     (disk_p0_kin_v500$sigma_m - disk_control_kin$sigma_m)) ) /
  c((disk_p0_kin_v50$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v100$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v150$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v200$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v250$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v300$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v350$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v400$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v450$lambda_re - disk_control_kin$lambda_re),
    (disk_p0_kin_v500$lambda_re - disk_control_kin$lambda_re)) 




