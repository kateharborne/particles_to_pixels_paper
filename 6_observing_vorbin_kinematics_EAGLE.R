# Kate Harborne 30/10/23
# ============================================================================ #
# Measure the kinematics within the half-light radius, using standard pixels 
# and vorbin pixels.
# This script provides an example of how the EAGLE galaxies have been observed 
# using SimSpin for completeness. These scripts were run using the Setonix 
# Supercomputer from the Pawsey Supercomputing Centre between 2023-24. For the 
# exact details used for the paper, please contact the corresponding author at 
# katherine.harborne@uwa.edu.au.
# ============================================================================ #

library(SimSpin)
library(magicaxis)

# Reading in the EAGLE galaxy data ---------------------------------------------
batch = commandArgs(trailingOnly = T)
batch_size = 25

# PATHS ========================================================================
setonix_root = ""

ss_file_loc = paste0(setonix_root, "/simspin_files/")
ss_files = list.files(path = ss_file_loc)
kin_out_loc = paste0(setonix_root, "/galaxy_info/")
functions_loc = paste0(setonix_root, "functions/")
cubes_loc = paste0(setonix_root, "/cubes/")
logs_loc  = paste0(setonix_root,"/logs/")

galaxy_centres_file = paste0(setonix_root, "/RefL0050N0752_subhalos_1e10.txt")

# RECORD FILES =================================================================
obs_prop_file = paste0(kin_out_loc, "observation_properties_mom2", batch,".txt")
ub_kin_prop_file = paste0(kin_out_loc, "kinematic_properties_mom2_unbinned_", batch,".txt")
vb_kin_prop_file = paste0(kin_out_loc, "kinematic_properties_mom2_voronoi_binned_", batch,".txt")
log_file      = paste0(logs_loc, "kinematic_observations_mom2_log_", batch,".txt")

# Writing the opening lines for output file obs_prop_file ----------------------
cat(paste0('Observation summaries of RefL0050N0752 galaxies from snap 28 using SimSpin ', packageVersion('SimSpin'), ' on date: ', Sys.Date()), '/ batch no: ', batch,
    sep='\n', file=obs_prop_file)
cat(paste0('-------------------------------------------------------------------------------------------------'),
    sep='\n', file=obs_prop_file, append=T)
cat(paste0("'GalaxyID','inc_deg','twist_deg','dist_z','seeing_fwhm','major_axis_a_kpc','minor_axis_b_kpc','angle','ellipticity','sersic_index'"),
    sep='\n', file=obs_prop_file, append=T)

# Writing the opening lines for output file ub_kin_prop_file -------------------
cat(paste0('Unbinned kinematic summaries of RefL0050N0752 galaxies from snap 28 using SimSpin ', packageVersion('SimSpin'), ' on date: ', Sys.Date()), '/ batch no: ', batch,
    sep='\n', file=ub_kin_prop_file)
cat(paste0('-------------------------------------------------------------------------------------------------'),
    sep='\n', file=ub_kin_prop_file, append=T)
cat(paste0("'GalaxyID','lambda_re','vsigma','j','sigma_m','rotational_velocity','unique_pixels'"),
    sep='\n', file=ub_kin_prop_file, append=T)

# Writing the opening lines for output file vb_kin_prop_file -------------------
cat(paste0('Voronoi-binned kinematic summaries of RefL0050N0752 galaxies from snap 28 using SimSpin ', packageVersion('SimSpin'), ' on date: ', Sys.Date()), '/ batch no: ', batch,
    sep='\n', file=vb_kin_prop_file)
cat(paste0('-------------------------------------------------------------------------------------------------'),
    sep='\n', file=vb_kin_prop_file, append=T)
cat(paste0("'GalaxyID','vorbin_number','lambda_re','vsigma','j','sigma_m','rotational_velocity','unique_pixels'"),
    sep='\n', file=vb_kin_prop_file, append=T)

# Log file ---------------------------------------------------------------------
cat(paste0("New observation set generation using SimSpin ", packageVersion("SimSpin"), " for batch ", batch, " on date: ", Sys.Date()),
    sep="\n", file=log_file)
cat(paste0("-------------------------------------------------------------------------------------------------"),
    sep="\n", file=log_file, append=T)
cat(paste0(" "),
    sep="\n", file=log_file, append=T)

# FUNCTIONS ====================================================================
source(paste0(functions_loc, "get_kinematic_properties.R"))
source(paste0(functions_loc, "get_photometric_properties.R"))

# MAIN =========================================================================

# Setting the observing telescope to be SAMI
ifu = telescope(type="SAMI", 
                signal_to_noise = NA)

# Initialising the photometric observation telescope - a KIDS equivalent
photo = telescope(type = "IFU",  
                  fov = 100, 
                  filter = "r",
                  aperture_shape = "square", 
                  spatial_res = 0.2, 
                  wave_res = 200,
                  signal_to_noise = NA)

galaxy_index = c(seq(1, length(ss_files), by=batch_size), length(ss_files))

eagle_data = read.csv(galaxy_centres_file, sep=",", comment.char = "#", header = T)

ob = 1

for (each in galaxy_index[as.numeric(batch)]:(galaxy_index[as.numeric(batch)+1]-1)){
  
  galaxy_id = stringr::str_remove(
    stringr::str_remove(ss_files[each], "SimSpin_EAGLE_snap28_50kpc_with_galaxyID_"),
    ".Rdata")
  
  ss_galaxy_file = paste0(ss_file_loc, ss_files[each])
  
  strategy = observing_strategy(dist_kpc_per_arcsec = 1, 
                                inc_deg = 70, 
                                twist_deg = 0,
                                blur = F)
  
  cat(paste0("BEGINNING OBSERVATION ", ob, " OF ", batch_size, " IN BATCH ", batch, " of galaxy ", galaxy_id),
      sep="\n", file=log_file, append=T)
  
  # Properties of the galaxy to be observed
  galaxy_id     = as.numeric(galaxy_id)
  index         = which(eagle_data$GalaxyID == galaxy_id)
  gal_mass      = eagle_data$MassType_Star[index]
  
  # Step 1: Construct an unblurred cube at KiDS resolution to compute photometric properties
  
  photo_mass_output = tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                                     telescope = photo, 
                                                     observing_strategy = strategy,
                                                     method = "velocity",
                                                     write_fits = F, verbose = F, cores = 1, 
                                                     mass_flag = T, moments = 2),
                               error = function(e){NA})
  
  mass_frac = sum(photo_mass_output$raw_images$mass_image, na.rm=T)/gal_mass
  
  mass_frac = c(mass_frac, 1)
  while (mass_frac[1] < 0.85 & mass_frac[2] != 0 & as.numeric(kpc_per_arcsec(strategy$distance)) < 8){ 
    # if the mass fraction in the image is too small, grow until sufficient to measure kinematics
    strategy = observing_strategy(dist_kpc_per_arcsec = as.numeric(kpc_per_arcsec(strategy$distance)+0.25), 
                                  inc_deg = 70, 
                                  twist_deg = 0,
                                  blur = F)
    
    photo_mass_output = tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                                       telescope = photo, 
                                                       observing_strategy = strategy,
                                                       method = "velocity",
                                                       write_fits = F, verbose = F, cores = 1, 
                                                       mass_flag = T, moments = 2),
                                 error = function(e){NA})
    
    
    mass_frac[2] = (sum(photo_mass_output$raw_images$mass_image, na.rm=T)/gal_mass) - mass_frac[1] 
    mass_frac[1] = sum(photo_mass_output$raw_images$mass_image, na.rm=T)/gal_mass
    
  } 
  
  if (mass_frac[1] < 0.85 & mass_frac[2] != 0 & as.numeric(kpc_per_arcsec(strategy$distance)) == 8){
    photo_output = NA
  } else {
    photo_output = tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                                  telescope = photo, 
                                                  observing_strategy = strategy,
                                                  method = "velocity",
                                                  write_fits = F, verbose = F, cores = 1, 
                                                  mass_flag = F, moments = 2),
                            error = function(e){NA})
  }
  
  if (all(is.na(photo_output))){
    # If that function returns NA's, the observation angle is not going to work for the analysis further on. 
    # Ditch this unbiased sample of observing angles and draw another set. 
    cat(paste0("ERROR observing angle ", ob, ". Cannot get 90% of the mass in field for GalaxyID ", galaxy_id), 
        sep= "\n", file = log_file, append = T)
    
  } else {
    cat(paste0(" -- Measuring photometry of galaxy -- "), 
        sep= "\n", file = log_file, append = T)
    
    # Step 2: Measuring the observed half-mass radius from the image 
    photometry_properties = get_photometry_properties(photo_flux_image = photo_output$observed_images$flux_image,
                                                      obs_summary = photo_output$observation,
                                                      rband_limit = 24.9, mass_frac=mass_frac,
                                                      plot=F)
    
    # Step 3: Observing the kinematics without binning
    
    unbinned_output = 
      build_datacube(simspin_file = ss_galaxy_file, 
                     telescope = ifu, 
                     observing_strategy = strategy,
                     method = "velocity",
                     write_fits = F, verbose = F, cores = 2, 
                     mass_flag = F, moments = 2)
    
    write_simspin_FITS(output_file = paste0(cubes_loc,
                                            "p2p_EAGLE_unbinned_stellar_velocity_galaxyID_",
                                            galaxy_id, ".fits"),     # where should the FITS file be written to?
                       simspin_datacube = unbinned_output,                     # output from `build_datacube`
                       object_name = paste0("p2p_EAGLE_unbinned_stellar_velocity_galaxyID_", galaxy_id),
                       telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                       observer_name = "K E Harborne",
                       input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                       split_save=F,              # save cube and images to different files?
                       mask=NA,                      # add a "mask" image to the output FITS
                       galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
    
    cat(paste0("> Finished observing unbinned ", ob, " / ", batch_size, " <"),
        sep= "\n", file = log_file, append = T)
    
    cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
    unbinned_kinematics =
      tryCatch(expr = get_kinematic_properties(ifu_output = unbinned_output,
                                               photometry_properties = photometry_properties, 
                                               particle_limit = 0, plot = F, voronoi=F),
               error = function(e){NA})
    
    
    # Step 4: Observing the kinematics with voronoi binning
    
    vbinned_output = 
      tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                     telescope = ifu, 
                                     observing_strategy = strategy,
                                     method = "velocity",
                                     write_fits = F,
                                     verbose = F, cores = 2, 
                                     mass_flag = F, voronoi_bin = T, 
                                     vorbin_limit = 200, moments = 2),
               error = function(e){NA})
    
    if (any(is.na(vbinned_output))){ 
      
      cat(paste0("ERROR observing voronoi 200 ", ob, "."), 
          sep= "\n", file = log_file, append = T)
      
      binned_kinematics = NA
      
    } else {
      
      write_simspin_FITS(output_file = paste0(cubes_loc,
                                              "p2p_EAGLE_vbinned_200_stellar_velocity_galaxyID_",
                                              galaxy_id, ".fits"),     # where should the FITS file be written to?
                         simspin_datacube = vbinned_output,                     # output from `build_datacube`
                         object_name = paste0("p2p_EAGLE_vbinned_200_stellar_velocity_galaxyID_", galaxy_id),
                         telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                         observer_name = "K E Harborne",
                         input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                         split_save=F,              # save cube and images to different files?
                         mask=NA,                      # add a "mask" image to the output FITS
                         galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
      
      cat(paste0("> Finished observing voronoi-binned ", ob, " / ", batch_size, " <"),
          sep= "\n", file = log_file, append = T)
      
      cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
      binned_kinematics =
        tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                 photometry_properties = photometry_properties, 
                                                 particle_limit = 0, plot = F, voronoi=T),
                 error = function(e){NA})
    }
    
    if (all(is.na(unbinned_kinematics)) | all(is.na(binned_kinematics))){
      
      cat(paste0("ERROR observing ", ob, "."), 
          sep= "\n", file = log_file, append = T)
      
    } else {
      
      # Record observation particulars
      cat(paste0(c(galaxy_id, 
                   vbinned_output$observation$inc_deg,
                   vbinned_output$observation$twist_deg,
                   vbinned_output$observation$z,
                   vbinned_output$observation$psf_fwhm,
                   photometry_properties$major_axis_a_kpc,
                   photometry_properties$minor_axis_b_kpc,
                   photometry_properties$primary_axis,
                   photometry_properties$ellipticity,
                   photometry_properties$sersic_index), collapse=","),
          sep='\n', file=obs_prop_file, append=T)
      
      # Record unbinned kinematics
      cat(paste0(c(galaxy_id,
                   unbinned_kinematics$lambda_re,
                   unbinned_kinematics$vsigma,
                   unbinned_kinematics$j,
                   unbinned_kinematics$sigma_m,
                   unbinned_kinematics$rotational_velocity,
                   unbinned_kinematics$unique_pixels), collapse=","),
          sep='\n', file=ub_kin_prop_file, append=T)
      
      # Record binned kinematics
      cat(paste0(c(galaxy_id,
                   200,
                   binned_kinematics$lambda_re,
                   binned_kinematics$vsigma,
                   binned_kinematics$j,
                   binned_kinematics$sigma_m,
                   binned_kinematics$rotational_velocity,
                   binned_kinematics$unique_pixels), collapse=","),
          sep='\n', file=vb_kin_prop_file, append=T)
      
      # Step 5 ----------------------------------------------------------------
      # Now with increased voronoi binning number
      
      vbinned_output = 
        tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                       telescope = ifu, 
                                       observing_strategy = strategy,
                                       method = "velocity",
                                       write_fits = F, 
                                       verbose = F, cores = 2, 
                                       mass_flag = F, voronoi_bin = T, 
                                       vorbin_limit = 300, moments = 2),
                 error = function(e){NA})
      
      if (any(is.na(vbinned_output))){ 
        
        cat(paste0("ERROR observing voronoi 300 ", ob, "."), 
            sep= "\n", file = log_file, append = T)
        
      } else {
        
        write_simspin_FITS(output_file = paste0(cubes_loc,
                                                "p2p_EAGLE_vbinned_300_stellar_velocity_galaxyID_",
                                                galaxy_id, ".fits"),     # where should the FITS file be written to?
                           simspin_datacube = vbinned_output,                     # output from `build_datacube`
                           object_name = paste0("p2p_EAGLE_vbinned_300_stellar_velocity_galaxyID_", galaxy_id),
                           telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                           observer_name = "K E Harborne",
                           input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                           split_save=F,              # save cube and images to different files?
                           mask=NA,                      # add a "mask" image to the output FITS
                           galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
        
        cat(paste0("> Finished observing 300 voronoi-binned ", ob, " / ", batch_size, " <"),
            sep= "\n", file = log_file, append = T)
        
        cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
        binned_kinematics =
          tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                   photometry_properties = photometry_properties, 
                                                   particle_limit = 0, plot = F, voronoi=T),
                   error = function(e){NA})
        
        # Record binned kinematics
        cat(paste0(c(galaxy_id,
                     300,
                     binned_kinematics$lambda_re,
                     binned_kinematics$vsigma,
                     binned_kinematics$j,
                     binned_kinematics$sigma_m,
                     binned_kinematics$rotational_velocity,
                     binned_kinematics$unique_pixels), collapse=","),
            sep='\n', file=vb_kin_prop_file, append=T)
      }
      
      
      vbinned_output = 
        tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                       telescope = ifu, 
                                       observing_strategy = strategy,
                                       method = "velocity",
                                       write_fits = F, 
                                       verbose = F, cores = 2, 
                                       mass_flag = F, voronoi_bin = T, 
                                       vorbin_limit = 400, moments = 2),
                 error = function(e){NA})
      
      if (any(is.na(vbinned_output))){ 
        
        cat(paste0("ERROR observing voronoi 400 ", ob, "."), 
            sep= "\n", file = log_file, append = T)
        
      } else {
        
        write_simspin_FITS(output_file = paste0(cubes_loc,
                                                "p2p_EAGLE_vbinned_400_stellar_velocity_galaxyID_",
                                                galaxy_id, ".fits"),     # where should the FITS file be written to?
                           simspin_datacube = vbinned_output,                     # output from `build_datacube`
                           object_name = paste0("p2p_EAGLE_vbinned_400_stellar_velocity_galaxyID_", galaxy_id),
                           telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                           observer_name = "K E Harborne",
                           input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                           split_save=F,              # save cube and images to different files?
                           mask=NA,                      # add a "mask" image to the output FITS
                           galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
        
        cat(paste0("> Finished observing 400 voronoi-binned ", ob, " / ", batch_size, " <"),
            sep= "\n", file = log_file, append = T)
        
        cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
        binned_kinematics =
          tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                   photometry_properties = photometry_properties, 
                                                   particle_limit = 0, plot = F, voronoi=T),
                   error = function(e){NA})
        
        # Record binned kinematics
        cat(paste0(c(galaxy_id,
                     400,
                     binned_kinematics$lambda_re,
                     binned_kinematics$vsigma,
                     binned_kinematics$j,
                     binned_kinematics$sigma_m,
                     binned_kinematics$rotational_velocity,
                     binned_kinematics$unique_pixels), collapse=","),
            sep='\n', file=vb_kin_prop_file, append=T)
      }
      
      
      vbinned_output = 
        tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                       telescope = ifu, 
                                       observing_strategy = strategy,
                                       method = "velocity",
                                       write_fits = F, 
                                       verbose = F, cores = 2, 
                                       mass_flag = F, voronoi_bin = T, 
                                       vorbin_limit = 500, moments = 2),
                 error = function(e){NA})
      
      if (any(is.na(vbinned_output))){ 
        
        cat(paste0("ERROR observing voronoi 500 ", ob, "."), 
            sep= "\n", file = log_file, append = T)
        
      } else {
        
        write_simspin_FITS(output_file = paste0(cubes_loc,
                                                "p2p_EAGLE_vbinned_500_stellar_velocity_galaxyID_",
                                                galaxy_id, ".fits"),     # where should the FITS file be written to?
                           simspin_datacube = vbinned_output,                     # output from `build_datacube`
                           object_name = paste0("p2p_EAGLE_vbinned_500_stellar_velocity_galaxyID_", galaxy_id),
                           telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                           observer_name = "K E Harborne",
                           input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                           split_save=F,              # save cube and images to different files?
                           mask=NA,                      # add a "mask" image to the output FITS
                           galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
        
        cat(paste0("> Finished observing 500 voronoi-binned ", ob, " / ", batch_size, " <"),
            sep= "\n", file = log_file, append = T)
        
        cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
        binned_kinematics =
          tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                   photometry_properties = photometry_properties, 
                                                   particle_limit = 0, plot = F, voronoi=T),
                   error = function(e){NA})
        
        # Record binned kinematics
        cat(paste0(c(galaxy_id,
                     500,
                     binned_kinematics$lambda_re,
                     binned_kinematics$vsigma,
                     binned_kinematics$j,
                     binned_kinematics$sigma_m,
                     binned_kinematics$rotational_velocity,
                     binned_kinematics$unique_pixels), collapse=","),
            sep='\n', file=vb_kin_prop_file, append=T)
        
      }
      
      vbinned_output = 
        tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                       telescope = ifu, 
                                       observing_strategy = strategy,
                                       method = "velocity",
                                       write_fits = F, 
                                       verbose = F, cores = 2, 
                                       mass_flag = F, voronoi_bin = T, 
                                       vorbin_limit = 600, moments = 2),
                 error = function(e){NA})
      
      if (any(is.na(vbinned_output))){ 
        
        cat(paste0("ERROR observing voronoi 600 ", ob, "."), 
            sep= "\n", file = log_file, append = T)
        
      } else {
        
        write_simspin_FITS(output_file = paste0(cubes_loc,
                                                "p2p_EAGLE_vbinned_600_stellar_velocity_galaxyID_",
                                                galaxy_id, ".fits"),     # where should the FITS file be written to?
                           simspin_datacube = vbinned_output,                     # output from `build_datacube`
                           object_name = paste0("p2p_EAGLE_vbinned_600_stellar_velocity_galaxyID_", galaxy_id),
                           telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                           observer_name = "K E Harborne",
                           input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                           split_save=F,              # save cube and images to different files?
                           mask=NA,                      # add a "mask" image to the output FITS
                           galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
        
        cat(paste0("> Finished observing 600 voronoi-binned ", ob, " / ", batch_size, " <"),
            sep= "\n", file = log_file, append = T)
        
        cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
        binned_kinematics =
          tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                   photometry_properties = photometry_properties, 
                                                   particle_limit = 0, plot = F, voronoi=T),
                   error = function(e){NA})
        
        # Record binned kinematics
        cat(paste0(c(galaxy_id,
                     600,
                     binned_kinematics$lambda_re,
                     binned_kinematics$vsigma,
                     binned_kinematics$j,
                     binned_kinematics$sigma_m,
                     binned_kinematics$rotational_velocity), collapse=","),
            sep='\n', file=vb_kin_prop_file, append=T)
      }
      
      vbinned_output = 
        tryCatch(expr = build_datacube(simspin_file = ss_galaxy_file, 
                                       telescope = ifu, 
                                       observing_strategy = strategy,
                                       method = "velocity",
                                       write_fits = F, 
                                       verbose = F, cores = 2, 
                                       mass_flag = F, voronoi_bin = T, 
                                       vorbin_limit = 700, moments = 2),
                 error = function(e){NA})
      
      if (any(is.na(vbinned_output))){ 
        
        cat(paste0("ERROR observing voronoi 700 ", ob, "."), 
            sep= "\n", file = log_file, append = T)
        
      } else {
        
        write_simspin_FITS(output_file = paste0(cubes_loc,
                                                "p2p_EAGLE_vbinned_700_stellar_velocity_galaxyID_",
                                                galaxy_id, ".fits"),     # where should the FITS file be written to?
                           simspin_datacube = vbinned_output,                     # output from `build_datacube`
                           object_name = paste0("p2p_EAGLE_vbinned_700_stellar_velocity_galaxyID_", galaxy_id),
                           telescope_name = "AAO", instrument_name = "SAMI",  # metadata information
                           observer_name = "K E Harborne",
                           input_simspin_file = paste0(ss_file_loc, ss_files[each]),
                           split_save=F,              # save cube and images to different files?
                           mask=NA,                      # add a "mask" image to the output FITS
                           galaxy_centre = c(0,0,0)) # centre to give mock RA and Dec in outputs
        
        cat(paste0("> Finished observing 700 voronoi-binned ", ob, " / ", batch_size, " <"),
            sep= "\n", file = log_file, append = T)
        
        cat("Measuring kinematics...", sep= "\n", file = log_file, append = T)
        binned_kinematics =
          tryCatch(expr = get_kinematic_properties(ifu_output = vbinned_output,
                                                   photometry_properties = photometry_properties, 
                                                   particle_limit = 0, plot = F, voronoi=T),
                   error = function(e){NA})
        
        # Record binned kinematics
        cat(paste0(c(galaxy_id,
                     700,
                     binned_kinematics$lambda_re,
                     binned_kinematics$vsigma,
                     binned_kinematics$j,
                     binned_kinematics$sigma_m,
                     binned_kinematics$rotational_velocity), collapse=","),
            sep='\n', file=vb_kin_prop_file, append=T)
      }
      
      cat(paste0(">>> Finished observing EAGLE galaxy ", galaxy_id, " ( ", ob, " / 100 ) <<<"),
          sep= "\n", file = log_file, append = T)
      cat(paste0(" "),
          sep= "\n", file = log_file, append = T)
      
    }
    
  }
  
  ob = ob + 1
}