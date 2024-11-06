# Kate Harborne 17/10/23
# ============================================================================ #
# Running bootstrap analysis: 
# Using the very high resolution model, we can examine the change in kinematics
# per pixel as the number of particles is bootstrap sampled. This is done for 
# the N-body bulge and disk models.
# Here we provide the code for plotting Figure 6.
# ============================================================================ #

library(SimSpin)
library(magicaxis)

# PATHS ========================================================================
loc = "simulations/hdf5/"
ss_loc = "simspin_files/"
sa_loc = "sim_analysis/"
data_loc = "data/"
plot_loc = "plots/"

ifu = telescope(type="MUSE", fov = 0.8)
sky_00 = observing_strategy(dist_kpc_per_arcsec = 0.75, 
                            inc_deg = 0, blur = F, 
                            pointing_kpc = c(1,0))
sky_90 = observing_strategy(dist_kpc_per_arcsec = 0.5, 
                            inc_deg = 90, blur = F, 
                            pointing_kpc = c(3,0))

# Assuming that the control models have already had the SimSpin files generated,
# the following line should fail:

disk_control = make_simspin_file(filename = paste0(loc, "res_disk_sm650_control.hdf5"),
                                 write_to_file = T, overwrite = T,
                                 cores = 3, centre = c(0,0,0),
                                 output = paste0(ss_loc, "disk_650_control_simspin_file.Rdata"))

bulge_control = make_simspin_file(filename = paste0(loc, "res_bulge_sm650_control.hdf5"),
                                  write_to_file = T, 
                                  cores = 3, centre = c(0,0,0),
                                  output = paste0(ss_loc, "bulge_650_control_simspin_file.Rdata"))

# FUNCTIONS ====================================================================

generate_files = function(simspin_file, telescope, observing_strategy, method){
  
  simspin_data = readRDS(simspin_file)
  observation = observation(telescope = telescope, observing_strategy = observing_strategy, method = method)
  
  galaxy_data = simspin_data$star_part
  
  # Twisting galaxy about the z-axis to look from an angle
  twisted_data = twist_galaxy(galaxy_data, twist_rad = observation$twist_rad)
  
  # Projecting the galaxy to given inclination
  obs_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  
  galaxy_data$x = (obs_data$x + observation$pointing_kpc[1])   # adjusting pointing of the aperture by x_kpc
  galaxy_data$y = obs_data$y                                   #   and y_kpc
  galaxy_data$z = (obs_data$z + observation$pointing_kpc[2])
  galaxy_data$vx = obs_data$vx; galaxy_data$vy = obs_data$vy; galaxy_data$vz = obs_data$vz
  
  remove(obs_data, twisted_data)
  
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z, breaks=observation$sbin_seq, labels=F)) - (observation$sbin)
  
  # Trimming particles that lie outside the aperture of the telescope
  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],]
  
  return(galaxy_data)
}

measure_kinematics = function(galaxy_data, telescope, observing_strategy, method){
  
  .meanwt = function(x,wt){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      val = sum(x*wt, na.rm=T)/sum(wt,na.rm=T)
    }
    return(val)
  } # weighted mean
  
  .varwt = function(x, wt, xcen){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      if (missing(xcen)){xcen = .meanwt(x,wt)}
      val = sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T)
    }
    return(val)
  } # weighted variance
  
  .sum_velocities = function(galaxy_sample, observation){
    vel_diff = function(lum, vy){diff((lum * pnorm(observation$vbin_edges, mean = vy, sd = 0)))}
    
    bins = mapply(vel_diff, galaxy_sample$luminosity, galaxy_sample$vy)
    
    return(rowSums(bins))
    
  }
  
  .losvd_fit_h3h4 = function(par, x, losvd){
    
    vel = par[1]
    sig = par[2]
    h3 = par[3]
    h4 = par[4]
    k  = 1
    
    w = (x - vel)/sig
    H3 = (1/sqrt(6))  * (((2*sqrt(2))* w^3) - ((3*sqrt(2)) * w))
    H4 = (1/sqrt(24)) * ((4* w^4) - (12 * w^2) + 3)
    
    measured_vlos = (k * exp(-0.5*(w^2))) * (1 + (h3*H3) + (h4*H4))
    return=sum((measured_vlos-losvd)^2)
  }
  
  .losvd_out_h3h4 = function(x, vel, sig, h3, h4){
    k=1
    w = (x - vel)/sig
    H3 = (1/sqrt(6))  * (((2*sqrt(2))* w^3) - ((3*sqrt(2)) * w))
    H4 = (1/sqrt(24)) * ((4* w^4) - (12 * w^2) + 3)
    
    measured_vlos = (k * exp(-0.5*(w^2))) * (1 + (h3*H3) + (h4*H4))
    return(measured_vlos)
  }
  
  
  .losvd_fit_vsig = function(par, x, losvd){
    
    vel = par[1]
    sig = par[2]
    k  = 1
    
    w = (x - vel)/sig
    measured_vlos = k * exp(-0.5*(w^2))
    return=sum((measured_vlos-losvd)^2)
  }
  
  .losvd_out_vsig = function(x, vel, sig, h3, h4){
    k=1
    w = (x - vel)/sig
    
    measured_vlos = k * exp(-0.5*(w^2)) 
    return(measured_vlos)
  }
  
  observation = SimSpin::observation(telescope = telescope, observing_strategy = observing_strategy, method = method)
  
  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]
  
  observation$vbin = 5*ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
  if (observation$vbin <= 2){observation$vbin = 3}
  
  observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
  observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2
  
  #vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = observation$sbin^2)
  part_map = array(data=0, dim = observation$sbin^2)
  vel_4_map  = array(data=0, dim = observation$sbin^2)
  disp_4_map = array(data=0, dim = observation$sbin^2)
  vel_2_map  = array(data=0, dim = observation$sbin^2)
  disp_2_map = array(data=0, dim = observation$sbin^2)
  h3_map   = array(data=0, dim = observation$sbin^2)
  h4_map   = array(data=0, dim = observation$sbin^2)
  res_4_map  = array(data=0, dim = observation$sbin^2)
  res_2_map  = array(data=0, dim = observation$sbin^2)
  pix_map  = array(data=0, dim = observation$sbin^2)
  
  for (i in 1:(dim(part_in_spaxel)[1])){
    
    num_part = part_in_spaxel$N[i] # number of particles in spaxel
    part_map[part_in_spaxel$pixel_pos[i]] = num_part
    pix_map[part_in_spaxel$pixel_pos[i]] = part_in_spaxel$pixel_pos[i]
    
    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    galaxy_sample[, luminosity := Mass, ]
    galaxy_sample[, band_lum := Mass, ]
    
    vel_spec = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)
    
    vel_ini = .meanwt(observation$vbin_seq, vel_spec)
    sd_ini  = sqrt(.varwt(observation$vbin_seq, vel_spec, vel_ini))
    
    kin  = tryCatch({stats::optim(par   = c(vel_ini,sd_ini,0,0),
                                  fn    = .losvd_fit_h3h4,
                                  x     = observation$vbin_seq,
                                  losvd = (vel_spec/(max(vel_spec, na.rm=T))),
                                  method="BFGS", control=list(reltol=1e-9))$par},
                    error = function(e){c(0,0,0,0)})
    
    vel_4_map[part_in_spaxel$pixel_pos[i]]  = kin[1]
    disp_4_map[part_in_spaxel$pixel_pos[i]] = kin[2]
    h3_map[part_in_spaxel$pixel_pos[i]]   = kin[3]
    h4_map[part_in_spaxel$pixel_pos[i]]   = kin[4]
    res_4_map[part_in_spaxel$pixel_pos[i]]  = mean(abs(.losvd_out_h3h4(x=observation$vbin_seq, vel=kin[1], sig=kin[2], 
                                                                       h3=kin[3], h4=kin[4]) -
                                                         vel_spec/(max(vel_spec, na.rm=T)) ), na.rm=T)
    
    kin2  = tryCatch({stats::optim(par   = c(vel_ini,sd_ini),
                                   fn    = .losvd_fit_vsig,
                                   x     = observation$vbin_seq,
                                   losvd = (vel_spec/(max(vel_spec, na.rm=T))),
                                   method="BFGS", control=list(reltol=1e-9))$par},
                     error = function(e){c(0,0)})
    
    vel_2_map[part_in_spaxel$pixel_pos[i]]  = kin2[1]
    disp_2_map[part_in_spaxel$pixel_pos[i]] = kin2[2]
    res_2_map[part_in_spaxel$pixel_pos[i]]  = mean(abs(.losvd_out_vsig(x=observation$vbin_seq, vel=kin2[1], sig=kin2[2]) -
                                                         vel_spec/(max(vel_spec, na.rm=T)) ), na.rm=T)
    
    
  }
  
  return(list(velocity_4_image = array(data = vel_4_map, dim = c(observation$sbin, observation$sbin)),
              dispersion_4_image = array(data = disp_4_map, dim = c(observation$sbin, observation$sbin)),
              velocity_2_image = array(data = vel_2_map, dim = c(observation$sbin, observation$sbin)),
              dispersion_2_image = array(data = disp_2_map, dim = c(observation$sbin, observation$sbin)),
              h3_image = array(data = h3_map, dim = c(observation$sbin, observation$sbin)), 
              h4_image = array(data = h4_map, dim = c(observation$sbin, observation$sbin)),
              residuals_4 = array(data = res_4_map, dim = c(observation$sbin, observation$sbin)),
              residuals_2 = array(data = res_2_map, dim = c(observation$sbin, observation$sbin)),
              particles = array(data = part_map, dim = c(observation$sbin, observation$sbin)),
              pixel_pos = array(data = pix_map, dim = c(observation$sbin, observation$sbin)))
  )
}

sample_galaxy = function(galaxy_data, ppp){
  
  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]
  new_galaxy_data = data.frame(array(data = 0.0, dim = c(dim(part_in_spaxel)[1]*ppp, 13)))
  names(new_galaxy_data) = names(galaxy_data)
  
  sample_seq = c(seq(1, nrow(new_galaxy_data), by=ppp), nrow(new_galaxy_data)+1)
  
  for (i in 1:(dim(part_in_spaxel)[1])){
    num_part = part_in_spaxel$N[i] # number of particles in spaxel
    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    new_galaxy_data[sample_seq[i]:(sample_seq[(i+1)]-1),] =
      galaxy_sample[sample(seq(1, num_part), size = ppp, replace = T),]
  }    
  
  return(data.table::as.data.table(new_galaxy_data))
}

# Bootstrap analysis (DISK at 90 deg) ==========================================

disk_base = build_datacube(simspin_file = paste0(ss_loc,"disk_650_control_simspin_file.Rdata"),
                           telescope = ifu, 
                           observing_strategy = sky_90, 
                           method = "velocity")

ppp = sort(seq(50, 500, by=10), decreasing = T)
samples_per_ppp = 500
number_of_obs = (samples_per_ppp*length(ppp)*16) + 16

bootstrap_disk =
  data.frame("number_of_ppp" = numeric(number_of_obs),
             "pixel_id" = numeric(number_of_obs),
             "velocity_4" = numeric(number_of_obs),
             "dispersion_4" = numeric(number_of_obs),
             "h3" = numeric(number_of_obs),
             "h4" = numeric(number_of_obs),
             "velocity_2" = numeric(number_of_obs),
             "dispersion_2" = numeric(number_of_obs))

bootstrap_seq = seq(1, number_of_obs, by=16)

disk_table = generate_files(simspin_file = paste0(ss_loc,"disk_650_control_simspin_file.Rdata"),
                            telescope = ifu, 
                            observing_strategy = sky_90, 
                            method = "velocity")

disk_obs = measure_kinematics(galaxy_data = disk_table, 
                              telescope = ifu, 
                              observing_strategy = sky_90, 
                              method = "velocity")

i = 1
bootstrap_disk$number_of_ppp[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$particles)
bootstrap_disk$pixel_id[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$pixel_pos)
bootstrap_disk$velocity_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$velocity_4_image)
bootstrap_disk$dispersion_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$dispersion_4_image)
bootstrap_disk$h3[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$h3_image)
bootstrap_disk$h4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$h4_image)
bootstrap_disk$velocity_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$velocity_2_image)
bootstrap_disk$dispersion_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$dispersion_2_image)

i = 2
for (p in ppp){
  for (sample in 1:samples_per_ppp){
    
    sampled_galaxy_data = sample_galaxy(disk_table, ppp = p)
    
    disk_obs = measure_kinematics(galaxy_data = sampled_galaxy_data, 
                                  telescope = ifu, 
                                  observing_strategy = sky_90, 
                                  method = "velocity")
    
    bootstrap_disk$number_of_ppp[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$particles)
    bootstrap_disk$pixel_id[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$pixel_pos)
    bootstrap_disk$velocity_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$velocity_4_image)
    bootstrap_disk$dispersion_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$dispersion_4_image)
    bootstrap_disk$h3[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$h3_image)
    bootstrap_disk$h4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$h4_image)
    bootstrap_disk$velocity_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$velocity_2_image)
    bootstrap_disk$dispersion_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(disk_obs$dispersion_2_image)
    i = i+1
    print(paste("Done with ", i, " of ", number_of_obs))
  }
}

saveRDS(bootstrap_disk, paste0(data_loc, "bootstrap_disk_inc90_sm650.Rdata"))

# Computing properties of averages and interquartile ranges for each bootstrap ====

no = length(rep(sort(seq(50, 500, by=10), decreasing = T), 16))

disk_velocity = data.frame("pixel_id" = rep(seq(1,16), each=46),
                           "ppp" = rep(sort(seq(50, 400, by=10), decreasing = T), 16),
                           "mean_v4" = numeric(no),
                           "p25_v4" = numeric(no),
                           "p75_v4" = numeric(no),
                           "mean_d4" = numeric(no),
                           "p25_d4" = numeric(no),
                           "p75_d4" = numeric(no),
                           "mean_v2" = numeric(no),
                           "p25_v2" = numeric(no),
                           "p75_v2" = numeric(no),
                           "mean_d2" = numeric(no),
                           "p25_d2" = numeric(no),
                           "p75_d2" = numeric(no),
                           "mean_3" = numeric(no),
                           "p25_3" = numeric(no),
                           "p75_3" = numeric(no),
                           "mean_4" = numeric(no),
                           "p25_4" = numeric(no),
                           "p75_4" = numeric(no))

for (group in 1:no){
  
  pixel_id = disk_velocity$pixel_id[group]
  ppp = disk_velocity$ppp[group]
  
  disk_group = bootstrap_disk[bootstrap_disk$number_of_ppp == ppp & bootstrap_disk$pixel_id == pixel_id,]
  
  disk_velocity$mean_v4[group] = mean(disk_group$velocity_4)
  disk_velocity$p25_v4[group] = quantile(disk_group$velocity_4, c(0.25))
  disk_velocity$p75_v4[group] = quantile(disk_group$velocity_4, c(0.75))
  
  disk_velocity$mean_d4[group] = mean(disk_group$dispersion_4)
  disk_velocity$p25_d4[group] = quantile(disk_group$dispersion_4, c(0.25))
  disk_velocity$p75_d4[group] = quantile(disk_group$dispersion_4, c(0.75))
  
  disk_velocity$mean_v2[group] = mean(disk_group$velocity_2)
  disk_velocity$p25_v2[group] = quantile(disk_group$velocity_2, c(0.25))
  disk_velocity$p75_v2[group] = quantile(disk_group$velocity_2, c(0.75))
  
  disk_velocity$mean_d2[group] = mean(disk_group$dispersion_2)
  disk_velocity$p25_d2[group] = quantile(disk_group$dispersion_2, c(0.25))
  disk_velocity$p75_d2[group] = quantile(disk_group$dispersion_2, c(0.75))
  
  disk_velocity$mean_3[group] = mean(disk_group$h3)
  disk_velocity$p25_3[group] = quantile(disk_group$h3, c(0.25))
  disk_velocity$p75_3[group] = quantile(disk_group$h3, c(0.75))
  
  disk_velocity$mean_4[group] = mean(disk_group$h4)
  disk_velocity$p25_4[group] = quantile(disk_group$h4, c(0.25))
  disk_velocity$p75_4[group] = quantile(disk_group$h4, c(0.75))
  
}

saveRDS(disk_velocity, paste0(data_loc,"/bootstrap_avgs_disk_inc90_sm650.Rdata"))

# Bootstrap analysis (BULGE at 0 deg) ==========================================
bulge_base = build_datacube(simspin_file = paste0(ss_loc,"bulge_650_control_simspin_file.Rdata"),
                            telescope = ifu, 
                            observing_strategy = sky_00, 
                            method = "velocity")

ppp = sort(seq(50, 500, by=10), decreasing = T)
samples_per_ppp = 500
number_of_obs = (samples_per_ppp*length(ppp)*16) + 16
bootstrap_seq = seq(1, number_of_obs, by=16)

bootstrap_bulge =
  data.frame("number_of_ppp" = numeric(number_of_obs),
             "pixel_id" = numeric(number_of_obs),
             "velocity_4" = numeric(number_of_obs),
             "dispersion_4" = numeric(number_of_obs),
             "h3" = numeric(number_of_obs),
             "h4" = numeric(number_of_obs),
             "velocity_2" = numeric(number_of_obs),
             "dispersion_2" = numeric(number_of_obs))

bulge_table = generate_files(simspin_file = "project_resolution/particles_to_pixels_paper/simspin_files/sm650/bulge_650_control_simspin_file.Rdata",
                             telescope = ifu, observing_strategy = sky_00, method = "velocity")

bulge_obs = measure_kinematics(galaxy_data = bulge_table, 
                               telescope = ifu, 
                               observing_strategy = sky_00, 
                               method = "velocity")

i = 1
bootstrap_bulge$number_of_ppp[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$particles)
bootstrap_bulge$pixel_id[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$pixel_pos)
bootstrap_bulge$velocity_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$velocity_4_image)
bootstrap_bulge$dispersion_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$dispersion_4_image)
bootstrap_bulge$h3[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$h3_image)
bootstrap_bulge$h4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$h4_image)
bootstrap_bulge$velocity_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$velocity_2_image)
bootstrap_bulge$dispersion_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$dispersion_2_image)

i = 2
for (p in ppp){
  for (sample in 1:samples_per_ppp){
    
    sampled_galaxy_data = sample_galaxy(bulge_table, ppp = p)
    
    bulge_obs = measure_kinematics(galaxy_data = sampled_galaxy_data, 
                                   telescope = ifu, 
                                   observing_strategy = sky_00, 
                                   method = "velocity")
    
    bootstrap_bulge$number_of_ppp[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$particles)
    bootstrap_bulge$pixel_id[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$pixel_pos)
    bootstrap_bulge$velocity_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$velocity_4_image)
    bootstrap_bulge$dispersion_4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$dispersion_4_image)
    bootstrap_bulge$h3[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$h3_image)
    bootstrap_bulge$h4[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$h4_image)
    bootstrap_bulge$velocity_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$velocity_2_image)
    bootstrap_bulge$dispersion_2[bootstrap_seq[i]:(bootstrap_seq[i]+15)] = c(bulge_obs$dispersion_2_image)
    
    i = i+1
    print(paste("Done with ", i, " of ", number_of_obs))
  }
}

saveRDS(bootstrap_bulge, paste0(data_loc,"bootstrap_bulge_inc00_sm650.Rdata"))

# Computing properties of averages and interquartile ranges for each bootstrap ====

no = length(rep(sort(seq(50, 500, by=10), decreasing = T), 16))

bulge_velocity = data.frame("pixel_id" = rep(seq(1,16), each=46),
                            "ppp" = rep(sort(seq(50, 500, by=10), decreasing = T), 16),
                            "mean_v4" = numeric(no),
                            "p25_v4" = numeric(no),
                            "p75_v4" = numeric(no),
                            "mean_d4" = numeric(no),
                            "p25_d4" = numeric(no),
                            "p75_d4" = numeric(no),
                            "mean_v2" = numeric(no),
                            "p25_v2" = numeric(no),
                            "p75_v2" = numeric(no),
                            "mean_d2" = numeric(no),
                            "p25_d2" = numeric(no),
                            "p75_d2" = numeric(no),
                            "mean_3" = numeric(no),
                            "p25_3" = numeric(no),
                            "p75_3" = numeric(no),
                            "mean_4" = numeric(no),
                            "p25_4" = numeric(no),
                            "p75_4" = numeric(no))


for (group in 1:no){
  
  pixel_id = bulge_velocity$pixel_id[group]
  ppp = bulge_velocity$ppp[group]
  
  bulge_group = bootstrap_bulge[bootstrap_bulge$number_of_ppp == ppp & bootstrap_bulge$pixel_id == pixel_id,]
  
  bulge_velocity$mean_v4[group] = mean(bulge_group$velocity_4)
  bulge_velocity$p25_v4[group] = quantile(bulge_group$velocity_4, c(0.25))
  bulge_velocity$p75_v4[group] = quantile(bulge_group$velocity_4, c(0.75))
  
  bulge_velocity$mean_d4[group] = mean(bulge_group$dispersion_4)
  bulge_velocity$p25_d4[group] = quantile(bulge_group$dispersion_4, c(0.25))
  bulge_velocity$p75_d4[group] = quantile(bulge_group$dispersion_4, c(0.75))
  
  bulge_velocity$mean_v2[group] = mean(bulge_group$velocity_2)
  bulge_velocity$p25_v2[group] = quantile(bulge_group$velocity_2, c(0.25))
  bulge_velocity$p75_v2[group] = quantile(bulge_group$velocity_2, c(0.75))
  
  bulge_velocity$mean_d2[group] = mean(bulge_group$dispersion_2)
  bulge_velocity$p25_d2[group] = quantile(bulge_group$dispersion_2, c(0.25))
  bulge_velocity$p75_d2[group] = quantile(bulge_group$dispersion_2, c(0.75))
  
  bulge_velocity$mean_3[group] = mean(bulge_group$h3)
  bulge_velocity$p25_3[group] = quantile(bulge_group$h3, c(0.25))
  bulge_velocity$p75_3[group] = quantile(bulge_group$h3, c(0.75))
  
  bulge_velocity$mean_4[group] = mean(bulge_group$h4)
  bulge_velocity$p25_4[group] = quantile(bulge_group$h4, c(0.25))
  bulge_velocity$p75_4[group] = quantile(bulge_group$h4, c(0.75))
  
}

saveRDS(bulge_velocity, paste0(data_loc,"bootstrap_avgs_bulge_inc00_sm650.Rdata"))

# Then getting data to provide the plot in Figure 6 ============================
# New reiteration of functions 

generate_files = function(simspin_file, telescope, observing_strategy, method){
  
  simspin_data = readRDS(simspin_file)
  observation = observation(telescope = telescope, observing_strategy = observing_strategy, method = method)
  
  galaxy_data = simspin_data$star_part
  
  # Twisting galaxy about the z-axis to look from an angle
  twisted_data = twist_galaxy(galaxy_data, twist_rad = observation$twist_rad)
  
  # Projecting the galaxy to given inclination
  obs_data = obs_galaxy(part_data = twisted_data, inc_rad = observation$inc_rad)
  
  galaxy_data$x = (obs_data$x + observation$pointing_kpc[1])   # adjusting pointing of the aperture by x_kpc
  galaxy_data$y = obs_data$y                                   #   and y_kpc
  galaxy_data$z = (obs_data$z + observation$pointing_kpc[2])
  galaxy_data$vx = obs_data$vx; galaxy_data$vy = obs_data$vy; galaxy_data$vz = obs_data$vz
  
  remove(obs_data, twisted_data)
  
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z, breaks=observation$sbin_seq, labels=F)) - (observation$sbin)
  
  # Trimming particles that lie outside the aperture of the telescope
  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$pixel_region[!is.na(observation$pixel_region)],]
  
  return(galaxy_data)
}

measure_kinematics = function(galaxy_data, telescope, observing_strategy, method){
  
  .meanwt = function(x,wt){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      val = sum(x*wt, na.rm=T)/sum(wt,na.rm=T)
    }
    return(val)
  } # weighted mean
  
  .varwt = function(x, wt, xcen){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      if (missing(xcen)){xcen = .meanwt(x,wt)}
      val = sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T)
    }
    return(val)
  } # weighted variance
  
  .sum_velocities = function(galaxy_sample, observation){
    vel_diff = function(lum, vy){diff((lum * pnorm(observation$vbin_edges, mean = vy, sd = 0)))}
    
    bins = mapply(vel_diff, galaxy_sample$luminosity, galaxy_sample$vy)
    
    return(rowSums(bins))
    
  }
  
  observation = SimSpin::observation(telescope = telescope, observing_strategy = observing_strategy, method = method)
  
  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]
  
  observation$vbin = 5*ceiling((max(abs(galaxy_data$vy))*2) / observation$vbin_size) # the number of velocity bins in the cube
  if (observation$vbin <= 2){observation$vbin = 3}
  
  observation$vbin_edges = seq(-(observation$vbin * observation$vbin_size)/2, (observation$vbin * observation$vbin_size)/2, by=observation$vbin_size)
  observation$vbin_seq   = observation$vbin_edges[1:observation$vbin] + diff(observation$vbin_edges)/2
  
  vel_spec = matrix(data = 0, ncol = observation$vbin, nrow = (dim(part_in_spaxel)[1]))
  part_map = array(data=0, dim = observation$sbin^2)
  pix_map  = array(data=0, dim = observation$sbin^2)
  
  for (i in 1:(dim(part_in_spaxel)[1])){
    
    num_part = part_in_spaxel$N[i] # number of particles in spaxel
    part_map[part_in_spaxel$pixel_pos[i]] = num_part
    pix_map[part_in_spaxel$pixel_pos[i]] = part_in_spaxel$pixel_pos[i]
    
    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    galaxy_sample[, luminosity := Mass, ]
    galaxy_sample[, band_lum := Mass, ]
    
    vel_spec[i,] = .sum_velocities(galaxy_sample = galaxy_sample, observation = observation)
    
  }
  
  return(list("particles"=array(data = part_map, dim = c(observation$sbin, observation$sbin)), 
              "pixel_pos"=array(data = pix_map, dim = c(observation$sbin, observation$sbin)), 
              "vel_spec"=vel_spec, "vel_seq"=observation$vbin_seq))
  
}

sample_galaxy = function(galaxy_data, ppp){
  
  # which particles sit in each spaxel?
  part_in_spaxel = galaxy_data[, list(val=list(ID), .N), by = "pixel_pos"]
  new_galaxy_data = data.frame(array(data = 0.0, dim = c(dim(part_in_spaxel)[1]*ppp, 13)))
  names(new_galaxy_data) = names(galaxy_data)
  
  sample_seq = c(seq(1, nrow(new_galaxy_data), by=ppp), nrow(new_galaxy_data)+1)
  
  for (i in 1:(dim(part_in_spaxel)[1])){
    num_part = part_in_spaxel$N[i] # number of particles in spaxel
    galaxy_sample = galaxy_data[ID %in% part_in_spaxel$val[[i]]]
    new_galaxy_data[sample_seq[i]:(sample_seq[(i+1)]-1),] =
      galaxy_sample[sample(seq(1, num_part), size = ppp, replace = T),]
  }    
  
  return(data.table::as.data.table(new_galaxy_data))
}

# Bootstrap analysis (DISK at 90 deg) ==========================================

disk_base = build_datacube(simspin_file = paste0(ss_loc,"disk_650_control_simspin_file.Rdata"),
                           telescope = ifu, 
                           observing_strategy = sky_90, 
                           method = "velocity")

ppp = c(500, 200, 50)#sort(seq(50, 500, by=10), decreasing = T)
samples_per_ppp = 500
number_of_obs = (samples_per_ppp*length(ppp))+1

bootstrap_disk =
  list("number_of_ppp" = numeric(number_of_obs),
       "vel_seq" = vector(mode = "list", length = number_of_obs),
       "vel_spec" = vector(mode = "list", length = number_of_obs))

bootstrap_seq = seq(1, number_of_obs)

disk_table = generate_files(simspin_file = paste0(ss_loc,"disk_650_control_simspin_file.Rdata"),
                            telescope = ifu, 
                            observing_strategy = sky_90, 
                            method = "velocity")

disk_obs = measure_kinematics(galaxy_data = disk_table, 
                              telescope = ifu, 
                              observing_strategy = sky_90, 
                              method = "velocity")

i = 1
bootstrap_disk$number_of_ppp[bootstrap_seq[i]] = disk_obs$particles
bootstrap_disk$pixel_id[bootstrap_seq[i]] = disk_obs$pixel_pos
bootstrap_disk$vel_spec[[bootstrap_seq[i]]] = disk_obs$vel_spec
bootstrap_disk$vel_seq[[bootstrap_seq[i]]] = disk_obs$vel_seq

i = 2
for (p in ppp){
  for (sample in 1:samples_per_ppp){
    
    sampled_galaxy_data = sample_galaxy(disk_table, ppp = p)
    
    disk_obs = measure_kinematics(galaxy_data = sampled_galaxy_data, 
                                  telescope = ifu, 
                                  observing_strategy = sky_90, 
                                  method = "velocity")
    
    bootstrap_disk$number_of_ppp[bootstrap_seq[i]] = disk_obs$particles
    bootstrap_disk$pixel_id[bootstrap_seq[i]] = disk_obs$pixel_pos
    bootstrap_disk$vel_spec[[bootstrap_seq[i]]] = disk_obs$vel_spec
    bootstrap_disk$vel_seq[[bootstrap_seq[i]]] = disk_obs$vel_seq
    print(paste("Done with ", i, " of ", number_of_obs))
    i = i+1
    
  }
}

saveRDS(bootstrap_disk, paste0(data_loc, "bootstrap_velspec_disk_inc90_sm650_forplot.Rdata"))

# Bootstrap analysis (BULGE at 0 deg) ==========================================
bulge_base = build_datacube(simspin_file = paste0(ss_loc,"bulge_650_control_simspin_file.Rdata"),
                            telescope = ifu, 
                            observing_strategy = sky_00, 
                            method = "velocity")

ppp = c(500, 200, 50)#sort(seq(50, 500, by=10), decreasing = T)
samples_per_ppp = 500
number_of_obs = (samples_per_ppp*length(ppp))+1

bootstrap_bulge =
  list("number_of_ppp" = numeric(number_of_obs),
       "vel_seq" = vector(mode = "list", length = number_of_obs),
       "vel_spec" = vector(mode = "list", length = number_of_obs))

bootstrap_seq = seq(1, number_of_obs)

bulge_table = generate_files(simspin_file = paste0(ss_loc,"bulge_650_control_simspin_file.Rdata"),
                             telescope = ifu, 
                             observing_strategy = sky_90, 
                             method = "velocity")

bulge_obs = measure_kinematics(galaxy_data = bulge_table, 
                               telescope = ifu, 
                               observing_strategy = sky_90, 
                               method = "velocity")

i = 1
bootstrap_bulge$number_of_ppp[bootstrap_seq[i]] = bulge_obs$particles
bootstrap_bulge$pixel_id[bootstrap_seq[i]] = bulge_obs$pixel_pos
bootstrap_bulge$vel_spec[[bootstrap_seq[i]]] = bulge_obs$vel_spec
bootstrap_bulge$vel_seq[[bootstrap_seq[i]]] = bulge_obs$vel_seq

i = 2
for (p in ppp){
  for (sample in 1:samples_per_ppp){
    
    sampled_galaxy_data = sample_galaxy(bulge_table, ppp = p)
    
    bulge_obs = measure_kinematics(galaxy_data = sampled_galaxy_data, 
                                   telescope = ifu, 
                                   observing_strategy = sky_90, 
                                   method = "velocity")
    
    bootstrap_bulge$number_of_ppp[bootstrap_seq[i]] = bulge_obs$particles
    bootstrap_bulge$pixel_id[bootstrap_seq[i]] = bulge_obs$pixel_pos
    bootstrap_bulge$vel_spec[[bootstrap_seq[i]]] = bulge_obs$vel_spec
    bootstrap_bulge$vel_seq[[bootstrap_seq[i]]] = bulge_obs$vel_seq
    print(paste("Done with ", i, " of ", number_of_obs))
    i = i+1
    
  }
}

saveRDS(bootstrap_bulge, paste0(data_loc, "bootstrap_velspec_bulge_inc90_sm650_forplot.Rdata"))


# Plotting sampled LOSVDs ======================================================
cex_val = 1.1

w = plot_array(0.40, 3)

bootstrap_disk = readRDS(paste0(data_loc, "bootstrap_velspec_disk_inc90_sm650_forplot.Rdata"))
bootstrap_bulge = readRDS(paste0(data_loc, "bootstrap_velspec_bulge_inc90_sm650_forplot.Rdata"))

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")
tol_pal = rep(tol_pal, length.out=1501)

plot_losvd = function(vel_seq, vel_spec, col, lwd){
  
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
  
  .meanwt = function(x,wt){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      val = sum(x*wt, na.rm=T)/sum(wt,na.rm=T)
    }
    return(val)
  } # weighted mean
  
  .varwt = function(x, wt, xcen){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      if (missing(xcen)){xcen = .meanwt(x,wt)}
      val = sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T)
    }
    return(val)
  } # weighted variance
  
  vel_ini = .meanwt(vel_seq, vel_spec)
  sd_ini  = sqrt(.varwt(vel_seq, vel_spec, vel_ini))
  
  kin_fit = stats::optim(par   = c(vel_ini,sd_ini),
                         fn    = .losvd_fit,
                         x     = vel_seq,
                         losvd = (vel_spec/(max(vel_spec, na.rm=T))),
                         method="BFGS", control=list(reltol=1e-9))$par
  
  v = seq(min(vel_seq), max(vel_seq))
  
  lines(v, .losvd_out(v, vel=kin_fit[1], sig = kin_fit[2]), lwd= lwd, col = col)
}

plot_residual = function(vel_seq, vel_spec, col, lwd, control_vseq, control_losvd){
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
  
  .meanwt = function(x,wt){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      val = sum(x*wt, na.rm=T)/sum(wt,na.rm=T)
    }
    return(val)
  } # weighted mean
  
  .varwt = function(x, wt, xcen){
    if (sum(wt,na.rm=T) == 0){
      val = 0
    } else {
      if (missing(xcen)){xcen = .meanwt(x,wt)}
      val = sum(wt*(x - xcen)^2, na.rm=T)/sum(wt, na.rm=T)
    }
    return(val)
  } # weighted variance
  
  vel_ini = .meanwt(vel_seq, vel_spec)
  sd_ini  = sqrt(.varwt(vel_seq, vel_spec, vel_ini))
  
  kin_fit = stats::optim(par   = c(vel_ini,sd_ini),
                         fn    = .losvd_fit,
                         x     = vel_seq,
                         losvd = (vel_spec/(max(vel_spec, na.rm=T))),
                         method="BFGS", control=list(reltol=1e-9))$par
  
  vel_ini_cont = .meanwt(control_vseq, control_losvd)
  sd_ini_cont  = sqrt(.varwt(control_vseq, control_losvd, vel_ini_cont))
  
  kin_fit_cont = stats::optim(par   = c(vel_ini_cont,sd_ini_cont),
                              fn    = .losvd_fit,
                              x     = control_vseq,
                              losvd = (control_losvd/(max(control_losvd, na.rm=T))),
                              method="BFGS", control=list(reltol=1e-9))$par
  
  v = seq(min(vel_seq), max(vel_seq))
  
  lines(v, 
        .losvd_out(v, vel=kin_fit[1], sig = kin_fit[2])-.losvd_out(v, vel=kin_fit_cont[1], sig = kin_fit_cont[2]),
        lwd = lwd, col = col)
  
}

# Plot of LOSVD ----------------------------------------------------------------
pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "LOSVD_bootstraps.pdf"), 
          width = 11, height = 9)

# Bulge ========================================================================
par(fig = c(w[1,1], w[1,2], 0.58, 1))
magplot(0, type="n",
        xlim = c(-600,600), 
        xaxt="n", ylim = c(0,1.2), yaxt="n",
        cex.lab = cex_val)
magaxis(side=1, labels = F)
magaxis(side=2, ylab="Normalised LOSVD", majorn=3, cex.lab = cex_val)
text(325,1.1, "500 particles/pixel", cex = cex_val, font=1)
text(-380,1.1, "Bulge Model", cex = cex_val, font=2)
for (i in which(bootstrap_bulge$number_of_ppp  == 500)){
  plot_losvd(bootstrap_bulge$vel_seq[[i]], bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3)
legend("topleft", legend = c("control"), lwd=3, col = tol_pal[1], bty="n", inset = c(0.01, 0.15))


par(fig = c(w[1,1], w[1,2], 0.43, 0.78), new=T)
magplot(0, type="n",
        xlim = c(-600,600), 
        xaxt="n", yaxt="n",
        cex.lab = cex_val)
magaxis(side=2, ylab=expression('\u0394'["X - control"]), majorn=3, cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
for (i in which(bootstrap_bulge$number_of_ppp  == 500)){
  plot_residual(bootstrap_bulge$vel_seq[[i]], 
                bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_bulge$vel_spec[[1]],
                control_vseq = bootstrap_bulge$vel_seq[[1]])
}
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])

par(fig = c(w[2,1], w[2,2], 0.58, 1), new=T)
magplot(0, type="n",
        xlim = c(-600,600), 
        xaxt="n", yaxt="n",  ylim = c(0,1.2),
        cex.lab = cex_val)
magaxis(side=1, labels = F)
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3)
text(325,1.1, "200 particles/pixel", cex = cex_val, font=1)
for (i in which(bootstrap_bulge$number_of_ppp  == 200)){
  plot_losvd(bootstrap_bulge$vel_seq[[i]], bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3)

par(fig = c(w[2,1], w[2,2], 0.43, 0.78), new=T)
magplot(0, type="n",
        xlim = c(-600,600), 
        xaxt="n",
        yaxt="n", cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
for (i in which(bootstrap_bulge$number_of_ppp  == 200)){
  plot_residual(bootstrap_bulge$vel_seq[[i]], 
                bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_bulge$vel_spec[[1]],
                control_vseq = bootstrap_bulge$vel_seq[[1]])
}
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])


par(fig = c(w[3,1], w[3,2], 0.58, 1), new=T)
magplot(0, type="n",
        xlim = c(-600,600),  ylim = c(0,1.2),
        xaxt="n", yaxt="n",
        cex.lab = cex_val)
magaxis(side=1, labels = F)
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3)
text(350,1.1, "50 particles/pixel", cex = cex_val, font=1)
for (i in which(bootstrap_bulge$number_of_ppp  == 50)){
  plot_losvd(bootstrap_bulge$vel_seq[[i]], bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_bulge$vel_seq[[1]], bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3)

par(fig = c(w[3,1], w[3,2], 0.43, 0.78), new=T)
magplot(0, type="n",
        xlim = c(-600,600), 
        xaxt="n",
        yaxt="n", cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
for (i in which(bootstrap_bulge$number_of_ppp  == 50)){
  plot_residual(bootstrap_bulge$vel_seq[[i]], 
                bootstrap_bulge$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_bulge$vel_spec[[1]],
                control_vseq = bootstrap_bulge$vel_seq[[1]])
}
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])
plot_residual(bootstrap_bulge$vel_seq[[1]], 
              bootstrap_bulge$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_bulge$vel_spec[[1]],
              control_vseq = bootstrap_bulge$vel_seq[[1]])

# Disk =========================================================================
par(fig = c(w[1,1], w[1,2], 0.15, 0.57), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        xaxt="n",  ylim = c(0,1.1), yaxt="n",
        cex.lab = cex_val)
magaxis(side=1, labels = F)
magaxis(side=2, ylab="Normalised LOSVD", majorn=3, cex.lab = cex_val)
text(-260,1, "Disc Model", cex = cex_val, font=2)
for (i in which(bootstrap_disk$number_of_ppp  == 500)){
  plot_losvd(bootstrap_disk$vel_seq[[i]], bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3)

par(fig = c(w[1,1], w[1,2], 0, 0.35), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        #xlab = "velocity, km/s", 
        xaxt="n",yaxt="n", cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
magaxis(side=2, ylab=expression('\u0394'["X - control"]), majorn=3, cex.lab = cex_val)
for (i in which(bootstrap_disk$number_of_ppp  == 500)){
  plot_residual(bootstrap_disk$vel_seq[[i]], 
                bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_disk$vel_spec[[1]],
                control_vseq = bootstrap_disk$vel_seq[[1]])
}
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])

par(fig = c(w[2,1], w[2,2], 0.15, 0.57), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        xaxt="n", yaxt="n",  ylim = c(0,1.1),
        cex.lab = cex_val)
magaxis(side=1, labels = F)
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3)
for (i in which(bootstrap_disk$number_of_ppp  == 200)){
  plot_losvd(bootstrap_disk$vel_seq[[i]], bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3)


par(fig = c(w[2,1], w[2,2], 0, 0.35), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        xaxt="n",
        yaxt="n", cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
for (i in which(bootstrap_disk$number_of_ppp  == 200)){
  plot_residual(bootstrap_disk$vel_seq[[i]], 
                bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_disk$vel_spec[[1]],
                control_vseq = bootstrap_disk$vel_seq[[1]])
}
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])


par(fig = c(w[3,1], w[3,2], 0.15, 0.57), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        xaxt="n", yaxt="n",  ylim = c(0,1.1),
        cex.lab = cex_val)
magaxis(side=1, labels = F)
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3)
for (i in which(bootstrap_disk$number_of_ppp  == 50)){
  plot_losvd(bootstrap_disk$vel_seq[[i]], bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1)
}
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6)
plot_losvd(bootstrap_disk$vel_seq[[1]], bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3)


par(fig = c(w[3,1], w[3,2], 0, 0.35), new=T)
magplot(0, type="n",
        xlim = c(-400,400), 
        yaxt="n", xaxt="n",
        cex.lab = cex_val)
magaxis(side=1, xlab = "velocity, km/s", mtline=1.5, majorn = 4)
for (i in which(bootstrap_disk$number_of_ppp  == 50)){
  plot_residual(bootstrap_disk$vel_seq[[i]], 
                bootstrap_disk$vel_spec[[i]], col = tol_pal[i], lwd = 1,
                control_losvd = bootstrap_disk$vel_spec[[1]],
                control_vseq = bootstrap_disk$vel_seq[[1]])
}
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = "white", lwd = 6,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])
plot_residual(bootstrap_disk$vel_seq[[1]], 
              bootstrap_disk$vel_spec[[1]], col = tol_pal[1], lwd = 3,
              control_losvd = bootstrap_disk$vel_spec[[1]],
              control_vseq = bootstrap_disk$vel_seq[[1]])

dev.off()
