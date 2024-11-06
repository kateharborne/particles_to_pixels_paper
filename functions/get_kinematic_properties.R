# find_kinematics functions with Voronoi support
# Kate Harborne 11/03/22

# Kinematic properties calculated from the mock observed IFU observations using
# the photometric properties for the measurement radius

get_kinematic_properties = function(ifu_output, photometry_properties,
                                    particle_limit, plot=F, voronoi=F){

  #' Function to return the common kinematic properties of the observed galaxy
  #' @param ifu_output SimSpin output from build_datacube.
  #' @param photometry_properties List output from get_photometry_properties.
  #' @param particle_limit Numeric describing the minimum number of particles
  #'                       per pixel for kinematic analysis
  #' @param plot Boolean to describe whether QC plots should be printed. Default
  #'             is FALSE.
  #' @param voronoi Boolean to describe whether Voronoi binning has been used
  #'                when constructing the mock images. Default is FALSE.

  properties = list("ellipticity" = photometry_properties$ellipticity,
                    "ellipticity_error"= photometry_properties$ellipticity_error,
                    "psi" = numeric(1),
                    "psi_error" = numeric(1),
                    "lambda_re" = numeric(1),
                    "vsigma" = numeric(1),
                    "j" = numeric(1),
                    "sigma_m" = numeric(1),
                    "rotational_velocity" = numeric(1),
                    "sersic_index" = photometry_properties$sersic_index,
                    "sersic_index_error" = photometry_properties$sersic_index_error,
                    "major_axis_a_kpc"= photometry_properties$major_axis_a_kpc,
                    "minor_axis_b_kpc"= photometry_properties$minor_axis_b_kpc,
                    "fill_factor" = numeric(1))

  # Compute which pixels contain enough particles for kinematic analysis
  particle_limit_mask = array(data = 0, dim = dim(ifu_output$raw_images$particle_image))
  particle_limit_mask[ifu_output$raw_images$particle_image >= particle_limit] = 1

  if (max(particle_limit_mask) == 0){
    properties = list("ellipticity" = NA,
                      "ellipticity_error"= NA,
                      "psi" = NA,
                      "psi_error" = NA,
                      "lambda_re" = NA,
                      "vsigma" = NA,
                      "j" = NA,
                      "sigma_m" = NA,
                      "rotational_velocity" = NA,
                      "sersic_index" = NA,
                      "sersic_index_error" = NA,
                      "major_axis_a_kpc" = NA,
                      "minor_axis_b_kpc" = NA,
                      "fill_factor" = NA)

    return(properties)

  } else {
    
    if (!"mass_image" %in% names(ifu_output$observed_images)){
      ifu_output$observed_images$mass_image = ifu_output$observed_images$flux_image
    }

    # mask out insufficient pixels
    masked_mass_image = ifu_output$observed_images$mass_image * particle_limit_mask
    masked_velocity_image = ifu_output$observed_images$velocity_image * particle_limit_mask
    masked_dispersion_image = ifu_output$observed_images$dispersion_image * particle_limit_mask

    if (plot){
      plot_flux(masked_mass_image, fig = c(0,0.5,0.5,1),
                units = ("Mass, Msol"),
                radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                             "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                             "ang" = photometry_properties$primary_axis+90))
      plot_age(ifu_output$raw_images$particle_image, fig = c(0.5,1,0.5,1), new = T, units = "Particle Number",
               radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                            "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                            "ang" = photometry_properties$primary_axis+90))
      plot_velocity(masked_velocity_image, fig = c(0,0.5,0,0.5), new = T,
                    radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                                 "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                                 "ang" = photometry_properties$primary_axis+90))
      plot_dispersion(masked_dispersion_image, fig=c(0.5,1,0,0.5), new = T,
               radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                            "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                            "ang" = photometry_properties$primary_axis+90))
    }

    # First recover the coordinates along x (returning central pixel values)
    # Because all images are centred and square, this axis also describes the
    # y axis coordinates (unnecessary to store at this stage)
    xcoord = (ifu_output$observation$sbin_seq +
                (ifu_output$observation$sbin_size/2))[1:ifu_output$observation$sbin]


    # Determine kinematic PA from images
    pa_kin = .find_kinematic_pa(masked_velocity_image,
                                masked_mass_image,
                                ifu_output$observation,
                                xcoord, plot=F)

    if (plot){
      centre = c(ifu_output$observation$sbin/2, ifu_output$observation$sbin/2) # the image centre
      theta_kin = pa_kin*pi/180; theta_pho = photometry_properties$primary_axis*pi/180
      n = ifu_output$observation$sbin # the length of the PA line indicator

      plot_flux(masked_mass_image, fig = c(0,1,0,1),
                units = ("Mass, Msol"),
                radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                             "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                             "ang" = photometry_properties$primary_axis+90))
      if (theta_kin %/% (pi/2) == 0){
        segments(centre[1], centre[2], (-n*sin(theta_kin))+centre[1], (n*cos(theta_kin))+centre[2], col ="red", lwd= 2)
      } else {
        segments(centre[1], centre[2], (-n*cos(theta_kin-(pi/2)))+centre[1], (-n*sin(theta_kin-(pi/2)))+centre[2], col ="red", lwd= 2)
      }
      if (theta_pho %/% (pi/2) == 0){
        segments(centre[1], centre[2], (-n*sin(theta_pho))+centre[1], (n*cos(theta_pho))+centre[2], col ="green", lwd= 2)
      } else {
        segments(centre[1], centre[2], (-n*cos(theta_pho-(pi/2)))+centre[1], (-n*sin(theta_pho-(pi/2)))+centre[2], col ="green", lwd= 2)
      }
      legend("topright", legend = c(expression("PA"["kin"]), expression("PA"["phot"])), lty=c(1,1),
             lwd=c(3,3), col=c("red", "green"), bty="n", inset=c(0.01, 0.02))

      plot_velocity(masked_velocity_image, fig = c(0,1,0,1))
      if (theta_kin %/% (pi/2) == 0){
        segments(centre[1], centre[2], (-n*sin(theta_kin))+centre[1], (n*cos(theta_kin))+centre[2], col ="red", lwd= 2)
      } else {
        segments(centre[1], centre[2], (-n*cos(theta_kin-(pi/2)))+centre[1], (-n*sin(theta_kin-(pi/2)))+centre[2], col ="red", lwd= 2)
      }
      if (theta_pho %/% (pi/2) == 0){
        segments(centre[1], centre[2], (-n*sin(theta_pho))+centre[1], (n*cos(theta_pho))+centre[2], col ="green", lwd= 2)
      } else {
        segments(centre[1], centre[2], (-n*cos(theta_pho-(pi/2)))+centre[1], (-n*sin(theta_pho-(pi/2)))+centre[2], col ="green", lwd= 2)
      }
      legend("topright", legend = c(expression("PA"["kin"]), expression("PA"["phot"])), lty=c(1,1),
             lwd=c(3,3), col=c("red", "green"), bty="n", inset=c(0.01, 0.02))
    }

    # Computing psi given the photometric major axis
    properties$psi = asin(abs(sin((photometry_properties$primary_axis - pa_kin)*pi/180)))*180/pi
    properties$psi_error = photometry_properties$primary_axis_error

    # Next computing the kinematics within the 1Re ellipse:
    # Finding observable kinematics within 1Re mask

    sbin = ifu_output$observation$sbin

    sbin_grid = list("x" = array(data = rep(xcoord, sbin), dim = c(sbin,sbin)),
                     "y" = array(data = rep(xcoord, each=sbin), dim = c(sbin,sbin)))

    mask = array(data = 0, dim = dim(ifu_output$raw_images$particle_image))

    # Semi-major axis in units of kpc
    a = photometry_properties$major_axis_a_kpc
    # Semi-minor axis in units of kpc
    b = photometry_properties$minor_axis_b_kpc

    for (xi in 1:sbin){
      for (yi in 1:sbin){
        xx = sbin_grid$x[xi,yi]
        yy = sbin_grid$y[xi,yi]
        rr = .ellipse(a = a, b = b, x = xx, y = yy, ang_deg = photometry_properties$primary_axis)
        if (rr <= 1){
          mask[xi,yi] = 1
        }
      }
    }
    
    if (plot){
      mask[mask==0] = NA
      plot_flux(masked_mass_image*mask, fig=c(0,0.33,0,1),
                radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                             "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                             "ang" = photometry_properties$primary_axis+90))
      plot_velocity(masked_velocity_image*mask, fig=c(0.33,0.66,0,1), new=T,
                    radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                                 "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                                 "ang" = photometry_properties$primary_axis+90))
      plot_dispersion(masked_dispersion_image*mask, fig=c(0.66,1,0,1), new=T,
                      radii = list("a" = properties$major_axis_a_kpc/ifu_output$observation$sbin_size,
                                   "b" = properties$minor_axis_b_kpc/ifu_output$observation$sbin_size,
                                   "ang" = photometry_properties$primary_axis+90))
      mask[is.na(mask)] = 0
    }

    properties$fill_factor = sum(mask)/(pi*(a/ifu_output$observation$sbin_size)*(b/ifu_output$observation$sbin_size))
    
    if (voronoi){properties$unique_pixels = sum(!is.na(unique(as.numeric(ifu_output$raw_images$voronoi_bins*mask))))}else{
      properties$unique_pixels = sum(!is.na(as.numeric(mask)))
    }
    
    # Measuring kinematics within the 1Re radius
    kin_props =
      .find_observed_kin(flux_image = masked_mass_image,
                         vel_image = masked_velocity_image,
                         disp_image = masked_dispersion_image,
                         mask = mask,
                         a = photometry_properties$major_axis_a_kpc,
                         b = photometry_properties$minor_axis_b_kpc,
                         ang = photometry_properties$primary_axis,
                         sbin_grid = sbin_grid,
                         plot=plot)

    properties$lambda_re = kin_props$lr
    properties$vsigma = kin_props$vsigma
    properties$j = kin_props$j
    properties$sigma_m = kin_props$sigma_m

    # Finally, fitting for the rotational velocity
    vrot = .get_vrot(vel_image = ifu_output$observed_images$velocity_image,
                     pa_kin = pa_kin,
                     obs_summary = ifu_output$observation,
                     x = xcoord,
                     plot = plot)

    properties$rotational_velocity = vrot

    return(properties)

  }
}

.find_kinematic_pa <- function(velocity_image, flux_image, obs_summary, x, plot=F){
  #' Function to determine the kinematic PA of the input image.
  #' @param velocity_image Numeric array representing the velocity image
  #' @param flux_image Numeric array representing the flux/mass image
  #' @param obs_summary List containing the details of the observation
  #'                            (i.e. object of the observation class)
  #' @param x Numeric array describing the coordinates at the centres of each
  #'                  pixel
  #' @param plot Boolean to describe whether QC plots should be printed. Default
  #'             is FALSE.

  sbin   = obs_summary$sbin
  nsteps = 361 # number of angles through which to rotate the images
  angles = seq(0, 180, length.out = nsteps)
  chi2   = numeric(nsteps) # chi^2 at each angle

  sbin_grid = list("x" = array(data = rep(x, sbin),
                               dim = c(sbin,sbin)),
                   "y" = array(data = rep(x, each=sbin),
                               dim = c(sbin,sbin)))

  # Trimming the segim so that the rotated image has a consistent shape
  r_segim = array(data = NA, dim = c(sbin, sbin))
  if (obs_summary$aperture_shape == "circular"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = .circular_ap((sbin-4))
  }
  if (obs_summary$aperture_shape == "square"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = matrix(data = 1, ncol = (sbin-4), nrow = (sbin-4))
  }
  if (obs_summary$aperture_shape == "hexagonal"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = .hexagonal_ap((sbin-4))
  }

  for (i in 1:nsteps){
    # Begin by rotating the image - transform coordinates, interpolate velocity and flux at new positions
    # need to compute for full sbin_grid
    new_coord = .rotate_points(as.numeric(sbin_grid$x), as.numeric(sbin_grid$y), angles[i])
    new_vel   = .interpolate_image(x, x, velocity_image, new_coord, r_segim, plot)
    new_flux  = .interpolate_image(x, x, flux_image, new_coord, r_segim, plot=F)

    # Next, compute the bi-anti-symmetric velocity map
    basv_map  = .make_biantisymmetric_maps(new_vel$z)
    basv_map_mean = (basv_map[[1]] + basv_map[[2]] + basv_map[[3]] + basv_map[[4]]) / 4

    if (plot){
      plot_velocity(basv_map[[1]], fig = c(0,0.5,0,0.5))
      plot_velocity(basv_map[[2]], fig = c(0.5,1,0,0.5), new=T)
      plot_velocity(basv_map[[3]], fig = c(0,0.5,0.5,1), new=T)
      plot_velocity(basv_map[[4]], fig = c(0.5,1,0.5,1), new=T)
      plot_velocity(basv_map_mean)
    }

    chi2[i] = sum((((basv_map_mean - new_vel$z)*new_flux$z)/new_flux$z)^2, na.rm=T)
  }

  if (plot){
    magicaxis::magplot(angles, chi2, type="l", col="blue", lwd=3)
  }

  pa_kin = angles[which(chi2 == min(chi2))][1]

  return(pa_kin) # adding 90 so that the definition of angles matches profound
}

.circular_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
  y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
  xx = x - xcentre; yy = y - ycentre
  rr = sqrt(xx^2 + yy^2)
  ap_region[rr<= sbin/2] = 1
  return(as.vector(ap_region))
}

.hexagonal_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = x - xcentre
      yy = y - ycentre
      rr = (2 * (sbin / 4) * (sbin * sqrt(3) / 4)) - ((sbin / 4) ) * abs(yy) - ((sbin * sqrt(3) / 4)) * abs(xx)
      if ((rr >= 0) && (abs(xx) < sbin/2) && (abs(yy) < (sbin  * sqrt(3) / 4))){
        ap_region[x,y] = 1
      }
    }
  }
  return(as.vector(ap_region))
}

.rotate_points <- function(x, y, ang){
  #' Function to compute new coordinates in a rotated reference frame.
  #' @param x Numeric array containing coordinates in x-axis
  #' @param y Numeric array containing corresponding coordinates in y-axis
  #' @param ang Numeric describing the angle by which to rotate the image

  theta = (ang) * pi/180
  x_new = x*cos(theta) - y*sin(theta)
  y_new = x*sin(theta) + y*cos(theta)
  return(list("x" = x_new, "y" = y_new))
}

.interpolate_image <- function(x, y, image, new_coord, segim, plot=F){
  #' Function to compute the image values in a modified reference frame
  #' @param x Numeric array describing the original x coordinates of the image
  #' @param y Numeric array describing the original y coordinates of the image
  #' @param image Numeric array describing the original image to be transformed
  #' @param new_coord List containing numeric arrays "x" and "y" to describe the
  #'                  new coordinates at which the image should be evaluated
  #' @param segim Numeric array describing the aperture of the image.
  #' @param plot Boolean to describe whether QC plots should be printed. Default
  #'             is FALSE.

  old_image = list("x" = x, "y" = y, "z" = image)
  new_image = list("x" = new_coord$x, "y" = new_coord$y, "z" = numeric(length(new_coord$x)))
  fields::interp.surface(old_image, cbind(new_coord$x, new_coord$y)) -> new_image$z

  new_image$z = array(new_image$z, dim = c(length(x), length(y))) * segim

  if (plot){
    plot_velocity(new_image$z)
  }

  return(new_image)
}

.make_biantisymmetric_maps<-function(image){
  #' Function to generate a biantisymmetric map from an input image
  #' @param image Numeric array corresponding to an image

  V1 = image
  V2 = t(apply(V1, 1, rev))*-1
  V3 = apply(V1, 2, rev)
  V4 = apply(V2, 2, rev)*-1

  return(list(V1, V2, V3, V4))

}

.ellipse = function(a, b, x, y, ang_deg){
  ang_rad = (ang_deg+90) * (pi/180)
  r = ((x*cos(ang_rad) + y*sin(ang_rad))^2 /a^2) + ((x*sin(ang_rad) - y*cos(ang_rad))^2 / b^2)
  return(r)
}

.find_observed_kin <- function(flux_image, vel_image, disp_image, mask,
                               a, b, ang, sbin_grid, plot=F){

  #' Function to compute common kinematic properties lambda_r, v/sigma,
  #' angular momentum (j) and the mass weighted dispersion.
  #' @param flux_image Numeric array containing the flux/mass map
  #' @param vel_image Numeric array containing the velocity map
  #' @param disp_image Numeric array containing the dispersion map
  #' @param mask Numeric array containing the area of the map to use to compute
  #'             kinematic quantities
  #' @param a Numeric describing the semi-major axis
  #' @param b Numeric describing the semi-minor axis
  #' @param ang Numeric describing the photometric position angle
  #' @param sbin_grid Numeric array describing the coordinates at the centres of each
  #'          pixel
  #' @param plot Boolean to describe whether QC plots should be printed. Default
  #'             is FALSE.

  masked_flux = mask * flux_image
  masked_vel = mask * vel_image
  masked_disp = mask * disp_image

  radius = .ellipse(a, b, sbin_grid$x, sbin_grid$y, ang_deg=ang)
  radius = radius*mask

  if (plot){
    plot_flux(masked_flux, fig = c(0,0.5,0.5,1), units = "Mass, Msol")
    plot_flux(radius, fig = c(0.5,1,0.5,1), new=T, units = "Radius")
    plot_velocity(masked_vel, fig = c(0,0.5,0,0.5), new=T)
    plot_dispersion(masked_disp, fig = c(0.5,1,0,0.5), new=T)
  }

  lr = sum(masked_flux * abs(masked_vel) * radius, na.rm = T) / sum(masked_flux * radius * sqrt(masked_vel^2 + masked_disp^2), na.rm=T)
  vsigma = sqrt(sum((masked_flux*masked_vel)^2, na.rm=T) / sum((masked_flux*masked_disp)^2, na.rm = T))
  j = sum(masked_flux * abs(masked_vel) * radius, na.rm=T) / sum(masked_flux, na.rm=T)
  sigma_m = sum(masked_flux * masked_disp, na.rm=T)/sum(masked_flux, na.rm=T)

  return(list("lr"=lr, "vsigma"=vsigma, "j"=j, "sigma_m"=sigma_m))
}

.get_vrot <- function(vel_image, pa_kin, obs_summary, x, plot){
  #' Function to compute common kinematic properties lambda_r, v/sigma,
  #' angular momentum (j) and the mass weighted dispersion.
  #' @param vel_image Numeric array containing the velocity map
  #' @param pa_kin Numeric describing the kinematic PA
  #' @param obs_summary List containing the details of the observation
  #'                            (i.e. object of the observation class)
  #' @param x Numeric array describing the coordinates at the centres of each
  #'          pixel

  sbin   = obs_summary$sbin
  sbin_grid = list("x" = array(data = rep(x, sbin), dim = c(sbin,sbin)),
                   "y" = array(data = rep(x, each=sbin), dim = c(sbin,sbin)))

  # Trimming the segim so that the rotated image has a consistent shape
  r_segim = array(data = NA, dim = c(sbin, sbin))
  if (obs_summary$aperture_shape == "circular"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = .circular_ap((sbin-4))
  }
  if (obs_summary$aperture_shape == "square"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = matrix(data = 1, ncol = (sbin-4), nrow = (sbin-4))
  }
  if (obs_summary$aperture_shape == "hexagonal"){
    r_segim[3:(sbin-2), 3:(sbin-2)] = .hexagonal_ap((sbin-4))
  }

  r_segim[r_segim==0] = NA
  new_coord = .rotate_points(as.numeric(sbin_grid$x), as.numeric(sbin_grid$y), (pa_kin-90))
  new_vel   = tryCatch(expr = .interpolate_image(x, x, vel_image, new_coord, r_segim, plot=plot),
                       error = function(e){NA})

  if (all(is.na(new_vel))){
    vrot = NA
    return(vrot)
  } else {
    mid_point = (sbin/2)
    end_point = sbin

    radii = sbin_grid$x[,mid_point]
    vel_slit = new_vel$z[,mid_point]
    vel_slit_l = rev(vel_slit[1:mid_point])
    vel_slit_r = vel_slit[(mid_point+1):end_point]

    if (length(vel_slit_l[-which(is.na(vel_slit_l))]) < 4 | length(vel_slit_r[-which(is.na(vel_slit_r))]) < 4 ){
      vrot = NA
      return(vrot)
    } else {

      ind_pos = which(abs(diff(diff(vel_slit_l))) == min(abs(diff(diff(vel_slit_l))), na.rm = T))[1]
      ind_neg = which(abs(diff(diff(vel_slit_r))) == min(abs(diff(diff(vel_slit_r))), na.rm = T))[1]

      if (is.null(ind_pos) | is.null(ind_neg)){
        vrot = NA
        return(vrot)
      } else {
        vrot = median(c(abs(vel_slit_l[ind_neg:mid_point]), abs(vel_slit_r[ind_pos:mid_point])), na.rm = T)
        if (plot){
          magicaxis::magplot(radii, vel_slit, type="p", pch=16, col="black",
                             xlab = "Radius, kpc", ylab = expression("LOS Velocity at PA"["kin"]*", km/s"))

          lines(radii[mid_point:1], vel_slit_l, lwd=1, lty=2, col="blue")
          lines(radii[(mid_point+1):end_point], vel_slit_r, lwd=1, lty=2, col="red")

          points(radii[mid_point:1][ind_neg:mid_point], vel_slit_l[ind_neg:mid_point], pch=16, col="blue")
          points(radii[(mid_point+1):end_point][ind_pos:mid_point], vel_slit_r[ind_pos:mid_point], pch=16, col="red")

          abline(h = vrot, col="green")
          abline(h= -vrot, col="green")
        }

        return(vrot)
      }
    }
  }
}
