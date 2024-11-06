# find_kinematics functions sorted into separate functions for Project Shapes
# Kate Harborne 09/03/22

# Properties such as the ellipticity, the principle axis, the measurement radius
# and sersic index would be computed from KIDS-like imaging rather than the SAMI
# maps.

get_photometry_properties = function(photo_flux_image, obs_summary,
                                     rband_limit=24.9, plot=F,
                                     mass_frac=1){

  #' Function to return the photometric properties of the observed galaxy
  #' @param photo_flux_image Numeric array containing the flux/mass map of the
  #'                         "photometric" resolution image
  #' @param obs_summary List containing the details of the observation
  #'                    (i.e. object of the observation class)
  #' @param rband_limit Integer describing minimum r-band magnitude in AB mag.
  #' @param plot Boolean to describe whether QC plots should be printed. Default
  #'             is FALSE.
  #' @returns A list of galaxy properties:
  #'
  #'   properties = list("ellipticity" = numeric(),        # fit half-light ellipticity
  #'                     "ellipticity_error" = numeric(),  # 1-sigma error on value
  #'                     "primary_axis" = numeric(),       # angle in degrees measured from vertical
  #'                     "primary_axis_error" = numeric(), # 1-sigma error on value
  #'                     "major_axis_a_kpc" = numeric(),       # major axis of half-light ellipse in kpc
  #'                     "minor_axis_b_kpc" = numeric(),       # minor axis of half-light ellipse in kpc
  #'                     "sersic_index" = numeric(),       # best fit one component sersic fit to SB profile
  #'                     "sersic_index_error" = numeric()) # 1-sigma error on this value

  # List of properties returned by function.
  properties = list("ellipticity" = numeric(),
                    "ellipticity_error" = numeric(),
                    "primary_axis" = numeric(),
                    "primary_axis_error" = numeric(),
                    "major_axis_a_kpc" = numeric(),
                    "minor_axis_b_kpc" = numeric(),
                    "sersic_index" = numeric(),
                    "sersic_index_error" = numeric())

  # Compute which pixels contain enough flux
  rband_limit_mask = array(data = 0, dim = dim(photo_flux_image))
  rband_limit_mask[photo_flux_image >= magAB2CGS(rband_limit)] = 1

  # mask out insufficient pixels
  masked_photo_flux_image = photo_flux_image * rband_limit_mask
  masked_photo_flux_image[is.na(masked_photo_flux_image)] = 0

  set.seed(42)
  noise = rnorm(length(masked_photo_flux_image), mean = 0, sd = min(masked_photo_flux_image[masked_photo_flux_image>0])/1000)
  noisy_photo_image = masked_photo_flux_image + noise
  
  # Have a look at the initial images given
  if (plot){
    plot_flux(photo_flux_image)
    plot_flux(masked_photo_flux_image)
    suppressWarnings(plot_flux(log10(noisy_photo_image)))
  }

  # Run isophote finder for measurement of half light radius and ellipticity
  lev = 20  # initial number of isophotes fitted
  segim = array(data = 1, dim = dim(photo_flux_image)) # initial segmentation map
  ellipse_data = tryCatch(expr = ProFound::profoundGetEllipses(image=noisy_photo_image,
                                                               pixscale = obs_summary$sbin_size,
                                                               segim=segim, levels = lev, plot = plot),
                          error = function(e){NA})

  while (is.na(ellipse_data[1]) & lev > 6){
    lev = lev - 1
    ellipse_data = tryCatch(expr = ProFound::profoundGetEllipses(image=noisy_photo_image,
                                                                 pixscale = obs_summary$sbin_size,
                                                                 segim=segim, levels = lev, plot = plot),
                            error = function(e){NA})
  }

  if (all(is.na(ellipse_data))){
    cat("Cannot fit minimum number of isophotes. \n")
    properties = list("ellipticity" = NA,
                      "ellipticity_error" = NA,
                      "primary_axis" = NA,
                      "primary_axis_error" = NA,
                      "major_axis_a_kpc" = NA,
                      "minor_axis_b_kpc" = NA,
                      "sersic_index" = NA,
                      "sersic_index_error" =  NA)
    return(properties)
  }  else {
    flux_frac = ellipse_data$ellipses$fluxfrac
    hm_id = which(abs(flux_frac-0.5*mass_frac) == min(abs(flux_frac-0.5*mass_frac)))
    um_id = which(abs(flux_frac-0.6*mass_frac) == min(abs(flux_frac-0.6*mass_frac)))
    lm_id = which(abs(flux_frac-0.4*mass_frac) == min(abs(flux_frac-0.4*mass_frac)))
    isops = lm_id:um_id

    properties$ellipticity = mean(1 - ellipse_data$ellipses$axrat[isops], na.rm=T)
    properties$ellipticity_error = sd(1 - ellipse_data$ellipses$axrat[isops], na.rm=T)
    properties$primary_axis = mean(ellipse_data$ellipses$ang[isops], na.rm=T)
    properties$primary_axis_error = sd(ellipse_data$ellipses$ang[isops], na.rm=T)
    properties$major_axis_a_kpc =  mean(ellipse_data$ellipses$radhi[isops], na.rm=T)#*obs_summary$sbin_size # size in kpc
    properties$minor_axis_b_kpc =  mean(ellipse_data$ellipses$radlo[isops], na.rm=T)#*obs_summary$sbin_size # size in kpc

    if (plot){
      par(pty="m", fig=c(0,1,0,1), family="serif", font=1)
      plot_flux(masked_photo_flux_image, radii = list("a" = properties$major_axis_a_kpc/obs_summary$sbin_size,
                                                      "b" = properties$minor_axis_b_kpc/obs_summary$sbin_size,
                                                      "ang" = properties$primary_axis+90))

      centre = c(obs_summary$sbin/2, obs_summary$sbin/2) # the image centre
      theta = properties$primary_axis * (pi/180) # the angle from vertical in radians
      theta_error = properties$primary_axis_error * (pi/180) # the angle from vertical in radians
      n = obs_summary$sbin # the length of the PA line indicator

      if (theta %/% 90 == 0){ # then we are in the first quadrant
        segments(centre[1], centre[2], (-n*sin(theta))+centre[1], (n*cos(theta))+centre[2], col ="red", lwd= 2)
        segments(centre[1], centre[2], (-n*sin(theta+theta_error))+centre[1], (n*cos(theta+theta_error))+centre[2], col =scales::alpha(rgb(1,0,0), alpha = 0.5), lwd= 2)
        segments(centre[1], centre[2], (-n*sin(theta-theta_error))+centre[1], (n*cos(theta-theta_error))+centre[2], col =scales::alpha(rgb(1,0,0), alpha = 0.5), lwd= 2)
      } else { # we are in the second quadrant
        segments(centre[1], centre[2], (-n*cos(theta-(pi/2)))+centre[1], (-n*sin(theta-(pi/2)))+centre[2], col ="red", lwd= 2)
        segments(centre[1], centre[2], (-n*cos(theta+theta_error-(pi/2)))+centre[1], (-n*sin(theta+theta_error-(pi/2)))+centre[2], col =scales::alpha(rgb(255,0,0), 0.5), lwd= 2)
        segments(centre[1], centre[2], (-n*cos(theta-theta_error-(pi/2)))+centre[1], (-n*sin(theta-theta_error-(pi/2)))+centre[2], col =scales::alpha(rgb(255,0,0), 0.5), lwd= 2)
      }

      legend("topright", legend = expression("PA with 1 "*paste(sigma)*" errors"), text.col = "red", bty="n")
    }

    # Computing the sersic index using ProFit
    start = 1
    model_fit = tryCatch(expr = suppressWarnings(
      optim(onecomponent, par=c(ellipse_data$ellipses$SB[1],#CGS2magAB(sum(masked_photo_flux_image)),
                                properties$major_axis_a_kpc, 4),
            rad=ellipse_data$ellipses$radhi[start:lev],
            SB=ellipse_data$ellipses$SB[start:lev],
            pixscale=obs_summary$spatial_res,
            method='BFGS')$par), #lower=c(-100, 0.5, 1), upper=c(rband_limit, 300, 20)
      error = function(e){NA})

    while (all(is.na(model_fit)) & lev >= 10){
      lev = lev - 1
      model_fit = tryCatch(expr = suppressWarnings(
        optim(onecomponent, par=c(ellipse_data$ellipses$SB[1],
                                  properties$major_axis_a_kpc, 4),
              rad=ellipse_data$ellipses$radhi[start:lev],
              SB=ellipse_data$ellipses$SB[start:lev],
              pixscale=obs_summary$spatial_res,
              method='BFGS')$par), #lower=c(-100, 0.5, 1), upper=c(rband_limit, 300, 20)
        error = function(e){NA})
    }

    residual = 100 * onecomponent_residuals(par=model_fit,
                                            rad = ellipse_data$ellipses$radhi,
                                            SB = ellipse_data$ellipses$SB,
                                            pixscale = obs_summary$spatial_res)/ellipse_data$ellipses$SB

    if (plot){
      rlocs=seq(0,500,by=0.1)

      par(pty="m", fig=c(0,1,0,0.5), family="serif", font=1)
      magplot(ellipse_data$ellipses$radhi, residual,
              xlab=' radius / pixels', ylab="residual, %", grid=TRUE,
              cex.axis = 1, cex.lab = 1, family="serif", font=1,
              #ylim=c(-9,9),
              col="black", pch=16)
      polygon(c(-10,max(ellipse_data$ellipses$radhi)+50,max(ellipse_data$ellipses$radhi)+50,-10,-10),
              c(sd(residual), sd(residual), -sd(residual), -sd(residual), sd(residual)), col=adjustcolor("grey", alpha.f = 0.5), border = NA)
      abline(h=0, lwd=2)
      legend("topright", inset=c(0.02, 0.03), legend=bquote(sigma==.(sd(residual)) ~ "%"), text.col = "black", bty="n")

      par(fig=c(0,1,0.28,1), new=TRUE, family="serif", font=1)
      magplot(ellipse_data$ellipses$radhi, ellipse_data$ellipses$SB, grid=TRUE,
              ylim=c(max(ellipse_data$ellipses$SB, na.rm = TRUE)+3,min(ellipse_data$ellipses$SB, na.rm = TRUE)-3),
              type='p', xaxt='n', ylab=expression('mag / arcsec'^2), lwd = 2, family="serif", font=1)
      galaxy=ProFit::profitRadialSersic(rlocs, mag=model_fit[1], re=model_fit[2], nser=model_fit[3])
      lines(rlocs, ProFound::profoundFlux2SB(galaxy, pixscale=obs_summary$spatial_res), col='magenta', lwd=2)
      legend("topright", inset=c(0.02, 0.1), legend=parse(text=sprintf('n == %f', model_fit[3])), text.col = "black", bty="n")

    }

    properties$sersic_index = model_fit[3]
    properties$sersic_index_error = (sd(residual)/100)*model_fit[3]
    return(properties)
  }

}

# Supporting functions for get_photometric_properties():

onecomponent=function(par=c(20, 3, 3), rad, SB, pixscale=1){
  bulge=ProFit::profitRadialSersic(rad, mag=par[1], re=par[2], nser=par[3])
  total=ProFound::profoundFlux2SB(bulge, pixscale=pixscale)
  return=sum((total-SB)^2)
} # fitting a 1D function to the surface brightness profile

onecomponent_residuals = function(par, rad, SB, pixscale=1){
  bulge=ProFit::profitRadialSersic(rad, mag=par[1], re=par[2], nser=par[3])
  total=ProFound::profoundFlux2SB(bulge, pixscale=pixscale)
  return(total-SB)
}

two_component=function(par=c(20, 3, 3), rad, SB, pixscale=1){
  bulge=ProFit::profitRadialSersic(rad, mag=par[1], re=par[2], nser=par[3])
  total=ProFound::profoundFlux2SB(bulge, pixscale=pixscale)
  return=sum((total-SB)^2)
} # fitting a 1D function to the surface brightness profile

magAB2CGS=function(x){10^(-0.4*(x+48.6))}
CGS2magAB=function(x){-2.5*log10(x)-48.6}
