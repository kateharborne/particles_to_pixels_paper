# Kate Harborne 01/11/23
# ============================================================================ #
# Taking the batched results from EAGLE and plotting comparison between the
# binned and unbinned observational measurements.
# Plotting figure 9 (and Fig B1) using data created from scripts run on Setonix. 
# ============================================================================ #

library(magicaxis)
library(SimSpin)

source("functions/shadowtext.R")
source("functions/plot_array.R")
source("functions/get_kinematic_properties.R")

# PATHS ========================================================================
kin_out_loc = "data/eagle_kinematics/"
plot_loc = "plots/"

# COMBINING BATCHED DATA RefL0050N0752 =========================================
# COMBINING BATCHED DATA =======================================================
observation_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752/observation_properties_mom2_1.txt"), skip = 5, header=F,
           col.names = c('GalaxyID','inc_deg','twist_deg', 'dist_z', 'seeing_fwhm',
                         'major_axis_a_kpc','minor_axis_b_kpc','angle',
                         'ellipticity','sersic_index'))

unbinned_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752/kinematic_properties_mom2_unbinned_1.txt"), skip=5, header=F,
           col.names = c('GalaxyID','lambda_re','vsigma','j','sigma_m','rotational_velocity', "unique_pixels"))

voronoi_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752/kinematic_properties_mom2_voronoi_binned_1.txt"), skip=4, header=T,
           col.names = c('GalaxyID','vorbin_no','lambda_re','vsigma','j','sigma_m','rotational_velocity', "unique_pixels"))

for (n in 2:20){
  observation_properties = 
    rbind(observation_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752/observation_properties_mom2_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','inc_deg','twist_deg',
                                 'dist_z','seeing_fwhm',
                                 'major_axis_a_kpc','minor_axis_b_kpc','angle',
                                 'ellipticity','sersic_index')))
  
  unbinned_properties =
    rbind(unbinned_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752/kinematic_properties_mom2_unbinned_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','lambda_re','vsigma','j',
                                 'sigma_m','rotational_velocity', "unique_pixels")))
  
  
  voronoi_properties =
    rbind(voronoi_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752/kinematic_properties_mom2_voronoi_binned_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','vorbin_no','lambda_re','vsigma','j',
                                 'sigma_m','rotational_velocity', "unique_pixels")))
  
}

observation_properties = na.omit(observation_properties)
unbinned_properties = na.omit(unbinned_properties)
voronoi_properties = na.omit(voronoi_properties)

# Adding stellar mass to observation properties --------------------------------
eagle = read.csv("data/RefL0050N0752_subhalos_1e10.txt",
                 comment.char = "#")
observation_properties$StellarMass = 0

for (each in 1:length(observation_properties$GalaxyID)){
  observation_properties$StellarMass[each] = eagle$StellarMass_Msol[eagle$GalaxyID == observation_properties$GalaxyID[each]]
}

remove(eagle)

voronoi_properties$ellipticity = 0
voronoi_properties$StellarMass = 0

for (ob in 1:length(voronoi_properties$GalaxyID)){
  id = which(observation_properties$GalaxyID %in% voronoi_properties$GalaxyID[ob])   
  voronoi_properties$ellipticity[ob] = observation_properties$ellipticity[id]
  voronoi_properties$StellarMass[ob] = observation_properties$StellarMass[id]
}

min_num_part = 50

sr_frac = (length(which(observation_properties$ellipticity < 0.4 & 
                          (0.08+(observation_properties$ellipticity/4)) > unbinned_properties$lambda_re)) / 2964)*100

no_systems = length(unbinned_properties$unique_pixels > min_num_part)
no_SRs = 0

for (ppp in c(200,300,400,500)){
  
  sr_frac = 
    rbind(sr_frac,
          (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                          (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                          voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part])) / 
             length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]))*100
    )
  
  no_systems =
    rbind(no_systems,
          length(voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]))
  
  no_SRs = 
    rbind(no_SRs,
          length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                         (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                         voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part])))
}


# Making the final plot --------------------------------------------------------

vpal = c(rgb(126, 41, 84, maxColorValue = 256),
         rgb(159, 74,150, maxColorValue = 256),
         rgb(194,106,119, maxColorValue = 256),
         rgb(220,205,125, maxColorValue = 256),
         rgb(221,221,221, maxColorValue = 256))

vpal_alpha = c(rgb(126/256, 41/256, 84/256, alpha = 0.5),
               rgb(159/256, 74/256,150/256, alpha = 0.5),
               rgb(194/256,106/256,119/256, alpha = 0.5),
               rgb(220/256,205/256,125/256, alpha = 0.5),
               rgb(221/256,221/256,221/256, alpha = 0.5))

el_in = seq(0.001, 0.95, length.out = 100)
beta = 0.7 * el_in
el_90 = sqrt(1 - ((1 - el_in)^2))
oe_90 = 0.5  * (((asin(el_90) / sqrt(1 - (el_90^2))) - el_90) / (el_90 - (asin(el_90) * sqrt(1 - (el_90^2)))))
vsig_90 = sqrt((((1 - beta) * oe_90) - 1) / ((0.15 * (1 - beta) * oe_90) + 1))
lr_90 = (1.1 * vsig_90) / sqrt(1 + ((1.1^2) * (vsig_90^2)))

el_sr = seq(0.001, 0.4, length.out = 50)
lr_sr = 0.08+(el_sr/4)

w = plot_array(0.22, 5)

pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "RefL0050N0752_binning_kinematics_new.pdf"), 
          width = 14, height = 10)

par(fig = c(w[1,1],w[1,2],0.66,0.99),
    mar = c(3.5,3.5,0.5,0.5))
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(observation_properties$ellipticity,
       unbinned_properties$lambda_re,
       pch=21, bg = vpal_alpha[5], col="black")
magcon(observation_properties$ellipticity,
       unbinned_properties$lambda_re,
       col = vpal[5],
       doim = F, lwd=2,
       add=T)
ub_SR_frac = (length(which(observation_properties$ellipticity < 0.4 & 
                             (0.08+(observation_properties$ellipticity/4)) > unbinned_properties$lambda_re)) / 2964)*100
text(0.17, 0.99, paste0(round(ub_SR_frac, 2), "% SR"), font = 2, cex=1.1)
text(0.80, 0.07,  "n = 0/474", font = 1, cex=1)


par(fig = c(w[2,1],w[2,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200],
       pch=16, col = vpal_alpha[4])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[4], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200],
       col =  vpal[4],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.87, 0.07,  "9/140", font = 1, cex=1)

par(fig = c(w[3,1],w[3,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "", #expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300],
       pch=16, col = vpal_alpha[3])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[3], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "4/76", font = 1, cex=1)

par(fig = c(w[4,1],w[4,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400],
       pch=16, col = vpal_alpha[2])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[2], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "3/48", font = 1, cex=1)

par(fig = c(w[5,1],w[5,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "", #expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500],
       pch=16, col = vpal_alpha[1])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[1], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "1/23", font = 1, cex=1)

# Vsigma

par(fig = c(w[1,1],w[1,2],0.33,0.66), new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(observation_properties$ellipticity,
       unbinned_properties$vsigma,
       pch=16, col = vpal_alpha[5])
points(observation_properties$ellipticity,
       unbinned_properties$vsigma,
       pch=21, bg = vpal_alpha[5], col="black")
magcon(observation_properties$ellipticity[unbinned_properties$vsigma <= 2],
       unbinned_properties$vsigma[unbinned_properties$vsigma <= 2],
       col = vpal[5],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[2,1],w[2,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200],
       pch=16, col = vpal_alpha[4])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[4])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200 & voronoi_properties$vsigma <= 2],
       col = vpal[4],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[3,1],w[3,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300],
       pch=16, col = vpal_alpha[3])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[3])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300 & voronoi_properties$vsigma <= 2],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[4,1],w[4,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400],
       pch=16, col = vpal_alpha[2])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[2])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400 & voronoi_properties$vsigma <= 2],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[5,1],w[5,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500],
       pch=16, col = vpal_alpha[1])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[1])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500 & voronoi_properties$vsigma <= 2],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)

# j m_star

par(fig = c(w[1,1],w[1,2],0,0.33), new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(observation_properties$StellarMass),
       log10(unbinned_properties$j),
       pch=16, col = vpal_alpha[5])
points(log10(observation_properties$StellarMass),
       log10(unbinned_properties$j),
       pch=21, bg = vpal_alpha[5])
magcon(log10(observation_properties$StellarMass)[log10(observation_properties$StellarMass) > 10],
       log10(unbinned_properties$j)[log10(observation_properties$StellarMass) > 10],
       col = vpal[5],
       doim = F, lwd=2,
       add=T)
shadowtext(11.1, 2.25, "No binning", font = 2, cex=1.1)

par(fig = c(w[2,1],w[2,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200]),
       pch=16, col = vpal_alpha[4])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[4])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[4],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "200 particles", font = 2, cex=1.1)

par(fig = c(w[3,1],w[3,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300]),
       pch=16, col = vpal_alpha[3])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[3])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "300 particles", font = 2, cex=1.1)

par(fig = c(w[4,1],w[4,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400]),
       pch=16, col = vpal_alpha[2])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[2])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "400 particles", font = 2, cex=1.1)

par(fig = c(w[5,1],w[5,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500]),
       pch=16, col = vpal_alpha[1])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[1])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "500 particles", font = 2, cex=1.1)
dev.off()

# The same, but with higher DM sampling ========================================
# COMBINING BATCHED DATA =======================================================
observation_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/observation_properties_mom2_1.txt"), skip = 5, header=F,
           col.names = c('GalaxyID','inc_deg','twist_deg', 'dist_z', 'seeing_fwhm',
                         'major_axis_a_kpc','minor_axis_b_kpc','angle',
                         'ellipticity','sersic_index'))

unbinned_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/kinematic_properties_mom2_unbinned_1.txt"), skip=5, header=F,
           col.names = c('GalaxyID','lambda_re','vsigma','j','sigma_m','rotational_velocity', "unique_pixels"))

voronoi_properties = 
  read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/kinematic_properties_mom2_voronoi_binned_1.txt"), skip=4, header=T,
           col.names = c('GalaxyID','vorbin_no','lambda_re','vsigma','j','sigma_m','rotational_velocity', "unique_pixels"))

for (n in 2:20){
  observation_properties = 
    rbind(observation_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/observation_properties_mom2_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','inc_deg','twist_deg',
                                 'dist_z','seeing_fwhm',
                                 'major_axis_a_kpc','minor_axis_b_kpc','angle',
                                 'ellipticity','sersic_index')))
  
  unbinned_properties =
    rbind(unbinned_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/kinematic_properties_mom2_unbinned_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','lambda_re','vsigma','j',
                                 'sigma_m','rotational_velocity', "unique_pixels")))
  
  
  voronoi_properties =
    rbind(voronoi_properties, 
          read.csv(paste0(kin_out_loc, "RefL0050N0752x7_752dm/kinematic_properties_mom2_voronoi_binned_", n,".txt"), 
                   skip = 5, header=F,
                   col.names = c('GalaxyID','vorbin_no','lambda_re','vsigma','j',
                                 'sigma_m','rotational_velocity', "unique_pixels")))
  
}

observation_properties = na.omit(observation_properties)
unbinned_properties = na.omit(unbinned_properties)
voronoi_properties = na.omit(voronoi_properties)

# Adding stellar mass to observation properties --------------------------------
eagle = read.csv("data/RefL0050N0752_7x752dm_subhalos_1e10.txt",
                 comment.char = "#")
observation_properties$StellarMass = 0

for (each in 1:length(observation_properties$GalaxyID)){
  observation_properties$StellarMass[each] = eagle$StellarMass_Msol[eagle$GalaxyID == observation_properties$GalaxyID[each]]
}

remove(eagle)

voronoi_properties$ellipticity = 0
voronoi_properties$StellarMass = 0

for (ob in 1:length(voronoi_properties$GalaxyID)){
  id = which(observation_properties$GalaxyID %in% voronoi_properties$GalaxyID[ob])   
  voronoi_properties$ellipticity[ob] = observation_properties$ellipticity[id]
  voronoi_properties$StellarMass[ob] = observation_properties$StellarMass[id]
}

min_num_part = 50

sr_frac = (length(which(observation_properties$ellipticity < 0.4 & 
                          (0.08+(observation_properties$ellipticity/4)) > unbinned_properties$lambda_re)) / 2964)*100

no_systems = length(unbinned_properties$unique_pixels > min_num_part)
no_SRs = 0

for (ppp in c(200,300,400,500)){
  
  sr_frac = 
    rbind(sr_frac,
          (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                          (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                          voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part])) / 
             length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]))*100
    )
  
  no_systems =
    rbind(no_systems,
          length(voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]))
  
  no_SRs = 
    rbind(no_SRs,
          length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                         (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                         voronoi_properties$lambda_re[voronoi_properties$vorbin_no == ppp & voronoi_properties$unique_pixels > min_num_part])))
}

# Making the final plot --------------------------------------------------------

vpal = c(rgb(126, 41, 84, maxColorValue = 256),
         rgb(159, 74,150, maxColorValue = 256),
         rgb(194,106,119, maxColorValue = 256),
         rgb(220,205,125, maxColorValue = 256),
         rgb(221,221,221, maxColorValue = 256))

vpal_alpha = c(rgb(126/256, 41/256, 84/256, alpha = 0.5),
               rgb(159/256, 74/256,150/256, alpha = 0.5),
               rgb(194/256,106/256,119/256, alpha = 0.5),
               rgb(220/256,205/256,125/256, alpha = 0.5),
               rgb(221/256,221/256,221/256, alpha = 0.5))

el_in = seq(0.001, 0.95, length.out = 100)
beta = 0.7 * el_in
el_90 = sqrt(1 - ((1 - el_in)^2))
oe_90 = 0.5  * (((asin(el_90) / sqrt(1 - (el_90^2))) - el_90) / (el_90 - (asin(el_90) * sqrt(1 - (el_90^2)))))
vsig_90 = sqrt((((1 - beta) * oe_90) - 1) / ((0.15 * (1 - beta) * oe_90) + 1))
lr_90 = (1.1 * vsig_90) / sqrt(1 + ((1.1^2) * (vsig_90^2)))

el_sr = seq(0.001, 0.4, length.out = 50)
lr_sr = 0.08+(el_sr/4)

w = plot_array(0.22, 5)

pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "RefL0050N0752x7_binning_kinematics.pdf"), 
          width = 14, height = 10)

par(fig = c(w[1,1],w[1,2],0.66,0.99),
    mar = c(3.5,3.5,0.5,0.5))
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(observation_properties$ellipticity,
       unbinned_properties$lambda_re,
       pch=21, bg = vpal_alpha[5], col="black")
magcon(observation_properties$ellipticity,
       unbinned_properties$lambda_re,
       col = vpal[5],
       doim = F, lwd=2,
       add=T)
ub_SR_frac = (length(which(observation_properties$ellipticity < 0.4 & 
                             (0.08+(observation_properties$ellipticity/4)) > unbinned_properties$lambda_re)) / 2964)*100
text(0.17, 0.99, paste0(round(ub_SR_frac, 2), "% SR"), font = 2, cex=1.1)
text(0.80, 0.07,  "n = 0/469", font = 1, cex=1)


par(fig = c(w[2,1],w[2,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200],
       pch=16, col = vpal_alpha[4])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[4], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200],
       col =  vpal[4],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.87, 0.07,  "7/149", font = 1, cex=1)

par(fig = c(w[3,1],w[3,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "", #expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300],
       pch=16, col = vpal_alpha[3])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[3], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "7/79", font = 1, cex=1)

par(fig = c(w[4,1],w[4,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400],
       pch=16, col = vpal_alpha[2])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[2], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 400  & voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "3/48", font = 1, cex=1)

par(fig = c(w[5,1],w[5,2],0.66,0.99),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "", #expression('\u03BB'['R']),
        xlim = c(0,1),
        ylim = c(0,1.05))
lines(el_in, lr_90, col=rgb(93/256,168/256,153/256), lwd=4)
lines(c(el_sr,0.4), c(lr_sr,0), col="black", lwd = 4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500],
       pch=16, col = vpal_alpha[1])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[1], col="black")
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500],
       voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)
vb_SR_frac = (length(which(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part] < 0.4 & 
                             (0.08+(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part]/4)) > 
                             voronoi_properties$lambda_re[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part])) / 
                length(voronoi_properties$ellipticity[voronoi_properties$vorbin_no == 500& voronoi_properties$unique_pixels > min_num_part]))*100
text(0.2, 0.99, paste0(round(vb_SR_frac, 0), "% SR"), font = 2, cex=1.1)
text(0.89, 0.07,  "2/28", font = 1, cex=1)

# Vsigma

par(fig = c(w[1,1],w[1,2],0.33,0.66), new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(observation_properties$ellipticity,
       unbinned_properties$vsigma,
       pch=16, col = vpal_alpha[5])
points(observation_properties$ellipticity,
       unbinned_properties$vsigma,
       pch=21, bg = vpal_alpha[5], col="black")
magcon(observation_properties$ellipticity[unbinned_properties$vsigma <= 2],
       unbinned_properties$vsigma[unbinned_properties$vsigma <= 2],
       col = vpal[5],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[2,1],w[2,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200],
       pch=16, col = vpal_alpha[4])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[4])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==200 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==200 & voronoi_properties$vsigma <= 2],
       col = vpal[4],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[3,1],w[3,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300],
       pch=16, col = vpal_alpha[3])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[3])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==300 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==300 & voronoi_properties$vsigma <= 2],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[4,1],w[4,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400],
       pch=16, col = vpal_alpha[2])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[2])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==400 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==400 & voronoi_properties$vsigma <= 2],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)

par(fig = c(w[5,1],w[5,2],0.33,0.66),
    new=T)
magplot(0, type="n", 
        xlab = expression('\u03B5'),
        ylab = "",#expression("V/"*'\u03C3'),
        xlim = c(0,1),
        ylim = c(0,2.05))
lines(el_in, vsig_90, col=rgb(93/256,168/256,153/256), lwd=4)
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500],
       pch=16, col = vpal_alpha[1])
points(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500 & voronoi_properties$unique_pixels > min_num_part],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500 & voronoi_properties$unique_pixels > min_num_part],
       pch=21, bg = vpal_alpha[1])
magcon(voronoi_properties$ellipticity[voronoi_properties$vorbin_no==500 & voronoi_properties$vsigma <= 2],
       voronoi_properties$vsigma[voronoi_properties$vorbin_no==500 & voronoi_properties$vsigma <= 2],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)

# j m_star

par(fig = c(w[1,1],w[1,2],0,0.33), new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(observation_properties$StellarMass),
       log10(unbinned_properties$j),
       pch=16, col = vpal_alpha[5])
points(log10(observation_properties$StellarMass),
       log10(unbinned_properties$j),
       pch=21, bg = vpal_alpha[5])
magcon(log10(observation_properties$StellarMass)[log10(observation_properties$StellarMass) > 10],
       log10(unbinned_properties$j)[log10(observation_properties$StellarMass) > 10],
       col = vpal[5],
       doim = F, lwd=2,
       add=T)
shadowtext(11.1, 2.25, "No binning", font = 2, cex=1.1)

par(fig = c(w[2,1],w[2,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200]),
       pch=16, col = vpal_alpha[4])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[4])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 200])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 200])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[4],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "200 particles", font = 2, cex=1.1)

par(fig = c(w[3,1],w[3,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300]),
       pch=16, col = vpal_alpha[3])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[3])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 300])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 300])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[3],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "300 particles", font = 2, cex=1.1)

par(fig = c(w[4,1],w[4,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400]),
       pch=16, col = vpal_alpha[2])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[2])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 400])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 400])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[2],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "400 particles", font = 2, cex=1.1)

par(fig = c(w[5,1],w[5,2],0,0.33),
    new=T)
magplot(0, type="n", 
        xlab = expression("log"[10]*" M"['\u002A']*" [ M"['\u0298']*" ] "),
        ylab = "",#expression("log"[10]*" j"['\u002A']*" [ kpc km/s ] "),
        xlim = c(9.9,11.5),
        ylim = c(-0.1,2.5))
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500]),
       pch=16, col = vpal_alpha[1])
points(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part]),
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500 & voronoi_properties$unique_pixels > min_num_part]),
       pch=21, bg = vpal_alpha[1])
magcon(log10(voronoi_properties$StellarMass[voronoi_properties$vorbin_no == 500])[log10(voronoi_properties$StellarMass) > 10],
       log10(voronoi_properties$j[voronoi_properties$vorbin_no == 500])[log10(voronoi_properties$StellarMass) > 10],
       col = vpal[1],
       doim = F, lwd=2,
       add=T)
shadowtext(11.05, 2.25, "500 particles", font = 2, cex=1.1)
dev.off()