# Kate Harborne 18/10/23
# ============================================================================ #
# Plotting the results of bootstrapping: 
# Generating the plots that show the bootstrap results in comparison with the
# resolution models. 
# This script produces Figures 3 and 4 in the manuscript.
# ============================================================================ #

library(SimSpin)
library(magicaxis)
source("functions/plot_array.R")

# PATHS ========================================================================
ss_loc = "simspin_files/"
data_loc = "data/"
plot_loc = "plots/"
cube_loc = paste0(data_loc, "/cubes/")

# Functions for plotting -------------------------------------------------------

plot_errors = function(x, mn, p25, p75, col, type="l"){
  if (type == "l"){
    lines(x, mn, lwd=3, col = col, lty=1)
    polygon(c(x, x[seq(length(x), 1)], x[1]), 
            c(p25, p75[seq(length(x), 1)], p25[1]),
            col = adjustcolor(col, alpha.f = 0.2), border = NA)
  }
  if (type == "l2"){
    lines(x, mn, lwd=3, col = col, lty=2)
    lines(x, p25, lwd=1, col = col, lty=2)
    lines(x, p75, lwd=1, col = col, lty=2)
  }
  if (type == "p"){
    for (each in 1:length(mn)){
      lines(c(x[each,2], x[each,2]), c(p25[each], p75[each]), lwd = 1, col = col[1])
      lines(c(x[each,1], x[each,3]), c(mn[each], mn[each]), lwd = 1, col = col[1])
    }
    
    if (length(col) == 2){
      points(x[,2], mn, pch = 21, col = col[1], bg = col[2])
    } else{
      points(x[,2], mn, pch = 16, col = col[1])
    }
    points(x[,2], p25, pch="-", col = col[1])
    points(x[,2], p75, pch="-", col = col[1])
    points(x[,1], mn, pch="|", col = col[1])
    points(x[,3], mn, pch = "|", col = col[1])
  }
}

tol_pal = c("#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255",
            "#332288", "#117733", "#44AA99", "#88CCEE", 
            "#DDCC77", "#CC6677", "#AA4499", "#882255")

w = plot_array(0.23, 6)

# Loading the data =============================================================

disk_res_data = readRDS(paste0(data_loc, "disk_650_res_data_mom4.Rdata"))
bulge_res_data = readRDS(paste0(data_loc, "bulge_650_res_data_mom4.Rdata"))

disk_res_data_mom2 = readRDS(paste0(data_loc, "disk_650_res_data_mom2.Rdata"))
bulge_res_data_mom2 = readRDS(paste0(data_loc, "bulge_650_res_data_mom2.Rdata"))

disk_velocity = readRDS(paste0(data_loc,"bootstrap_avgs_disk_inc90_sm650.Rdata"))
bulge_velocity = readRDS(paste0(data_loc,"bootstrap_avgs_bulge_inc00_sm650.Rdata"))

bootstrap_disk = readRDS(paste0(data_loc, "bootstrap_disk_inc90_sm650.Rdata"))
bootstrap_bulge = readRDS(paste0(data_loc, "bootstrap_bulge_inc00_sm650.Rdata"))

ppp = sort(seq(50, 500, by=10), decreasing = T)
true_disk = bootstrap_disk[!bootstrap_disk$number_of_ppp %in% ppp,]
true_bulge = bootstrap_bulge[!bootstrap_bulge$number_of_ppp %in% ppp,]

# Compute "average" residual across all pixels for mom=2 mom=4 bulge comparison ------

mean_bulge_v2 = (bulge_velocity$mean_v2[bulge_velocity$pixel_id == 1] - 
                   true_bulge$velocity_2[true_bulge$pixel_id == 1]) / true_bulge$velocity_2[true_bulge$pixel_id == 1]

p25_bulge_v2 =  (bulge_velocity$p25_v2[bulge_velocity$pixel_id == 1]- 
                   true_bulge$velocity_2[true_bulge$pixel_id == 1]) / true_bulge$velocity_2[true_bulge$pixel_id == 1]

p75_bulge_v2 = (bulge_velocity$p75_v2[bulge_velocity$pixel_id == 1]- 
                  true_bulge$velocity_2[true_bulge$pixel_id == 1]) / true_bulge$velocity_2[true_bulge$pixel_id == 1]

mean_bulge_v4 = (bulge_velocity$mean_v4[bulge_velocity$pixel_id == 1] - 
                   true_bulge$velocity_4[true_bulge$pixel_id == 1]) / true_bulge$velocity_4[true_bulge$pixel_id == 1]

p25_bulge_v4 =  (bulge_velocity$p25_v4[bulge_velocity$pixel_id == 1]- 
                   true_bulge$velocity_4[true_bulge$pixel_id == 1]) / true_bulge$velocity_4[true_bulge$pixel_id == 1]

p75_bulge_v4 = (bulge_velocity$p75_v4[bulge_velocity$pixel_id == 1]- 
                  true_bulge$velocity_4[true_bulge$pixel_id == 1]) / true_bulge$velocity_4[true_bulge$pixel_id == 1]

mean_bulge_d2 = (bulge_velocity$mean_d2[bulge_velocity$pixel_id == 1] - 
                   true_bulge$dispersion_2[true_bulge$pixel_id == 1]) / true_bulge$dispersion_2[true_bulge$pixel_id == 1]

p25_bulge_d2 = (bulge_velocity$p25_d2[bulge_velocity$pixel_id == 1]- 
                  true_bulge$dispersion_2[true_bulge$pixel_id == 1]) / true_bulge$dispersion_2[true_bulge$pixel_id == 1]

p75_bulge_d2 = (bulge_velocity$p75_d2[bulge_velocity$pixel_id == 1]- 
                  true_bulge$dispersion_2[true_bulge$pixel_id == 1]) / true_bulge$dispersion_2[true_bulge$pixel_id == 1]

mean_bulge_d4 = (bulge_velocity$mean_d4[bulge_velocity$pixel_id == 1] - 
                   true_bulge$dispersion_4[true_bulge$pixel_id == 1]) / true_bulge$dispersion_4[true_bulge$pixel_id == 1]

p25_bulge_d4 = (bulge_velocity$p25_d4[bulge_velocity$pixel_id == 1]- 
                  true_bulge$dispersion_4[true_bulge$pixel_id == 1]) / true_bulge$dispersion_4[true_bulge$pixel_id == 1]

p75_bulge_d4 = (bulge_velocity$p75_d4[bulge_velocity$pixel_id == 1]- 
                  true_bulge$dispersion_4[true_bulge$pixel_id == 1]) / true_bulge$dispersion_4[true_bulge$pixel_id == 1]

for (n in 2:16){
  
  mean_bulge_v2 = mean_bulge_v2 +
    (bulge_velocity$mean_v2[bulge_velocity$pixel_id == n] - 
       true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) 
       )
  
  p25_bulge_v2 = p25_bulge_v2 +
    (bulge_velocity$p25_v2[bulge_velocity$pixel_id == n]- 
       true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) 
       )
  
  p75_bulge_v2 = p75_bulge_v2 +
    (bulge_velocity$p75_v2[bulge_velocity$pixel_id == n]- 
       true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) 
       )
  
  mean_bulge_v4 = mean_bulge_v4 +
    (bulge_velocity$mean_v4[bulge_velocity$pixel_id == n] - 
       true_bulge$velocity_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_4[true_bulge$pixel_id == n]) 
       )
  
  p25_bulge_v4 = p25_bulge_v4 +
    (bulge_velocity$p25_v4[bulge_velocity$pixel_id == n]- 
       true_bulge$velocity_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_4[true_bulge$pixel_id == n]) 
       )
  
  p75_bulge_v4 = p75_bulge_v4 +
    (bulge_velocity$p75_v4[bulge_velocity$pixel_id == n]- 
       true_bulge$velocity_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$velocity_4[true_bulge$pixel_id == n])
       )
  
  mean_bulge_d2 = mean_bulge_d2 +
    (bulge_velocity$mean_d2[bulge_velocity$pixel_id == n] - 
       true_bulge$dispersion_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_2[true_bulge$pixel_id == n]) 
       )
  
  p25_bulge_d2 = p25_bulge_d2 +
    (bulge_velocity$p25_d2[bulge_velocity$pixel_id == n]- 
       true_bulge$dispersion_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_2[true_bulge$pixel_id == n]) 
       )
  
  p75_bulge_d2 = p75_bulge_d2 +
    (bulge_velocity$p75_d2[bulge_velocity$pixel_id == n]- 
       true_bulge$dispersion_2[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_2[true_bulge$pixel_id == n]) 
       )
  
  mean_bulge_d4 = mean_bulge_d4 +
    (bulge_velocity$mean_d4[bulge_velocity$pixel_id == n] - 
       true_bulge$dispersion_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_4[true_bulge$pixel_id == n]) 
       )
  
  p25_bulge_d4 = p25_bulge_d4 +
    (bulge_velocity$p25_d4[bulge_velocity$pixel_id == n]- 
       true_bulge$dispersion_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_4[true_bulge$pixel_id == n]) 
       )
  
  p75_bulge_d4 = p75_bulge_d4 +
    (bulge_velocity$p75_d4[bulge_velocity$pixel_id == n]- 
       true_bulge$dispersion_4[true_bulge$pixel_id == n]) / (
         abs(true_bulge$dispersion_4[true_bulge$pixel_id == n])
       )
  
}

mean_bulge_v2 = mean_bulge_v2/16
p25_bulge_v2 = p25_bulge_v2/16
p75_bulge_v2 = p75_bulge_v2/16
mean_bulge_v4 = mean_bulge_v4/16
p25_bulge_v4 = p25_bulge_v4/16
p75_bulge_v4 = p75_bulge_v4/16

mean_bulge_d2 = mean_bulge_d2/16
p25_bulge_d2 = p25_bulge_d2/16
p75_bulge_d2 = p75_bulge_d2/16
mean_bulge_d4 = mean_bulge_d4/16
p25_bulge_d4 = p25_bulge_d4/16
p75_bulge_d4 = p75_bulge_d4/16

# Plotting the bulge results ===================================================
pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "bulge_bootstrap_model.pdf"), 
          width = 5.5, height = 10)

par(fig = c(0,1,w[6,1],w[6,2]), 
    mar = c(3.5,3.5,0.5,0.5))
magplot(0, type = "n", grid = T,
        ylim = c(-15,15), xlim = c(50,500),
        xaxt = "n", 
        ylab =  expression('\u0394'*'V / V'["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_v2[bulge_velocity$pixel_id == n] - 
                 true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
                   (abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) + abs(bulge_velocity$mean_v2[bulge_velocity$pixel_id == n]))/2 
                 ),
              (bulge_velocity$p25_v2[bulge_velocity$pixel_id == n]- 
                 true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
                   (abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) + abs(bulge_velocity$mean_v2[bulge_velocity$pixel_id == n]))/2 
                 ),
              (bulge_velocity$p75_v2[bulge_velocity$pixel_id == n]- 
                 true_bulge$velocity_2[true_bulge$pixel_id == n]) / (
                   (abs(true_bulge$velocity_2[true_bulge$pixel_id == n]) + abs(bulge_velocity$mean_v2[bulge_velocity$pixel_id == n]))/2 
                 ),
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

text(450, 97, "Bulge Models", font = 2, cex=1.1)
text(470, 72, "Gaussian", font = 1, cex=1.1)

plot_errors(bulge_velocity$ppp[1:46],
            mean_bulge_v4,
            p25_bulge_v4,
            p75_bulge_v4,
            col = "grey", type="l2")

plot_errors(x = rbind(bulge_res_data_mom2$InnerN$p0, 
                      bulge_res_data_mom2$InnerN$p16, 
                      bulge_res_data_mom2$InnerN$p50, 
                      bulge_res_data_mom2$InnerN$p84,
                      bulge_res_data_mom2$InnerN$p100), 
            mn = rbind(bulge_res_data_mom2$InnerVelocity$p0[2], 
                       bulge_res_data_mom2$InnerVelocity$p16[2], 
                       bulge_res_data_mom2$InnerVelocity$p50[2], 
                       bulge_res_data_mom2$InnerVelocity$p84[2], 
                       bulge_res_data_mom2$InnerVelocity$p100[2]),
            p25 = rbind(bulge_res_data_mom2$InnerVelocity$p0[1], 
                        bulge_res_data_mom2$InnerVelocity$p16[1], 
                        bulge_res_data_mom2$InnerVelocity$p50[1], 
                        bulge_res_data_mom2$InnerVelocity$p84[1], 
                        bulge_res_data_mom2$InnerVelocity$p100[1]),
            p75 = rbind(bulge_res_data_mom2$InnerVelocity$p0[3], 
                        bulge_res_data_mom2$InnerVelocity$p16[3], 
                        bulge_res_data_mom2$InnerVelocity$p50[3], 
                        bulge_res_data_mom2$InnerVelocity$p84[3], 
                        bulge_res_data_mom2$InnerVelocity$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(bulge_res_data_mom2$OuterN$p0, 
                      bulge_res_data_mom2$OuterN$p16, 
                      bulge_res_data_mom2$OuterN$p50, 
                      bulge_res_data_mom2$OuterN$p84,
                      bulge_res_data_mom2$OuterN$p100), 
            mn = rbind(bulge_res_data_mom2$OuterVelocity$p0[2], 
                       bulge_res_data_mom2$OuterVelocity$p16[2], 
                       bulge_res_data_mom2$OuterVelocity$p50[2], 
                       bulge_res_data_mom2$OuterVelocity$p84[2], 
                       bulge_res_data_mom2$OuterVelocity$p100[2]),
            p25 = rbind(bulge_res_data_mom2$OuterVelocity$p0[1], 
                        bulge_res_data_mom2$OuterVelocity$p16[1], 
                        bulge_res_data_mom2$OuterVelocity$p50[1], 
                        bulge_res_data_mom2$OuterVelocity$p84[1], 
                        bulge_res_data_mom2$OuterVelocity$p100[1]),
            p75 = rbind(bulge_res_data_mom2$OuterVelocity$p0[3], 
                        bulge_res_data_mom2$OuterVelocity$p16[3], 
                        bulge_res_data_mom2$OuterVelocity$p50[3], 
                        bulge_res_data_mom2$OuterVelocity$p84[3], 
                        bulge_res_data_mom2$OuterVelocity$p100[3]),
            col= c("black","white"), 
            type = "p")

par(fig = c(0,1,w[5,1], w[5,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-0.5,0.15), xlim = c(50,500),
        xaxt = "n",
        ylab = expression('\u0394'*'\u03C3'*" / \u03C3"["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_d2[bulge_velocity$pixel_id == n] - 
                 true_bulge$dispersion_2[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_2[true_bulge$pixel_id == n],
              (bulge_velocity$p25_d2[bulge_velocity$pixel_id == n]- 
                 true_bulge$dispersion_2[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_2[true_bulge$pixel_id == n],
              (bulge_velocity$p75_d2[bulge_velocity$pixel_id == n]- 
                 true_bulge$dispersion_2[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_2[true_bulge$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}

plot_errors(bulge_velocity$ppp[1:46],
            mean_bulge_d4,
            p25_bulge_d4,
            p75_bulge_d4,
            col = "grey", type="l2")

plot_errors(x = rbind(bulge_res_data_mom2$InnerN$p0, 
                      bulge_res_data_mom2$InnerN$p16, 
                      bulge_res_data_mom2$InnerN$p50, 
                      bulge_res_data_mom2$InnerN$p84,
                      bulge_res_data_mom2$InnerN$p100), 
            mn = rbind(bulge_res_data_mom2$InnerDispersion$p0[2], 
                       bulge_res_data_mom2$InnerDispersion$p16[2], 
                       bulge_res_data_mom2$InnerDispersion$p50[2], 
                       bulge_res_data_mom2$InnerDispersion$p84[2], 
                       bulge_res_data_mom2$InnerDispersion$p100[2]),
            p25 = rbind(bulge_res_data_mom2$InnerDispersion$p0[1], 
                        bulge_res_data_mom2$InnerDispersion$p16[1], 
                        bulge_res_data_mom2$InnerDispersion$p50[1], 
                        bulge_res_data_mom2$InnerDispersion$p84[1], 
                        bulge_res_data_mom2$InnerDispersion$p100[1]),
            p75 = rbind(bulge_res_data_mom2$InnerDispersion$p0[3], 
                        bulge_res_data_mom2$InnerDispersion$p16[3], 
                        bulge_res_data_mom2$InnerDispersion$p50[3], 
                        bulge_res_data_mom2$InnerDispersion$p84[3], 
                        bulge_res_data_mom2$InnerDispersion$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(bulge_res_data_mom2$OuterN$p0, 
                      bulge_res_data_mom2$OuterN$p16, 
                      bulge_res_data_mom2$OuterN$p50, 
                      bulge_res_data_mom2$OuterN$p84,
                      bulge_res_data_mom2$OuterN$p100), 
            mn = rbind(bulge_res_data_mom2$OuterDispersion$p0[2], 
                       bulge_res_data_mom2$OuterDispersion$p16[2], 
                       bulge_res_data_mom2$OuterDispersion$p50[2], 
                       bulge_res_data_mom2$OuterDispersion$p84[2], 
                       bulge_res_data_mom2$OuterDispersion$p100[2]),
            p25 = rbind(bulge_res_data_mom2$OuterDispersion$p0[1], 
                        bulge_res_data_mom2$OuterDispersion$p16[1], 
                        bulge_res_data_mom2$OuterDispersion$p50[1], 
                        bulge_res_data_mom2$OuterDispersion$p84[1], 
                        bulge_res_data_mom2$OuterDispersion$p100[1]),
            p75 = rbind(bulge_res_data_mom2$OuterDispersion$p0[3], 
                        bulge_res_data_mom2$OuterDispersion$p16[3], 
                        bulge_res_data_mom2$OuterDispersion$p50[3], 
                        bulge_res_data_mom2$OuterDispersion$p84[3], 
                        bulge_res_data_mom2$OuterDispersion$p100[3]),
            col= c("black","white"), 
            type = "p")

legend("bottomright", legend=c("Average Gauss-Hermite profile fit"), lty=2, lwd=3, col="grey", bty="n")


par(fig = c(0,1,w[4,1],w[4,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-70,120), xlim = c(50,500),
        xaxt = "n", 
        ylab = expression('\u0394'*'V / V'["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_v4[bulge_velocity$pixel_id == n] - 
                 true_bulge$velocity_4[true_bulge$pixel_id == n]) /
                true_bulge$velocity_4[true_bulge$pixel_id == n],
              (bulge_velocity$p25_v4[bulge_velocity$pixel_id == n]- 
                 true_bulge$velocity_4[true_bulge$pixel_id == n]) /
                true_bulge$velocity_4[true_bulge$pixel_id == n],
              (bulge_velocity$p75_v4[bulge_velocity$pixel_id == n]- 
                 true_bulge$velocity_4[true_bulge$pixel_id == n]) /
                true_bulge$velocity_4[true_bulge$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

text(445, 97, "Gauss-Hermite", font = 1, cex=1.1)

plot_errors(bulge_velocity$ppp[1:46],
            mean_bulge_v2,
            p25_bulge_v2,
            p75_bulge_v2,
            col = "grey", type="l2")

plot_errors(x = rbind(bulge_res_data$InnerN$p0, 
                      bulge_res_data$InnerN$p16, 
                      bulge_res_data$InnerN$p50, 
                      bulge_res_data$InnerN$p84,
                      bulge_res_data$InnerN$p100), 
            mn = rbind(bulge_res_data$InnerVelocity$p0[2], 
                       bulge_res_data$InnerVelocity$p16[2], 
                       bulge_res_data$InnerVelocity$p50[2], 
                       bulge_res_data$InnerVelocity$p84[2], 
                       bulge_res_data$InnerVelocity$p100[2]),
            p25 = rbind(bulge_res_data$InnerVelocity$p0[1], 
                        bulge_res_data$InnerVelocity$p16[1], 
                        bulge_res_data$InnerVelocity$p50[1], 
                        bulge_res_data$InnerVelocity$p84[1], 
                        bulge_res_data$InnerVelocity$p100[1]),
            p75 = rbind(bulge_res_data$InnerVelocity$p0[3], 
                        bulge_res_data$InnerVelocity$p16[3], 
                        bulge_res_data$InnerVelocity$p50[3], 
                        bulge_res_data$InnerVelocity$p84[3], 
                        bulge_res_data$InnerVelocity$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(bulge_res_data$OuterN$p0, 
                      bulge_res_data$OuterN$p16, 
                      bulge_res_data$OuterN$p50, 
                      bulge_res_data$OuterN$p84,
                      bulge_res_data$OuterN$p100), 
            mn = rbind(bulge_res_data$OuterVelocity$p0[2], 
                       bulge_res_data$OuterVelocity$p16[2], 
                       bulge_res_data$OuterVelocity$p50[2], 
                       bulge_res_data$OuterVelocity$p84[2], 
                       bulge_res_data$OuterVelocity$p100[2]),
            p25 = rbind(bulge_res_data$OuterVelocity$p0[1], 
                        bulge_res_data$OuterVelocity$p16[1], 
                        bulge_res_data$OuterVelocity$p50[1], 
                        bulge_res_data$OuterVelocity$p84[1], 
                        bulge_res_data$OuterVelocity$p100[1]),
            p75 = rbind(bulge_res_data$OuterVelocity$p0[3], 
                        bulge_res_data$OuterVelocity$p16[3], 
                        bulge_res_data$OuterVelocity$p50[3], 
                        bulge_res_data$OuterVelocity$p84[3], 
                        bulge_res_data$OuterVelocity$p100[3]),
            col= c("black","white"), 
            type = "p")


par(fig = c(0,1,w[3,1], w[3,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-0.5,0.15), xlim = c(50,500),
        xaxt = "n",
        ylab = expression('\u0394'*'\u03C3'*" / \u03C3"["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_d4[bulge_velocity$pixel_id == n] - 
                 true_bulge$dispersion_4[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_4[true_bulge$pixel_id == n],
              (bulge_velocity$p25_d4[bulge_velocity$pixel_id == n]- 
                 true_bulge$dispersion_4[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_4[true_bulge$pixel_id == n],
              (bulge_velocity$p75_d4[bulge_velocity$pixel_id == n]- 
                 true_bulge$dispersion_4[true_bulge$pixel_id == n]) /
                true_bulge$dispersion_4[true_bulge$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}

plot_errors(bulge_velocity$ppp[1:46],
            mean_bulge_d2,
            p25_bulge_d2,
            p75_bulge_d2,
            col = "grey", type="l2")

plot_errors(x = rbind(bulge_res_data$InnerN$p0, 
                      bulge_res_data$InnerN$p16, 
                      bulge_res_data$InnerN$p50, 
                      bulge_res_data$InnerN$p84,
                      bulge_res_data$InnerN$p100), 
            mn = rbind(bulge_res_data$InnerDispersion$p0[2], 
                       bulge_res_data$InnerDispersion$p16[2], 
                       bulge_res_data$InnerDispersion$p50[2], 
                       bulge_res_data$InnerDispersion$p84[2], 
                       bulge_res_data$InnerDispersion$p100[2]),
            p25 = rbind(bulge_res_data$InnerDispersion$p0[1], 
                        bulge_res_data$InnerDispersion$p16[1], 
                        bulge_res_data$InnerDispersion$p50[1], 
                        bulge_res_data$InnerDispersion$p84[1], 
                        bulge_res_data$InnerDispersion$p100[1]),
            p75 = rbind(bulge_res_data$InnerDispersion$p0[3], 
                        bulge_res_data$InnerDispersion$p16[3], 
                        bulge_res_data$InnerDispersion$p50[3], 
                        bulge_res_data$InnerDispersion$p84[3], 
                        bulge_res_data$InnerDispersion$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(bulge_res_data$OuterN$p0, 
                      bulge_res_data$OuterN$p16, 
                      bulge_res_data$OuterN$p50, 
                      bulge_res_data$OuterN$p84,
                      bulge_res_data$OuterN$p100), 
            mn = rbind(bulge_res_data$OuterDispersion$p0[2], 
                       bulge_res_data$OuterDispersion$p16[2], 
                       bulge_res_data$OuterDispersion$p50[2], 
                       bulge_res_data$OuterDispersion$p84[2], 
                       bulge_res_data$OuterDispersion$p100[2]),
            p25 = rbind(bulge_res_data$OuterDispersion$p0[1], 
                        bulge_res_data$OuterDispersion$p16[1], 
                        bulge_res_data$OuterDispersion$p50[1], 
                        bulge_res_data$OuterDispersion$p84[1], 
                        bulge_res_data$OuterDispersion$p100[1]),
            p75 = rbind(bulge_res_data$OuterDispersion$p0[3], 
                        bulge_res_data$OuterDispersion$p16[3], 
                        bulge_res_data$OuterDispersion$p50[3], 
                        bulge_res_data$OuterDispersion$p84[3], 
                        bulge_res_data$OuterDispersion$p100[3]),
            col= c("black","white"), 
            type = "p")

legend("bottomright", legend=c("Average Gaussian profile fit"), lty=2, lwd=3, col="grey", bty="n")

par(fig = c(0,1,w[2,1], w[2,2]), new=T)
magplot(0, type = "n", grid = T,
        xaxt = "n",
        ylim = c(-36,36), xlim = c(50,500),
        ylab = expression('\u0394'*"h"[3]*" / h"["3, control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_3[bulge_velocity$pixel_id == n] - 
                 true_bulge$h3[true_bulge$pixel_id == n]) /
                true_bulge$h3[true_bulge$pixel_id == n],
              (bulge_velocity$p25_3[bulge_velocity$pixel_id == n] - 
                 true_bulge$h3[true_bulge$pixel_id == n]) /
                true_bulge$h3[true_bulge$pixel_id == n],
              (bulge_velocity$p75_3[bulge_velocity$pixel_id == n] -
                 true_bulge$h3[true_bulge$pixel_id == n]) / 
                true_bulge$h3[true_bulge$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}

plot_errors(x = rbind(bulge_res_data$InnerN$p0, 
                      bulge_res_data$InnerN$p16, 
                      bulge_res_data$InnerN$p50, 
                      bulge_res_data$InnerN$p84,
                      bulge_res_data$InnerN$p100), 
            mn = rbind(bulge_res_data$InnerH3$p0[2], 
                       bulge_res_data$InnerH3$p16[2], 
                       bulge_res_data$InnerH3$p50[2], 
                       bulge_res_data$InnerH3$p84[2], 
                       bulge_res_data$InnerH3$p100[2]),
            p25 = rbind(bulge_res_data$InnerH3$p0[1], 
                        bulge_res_data$InnerH3$p16[1], 
                        bulge_res_data$InnerH3$p50[1], 
                        bulge_res_data$InnerH3$p84[1], 
                        bulge_res_data$InnerH3$p100[1]),
            p75 = rbind(bulge_res_data$InnerH3$p0[3], 
                        bulge_res_data$InnerH3$p16[3], 
                        bulge_res_data$InnerH3$p50[3], 
                        bulge_res_data$InnerH3$p84[3], 
                        bulge_res_data$InnerH3$p100[3]),
            col="black", 
            type = "p")

plot_errors(x = rbind(bulge_res_data$OuterN$p0, 
                      bulge_res_data$OuterN$p16, 
                      bulge_res_data$OuterN$p50, 
                      bulge_res_data$OuterN$p84,
                      bulge_res_data$OuterN$p100), 
            mn = rbind(bulge_res_data$OuterH3$p0[2], 
                       bulge_res_data$OuterH3$p16[2], 
                       bulge_res_data$OuterH3$p50[2], 
                       bulge_res_data$OuterH3$p84[2], 
                       bulge_res_data$OuterH3$p100[2]),
            p25 = rbind(bulge_res_data$OuterH3$p0[1], 
                        bulge_res_data$OuterH3$p16[1], 
                        bulge_res_data$OuterH3$p50[1], 
                        bulge_res_data$OuterH3$p84[1], 
                        bulge_res_data$OuterH3$p100[1]),
            p75 = rbind(bulge_res_data$OuterH3$p0[3], 
                        bulge_res_data$OuterH3$p16[3], 
                        bulge_res_data$OuterH3$p50[3], 
                        bulge_res_data$OuterH3$p84[3], 
                        bulge_res_data$OuterH3$p100[3]),
            col= c("black","white"), 
            type = "p")

par(fig = c(0,1,w[1,1], w[1,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-10,40), xlim = c(50,500),
        xlab = "Particle number per pixel",
        ylab = expression('\u0394'*"h"[4]*" / h"["4, control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(bulge_velocity$ppp[bulge_velocity$pixel_id == n],
              (bulge_velocity$mean_4[bulge_velocity$pixel_id == n] -
                 true_bulge$h4[true_bulge$pixel_id == n]) / 
                true_bulge$h4[true_bulge$pixel_id == n],
              (bulge_velocity$p25_4[bulge_velocity$pixel_id == n] -
                 true_bulge$h4[true_bulge$pixel_id == n]) /
                true_bulge$h4[true_bulge$pixel_id == n],
              (bulge_velocity$p75_4[bulge_velocity$pixel_id == n] -
                 true_bulge$h4[true_bulge$pixel_id == n]) /
                true_bulge$h4[true_bulge$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}

plot_errors(x = rbind(bulge_res_data$InnerN$p0, 
                      bulge_res_data$InnerN$p16, 
                      bulge_res_data$InnerN$p50, 
                      bulge_res_data$InnerN$p84,
                      bulge_res_data$InnerN$p100), 
            mn = rbind(bulge_res_data$InnerH4$p0[2], 
                       bulge_res_data$InnerH4$p16[2], 
                       bulge_res_data$InnerH4$p50[2], 
                       bulge_res_data$InnerH4$p84[2], 
                       bulge_res_data$InnerH4$p100[2]),
            p25 = rbind(bulge_res_data$InnerH4$p0[1], 
                        bulge_res_data$InnerH4$p16[1], 
                        bulge_res_data$InnerH4$p50[1], 
                        bulge_res_data$InnerH4$p84[1], 
                        bulge_res_data$InnerH4$p100[1]),
            p75 = rbind(bulge_res_data$InnerH4$p0[3], 
                        bulge_res_data$InnerH4$p16[3], 
                        bulge_res_data$InnerH4$p50[3], 
                        bulge_res_data$InnerH4$p84[3], 
                        bulge_res_data$InnerH4$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(bulge_res_data$OuterN$p0, 
                      bulge_res_data$OuterN$p16, 
                      bulge_res_data$OuterN$p50, 
                      bulge_res_data$OuterN$p84,
                      bulge_res_data$OuterN$p100), 
            mn = rbind(bulge_res_data$OuterH4$p0[2], 
                       bulge_res_data$OuterH4$p16[2], 
                       bulge_res_data$OuterH4$p50[2], 
                       bulge_res_data$OuterH4$p84[2], 
                       bulge_res_data$OuterH4$p100[2]),
            p25 = rbind(bulge_res_data$OuterH4$p0[1], 
                        bulge_res_data$OuterH4$p16[1], 
                        bulge_res_data$OuterH4$p50[1], 
                        bulge_res_data$OuterH4$p84[1], 
                        bulge_res_data$OuterH4$p100[1]),
            p75 = rbind(bulge_res_data$OuterH4$p0[3], 
                        bulge_res_data$OuterH4$p16[3], 
                        bulge_res_data$OuterH4$p50[3], 
                        bulge_res_data$OuterH4$p84[3], 
                        bulge_res_data$OuterH4$p100[3]),
            col= c("black","white"), 
            type = "p")

dev.off()

# Compute "average" residual across all pixels for mom=2 mom=4 disk comparison ------

mean_disk_v2 = (disk_velocity$mean_v2[disk_velocity$pixel_id == 1] - 
                  true_disk$velocity_2[true_disk$pixel_id == 1]) /
  true_disk$velocity_2[true_disk$pixel_id == 1]

p25_disk_v2 =  (disk_velocity$p25_v2[disk_velocity$pixel_id == 1]- 
                  true_disk$velocity_2[true_disk$pixel_id == 1]) /
  true_disk$velocity_2[true_disk$pixel_id == 1]

p75_disk_v2 = (disk_velocity$p75_v2[disk_velocity$pixel_id == 1]- 
                 true_disk$velocity_2[true_disk$pixel_id == 1]) /
  true_disk$velocity_2[true_disk$pixel_id == 1]

mean_disk_v4 = (disk_velocity$mean_v4[disk_velocity$pixel_id == 1] - 
                  true_disk$velocity_4[true_disk$pixel_id == 1]) /
  true_disk$velocity_4[true_disk$pixel_id == 1]

p25_disk_v4 =  (disk_velocity$p25_v4[disk_velocity$pixel_id == 1]- 
                  true_disk$velocity_4[true_disk$pixel_id == 1]) /
  true_disk$velocity_4[true_disk$pixel_id == 1]

p75_disk_v4 = (disk_velocity$p75_v4[disk_velocity$pixel_id == 1]- 
                 true_disk$velocity_4[true_disk$pixel_id == 1]) /
  true_disk$velocity_4[true_disk$pixel_id == 1]


mean_disk_d2 = (disk_velocity$mean_d2[disk_velocity$pixel_id == 1] - 
                  true_disk$dispersion_2[true_disk$pixel_id == 1]) /
  true_disk$dispersion_2[true_disk$pixel_id == 1]

p25_disk_d2 = (disk_velocity$p25_d2[disk_velocity$pixel_id == 1]- 
                 true_disk$dispersion_2[true_disk$pixel_id == 1]) /
  true_disk$dispersion_2[true_disk$pixel_id == 1]

p75_disk_d2 = (disk_velocity$p75_d2[disk_velocity$pixel_id == 1]- 
                 true_disk$dispersion_2[true_disk$pixel_id == 1]) /
  true_disk$dispersion_2[true_disk$pixel_id == 1]

mean_disk_d4 = (disk_velocity$mean_d4[disk_velocity$pixel_id == 1] - 
                  true_disk$dispersion_4[true_disk$pixel_id == 1]) /
  true_disk$dispersion_4[true_disk$pixel_id == 1]

p25_disk_d4 = (disk_velocity$p25_d4[disk_velocity$pixel_id == 1]- 
                 true_disk$dispersion_4[true_disk$pixel_id == 1]) /
  true_disk$dispersion_4[true_disk$pixel_id == 1]

p75_disk_d4 = (disk_velocity$p75_d4[disk_velocity$pixel_id == 1]- 
                 true_disk$dispersion_4[true_disk$pixel_id == 1]) /
  true_disk$dispersion_4[true_disk$pixel_id == 1]

for (n in 2:16){
  
  mean_disk_v2 = mean_disk_v2 +
    (disk_velocity$mean_v2[disk_velocity$pixel_id == n] - 
       true_disk$velocity_2[true_disk$pixel_id == n]) /
    true_disk$velocity_2[true_disk$pixel_id == n]
  
  p25_disk_v2 = p25_disk_v2 +
    (disk_velocity$p25_v2[disk_velocity$pixel_id == n]- 
       true_disk$velocity_2[true_disk$pixel_id == n]) /
    true_disk$velocity_2[true_disk$pixel_id == n]
  
  p75_disk_v2 = p75_disk_v2 +
    (disk_velocity$p75_v2[disk_velocity$pixel_id == n]- 
       true_disk$velocity_2[true_disk$pixel_id == n]) /
    true_disk$velocity_2[true_disk$pixel_id == n]
  
  mean_disk_v4 = mean_disk_v4 +
    (disk_velocity$mean_v4[disk_velocity$pixel_id == n] - 
       true_disk$velocity_4[true_disk$pixel_id == n]) /
    true_disk$velocity_4[true_disk$pixel_id == n]
  
  p25_disk_v4 = p25_disk_v4 +
    (disk_velocity$p25_v4[disk_velocity$pixel_id == n]- 
       true_disk$velocity_4[true_disk$pixel_id == n]) /
    true_disk$velocity_4[true_disk$pixel_id == n]
  
  p75_disk_v4 = p75_disk_v4 +
    (disk_velocity$p75_v4[disk_velocity$pixel_id == n]- 
       true_disk$velocity_4[true_disk$pixel_id == n]) /
    true_disk$velocity_4[true_disk$pixel_id == n]
  
  mean_disk_d2 = mean_disk_d2 +
    (disk_velocity$mean_d2[disk_velocity$pixel_id == n] - 
       true_disk$dispersion_2[true_disk$pixel_id == n]) /
    true_disk$dispersion_2[true_disk$pixel_id == n]
  
  p25_disk_d2 = p25_disk_d2 +
    (disk_velocity$p25_d2[disk_velocity$pixel_id == n]- 
       true_disk$dispersion_2[true_disk$pixel_id == n]) /
    true_disk$dispersion_2[true_disk$pixel_id == n]
  
  p75_disk_d2 = p75_disk_d2 +
    (disk_velocity$p75_d2[disk_velocity$pixel_id == n]- 
       true_disk$dispersion_2[true_disk$pixel_id == n]) /
    true_disk$dispersion_2[true_disk$pixel_id == n]
  
  mean_disk_d4 = mean_disk_d4 +
    (disk_velocity$mean_d4[disk_velocity$pixel_id == n] - 
       true_disk$dispersion_4[true_disk$pixel_id == n]) /
    true_disk$dispersion_4[true_disk$pixel_id == n]
  
  p25_disk_d4 = p25_disk_d4 +
    (disk_velocity$p25_d4[disk_velocity$pixel_id == n]- 
       true_disk$dispersion_4[true_disk$pixel_id == n]) /
    true_disk$dispersion_4[true_disk$pixel_id == n]
  
  p75_disk_d4 = p75_disk_d4 +
    (disk_velocity$p75_d4[disk_velocity$pixel_id == n]- 
       true_disk$dispersion_4[true_disk$pixel_id == n]) /
    true_disk$dispersion_4[true_disk$pixel_id == n]
  
}

mean_disk_v2 = mean_disk_v2/16
p25_disk_v2 = p25_disk_v2/16
p75_disk_v2 = p75_disk_v2/16
mean_disk_v4 = mean_disk_v4/16
p25_disk_v4 = p25_disk_v4/16
p75_disk_v4 = p75_disk_v4/16

mean_disk_d2 = mean_disk_d2/16
p25_disk_d2 = p25_disk_d2/16
p75_disk_d2 = p75_disk_d2/16
mean_disk_d4 = mean_disk_d4/16
p25_disk_d4 = p25_disk_d4/16
p75_disk_d4 = p75_disk_d4/16


# Plotting the disk results ===================================================
pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "disk_bootstrap_model.pdf"), 
          width = 5.5, height = 10)

par(fig = c(0,1,w[6,1],w[6,2]),
    mar = c(3.5,3.5,0.5,0.5))
magplot(0, type = "n", grid = T,
        ylim = c(-0.09,0.09), xlim = c(50,500),
        xaxt = "n", 
        ylab = expression('\u0394'*'V / V'["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_v2[disk_velocity$pixel_id == n] - 
                 true_disk$velocity_2[true_disk$pixel_id == n]) /
                true_disk$velocity_2[true_disk$pixel_id == n],
              (disk_velocity$p25_v2[disk_velocity$pixel_id == n]- 
                 true_disk$velocity_2[true_disk$pixel_id == n]) /
                true_disk$velocity_2[true_disk$pixel_id == n],
              (disk_velocity$p75_v2[disk_velocity$pixel_id == n]- 
                 true_disk$velocity_2[true_disk$pixel_id == n]) /
                true_disk$velocity_2[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

text(455, 0.07, "Disc Models", font = 2, cex=1.1)
text(470, 0.045, "Gaussian", font = 1, cex=1.1)

plot_errors(disk_velocity$ppp[1:46],
            mean_disk_v4,
            p25_disk_v4,
            p75_disk_v4,
            col = "grey", type="l2")

plot_errors(x = rbind(disk_res_data_mom2$InnerN$p0, 
                      disk_res_data_mom2$InnerN$p16, 
                      disk_res_data_mom2$InnerN$p50, 
                      disk_res_data_mom2$InnerN$p84,
                      disk_res_data_mom2$InnerN$p100), 
            mn = rbind(disk_res_data_mom2$InnerVelocity$p0[2], 
                       disk_res_data_mom2$InnerVelocity$p16[2], 
                       disk_res_data_mom2$InnerVelocity$p50[2], 
                       disk_res_data_mom2$InnerVelocity$p84[2], 
                       disk_res_data_mom2$InnerVelocity$p100[2]),
            p25 = rbind(disk_res_data_mom2$InnerVelocity$p0[1], 
                        disk_res_data_mom2$InnerVelocity$p16[1], 
                        disk_res_data_mom2$InnerVelocity$p50[1], 
                        disk_res_data_mom2$InnerVelocity$p84[1], 
                        disk_res_data_mom2$InnerVelocity$p100[1]),
            p75 = rbind(disk_res_data_mom2$InnerVelocity$p0[3], 
                        disk_res_data_mom2$InnerVelocity$p16[3], 
                        disk_res_data_mom2$InnerVelocity$p50[3], 
                        disk_res_data_mom2$InnerVelocity$p84[3], 
                        disk_res_data_mom2$InnerVelocity$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data_mom2$OuterN$p0, 
                      disk_res_data_mom2$OuterN$p16, 
                      disk_res_data_mom2$OuterN$p50, 
                      disk_res_data_mom2$OuterN$p84,
                      disk_res_data_mom2$OuterN$p100), 
            mn = rbind(disk_res_data_mom2$OuterVelocity$p0[2], 
                       disk_res_data_mom2$OuterVelocity$p16[2], 
                       disk_res_data_mom2$OuterVelocity$p50[2], 
                       disk_res_data_mom2$OuterVelocity$p84[2], 
                       disk_res_data_mom2$OuterVelocity$p100[2]),
            p25 = rbind(disk_res_data_mom2$OuterVelocity$p0[1], 
                        disk_res_data_mom2$OuterVelocity$p16[1], 
                        disk_res_data_mom2$OuterVelocity$p50[1], 
                        disk_res_data_mom2$OuterVelocity$p84[1], 
                        disk_res_data_mom2$OuterVelocity$p100[1]),
            p75 = rbind(disk_res_data_mom2$OuterVelocity$p0[3], 
                        disk_res_data_mom2$OuterVelocity$p16[3], 
                        disk_res_data_mom2$OuterVelocity$p50[3], 
                        disk_res_data_mom2$OuterVelocity$p84[3], 
                        disk_res_data_mom2$OuterVelocity$p100[3]),
            col= c("black","white"), 
            type = "p")


par(fig = c(0,1,w[5,1], w[5,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-0.3,0.045), xlim = c(50,500),
        xaxt = "n",
        ylab = expression('\u0394'*'\u03C3'*" / \u03C3"["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_d2[disk_velocity$pixel_id == n] - 
                 true_disk$dispersion_2[true_disk$pixel_id == n]) /
                true_disk$dispersion_2[true_disk$pixel_id == n],
              (disk_velocity$p25_d2[disk_velocity$pixel_id == n]- 
                 true_disk$dispersion_2[true_disk$pixel_id == n]) /
                true_disk$dispersion_2[true_disk$pixel_id == n],
              (disk_velocity$p75_d2[disk_velocity$pixel_id == n]- 
                 true_disk$dispersion_2[true_disk$pixel_id == n]) /
                true_disk$dispersion_2[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}

plot_errors(disk_velocity$ppp[1:46],
            mean_disk_d4,
            p25_disk_d4,
            p75_disk_d4,
            col = "grey", type="l2")

plot_errors(x = rbind(disk_res_data_mom2$InnerN$p0, 
                      disk_res_data_mom2$InnerN$p16, 
                      disk_res_data_mom2$InnerN$p50, 
                      disk_res_data_mom2$InnerN$p84,
                      disk_res_data_mom2$InnerN$p100), 
            mn = rbind(disk_res_data_mom2$InnerDispersion$p0[2], 
                       disk_res_data_mom2$InnerDispersion$p16[2], 
                       disk_res_data_mom2$InnerDispersion$p50[2], 
                       disk_res_data_mom2$InnerDispersion$p84[2], 
                       disk_res_data_mom2$InnerDispersion$p100[2]),
            p25 = rbind(disk_res_data_mom2$InnerDispersion$p0[1], 
                        disk_res_data_mom2$InnerDispersion$p16[1], 
                        disk_res_data_mom2$InnerDispersion$p50[1], 
                        disk_res_data_mom2$InnerDispersion$p84[1], 
                        disk_res_data_mom2$InnerDispersion$p100[1]),
            p75 = rbind(disk_res_data_mom2$InnerDispersion$p0[3], 
                        disk_res_data_mom2$InnerDispersion$p16[3], 
                        disk_res_data_mom2$InnerDispersion$p50[3], 
                        disk_res_data_mom2$InnerDispersion$p84[3], 
                        disk_res_data_mom2$InnerDispersion$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data_mom2$OuterN$p0, 
                      disk_res_data_mom2$OuterN$p16, 
                      disk_res_data_mom2$OuterN$p50, 
                      disk_res_data_mom2$OuterN$p84,
                      disk_res_data_mom2$OuterN$p100), 
            mn = rbind(disk_res_data_mom2$OuterDispersion$p0[2], 
                       disk_res_data_mom2$OuterDispersion$p16[2], 
                       disk_res_data_mom2$OuterDispersion$p50[2], 
                       disk_res_data_mom2$OuterDispersion$p84[2], 
                       disk_res_data_mom2$OuterDispersion$p100[2]),
            p25 = rbind(disk_res_data_mom2$OuterDispersion$p0[1], 
                        disk_res_data_mom2$OuterDispersion$p16[1], 
                        disk_res_data_mom2$OuterDispersion$p50[1], 
                        disk_res_data_mom2$OuterDispersion$p84[1], 
                        disk_res_data_mom2$OuterDispersion$p100[1]),
            p75 = rbind(disk_res_data_mom2$OuterDispersion$p0[3], 
                        disk_res_data_mom2$OuterDispersion$p16[3], 
                        disk_res_data_mom2$OuterDispersion$p50[3], 
                        disk_res_data_mom2$OuterDispersion$p84[3], 
                        disk_res_data_mom2$OuterDispersion$p100[3]),
            col= c("black","white"), 
            type = "p")

legend("bottomright", legend=c("Average Gauss-Hermite profile fit"), lty=2, lwd=3, col="grey", bty="n")

par(fig = c(0,1,w[4,1],w[4,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-0.09,0.09), xlim = c(50,500),
        xaxt = "n", 
        ylab = expression('\u0394'*'V / V'["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_v4[disk_velocity$pixel_id == n] - 
                 true_disk$velocity_4[true_disk$pixel_id == n]) /
                true_disk$velocity_4[true_disk$pixel_id == n],
              (disk_velocity$p25_v4[disk_velocity$pixel_id == n]- 
                 true_disk$velocity_4[true_disk$pixel_id == n]) /
                true_disk$velocity_4[true_disk$pixel_id == n],
              (disk_velocity$p75_v4[disk_velocity$pixel_id == n]- 
                 true_disk$velocity_4[true_disk$pixel_id == n]) /
                true_disk$velocity_4[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

plot_errors(disk_velocity$ppp[1:46],
            mean_disk_v2,
            p25_disk_v2,
            p75_disk_v2,
            col = "grey", type="l2")

text(445, 0.07, "Gauss-Hermite", font = 1, cex=1.1)

plot_errors(x = rbind(disk_res_data$InnerN$p0, 
                      disk_res_data$InnerN$p16, 
                      disk_res_data$InnerN$p50, 
                      disk_res_data$InnerN$p84,
                      disk_res_data$InnerN$p100), 
            mn = rbind(disk_res_data$InnerVelocity$p0[2], 
                       disk_res_data$InnerVelocity$p16[2], 
                       disk_res_data$InnerVelocity$p50[2], 
                       disk_res_data$InnerVelocity$p84[2], 
                       disk_res_data$InnerVelocity$p100[2]),
            p25 = rbind(disk_res_data$InnerVelocity$p0[1], 
                        disk_res_data$InnerVelocity$p16[1], 
                        disk_res_data$InnerVelocity$p50[1], 
                        disk_res_data$InnerVelocity$p84[1], 
                        disk_res_data$InnerVelocity$p100[1]),
            p75 = rbind(disk_res_data$InnerVelocity$p0[3], 
                        disk_res_data$InnerVelocity$p16[3], 
                        disk_res_data$InnerVelocity$p50[3], 
                        disk_res_data$InnerVelocity$p84[3], 
                        disk_res_data$InnerVelocity$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data$OuterN$p0, 
                      disk_res_data$OuterN$p16, 
                      disk_res_data$OuterN$p50, 
                      disk_res_data$OuterN$p84,
                      disk_res_data$OuterN$p100), 
            mn = rbind(disk_res_data$OuterVelocity$p0[2], 
                       disk_res_data$OuterVelocity$p16[2], 
                       disk_res_data$OuterVelocity$p50[2], 
                       disk_res_data$OuterVelocity$p84[2], 
                       disk_res_data$OuterVelocity$p100[2]),
            p25 = rbind(disk_res_data$OuterVelocity$p0[1], 
                        disk_res_data$OuterVelocity$p16[1], 
                        disk_res_data$OuterVelocity$p50[1], 
                        disk_res_data$OuterVelocity$p84[1], 
                        disk_res_data$OuterVelocity$p100[1]),
            p75 = rbind(disk_res_data$OuterVelocity$p0[3], 
                        disk_res_data$OuterVelocity$p16[3], 
                        disk_res_data$OuterVelocity$p50[3], 
                        disk_res_data$OuterVelocity$p84[3], 
                        disk_res_data$OuterVelocity$p100[3]),
            col= c("black","white"), 
            type = "p")


par(fig = c(0,1,w[3,1], w[3,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-0.3,0.045), xlim = c(50,500),
        xaxt = "n",
        ylab = expression('\u0394'*'\u03C3'*" / \u03C3"["control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_d4[disk_velocity$pixel_id == n] - 
                 true_disk$dispersion_4[true_disk$pixel_id == n]) /
                true_disk$dispersion_4[true_disk$pixel_id == n],
              (disk_velocity$p25_d4[disk_velocity$pixel_id == n]- 
                 true_disk$dispersion_4[true_disk$pixel_id == n]) /
                true_disk$dispersion_4[true_disk$pixel_id == n],
              (disk_velocity$p75_d4[disk_velocity$pixel_id == n]- 
                 true_disk$dispersion_4[true_disk$pixel_id == n]) /
                true_disk$dispersion_4[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
  
}


plot_errors(disk_velocity$ppp[1:46],
            mean_disk_d2,
            p25_disk_d2,
            p75_disk_d2,
            col = "grey", type="l2")

plot_errors(x = rbind(disk_res_data$InnerN$p0, 
                      disk_res_data$InnerN$p16, 
                      disk_res_data$InnerN$p50, 
                      disk_res_data$InnerN$p84,
                      disk_res_data$InnerN$p100), 
            mn = rbind(disk_res_data$InnerDispersion$p0[2], 
                       disk_res_data$InnerDispersion$p16[2], 
                       disk_res_data$InnerDispersion$p50[2], 
                       disk_res_data$InnerDispersion$p84[2], 
                       disk_res_data$InnerDispersion$p100[2]),
            p25 = rbind(disk_res_data$InnerDispersion$p0[1], 
                        disk_res_data$InnerDispersion$p16[1], 
                        disk_res_data$InnerDispersion$p50[1], 
                        disk_res_data$InnerDispersion$p84[1], 
                        disk_res_data$InnerDispersion$p100[1]),
            p75 = rbind(disk_res_data$InnerDispersion$p0[3], 
                        disk_res_data$InnerDispersion$p16[3], 
                        disk_res_data$InnerDispersion$p50[3], 
                        disk_res_data$InnerDispersion$p84[3], 
                        disk_res_data$InnerDispersion$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data$OuterN$p0, 
                      disk_res_data$OuterN$p16, 
                      disk_res_data$OuterN$p50, 
                      disk_res_data$OuterN$p84,
                      disk_res_data$OuterN$p100), 
            mn = rbind(disk_res_data$OuterDispersion$p0[2], 
                       disk_res_data$OuterDispersion$p16[2], 
                       disk_res_data$OuterDispersion$p50[2], 
                       disk_res_data$OuterDispersion$p84[2], 
                       disk_res_data$OuterDispersion$p100[2]),
            p25 = rbind(disk_res_data$OuterDispersion$p0[1], 
                        disk_res_data$OuterDispersion$p16[1], 
                        disk_res_data$OuterDispersion$p50[1], 
                        disk_res_data$OuterDispersion$p84[1], 
                        disk_res_data$OuterDispersion$p100[1]),
            p75 = rbind(disk_res_data$OuterDispersion$p0[3], 
                        disk_res_data$OuterDispersion$p16[3], 
                        disk_res_data$OuterDispersion$p50[3], 
                        disk_res_data$OuterDispersion$p84[3], 
                        disk_res_data$OuterDispersion$p100[3]),
            col= c("black","white"), 
            type = "p")


legend("bottomright", legend=c("Average Gaussian profile fit"), lty=2, lwd=3, col="grey", bty="n")


par(fig = c(0,1,w[2,1], w[2,2]), new=T)
magplot(0, type = "n", grid = T,
        xaxt = "n",
        ylim = c(-1.8,1.8), xlim = c(50,500),
        ylab = expression('\u0394'*"h"[3]*" / h"["3, control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_3[disk_velocity$pixel_id == n] - 
                 true_disk$h3[true_disk$pixel_id == n]) /
                true_disk$h3[true_disk$pixel_id == n],
              (disk_velocity$p25_3[disk_velocity$pixel_id == n] - 
                 true_disk$h3[true_disk$pixel_id == n]) /
                true_disk$h3[true_disk$pixel_id == n],
              (disk_velocity$p75_3[disk_velocity$pixel_id == n] -
                 true_disk$h3[true_disk$pixel_id == n]) / 
                true_disk$h3[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

plot_errors(x = rbind(disk_res_data$InnerN$p0, 
                      disk_res_data$InnerN$p16, 
                      disk_res_data$InnerN$p50, 
                      disk_res_data$InnerN$p84,
                      disk_res_data$InnerN$p100), 
            mn = rbind(disk_res_data$InnerH3$p0[2], 
                       disk_res_data$InnerH3$p16[2], 
                       disk_res_data$InnerH3$p50[2], 
                       disk_res_data$InnerH3$p84[2], 
                       disk_res_data$InnerH3$p100[2]),
            p25 = rbind(disk_res_data$InnerH3$p0[1], 
                        disk_res_data$InnerH3$p16[1], 
                        disk_res_data$InnerH3$p50[1], 
                        disk_res_data$InnerH3$p84[1], 
                        disk_res_data$InnerH3$p100[1]),
            p75 = rbind(disk_res_data$InnerH3$p0[3], 
                        disk_res_data$InnerH3$p16[3], 
                        disk_res_data$InnerH3$p50[3], 
                        disk_res_data$InnerH3$p84[3], 
                        disk_res_data$InnerH3$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data$OuterN$p0, 
                      disk_res_data$OuterN$p16, 
                      disk_res_data$OuterN$p50, 
                      disk_res_data$OuterN$p84,
                      disk_res_data$OuterN$p100), 
            mn = rbind(disk_res_data$OuterH3$p0[2], 
                       disk_res_data$OuterH3$p16[2], 
                       disk_res_data$OuterH3$p50[2], 
                       disk_res_data$OuterH3$p84[2], 
                       disk_res_data$OuterH3$p100[2]),
            p25 = rbind(disk_res_data$OuterH3$p0[1], 
                        disk_res_data$OuterH3$p16[1], 
                        disk_res_data$OuterH3$p50[1], 
                        disk_res_data$OuterH3$p84[1], 
                        disk_res_data$OuterH3$p100[1]),
            p75 = rbind(disk_res_data$OuterH3$p0[3], 
                        disk_res_data$OuterH3$p16[3], 
                        disk_res_data$OuterH3$p50[3], 
                        disk_res_data$OuterH3$p84[3], 
                        disk_res_data$OuterH3$p100[3]),
            col= c("black","white"), 
            type = "p")

par(fig = c(0,1,w[1,1], w[1,2]), new=T)
magplot(0, type = "n", grid = T,
        ylim = c(-20,20), xlim = c(50,500),
        xlab = "Particle number per pixel",
        ylab = expression('\u0394'*"h"[4]*" / h"["4, control"]))
abline(h=0, col = "black")
for (n in 1:16){
  plot_errors(disk_velocity$ppp[disk_velocity$pixel_id == n],
              (disk_velocity$mean_4[disk_velocity$pixel_id == n] -
                 true_disk$h4[true_disk$pixel_id == n]) / 
                true_disk$h4[true_disk$pixel_id == n],
              (disk_velocity$p25_4[disk_velocity$pixel_id == n] -
                 true_disk$h4[true_disk$pixel_id == n]) /
                true_disk$h4[true_disk$pixel_id == n],
              (disk_velocity$p75_4[disk_velocity$pixel_id == n] -
                 true_disk$h4[true_disk$pixel_id == n]) /
                true_disk$h4[true_disk$pixel_id == n],
              col = tol_pal[n])
  #  readline(prompt = "Press enter to continue")
}

plot_errors(x = rbind(disk_res_data$InnerN$p0, 
                      disk_res_data$InnerN$p16, 
                      disk_res_data$InnerN$p50, 
                      disk_res_data$InnerN$p84,
                      disk_res_data$InnerN$p100), 
            mn = rbind(disk_res_data$InnerH4$p0[2], 
                       disk_res_data$InnerH4$p16[2], 
                       disk_res_data$InnerH4$p50[2], 
                       disk_res_data$InnerH4$p84[2], 
                       disk_res_data$InnerH4$p100[2]),
            p25 = rbind(disk_res_data$InnerH4$p0[1], 
                        disk_res_data$InnerH4$p16[1], 
                        disk_res_data$InnerH4$p50[1], 
                        disk_res_data$InnerH4$p84[1], 
                        disk_res_data$InnerH4$p100[1]),
            p75 = rbind(disk_res_data$InnerH4$p0[3], 
                        disk_res_data$InnerH4$p16[3], 
                        disk_res_data$InnerH4$p50[3], 
                        disk_res_data$InnerH4$p84[3], 
                        disk_res_data$InnerH4$p100[3]),
            col="black", 
            type = "p")


plot_errors(x = rbind(disk_res_data$OuterN$p0, 
                      disk_res_data$OuterN$p16, 
                      disk_res_data$OuterN$p50, 
                      disk_res_data$OuterN$p84,
                      disk_res_data$OuterN$p100), 
            mn = rbind(disk_res_data$OuterH4$p0[2], 
                       disk_res_data$OuterH4$p16[2], 
                       disk_res_data$OuterH4$p50[2], 
                       disk_res_data$OuterH4$p84[2], 
                       disk_res_data$OuterH4$p100[2]),
            p25 = rbind(disk_res_data$OuterH4$p0[1], 
                        disk_res_data$OuterH4$p16[1], 
                        disk_res_data$OuterH4$p50[1], 
                        disk_res_data$OuterH4$p84[1], 
                        disk_res_data$OuterH4$p100[1]),
            p75 = rbind(disk_res_data$OuterH4$p0[3], 
                        disk_res_data$OuterH4$p16[3], 
                        disk_res_data$OuterH4$p50[3], 
                        disk_res_data$OuterH4$p84[3], 
                        disk_res_data$OuterH4$p100[3]),
            col= c("black","white"), 
            type = "p")

dev.off()

# For the table ================================================================

bootstrap_bulge$res_vel_2 = 0
bootstrap_bulge$res_sig_2 = 0
bootstrap_bulge$res_vel_4 = 0
bootstrap_bulge$res_sig_4 = 0
bootstrap_bulge$res_h3  = 0
bootstrap_bulge$res_h4  = 0
bootstrap_bulge$rel_res_vel_2 = 0
bootstrap_bulge$rel_res_sig_2 = 0
bootstrap_bulge$rel_res_vel_4 = 0
bootstrap_bulge$rel_res_sig_4 = 0
bootstrap_bulge$rel_res_h3  = 0
bootstrap_bulge$rel_res_h4  = 0

bootstrap_disk$res_vel_2  = 0
bootstrap_disk$res_sig_2  = 0
bootstrap_disk$res_vel_4  = 0
bootstrap_disk$res_sig_4  = 0
bootstrap_disk$res_h3   = 0
bootstrap_disk$res_h4   = 0
bootstrap_disk$rel_res_vel_2  = 0
bootstrap_disk$rel_res_sig_2  = 0
bootstrap_disk$rel_res_vel_4  = 0
bootstrap_disk$rel_res_sig_4  = 0
bootstrap_disk$rel_res_h3   = 0
bootstrap_disk$rel_res_h4   = 0
no_obs = length(bootstrap_bulge$velocity_2[bootstrap_bulge$pixel_id == 1])

for (n in 1:16){
  true_pixel_bulge = true_bulge[true_bulge$pixel_id == n,] 
  
  bootstrap_bulge$res_vel_2[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$velocity_2[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$velocity_2, no_obs))
  
  bootstrap_bulge$res_sig_2[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$dispersion_2[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$dispersion_2, no_obs))
  
  bootstrap_bulge$res_vel_4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$velocity_4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$velocity_4, no_obs))
  
  bootstrap_bulge$res_sig_4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$dispersion_4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$dispersion_4, no_obs))
  
  bootstrap_bulge$res_h3[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$h3[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$h3, no_obs))
  
  bootstrap_bulge$res_h4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$h4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$h4, no_obs))
  
  bootstrap_bulge$rel_res_vel_2[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$velocity_2[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$velocity_2, no_obs))/rep(true_pixel_bulge$velocity_2, no_obs)
  
  bootstrap_bulge$rel_res_sig_2[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$dispersion_2[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$dispersion_2, no_obs))/rep(true_pixel_bulge$dispersion_2, no_obs)
  
  bootstrap_bulge$rel_res_vel_4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$velocity_4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$velocity_4, no_obs))/rep(true_pixel_bulge$velocity_4, no_obs)
  
  bootstrap_bulge$rel_res_sig_4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$dispersion_4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$dispersion_4, no_obs))/rep(true_pixel_bulge$dispersion_4, no_obs)
  
  bootstrap_bulge$rel_res_h3[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$h3[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$h3, no_obs))/rep(true_pixel_bulge$h3, no_obs)
  
  bootstrap_bulge$rel_res_h4[bootstrap_bulge$pixel_id == n] =
    (bootstrap_bulge$h4[bootstrap_bulge$pixel_id == n] -
       rep(true_pixel_bulge$h4, no_obs))/rep(true_pixel_bulge$h4, no_obs)
  
  
  true_pixel_disk = true_disk[true_disk$pixel_id == n,] 
  
  bootstrap_disk$res_vel_2[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$velocity_2[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$velocity_2, no_obs))
  
  bootstrap_disk$res_sig_2[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$dispersion_2[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$dispersion_2, no_obs))
  
  bootstrap_disk$res_vel_4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$velocity_4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$velocity_4, no_obs))
  
  bootstrap_disk$res_sig_4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$dispersion_4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$dispersion_4, no_obs))
  
  bootstrap_disk$res_h3[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$h3[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$h3, no_obs))
  
  bootstrap_disk$res_h4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$h4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$h4, no_obs))
  
  bootstrap_disk$rel_res_vel_2[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$velocity_2[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$velocity_2, no_obs))/rep(true_pixel_disk$velocity_2, no_obs)
  
  bootstrap_disk$rel_res_sig_2[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$dispersion_2[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$dispersion_2, no_obs))/rep(true_pixel_disk$dispersion_2, no_obs)
  
  bootstrap_disk$rel_res_vel_4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$velocity_4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$velocity_4, no_obs))/rep(true_pixel_disk$velocity_4, no_obs)
  
  bootstrap_disk$rel_res_sig_4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$dispersion_4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$dispersion_4, no_obs))/rep(true_pixel_disk$dispersion_4, no_obs)
  
  bootstrap_disk$rel_res_h3[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$h3[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$h3, no_obs))/rep(true_pixel_disk$h3, no_obs)
  
  bootstrap_disk$rel_res_h4[bootstrap_disk$pixel_id == n] =
    (bootstrap_disk$h4[bootstrap_disk$pixel_id == n] -
       rep(true_pixel_disk$h4, no_obs))/rep(true_pixel_disk$h4, no_obs)
  
}

# At each number of particles, compute the km/s SD
error_table_bulge = data.frame("ppp" = ppp,
                               "velocity_sd_2" = numeric(length(ppp)),
                               "dispersion_sd_2" = numeric(length(ppp)),
                               "velocity_sd_4" = numeric(length(ppp)),
                               "dispersion_sd_4" = numeric(length(ppp)),
                               "h3_sd" = numeric(length(ppp)),
                               "h4_sd" = numeric(length(ppp)),
                               "rel_velocity_sd_2" = numeric(length(ppp)),
                               "rel_dispersion_sd_2" = numeric(length(ppp)),
                               "rel_velocity_sd_4" = numeric(length(ppp)),
                               "rel_dispersion_sd_4" = numeric(length(ppp)),
                               "rel_h3_sd" = numeric(length(ppp)),
                               "rel_h4_sd" = numeric(length(ppp))) 

error_table_disk  = data.frame("ppp" = ppp,
                               "velocity_sd_2" = numeric(length(ppp)),
                               "dispersion_sd_2" = numeric(length(ppp)),
                               "velocity_sd_4" = numeric(length(ppp)),
                               "dispersion_sd_4" = numeric(length(ppp)),
                               "h3_sd" = numeric(length(ppp)),
                               "h4_sd" = numeric(length(ppp)),
                               "rel_velocity_sd_2" = numeric(length(ppp)),
                               "rel_dispersion_sd_2" = numeric(length(ppp)),
                               "rel_velocity_sd_4" = numeric(length(ppp)),
                               "rel_dispersion_sd_4" = numeric(length(ppp)),
                               "rel_h3_sd" = numeric(length(ppp)),
                               "rel_h4_sd" = numeric(length(ppp))) 

for (p in ppp){
  error_table_bulge$velocity_sd_2[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_vel_2[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$dispersion_sd_2[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_sig_2[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$velocity_sd_4[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_vel_4[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$dispersion_sd_4[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_sig_4[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$h3_sd[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_h3[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$h4_sd[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_h4[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_velocity_sd_2[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$rel_res_vel_2[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_dispersion_sd_2[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$rel_res_sig_2[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_velocity_sd_4[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$rel_res_vel_4[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_dispersion_sd_4[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$rel_res_sig_4[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_h3_sd[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$res_h3[bootstrap_bulge$number_of_ppp == p])
  
  error_table_bulge$rel_h4_sd[error_table_bulge$ppp == p] =
    2*sd(bootstrap_bulge$rel_res_h4[bootstrap_bulge$number_of_ppp == p])
  
  
  error_table_disk$velocity_sd_2[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_vel_2[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$dispersion_sd_2[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_sig_2[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$velocity_sd_4[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_vel_4[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$dispersion_sd_4[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_sig_4[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$h3_sd[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_h3[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$h4_sd[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$res_h4[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_velocity_sd_2[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_vel_2[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_dispersion_sd_2[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_sig_2[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_velocity_sd_4[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_vel_4[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_dispersion_sd_4[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_sig_4[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_h3_sd[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_h3[bootstrap_disk$number_of_ppp == p])
  
  error_table_disk$rel_h4_sd[error_table_disk$ppp == p] =
    2*sd(bootstrap_disk$rel_res_h4[bootstrap_disk$number_of_ppp == p])
}

assymptote_fun = function(par= c(1,1,1), x = ppp, y){
  fit_y = par[1]*log10(x)^(par[2]) + par[3]
  return = sum((fit_y - y)^2)
}

assymptote_plot = function(par, x = ppp){
  y =   par[1]*log10(x)^(par[2]) + par[3]
  lines(x, y, col = tol_pal[5], lwd= 2)
}

assymptote_y = function(par, y){
  x = 10^(((y - par[3])/par[1])^(1/par[2]))
  return(x)
}


# Plots -----------------------------------------------------------------------
w = plot_array(0.23, 6)

assymptote_y(bulge_v_fit, 20)
assymptote_y(disk_v_fit, 20)
assymptote_y(bulge_v_fit, 15)
assymptote_y(disk_v_fit, 15)
assymptote_y(bulge_v_fit, 10)
assymptote_y(disk_v_fit, 10)

assymptote_y(bulge_s_fit, 30)
assymptote_y(disk_s_fit, 30)
assymptote_y(bulge_s_fit, 20)
assymptote_y(disk_s_fit, 20)
assymptote_y(bulge_s_fit, 10)
assymptote_y(disk_s_fit, 10)

assymptote_y(bulge_3_fit, 0.15)
assymptote_y(disk_3_fit, 0.15)
assymptote_y(bulge_3_fit, 0.1)
assymptote_y(disk_3_fit, 0.1)
assymptote_y(bulge_3_fit, 0.05)
assymptote_y(disk_3_fit, 0.05)

assymptote_y(bulge_4_fit, 0.15)
assymptote_y(disk_4_fit, 0.15)
assymptote_y(bulge_4_fit, 0.1)
assymptote_y(disk_4_fit, 0.1)
assymptote_y(bulge_4_fit, 0.05)
assymptote_y(disk_4_fit, 0.05)


pdf.options(fonts="Arial")
cairo_pdf(file=paste0(plot_loc, "p2p_raw_estimation.pdf"), 
          width = 7, height = 10)

par(fig = c(0,0.53,w[6,1],w[6,2]),
    mar = c(3.5,3.5,0.5,0.5))
magplot(error_table_bulge$ppp, error_table_bulge$velocity_sd_2, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = expression('2 \u03C3'["raw 2-mom velocity"]),
        pch=16, col = tol_pal[1], xlim = c(50,500))
text(x = 455, y = 41, "Bulge", font = 2)
text(x = 440, y = 37, "Gaussian", font = 1)
abline(h=0, col="black")

bulge_v_fit = optim(par=c(44,-5,-0.2),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$velocity_sd_2,
                    method="BFGS")$par

assymptote_plot(par = bulge_v_fit)

par(fig = c(0,0.53,w[5,1],w[5,2]),
    new=T)
magplot(error_table_bulge$ppp, error_table_bulge$dispersion_sd_2, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = expression('2 \u03C3'["raw 2-mom dispersion"]),
        pch=16, col=tol_pal[1], xlim = c(50,500))
abline(h=0, col="black")

bulge_s_fit = optim(par=c(5,-0.5,-3),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$dispersion_sd_2,
                    method="BFGS")$par

assymptote_plot(par = bulge_s_fit)

par(fig = c(0,0.53,w[4,1],w[4,2]),
    new=T)
magplot(error_table_bulge$ppp, error_table_bulge$velocity_sd_4, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = expression('2 \u03C3'["raw 4-mom velocity"]),
        pch=16, col = tol_pal[1], xlim = c(50,500))
text(x = 415, y = 48, "Gauss-Hermite", font = 1)
abline(h=0, col="black")

bulge_v_fit = optim(par=c(44,-5,-0.2),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$velocity_sd_4,
                    method="BFGS")$par

assymptote_plot(par = bulge_v_fit)

par(fig = c(0,0.53,w[3,1],w[3,2]),
    new=T)
magplot(error_table_bulge$ppp, error_table_bulge$dispersion_sd_4, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = expression('2 \u03C3'["raw 4-mom dispersion"]),
        pch=16, col=tol_pal[1], xlim = c(50,500))
abline(h=0, col="black")

bulge_s_fit = optim(par=c(5,-0.5,-3),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$dispersion_sd_4,
                    method="BFGS")$par

assymptote_plot(par = bulge_s_fit)


par(fig = c(0,0.53,w[2,1],w[2,2]),
    new=T)
magplot(error_table_bulge$ppp, error_table_bulge$h3_sd,  
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = expression('2 \u03C3'["raw h3"]),
        pch=16, col=tol_pal[1], ylim = c(0.05,0.38), xlim = c(50,500))
abline(h=0, col="black")

bulge_3_fit = optim(par=c(1,-3,-0.01),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$h3_sd,
                    method="BFGS")$par

assymptote_plot(par = bulge_3_fit)

par(fig = c(0,0.53,w[1,1],w[1,2]),
    new=T)
magplot(error_table_bulge$ppp, error_table_bulge$h4_sd, 
        xlab = "Number of particles per pixel",
        ylab = expression('2 \u03C3'["raw h4"]),
        pch=16, col=tol_pal[1], xlim = c(50,500))
abline(h=0, col="black")

bulge_4_fit = optim(par=c(1,-3,-0.01),
                    fn = assymptote_fun,
                    x = ppp,
                    y = error_table_bulge$h4_sd,
                    method="BFGS")$par

assymptote_plot(par = bulge_4_fit)


par(fig = c(0.47,1,w[6,1],w[6,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$velocity_sd_2, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = "",#expression('2 \u03C3'["raw velocity"]),
        pch=16, col = tol_pal[1], xlim = c(50,500))
text(x = 460, y = 22, "Disc", font = 2)
text(x = 440, y = 20, "Gaussian", font = 1)
abline(h=0, col="black")

disk_v_fit = optim(par=c(44,-5,-0.2),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$velocity_sd_2,
                   method="BFGS")$par

assymptote_plot(par = disk_v_fit)

par(fig = c(0.47,1,w[5,1],w[5,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$dispersion_sd_2, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = "",#expression('2 \u03C3'["raw dispersion"]),
        pch=16, col=tol_pal[1], xlim = c(50,500))
abline(h=0, col="black")

disk_s_fit = optim(par=c(5,-0.5,-3),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$dispersion_sd_2,
                   method="BFGS")$par

assymptote_plot(par = disk_s_fit)

par(fig = c(0.47,1,w[4,1],w[4,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$velocity_sd_4, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = "",#expression('2 \u03C3'["raw velocity"]),
        pch=16, col = tol_pal[1], xlim = c(50,500))
text(x = 415, y = 24, "Gauss-Hermite", font = 1)
abline(h=0, col="black")

disk_v_fit = optim(par=c(44,-5,-0.2),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$velocity_sd_4,
                   method="BFGS")$par

assymptote_plot(par = disk_v_fit)

par(fig = c(0.47,1,w[3,1],w[3,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$dispersion_sd_4, 
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = "",#expression('2 \u03C3'["raw dispersion"]),
        pch=16, col=tol_pal[1], xlim = c(50,500))
abline(h=0, col="black")

disk_s_fit = optim(par=c(5,-0.5,-3),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$dispersion_sd_4,
                   method="BFGS")$par

assymptote_plot(par = disk_s_fit)


par(fig = c(0.47,1,w[2,1],w[2,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$h3_sd,  
        #xlab = "Number of particles per pixel",
        xaxt="n",
        ylab = "",#expression('2 \u03C3'["raw h3"]),
        pch=16, col=tol_pal[1], ylim = c(0.07,0.28), xlim = c(50,500))
abline(h=0, col="black")

disk_3_fit = optim(par=c(1,-3,-0.01),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$h3_sd,
                   method="BFGS")$par

assymptote_plot(par = disk_3_fit)

par(fig = c(0.47,1,w[1,1],w[1,2]),
    new=T)
magplot(error_table_disk$ppp, error_table_disk$h4_sd, 
        xlab = "Number of particles per pixel",
        ylab = "",#expression('2 \u03C3'["raw h4"]),
        pch=16, col=tol_pal[1], ylim = c(0.05,0.38), xlim = c(50,500))
abline(h=0, col="black")

disk_4_fit = optim(par=c(1,-3,-0.01),
                   fn = assymptote_fun,
                   x = ppp,
                   y = error_table_disk$h4_sd,
                   method="BFGS")$par

assymptote_plot(par = disk_4_fit)

dev.off()

