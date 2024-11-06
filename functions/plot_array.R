plot_array = function(plot_size, n_plots){
  bounds_array = matrix(data = NA, nrow = n_plots, ncol = 2)
  bounds_array[1,] = c(0, plot_size)
  bounds_array[n_plots,] = c(1-plot_size, 1)
  
  overlap = ((plot_size * n_plots) - 1)/(n_plots - 1)
  if (n_plots > 2){end = (n_plots-1)} else {end = 2}
  for (i in 2:end){
    bounds_array[i,] = c((bounds_array[i-1,2] - overlap), (bounds_array[i-1,2] - overlap) + plot_size)
  }
  return(bounds_array)
}
