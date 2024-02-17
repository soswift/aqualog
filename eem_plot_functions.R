# functions for plotting EEMs
# ultra-simplified versions of staRdom/eemR reading and plotting
require(data.table)
require(ggplot2)
require(viridisLite)

# plot_eem() creates a binned raster plot of the eem
# em = emission wavelength
# ex = excitation wavelength
# Fluorescence = fluorescence measurement
# sample_name = character string name for plot tile

plot_eem <- function(eem, sample_name=NULL, contour = T, color_breaks = NULL, rows_as_names = F){
  
  require(data.table)
  
  # long format
  eem_dt <- as.data.table(eem, keep.rownames = rows_as_names)
  names(eem_dt)[1] <- "em"
  eem_long <- suppressWarnings(melt.data.table(eem_dt, id.vars = "em", variable.factor = F))
  names(eem_long) <- c("em","ex","Value")
  eem_long[ , ex:= as.numeric(sub("X","",ex))]
  
  # create plots with uniform scale if color breaks have been specified
  if(!is.null(color_breaks)) {
    
    ncolors = length(color_breaks[-1])
    
    eem_long[, color_group := cut(Value,
                                  breaks = color_breaks,
                                  labels = color_breaks[-1])] 
    
  # uniform scale plot
    p = ggplot(eem_long, aes(y = as.numeric(em),
                             x = ex,
                             color = color_group,
                             fill = color_group,
                             z = Value))+
      geom_tile()+
      scale_fill_manual(name = "Value",breaks = color_breaks[-1], values = viridis(ncolors), na.value = "white")+
      scale_color_manual(name = "Value",breaks = color_breaks[-1], values = viridis(ncolors), na.value = "white")+
      labs(title = paste("Sample:", sample_name),
           x = "Excitation (nm)",
           y = "Emission (nm)")+
      theme_minimal()
      theme()
    
    
 }else{
  
  # create free-scale plot
  p = ggplot(eem_long, aes(y = as.numeric(em),
                           x = ex,
                           color = Value,
                           fill = Value,
                           z = Value))+
          geom_tile()+
          scale_fill_stepsn(n.breaks = 10, colors = viridis(10), na.value = "white")+
          scale_color_stepsn(n.breaks = 10, colors = viridis(10), na.value = "white")+
          labs(title = paste("Sample:", sample_name),
               x = "Excitation (nm)",
               y = "Emission (nm)")+
          theme_minimal()
 }
        
  if(contour == T){
    p = p + geom_contour(color = "white",bins = 3)
  }
  return(p)
}

# Annotate an fdom plot with JP5 box and value
jp5_annotation <- function(eem_plot, eem_indx, sample_id = paste("Sample:", eem_indx$UniqueID, eem_indx$collection_code, eem_indx$received_from)){
  
  jp5_value = eem_indx$JP5_normalized
  
  eem_plot+
  annotate("path",
           y=c( 308, 344, 344, 308, 308),
           x=c(260, 260,290,290,260),
           color = "gray",
           alpha = 0.8)+
    labs(title = sample_id,
         subtitle = paste("Normalized JP5 = ", round(jp5_value, 3), "Empirical JP5 = ", round(eem_indx$JP5_empirical, 3)))
}

# Plot JP5 contamination 
# c_range example = c(-Inf, 0.5, 1.0, 3,   Inf)
# c_level example = contam_levels <- c("Non-Detect" = "gray", "Low (8<ppb)" = "yellow", "Moderate (40-200ppb)" = "orange", "High (> 200ppb)" = "red")
jp5_severity_plot <- function(jp5_value,
                              jp5_level,
                              c_range,
                              c_level){
  
  # make plot data structure with empty columns on either side
  df <- data.frame(jp5 = c(NA,jp5_value,NA),
                   sample = factor( c("a_null","sample","z_null")),
                   severity = c(NA,jp5_level,NA))
  
  
  ggplot(df, aes(y = jp5, x = sample, color = jp5_level)) + 
  geom_point( size = 8, pch = 17)+
  labs(x = "", y = "JP5 Normalized Fluor.")+
  scale_color_manual(name='JP5 Detection',
                       breaks=names(c_level),
                       values=c_level)+
  annotate("rect",
           ymin = c_range[1], ymax = c_range[2],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = c_level[[1]])+
  annotate("rect",
           ymin = c_range[2], ymax = c_range[3],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = c_level[[2]])+
  annotate("rect",
           ymin = c_range[3], ymax = c_range[4],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = c_level[[3]])+
  annotate("rect",
           ymin = c_range[4], ymax = c_range[5],
           xmin = "a_null", xmax = "z_null",
         alpha = 0.2,
         fill = c_level[[4]])+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
}

# plot_ex() creates a line plot of emission values for a specified excitation value
# excitation_nm = integer, specified excitation value (e.g. 266)
# eem = matrix, standard eem format with first column "X" containing emission wavelengths
# sample_name = character string name for plot tile
plot_ex = function(excitation_nm, eem, sample_name){
  
  # pull excitation and emission values from data
  exs = as.integer(sub("X","",names(eem)))
  ems = eem[[1]]
  
  
  # find closest excitation value
  ex_dif = abs(exs - excitation_nm)
  closest_ex = which(ex_dif == min(ex_dif, na.rm = T))
  
  ex_dat = data.frame(ems = ems, value = eem[,closest_ex])
  ex_dat = na.omit(ex_dat)
  
  ggplot(ex_dat, aes(x = ems, y = value))+
    geom_line()+
    labs(title = sample_name,
         subtitle = paste("Excitation =", excitation_nm, "nm"),
         x = "Emission (nm)",
         y = "Fluorescence")
}

# Testing ----
#  eem_file = "~/Documents/Bioinformatics/Projects/Redhill/data/fdom_runs/20220114/processed_data/processed_matrices/20220114_RH_100_Group002Sample0007_clean.csv"
#  indx_file = "~/Documents/Bioinformatics/Projects/Redhill/data/fdom_compiled/all_sample_indices.csv"
#  EEM = fread(eem_file, header = T)
#  indx = fread(indx_file, header = T)
#  p = plot_eem(EEM)
#  p
# 
#  jp5_annotation(p, indx)
# 
# ggsave(p, "test.pdf")
#

