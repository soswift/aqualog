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

plot_eem <- function(eem, sample_name=NULL, contour = T){
  
  require(data.table)
  
  # long format
  eem_dt <- as.data.table(eem, keep.rownames=F)
  names(eem_dt)[1] <- "em"
  eem_long <- melt.data.table(eem_dt, id.vars = "em", variable.factor = F)
  names(eem_long) <- c("em","ex","Value")
  eem_long[ , ex:= as.numeric(sub("X","",ex))]
  
  # plot
  p = ggplot(eem_long, aes(y = em,
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
        
  if(contour == T){
    p = p + geom_contour(color = "white",bins = 3)
  }
  return(p)
}

# Annotate an fdom plot with JP5 box and value
jp5_annotation <- function(eem_plot, eem_indx){
  
  jp5_value = eem_indx$JP5_normalized
  
  eem_plot+
  annotate("path",
           y=c( 308, 344, 344, 308, 308),
           x=c(260, 260,290,290,260),
           color = "gray",
           alpha = 0.8)+
    labs(title = paste("Sample:", eem_indx$UniqueID, eem_indx$collection_code, eem_indx$received_from),
         subtitle = paste("Normalized JP5 = ", round(jp5_value, 3), "Empirical JP5 = ", round(eem_indx$JP5_empirical, 3)))
}

# Plot JP5 contamination 
jp5_severity_plot <- function(jp5_value,
                              jp5_level,
                              c_range = contam_ranges,
                              c_level = contam_levels){
  
  # make plot data structure with empty columns on either side
  df <- data.frame(jp5 = c(NA,jp5_value,NA),
                   sample = factor( c("a_null","sample","z_null")),
                   severity = c(NA,jp5_level,NA))
  
  
  ggplot(df, aes(y = jp5, x = sample, color = jp5_level)) + 
  geom_point( size = 8, pch = 17)+
  labs(x = "", y = "JP5 Normalized Fluor.")+
  scale_color_manual(name='JP5 Detection',
                       breaks=names(contam_levels),
                       values=contam_levels)+
  annotate("rect",
           ymin = contam_ranges[1], ymax = contam_ranges[2],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = contam_levels[[1]])+
  annotate("rect",
           ymin = contam_ranges[2], ymax = contam_ranges[3],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = contam_levels[[2]])+
  annotate("rect",
           ymin = contam_ranges[3], ymax = contam_ranges[4],
           xmin = "a_null", xmax = "z_null",
           alpha = 0.2,
           fill = contam_levels[[3]])+
  annotate("rect",
           ymin = contam_ranges[4], ymax = 0.4,
           xmin = "a_null", xmax = "z_null",
         alpha = 0.2,
         fill = contam_levels[[4]])+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
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

