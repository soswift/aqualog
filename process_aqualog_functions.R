#!/usr/bin/Rscript
  # first attempt to translate aqualog raw data processing from matlab to R
  # this script is executable from the command line. Tested on linux. 

# Set Global Parameters --------------------------
library(eemR)
library(dplyr)

# Aqualog parameters
# by default, in the aqualog output file the x-axis is emission, y-axis is excitation
# EmStart & EmEnd = estimate the range of emissions matching the number of Aqualog export rows
# ExStart & ExEnd = estimate the range of excitations matching the number of aqualog export columns
# EMstep = There is an EmStep distance between each emission value (rows). This is dictated by the number of pixels assigned to integrate 4.65nm = 8 pixels.
# EXStep = There is an ExStep distance between each excitation value (columns).

Slitwidth = 5
ExStart   = 240
ExEnd     = 500
ExStep    = 5
EmStart   = 248   
EmEnd     = 824.6 
EmStep    = 4.65  

# Based on start, end, and step, determine the Em/Ex values measured at each cell in the matrix
# Ex and Em identify the axis values for Ex/Em, example Excitation: 500nm 495nm  ... 245nm 240nm
Ex = seq(from = ExEnd, to = ExStart, by = -ExStep) 
Em = seq(from = EmStart, to = EmEnd, by = EmStep)

# Select the appropriate range of the scan to read based on the number of steps
# Typical Aqualog scans for Nelson Lab are 53 columns (Ex) and 125 rows (Em)
Ex_range <- length(Ex)
Em_range <- length(Em)

params = c(Slitwidth = Slitwidth,
           ExStart   = ExStart,
           ExEnd     = ExEnd,
           ExStep    = ExStep,
           EmStart   = EmStart,   
           EmEnd     = EmEnd,
           EmStep    = EmStep)

print("Aqualog Parameters")
writeLines(paste(names(params),params, sep = "\n"))

# Check Blanks -------------------

# check_blanks() imports blanks from each group, removes scattering, and plots the blank EEM. 
  # plots are returned as a list and are written out to a pdf file
  # this function uses the package eemR for scatter removal and plotting.

# blank checks
check_blanks = function(run_name,
                        data_directory,
                        out_dir) {
  require(eemR)
  # organize data files
  all_files = list.files(data_directory, full.names = T)
  
  BEMs = all_files[grepl("BEM", all_files)]
  
  BEM_dir = file.path(out_dir,"BEM")
  
  dir.create(BEM_dir, showWarnings = F)
  file.copy(BEMs, BEM_dir, overwrite = T)
  
  # read in BEM files (blank data)
  blanks = eem_read(file = BEM_dir,
                    recursive = T,
                    import_function = "aqualog")
  
  group_blanks = eem_extract(blanks,
                             sample = "Sample0001",
                             keep = T)
  
  print(paste("Number of Blank groups =", length(group_blanks)))
  
  # remove scattering from blanks
  blanks_no_scat = eem_remove_scattering(group_blanks, "rayleigh", 1, 11) %>%
    eem_remove_scattering("rayleigh", 2, 11) %>%
    eem_remove_scattering("raman", 1, 11) %>%
    eem_remove_scattering("raman", 2, 11)
  
  
  
  # write to PDF
  outfile = paste0("blanks_", run_name,".pdf")
  
  print(paste("writing blank plots to file", outfile))
  
  blank_plots = list()
  
  pdf(file = file.path(out_dir, outfile),
      width = 10,
      height = 8)
  for (i in 1:length(blanks_no_scat)) {
    blank_plots[[i]] = plot(blanks_no_scat, which = i)
  }
  dev.off()
  
  # return plots as list
  return(blank_plots)
}

# Data file organization --------------
  
# organize_aqualog() organizes and summarizes data from an aqualog run.
  #Assumes default aqualog naming scheme. 
  #Automates association of samples, groups, filepaths, and blanks

  # data_directory = Identify the input directory containing the raw aqualog files
  # run_name = Provide a descriptive run name that will be prepended to files
  # sample_key_file = Provide a .tsv file with column "UniqueID" with samples listed in the order they were loaded

  # In this section, relevant fdom data files are identified and organized
  # Assumptions:
  # 1) all of the data files are in a single directory
  # 2) Blanks can be identified by *.ogw files
  # 3) Within a sample group, the first sample's BEM can be used as a blank for all samples in the group
  # 4) For each group, there are two blanks: 1) blank run before the group; 2) blank run after the group
  # 5) For the last group, the two blanks are a little different: 1) blank run before the group; 2) the blank from the first group

organize_aqualog <- function(data_directory,
                             run_name,
                             sample_key_file,
                             out_dir = data_directory){  
    
    # aqualog outputs data files (*.dat) and blank files (*.ogw)
    # *.dat files include raw data BEM, SEM, and ABS
    writeLines(paste0("Reading data files in directory: ", data_directory))
    
    all_files <- list.files(data_directory)
    dat_files <- grep("\\.dat", all_files, value = T)
    blank_files <- grep("Sample0001BEM", dat_files, value = T)
    SEM_files <- grep("SEM", dat_files, value = T)
    ABS_files <- grep("ABS", dat_files, value = T)
  
    # summarize files and samples
    samples      <- gsub("(Sample\\d+).+", "\\1", SEM_files, perl = T)
    sample_names <- unique(samples)
   
    N_samples    <- length(sample_names)
    N_blanks     <- length(blank_files)
  
    if(N_samples < 1)  stop("No sample data files (*.SEM) in the directory!") 
    if(N_blanks  < 1)  stop("No blank files (*.ogw) in the directory!")
      
    groups       <- gsub("(Group\\d+).+", "\\1", sample_names, perl = T)
    group_names  <- unique(groups)
    N_groups     <- length(group_names)

    samples_per_group <- sapply(group_names, function(x) length(groups[groups == x]))
    typical_group_size <- max(samples_per_group)
    last_group_size    <- tail(samples_per_group, 1)
   
    if(any(samples_per_group > 7)){
      warning("WARNING: there are a lot of samples per blank, maybe double check that"); print(samples_per_group)
    }
    
    if(N_groups != N_blanks){
       warning("WARNING: The number of groups does not match the number of blanks!")
    }
    

    # match data files to sample names
    writeLines("Reading in sample key")
    sample_key <- read.table( file.path(data_directory,sample_key_file),
                             stringsAsFactors = F,
                             header = T,
                             fill = T,
                             sep = "\t")
    
    load_order <- sample_key$UniqueID
    
    if(!("UniqueID" %in% colnames(sample_key))){
      stop("ERROR Column 'UniqueID' missing from sample key")
    }
    
    key_blanks <- grepl("blank", load_order, ignore.case = T)

    if(any(key_blanks)){
      writeLines(paste("Blanks detected in sample key =", sum(key_blanks)))
      load_order <- load_order[!grepl("blank", load_order, ignore.case = T)]
    }

    if(N_samples != length(load_order)){
       stop(paste("ERROR Sample key mismatch",
       "samples in key =", length(load_order),
       "sample data files =", N_samples))
    }
    
    # print summary 
    writeLines(paste0("Total N samples = ", N_samples , "\n",
                      "Total N groups = ", N_groups, "\n",
                      "Total N blanks  = ", N_blanks, "\n",
                      "N samples in each group = ", typical_group_size, "\n",
                      "N samples in last group = ", last_group_size
    )
    )
  
  # Assign blanks for each group
    # for blank1, use the blank corresponding to the group name
    # for blank2, use the blank of the next group in a circular fashion, so the last group gets the first as blank2
    blank1 <- sapply(group_names, function(x) grep(x, blank_files, value = T))
    blank2 <- c(tail(blank1, -1), head(blank1, 1))
    names(blank2) <- group_names
    
  # compile sample file summary into a spreadsheet
    sample_sheet <- data.frame(
                          UniqueID = load_order,
                          run_name = run_name,
                          run_order = 1:N_samples,
                          sample = sample_names,
                          group  = groups,
                          SEM_file = SEM_files,
                          blank1 = blank1[groups],
                          blank2 = blank2[groups],
                          stringsAsFactors = F
    )
   # if there are absorbance files, add that information, otherwise skip
    ABS_collected <- FALSE
    
    if(length(ABS_files) > 1){
      ABS_collected <- TRUE
      writeLines("Absorbance files detected, including absorbance data and IFE correction")
      sample_sheet$ABS_file <- ABS_files
    } else( warning( "No absorbance data files (*ABS.dat) found in directory, skipping absorbance analysis" ) )
    

  # set up output directory 'processed_data' and write out the complied sample information
    proc_dir <- file.path(out_dir,"processed_data")
    if(!dir.exists(proc_dir)){
      dir.create(proc_dir)
    }

    # write out sample sheet to a CSV file
    sample_sheet_file = paste0(run_name, "_sample_sheet.csv")
    writeLines(paste0("writing out sample_sheet as: ", sample_sheet_file))
    write.csv(sample_sheet, file.path(out_dir,"processed_data",sample_sheet_file), row.names = F)
    
    return(sample_sheet)
}

# Correct EEMs ------------------
# correct_EEMs() performs corrections on the raw fluorescence data
  # Read in the excitation/emission matrices for each sample in a run.
  # If correcting multiple runs, make sure thye're in separate directories.
  # Relies on information in sample sheet from organize_aqualog().
  # Alternatively, provide a manually generated sample sheet.
  # see organize_aqualog() function for more details.

  # specify a sample sheet and a run name
 
 # required sample sheet columns:
  # UniqueID = character, identifies each unique sample
  # run_order = intiger, order in which sample was run (1: N samples)
  # sample = character, aqualog data file name(e.g. Group001Sample001)
  # group = chracter, blank group for sample (e.g. Group001)
  # SEM_file = character, filepath of aqualog data file (e.g. Group001Sample0001SEM.dat)
  # blank1 = chracter, filepath of first blank for group (e.g. Group001Sample0001BEM.dat)
  # blank2 = chracter, filepath of first blankf for group (e.g. Group002Sample0001BEM.dat)


correct_EEMs <- function(data_directory, 
                         sample_sheet,
                         run_name,
                         out_dir = data_directory,
                         blank_subtract = T){
  
  # SEM_files = filepaths for raw Aqualog SEM files
  # ABS_files = filepaths of ABS files,
  # blank_files = filepaths for blanks
  # groups = group id for each SEM file
  
  # TODO: process each blank once, associate blanks w samples
    
  # pull information from sample sheet
    groups = sample_sheet$group
    
    SEM_files = file.path(data_directory, as.character(sample_sheet$SEM_file))
    ABS_files = file.path(data_directory, as.character(sample_sheet$ABS_file))

    blank1 = file.path(data_directory, as.character(sample_sheet$blank1))
    blank2 = file.path(data_directory, as.character(sample_sheet$blank2))
    
    blank1_names = as.character(sample_sheet$blank1)
    blank2_names = as.character(sample_sheet$blank2)
    sample_names = sample_sheet$sample
    UniqueID = sample_sheet$UniqueID

    N_samples = nrow(sample_sheet)
    
  
  # get SEM data from files as matrix
    writeLines("extracting EEM data from files")
    if(all(file.exists(SEM_files))){
      print("Good. All SEM files exist.")
    }else(stop("Not all SEM files exist!"))

    read_SEM <- function(sem_file){
      sem_file <- as.character(sem_file)
      if(!file.exists(sem_file)) stop(paste("File not readable: ", sem_file))
      # read specified region of the emission/excitation matrix
      as.matrix(read.table(sem_file, header = T, row.names = 1)[1:Em_range, 1:Ex_range])
    }
    
    Sample_EEMs <- lapply(SEM_files, read_SEM) # sample EEM data
    Blank1_EEMs <- lapply(blank1, read_SEM) # blank1 EEM data
    Blank2_EEMs <- lapply(blank2, read_SEM) # blank2 EEM data
    
    # check blanks
    
    
  
  # 1. Read in absorbance data ------
    # If there was absorbance data recorded, read the ABS files
    # Note: absorbance measurements are used to correct excitation/emission values
    # The equation and calculations to perform the correction are detailed below
  
    if(!is.null(ABS_files)){
    writeLines("Absorbance (ABS) files detected, reading absorbance data")
    # absorbance values are collected at the same wavelengths as excitation
    Abs = Ex
    Ab_range <- length(Abs)
    
    # read in specified region of the absorbance matrix
    read_ABS <- function(abs_file){
      abs_table <- read.table(abs_file, row.names = 1)
      abs_vec <- abs_table[[1]]
      names(abs_vec) <- row.names(abs_table)
      return(abs_vec)
    }
  
   Sample_ABSs <- lapply(ABS_files, read_ABS)
     
    # Generate a correction factor for each EEM value (Kothwala et al. 2013, 10.4319/lom.2013.11.616, equation 2)
    # this is the 'Inner Filter Effect' correction
    # To perform the correction, identify the absorbance value corresponding to each emission and excitation wavelength
    # For excitation, this is easy because absorbance is measured at the same wavelengths as excitation
    # To pair up absorbance readings with emissions, find the nearest absorbance measurement that is < the emission wavelength
  
    Abs_equals_Em  <- sapply(Em, function(x) as.character(Abs[Abs < x][1]))
    Abs_equals_Ex  <- as.character(Ex)
    
    # define formula for calculating correction at each combination of emission and excitation
    # x = the absorption value at a given emission wavelength
    # y = the absorption value at a given excitation wavelength
  
    ABA_IFE_corection <- function(x, y) 10^(0.5 *(x + y))
    
    
    # for each sample:
    # calculate a corrected fluorescence for each pair of Excitation/Emission values
    Sample_IFE_cors <- lapply(Sample_ABSs, function(abs_vec){
    
      x = abs_vec[Abs_equals_Em]
      y = abs_vec[Abs_equals_Ex]
      
      IFE_cor <- outer(x, y, ABA_IFE_corection)
      return(IFE_cor)
    })
    } else writeLines("No absorbance (.ABS) files detected, skipping absorbance correction")
  
  # 2. Perform EEM corrections -------------------------------
    # Each sample's EEM matrix needs to be corrected in order to remove known sources of noise
    # the corrections are:
    # 1) Inner Filter Correction
    # 2) Raman scaling
    # 3) Blank subtraction
    # 4) Rayleigh scatter removal
    
    # Corrections follow 
    # Kothwala et al. 2013, DOI: 10.4319/lom.2013.11.616
    writeLines("calculating Raman Areas")
    
    # Calculate Raman areas
    # Kothawala et al. 2013 pg 619 'Fluorescence measurements'
    # Raman area of pure water integrated over emission range of 381 nm - 426 nm at excitation 350 nm
    # There might not be measurements at exactly 381nm and 426nm, so take the conservative approach:
    # allow the range of emission values to expand a little bit beyond the precise limits, just to be safe
    
    # define ranges to correct
    Ram_em_range <- list(low  = tail(Em[Em<=380], 1),
    	  	               high = head(Em[Em>=426], 1))
    
    Ram_em_rows <- ( Em >= Ram_em_range$low &
  		               Em <= Ram_em_range$high ) 
    
    Ram_ex <- 350
    Ram_ex_col <- which(Ex == Ram_ex)
    
    # function to integrate over Raman area of an EEM matrix
    Raman_integration <- function(blank_avg){
      EmStep * sum(blank_avg[ Ram_em_rows, Ram_ex_col ])
    }
    
    # average blanks 1 and 2 for each group and apply integration function
    BlankAverages <- lapply(1:length(Blank1_EEMs), 
                            function(i) (Blank1_EEMs[[i]] + Blank2_EEMs[[i]]) / 2 ) 
    

    BlankRamanAreas <- lapply(BlankAverages, Raman_integration)
  
    # a) Apply Inner Filter Correction
    # multiply sample EEM values by the IFE correction values calculated from absorbances (see section 3. above)
    if(!is.null(ABS_files)){
    
       writeLines("performing Inner Filter Correction")
       Sample_EEMs_IFE <- lapply(1:N_samples,
                               function(i) Sample_EEMs[[i]] * Sample_IFE_cors[[i]] )
    }else{ 

	writeLines("skipping Inner Filter Correction")
        Sample_EEMs_IFE = Sample_EEMs
    }
    
    # b) Apply Raman Scaling
    writeLines("performing Raman scaling")
    # use the Raman area derived from the appropriate pair of blanks for the sample
    Sample_EEMs_Ram <- lapply(1:N_samples,
                              function(i) Sample_EEMs_IFE[[i]] / BlankRamanAreas[[i]])
    # also correct the blanks
    Blank1_EEMs_Ram <- lapply(1:N_samples,
                         function(i) Blank1_EEMs[[i]] / BlankRamanAreas[[i]])
    Blank2_EEMs_Ram <- lapply(1:N_samples,
                         function(i) Blank2_EEMs[[i]] / BlankRamanAreas[[i]])
  
    # c) Subtract Blanks (Raman Scaling of blanks by average)
    if(isTRUE(blank_subtract)){
      # divide the averaged blanks by their Raman area, then subtract that from the samples in each group
      BlankAverages_Ram <- lapply(1:N_samples,
                                  function(i) BlankAverages[[i]] / BlankRamanAreas[[i]] )
      # subtract the Raman scaled blanks from the Raman scaled samples in the appropriate group
      Sample_EEMs_BlkSub <- lapply(1:N_samples,
                                   function(i) Sample_EEMs_Ram[[i]] - BlankAverages_Ram[[i]])
      
      Blank1_BlkSub <- lapply(1:N_samples,
                              function(i) Blank1_EEMs_Ram[[i]] - Blank2_EEMs_Ram[[i]])
      Blank2_BlkSub <- lapply(1:N_samples,
                              function(i) Blank2_EEMs_Ram[[i]] - Blank1_EEMs_Ram[[i]])
    }else{
      # skip blank subtraction (e.g. for diagnosing blank issues)
      Sample_EEMs_BlkSub <- Sample_EEMs_Ram
      Blank1_BlkSub <- Blank1_EEMs_Ram
      Blank2_BlkSub <- Blank2_EEMs_Ram
      
    }
    # d) Remove Rayleigh scatter
    writeLines("performing Rayleigh scatter subtraction")
    
    #subtract Rayleigh scatter where excitation is > emission
    # remove the rayleigh scatter by cutting based on slit widths
    w = Em[1] - Ex - (2 * Slitwidth)
    s = w * -1 / EmStep
    
    s_remove = s[w <0]
    
    # function to remove Rayleigh scattering from an EEM matrix
    Remove_Rayleigh <- function(eem){
      # for each extraction wavelength, if w is negative, subtract all values less than s
      for(i in which(w < 0)){
        eem[ 1:s[i] , i] <- 0
      }
      return(eem)
    }
    
    # remove Rayleigh scattering from samples and blanks
    Sample_EEMs_Ray <- lapply(Sample_EEMs_BlkSub, Remove_Rayleigh)
    Blank1_Ray <- lapply(Blank1_BlkSub, Remove_Rayleigh)
    Blank2_Ray <- lapply(Blank2_BlkSub, Remove_Rayleigh)
    
    # e) Crop out scattering regions.
    
    # define cropping functions:
    # Calc_Crops() generates excitation values that should be cropped to across a range of emissions values.
    # em_range = the range of emission values to be cropped out of the matrix
    # ex_range = the range of excitation vales to be cropped out of the matrix
    # together, these two ranges define a line
    # the triangle of values that fall above or below the line will be cropped out by Crop_EEM()
   
    Calculate_Crop <- function(em_range, ex_range){
      # model the crop line
      crop_line <- lm(ex_range ~ em_range)
      
      # calculate the excitation cutoff point for each emission value 
      ex_crops <- coef(crop_line)[2] * Em + coef(crop_line)[1]
      return(ex_crops)
    }
      
    # now that the crop values are defined, pass them on to Crop_EEM()
    # specify direction to crop:
       # "left" =  crop excitation values that less than the line
       # "right" = crop excitation values that are greater than the line
    
    Crop_EEM <- function(eem,
                         ex_crops,
                         crop_direction = c("left","right")){
      
      if(!crop_direction %in% c("left","right")){
        stop("Crop direction not specified")
      }
  
      if (crop_direction == "left") {
        # for each emission value, check if excitation values are below the crop line
        # if they are, set to NA
        for (i in 1:nrow(eem)) {
          eem[i, Ex <= ex_crops[i]] <- NA
        }
      }
      
      if (crop_direction == "right") {
        # for each emission value, check if excitation values are above the crop line
        # if they are, set to NA
        for (i in 1:nrow(eem))
          eem[i, Ex >= ex_crops[i]] <- NA
      } 
      return(eem)
    }
    
    # Crop1: crop the area corresponding to the upper left of EEM plots
    top_left <- Calculate_Crop(em_range = c(450, EmEnd),
                                 ex_range = c(ExStart, 425)) 
    
    Sample_EEMs_Crop1 <- lapply(Sample_EEMs_Ray,
                                Crop_EEM,
                                ex_crops = top_left,
                                crop_direction = "left")
    
    # Crop2: crop the area corresponding to the lower right of EEM plots
    bottom_right <- Calculate_Crop(em_range = c(255, 520),
                                  ex_range = c(ExStart, ExEnd))
    
    Sample_EEMs_Crop2 <- lapply(Sample_EEMs_Crop1,
                                Crop_EEM,
                                ex_crops = bottom_right,
                                crop_direction = "right")
    
  # 3. Write out processed sample EEM matrices -----------------------------
    # write out the proccessed data as .csv files
    
    # set negative values in the corrected EEM to zero
    Sample_EEMs_NoNeg<- lapply(Sample_EEMs_Crop2, function(eem){
      eem[eem < 0] <- 0
      return(eem)
    })
    
    # clean names of excitation spectra (X500 -> 500)
    Clean_Colnames <- function(eem){
      colnames(eem) <- sub("X","",colnames(eem))
      return(eem)
    }
    
    Sample_EEMs_Final <- lapply(Sample_EEMs_NoNeg, Clean_Colnames)
    Blank1_Final <- lapply(Blank1_Ray, Clean_Colnames)
    Blank2_Final <- lapply(Blank2_Ray, Clean_Colnames)
    
    names(Sample_EEMs_Final) <- paste(UniqueID, sample_names, sep = "_")
    names(Blank1_Final) <- blank1
    names(Blank2_Final) <- blank2
    
    
    # write out
    proc_dir   = file.path(out_dir, "processed_data")
    out_path   = file.path(proc_dir, "processed_matrices")
    blank_path = file.path(proc_dir,"blanks")
    
    if(!dir.exists(out_path)){
      dir.create(out_path, recursive = T)
    }
    
    if(!dir.exists(blank_path)){
      dir.create(blank_path)
    }
    
    # blanks
    blank1_output_files <- file.path(blank_path, paste(run_name, blank1_names, "clean.csv",sep="_" ))
    
    invisible( lapply(1:length(Blank1_Final),
                      function(i)   write.csv(Blank1_Final[[i]], blank1_output_files[i])))
    writeLines(unique(paste(blank1_names, round(sapply(Blank1_Final, sum)))), con = file.path(blank_path,"blank_sums.txt"))
    
    
    # samples
    sample_output_files <- file.path(out_path, paste(run_name, names(Sample_EEMs_Final), "clean.csv", sep = "_"))

    writeLines("writing out cleaned EEM files to directory processed_data > processed_matrices")
    saveRDS(Sample_EEMs_Final,
            file.path(out_dir,"processed_data", paste0(run_name, "_EEM.rds")))
    
    invisible( lapply(1:N_samples,
           		function(i)   write.csv(Sample_EEMs_Final[[i]], sample_output_files[i])))
    
    
    
  return(Sample_EEMs_Final)
}


# Calculate indices ----------------------------------
# This function takes a list of cleaned sample EEMs and calculates various indices for each.
# These indices include published, single point indices (E.g. Coble A,B,C,M,T).
# They also include custom indices based on empirical data, such as JP5.
# It's easy to add more indices. Follow the format of previous indices and change em/ex values.
# The data is technically discrete measurements at exictation/emission steps.

fdom_indices <- function(Sample_EEMs,
                         run_name,
                         sample_sheet,
                         out_dir = data_dir){
  
  # using the cleaned EEM matrices, calculate indices for various compounds/groups of compounds
  writeLines("calculating indices")
  
  # define a function to find the emission wavelenthg measured in the data that is closest
  # to a specified emission value
  Get_Nearest_Em <- function(Em_val){
    Em_dif <- abs(Em - Em_val)
    which(Em_dif == min(Em_dif))
  }
  
  # for the purpose of summing indices, set cropped regions to 0
  Sample_EEMs <- lapply(Sample_EEMs, function(eem){
      eem[is.na(eem)] <- 0
      return(eem)
      }
    )
  
  
  # List emission and excitation values of each index
  # simple indeces can be defined as a vector of c(emission value, excitation value)
  # more complicated indices relate multiple emission/excitation values to each other

   
    # simple indices involving a single peak
    # vectors of c(emission, excitation)
    simple_indices <- list(
                  CobleA = c(450,320),
                  CobleB = c(305,275),
                  CobleC = c(445, 345),
                  CobleM = c(410, 310),
                  CobleT = c(340,275),
                  
                  Fpeak = c(299,240),
                  Stedmon_D = c(509,390),
                  Optical_Brighteners = c(360, 435),
                  dieselBandII = c(510,410),
                  Petroleum = c(510,270),
                  Lignin = c(360,240)
    )
    
    indices_out <- list()
    
  # Calculate simple indices for each sample
    indices_out[["simple"]] <- sapply(simple_indices,
                                      function(index) {
                                        nearest_em <- Get_Nearest_Em(index[1])
                                        exact_ex <- which(Ex == index[2])
                                        
                                        # generate the index for each sample
                                        sapply(Sample_EEMs,
                                               function(sample_eem)
                                                 sample_eem[nearest_em, exact_ex])
                                      })

  # Calculate complex indices for each sample
    # these indices involve ranges of emissions and/or ratios of emissions
    
    # BIX, divide 380 nm by 430 nm emission @ 310 nm excitation
    indices_out[["BIX"]] <- sapply(Sample_EEMs,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(380, 430), Get_Nearest_Em)
                                     exact_ex <- which(Ex == 310)
                                     
                                     sample_eem[nearest_ems[[1]], exact_ex] / sample_eem[nearest_ems[[2]], exact_ex]
                                   })
    
    # HIX, divide the sum of 434 nm - 480 nm emissions by the sum of 300 nm - 346 nm emission @ 255 nm excitation
    indices_out[["HIX"]] <- sapply(Sample_EEMs,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(434, 480, 300, 346), Get_Nearest_Em)
                                     exact_ex <- which(Ex == 255)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex]) /
                                       sum(sample_eem[nearest_ems[[3]]:nearest_ems[[4]], exact_ex])
                                   })
    
    # FI, divide 470 nm by 520 nm emission @ 370 nm excitation
    indices_out[["FI"]] <- sapply(Sample_EEMs,
                                  function(sample_eem) {
                                    nearest_ems <- lapply(c(470, 520), Get_Nearest_Em)
                                    exact_ex <- which(Ex == 310)
                                    
                                    sample_eem[nearest_ems[[1]], exact_ex] / sample_eem[nearest_ems[[2]], exact_ex]
                                  })
    
    # M to C ratio
    indices_out[["M_to_C"]] <- indices_out[["simple"]][ , "CobleC"] / 
                                indices_out[["simple"]][ , "CobleM"]
    
    # JP5 Jet fuel, based on ordination of original community samples. This module identifies contaminated samples.
    indices_out[["JP5_empirical"]] <- sapply(Sample_EEMs,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(308, 344), Get_Nearest_Em)
                                     exact_ex <- which(Ex %in% 260:290)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex])
                                   })
    
    # JP5 Jet fuel, normalized to background fluorescence
    # integration of the range Ex260-290 and Em308-344 inclusive AFTER cropping divided by the remaining sum of the scan
    indices_out[["JP5_normalized"]] <- sapply(Sample_EEMs,
                                             function(sample_eem) {
                                               nearest_ems <- lapply(c(308, 344), Get_Nearest_Em)
                                               exact_ex <- which(Ex %in% 260:290)
                                               
                                               jp5 <- sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex])
                                               not_jp5 <- sum(sample_eem) - jp5
                                               return(jp5/not_jp5)
                                            
                                             })
    
    # Diesel, fairly wide range (Chen et al 2021 0.2166/aqua.2021.120)
    indices_out[["Diesel_Chen"]] <- sapply(Sample_EEMs,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(280, 445), Get_Nearest_Em)
                                     exact_ex <- which(Ex %in% 220:275)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex])
                                   })

    writeLines(paste("index calculated: ", c(names(simple_indices), names(indices_out)[-1]), collapse = "\n"))

    # write out calculated indices as a CSV file
    writeLines("writing out indices to processed_data > fDOM_indices_out.csv")
    indices_dat <- as.data.frame(do.call("cbind", indices_out))
    indices_dat$sample   <- gsub(".+(Group.+)","\\1",row.names(indices_dat))
    indices_dat$run_name <- rep(run_name, nrow(indices_dat))
    indices_dat$UniqueID <- sample_sheet$UniqueID[match(indices_dat$sample, sample_sheet$sample)]
    
    out_file <- paste0(run_name,"_fDOM_indices_out.csv")
    write.csv(indices_dat, file.path(out_dir, "processed_data", out_file), row.names = F)
  
    writeLines(paste0("DONE: indices written out to ", out_file))
    return(indices_dat)
}

# Process Aqualog -----------
# Main function process_aqualog() wraps up the three modules above
# If a sample organization sheet has already been created that indicates filenames, blank files, etc. 
# that custom ordering can be indicated using the org_file variable.
process_aqualog <- function(data_directory,
                            run_name,
                            org_file = NULL,
                            sample_key_file,
                            out_dir = data_directory){
  
  # start log
  logtime <- format(Sys.time(), "%Y%m%d_%H:%M:%S")
  log_file = file.path(out_dir, paste0(run_name,"_",logtime,".log"))
  
  sink(log_file)
  on.exit(sink(file = NULL))
  
  # generate file organization map if needed
  if(is.null(org_file)){
    org <- organize_aqualog(data_directory = data_directory,
                          run_name = run_name,
                          sample_key_file = sample_key_file,
                          out_dir = out_dir)
  }else{
    org <- read.csv(file.path(data_directory, org_file), header = T, stringsAsFactors = F)
  }
  
  # process eems
  corrected_eems <- correct_EEMs(data_directory = data_directory,
                                 sample_sheet = org,
                                 run_name = run_name,
                                 out_dir = out_dir)
  # calculate indices
  inds <- fdom_indices(corrected_eems,
                       sample_sheet = org,
                       run_name = run_name,
                       out_dir = out_dir)
  # check blanks
  blanks <- check_blanks(run_name = run_name,
                         data_directory = data_directory,
                         out_dir = out_dir)
                         
                         
  sink()
  
  print(readLines(log_file))
  
  return(list(sample_sheet = org, 
              EEMs = corrected_eems,
              indices = inds,
              blanks = blanks))
}
  
# Compile data across runs --------------------------

# helper function for reading flat eem matrices from directory
# by default, looks for naming pattern "clean.csv" to identify eem files
# assumes file is a csv with a header
# returns a list of matrices named with the filenames

read_eems = function(eem_dir, pattern = "clean.csv"){
  eem_files = list.files(eem_dir)
  eem_files = eem_files[grepl(pattern,eem_files)]
  
  run_eems = lapply(eem_files, function(x){
    read.csv(file.path(eem_dir, x), header = T, stringsAsFactors = F)
  })
  
  names(run_eems) = eem_files
  return(run_eems)
}


# read_processed_data() reads in the processed data from a run.
# it also double checks that run sheets match the data files.

read_processed_data = function(run_out_dir){
  
  # sample sheet and indices
  proc_data = file.path(run_out_dir, "processed_data/")
  
  sample_file  = grep("sample_sheet.csv", list.files(proc_data), value = T)
  indices_file = grep( "indices_out.csv", list.files(proc_data), value = T)
  
  sample_path  = file.path(proc_data, sample_file)
  indices_path = file.path(proc_data, indices_file)
  
  
  if(!(length(indices_path) == 1)) stop(paste("ERROR check for single indices file", run_out_dir))
  
  run_sheet   = read.csv(sample_path, header = T, stringsAsFactors = F)
  run_indices = read.csv(indices_path, header = T, stringsAsFactors = F)
  
  if(!any(run_indices$UniqueID %in% sample_sheet$UniqueID)){
    stop(
      paste("ERROR Samples in indices don't match run sample sheet", run_out_dir),
      " Samples: ",
      paste(run_indices$UniqueID[!(run_indices$UniqueID %in% sample_sheet$nelson_lab_id)])
    )
  }

  
  # sample EEMs
  run_eems = read_eems(eem_dir = file.path(proc_data, "processed_matrices"))
  eem_samples = gsub(".+(Group.+)_.+", "\\1", names(run_eems))
  eem_uniqID  = run_sheet$UniqueID[match(eem_samples, run_sheet$sample)]
  names(run_eems) = eem_uniqID
  
  # processed blank EEMs
  run_blanks = read_eems(eem_dir = file.path(proc_data,"blanks"))
  
  # run checks
  n_samples = list(sheet = nrow(run_sheet),
                   indices = nrow(run_indices),
                   eems = length(run_eems))
  
  writeLines(paste( "Samples in", 
                    "sheet = ", n_samples$sheet,
                    "indices = ", n_samples$indices,
                    "eems = ", n_samples$eems,
                    "for run", unique(run_sheet$run_name)))
  
  if(length(unique(n_samples)) != 1) stop(
    paste("Sample N mismatch in run ",
          run_out_dir,
          paste(names(n_samples),
                n_samples,
                collapse = " "))
    )
  
  
  #output sample sheet, indices, eems
  out = list(samples = run_sheet,
             indices = run_indices,
             eems = run_eems,
             blanks = run_blanks)
  return(out)
}

# compile_runs() compiles processed data across multiple runs
# it locates file for each run and calls read_processed_data()
# indices and sample metadata are merged into a single data.frame
# EEMs are merged into a concatenated matrix, which is saved as a csv
# EEMS are also as separate matrices in an R list.

compile_runs <- function(run_dirs,
                         out_dir,
                         data_dir){
  
  processed_runs = lapply(run_dirs, read_processed_data)
  
  
  # indices
  run_indices <- lapply(processed_runs, function(x) x$indices)
  compiled_indices  = do.call("rbind", run_indices)
  
  dupd_sams = duplicated(compiled_indices$UniqueID)
  if(any(dupd_sams)){
    stop( paste("ERROR duplicated UniqueID:",
                compiled_indices$UniqueID[dupd_sams],
                compiled_indices$run_name)
    )
  }
  
  # sample metadata
  all_samples = merge(sample_sheet, compiled_indices, by = "UniqueID", all.y = T,)
  if(nrow(compiled_indices) != nrow(all_samples)) stop("ERROR sample indices merge error")
  
  # EEM matrices
  run_eems = lapply(processed_runs, function(x) x$eems)
  run_eems = unlist(run_eems, recursive = F)
  all_eems = do.call("rbind", run_eems)
  all_eems$UniqueID = gsub("\\..+","",row.names(all_eems))
  
  # blanks
  run_blanks = lapply(processed_runs, function(x) x$blanks)
  run_blanks = unlist(run_blanks, recursive = F)
  
  blank_sums = sapply(run_dirs, function(x) read.table(file.path(x,"processed_data","blanks","blank_sums.txt")))
  
  # write out
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  write.csv(all_samples, file.path(out_dir,"all_sample_indices.csv"), row.names = F)
  write.csv(all_eems, file.path(out_dir,"all_sample_processed_EEMs.csv"), row.names = F)
  write.csv(blank_sums, file.path(out_dir,"all_blank_sums.csv"))
             
  saveRDS(run_eems, file.path(out_dir,"EEMs.rds"))
  saveRDS(run_blanks, file.path(out_dir, "blanks.rds"))

  return(list(indices = all_samples,
              eems = run_eems))
}


