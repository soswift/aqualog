#!/usr/bin/Rscript
  # first attempt to translate aqualog raw data processing from matlab to R
  # this script is executable from the command line. Tested on linux. 


  # 1. Organize the data files --------------------
    # In this section, relevant fdom data files are identified and organized
    # Assumptions:
    # 1) all of the data files are in a single directory
    # 2) Blanks can be identified by *.ogw files
    # 3) Within a sample group, the first sample's BEM can be used as a blank for all samples in the group
    # 4) For each group, there are two blanks: 1) blank run before the group; 2) blank run after the group
    # 5) For the last group, the two blanks are a little different: 1) blank run before the group; 2) the blank from the first groupi
  
    # parse arguments
    library(optparse)
  
    option_list <- list(
    make_option(c("-i", "--input_dir"),
                action="store", default=NA, type = 'character',
                help="Identify the input directory containing the raw aqualog files"),

    make_option(c("-d", "--description"),
                action="store", default=NA, type='character',
			          help="Provide a descriptive run name that will be append to files"),
    
    make_option(c("-k", "--key_to_samples"),
                action="store", default=NA, type = 'character',
                help = "Provide a .tsv file with column 'UniqueID' indicating order in which samples were loaded")
    )
    
    opt = parse_args(OptionParser(option_list=option_list))
    if(any(is.na(opt))){
      stop(paste0("ERROR missing argument: ", names(opt[is.na(opt)])))
    }
    
    # organize command line inputs
    data_directory = opt$input_dir
    run_name = opt$description
    sample_key_file = opt$key_to_samples
    
    # set working directory
    setwd(data_directory)
    
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
    sample_key <- read.table(sample_key_file,
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
                          run_order = 1:N_samples,
                          sample = sample_names,
                          group  = groups,
                          SEM_file = SEM_files,
                          blank1 = blank1[groups],
                          blank2 = blank2[groups]
    )

   # if there are absorbance files, add that information in
    ABS_collected <- FALSE
    
    if(length(ABS_files) > 1){
      ABS_collected <- TRUE
      writeLines("Absorbance files detected, including absorbance data and IFE correction")
      sample_sheet$ABS_file <- ABS_files
    } else( warning( "No absorbance data files (*ABS.dat) found in directory, skipping absorbance analysis" ) )
    

  # set up output directory 'processed_data' and write out the complied sample information
    if(!dir.exists("processed_data")){
      dir.create("processed_data")
    }

    # write out sample sheet to a CSV file
    sample_sheet_file = paste0(run_name, "_sample_sheet.csv")
    writeLines(paste0("writing out sample_sheet as: ", sample_sheet_file))
    write.csv(sample_sheet, file.path("processed_data",sample_sheet_file), row.names = F)
  
  # 2. Read in EEM data --------------------------------------------------
    # Using the files organized above, read in the excitation/emission matrices for each sample
  
  # Aqualog parameters
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
  
  
  # calculate array size
    # these develop axis values for Ex/Em, example 500 495 490 ... 245 240
    Ex = seq(from = ExEnd, to = ExStart, by = -ExStep) 
    Em = seq(from = EmStart, to = EmEnd, by = EmStep)
  
    # this selects the appropriate range of the scans to be read in by dlmread
    # Typical Aqualog scans for Nelson Lab are 53 columns (Ex) and 125 rows (Em)
    Ex_range <- length(Ex)
    Em_range <- length(Em)
  
  
  # get SEM data from files as matrix
    writeLines("extracting EEM data from files")

    read_SEM <- function(sem_file){
      # read specified region of the emission/excitation matrix
      as.matrix(read.table(sem_file, header = T, row.names = 1)[1:Em_range, 1:Ex_range])
    }
    
    Sample_EEMs <- lapply(SEM_files, read_SEM) # sample EEM data
    Blank1_EEMs <- lapply(blank1, read_SEM) # blank1 EEM data
    Blank2_EEMs <- lapply(blank2, read_SEM) # blank2 EEM data
    
    # average blanks 1 and 2
    BlankAverages <- lapply(1:length(Blank1_EEMs), function(i) (Blank1_EEMs[[i]] + Blank2_EEMs[[i]]) / 2 ) 
  
  # 3. Read in absorbance data ------------------------------------------
    # If there was absorbance data recorded, read the ABS files
    # Note: absorbance measurements are used to correct excitation/emission values
    # The equation and calculations to perform the correction are detailed below
  
    if(isTRUE(ABS_collected)){
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
  
  # 4. Perform EEM corrections -------------------------------
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
    
    BlankRamanAreas <- lapply(BlankAverages, Raman_integration)
  
  # a) Apply Inner Filter Correction
    # multiply sample EEM values by the IFE correction values calculated from absorbances (see section 3. above)
    if(isTRUE(ABS_collected)){
    
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
                              function(i) Sample_EEMs_IFE[[i]] / BlankRamanAreas[[ which(group_names == groups[i]) ]])
    # also correct the blanks
    Blank1_EEMs_Ram <- lapply(1:N_blanks,
                         function(i) Blank1_EEMs[[i]] / BlankRamanAreas[[i]])
    Blank2_EEMs_Ram <- lapply(1:N_blanks,
                         function(i) Blank2_EEMs[[i]] / BlankRamanAreas[[i]])
  
  # c) Subtract Blanks (Raman Scaling of blanks by average)
    # divide the averaged blanks by their Raman area, then subtract that from the samples in each group
    BlankAverages_Ram <- lapply(1:N_blanks,
                                function(i) BlankAverages[[i]] / BlankRamanAreas[[i]] )
    # subtract the Raman scaled blanks from the Raman scaled samples in the appropriate group
    Sample_EEMs_BlkSub <- lapply(1:N_samples,
                                 function(i) Sample_EEMs_Ram[[i]] - BlankAverages_Ram[[ which(group_names == groups[i]) ]])
  
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
    Blank1_Ray <- lapply(Blank1_EEMs_Ram, Remove_Rayleigh)
    Blank2_Ray <- lapply(Blank2_EEMs_Ram, Remove_Rayleigh)
    
   
    
  # 5. Write out processed sample EEM matrices -----------------------------
    
    # clean names of excitation spectra (X500 -> 500)
    Clean_Colnames <- function(eem){
      colnames(eem) <- sub("X","",colnames(eem))
      return(eem)
    }
    
    Sample_EEMs_Final <- lapply(Sample_EEMs_Ray, Clean_Colnames)
    Blank1_Final <- lapply(Blank1_Ray, Clean_Colnames)
    Blank2_Final <- lapply(Blank2_Ray, Clean_Colnames)
    
    names(Sample_EEMs_Final) <- sample_names
    names(Blank1_Final) <- blank1
    names(Blank2_Final) <- blank2
    
  
    if(!dir.exists(file.path("processed_data","processed_matrices"))){
      dir.create(file.path("processed_data","processed_matrices"))
    }
    
    
    sample_output_files <- file.path("processed_data","processed_matrices",
                                     paste(run_name, names(Sample_EEMs_Final), "clean.csv", sep = "_"))
    
    writeLines("writing out cleaned EEM files to directory processed_data > processed_matrices")
    invisible( lapply(1:N_samples,
           		function(i)   write.csv(Sample_EEMs_Final[[i]], sample_output_files[i]) ))
    
  # 6. Calculate indices ----------------------------------
    # using the cleaned EEM matrices, calculate indices for various compounds/groups of compounds
    writeLines("calculating indices")
    # define a function to find the emission value in the data that is closest to the emission value specified by the index
    Get_Nearest_Em <- function(Em_val){
      Em_dif <- abs(Em - Em_val)
      which(Em_dif == min(Em_dif))
    }
    
  # List emission and excitation values of each index
  # simple indeces can be defined as a vector of c(emission value, excitation value)
  # more complicated indices relate multiple emission/excitation values to each other
  
    # simple indices involving a single peak
    simple_indices <- list(
                  CobleA = c(450,320),
                  CobleB = c(305,275),
                  CobleC = c(445, 345),
                  CobleM = c(410, 310),
                  CobleT = c(340,275),
                  
                  Fpeak = c(299,240),
                  Stedmon_D = c(509,390),
                  Optical_Brighteners = c(435,360),
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
                                        sapply(Sample_EEMs_Final,
                                               function(sample_eem)
                                                 sample_eem[nearest_em, exact_ex])
                                      })
    
  # Calculate complex indices for each sample
    # these indices involve ranges of emissions and/or ratios of emissions
    
    # BIX, divide 380 nm by 430 nm emission @ 310 nm excitation
    indices_out[["BIX"]] <- sapply(Sample_EEMs_Final,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(380, 430), Get_Nearest_Em)
                                     exact_ex <- which(Ex == 310)
                                     
                                     sample_eem[nearest_ems[[1]], exact_ex] / sample_eem[nearest_ems[[2]], exact_ex]
                                   })
    
    # HIX, divide the sum of 434 nm - 480 nm emissions by the sum of 300 nm - 346 nm emission @ 255 nm excitation
    indices_out[["HIX"]] <- sapply(Sample_EEMs_Final,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(434, 480, 300, 346), Get_Nearest_Em)
                                     exact_ex <- which(Ex == 255)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex]) /
                                       sum(sample_eem[nearest_ems[[3]]:nearest_ems[[4]], exact_ex])
                                   })
    
    # FI, divide 470 nm by 520   nm emission @ 370 nm excitation
    indices_out[["FI"]] <- sapply(Sample_EEMs_Final,
                                  function(sample_eem) {
                                    nearest_ems <- lapply(c(470, 520), Get_Nearest_Em)
                                    exact_ex <- which(Ex == 310)
                                    
                                    sample_eem[nearest_ems[[1]], exact_ex] / sample_eem[nearest_ems[[2]], exact_ex]
                                  })
    
    # M to C ratio
    indices_out[["M_to_C"]] <- indices_out[["simple"]][ , "CobleC"] / 
                                indices_out[["simple"]][ , "CobleM"]
    # JP5 Jet fuel, based on ordination of original community samples. This module identifies contaminated samples.
    indices_out[["JP5_empirical"]] <- sapply(Sample_EEMs_Final,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(308, 344), Get_Nearest_Em)
                                     exact_ex <- which(Ex %in% 260:290)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex])
                                   })
    # Diesel, fairly wide range (Chen et al 2021 0.2166/aqua.2021.120)
    indices_out[["Diesel_Chen"]] <- sapply(Sample_EEMs_Final,
                                   function(sample_eem) {
                                     nearest_ems <- lapply(c(280, 445), Get_Nearest_Em)
                                     exact_ex <- which(Ex %in% 220:275)
                                     
                                     sum(sample_eem[nearest_ems[[1]]:nearest_ems[[2]], exact_ex])
                                   })

    writeLines(paste("index calculated: ", c(names(simple_indices), names(indices_out)[-1]), collapse = "\n"))

    # write out calculated indices as a CSV file
    writeLines("writing out indices to processed_data > fDOM_indices_out.csv")
    indices_dat <- as.data.frame(do.call("cbind", indices_out))
    indices_dat$sample   <- row.names(indices_dat)
    indices_dat$run_name <- rep(run_name, nrow(indices_dat))
    indices_dat$UniqueID <- sample_sheet$UniqueID[match(indices_dat$sample, sample_sheet$sample)]
    
    write.csv(indices_dat, file.path("processed_data", paste0(run_name,"_fDOM_indices_out.csv")), row.names = F)
  
    writeLines("DONE")
  
    
