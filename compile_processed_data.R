# This script compiles the outputs of the process_aqualog.R script across multiple runs
# processed EEM matrices are saved as an R object
# indices are compiled into a master table with all samples and all indices
  setwd("~/Documents/Bioinformatics/Projects/Redhill/aqualog_script/script/")
  
  # output directory for compiled data
  out_dir <- "../../data/compiled/"

  # identify master sample sheet
  sample_sheet_file <- "../../data/RH_Nelson_Lab_Sample_Sheet - all_samples.csv"
  sample_sheet <- read.csv(sample_sheet_file, stringsAsFactors = F, header = T)

  # identify run directories in the 'data' directory
  data_dir <- "../../data"
  
  runs <- c("11Dec2021",
            "20211217")
  
  # read function -----------------------------------------
  # This function reads in the all the relevant data for a run
  # As it reads in, it checks that the sample numbers all match up
  
  read_processed_data <- function(run){
    run_path <- file.path(data_dir, run, "processed_data")
  
    # sample sheet
    run_file <- file.path(run_path,
                          paste0(run, "_sample_sheet.csv"))
     if (!file.exists(run_file)) stop(
       paste("ERROR Can't find sample sheet", run)
       )
    run_sheet <- read.csv(run_file, header = T, stringsAsFactors = F)

    # indices
    indices_file <- file.path(run_path,
                              paste0(run, "_fDOM_indices_out.csv"))
     if(!file.exists(indices_file)) stop(
       paste("ERROR Can't find calculated indices for run", run)
       )
    
    run_indices <- read.csv(indices_file, header = T, stringsAsFactors = F)
    if (!any(run_indices$UniqueID %in% sample_sheet$nelson_lab_id)) stop(
      paste( "ERROR Samples in indices don't match run sample sheet", run)
      )

    # processed EEM matrices
    eem_dir   <- file.path(run_path, "processed_matrices")
    eem_files <-  list.files(eem_dir)
    
    run_eems <- lapply(eem_files,
                       function(x)
                         read.csv(
                           file.path(eem_dir, x),
                           header = T,
                           stringsAsFactors = F
                         ))
  
    eem_samples <- gsub(".+(Group.+)_.+", "\\1", eem_files)
    eem_uniqID <- run_sheet$UniqueID[match(eem_samples, run_sheet$sample)]
    
    names(run_eems) <- eem_uniqID
    
    # check that numbers match
    n_samples = list(sheet = nrow(run_sheet),
                  indices = nrow(run_indices),
                  eems = length(run_eems))
    
    writeLines(paste( "Samples in", 
                 "sheet = ", n_samples$sheet,
                 "indices = ", n_samples$indices,
                 "eems = ", n_samples$eems,
                 "for run", run))
    
    if(length(unique(n_samples)) != 1) stop(
      paste("Sample N mismatch in run ",run)
    )
    
    #output sample sheet, indices, eems
    out <- list(samples = run_sheet,
                indices = run_indices,
                eems = run_eems)
    return(out)
  }
  

# Compile data -------------------------------------
  
# compile output
  data_list <- lapply(runs, read_processed_data)

  # indices
  indices_list <- lapply(data_list, function(x) x$indices)
  all_indices <- do.call("rbind", indices_list)

  dupd_sams <- duplicated(all_indices$UniqueID)
  if(any(dupd_sams)) stop(
    paste("ERROR duplicated UniqueID:", 
          all_indices$UniqueID[dupd_sams],
          all_indices$run_name)
  )

  all_samples<- merge(sample_sheet, all_indices, by = "UniqueID", all.y = T,)
  if(nrow(all_indices) != nrow(all_samples)) stop("ERROR sample indices merge error")

  # EEMs
  eem_run_list <- lapply(data_list, function(x) x$eems)
  eem_list <- unlist(eem_run_list, recursive = F)
  all_eems <- do.call("rbind", eem_list)
  all_eems$UniqueID <- gsub("\\..+","",row.names(all_eems))

# write out
  write.csv(all_samples, file.path(out_dir,"all_sample_indices.csv"), row.names = F)
  write.csv(all_eems, file.path(out_dir,"all_sample_processed_EEMs.csv"), row.names = F)
  saveRDS(eem_list, file.path(out_dir,"all_EEMs.rds"))
                                 