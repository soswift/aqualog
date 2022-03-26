# This script compiles the outputs of the process_aqualog.R script across multiple runs
# processed EEM matrices are saved as an R object
# indices are compiled into a master table with all samples and all indices
  setwd("~/Documents/Bioinformatics/Projects/Redhill/processing/scripts/")
  
  out_dir  = "../../data/fdom_compiled/"
  data_dir = "../../data/fdom_runs"
  
  runs = c("11Dec2021",
           "20211217",
           "20220114")

  # master sample sheet
  sample_sheet_file = "../../data/metadata/RH_Nelson_Lab_Sample_Sheet - all_samples.csv"
  sample_sheet = read.csv(sample_sheet_file, stringsAsFactors = F, header = T)

  head(sample_sheet)
  

# Compile data -------------------------------------
  
  # read data from runs
  data_list = lapply(runs, read_processed_data)

  compiled_indices  = do.call("rbind", indices_list)

  # check
  dupd_sams = duplicated(compiled_indices$UniqueID)
  if(any(dupd_sams)){
    stop( paste("ERROR duplicated UniqueID:",
                compiled_indices$UniqueID[dupd_sams],
                compiled_indices$run_name)
    )
  }

  # samples and indices
  all_samples= merge(sample_sheet, compiled_indices, by = "UniqueID", all.y = T,)
  if(nrow(compiled_indices) != nrow(all_samples)) stop("ERROR sample indices merge error")

  # EEMs
  eem_run_list = lapply(data_list, function(x) x$eems)
  eem_list = unlist(eem_run_list, recursive = F)
  all_eems = do.call("rbind", eem_list)
  all_eems$UniqueID = gsub("\\..+","",row.names(all_eems))

# write out
  write.csv(all_samples, file.path(out_dir,"all_sample_indices.csv"), row.names = F)
  write.csv(all_eems, file.path(out_dir,"all_sample_processed_EEMs.csv"), row.names = F)
  saveRDS(eem_list, file.path(out_dir,"all_EEMs.rds"))
                                 