

### Testing -------------------------------------------------
test_dir = "~/Documents/Bioinformatics/Projects/Redhill/processing/scripts/testing/"

setwd(test_dir)
source("../process_aqualog_functions.R")
source("../eem_plot_functions.R")

data_dir = "~/Documents/Bioinformatics/Projects/Redhill/data/fdom_runs/20211217/"
run = "test"
key = "RH_round2_run_record.tsv"

test_out = file.path(test_dir, "test_out")

# processing functions -----------------
# organize data
org_d = organize_aqualog(
                           data_directory = data_dir,
                           run_name = run,
                           sample_key_file = key,
                           out_dir = test_out
                         )
# process EEM data
cor_em = correct_EEMs(
                       data_directory = data_dir,
                       sample_sheet = org_d,
                       run_name = run,
                       out_dir = test_out
                       )

# calculate indices
inds = fdom_indices(
                      Sample_EEMs = cor_em,
                      sample_sheet = org_d,
                      run_name = run,
                      out_dir = test_out
                    )

# main function (wraps the three previous functions)
main = process_aqualog(
                        data_directory = data_dir,
                        run_name = run,
                        sample_key_file = key,
                        out_dir = test_out
                      )


# Compile functions --------------
runs = c("test_out")

compiled_out_dir  = file.path(test_out, "compiled")
procd_data_dir = test_dir


all_data = compile_runs(run_dirs = runs,
                        data_dir = procd_data_dir,
                        out_dir = compiled_out_dir)



# plotting functions --------------------
indx = main$indices
eems = readRDS("test_out/processed_data/test_EEM.rds")
samples = main$sample_sheet

sample_sheet_file = "../../../data/metadata/RH_Nelson_Lab_Sample_Sheet - all_samples.csv"
std_curve_file = "~/Documents/Bioinformatics/Projects/Redhill/data/fdom_compiled/dilution_indices.csv"

sample_sheet = fread(sample_sheet_file)
std_curve = fread(std_curve_file)

# TODO: troubleshoot eems list plotting huge y axis


plots = lapply(1:length(eems), function(i) {
  sample_id = samples$UniqueID[[i]]
  collection_code = sample_sheet[UniqueID == sample_id, collection_code]
  
  jp5_value = indx$JP5_normalized[indx$UniqueID == sample_id]
  
  
  # eem plot
  eem_p = plot_eem(eems[[i]], sample_name = sample_id)
  eem_p = jp5_annotation(eem_p, round(jp5_value, 3))
  
  # severity plot
  sev_p = jp5_severity_plot(jp5_value = jp5_value)
  
  # arrange
  p = plot_grid(eem_p, sev_p, rel_widths = c(2, 1))
  
  return(p)
})

names(plots) = names(eems)

pdf(file.path("test_out","plots.pdf"), width = 10)
plots
dev.off()
