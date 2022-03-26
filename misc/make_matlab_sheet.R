# file directory
path = "path"
run_name = "round_3"

# set working directory
setwd(path)

# aqualog outputs data files (*.dat) and blank files (*.ogw)
# *.dat files include raw data BEM, SEM, and ABS
writeLines(paste0("Reading data files in directory: ", path))

all_files <- list.files(path)
dat_files <- grep("\\.dat", all_files, value = T)
blank_files <- grep("Sample0001BEM", dat_files, value = T)
SEM_files <- grep("SEM", dat_files, value = T)
ABS_files <- grep("ABS", dat_files, value = T)

# summarize files and samples
samples      <- gsub("(Sample\\d+).+", "\\1", dat_files, perl = T)
sample_names <- unique(samples)

N_samples    <- length(sample_names)
N_blanks     <- length(blank_files)

if(N_samples < 1)  stop("No sample data files (*.dat) in the directory!") 
if(N_blanks  < 1)  stop("No blank files (*.ogw) in the directory!")

groups       <- gsub("(Group\\d+).+", "\\1", SEM_files, perl = T)
group_names  <- unique(groups) 
N_groups     <- length(group_names)

samples_per_group <- sapply(group_names, function(x) length(groups[groups == x]))

if(any(samples_per_group > 7)){ warning("There are a lot of samples per blank, maybe double check that"); print(samples_per_group)}
if(N_groups != N_blanks) warning("The number of groups does not match the number of blanks!")

typical_group_size <- names(sort(table(samples_per_group), decreasing = TRUE)[1])
last_group_size    <- tail(samples_per_group, 1)


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

# Quick script to make matlab sample sheet
sample_sheet <- data.frame(
  unique_id = paste0("ID_", 1:length(sample_names)),
  sample = sample_names,
  group  = groups,
  SEM_file = SEM_files,
  blank1 = blank1[groups],
  blank2 = blank2[groups]
)

# add absorbance files
sample_sheet$ABS_file <- ABS_files

# insert blanks between samples
Blank <-rep("Sample", N_samples + N_blanks)
Blank[ seq(1, 
           N_samples + N_blanks - last_group_size,
           by = as.numeric(typical_group_size)+1)] <- "Blank"

matlab_blanks <- data.frame( UniqueID = seq_along(Blank),
                             Blank = Blank)


# matlab sample sheet
matlab_sheet <- data.frame( UniqueID = matlab_blanks[matlab_blanks$Blank == "Sample", "UniqueID"],
                            Samplepath = file.path(path,sample_sheet$SEM_file),
                            Blank1 = file.path(path,sample_sheet$blank1),
                            Blank2 = file.path(path, sample_sheet$blank2),
                            abs = file.path(path, sample_sheet$ABS_file),
                            stringsAsFactors = F
)

# merge sample sheet with blanks
matlab_sheet <- merge(matlab_blanks, matlab_sheet,  by = "UniqueID", all.x = T)

# add file paths for Blanks
blank_paths <- unique(matlab_sheet$Blank1)[-1]

matlab_sheet[matlab_sheet$Blank == "Blank", "Samplepath"] <- blank_paths
matlab_sheet[matlab_sheet$Blank == "Blank", "Blank1"] <- c(tail(blank_paths, -1), head(blank_paths, 1))
matlab_sheet[matlab_sheet$Blank == "Blank", "Blank2"] <- c(tail(blank_paths, -2), head(blank_paths, 2))
matlab_sheet[matlab_sheet$Blank == "Blank", "abs"] <- grep("Sample0001",matlab_sheet$abs, value = T)

# check for whitespace
if(any(grepl(" ", matlab_sheet))) warning("!!! There's whitespace in your sample sheet, will cause error !!!")

# write out
write.table(matlab_sheet, paste0(run_name,"_matlab_sheet.tsv"), sep = "\t", row.names = F)

            