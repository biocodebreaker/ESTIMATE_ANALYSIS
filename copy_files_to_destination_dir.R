# set the directory where the TSV files are located
dir_path <- "/home/r1816512/.cache/GenomicDataCommons"


# set the destination directory where the copied files will be saved
destination_dir <- "/home/r1816512/TSV_subset"

# get a list of all TSV files in the source directory
filelist <- list.files(dir_path, pattern="*.tsv$", recursive=TRUE, full.names=TRUE)

# copy the first 10 TSV files to the destination directory
for (i in 1:10) {
  file.copy(from = filelist[i], to = destination_dir)
}