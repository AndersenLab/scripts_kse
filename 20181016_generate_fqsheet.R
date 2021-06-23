library(tidyverse)

# make sure you are in the same directory as the files, or provide a directory here
files = list.files()

# filter only the files you are interested in, in my case only "ECA" names that were over 1000
files2 = grep("ECA10", files, value = T)

# make a dummy dataframe
df <- data.frame(strain = NA, fastq_pair_id = NA, library = NA, fastq_1_path = NA, fastq_2_path = NA)

# get the unique sample IDs
samples <- unique(stringr::str_split_fixed(files2, "_", 2)[,1])

# make the fq_sheet
for(i in samples) {
    strain = i
    fastq1 = grep("R1", grep(strain, files2, value = T), value = T)
    fastq2 = grep("R2", grep(strain, files2, value = T), value = T)
    id = stringr::str_split_fixed(fastq1, "_R", 2)[,1]
    library = stringr::str_split_fixed(id, "_", 2)[,2]
    df <- rbind(df, c(strain, id, library, fastq1, fastq2))
}

df <- df %>% 
    na.omit()

readr::write_tsv(df, "fq_sheet.tsv", col_names = F)