# Create a data frame for 4 samples with paired-end fastqs
sample_sheet <- data.frame(
  sample_id = c("statocyst_1", "statocyst_2", "gill_1", "gill_2"),
  tissue    = c("statocyst", "statocyst", "gill", "gill"),
  replicate = c(1, 2, 1, 2),
  fastq_1 = c("1_S1_R1.fastq.gz",
              "2_S5_R1.fastq.gz",
              "DP-Gill-1_S1_R1.fastq.gz",
              "DP-Gill-2_S2_R1.fastq.gz"),
  fastq_2 = c("1_S1_R2.fastq.gz",
              "2_S5_R2.fastq.gz",
              "DP-Gill-1_S1_R2.fastq.gz",
              "DP-Gill-2_S2_R2.fastq.gz"),
  condition = c("Statocyst", "Statocyst", "Gill", "Gill")
)

sample_sheet

# Save it as CSV
write.csv(sample_sheet, "data/metadata/sample_sheet.csv", row.names = FALSE)