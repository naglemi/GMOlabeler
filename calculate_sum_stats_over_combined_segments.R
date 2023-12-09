library(data.table)
library(optparse)
`%notin%` <- Negate(`%in%`)

# Read arguments from command line  ---------------------------------------

option_list = list(
  make_option(c("-d", "--datapath"),
              type="character",
              default=NULL,
              help="data output from GMOlabeler",
              metavar="character"),
  make_option(c("-o", "--output_dir"),
              type="character",
              default=NULL,
              help="Output directory with results for GMOlabeler",
              metavar="character"),
  make_option(c("-k", "--keypath"),
              type="character",
              default=NULL,
              help="Path to the segmentation model key",
              metavar="character"),
  make_option(c("-e", "--exclude_tissues"),
              type="character",
              default="",
              help="String of tissue names to exclude, separated by spaces",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

arg_out_path <- paste0(opt$output_dir, opt$datapath, "combine_args.rds")

saveRDS(opt, file = arg_out_path)

### IF DEBUGGING IN RSTUDIO, UNCOMMENT THIS LINE INSTEAD OF USING OptParser
#opt <- readRDS("/home/gmobot/GMOlabeler/combine_args.rds")
#opt$output_dir <- "/home/gmobot/GMOGUI/output/gmodetector_out/gmolabeler_logic_outputs/output/"

datapath <- paste0(opt$output_dir, opt$datapath, "stats.csv")
output <- fread(datapath)

#### First do for all segments, combine all
all_segment_npixel_table <- aggregate(n_pixels_passing_threshold ~ filename + grid_item, data = output, FUN=sum)
all_segment_maxsignal_table <- aggregate(max_signal ~ filename + grid_item, data = output, FUN=max)
all_segment_totalsignal_table <- aggregate(total_signal ~ filename + grid_item, data = output, FUN=sum)
all_segment_totalpixel_table <- aggregate(total_pixels ~ filename + grid_item, data = output, FUN=sum)

all_segment_table_combined <- merge(merge(merge(all_segment_npixel_table,
                                                all_segment_maxsignal_table),
                                          all_segment_totalsignal_table),
                                    all_segment_totalpixel_table)

all_segment_table_combined$segment_hex <-
  rep("All_tissue", n=nrow(all_segment_table_combined))

print("Completed for whole-tissue batch")

# Load the segmentation model key if provided
if (!is.null(opt$keypath) && nzchar(opt$keypath)) {
  segmentation_key <- fread(opt$keypath)
  exclude_tissues_vector <- unlist(strsplit(opt$exclude_tissues, " "))
  hexes_to_exclude <- segmentation_key$hex_code[
    which(segmentation_key$Tissue %in% exclude_tissues_vector)]
} else {
  segmentation_key <- data.table(Tissue=character(), hex_code=character(), Integer=numeric())
  hexes_to_exclude <- c("CC0000")
}

# Subset the data.table using this logical mask
output_sans_excluded <- output[output$segment_hex %notin% hexes_to_exclude, ]

print("Now to compute only for regenerated tissues")
print("Rows before excluding: ")
print(nrow(output))

print("Rows after excluding: ")
print(nrow(output_sans_excluded))

### Do the same thing as before with the excluded tissues
excluded_npixel_table <- aggregate(n_pixels_passing_threshold ~ filename + grid_item, data = output_sans_excluded, FUN=sum)
excluded_maxsignal_table <- aggregate(max_signal ~ filename + grid_item, data = output_sans_excluded, FUN=max)
excluded_totalsignal_table <- aggregate(total_signal ~ filename + grid_item, data = output_sans_excluded, FUN=sum)
excluded_totalpixel_table <- aggregate(total_pixels ~ filename + grid_item, data = output_sans_excluded, FUN=sum)

excluded_table_combined <- merge(merge(merge(excluded_npixel_table,
                                             excluded_maxsignal_table),
                                       excluded_totalsignal_table),
                                 excluded_totalpixel_table)

excluded_table_combined$segment_hex <- rep("All_regenerated_tissue", n=nrow(excluded_table_combined))

##### Combine these two from above
sum_segment_tables_combined <- rbind(all_segment_table_combined,
                                     excluded_table_combined)

# Fill in other columns
sum_segment_tables_combined$V1 <- rep("observation", n=nrow(sum_segment_tables_combined))
sum_segment_tables_combined$mean_signal <- rep(NA, n=nrow(sum_segment_tables_combined))
sum_segment_tables_combined$intensity_threshold <- rep(NA, n=nrow(sum_segment_tables_combined))

# Put columns in same order so we can rbind combined segment stats with original data (separate segments)
setcolorder(sum_segment_tables_combined, neworder=colnames(output))

final_combined_table <- rbind(output, sum_segment_tables_combined)
new_datapath <- gsub("stats", "stats_with_sums_over_tissues", datapath)

print(paste0("Writing to ", new_datapath))
fwrite(final_combined_table, new_datapath)

print("Script completed successfully.")
