# Take sums of everything in the outputs from GMOlabeler, creating sum categories.

library(data.table)
library(optparse)

# Read arguments from command line  ---------------------------------------


option_list = list(
  make_option(c("-d", "--datapath"),
              type="character",
              default=NULL,
              help="data output from GMOlabeleR",
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#datapath <- "T16_DEV_genes/EA/wk7/"
datapath <- paste0("/scratch2/NSF_GWAS/GMOlabeler/output/", opt$datapath, "stats.csv")
output <-fread(datapath)

#### First do for all segments, combine all
all_segment_npixel_table <- aggregate(n_pixels_passing_threshold ~ filename + grid_item, data = output, FUN=sum)
all_segment_maxsignal_table <- aggregate(max_signal ~ filename + grid_item, data = output, FUN=max)
all_segment_totalsignal_table <- aggregate(total_signal ~ filename + grid_item, data = output, FUN=sum)
all_segment_totalpixel_table <- aggregate(total_pixels ~ filename + grid_item, data = output, FUN=sum)

all_segment_table_combined <- merge(merge(merge(all_segment_npixel_table,
                                                all_segment_maxsignal_table),
                                          all_segment_totalsignal_table),
                                    all_segment_totalpixel_table)


all_segment_table_combined$segment_hex <- rep("All_tissue", n=nrow(all_segment_table_combined))


#### Now only combine for callus and shoot (or more specifically, everything except stem)
output_sans_stem <- output[which(output$segment_hex!="CC0000"),]

### Do the same thing as before
stemless_npixel_table <- aggregate(n_pixels_passing_threshold ~ filename + grid_item, data = output_sans_stem, FUN=sum)
stemless_maxsignal_table <- aggregate(max_signal ~ filename + grid_item, data = output_sans_stem, FUN=max)
stemless_totalsignal_table <- aggregate(total_signal ~ filename + grid_item, data = output_sans_stem, FUN=sum)
stemless_totalpixel_table <- aggregate(total_pixels ~ filename + grid_item, data = output, FUN=sum)

stemless_table_combined <- merge(merge(merge(stemless_npixel_table,
                                                stemless_maxsignal_table),
                                          stemless_totalsignal_table),
                                    stemless_totalpixel_table)

stemless_table_combined$segment_hex <- rep("All_regenerated_tissue", n=nrow(stemless_table_combined))

##### Combine these two from above
sum_segment_tables_combined <- rbind(all_segment_table_combined,
                                     stemless_table_combined)

# Fill in other columns
sum_segment_tables_combined$V1 <- rep("observation", n=nrow(sum_segment_tables_combined))
sum_segment_tables_combined$mean_signal <- rep(NA, n=nrow(sum_segment_tables_combined))
sum_segment_tables_combined$intensity_threshold <- rep(NA, n=nrow(sum_segment_tables_combined))

# Put columns in same order so we can rbind combined segment stats with original data (separate segments)
setcolorder(sum_segment_tables_combined, neworder=colnames(output))

final_combined_table <- rbind(output, sum_segment_tables_combined)
new_datapath <- gsub("stats", "stats_with_sums_over_tissues", datapath)

print(paste0("Writing output with sums statistics calculated over combined tissue segments to: ",
             new_datapath))
fwrite(final_combined_table, new_datapath, row.names = FALSE, quote = FALSE)
