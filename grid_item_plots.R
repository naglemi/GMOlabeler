# Warning: This script is not my most elegantly-written work.
#  There remain a plethora of debugging statements, which are being left here
#  in case any bugs resulting from corner cases appear again. Ideally, these
#  should be replaced with error codes that catch and describe the errors
#  and solutions for the user to fix them.
#  Also, plots have not been functionalized and there is a lot of redundant
#  code for making plots. I'm keeping it this way for now, as these plots
#  are tentative and we may wish to add more soon.
#  To navigate this code, it is suggested to use the RStudio document
#  outline panel (Ctrl+Shift+O)

library(data.table)
library(scales)
library(readxl)
library(ggplot2)
library(randomcoloR)
library(optparse)
library(tools)
library(stringr)
library(permutes)
library(buildmer)
library(tidyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

# Read arguments from command line  ---------------------------------------

option_list = list(
  make_option(c("-d", "--datapath1"),
              type="character",
              default=NULL,
              help="data output from GMOlabeleR",
              metavar="character"),
  make_option(c("-r", "--randomization_datasheet_path"),
              type="character",
              default=NULL,
              help="path to randomization datasheet",
              metavar="character"),
  make_option(c("-S", "--samples-pre-labeling"),
              type="character",
              default=NULL,
              help=paste("path to output from pre_label.R, used",
              "again here to filter out from missing",
              "explant data plates that were imaged twice or more",
              "by removing all images except the last from analysis."),
              metavar="character"),
  make_option(c("-p", "--pixel_threshold"),
              type="numeric",
              default=5,
              help=paste("Number of pixels passing intensity threshold for",
                         "binary classifications of transgenic or not"),
              metavar="numeric"),
  make_option(c("-v", "--variable"),
              type="character",
              default="categorical",
              help="categorical or continuous",
              metavar="numeric"),
  make_option(c("-m", "--missing"),
              type="numeric",
              default="FALSE",
              help=paste("TRUE (1) if some explants are missing.",
                         "Must also specify list with -M"),
              metavar="numeric"),
  make_option(c("-M", "--MissingList"),
              type="character",
              default=FALSE,
              help="A file output from automated detection of missing explants",
              metavar="numeric"),
  make_option(c("-g", "--grid_type"),
              type="numeric",
              default="12",
              help="Grid type (currently supporting 12 and 20)",
              metavar="numeric"),
  make_option(c("-s", "--sort"),
              type="numeric",
              default="1",
              help=paste("Whether to sort genotypes according to difference in",
                         "effects of treatment on trangenic callus -",
                         "numeric 0 (False) or 1 (True)"),
              metavar="numeric"),
  make_option(c("-H", "--height"),
              type="numeric",
              default="7",
              help="Plot height (inches)",
              metavar="numeric"),
  make_option(c("-w", "--width"),
              type="numeric",
              default=-9,
              help="Plot width (inches)",
              metavar="numeric"),
  make_option(c("-e", "--debug"),
              type="numeric",
              default=1,
              help="Print debugging output",
              metavar="numeric"),
  make_option(c("-R", "--Reporter"),
              type="character",
              default="",
              help=paste("Specify reporter so when",
                         "running with multiple reporters, results",
                         "for each are collected and saved separately."),
              metavar="character"),
  make_option(c("-k", "--keypath"),
              type="character",
              default=FALSE,
              help="Path to the segmentation model key",
              metavar="character"),
  make_option(c("-u", "--outdir"),
              type="character",
              default=getwd(),
              help="Directory to save outputs in (default is the current working directory)",
              metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

# opt <- readRDS("/home/gmobot/GMOGUI/output/gmodetector_out//gmolabeler_stats_plots/ETFB/wk26///DsRed/gridplot_args.rds")
# opt$outdir <- "/home/gmobot/GMOGUI/output/gmodetector_out/"
# opt$randomization_datasheet_path <- gsub("/mnt/Elements_24/", "/media/gmobot/Elements_24/Elements_24/", opt$randomization_datasheet_path)
# opt$`samples-pre-labeling` <- gsub("/mnt/Elements_24/", "/media/gmobot/Elements_24/Elements_24/", opt$`samples-pre-labeling`)
# opt$keypath <- "/home/gmobot/GMOGUI/euc_cubeschool_a1-v4_cpua2.key.csv"

#setwd("/home/michael/GMOlabeler/")
#opt$debug <- 1

# Import and preprocess GMOlabeler data ----------------------------------------------

wd <- paste0(opt$outdir, "/gmolabeler_stats_plots/", opt$datapath1, "/", opt$Reporter, "/")
if(!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

arg_out_path <- paste0(wd, "gridplot_args.rds")

if(opt$debug == TRUE){
  print(paste0("Saving list of input arguments to : ", arg_out_path))
}

saveRDS(opt, file = arg_out_path)

datapath <- paste0(opt$outdir, "/gmolabeler_logic_outputs/", opt$datapath1, opt$Reporter,
                   "/stats_with_sums_over_tissues.csv")

if(opt$debug == TRUE){
  print(paste0("Reading in output from GMOlabeler at path: ", datapath))
}

if(opt$keypath != FALSE){
  pixel_demographics <- fread(opt$keypath)
  pixel_demographics$V1 <- NULL
  pixel_demographics$Integer <- NULL
  morep_demographics = data.frame(cbind(c('All_tissue',
                                          'All_regenerated_tissue'),
                                        c('None',
                                          'None')))
  colnames(pixel_demographics) <- colnames(morep_demographics) <-
    c('Tissue', 'hex_code')
  pixel_demographics <- rbind(pixel_demographics, morep_demographics)

} else {
    print("Assuming default tissues for original poplar model because not specified with `keypath`")
    pixel_demographics = data.frame(cbind(c('Shoot', 'Callus', 'Stem', 'All_tissue',
                                            'All_regenerated_tissue', 'Background'),
                                          c('00CC11', '0006CC', 'CC0000',
                                            'All_tissue',
                                            'All_regenerated_tissue', '000000'),
                                          c('green', 'blue', 'red',
                                            NA, NA, 'black')))
    colnames(pixel_demographics) <- c('Tissue', 'hex_code', 'color')
}

output <- fread(datapath)
output <- output[which(!is.na(grid_item)), ]

if(opt$debug == TRUE){
  cat("\n")
  print(paste0("Rows in output from GMOlabeler: ", nrow(output)))

  cat("\n")
  print(paste0("Max n_pixels_passing_threshold in output from GMOlabeler: ",
               max(na.omit(output$n_pixels_passing_threshold))))

  cat("\n")
  print(paste0("Max total_signal in output from GMOlabeler: ",
               max(na.omit(output$total_signal))))
}

output$filename <- gsub("_segment_uncropped_processed", "", output$filename)

# Saw this naming error in T19 (e.g. TAP1_I5.0_F1.9__L100_105739_1_0_1)
output$filename <- gsub("__", "_", output$filename)

# Merge in missing explant data -------------------------------------------

if (opt$MissingList %in% c(FALSE, 0, NA, "None", "none")) {
  # Create column names
  col_names <- c("image_name", as.character(1:opt$grid_type))
  # Create a row of data
  filler_data <- c("dummy", rep("NC", opt$grid_type))
  # Create dataframe
  missing_explant_data <-
    data.frame(matrix(filler_data, nrow = 1, ncol = length(filler_data)))
  names(missing_explant_data) <- col_names
} else {
  missing_explant_datapath <- opt$MissingList

  # character class should be specified explicitly to stop columns from
  #  being read as numeric when missing data comes from score data
  missing_explant_data <- fread(missing_explant_datapath,
                                header=TRUE,
				colClasses=c("character"))
}
#print("preparing to pivot missing explant data table")

# Need to exclude missing explant data for non-final imags of plates that
#  were re-imaged.

samples_pre_label <- fread(opt$`samples-pre-labeling`)
nonredundant_rgb_list <- basename(
  gsub("_processed", "", samples_pre_label$rgb))

missing_explant_data <- missing_explant_data[which(
  missing_explant_data$image_name %in% nonredundant_rgb_list), ]

missing_explant_data_tidy <- pivot_longer(
  data= missing_explant_data,
  cols = colnames(missing_explant_data)[-1],
  names_to = "grid_item",
  values_to = "present",
  values_drop_na = TRUE
)

colnames(missing_explant_data_tidy)[1] <- "filename"

# Saw this naming error in T19 (e.g. TAP1_I5.0_F1.9__L100_105739_1_0_1)
missing_explant_data_tidy$filename <- gsub("__", "_",
                                           missing_explant_data_tidy$filename)

# Get rid of file extension to be consistent with outputs
#  from automated missing explant script
output$filename <- gsub(".jpg", "", output$filename)

if(opt$debug == TRUE){
  cat("\n")
  print("Look at the top of output from GMOlabeler")
  print(head(output$filename))

  cat("\n")
  print("Look at the top of missing explant data")
  print(head(missing_explant_data_tidy))
}

output$filename <- basename(output$filename)
output$grid_item <- as.character(output$grid_item)

if(opt$debug == TRUE){
  cat("\n")
  print("Look at the head of filenames for missing explant data")
  print(head(missing_explant_data_tidy$filename))

  cat("\n")
  print("Look at the colnames of output from GMOlabeler")
  print(colnames(output))

  cat("\n")
  print("Look at the colnames of missing explant data")
  print(colnames(missing_explant_data_tidy))

  # We merge like this since we only import missing explant data for explant
  # that are missing. The rest is to be imputed as not being missing.

  cat("\n")
  print(paste0(
    "Before merging the GMOlabeler and missing output data,",
    "they have rows of ", nrow(output), " and ",
    nrow(missing_explant_data_tidy), " respectively."))
}

output <- merge(x = output,
                y = missing_explant_data_tidy,
                by = c("filename", "grid_item"),
                all.x = TRUE,
                all.y = FALSE)

# Remove explants labeled as missing or contaminated

if(opt$debug == TRUE){
  cat("\n")
  print(paste0("After merging but before removing missing/contminated",
               " explants, merged output has ",
               nrow(output), " rows."))
}

output$present[is.na(output$present)] <- "Y"

# cat("\n")
# print(paste0("All of these ",
#              nrow(output[which(output$present=="N"),]),
#              " rows are for missing explant"))
# print(output[which(output$present=="N"),])
#
# cat("\n")
# print(paste0("All of these ",
#              nrow(output[which(output$present=="C"),]),
#              " rows are for contaminated explant"))
# output[which(output$present=="C"),]

missing_explant_data[missing_explant_data=='NC'] <- 1
missing_explant_data[missing_explant_data=='Y'] <- 1
missing_explant_data[missing_explant_data=='N'] <- 0
missing_explant_data[missing_explant_data=='C'] <- 0
missing_explant_data[missing_explant_data=='P'] <- 0
missing_explant_data[missing_explant_data=='M'] <- 0

if(opt$debug == TRUE){
  cat("\n")
  print(paste0(
    "After removing the missing/contaminated explants, merged output has ",
    nrow(output)))
}

if(opt$debug == TRUE){
  cat("\n")
  print(paste0("After processing GMOlabeler output to merge with missing explant",
               "data, Head of data (w/ first 5 col): "))
  output[1:5,1:5]
}

#' Given macroPhor Array output filename, parse out tray and plate IDs
#'
#' @param filename for macroPhor Array output, with standardized naming
#'
#' @return A character string with tray ID and plate ID delimited by "_"
#' @export
#'
#' @examples
parse_trayplateID <- function(name_being_parsed){
  pass_to_dodge_error <- name_being_parsed
  imgpath_stripped <- file_path_sans_ext(basename(pass_to_dodge_error))
  trayID <- str_split_fixed(imgpath_stripped, "_", 2)[1]

  assign_ID_index_from_row_column_on_tray <- function(data_to_parse = filename,
                                                      components_list,
                                                      mode="table",
                                                      verbose=FALSE){
    dictionary <- cbind(c(0,0,0,0,0,0,0,
                          1,1,1,1,1,1,1,
                          2,2,2,2,2,2,2),
                        c(0,1,2,3,4,5,6,
                          0,1,2,3,4,5,6,
                          0,1,2,3,4,5,6),
                        c(1:21))
    dictionary <- as.data.table(dictionary)
    colnames(dictionary) <- c("row", "column", "ID")
    # Get the ID of position in tray in according to row and column
    dictionary$row_column <- paste0(dictionary$row, "_", dictionary$column)
    dictionary[,1:2] <- NULL
    if(mode=="table"){
      # Set colnames for spectral components if multiple are same
      #colnames(data_to_parse)[1:length(components_list)] <- components_list
      data_merged <- merge(data_to_parse, dictionary, by="row_column",
                           all.x = TRUE, all.y = TRUE)
      return(data_merged)
    }
    if(mode=="filename"){

      # Patch added in v0.19 for compatibility regardless of whether
      #  "_cyan" is at end of filename
      if(grepl("cyan", data_to_parse)==1){
        ndelimiters=9
      }else{
        ndelimiters=8
      }

      row <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)),
                             "_", 9)[ndelimiters-1]
      # Changed in v0.19 along with patch above
      col <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)),
                             "_", ndelimiters)[ndelimiters]
      row_col <- paste0(row, "_", col)
      ID <- dictionary[which(dictionary$row_column == row_col),]$ID
      # Debugging lines added in v0.19
      if(verbose==TRUE){
        print(paste0("This row is ", row))
        print(paste0("This col is ", col))
        print(paste0("This row_col is ", row_col))
        print(paste0("This filename (stripped) is ",
                     basename(file_path_sans_ext(data_to_parse))))
        print(paste0("This ID about to be returned from ",
                     "assign_ID_index_from_roW_column_on_tray is ", ID))
      }

      return(ID)
    }
  }

  plateID <- assign_ID_index_from_row_column_on_tray(
    data_to_parse = imgpath_stripped, mode="filename")

  trayplateID <- paste0(trayID, "_", plateID)
  return(trayplateID)
}

for(i in 1:nrow(output)){

  output$ID[i] <- parse_trayplateID(name_being_parsed = output$filename[i])
}


# Read randomization datasheet, clean if needed ---------------------------

if(grepl("\\.csv", opt$randomization_datasheet_path)){
  randomization_datasheet <- fread(opt$randomization_datasheet_path)
}
if(grepl("\\.xls", opt$randomization_datasheet_path)){
  randomization_datasheet <- as.data.table(
    read_excel(opt$randomization_datasheet_path))
}

colnames(randomization_datasheet) <-
  gsub("Tissue Type", "Tissue type", colnames(randomization_datasheet))

# If there are leaf explants, note that we can only
#  analyze stem explants in the current version
if (sum(grepl('Tissue type', colnames(randomization_datasheet))) >= 1){
  message(paste("Warning! This dataset contains multiple explant types. ",
                "We will only analyze stem/petiole explants."))

  if(opt$debug == TRUE){
    print("Randomization datasheet contains levels of tissue type: ")
    print(levels(factor(randomization_datasheet$`Tissue type`)))
    print(paste0("Nrow before subsetting to stem only: ",
                 nrow(randomization_datasheet)))
  }

  randomization_datasheet <- randomization_datasheet[which(
    randomization_datasheet$`Tissue type` == "S" |
    randomization_datasheet$`Tissue type` == "P" |
    randomization_datasheet$`Tissue type` == "S/P"), ]

  if(opt$debug == TRUE){
    print(paste0("Nrow after subsetting to stem/petiole only: ",
                 nrow(randomization_datasheet)))
  }
}

colnames(randomization_datasheet) <- gsub("total_explants",
                                          "total_explants_from_master_data",
                                          colnames(randomization_datasheet))

## Let's deal with inconsistently named columns here
colnames(randomization_datasheet) <- gsub("Genotype ID", "Genotype_ID",
                                          colnames(randomization_datasheet))
# If we don't already have a column named Treatment name,
#  convert the Treatment column to this.
if(sum(grepl("Treatment name", colnames(randomization_datasheet))) == 0){
  colnames(randomization_datasheet) <- gsub("Treatment", "Treatment name",
                                            colnames(randomization_datasheet))
}

if(opt$debug == TRUE){
  print(paste0("Upon loading randomization_datasheet, it has how many rows? ",
               nrow(randomization_datasheet)))
}

colnames(randomization_datasheet)[1:2] <- c("Image #", "Tray ID")
randomization_datasheet$ID <- paste0(randomization_datasheet$`Tray ID`,
                             "_",
                             randomization_datasheet$`Image #`)

output$transgenic <- rep(0, nrow(output))

output$transgenic[which(
  output$n_pixels_passing_threshold > opt$pixel_threshold)] <- 1

output$segment_present <- 0
output$segment_present[!is.na(output$total_signal)] <- 1

for (i in 1:nrow(pixel_demographics)){
    output$segment_hex <- gsub(pixel_demographics$hex_code[i],
                                      pixel_demographics$Tissue[i],
                                      output$segment_hex)
}

## Calculate total number of missing explants per plate
if (opt$MissingList %in% c(FALSE, 0, NA, "None", "none")) {
  randomization_datasheet$total_explants <- opt$grid_type
} else {

  missing_explant_data[missing_explant_data=='Y'] <- 1
  missing_explant_data[missing_explant_data=='N'] <- 0
  missing_explant_data[missing_explant_data=='C'] <- 0
  missing_explant_data[missing_explant_data=='P'] <- 0

  if(opt$debug == TRUE){
    cat("\n")
    print("Calculating total numbers of missing explants per plate")
  }

  missing_explant_data_with_totals <- missing_explant_data

  missing_explant_data_matrix <- as.matrix(
    missing_explant_data[,2:ncol(missing_explant_data)])

  class(missing_explant_data_matrix) <- "numeric"
  missing_explant_data_matrix[is.na(missing_explant_data_matrix)] <- 0

  missing_explant_data_with_totals$total_explants <-
    rowSums(missing_explant_data_matrix)

  if(opt$debug == TRUE){
    print("Dimensions of missing explant data with totals")
    print(dim(missing_explant_data_with_totals))
    print(paste("Minimum number of explants on any plate after",
                " removing missing and contaminated explants is..."))
    print(min(na.omit(missing_explant_data_with_totals$total_explants)))

    print(head(missing_explant_data_with_totals))
  }

  # Format to get ready to merge with randomization datasheet
  #  (which contains totals values per plate)
  colnames(missing_explant_data_with_totals)[1] <- "filename"

  missing_explant_data_with_totals <- cbind(
    missing_explant_data_with_totals,
    rep(NA, nrow(missing_explant_data_with_totals)))

  colnames(missing_explant_data_with_totals)[length(
    missing_explant_data_with_totals)] <- "ID"

  missing_explant_data_with_totals$filename <- gsub(
    "_rgb", "",missing_explant_data_with_totals$filename)

  missing_explant_data_with_totals$filename <- gsub(
    ".jpg", "", missing_explant_data_with_totals$filename)

  if(opt$debug == TRUE){
    cat("\n")
    print(paste("About to add IDs to missing explant data",
                "based on filenames. Here is head of filenames."))

    print(head(missing_explant_data_with_totals$filename))
  }

  cat("\n")
  if(nrow(missing_explant_data_with_totals) > 0){ # wont be case if "None" option
    for(i in 1:nrow(missing_explant_data_with_totals)){
      missing_explant_data_with_totals$ID[i] <- parse_trayplateID(
        name_being_parsed = missing_explant_data_with_totals$filename[i])
    }
  }

  if(opt$debug == TRUE){
    print(paste("We have added IDs to missing explant",
                "data based on filenames. Here is head of IDs."))
    print(head(missing_explant_data_with_totals$ID))
    cat("\n")

    print(paste("N rows of randomization datasheet",
                "before merging with missing explant data: "))
    print(nrow(randomization_datasheet))
    print("Head of IDs in missing explant data")
    print(head(missing_explant_data_with_totals$ID))
    print("Head of IDs in randomization datasheet")
    print(head(randomization_datasheet$ID))

    cat("\n")
    print(paste("Now about to merge randomization_datasheet",
                "and missing explant data."))
    print(paste0("Rows in randomization_datasheet: ",
                 nrow(randomization_datasheet)))
    print(paste0("Rows in missing explant data: ",
                 nrow(missing_explant_data_with_totals)))

    cat("\n")
    print(paste("Now about to merge randomization datasheet ",
                "and missing explant data, by ID. Look at head of IDs for both."))
    print(head(randomization_datasheet$ID))
    print(head(missing_explant_data_with_totals$ID))
  }

  randomization_datasheet <- merge(x = randomization_datasheet,
                                   y = missing_explant_data_with_totals,
                                   by = "ID",
                                   all.x = TRUE,
                                   all.y = FALSE)

  if(opt$debug == TRUE){
    print(head(randomization_datasheet))
    print(randomization_datasheet$total_explants)

    cat("\n")
    print("Dimensions of missing explant data with totals")
    print(dim(missing_explant_data_with_totals))
    print(paste("Levels of total explants on any plate",
                "after removing missing and contaminated explants is..."))
    print(levels(factor(randomization_datasheet$total_explants)))
    print(paste("Minimum number of explants on any plate",
                "after removing missing and contaminated explants is..."))
    print(min(na.omit(randomization_datasheet$total_explants)))
    cat("\n")

    print(paste0("Rows in merged output: ", nrow(randomization_datasheet)))
    print(randomization_datasheet$total_explants)

    cat("\n")

  }
}



#randomization_datasheet$total_explants[is.na(
#randomization_datasheet$total_explants)] <- opt$grid_type

randomization_datasheet$total_explants[is.na(
  randomization_datasheet$total_explants)] <- 0

if(opt$debug == TRUE){

  print(paste("N rows of randomization datasheet",
              "after merging with missing explant data: "))
  print(nrow(randomization_datasheet))

  print(paste("Column names of randomization datasheet",
              "after merging with missing explant data: "))
  print(colnames(randomization_datasheet))

  cat("\n")
  print(paste0("In output, Maximum observations for a segment",
               "in a grid item in a filename (should be 1): ",
               max(table(output$filename, output$segment_hex, output$grid_item))))
  print(paste0("In randomization datasheet,",
               "Maximum observations for an ID (should be 1): ",
               max(table(randomization_datasheet$ID))))

  cat("\n")
  print(paste("Before calculating totals, let's look at",
              "the column total_explants in randomization datasheet"))
  print(head(randomization_datasheet$total_explants))
  print(levels(factor(randomization_datasheet$total_explants)))
  cat("\n")
}

# Remove missing grid positions from output
## This specific step is needed for the situation where a grid item is
## removed from total explants, but there's still something there
## (due to overlap) so total % phenotype can be
## incorrectly counted as 1.

if (opt$MissingList %!in% c(FALSE, 0, NA, "None", "none")) {

  missing_explant_data_tidy$filename <- gsub("_rgb", "",
                                             missing_explant_data_tidy$filename)

  missing_explant_data_tidy$ID <- rep(NA, nrow(missing_explant_data_tidy))

  for(i in 1:nrow(missing_explant_data_tidy)){
    missing_explant_data_tidy$ID[i] <- parse_trayplateID(
      name_being_parsed = missing_explant_data_tidy$filename[i])
  }

  IDs_to_drop <- missing_explant_data_tidy[which(
    missing_explant_data_tidy$present!="NC"),]

  IDs_to_drop$ID_exp <- paste0(IDs_to_drop$ID, "_", IDs_to_drop$grid_item)
  output$ID_exp <- paste0(output$ID, "_", output$grid_item)

  if(opt$debug == TRUE){
    cat("\n")
    print(paste0("Nrow of output before dropping missing explants is: ",
                 nrow(output)))
    print(paste("Head of ID_exp columns from both output",
                "data table and IDs_to_drop data tbales, respectively:"))
    print(head(output$ID_exp))
    print(head(IDs_to_drop$ID_exp))
  }

  output <- output[output$ID_exp %!in% IDs_to_drop$ID_exp,]
  if(opt$debug == TRUE){
    print(paste0("Nrow of output AFTER dropping missing explants is: ",
                 nrow(output)))
    print("Head of output")
    print(head(output, 1))
    cat("\n")

    print(head(IDs_to_drop))
    print(colnames(output))
  }

  if(sum(is.na(randomization_datasheet$filename)) ==
     length(randomization_datasheet$filename)){
    stop(paste("Error: NA for all file names. Do the tray IDs on",
               "the randomization datasheet match the filenames of images?"))
  }

}


# Calculate transgenic ----------------------------------------------------

## First calculate all transgenic
# Initialize new columns
for (j in 1:nrow(pixel_demographics)){
    new_column_name <- paste0('n_transgenic_',
                              pixel_demographics$Tissue[j])

    randomization_datasheet <- cbind(
      randomization_datasheet,
      as.numeric(rep(NA,nrow(randomization_datasheet))))

    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <-
      new_column_name

    new_column_name <- paste0('portion_transgenic_',
                              pixel_demographics$Tissue[j])

    randomization_datasheet <- cbind(
      randomization_datasheet,
      as.numeric(rep(NA,nrow(randomization_datasheet))))

    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <-
      new_column_name
}

for (i in 1:nrow(randomization_datasheet)){
    for (j in 1:nrow(pixel_demographics)){
        column_of_interest <- paste0('n_transgenic_',
                                     pixel_demographics$Tissue[j])

        total_transgenic <- sum(na.omit(output$transgenic[which(
          output$ID==randomization_datasheet$ID[i] &
            output$segment_hex == pixel_demographics$Tissue[j])]))

        # print(paste0('Total transgenic explants on plate ',
        #             randomization_datasheet$ID[i],
        #             ' and tissue ',
        #             pixel_demographics$Tissue[j],
        #             ' is ',
        #             total_transgenic))

        randomization_datasheet[i, eval(column_of_interest)] <- total_transgenic

        column_of_interest <- paste0('portion_transgenic_',
                                     pixel_demographics$Tissue[j])

        # print(paste0('Total transgenic explants on plate ',
        #             randomization_datasheet$ID[i],
        #             ' and tissue ',
        #             pixel_demographics$Tissue[j],
        #             ' is ',
        #             total_transgenic))

        # Adding this patch because otherwise we get all NaN for min and
        #  max of each trait in debugging stdout (debugging needs only...?)
        if(randomization_datasheet$total_explants[i]==0){
          randomization_datasheet[i, eval(column_of_interest)] <- NA
        }
        if(randomization_datasheet$total_explants[i]>0){
          randomization_datasheet[i, eval(column_of_interest)] <-
            total_transgenic /
            randomization_datasheet$total_explants[i]
        }

    }
}

if(opt$debug == TRUE){
  print(paste0("Maximum # grid positions with transgenic stem in any plate: ",
               max(na.omit(randomization_datasheet$n_transgenic_Stem))))
  print(paste0("Maximum # grid positions with transgenic callus in any plate: ",
               max(na.omit(randomization_datasheet$n_transgenic_Callus))))
  print(paste0("Maximum # grid positions with transgenic shoot in any plate: ",
               max(na.omit(randomization_datasheet$n_transgenic_Shoot))))
}

# Calculate total ---------------------------------------------------------
## Now calculate all of each tissue whether transgenic or not
# Initialize new columns
for (j in 1:nrow(pixel_demographics)){
  if(opt$debug == TRUE) print(pixel_demographics[j, ])
    new_column_name <- paste0('n_', pixel_demographics$Tissue[j])

    randomization_datasheet <- cbind(
      randomization_datasheet,
      as.numeric(rep(NA, nrow(randomization_datasheet))))

    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <-
      new_column_name

    new_column_name <- paste0('portion_', pixel_demographics$Tissue[j])

    randomization_datasheet <- cbind(
      randomization_datasheet,
      as.numeric(rep(NA, nrow(randomization_datasheet))))

    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <-
      new_column_name
}

for (i in 1:nrow(randomization_datasheet)){
    for (j in 1:nrow(pixel_demographics)){
        column_of_interest <- paste0('n_', pixel_demographics$Tissue[j])

        total_tissue <- sum(
          na.omit(output$segment_present[which(
            output$ID==randomization_datasheet$ID[i] &
              output$segment_hex == pixel_demographics$Tissue[j])]))

        randomization_datasheet[i, eval(column_of_interest)] <- total_tissue
        if(total_tissue>20){

          if(opt$debug == TRUE){
            cat("\n")
            print(paste0("This ID is : ", randomization_datasheet$ID[i]))
            print(paste0("This segment_hex is : ", pixel_demographics$Tissue[j]))
            print(paste0("The rows we are summing over from output",
                         "to calculate this obviously wrong statistic:"))

            print(output$segment_present[which(
              output$ID==randomization_datasheet$ID[i] &
                output$segment_hex == pixel_demographics$Tissue[j])])
          }

          stop(paste("There are more than 20 grid positions with",
                     "this tissue on this plate? This can't be right."))
        }

        column_of_interest <- paste0('portion_', pixel_demographics$Tissue[j])

        total_tissue <- sum(na.omit(output$segment_present[which(
          output$ID==randomization_datasheet$ID[i] &
            output$segment_hex == pixel_demographics$Tissue[j])]))

        randomization_datasheet[i, eval(column_of_interest)] <-
          total_tissue / randomization_datasheet$total_explants[i]


    }
}

if(opt$debug == TRUE){
  print(paste0("Minimum # grid positions with explants in any plate: ",
               min(na.omit(randomization_datasheet$total_explants))))
  print(paste0("Maximum # grid positions with explants in any plate: ",
               max(na.omit(randomization_datasheet$total_explants))))
  print(paste0("Maximum # grid positions with callus in any plate: ",
               max(na.omit(randomization_datasheet$n_Callus))))
  print(paste0("Maximum # grid positions with shoot in any plate: ",
               max(na.omit(randomization_datasheet$n_Shoot))))
}

if(max(na.omit(randomization_datasheet$n_Callus))>20){
  stop(paste("Error! Why are there more than 20 calli",
             "on a plate? Are there two images per sample?"))
}

if(opt$debug == TRUE){
  for (j in 1:nrow(pixel_demographics)){
    column_of_interest <- paste0('portion_transgenic_',
                                 pixel_demographics$Tissue[j])

    print(paste0("Column of interest is ", column_of_interest))
    print(paste0("For transgenic tissue,",
                 "What is the range of values for portion for ",
                 pixel_demographics$Tissue[j]))

    print(min(randomization_datasheet[, get(column_of_interest)],
              na.rm = TRUE))

    print(max(randomization_datasheet[, get(column_of_interest)],
              na.rm = TRUE))

    print(paste0("What are factor levels?"))
    print(levels(factor(randomization_datasheet[, get(column_of_interest)])))
    cat("\n")

    column_of_interest <- paste0('portion_', pixel_demographics$Tissue[j])
    print(paste0("Column of interest is ", column_of_interest))

    print(paste0("For total tissue, What is the",
                 "range of values for portion for ",
                 pixel_demographics$Tissue[j]))

    print(min(randomization_datasheet[, get(column_of_interest)],
              na.rm = TRUE))
    print(max(randomization_datasheet[, get(column_of_interest)],
              na.rm = TRUE))

    print(paste0("What are factor levels?"))
    print(levels(factor(randomization_datasheet[, get(column_of_interest)])))
    cat("\n")
  }
}

## Make plots
if(opt$debug == TRUE){
  print("N rows of randomization datasheet: ")
  print(nrow(randomization_datasheet))
}

randomization_datasheet <- randomization_datasheet[!is.na(
  randomization_datasheet$Genotype_ID),]

if(opt$debug == TRUE){
  print("N rows of randomization datasheet: ")
  print(nrow(randomization_datasheet))

  cat("\n")
  print("Plots will be saved to: ")
  print(wd)
  cat("\n")
}

# Calculate size of plot output -------------------------------------------

## Before we make any plots, let's make sure their size will be ok
## We want a plot size that gives 1/4in per genotype*treatment

if(opt$width == -9){
  n_genotypes <- length(levels(factor(randomization_datasheet$Genotype_ID)))
  n_treatments <- length(levels(factor(randomization_datasheet$Treatment)))

  # We want the plots to be no smaller than 7 in wide to allow for titles
  plot_horz_in <- max(c(n_genotypes*n_treatments/3),
                      7)
}
if(opt$width != -9){
  plot_horz_in <- opt$width
}

if(plot_horz_in > 10){
  angle <- 75
}
if(plot_horz_in <= 10){
  angle <- 0
}

# Sort genotypes for plots ------------------------------------------------

# if(opt$sort==1){
#   colnames(randomization_datasheet) <- gsub("Treatment name", "Treatment",
#                                             colnames(randomization_datasheet))
#
#   rd_to_melt <- data.frame(
#     cbind(randomization_datasheet$Genotype_ID,
#           randomization_datasheet$Treatment,
#           randomization_datasheet$portion_transgenic_Callus))
#
#   colnames(rd_to_melt) <- c("Genotype_ID", "Treatment",
#                             "portion_transgenic_Callus")
#
#   rd_to_melt$portion_transgenic_Callus <-
#     as.numeric(as.character(rd_to_melt$portion_transgenic_Callus))
#
#   rd_to_melt_aggregated <- aggregate(
#     portion_transgenic_Callus ~ Genotype_ID + Treatment,
#     data = rd_to_melt, FUN = mean, drop = TRUE)
#
#   rd_to_melt_aggregated_spread <- spread(rd_to_melt_aggregated,
#                                          key = "Treatment",
#                                          value = "portion_transgenic_Callus")
#
#   rd_to_melt_aggregated_spread$difference <-
#     rd_to_melt_aggregated_spread[,ncol(rd_to_melt_aggregated_spread)] -
#     rd_to_melt_aggregated_spread[,2]
#
#   ordered_genotypes <- rd_to_melt_aggregated_spread[order(
#     rd_to_melt_aggregated_spread$difference),]$Genotype_ID
#
#   randomization_datasheet$Genotype_ID <-
#     factor(randomization_datasheet$Genotype_ID,
#            levels = ordered_genotypes)
#
#   colnames(randomization_datasheet) <- gsub("Treatment", "Treatment name",
#                                             colnames(randomization_datasheet))
# }


if(!dir.exists("Perc of explants with X")){
  dir.create("Perc of explants with X")
}

setwd("Perc of explants with X")

# By-plate plots  ----------------------------------------------------
# Initialize a counter for plot IDs
plot_counter <- 1

generate_plot <- function(data, tissue_type, variant, plot_type, plot_horz_in, plot_height, reporter, plot_id) {
  y_value <- ifelse(variant == "total",
                    paste0("portion_", tissue_type),
                    paste0("portion_transgenic_", tissue_type))

  y_title <- gsub("_", " ", paste0("Grid positions with ", tissue_type))
  main_title <- gsub("_", " ", paste0("Rates of ", variant, " ", tissue_type, " regeneration"))

  base_plot <- ggplot(data, aes_string(x = "`Treatment name`", y = y_value, group = "`Treatment name`")) +
    theme_dark() +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    stat_summary(fun = mean, geom = "point", shape = 1, size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(size = rel(1.3)),
          plot.title = element_text(size = rel(1.7))) +
    ylab(tools::toTitleCase(y_title)) +
    ggtitle(tools::toTitleCase(main_title))

  if (plot_type == "boxplot") {
    plot <- base_plot +
      geom_boxplot(outlier.shape = NA, fill = randomColor(luminosity = c("bright")), color = "white", alpha = 0.2) +
      facet_grid(~`Genotype_ID`) +
      geom_jitter(width = 0.25, height = 0.001, size = 1)
  } else {
    plot <- base_plot +
      geom_violin(fill = randomColor(luminosity = c("bright")), color = "white", alpha = 0.2) +
      geom_jitter(width = 0.25, height = 0.001, size = 1, aes_string(color = "`Genotype_ID`")) +
      scale_color_brewer(palette = "Set1") +
      theme(legend.position = "right")
  }

  folder_name <- if (plot_type == "boxplot") "By_Genotype" else "Genotypes_Combined"
  file_name_suffix <- ifelse(variant == "total", "_total", "_transgenic")
  file_name <- paste0("./", folder_name, "/", plot_id, "_", tissue_type, file_name_suffix, "_", plot_type, ".png")

  if (!dir.exists(folder_name)) dir.create(folder_name)
  ggsave(file_name, plot = plot, width = plot_horz_in, height = plot_height, units = "in")
}

# Loop over tissues, variants, and plot types
variants <- c("total", "transgenic")
plot_types <- c("boxplot", "violin")
plot_height <- opt$height
reporter <- opt$Reporter

for (tissue in pixel_demographics$Tissue) {
  for (variant in variants) {
    for (plot_type in plot_types) {
      # Construct a sequential and meaningful plot ID
      plot_id <- paste0(plot_counter, ifelse(plot_type == "boxplot", "A", "B"), ifelse(variant == "total", "i", "ii"))
      generate_plot(randomization_datasheet, tissue, variant, plot_type, plot_horz_in,
                    plot_height, reporter, plot_id)
      plot_counter <- plot_counter + 1
    }
  }
}



# Change dir for raw signal plots --------------------------------------------------------------

setwd("../")

if(!dir.exists("Raw signal stats")){
  dir.create("Raw signal stats")
}

setwd("Raw signal stats")

# Per explant signal plots ------------------------------------------------

if(opt$debug == TRUE){
  cat("\n")
  print("List of treatments in this randomization_datasheet: ")
  print(levels(factor(randomization_datasheet$`Treatment name`)))
  print("Now about to merge randomization_datasheet and GMOlabeler output.")
  print(paste0("Rows in randomization_datasheet: ",
               nrow(randomization_datasheet)))
  print(paste0("Rows in GMOlabeler output: ", nrow(output)))

  print("Before merging, look at table of ID in both. ")
  print(table(randomization_datasheet$ID))
  print(table(output$ID))
}

output <- output[!is.na(output$ID), ]
output <- output[!is.na(output$grid_item), ]

# The cartesian merge error here can result if we include missing explant
#  analysis results for a given plate twice. In those cases, the `pre_label.R`
#  output needs to be loaded to make sure redundant samples or problematic imgs
#  are excluded.
combined_data <- merge(output, randomization_datasheet, by="ID",
                       all.x=TRUE, allow.cartesian = FALSE)

if(opt$debug == TRUE){
  print(paste0("Rows in merged output: ", nrow(combined_data)))
  print("List of treatments in this combined data: ")
  print(levels(factor(combined_data$`Treatment name`)))
  cat("\n")
}

combined_data$Genotype_ID <- as.factor(combined_data$Genotype_ID)

# When not all images are in randomization datasheet
#  (e.g. some excluded due to being certain explant type),
#  this line needed to stop NA from being plotted.
combined_data <- combined_data[which(!is.na(combined_data$`Treatment name`)), ]
for(i in 1:(nrow(pixel_demographics)-1)){

  if(opt$debug == TRUE){
    print(paste0("Making final plots for tissue: "))
    print(pixel_demographics[i, ])
    print("Dim before and after subsetting combined data to this tissue only")
    print(dim(combined_data))
  }

  data_subset <- combined_data[which(
    combined_data$segment_hex == pixel_demographics$Tissue[i]),]

  if(opt$debug == TRUE){
    print(dim(data_subset))
    print("Max mean_signal")
    print(max(na.omit(data_subset$mean_signal)))
  }

  if (is.infinite(max(na.omit(data_subset$mean_signal)))){
    next
    }

  if(pixel_demographics$Tissue[i]!="All_regenerated_tissue" &
     pixel_demographics$Tissue[i]!="All_tissue"){
    # We do not have means calculated for these.
    if(opt$debug == TRUE) print("Mean signal plots")
    p1 <- ggplot(data_subset, aes(x=`Treatment name`, y=mean_signal,
                                  group=`Treatment name`)) +
      geom_boxplot(outlier.shape=NA,
                   fill=randomColor(luminosity=c("bright")),
                   color="black",
                   alpha=0.2) +
      theme_gray() +
      geom_jitter(width=0.05, height=0.01, size = 1) +
      facet_grid(~`Genotype_ID`) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold",
                                       size = rel(1.4)),
            axis.title.x = element_text(size = rel(1.3)),
            axis.text.y = element_text(size = rel(1.5), vjust=-1),
            axis.title.y = element_text(size = rel(1.3)),
            strip.text = element_text(angle = angle, size=rel(1.3)),
            plot.title = element_text(size=rel(1.7))) +
      ylab(paste0("Mean ", opt$Reporter, " signal")) +
      #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
      stat_summary(fun.y=mean, geom="point", shape = 1, size = 3) +
      ggtitle(paste0("Mean ", opt$Reporter, " signal in ",
                     tolower(pixel_demographics$Tissue[i]))) +
      scale_y_continuous(labels = scales::scientific)

    if(opt$debug == TRUE){
      print("Finished making mean signal plot... Now to save and make the rest.")
      print(p1)
    }

    ggsave(
      paste0("C",i,"_", pixel_demographics$Tissue[i],"_Mean_signal.png"),
      plot = last_plot(), width = plot_horz_in, height = opt$height,
      units = "in")
  }

  if(opt$debug == TRUE) print("Max signal plots")

  p2 <- ggplot(data_subset, aes(x=`Treatment name`, y=max_signal,
                                group=`Treatment name`)) +
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold",
                                     size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab(paste0("Max ", opt$Reporter, " signal")) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    stat_summary(fun.y=mean, geom="point", shape = 1, size = 3) +
    ggtitle(paste0("Max ", opt$Reporter, " signal in ",
                   tolower(pixel_demographics$Tissue[i])))+
    scale_y_continuous(labels = scales::scientific)

  if(opt$debug == TRUE) print(p2)

  ggsave(
    paste0("B",i,"_", pixel_demographics$Tissue[i],"_Max_signal.png"),
      plot = last_plot(),   width = plot_horz_in,   height = opt$height,
    units = "in")

  if(opt$debug == TRUE) print("Total signal plots")

  p3 <- ggplot(data_subset, aes(x=`Treatment name`, y=total_signal,
                                group=`Treatment name`)) +
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold",
                                     size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab(paste0("Total ", opt$Reporter, " signal")) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    stat_summary(fun.y=mean, geom="point", shape = 1, size = 3) +
    ggtitle(paste0("Total ", opt$Reporter, " signal in ",
                   tolower(pixel_demographics$Tissue[i]))) +
    scale_y_continuous(labels = scales::scientific)

  if(opt$debug == TRUE) print(p3)

  ggsave(
    paste0("A",i,"_", pixel_demographics$Tissue[i],"_Total_signal.png"),
      plot = last_plot(), width = plot_horz_in, height = opt$height,
    units = "in")

  p4 <- ggplot(data_subset, aes(x=`Treatment name`, y=total_pixels,
                                group=`Treatment name`)) +
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold",
                                     size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab(paste0("Total ", tolower(pixel_demographics$Tissue[i]), " pixels")) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    stat_summary(fun.y=mean, geom="point", shape = 1, size = 3) +
    ggtitle(paste0("Total pixels of ", tolower(pixel_demographics$Tissue[i]))) +
    scale_y_continuous(labels = scales::scientific)

  if(opt$debug == TRUE) print(p4)

  ggsave(
    paste0("A",i,"_", pixel_demographics$Tissue[i],"_Total_pixels.png"),
    plot = last_plot(), width = plot_horz_in, height = opt$height, units = "in")

  p5 <- ggplot(data_subset, aes(x=`Treatment name`,
                                y=n_pixels_passing_threshold,
                                group=`Treatment name`)) +
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold",
                                     size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab(paste0("Total ", opt$Reporter, "+ ",
                tolower(pixel_demographics$Tissue[i]), " pixels")) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    stat_summary(fun.y=mean, geom="point", shape = 1, size = 3) +
    ggtitle(paste0("Total ", opt$Reporter, "+ pixels in ",
                   tolower(pixel_demographics$Tissue[i]))) +
    scale_y_continuous(labels = scales::scientific)

  if(opt$debug == TRUE) print(p5)

  ggsave(
    paste0("A",i,"_", pixel_demographics$Tissue[i],
           "_pixels_w_strong_signal.png"),
    plot = last_plot(), width = plot_horz_in, height = opt$height, units = "in")
}



# Write out summary statistics ------------------

# Transformation/regeneration frequency summary stats

setwd("../")

if(opt$debug == TRUE) print(paste0("Saving summary stats in folder: ", getwd()))

fwrite(randomization_datasheet, "summary_stats_by_plate.csv")

if(opt$debug == TRUE) print("Saved summary_stats_by_plate.csv")

# Signal summary stats

colnames(combined_data) <- gsub("filename", "filename.x",
                                colnames(combined_data))

explant_data <- cbind(combined_data$ID,
                      combined_data$grid_item,
                      combined_data$Genotype_ID,
                      combined_data$`Treatment name`,
                      combined_data$filename.x,
                      combined_data$segment_hex,
                      combined_data$segment_present,
                      combined_data$mean_signal,
                      combined_data$max_signal,
                      combined_data$total_signal,
                      combined_data$total_pixels,
                      combined_data$n_pixels_passing_threshold)

colnames(explant_data) <- c("ID",
                            "Grid position",
                            "Genotype ID",
                            "Treatment",
                            "Filename",
                            "Tissue segment",
                            "Segment present",
                            "Mean signal",
                            "Max signal",
                            "Total signal",
                            "Total pixels",
                            "Total pixels with signif. signal")

if(sum(grepl("Tissue type", colnames(combined_data))) == 1){
  explant_data <- cbind(combined_data$`Tissue type`,
                        explant_data)

  colnames(explant_data)[1] <- c("Explant type")

}

if(sum(grepl("Block", colnames(combined_data))) == 1){
  explant_data <- cbind(combined_data$Block,
                        explant_data)

  colnames(explant_data)[1] <- c("Block")

}



explant_data <- as.data.frame(explant_data)

explant_data$`Segment present` <- as.logical(explant_data$`Segment present`)

fwrite(explant_data, "summary_stats_by_explant.csv")
if(opt$debug == TRUE) print("Saved summary_stats_by_explant.csv")


# Perform statistical tests -----------------------------------------------
# New section being added 4/19/22
# We will do tests over whole-plate statistics (as % of grid items)
#   in `randomization_datasheet` as well as over grid-item statistics
#   (fluorophore signal raw stats) in `combined_data`

if(opt$debug == TRUE) print("Completed workflow up until statistical tests")

library(nlme)
library(emmeans)
library(plyr)

# First, for whole-plate stats.
traits <- colnames(randomization_datasheet)[
  grepl("portion",
        colnames(randomization_datasheet))]

colnames(randomization_datasheet) <- gsub(" ", "_", colnames(
  randomization_datasheet
))

# randomization_datasheet_dummy <- randomization_datasheet
# randomization_datasheet_dummy$Block <- 1
#
# emmeans_table <- function(lme_object){
#   emm1 = emmeans(lme_object, specs = pairwise ~ Treatment_name)#,
#                  #regrid = "response") # doesn't work with regrid or transform
#   contrasts <- emm1$contrasts
#   contrasts_w_CI <- confint(contrasts)
#   contrasts_w_CI$p <- summary(contrasts, infer = TRUE)$p.value
#   contrasts_w_CI
# }

# i <- 1
# trait <- traits[i]
#
# if(length(levels(factor(randomization_datasheet$Block))) > 1){
#   formula <- as.formula(paste(traits[i], "~ Treatment_name + Block"))
# }
#
# linreg_perms <- perm.lm(formula = formula,
#                         data = randomization_datasheet,
#                         nperm = 1000,
#                         type = "regression",
#                         progress = TRUE)


### Deploy over plate traits (version if only one genotype)

ptm <- proc.time()

all_anova_plate_level <- data.frame()
all_linreg_plate_level <- data.frame()

for(i in 1:length(traits)){

  if(opt$debug == TRUE){
    print(traits[i])
    if(length(levels(factor(randomization_datasheet[, get(traits[i])]))) == 1){

      print("Single value")
      print(length(levels(factor(randomization_datasheet[, get(traits[i])]))))
      next # Skip trait if everything has same value
    }
    print(paste0("Trait has values: "))
    print(levels(factor(randomization_datasheet[, get(traits[i])])))
  }

  if(length(levels(factor(randomization_datasheet$Block))) > 1){
    formula <- paste(traits[i], "~ Treatment_name + Block")
  }

  if(length(levels(factor(randomization_datasheet$Block))) == 1){
    formula <- paste(traits[i], "~ Treatment_name")
  }

  if(length(levels(factor(randomization_datasheet$Genotype_ID))) > 1){
    formula <- paste(formula, "+ (1 | Genotype_ID)")
  }

  formula <- as.formula(formula)

  if(opt$debug == TRUE) print(formula)

  if(length(levels(factor(randomization_datasheet$Genotype_ID))) > 1){
    linreg_perms <- perm.lmer(formula = formula,
                              data = randomization_datasheet,
                              nperm = 1000,
                              type = "regression",
                              progress = TRUE)

    anova_perms <- perm.lmer(formula = formula,
                              data = randomization_datasheet,
                              nperm = 1000,
                              type = "anova",
                              progress = TRUE)
  }

  if(length(levels(factor(randomization_datasheet$Genotype_ID))) == 1){
    linreg_perms <- perm.lm(formula = formula,
                            data = randomization_datasheet,
                            nperm = 1000,
                            type = "regression",
                            progress = TRUE)

    anova_perms <- perm.lm(formula = formula,
                           data = randomization_datasheet,
                           nperm = 1000,
                           type = "anova",
                           progress = TRUE)
  }

  linreg_perms$trait <- traits[i]
  anova_perms$trait <- traits[i]

  if(opt$debug == TRUE) print(linreg_perms)

  all_anova_plate_level <- rbind.fill(all_anova_plate_level, anova_perms)
  all_linreg_plate_level <- rbind.fill(all_linreg_plate_level, linreg_perms)
}

traits <- c("mean_signal",
            "max_signal",
            "total_signal",
            "total_pixels",
            "n_pixels_passing_threshold")

all_anova_explant_level <- data.frame()
all_linreg_explant_level <- data.frame()

colnames(combined_data) <- gsub(" ", "_",
                                colnames(combined_data))

for(i in 1:(nrow(pixel_demographics)-1)){

  data_subset <- combined_data[which(
    combined_data$segment_hex == pixel_demographics$Tissue[i]),]

  if(opt$debug == TRUE) print(pixel_demographics$Tissue[i])

  if (is.infinite(max(na.omit(data_subset$mean_signal)))){
    next
  }

  if(pixel_demographics$Tissue[i]!="All_regenerated_tissue" &
     pixel_demographics$Tissue[i]!="All_tissue"){

    for(j in 1:length(traits)){
      print(traits[j])

      # If only one treatment has non-NA data for trait x tissue, skip trait
      if(length(levels(factor(data_subset$Treatment_name[which(
        !is.na(data_subset[, get(traits[j])]))]))) == 1){
        next
      }

      if(length(levels(factor(randomization_datasheet$Block))) > 1){
        formula <- paste(traits[j], "~ Treatment_name + Block")
      }

      if(length(levels(factor(randomization_datasheet$Block))) == 1){
        formula <- paste(traits[j], "~ Treatment_name")
      }

      if(length(levels(factor(randomization_datasheet$Genotype_ID))) > 1){
        formula <- paste(formula, "+ (1 | Genotype_ID)")
      }

      formula <- as.formula(formula)

      if(length(levels(factor(randomization_datasheet$Genotype_ID))) > 1){
        linreg_perms <- perm.lmer(formula = formula,
                                  data = data_subset,
                                  nperm = 1000,
                                  type = "regression",
                                  progress = TRUE)

        anova_perms <- perm.lmer(formula = formula,
                                 data = data_subset,
                                 nperm = 1000,
                                 type = "anova",
                                 progress = TRUE)
      }

      if(length(levels(factor(randomization_datasheet$Genotype_ID))) == 1){
        linreg_perms <- perm.lm(formula = formula,
                                data = data_subset,
                                nperm = 1000,
                                type = "regression",
                                progress = TRUE)

        anova_perms <- perm.lm(formula = formula,
                               data = data_subset,
                               nperm = 1000,
                               type = "anova",
                               progress = TRUE)
      }

      linreg_perms$trait <- traits[j]
      anova_perms$trait <- traits[j]
      linreg_perms$tissue <- pixel_demographics$Tissue[i]
      anova_perms$tissue <- pixel_demographics$Tissue[i]

      all_anova_explant_level <- rbind.fill(all_anova_explant_level,
                                            anova_perms)
      all_linreg_explant_level <- rbind.fill(all_linreg_explant_level,
                                             linreg_perms)

    }
  }
}

finish_time <- proc.time() - ptm
if(opt$debug == TRUE) print(finish_time)

fwrite(all_anova_plate_level, "ANOVA_plate_stats.csv")
fwrite(all_linreg_plate_level, "Linreg_plate_stats.csv")

fwrite(all_anova_explant_level, "ANOVA_explant_statistics.csv")
fwrite(all_linreg_explant_level, "Linreg_explant_statistics.csv")


