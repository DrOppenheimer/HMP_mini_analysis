# Notes
# The easiest way to use this tutorial is probably to load it into RStudio
# Anything following a "#" symbol will be interpreted as comments.
# If you runto into trouble -- get multiple R errors, I highly recommand
# that you switch to R gui without R studio or R from the system prompt/terminal
# This is a convenient way to avoid limitations imposed by R Studio, and is frequently
# necessary if you are analyzing a large number of samples (100s to 1,000s)

############################################################################################################################
############################################################################################################################
### Load matR and accessorty functions
############################################################################################################################
############################################################################################################################
# First, check your R version
R.Version()
# You should always try to keep your R up to date. This tutorial will require version 3 or later.
############################################################################################################################
# uninstall current matR (This is only necessary if your current version was installed from CRAN 
# or any source other than the early release candidate hosted on github)
remove.packages("matR")
# You'll see an error message if you don't already have it installed
############################################################################################################################
# load some additional packages - you may have to install them first
install.packages("matlab")
library(matlab)
############################################################################################################################
# install the devtols package
install.packages("devtools")
library(devtools)
# use devtools to install the early releas candidate of matR from its git repository
# install_github(repo="MG-RAST/matR", dependencies=FALSE, ref="3d068d0c4c644083f588379ca111a575a409c797")
install_github(repo="MG-RAST/matR", dependencies=FALSE, ref="early-release")
library(matR)
dependencies()
############################################################################################################################
# Source this function to donwload accessory functions - in this example they are all hosted on github
source_https <- function(url, ...) {
  require(RCurl)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
############################################################################################################################
# a function to perofrm batch download of annotation abundance data from MG-RAST
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/matR_batch_dl.r")
# we do not use this function in this example but include it in case you want to download data for
# for than ~50 metagenomes at a time
############################################################################################################################
# a function to download metadata from MG-RAST
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/get_metadata.R")
############################################################################################################################
# a function to import metadata that has already been downloaded and written to a flat file (using get_metadata and export_data)
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/import_metadata.r")
############################################################################################################################
# function to produce raw pco (performs the calculations and produces a flat file output but does not create images/vizualizations)
# PCoAs can be calculated from raw or normalized annotation abundance data with a variety of distance/dissimilarity metrics.
source_https("https://raw.githubusercontent.com/MG-RAST/AMETHST/master/plot_pco.r")
############################################################################################################################
# function to generate metadata colored PCoAs using the outputs from the metadata and pco function above
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/plot_fun.12-10-13/pcoa/render_calculated_pcoa.dev.v14.r")
############################################################################################################################
# Function to perform filtering (removal of singletons and extreme low abundance entries) and normalization (with a number
# different methods including log transformation and standardization, negative bionomial based variance stabiliation(with
# the DESeq package), and quantile normalization )
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/normalize_fun.2-27-14/norm_redux.v5.r")
# MGRAST_preprocessing
############################################################################################################################
# function to perform calculations for and to generate heatmaps from raw or normalized annotation abundance data.
# The function can also use metadata to annotate the image automatically
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/heatmap_fun/heatmap_dendrogram.from_file.3-18-14.r")
############################################################################################################################
# function to perform statistical tests on normalized or raw annotation abundance data
# The function can use metadata to automatically segregate samples into groups for testing
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/stats_fun.2-27-14/matR_stats_from_files.r")
############################################################################################################################
# simple function to export annotation abundance data in a format that can easily be used to transport and re-import the data
# for more advanced R users -- this is just a wrapper for write.table with appropriate options selected
export_data <- function(data_object, file_name){
  write.table(data_object, file=file_name, sep="\t", col.names = NA, row.names = TRUE, quote = FALSE, eol="\n")
}
############################################################################################################################
# simple function to imprt annotation abundance data from tab delimited text
# function assumes R appropriate formatting:
#     The first field is empty
#     The remainder of the first row contains column headers (samples)
#     The remainder of the first column contains column headers (depending on the dataset, rows can represent taxonomic or function categories)
#     All other fields contain numerical data (float or integer) that correspond to annotation abundance measures for the indicated
#          sample(column) and category (row; depending on the dataset, rows can represent taxonomic or function categories)
# for more advanced R users -- this is just a wrapper for write.table with appropriate options selected
import_data <- function(file_name)
{
  data.matrix(read.table(file_name, row.names=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE))
}
############################################################################################################################
# Function to remove singleton and low abundance values
remove.singletons <- function (x, lim.entry, lim.row, debug=FALSE) {
  x <- as.matrix (x)
  x [is.na (x)] <- 0
  x [x < lim.entry] <- 0 # less than limit changed to 0
  #x [ apply(x, MARGIN = 1, sum) >= lim.row, ] # THIS DOES NOT WORK - KEEPS ORIGINAL MATRIX
  x <- x [ apply(x, MARGIN = 1, sum) >= lim.row, ] # row sum equal to or greater than limit is retained
  if (debug==TRUE){write.table(x, file="sg_removed.txt", sep="\t", col.names = NA, row.names = TRUE, quote = FALSE, eol="\n")}
  x  
}
############################################################################################################################

############################################################################################################################
############################################################################################################################
### Download and filterfunctional annotation abundance data for 30 samples selected from the HMP
### Data will be downloaded and then filtered with respect to average annotation e-value, average match length,
### and average match percent id. 
############################################################################################################################
############################################################################################################################

# In this example we will explore functional annotation abundance data for a subset of samples in the HMP project
# hosted on MG-RAST (URL)
# These data are public -- to use private data, you have to have access to the datasets with your account and
# need to make sure that you are using a valid key -- this is covered in other examples.

# got to the directory where we want to perform the analysis (This will vary by your machine)
# if you are using RStudio you can got to "Session" and then "Set Working Directory" and "Choose Directory"
# I would recommend copying that value here (it will be displayed in the R terminal once selected)
setwd("~/Documents/Projects/HMP/tutorial.HMP_subset.9-1-15")

# import a list of mgrast ids; note that ids should have the "mgm" prefix
my_ids <- readIDs("HMP.30_mgrast_ids.9-1-15")
# This will also work
# my_ids <- scan(file="HMP.30_mgrast_ids.9-1-15", what="character")

# Establish some variables that will be used below to define the data that we download
my_annot = "function"    # c("function, "organism") # i.e. functional or taxonomic abundance data
my_level = "level3"      # many choices - depends function and source
my_source = "Subsystems" # many possibilities - we receommend Subsystems if you aren't sure what to choose

# create view(s) for the annotation abundance counts as well as the additional values you will use to sort them.
my_views <- list(
  my_counts=     c(entry="counts",     annot=my_annot, level=my_level, source=my_source),
  my_avg_evalues=c(entry="evalue",     annot=my_annot, level=my_level, source=my_source),
  my_avg_length= c(entry="length",     annot=my_annot, level=my_level, source=my_source),
  my_avg_pid=    c(entry="percentid",  annot=my_annot, level=my_level, source=my_source)
)
# This view will download the annotation abundances (counts) as well as the corresponding
# average  e-values, average hit length, and average hit percent identity.

# perform download, creating a matR collection that will contain all 4 types of data for the selected metagenomes
my_collection <- collection(my_ids, my_views)
# Note that this function will take a minute or two to complete -- it will download the abundance, evalue, length, and
# percentd tables for all metagenomes in the list. It will not work for more than ~ 50 metagenomes at a time. For large
# groups of metagenomes you can use the batch downloader provided above

# Before proceeding to the next step, make sure that all of the data are downloaded.
# Issue these commands to check.
my_collection$my_counts
my_collection$my_avg_evalues
my_collection$my_avg_length # note that this is the length of the hit of the in silco translated protein
my_collection$my_avg_pid
# If you get a ..."data is pending" ... for any of the above, wait a few minutes and try again
# When you see a matrix of values -- the data have been successfully downloaded

# export the 4 matrix objects (not completely necessary - but makes things a little easier later)
my_counts.matrix      <- my_collection$my_counts
my_avg_evalues.matrix <- my_collection$my_avg_evalues
my_avg_length.matrix  <- my_collection$my_avg_length
my_avg_pid.matrix     <- my_collection$my_avg_pid

# you can also write the raw data to file (I HIGHLY recommend that you do)
# you can always use these files to perform the analyses below without having
# to wait for the full download from MG-RAST
export_data(my_counts.matrix, "my_counts.txt")
export_data(my_avg_evalues.matrix, "my_evalues.txt")
export_data(my_avg_length.matrix, "my_lengths.txt") 
export_data(my_avg_pid.matrix, "my_pids.txt")

# Before we produce any plots - we'll use a simple trick to get R's graph parameters back to their defaults
dev.off()
dev.new()
original.par <- par()

# for fun you can look at the distribution of each type of value in your data
# All values in all samples at the same time
split.screen(c(2,2))
screen(1)
hist(my_counts.matrix)
screen(2)
hist(my_avg_evalues.matrix)
screen(3)
hist(my_avg_length.matrix) # note that this is the length of the hit of the in silco translated protein
screen(4)
hist(my_avg_pid.matrix)

# sometimes looking at the values can give you a better idea of where you should place
# your cutoffs --- alternatively, if to many/few values are removed, this will allog you
# to see why
# Values separated by sample (boxplots)

par(original.par) # reload the default graphical parameters we defined above
split.screen(c(2,2))
screen(1)
boxplot(my_counts.matrix, las=2, main="counts")
screen(2)
boxplot(my_avg_evalues.matrix, las=2, main="avg_evalue")
screen(3)
boxplot(my_avg_length.matrix, las=2, main="avg_length(in silico translated protein hits)")
screen(4)
boxplot(my_avg_pid.matrix, las=2, main="avg_pid")
# if you see a "Error in plot.new() : figure margins too large" error, just increase the 
# RStudio window for plots to a larger size. Alternative, you can run all of this code in R ui
# or system command line to avoid limitations imposed by R Studio


# place my filter values in variables -- filters will remove all values above(above *_min) or below(below *_max) the indicated value
# Values are currently set to MG-RAST defaults; these will obviously need to be edited per your requirements.


# use a nested loop to go through all of the count values, replacing them with 0 if they do not meet the filter
# In each of the loops, abundance values with a filtered value (e.g. avg evalue) that does not pass the filter
# are replcaed with 0
# Then all of the abundance rows that contain only zeros are removed

# make a copy of the abundance values that will be modified by filtering

my_counts.filtered.matrix <- my_counts.matrix

# filter abundances by evalue (replace values with e-value > my_avg_evalue_exp_max with 0)
# first look at the distribution of e-value values to determine a reasonable cutoff
par(original.par) # reload the default graphical parameters we defined above
hist(my_avg_evalues.matrix)
# E-06 looks like a reasonable place to cut
my_avg_evalue_exp_max <- -06 # upper limit for evalue exponent
for (i in 1:nrow(my_counts.filtered.matrix)){
  for (j in 1:ncol(my_counts.filtered.matrix)){
    if( my_avg_evalues.matrix[i,j] > my_avg_evalue_max ){
      my_counts.filtered.matrix[i,j] <- 0
    }
  }
}

# filter abundances avg_length (replace values with average length < my_avg_length_min with 0)
par(original.par) # reload the default graphical parameters we defined above
hist(my_avg_length.matrix)
# Lengths (average database hit length(translated protein length)) are a little on the short side
# we will use a suitably liberal minimum value
my_avg_length_min <- 20 # avg length of hits
for (i in  1:nrow(my_counts.filtered.matrix)){
  for (j in 1:ncol(my_counts.filtered.matrix)){
    if( my_avg_length.matrix[i,j] < my_avg_length_min ){
      my_counts.filtered.matrix[i,j] <- 0
    }
  }
}

# filter abundances avg_pid (replace values with average < my_avg_pid_min with 0 )
par(original.par) # reload the default graphical parameters we defined above
hist(my_avg_pid.matrix)
# 75 appears to be a reasonable minimum value
my_avg_pid_min    <- 75 # MG-RAST does not have a pid filter by default
for (i in 1:nrow(my_counts.filtered.matrix)){
  for (j in 1:ncol(my_counts.filtered.matrix)){
    if( my_avg_pid.matrix[i,j] < my_avg_pid_min ){
      my_counts.filtered.matrix[i,j] <- 0
    }
  }
}

# Now that we have filtered the data - replacing undesirable values with 0,
# It is possible that entire rows (functions(level3 subsystems))
# have been replaced entirely with zeroes or small number of extremely small values.
# Here we filter out such rows and columns - then check to see how many were lost
# from the original raw data.
removeSg_valueMin     = 2 # lowest retained individual value (lower converted to 0)
removeSg_rowMin       = 4 # lowest retained row sum (lower, row is removed)
# remove singletons
my_counts.filtered.matrix <- remove.singletons(x=my_counts.filtered.matrix, lim.entry=removeSg_valueMin, lim.row=removeSg_rowMin)

# see if any rows or columns have been lost
dim(my_counts.matrix) ; dim(my_counts.filtered.matrix)
# we can see that 102 rows were removed due to presence of extremely low values

# sort the filtered data by row (not necessary, but I usually like to
my_counts.filtered.matrix <- my_counts.filtered.matrix[ order(rownames(my_counts.filtered.matrix)), ]

# sort the unfiltered data the same way for export below
my_counts.matrix  <- my_counts.matrix[ order(rownames(my_counts.matrix )), ]

# If you want to export any of the values (say as tab delimited text for excel etc, this is an easy way to do it
# load this function

export_data(data_object=my_counts.matrix, file_name="unfiltered_counts.txt")
export_data(data_object=my_counts.filtered.matrix, file_name="filtered_counts.txt")
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
### Normalize annotation abundance data
############################################################################################################################
############################################################################################################################
# Note that the normalization function has some overlap with the screening methods described in the previous 
# section. We ignore there here, and use the default noramlization procedure(DESeq_blind) for the sake of simplicity.
# You can normalize annotation abundance data from the file created above:
MGRAST_preprocessing(data_in=my_counts.filtered.matrix, data_type="r_matrix")
# or, if it is still loaded, from the R object that contains the same data
MGRAST_preprocessing("filtered_counts.txt", data_type="file")
# Note that the normalization function returns a lot of information -- this information is 
# also saved to a log ( my_counts.filtered.matrix.DESeq_blind.PREPROCESSED.txt.log )
# It also saves the regression used by DESeq to normalize the data ( my_counts.filtered.matrix.DESeq_regression.png )
# A matrix that contains the normalized values ( my_counts.filtered.matrix.DESeq_blind.PREPROCESSED )
# A file that contains the same normalized values ( filtered_counts.txt.DESeq_blind.PREPROCESSED.txt )
# It is always good practise to see of normalization eliminated any of rows or columns
dim(my_counts.filtered.matrix); dim(my_counts.filtered.matrix.DESeq_blind.PREPROCESSED)
# It did not
# Now we can take a look at the distribution of the original raw data and along
# with the data that have undergone filtering and normalization
par(original.par) # reload the default graphical parameters we defined above
split.screen(c(2,1))
screen(1)
hist(my_counts.matrix, main="distribution of original raw values")
screen(2)
hist(my_counts.filtered.matrix.DESeq_blind.PREPROCESSED, main="distribution of values after filter and norm")
# Looks like we have a bimodal ditribution of abundance values --- but that's a story for another day ...
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
### Download sample metadata
### Filter metadata by data retained  in the metadata filtering
############################################################################################################################
############################################################################################################################
# You can download the metadata for your samples like this:
my_metadata <- get_metadata(mgid_list="HMP.30_mgrast_ids.9-1-15", output_file="HMP_30_samples.metadata")
# This function can take a minute or two to run -- it performs an independent querry to retrieve the metadata for each sample
# in the id list, and then formats the metadata to make it more palatable for R and/or human viewing.
# For more advanced users -- this function retrieves the fully nested metadata object from the API and then 
# recursively flatens it into a simple list. Successive lists are merged by rownames (name of a metadata variable)
# If you leave off the "output_file" option, the metadata will be written to the my_metadata matrix object, but not
# to a file.

# When data are filtered like those above - it is possible to remove samples (values are all 0, or too low to be retained)
# In addition - it is usually convenient to make sure that data and metadata columns have the same
# ordering. Both problems can be addressed with a command like the following
my_metadata.reordered <- my_metadata[ (colnames(my_counts.filtered.matrix)), ] # Note that in the abundance table, samples
# are per column; in the metadata, samples are per row
# I'd recommend exporting this reordered metadata
export_data(my_metadata.reordered, "filtered_counts.metadata.txt")
############################################################################################################################
############################################################################################################################

# NOTE: From this point on, images will be created as files in your working directory,
# while RStudio is a great tool I've found that it does not handle production of multiple
# figures in a consistent way. 

# You'll be able to see the figures by opening them from your working diretory(use the "getwd()" command to remind
# yourself where that is)

############################################################################################################################
############################################################################################################################
### Use abundance data and metadata to create PCoAs
############################################################################################################################
############################################################################################################################
# PCoAs are a great tool to present an entire dataset in a dimensionally reduced form
# One of the problems we frequently run into is generating multiple PCoAs for the same data.
# (i.e.) To produce a single PCoA and then create multiple images from it, coloring the points
# with different metadata values. Here we use a function that automates this labrious process

# First - calculate the raw pcoa, we'll use euclidean distance and Bray-Curtis Similarity
# for the sake of comparison.

# With Eucludean:
MGRAST_plot_pco(file_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt", dist_method="euclidean")
# creates files:
# filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA # the PCoA results
# filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.DIST # the distance/dissimilarity matrix used to compute the PCoA
# With Bray_curtis
MGRAST_plot_pco(file_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt", dist_method="euclidean")
# creates file: 
# filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA # the PCoA results
# filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.DIST # the distance/dissimilarity matrix used to compute the PCoA
# You can see the other available metrics if you call the function without arguments
MGRAST_plot_pco()

# We can create visualizations of a PCoA without metadata
render_pcoa.v14(PCoA_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", color_list="red")

# We can create a visualization that uses a single column of the metadata
# reference by column name
render_pcoa.v14(PCoA_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", metadata_table="filtered_counts.metadata.txt", metadata_column_index="env_package.data.body_site")
# or column number (base 1 index)
render_pcoa.v14(PCoA_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", metadata_table="filtered_counts.metadata.txt", metadata_column_index=1)

# We can also automatically produce PCoAs colored for every metadata value we have
render_pcoa.v14(PCoA_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", metadata_table="filtered_counts.metadata.txt", use_all_metadata_columns=TRUE)
# Note all outputs will be in your working directory
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
### Use abundance data and metadata to create heatmap dendrograms
############################################################################################################################
# A heatmap dendrogram can be generated from raw or normalized data
# There are a number of options for both the distance/dissimilarity metric
# and the clustering algorithm

# We use defaults for distance/dissimilarity (euclidean) and clustering (ward)
heatmap_dendrogram.from_file(file_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt")

# We can also add a color bar generated from any of the metadata values
heatmap_dendrogram.from_file(file_in="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt",
                             metadata_table="filtered_counts.metadata.txt", 
                             metadata_column_index=1
)
# Note all outputs will be in your working directory ( filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.HD.png  )
# Note: that the heatmap dendrogram creates a text file that contains the data reordered as it appears in the
# heatmap dendrogram ( filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.HD_sorted.txt  )
############################################################################################################################
############################################################################################################################

############################################################################################################################
### Statistical subselection
############################################################################################################################
# We can use any metadata to generate sample groupings for which statistics can be calculated
# a number of statistical tests are available - use (  ) to see the current list.

# Here we will perform a Kruskal-Wallis (non parametric ANOVA) to identify the level 3 subsystems that exhibit the most 
# significant differences among the three body sampling locations

stats_from_files(data_table ="filtered_counts.txt.DESeq_blind.PREPROCESSED.txt",
                 metadata_table="filtered_counts.metadata.txt", 
                 metadata_column=1
                 )
# Note all outputs will be in your working directory
# Stat tests are sorted by ascending Benjamini & Hochberg FDR

# You can easily use statistical results to subselect the data -- here is an example

# import the results of the statistical test
my_stats <- import_data("filtered_counts.txt.DESeq_blind.PREPROCESSED.txt.STATS_RESULTS.txt")
# convert to a dta frame
my_stats.dataframe <- data.frame(my_stats)

# subselect the subsystems with FDR <= 0.001
my_stats.lt_0p001 <- my_stats.dataframe[ my_stats.dataframe$Kruskal.Wallis..fdr <= 0.01, ]

# convert back to matrix and export just the coloumns that have the data values 
# import and use dim on the original data to check the number of columns added by stats
my_original_data <- import_data("filtered_counts.txt.DESeq_blind.PREPROCESSED.txt")
dim(my_original_data)
# 30 columns so:
subselected_data <- as.matrix(my_stats.lt_0p001)[,1:30]
export_data(subselected_data, "FDR_subselected_data.txt")

# Now generate a heatmap from the statistically subselected data
heatmap_dendrogram.from_file(file_in="FDR_subselected_data.txt",
                             metadata_table="filtered_counts.metadata.txt", 
                             metadata_column_index=1
)
# FDR_subselected_data.txt.HD.png
# Note: that the heatmap dendrogram creates a text file that contains the data reordered as it appears in the
# heatmap dendrogram ( FDR_subselected_data.txt.HD_sorted.txt )
############################################################################################################################
############################################################################################################################





############################################################################################################################
############################################################################################################################




































