# I've previously downloaded functional abundance data and metadata
# and precalculated a PCoA so you can quickly generate visualizations
# from it for the 1,600+ samples in the HMP jumstpart project:

# Here we use a few simple scripts to produce metadata colored PCoAs for this data
# Set the "HMP_large_example" subdirectory of this repsoitory as your working directory
setwd("~/git/HMP_mini_analysis/HMP_large_example")

# I downloaded, filtered, and normalised the abundance data.
# The abundance data were used to generat a flat PCoA
# The result is in the "  " file in this directory.

# I also downloaded the metadata like this. Note that this command is commented out.
# It would take a few minutes to run, depending on your internet connection and the 
# current workload on the MG-RAST servers
#      get_metadata(mgid_list="HMP_jumstart_ids.txt")
# The result is in the "" file in this directory

# Now we can use the abundance data and samples metadata to create a PCoA
# Of a single metadata value

# or the column header for that column
render_pcoa.v14(PCoA_in="HMP.Jumpstart.DESeq_normed.euclidean.PCoA", metadata_table="HMP_jumpstart_metadata.txt", metadata_column_index="env_package.data.body_site")
# another example
render_pcoa.v14(PCoA_in="HMP.Jumpstart.DESeq_normed.euclidean.PCoA", metadata_table="HMP_jumpstart_metadata.txt", metadata_column_index="sample.data.host_disease")

# or one for each metdata value --- note this is commented out -- it will take a while to run
#       render_pcoa.v14(PCoA_in="HMP.Jumpstart.DESeq_normed.euclidean.PCoA", metadata_table="HMP_jumstart_metadata.txt", use_all_metadata_columns=TRUE)
# and that your computer may not have enough memory for te task - I have found
# that it requires about 1-2Gb of RAM (True on Windows, OSX, and Unix/Linux)











############################################################################################################################
# Source this function to donwload accessory functions - in this example they are all hosted on github
source_https <- function(url, ...) {
  require(RCurl)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}
############################################################################################################################


############################################################################################################################
# a function to perform batch download of annotation abundance data from MG-RAST
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/matR_batch_dl.r")
# we do not use this function in this example but include it in case you want to download data for
# for than ~50 metagenomes at a time
# simple example - using all defaults (level 3 subsystem download)
#     matR_batch_dl(mgid_list="HMP_jumstart_ids.txt")
# The google group has a much more detailed example:
# https://groups.google.com/forum/#!topic/matr-forum/T-9q1evhlJ8
# Note that we have explicitly left out the steps needed to produce the calculated PCoA
# from the raw data -- you should be able to fill in the blanks with the 30 metagenome example.
# Still having trouble - use the google group for help:
# https://groups.google.com/forum/#!topic/matr-forum
############################################################################################################################

############################################################################################################################
# a function to download metadata from MG-RAST
source_https("https://raw.githubusercontent.com/DrOppenheimer/matR-apps/master/get_metadata.R")
# This function will take minutes to run. Simple example:
#     get_metadata(mgid_list="HMP_jumstart_ids.txt", output_name="HMP_jumpstart_metadata")
############################################################################################################################