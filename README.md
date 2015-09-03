# HMP_mini_analysis
Code and results for an analysis of a subset of HMP data using R and MG-RAST

The only files you need are:

    tutorial.HMP_subset.9-1-15.R

    HMP.30_mgrast_ids.9-1-15

All other files are just examples of the outputs that are generated when you 
perform the analysis outlined in the R file applied to the mgids in the 
second. These are included as examples - or can be used as inputs if the API
is acting up or if you are conducting a tutorial with a large number of users
at the same time.

This repository also contains a folder "HMP_large_example"
It contains a precomputed functional PCoA (level 3 Subsystems) for 
1,606 samples in the HMP jumstart project.
See this link for futher details about the data: http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=385 

I used methods identical to those outlined in the "Mini" example to produce a PCoA
of the complete set of samples. You can easily generate visualizations from the 
precomputed PCoA and pre-downloaded metadata with the R script contained in the
directory.

Note that further examples and instructions can be found on the google group for matR:
https://groups.google.com/forum/#!forum/matr-forum

Cheers
Kevin P. Keegan
9-1-15
