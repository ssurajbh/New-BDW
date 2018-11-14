##### Packages to Install #####
install.packages("tidyverse")
install.packages("readr")
install.packages("dplyr")
install.packages("fansi")
install.packages("stringi")
install.packages("ape")
install.packages("RSQLite")
install.packages("DECIPHER")
install.packages("stringr")
source("https://bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
abiocLite("muscle")
biocLite("msa")
biocLite("DECIPHER")

##### Data Analysis ####

############### Importing the data file



#bold_data <- read.delim(file.choose(), header = T)
#bold_data



library(tidyverse)

bold_data <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Pseudomonas&format=tsv") 

#N Loading selected data set of Daphnia specimen using Bold Api


#N The Data set checking and analysis works well but coild be expanded upon


############### Checking data-sets

############### Class and summary
class(bold_data) #tbl data.frame
summary(bold_data) #540 entries
str(bold_data)

############### Variable names.
names(bold_data)

############### Species richness of data 

length(unique(bold_data$species_name)) #40 

#N The BIN URI markercode was a good option for the original dataset, however the genus Pseudomonas lacks BIN URIs, they are NA in the data. As an alternative we can pick 

############### No. of BINs 

length(unique(bold_data$bin_uri))

############### Specimen bearing BINS.

length(grep("[:]", bold_data$bin_uri)) 

############### Data Reduction

# data3 <- bold_data[c(grep(":", bold_data$bin_uri)),]
# 
# length(data3$bin_uri)
# 
# sum(is.na(data3$bin_uri))      
# 
# sum(is.na(bold_data$bin_uri))


#N selecting data with bin_uris 

library(dplyr)
bold_data %>%
  select(bin_uri) -> data_bold #select function to extract wanted marker columns


############## Country with most barcode data 


library(fansi)

bold_data %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))


#United States 54



############### Biodiversity analysis

install.packages("vegan")
library(vegan)

### Groupind data and counting records in each bin

data_bold <- data_bold %>%
  group_by(bin_uri) %>%
  count(bin_uri)

############## changing the data set format

library(dplyr)

library(magrittr)

library(tidyverse)

spread_data <- spread(data_bold, bin_uri, n)


############## checking unique markers

unique(bold_data$markercode)





library(stringi)
library(ape)
library(RSQLite)
library(stringr)
library(Biostrings)
library(muscle)
library(msa)
library(DECIPHER)



############### Database, with list of markercodes

MarkerCodeList <- bold_data %>%
  group_by(markercode) %>%
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

############## Converting to data string

bold_dataString <- DNAStringSet(bold_data$nucleotides)


############### Subseting the dataframe bold_data to retain those records having a COI-5P markercode.

DataCOI <- bold_data %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect('nucleotides', "[ACGT]"))
DataCOI

############### Calculating Nucleotide sequences.

######## checking data type
class(DataCOI$nucleotides)

######## Converting to 18s markerset
Data18snucleotides <- DNAStringSet(DataCOI$nucleotides)
Data18snucleotides
