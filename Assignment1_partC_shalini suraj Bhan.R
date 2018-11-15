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
install.packages("vegan")
source("https://bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
abiocLite("muscle")
biocLite("msa")
biocLite("DECIPHER")
install.packages('rworldmap')

#### Data Acquisition ####

#N Loading selected data set of Daphnia specimen using Bold Api
library(tidyverse) %>% bold_data <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Pseudomonas&format=tsv")

##### Data Analysis ####

##### Part 1: Data Filtering ####

library(stringi) 
library(ape)
library(RSQLite)
library(stringr)
library(Biostrings)
library(muscle)
library(msa)
library(DECIPHER)

############### Checking data-sets

class(bold_data) #tbl data.frame
summary(bold_data) #540 entries
str(bold_data)

############### Variable names.
names(bold_data)

############### Species richness of data 

length(unique(bold_data$species_name)) #40 

#Using the select function instead of grepping 

#N looking at NA values

sum(is.na(bold_data$bin_uri))


#N selecting data 

library(dplyr)
bold_data %>%
  select(processid, species_name, nucleotides, country) -> bold_data2 #select function to extract wanted marker columns

length(is.na(bold_data2)) #identify NA values

bold_data2 <- na.omit(bold_data2) #N Option to remove NA data

############## checking unique markers

unique(bold_data$markercode)

############### Database, with list of markercodes

MarkerCodeList <- bold_data %>%
  group_by(markercode) %>%
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

############### Subseting the dataframe bold_data to retain those records having a COI-5P markercode.

DataCOI <- bold_data %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect('nucleotides', "[ACGT]"))
DataCOI


##### Part 2 Biodiversity ####

#N Country with most barcode data gives top 10

bold_data2 %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))


#N Possible visualisation

library(rworldmap)
bold_dataframe <- as.data.frame(bold_data2)
countries <- as.data.frame(table(bold_dataframe$country))
colnames(countries) <- c("country", "value")
matched <- joinCountryData2Map(countries, joinCode="NAME", nameJoinColumn="country") #matching countries from data to those in the Rworldmap package: 26 matched, 217 did not 
mapCountryData(matched, nameColumnToPlot="value", mapTitle="Pseudomonas samples per country ", catMethod = "pretty", colourPalette = "white2Black")

library(vegan)

############## changing the data set format

library(dplyr)
library(magrittr)
library(tidyverse)


##### Part 3 Phylogeny ####

library(dplyr); bold_data2 %>% group_by(country) %>% sample_n(1) -> bold_sample #sample data by countries, to reduce dataset 



#N Converted to DNA string but no further analysis performed 

############## Converting to data string

bold_data_string <- DNAStringSet(bold_sample$nucleotides) #N No NA values accepted

#NPicking up from where the anaylsis left off using the DNA stringset to perform a muscle alignment 
#N create dendrogram distance matrix
#N Possibly Specaccum for biodiversity 

bold_data_alignment <- DNAStringSet(muscle::muscle(bold_data_string, maxiters = 2,diags = TRUE)) #running muscle MSA 

bold_data_DNAbin<- as.DNAbin(bold_data_alignment)

distanceMatrix1 <- dist.dna(bold_data_DNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE) #creating a a pairwise distance matrix using the dist.dna function and selecting the Tamura-Nei 93 (TN93) model which offers a substitution model of nucleotide evolution

clusters.Daphnia16S <- IdClusters(distanceMatrix1,method = "single", cutoff= 0.02, showPlot = TRUE, type = "both", verbose = TRUE) #plotting the clusters using IDclusters creating a distnace similarity matrix at a 0.02 cutoff to provide enough clusters without overdoing it. 


