<<<<<<< HEAD
<<<<<<< HEAD
##### Packages to Install #####
=======
#N To tidy up the script, I created a list of all the packages that the user is asked to install. 
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
=======
#N To tidy up the script, I created a list of all the packages that the user is asked to install. 
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f

install.packages("tidyverse")
install.packages("readr")
install.packages("dplyr")
<<<<<<< HEAD
<<<<<<< HEAD
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

#### Importing the data file ####

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

=======
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f

############### Importing the data file

#N renamde data2 into bold_data for clarity 

#bold_data <- read.delim(file.choose(), header = T)
#bold_data

#N to avoid having to download a data file and import it, we can simply use the read_tsv function and read one directly from a bold URL 

#N Shalini's plan for her resarch question was to analyse gut microbe data. She did so using 3784 bold samples including 127 unique species of the mammalian class, including 9 orders: Artiodactyla Carnivora Cetacea Chiroptera Didelphimorphia Lagomorpha Perissodactyla Primates Rodentia Soricomorpha). For the improvements I picked Pseudomonas as the chosen genus for analysis to help narrow down the results and focus on the original research question, though we could substitute it for any other bacterium present in the gut flora. 

library(tidyverse)

bold_data <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Pseudomonas&format=tsv") 

#NLoading selected data set of Daphnia specimen using Bold Api


#N The Data set checking and analysis works well but coild be expanded upon


############### Checking data-sets

############### Class and summary
<<<<<<< HEAD
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
class(bold_data) #tbl data.frame
summary(bold_data) #540 entries
str(bold_data)

############### Variable names.
names(bold_data)

############### Species richness of data 

length(unique(bold_data$species_name)) #40 

<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
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

install.packages("fansi")
library(fansi)

bold_data %>%
<<<<<<< HEAD
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))


<<<<<<< HEAD
<<<<<<< HEAD
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
#N Specaccum for biodiversity 


bold_data_alignment <- DNAStringSet(muscle::muscle(bold_data_string, maxiters = 2,diags = TRUE)) #running muscle MSA 


bold_data_DNAbin<- as.DNAbin(bold_data_alignment)

distanceMatrix1 <- dist.dna(bold_data_DNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE) #creating a a pairwise distance matrix using the dist.dna function and selecting the Tamura-Nei 93 (TN93) model which offers a substitution model of nucleotide evolution

clusters.Daphnia16S <- IdClusters(distanceMatrix1,method = "single", cutoff= 0.02, showPlot = TRUE, type = "both", verbose = TRUE) #plotting the clusters using IDclusters creating a distnace similarity matrix at a 0.02 cutoff to provide enough clusters without overdoing it. 


=======
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
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



install.packages("stringi")
library(stringi)

install.packages("ape")
library(ape)

install.packages("RSQLite")
library(RSQLite)

install.packages("DECIPHER")
install.packages("stringr")
library(stringr)

source("https://bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
abiocLite("muscle")
biocLite("msa")
biocLite("DECIPHER")
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
<<<<<<< HEAD
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
=======
>>>>>>> 0e05fd4172b572555e990a191b4d045590da1a3f
