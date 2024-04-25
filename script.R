heights <- c(63, 69, 60, 65, NA, 68, 61, 70, 61, 59, 64, 69, 63, 63, NA, 72, 65, 64, 70, 63, 65)
(heightsNoNA <- na.omit(heights))
median(heights, na.rm = TRUE)
length(heightsNoNA[heightsNoNA > 67])
numeric(3)
logical(3)
character(3)

download.file(url = "https://github.com/carpentries-incubator/bioc-intro/raw/main/episodes/data/rnaseq.csv",
              destfile = "data/rnaseq.csv")
rna <- read.csv("data/rnaseq.csv")
rna
head(rna)
View(rna)
str(rna)
str(rna["gene"])
str(rna[,"gene"])

rna_200 <- rna[200, ]
rna[nrow(rna), ]
tail(rna, 1)
rna_middle <- rna[nrow(rna)/2, ]
rna[-(7:nrow(rna)), ]

country_climate <- data.frame(
  country = c("Canada", "Panama", "South Africa", "Australia"),
  climate = c("cold", "hot", "temperate", "hot/temperate"),
  temperature = c(10, 30, 18, "15"),
  northern_hemisphere = c(TRUE, TRUE, FALSE, "FALSE"),
  has_kangaroo = c(FALSE, FALSE, FALSE, 1)
)
str(country_climate)
installed.packages()
instPacks <- installed.packages()
str(instPacks)
## create the matrix
ip <- installed.packages()
head(ip)
## try also View(ip)
## number of package
nrow(ip)
## names of all installed packages
rownames(ip)
## type of information we have about each package
colnames(ip)

set.seed(123)
rmat <- matrix(rnorm(3000), 1000, 3)
dim(rmat)
head(rmat)
library(tidyverse)
library(lubridate)

date1 <- ymd(20240712)
date1
date2 <- ymd("23 Oct 23")
date2
date3 <- dmy("23rd October 1945")

rna <- read_csv("data/rnaseq.csv")
## view the data
rna

rna <- read.csv("data/rnaseq.csv")
rna

library(magrittr)

rna3 <- rna %>% filter(sex == "Female" & time == 0 & expression > 50000) %>% select(gene, sample, time, expression, age)
rna3


rna4 <- rna %>% mutate(expression = log(expression)) %>% 
  filter(chromosome_name %in% c("X", "Y"), 
         !is.na(phenotype_description), 
         expression > 5) %>% 
  select(gene, chromosome_name, phenotype_description, sample, expression)
rna4
