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

rna5 <- rna %>% filter(gene == "Dok3") %>% 
  group_by(time) %>% 
  summarise(mean_expression = mean(expression))
rna5

rna %>% count(sample)
rna

rna %>% group_by(sample) %>% summarise(n = sum(expression)) %>% arrange(desc(n))

rna %>% filter(sample == "GSM2545336") %>% count(gene_biotype) %>% arrange(desc(n))

rna %>% filter(phenotype_description == "abnormal DNA methylation") %>% 
  mutate(expression = log(expression)) %>% 
  group_by(gene, time) %>% summarise(mean = mean(expression)) %>% 
  arrange(gene, time)

(rna_wider <- rna %>% select(gene, mouse, expression) %>% 
  pivot_wider(names_from = mouse, values_from = expression))

(rna_longer <- rna_wider %>% 
    pivot_longer(names_to = "mouse", values_to = "expression", -gene))

(q2 <- rna %>% filter(chromosome_name %in% c("X", "Y")) %>% 
    select(sex, chromosome_name, expression) %>% 
    group_by(sex, chromosome_name) %>% 
    summarise(mean_expression = mean(expression)) %>% 
    pivot_wider(names_from = sex, values_from = mean_expression))

(q3 <- rna %>% select(gene, expression, time) %>% 
    group_by(gene, time) %>% 
    summarise(mean_expression = mean(expression)) %>% 
    pivot_wider(names_from = time, values_from = mean_expression))

(q4 <- q3 %>% mutate(fold_0to8 = `8`/`0`, fold_8to4 = `8`/`4`) %>% 
    pivot_longer(names_to = "comparison", values_to = "foldChanges", fold_0to8:fold_8to4))

(q4a <- q4 %>% select(gene, comparison, foldChanges))
q4a %>% pull(foldChanges)

library(datasauRus)
install.packages("datasauRus")
if(requireNamespace("ggplot2")){
  library(ggplot2)
  ggplot(datasaurus_dozen, aes(x = x, y = y, colour = dataset))+
    geom_point()+
    theme_void()+
    theme(legend.position = "none")+
    facet_wrap(~dataset, ncol = 3)
}

download.file(url = "https://carpentries-incubator.github.io/bioc-intro/data/annot1.csv",
              destfile = "data/annot1.csv")
annot1 <- read_csv(file = "data/annot1.csv")
annot1

rna_mini <- rna %>%
  select(gene, sample, expression) %>%
  head(10)
rna_mini

full_join(rna_mini, annot1)

download.file(url = "https://carpentries-incubator.github.io/bioc-intro/data/annot2.csv",
              destfile = "data/annot2.csv")
annot2 <- read_csv(file = "data/annot2.csv")
annot2

download.file(url = "https://carpentries-incubator.github.io/bioc-intro/data/annot3.csv",
              destfile = "data/annot3.csv")
annot3 <- read_csv(file = "data/annot3.csv")
annot3

full_join(rna_mini, annot3)
write_csv(q3, file = "data_output/rna_wide.csv")

ggplot(data = rna)
library(ggplot2
        )
ggplot(data = rna, mapping = aes(x = expression))
ggplot(data = rna, mapping = aes(x = expression)) +
  geom_histogram(
                 )
rna <- rna %>%
  mutate(expression_log = log2(expression + 1))
ggplot(rna, aes(x = expression_log)) + geom_histogram()

ggplot(data = rna, mapping = aes(x = expression +1)) +
  geom_histogram(
                 ) + scale_x_log10()

rna_fc <- rna %>% select(gene, time,
                         gene_biotype, expression_log) %>%
  group_by(gene, time, gene_biotype) %>%
  summarize(mean_exp = mean(expression_log)) %>%
  pivot_wider(names_from = time,
              values_from = mean_exp) %>%
  mutate(time_8_vs_0 = `8` - `0`, time_4_vs_0 = `4` - `0`)

ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point()
ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point(alpha = 0.3) +
  geom_abline(0, 2)

ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point(alpha = 0.3, aes(color = gene_biotype)) 

?abline
library(hexbin)
ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_hex() 

ggplot(data = rna, mapping = aes(x = expression_log, y = sample)) +
  geom_point(alpha = 0.3, aes(color = time)) 

ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_boxplot()

ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_jitter(alpha = 0.2, color = "tomato") +
  geom_boxplot(alpha = 0)

ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_jitter(alpha = 0.2, color = "tomato") +
  geom_boxplot(alpha = 0) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))

ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample, color = time)) +
  geom_jitter(alpha = 0.2) +
  geom_boxplot(alpha = 0) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))

# time as factor
ggplot(data = rna,
       mapping = aes(y = expression_log,
                     x = sample)) +
  geom_jitter(alpha = 0.2, aes(color = as.factor(time))) +
  geom_boxplot(alpha = 0) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))

ggplot(data = rna,
       mapping = aes(y = expression_log,
                     x = sample)) +
  geom_violin(aes(fill = as.factor(time))) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))

ggplot(data = rna,
       mapping = aes(y = expression_log,
                     x = sample)) +
  geom_violin(aes(fill = sex)) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))

rna_fc <- rna_fc %>% arrange(desc(time_8_vs_0))

genes_selected <- rna_fc$gene[1:10]

sub_rna <- rna %>%
  filter(gene %in% genes_selected)

mean_exp_by_time <- sub_rna %>%
  group_by(gene,time) %>%
  summarize(mean_exp = mean(expression_log))

ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp, group = gene)) +
  geom_line()
ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp, color = gene)) +
  geom_line()

ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp)) + geom_line() +
  facet_wrap(~ gene)

ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y")

mean_exp_by_time_sex <- sub_rna %>%
  group_by(gene, time, sex) %>%
  summarize(mean_exp = mean(expression_log))
ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y")

rna6 <- rna %>% select(expression_log, chromosome_name, time) %>% 
  group_by(chromosome_name, time) %>% 
  summarise(mean_expression = mean(expression_log))
rna6
ggplot(data = rna6, mapping = aes(x = time, y = mean_expression)) +
  geom_line() +
  facet_wrap(~ chromosome_name, scales = "free_y")

# One column, facet by rows
ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = gene)) +
  geom_line() +
  facet_grid(sex ~ .)

ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression")
library(patchwork)
rna$chromosome_name <- factor(rna$chromosome_name,
                              levels = c(1:19,"X","Y"))

count_gene_chromosome <- rna %>% select(chromosome_name, gene) %>%
  distinct() %>% ggplot() +
  geom_bar(aes(x = chromosome_name), fill = "seagreen",
           position = "dodge", stat = "count") +
  labs(y = "log10(n genes)", x = "chromosome") +
  scale_y_log10()

count_gene_chromosome

exp_boxplot_sex <- ggplot(rna, aes(y=expression_log, x = as.factor(time),
                                   color=sex)) +
  geom_boxplot(alpha = 0) +
  labs(y = "Mean gene exp",
       x = "time") + theme(legend.position = "none")

exp_boxplot_sex

library("patchwork")
count_gene_chromosome + exp_boxplot_sex
count_gene_chromosome / exp_boxplot_sex
count_gene_chromosome | exp_boxplot_sex
library("gridExtra")

my_plot <- ggplot(data = mean_exp_by_time_sex,
                  mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression") +
  guides(color=guide_legend(title="Gender")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "royalblue4", size = 12),
        axis.text.y = element_text(colour = "royalblue4", size = 12),
        text = element_text(size = 16),
        panel.grid = element_line(colour="lightsteelblue1"),
        legend.position = "top")
ggsave("fig_output/mean_exp_by_time_sex.png", my_plot, width = 15,
       height = 10)

count_matrix <- read.csv("data/count_matrix.csv",
                         row.names = 1) %>%
  as.matrix()

count_matrix[1:5, ]

sample_metadata <- read.csv("data/sample_metadata.csv")
sample_metadata

gene_metadata <- read.csv("data/gene_metadata.csv")
gene_metadata[1:10, 1:4]

## BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

stopifnot(rownames(count_matrix) == gene_metadata$gene)
stopifnot(colnames(count_matrix) == sample_metadata$sample)
se <- SummarizedExperiment(assays = list(counts = count_matrix),
                           colData = sample_metadata,
                           rowData = gene_metadata)
se

saveRDS(se, file = "data_output/se.rds")
rm(se)
se <- readRDS("data_output/se.rds")
head(se)

head(assay(se))
colData(se)
head(rowData(se))

assay(se)[1:3, colData(se)$time != 4]
# Equivalent to
assay(se)[1:3, colData(se)$time == 0 | colData(se)$time == 8]
# Equivalent to
rna |>
  filter(gene %in% c("Asl", "Apod", "Cyd2d22")) |>
  filter(time != 4) |> select(expression)

colData(se)$center <- rep("University of Illinois", nrow(colData(se)))
colData(se)

#BiocManager::install("tidySummarizedExperiment")
library("tidySummarizedExperiment")

se
