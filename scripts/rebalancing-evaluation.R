library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

dataDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/rebalance-evaluation"

# Dataset to verbose label
dataset2Label <- function(s) {
  datatype <- str_extract(s, "^(dna|aa)")
  dataset <- str_extract(s, "(?<=_)[:alnum:]+(?=_)")
  ranksPerNode <- str_extract(s, "(?<=_)[:digit:]+(?=@)")
  nNodes <- str_extract(s, "(?<=@)[:digit:]+$")
  description <- paste(toupper(datatype), " dataset ", dataset, "\n", nNodes, " nodes with ", ranksPerNode, " ranks each", sep = "")
  return(description)
}

### Data loading ###
# load profiling data from single csv file and add `dataset` column
readFile <- function(flnm) {
  read_csv(
    flnm,
    col_types = cols(
      time = col_integer(),
      rank = col_integer(),
      workMs = col_double()
    )) %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
stats <-
  list.files(
    path = dataDir,  
    pattern = "*.workByRank.csv", 
    full.names = T) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset))

### Plotting ###
datasetLabels <- unique(stats$dataset)
names(datasetLabels) <- datasetLabels
datasetLabels <- map_chr(datasetLabels, dataset2Label)

stats %>%
  #filter(dataset == "aa_rokasA1_10@4") %>%
  group_by(dataset, time) %>%
  summarize(
    meanWorkMs = mean(workMs),
    maxWorkMs = max(workMs)
  ) %>%
  mutate(
    worseThanMeanRel = maxWorkMs / meanWorkMs,
    worseThanMeanAbs = maxWorkMs - meanWorkMs
  ) %>%
ggplot() +
  geom_point(aes(x = time, y = worseThanMeanAbs)) +
  facet_wrap(
    ~dataset,
    labeller = labeller(dataset = datasetLabels),
    scale = "free"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 60, hjust = 1),
    #legend.position = c(0.23, 0.75)
    #legend.title = element_blank()
  )
