library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(lemon)
library(functional)
library(forcats)

### Helpers ###
# Bin factor levels and helper funcitons
binFactorLevels <- map_chr(0:63, function(bit) {
  sprintf("[2^%02d,2^%02d) ns", bit, bit + 1)
})

binFactor2num <- function (factors) map_int(factors, ~which(. == binFactorLevels))
binFactorMinNum <- function (factors) min(binFactor2num(factors))
binFactorMaxNum <- function (factors) max(binFactor2num(factors))
binFactorRangeNum <- function (factors) binFactorMinNum(factors):binFactorMaxNum(factors)
binFactorRange <- function (factors) binFactorLevels[binFactorRangeNum(factors)]

# Dataset to verbose label
dataset2Label <- function(s) {
  datatype <- str_extract(s, "^(dna|aa)")
  dataset <- str_extract(s, "(?<=_)[:alnum:]+(?=_)")
  ranksPerNode <- str_extract(s, "(?<=_)[:digit:]+(?=@)")
  nNodes <- str_extract(s, "(?<=@)[:digit:]+$")
  description <- paste(toupper(datatype), " dataset ", dataset, "\n", nNodes, " nodes with ", ranksPerNode, " ranks each", sep = "")
  return(description)
}

### Setup variables ###
csvDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/"

### Data loading ###
# load profiling data from single csv file and add `dataset` column
readFile <- function(flnm) {
  read_csv(flnm, col_types = "iciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii") %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
data <-
  list.files(
    path = csvDir,  
    pattern = "*.csv", 
    full.names = T) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset))

### Secondary metrices ###
data$eventCount <- data %>% select(ends_with("ns")) %>% rowSums

hist_quantile_bin <- function(row, quantile) {
  runningSum <- 0
  for (bin in colnames(row)) {
    if (endsWith(bin, ") ns")) {
      runningSum <- runningSum + row[[bin]];
      if (runningSum >= row$eventCount * quantile) {
        return(bin);
      }
    }
  }
}

data$medianBin <- by(data, 1:nrow(data), Curry(hist_quantile_bin, quantile = 0.5))
data$medianBin <- factor(data$medianBin, levels = binLevels)
data$q05 <- by(data, 1:nrow(data), Curry(hist_quantile_bin, quantile = 0.05))
data$q95 <- by(data, 1:nrow(data), Curry(hist_quantile_bin, quantile = 0.95))

summary <- data %>%
  group_by(dataset, timer, medianBin) %>%
  summarise(count = n()) %>%
  mutate(fraction = count / sum (count))

### Plotting ###

# Median time spent bin
datasetLabels <- unique(summary$dataset)
names(datasetLabels) <- datasetLabels
datasetLabels <- map_chr(datasetLabels, dataset2Label)

ggplot() +
  geom_bar(
    data = filter(summary, timer == "Work"),
    mapping = aes(x = medianBin, y = fraction, fill = "Work"),
    stat = "identity",
    position = position_nudge(x = 0.2),
    width = 0.4) +
  geom_bar(
    data = filter(summary, timer == "MPI_Allreduce"),
    mapping = aes(x = medianBin, y = fraction, fill = "MPI_Allreduce"),
    stat = "identity",
    position = position_nudge(x = -0.2),
    width = 0.4) +
  scale_x_discrete(limits = binFactorRange(data$medianBin)) +
  facet_rep_wrap(
    ~dataset,
    repeat.tick.labels = TRUE,
    labeller = labeller(dataset = datasetLabels)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "time spent",
    y = "fraction of ranks with their median time spent equal to \"time spent\"",
    fill = "Code Segment")

# Variance in time spent in different code segments
ggplot() +
  geom_linerange(
    data = filter(data, timer == "MPI_Allreduce"),
    mapping = aes(ymin = q05, ymax = q95, x = rank, color = "MPI_Allreduce"),
    position = position_dodge(.2)
  ) +
  # geom_linerange(
  #   data = filter(data, timer == "Work"),
  #   mapping = aes(ymin = q05, ymax = q90, x = rank, color = "Work"),
  #   position = position_dodge(.2)
  # ) +
  #coord_flip() +
  facet_wrap(~dataset, scale = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2") #+
  #scale_y_discrete(limits = binFactorRange(union(data$q05, data$q90)))
