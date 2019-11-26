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
binFactorLevels <- c("0 ns", map_chr(0:63, function(bit) {
  sprintf("[2^%02d,2^%02d) ns", bit, bit + 1)
}))

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
csvDir <- "/home/lukas/Documents/Uni/Masterarbeit/raxml-run/collection3"

### Data loading ###
# load profiling data from single csv file and add `dataset` column
readFile <- function(flnm) {
  read_csv(flnm, col_types = "icciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii") %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileData <-
  list.files(
    path = csvDir,  
    pattern = "*.proFile.csv", 
    full.names = T) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset))

# Clean-up
proFileData <- filter(proFileData, secondsPassed > 0) %>% # For secondsPassed = 0, there are to few measurements
  mutate(processor = str_extract(processor, "[:alnum:]+(?=.localdomain)")) # Remove .localdomain at end of hostnames

# Conversion to double is needed to avoid integer overflow
proFileData_timeless <- proFileData %>%
  group_by(rank, processor, timer, dataset) %>%
  select(-secondsPassed) %>%
  summarise_each(funs(sum(as.numeric(.)))) %>%
  ungroup()

### Secondary metrices ###
proFileData_timeless$eventCount <- proFileData_timeless %>%
  select(-rank, -processor, -timer, -dataset) %>%
  rowSums

hist_quantile_bin <- function(row, quantile) {
  runningSum <- 0
  for (bin in colnames(row)) {
    if (endsWith(bin, " ns")) {
      runningSum <- runningSum + row[[bin]];
      if (runningSum >= row$eventCount * quantile) {
        return(bin);
      }
    }
  }
}

proFileData_timeless$medianBin <- by(proFileData_timeless, 1:nrow(proFileData_timeless), Curry(hist_quantile_bin, quantile = 0.5))
proFileData_timeless$medianBin <- factor(proFileData_timeless$medianBin, levels = binFactorLevels)
proFileData_timeless$q05 <- by(proFileData_timeless, 1:nrow(proFileData_timeless), Curry(hist_quantile_bin, quantile = 0.05))
proFileData_timeless$q95 <- by(proFileData_timeless, 1:nrow(proFileData_timeless), Curry(hist_quantile_bin, quantile = 0.95))

### Plotting ###

datasetLabels <- unique(proFileData_timeless$dataset)
names(datasetLabels) <- datasetLabels
datasetLabels <- map_chr(datasetLabels, dataset2Label)

datasets <- unique(proFileData_timeless$dataset)
maxRanks <- map_int(datasets, function(ds) {
  max(filter(proFileData_timeless, dataset == ds)$rank)
})
names(maxRanks) <- datasets

lowerBinBorder <- function(bins) {
  sapply(bins, function(bin) {
    if (bin == "0 ns") {
      return(1)
    }
    pow <- as.numeric(str_extract_all(bin, "(?<=\\^)[:digit:]{1,2}")[[1]][1])
    return(2**pow)
})}

upperBinBorder <- function(bins) {
  sapply(bins, function(bin) {
    if (bin == "0 ns") {
      return(1)
    }
    pow <- as.numeric(str_extract_all(bin, "(?<=\\^)[:digit:]{1,2}")[[1]][2])
    return(2**pow - 1)
  })}

midBin <- function(bins) {
  sapply(bins, function(bin) {
    if (bin == "0 ns") {
      return(1)
    }
    powL <- as.numeric(str_extract_all(bin, "(?<=\\^)[:digit:]{1,2}")[[1]][1])
    return(1.5 * 2**powL)
})}

ns2us <- function(timeInNs) {
  return(timeInNs / 1000)
}

yBreaksGenerator <- function(limits) {
  lower <- limits[1]
  upper <- limits[2]
  
  breaks <- sapply(seq(from = 0, to = log10(upper)), function(x) 10**x)
  return(breaks)
}

yLabelGenerator <- function(breaks) {
  sapply(breaks, function(time) {
    if (is.na(time)) {
      return("")
    }
    unit <- "ns"
    if (time >= 1000) {
      time = time / 1000
      unit <- "µs"
    }
    if (time >= 1000) {
      time = time / 1000
      unit <- "ms"
    }
    if (time >= 1000) {
      time = time / 1000
      unit <- " s"
    }
    return(paste(time, unit))
})}

xBreaksGenerator <- function(limits) {
  lower <- limits[1]
  upper <- limits[2]
  range <- upper - lower
  withoutPadding <- range / 1.02
  interSegementSpace <- withoutPadding - withoutPadding / 1.03
  dataRange <- withoutPadding - interSegementSpace
  rankRange <- dataRange / 2
  midOfLeft <- rankRange / 2
  midOfRight <- rankRange + interSegementSpace + rankRange / 2
  
  breaks <- c(midOfLeft, midOfRight)
  return(breaks)
}

xLabelsGenerator <- function(breaks) {
  return(c("MPI_Allreduce", "Work"))
}

twentyColors = c("#023fa5", "#bec1d4", "#7d87b9", "#d6bcc0", "#8e063b", "#bb7784", "#4a6fe3", "#b5bbe3", "#8595e1", "#e07b91", "#e6afb9", 
                 "#11c638", "#d33f6a", "#8dd593", "#ead3c6", "#c6dec7", "#f0b98d", "#0fcfc0", "#ef9708", "#9cded6", "#f3e1eb", "#d5eae7",
                 "#f6c4e1", "#f79cd4")

ggplot() +
  geom_linerange(
    data = filter(proFileData_timeless, timer == "MPI_Allreduce"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank, color = as.character(as.integer(rank / 20)))
  ) +
  geom_linerange(
    data = filter(proFileData_timeless, timer == "Work"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1, color = as.character(as.integer(rank / 20)))
  ) +
  geom_point(
    data = filter(proFileData_timeless, timer == "MPI_Allreduce"),
    mapping = aes(y = midBin(medianBin), x = rank),
    color = "black",
    size = 0.5
  ) +
  geom_point(
    data = filter(proFileData_timeless, timer == "Work"),
    mapping = aes(y = midBin(medianBin), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1),
    color = "black",
    size = 0.5
  ) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) +
  #scale_y_discrete(
  #  limits = binFactorRange(union(proFileData_timeless$q95, proFileData_timeless$q05))
  #) +
  scale_y_log10(
    breaks = yBreaksGenerator,
    labels = yLabelGenerator
  ) +
  scale_x_continuous(
    breaks = xBreaksGenerator,
    labels = xLabelsGenerator,
    expand = c(0.02, 0)
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_color_manual(values = twentyColors) +
  guides(color = FALSE) +
  labs(
    x = "rank",
    y = "0.05 to 0.95 quantiles of time difference to fastest rank",
    colour = "Code Segment"
  )
  