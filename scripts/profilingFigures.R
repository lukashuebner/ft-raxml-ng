library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(lemon)
library(functional)
library(forcats)
library(svglite)

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
csvDirRelative <- "/home/lukas/Documents/Uni/Masterarbeit/raxml-run/relative"

### Data loading ###
# load profiling data from single csv file and add `dataset` column
readFileRelative <- function(flnm) {
  read_csv(flnm, col_types = "icciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii") %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileDataRelative <-
  list.files(
    path = csvDirRelative,  
    pattern = "*.proFile.csv", 
    full.names = T) %>% 
  map_df(~readFileRelative(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset))

# Clean-up
proFileDataRelative <- filter(proFileDataRelative, secondsPassed > 0) %>% # For secondsPassed = 0, there are to few measurements
  mutate(processor = str_extract(processor, "[:alnum:]+(?=.localdomain)")) # Remove .localdomain at end of hostnames

# Conversion to double is needed to avoid integer overflow
proFileDataRelative_timeless <- proFileDataRelative %>%
  group_by(rank, processor, timer, dataset) %>%
  select(-secondsPassed) %>%
  summarise_each(funs(sum(as.numeric(.)))) %>%
  ungroup()

### Secondary metrices ###
proFileDataRelative_timeless$eventCount <- proFileDataRelative_timeless %>%
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

proFileDataRelative_timeless$medianBin <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.5))
proFileDataRelative_timeless$medianBin <- factor(proFileDataRelative_timeless$medianBin, levels = binFactorLevels)
proFileDataRelative_timeless$q00 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.00))
proFileDataRelative_timeless$q05 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.05))
proFileDataRelative_timeless$q25 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.25))
proFileDataRelative_timeless$q75 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.75))
proFileDataRelative_timeless$q95 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 0.95))
proFileDataRelative_timeless$q100 <- by(proFileDataRelative_timeless, 1:nrow(proFileDataRelative_timeless), Curry(hist_quantile_bin, quantile = 1.00))


### Plotting ###

datasetLabels <- unique(proFileDataRelative_timeless$dataset)
names(datasetLabels) <- datasetLabels
datasetLabels <- map_chr(datasetLabels, dataset2Label)

datasets <- unique(proFileDataRelative_timeless$dataset)
maxRanks <- map_int(datasets, function(ds) {
  max(filter(proFileDataRelative_timeless, dataset == ds)$rank)
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
    data = filter(proFileDataRelative_timeless, timer == "MPI_Allreduce"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank, color = as.character(as.integer(rank / 20)))
  ) +
  geom_linerange(
    data = filter(proFileDataRelative_timeless, timer == "Work"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1, color = as.character(as.integer(rank / 20)))
  ) +
  geom_point(
    data = filter(proFileDataRelative_timeless, timer == "MPI_Allreduce"),
    mapping = aes(y = midBin(medianBin), x = rank),
    color = "black",
    size = 0.5
  ) +
  geom_point(
    data = filter(proFileDataRelative_timeless, timer == "Work"),
    mapping = aes(y = midBin(medianBin), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1),
    color = "black",
    size = 0.5
  ) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) +
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

# Boxplot of range of difference to fastest rank
proFileDataRelative_timeless_long <- proFileDataRelative_timeless %>%
  filter(dataset == "dna_rokasD1_20@4", timer == "MPI_Allreduce") %>%
  gather(
    key = "bin",
    value = "count",
    ends_with(" ns")
  )

ggplot(
    data = filter(proFileDataRelative_timeless_long),
    aes(x = rank)
  ) +
  geom_crossbar(aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), y = midBin(medianBin))) +
  geom_text(
    data = filter(proFileDataRelative_timeless_long,
      midBin(bin) > upperBinBorder(q95) | midBin(bin) < lowerBinBorder(q05),
      count != 0
    ),
    aes(
      y = midBin(bin),
      label = signif(count / eventCount, digits = 1)
    ),
    size = 2.5
  ) + 
  facet_wrap(~processor, scale = "free") +
  scale_y_log10(
    breaks = yBreaksGenerator,
    labels = yLabelGenerator
  ) +
  # scale_x_continuous(
  #   breaks = xBreaksGenerator,
  #   labels = xLabelsGenerator,
  #   expand = c(0.02, 0)
  # ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_color_manual(values = rep(twentyColors, 4)) +
  guides(color = FALSE) +
  labs(
    x = "rank",
    y = "range of time difference to fastest rank",
    colour = "Code Segment"
  )

### Behaviour of plots over time
# proFileDataRelative_long <- proFileDataRelative %>%
#   gather(key = bin, value = count, ends_with(" ns"))

proFileDataRelative$eventCount <- proFileDataRelative %>%
  select(ends_with(" ns")) %>%
  rowSums

proFileDataRelative$slowest <- by(proFileDataRelative, 1:nrow(proFileDataRelative), Curry(hist_quantile_bin, quantile = 1))
proFileDataRelative$slowest <- factor(proFileDataRelative$slowest, levels = binFactorLevels)

proFileDataRelative$q95 <- by(proFileDataRelative, 1:nrow(proFileDataRelative), Curry(hist_quantile_bin, quantile = 0.95))
proFileDataRelative$q95 <- factor(proFileDataRelative$q95, levels = binFactorLevels)

undominatedRanks <- proFileDataRelative %>%
  filter(timer == "Work") %>%
  group_by(dataset, secondsPassed) %>%
  arrange() %>%
  summarise(
    max = max(midBin(slowest)),
    rank = which.max(midBin(slowest)))

undominatedRanks <-
  left_join(
    undominatedRanks,
    filter(proFileDataRelative, timer == "Work"),
    by = c("dataset", "secondsPassed", "rank")
  ) %>%
  select(dataset, rank) %>%
  unique

getUndominatedRanks <- function(ds) {
  return(unique(filter(undominatedRanks, dataset == ds)$rank))
}

undominatedRanks_list = map(unique(undominatedRanks$dataset), getUndominatedRanks)
names(undominatedRanks_list) = unique(undominatedRanks$dataset)

proFileDataRelative$color = map2_int(proFileDataRelative$rank, proFileDataRelative$dataset, function(rank, dataset) {
  #return(rank %in% undominatedRanks_list[[dataset]])
  match(rank, undominatedRanks_list[[dataset]]) # Will be NA if not found
})

# Bar chart difference slowest, fastest
ggplot(
    data = filter(proFileDataRelative, !is.na(color))
  ) +
  geom_col(aes(x = as.factor(secondsPassed), y = midBin(slowest), fill = as.factor(color)), position = position_dodge()) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels), ncol = 1) +
  scale_y_log10(
    breaks = yBreaksGenerator,
    labels = yLabelGenerator
  ) +
  # scale_x_continuous(
  #   breaks = xBreaksGenerator,
  #   labels = xLabelsGenerator,
  #   expand = c(0.02, 0)
  # ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = FALSE) +
  labs(
    x = "checkpoint",
    y = "time difference of slowest and fastest rank"
  )

# Slowest rank over time
datasetLabelsOneLine = map_chr(datasetLabels, function(s) { str_replace(s, "\n", " ")})

slowestData <- filter(proFileDataRelative, !is.na(color)) %>%
  group_by(secondsPassed, dataset) %>%
  summarise(slowest = max(midBin(slowest)))
  
ggplot(data = slowestData) +
  geom_point(aes(x = secondsPassed, y = slowest)) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabelsOneLine), ncol = 1) +
  scale_y_continuous(
    labels = yLabelGenerator
  ) +
  # scale_x_continuous(
  #   breaks = xBreaksGenerator,
  #   labels = xLabelsGenerator,
  #   expand = c(0.02, 0)
  # ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = FALSE) +
  labs(
    x = "runtime",
    y = "time difference of slowest and fastest rank"
  )
#ggsave(filename = "/slowest_over_time_small.svg", device = "svg", path = "~", 
#       height = 6.86, width = 12.2, units = "in", dpi = 200)

### Absolute plots
csvDirAbsolute <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/absolute"

# load profiling data from single csv file and add `dataset` column
readFileAbsolute <- function(flnm) {
  read_csv(flnm, col_types = "iciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii") %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileDataAbsolute <-
  list.files(
    path = csvDirAbsolute,  
    pattern = "*.csv", 
    full.names = T) %>% 
  map_df(~readFileAbsolute(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset))

proFileDataAbsolute$eventCount <- proFileDataAbsolute %>%
  select(-rank, -timer, -dataset) %>%
  rowSums

proFileDataAbsolute$medianBin <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.5))
proFileDataAbsolute$medianBin <- factor(proFileDataAbsolute$medianBin, levels = binFactorLevels)
proFileDataAbsolute$q00 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.00))
proFileDataAbsolute$q05 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.05))
proFileDataAbsolute$q25 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.25))
proFileDataAbsolute$q75 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.75))
proFileDataAbsolute$q95 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 0.95))
proFileDataAbsolute$q100 <- by(proFileDataAbsolute, 1:nrow(proFileDataAbsolute), Curry(hist_quantile_bin, quantile = 1.00))


# How many of the difference-to-fastest measurements were above the worst absolute-time measurement?
absoluteWorstTiming <- proFileDataAbsolute %>% select(rank, timer, dataset, q100) %>% rename(absoluteWorst = q100)
relatativeWorstTiming <- proFileDataRelative_timeless %>% select(rank, timer, dataset, q100) %>% rename(relativeWorst = q100)
worstTiming <- inner_join(absoluteWorstTiming, relatativeWorstTiming, by = c("rank", "timer", "dataset")) %>% filter(midBin(relativeWorst) > midBin(absoluteWorst))
left_join(worstTiming, proFileDataRelative_timeless, by = c("rank", "timer", "dataset")) %>%
  select(-q00, -q05, -q25, -q75, -q95, -q100, -medianBin) %>%
  gather(key = "bin", value = "count", ends_with(" ns")) %>%
  filter(midBin(bin) > midBin(absoluteWorst)) %>%
  group_by(timer, dataset, eventCount) %>%
  summarise(worseThanWorstCnt = sum(count)) %>%
  mutate(worseThanWorstFreq = worseThanWorstCnt / eventCount)
  
  
### Plotting ###

datasetLabels <- unique(proFileDataAbsolute$dataset)
names(datasetLabels) <- datasetLabels
datasetLabels <- map_chr(datasetLabels, dataset2Label)

datasets <- unique(proFileDataAbsolute$dataset)
maxRanks <- map_int(datasets, function(ds) {
  max(filter(proFileDataAbsolute, dataset == ds)$rank)
})
names(maxRanks) <- datasets

ggplot() +
  geom_linerange(
    data = filter(proFileDataAbsolute, timer == "MPI_Allreduce"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank, color = as.character(as.integer(rank / 20)))
  ) +
  geom_linerange(
    data = filter(proFileDataAbsolute, timer == "Work"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1, color = as.character(as.integer(rank / 20)))
  ) +
  geom_point(
    data = filter(proFileDataAbsolute, timer == "MPI_Allreduce"),
    mapping = aes(y = midBin(medianBin), x = rank),
    color = "black",
    size = 0.5
  ) +
  geom_point(
    data = filter(proFileDataAbsolute, timer == "Work"),
    mapping = aes(y = midBin(medianBin), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1),
    color = "black",
    size = 0.5
  ) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) +
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
    y = "range of absolute time spent in code segment (0.05 to 0.95 quantiles)",
    colour = "Code Segment"
  )

### Different MPI Implementations ###
csvDirDifferentMPI <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/differentMpiImplementations"

# load profiling data from single csv file and add `dataset` column
readFileDiffMPI <- function(flnm) {
  read_csv(flnm, col_types = "icciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii") %>%
    mutate(dataset = str_extract(basename(flnm), "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}")) %>%
    mutate(mpi = str_extract(basename(flnm), "^([:alnum:]|\\.)+"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileData_diffMPI <-
  list.files(
    path = csvDirDifferentMPI,  
    pattern = "*.proFile.csv", 
    full.names = T) %>% 
  map_df(~readFileDiffMPI(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset))

# Clean-up
proFileData_diffMPI <- filter(proFileData_diffMPI, secondsPassed > 0) %>% # For secondsPassed = 0, there are to few measurements
  mutate(processor = str_extract(processor, "[:alnum:]+(?=.localdomain)")) # Remove .localdomain at end of hostnames

# Conversion to double is needed to avoid integer overflow
proFileData_diffMPI_timeless <- proFileData_diffMPI %>%
  group_by(rank, processor, timer, dataset, mpi) %>%
  select(-secondsPassed) %>%
  summarise_each(funs(sum(as.numeric(.)))) %>%
  ungroup()

### Secondary metrices ###
proFileData_diffMPI_timeless$eventCount <- proFileData_diffMPI_timeless %>%
  select(-rank, -processor, -timer, -dataset, -mpi) %>%
  rowSums

proFileData_diffMPI_timeless$medianBin <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.5))
proFileData_diffMPI_timeless$medianBin <- factor(proFileData_diffMPI_timeless$medianBin, levels = binFactorLevels)
proFileData_diffMPI_timeless$q00 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.00))
proFileData_diffMPI_timeless$q05 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.05))
proFileData_diffMPI_timeless$q25 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.25))
proFileData_diffMPI_timeless$q75 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.75))
proFileData_diffMPI_timeless$q95 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 0.95))
proFileData_diffMPI_timeless$q100 <- by(proFileData_diffMPI_timeless, 1:nrow(proFileData_diffMPI_timeless), Curry(hist_quantile_bin, quantile = 1.0))

ggplot() +
  geom_linerange(
    data = filter(proFileData_diffMPI_timeless, timer == "MPI_Allreduce"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank, color = as.character(as.integer(rank / 20)))
  ) +
  geom_linerange(
    data = filter(proFileData_diffMPI_timeless, timer == "Work"),
    mapping = aes(ymin = lowerBinBorder(q05), ymax = upperBinBorder(q95), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1, color = as.character(as.integer(rank / 20)))
  ) +
  geom_point(
    data = filter(proFileData_diffMPI_timeless, timer == "MPI_Allreduce"),
    mapping = aes(y = midBin(medianBin), x = rank),
    color = "black",
    size = 0.5
  ) +
  geom_point(
    data = filter(proFileData_diffMPI_timeless, timer == "Work"),
    mapping = aes(y = midBin(medianBin), x = rank + maxRanks[as.character(dataset)] * 1.03 + 1),
    color = "black",
    size = 0.5
  ) +
  facet_grid(vars(mpi), vars(dataset), scale = "free_x", labeller = labeller(dataset = datasetLabels)) +
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

### Overall statistics ###
csvDirOverallStats <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/overall"

# load profiling data from single csv file and add `dataset` column
readFileOverallStats<- function(flnm) {
  read_csv(flnm, col_types = cols(
      rank = col_integer(),
      processor = col_character(),
      nsSumInsideMPI = col_double(),
      nsSumOutsideMPI = col_double(),
      nsSumWait = col_double(),
      nsSumWork = col_double(),
      numIterations = col_double(),
      timesIWasSlowest = col_double()
    )) %>%
    mutate(dataset = str_extract(basename(flnm), "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileData_overallStats <-
  list.files(
    path = csvDirOverallStats,  
    pattern = "*.overallStats.csv", 
    full.names = T) %>% 
  map_df(~readFileOverallStats(.)) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset))

proFileData_overallStats <- inner_join(
  proFileData_overallStats,
  proFileData_overallStats %>%
    group_by(dataset) %>%
    summarise(avgWork = mean(nsSumWork)) %>%
    ungroup(),
  by = "dataset"
)

# Compute number of ranks
proFileData_overallStats <- inner_join(
  by = "dataset",
  proFileData_overallStats,
  proFileData_overallStats %>%
    group_by(dataset) %>%
    summarise(nRanks = max(rank) + 1)
)

# Imbalance of work across ranks
ggplot(data = proFileData_overallStats) +
  geom_histogram(aes(nsSumWork / avgWork), stat = "bin", binwidth = 0.01) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) + 
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "work on this rank / average work across all ranks",
    y = "rank count"
  )

# How often was a rank the slowest one?
ggplot(data = proFileData_overallStats) +
  geom_histogram(aes(timesIWasSlowest / numIterations), stat = "bin", binwidth = 0.01) +
  geom_vline(aes(xintercept = 1 / nRanks)) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) + 
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "fraction of iterations a rank took longest to finish work",
    y = "rank count"
  )

# Imbalance of work and communication
ggplot(data = proFileData_overallStats) +
  geom_histogram(aes(nsSumOutsideMPI / (nsSumInsideMPI + nsSumOutsideMPI)), stat = "bin", binwidth = 0.05) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) + 
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "time spent doing work / total runtime",
    y = "rank count"
  )

## Correlation of work/runtime ~ nRanks
proFileData_overallStats_summary <- proFileData_overallStats %>%
  group_by(dataset, nRanks) %>%
  summarise(fractionDoingWork = mean(nsSumOutsideMPI / (nsSumInsideMPI + nsSumOutsideMPI)))

dataDistribution <- tribble(
  ~dataset           , ~partitions, ~sites, ~weightPerThread,
  "dna_rokasD6_20@10", 1         ,  1133    ,  18128,
  "dna_rokasD2a_20@1", 1         , 57134    , 914144,
  "dna_rokasD4_20@8" , 1         ,  1037    ,  16592,
  "dna_rokasD1_20@4" , 1         ,  9331    , 149296,
  "dna_rokasD1_20@20", 1         ,  1867    ,  29872,
  "dna_rokasD1_20@18", 1         ,  2074    ,  33184,
  "dna_ShiD9_20@1"   , 1         ,   666    ,  10656,
  "aa_rokasA4_20@8"  , 1         ,  9675    , 193500
) %>% mutate(dataset = as.factor(dataset))

proFileData_overallStats_summary <- inner_join(
  by = "dataset",
  proFileData_overallStats_summary,
  dataDistribution
)

### Fractional Profiling ###
csvDirFractional <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/fractional"

# load profiling data from single csv file and add `dataset` column
# readFileFractional<- function(flnm) {
#   read_csv(flnm, col_types = cols(
#     rank = col_integer(),
#     processor = col_character(),
#     timer = col_character(),
#     secondsPassed = col_double(),
#     bin = col_character(),
#     count = col_double()
#   )) %>%
#     mutate(dataset = str_extract(basename(flnm), "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}")) %>%
#     group_by(rank, processor, timer, dataset, bin) %>%
#     select(-secondsPassed) %>%
#     summarise(count = sum(count))
# }
# 
# # Load all profiling data from multiple runs and aggregate it into `data`
# proFileData_fractionalStats <-
#   list.files(
#     path = csvDirFractional,  
#     pattern = "*.proFile.csv", 
#     full.names = T) %>% 
#   map_df(~readFileFractional(.)) %>%
#   as_tibble()
# 
# write_csv(proFileData_fractionalStats, file.path(csvDirFractional, "fractional_timeless.csv.bz2"))

proFileData_fractionalStats <- read_csv(
  file.path(csvDirFractional, "fractional_timeless.csv.bz2"),
  col_types = cols(
    rank = col_integer(),
    processor = col_character(),
    timer = col_character(),
    bin = col_character(),
    count = col_double()
))

fractionalMidBin <- function(bin) {
  from <- as.double(str_extract(bin, "^[:digit:]{1,2}\\.[:digit:]{3}"))
  to <- as.double(str_extract(bin, "[:digit:]{1,2}\\.[:digit:]{3}$"))
  
  return((from + to) / 2)
}

proFileData_fractionalStats$midBin <- fractionalMidBin(proFileData_fractionalStats$bin)

freq_quantiles <- function(data, q) {
  data %>%
    group_by(rank, processor, timer, dataset) %>%
    arrange(fractionalMidBin(bin)) %>%
    mutate(quantile = cumsum(count / sum(count))) %>%
    filter(quantile >= q) %>%
    slice(1) %>%
    select(-quantile, -count, -bin)
}

proFileData_fractionalStats_summary <-
  inner_join(
    inner_join(
      inner_join(
        proFileData_fractionalStats %>%
          group_by(rank, processor, timer, dataset) %>%
          summarise(
            min = min(midBin),
            max = max(midBin),
        ),
        freq_quantiles(proFileData_fractionalStats, 0.00) %>% rename(q01 = midBin),
        by = c("rank", "processor", "timer", "dataset")
      ),
      freq_quantiles(proFileData_fractionalStats, 1.00) %>% rename(q99 = midBin),
      by = c("rank", "processor", "timer", "dataset")
    ),
    freq_quantiles(proFileData_fractionalStats, 0.5) %>% rename(median = midBin),
    by = c("rank", "processor", "timer", "dataset")
  )

ggplot() +
  geom_linerange(
    data = filter(proFileData_fractionalStats_summary, timer == "MPI_Allreduce"),
    mapping = aes(ymin = q01, ymax = q99, x = rank, color = as.character(as.integer(rank / 20)))
  ) +
  geom_linerange(
    data = filter(proFileData_fractionalStats_summary, timer == "Work"),
    mapping = aes(ymin = q01, ymax = q99, x = rank + maxRanks[as.character(dataset)] * 1.03 + 1, color = as.character(as.integer(rank / 20)))
  ) +
  geom_point(
    data = filter(proFileData_fractionalStats_summary, timer == "MPI_Allreduce"),
    mapping = aes(y = median, x = rank),
    color = "black",
    size = 0.5
  ) +
  geom_point(
    data = filter(proFileData_fractionalStats_summary, timer == "Work"),
    mapping = aes(y = median, x = rank + maxRanks[as.character(dataset)] * 1.03 + 1),
    color = "black",
    size = 0.5
  ) +
  facet_wrap(~dataset, scale = "free", labeller = labeller(dataset = datasetLabels)) +
  scale_y_log10(breaks = scales::pretty_breaks(n = 6)) +
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
    y = "(time spent on work/comm. package) / (avg. time across all ranks)",
    colour = "Code Segment"
  )

### Are site-repeats the cause of imbalance? ###
csvDirSiteRepeats <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/site-repeats"

# load profiling data from single csv file and add `dataset` column
readFileSiteRepeats<- function(flnm) {
  read_csv(flnm, col_types = cols(
    rank = col_integer(),
    processor = col_character(),
    nsSumInsideMPI = col_double(),
    nsSumOutsideMPI = col_double(),
    nsSumWait = col_double(),
    nsSumWork = col_double(),
    numIterations = col_double(),
    timesIWasSlowest = col_double()
  )) %>%
    mutate(dataset = str_extract(basename(flnm), "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}")) %>%
    mutate(sr = str_extract(basename(flnm), "^(no-sr|sr)"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
proFileData_siteRepeats <-
  list.files(
    path = csvDirSiteRepeats,  
    pattern = "*.overallStats.csv", 
    full.names = T) %>% 
  map_df(~readFileSiteRepeats(.)) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset))

proFileData_siteRepeats <- inner_join(
  proFileData_siteRepeats,
  proFileData_siteRepeats %>%
    group_by(dataset, sr) %>%
    summarise(avgWork = mean(nsSumWork)) %>%
    ungroup(),
  by = c("dataset", "sr")
)

# Compute number of ranks
proFileData_siteRepeats <- inner_join(
  by = c("dataset", "sr"),
  proFileData_siteRepeats,
  proFileData_siteRepeats %>%
    group_by(dataset, sr) %>%
    summarise(nRanks = max(rank) + 1)
)

# Imbalance of work across ranks
ggplot(data = proFileData_siteRepeats) +
  geom_histogram(aes(nsSumWork / avgWork), stat = "bin", binwidth = 0.01) +
  facet_grid(vars(sr), vars(dataset), scale = "free_y", labeller = labeller(dataset = datasetLabels)) + 
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "work on this rank / average work across all ranks",
    y = "rank count"
  )

# Imbalance of work and communication
ggplot(data = proFileData_siteRepeats) +
  geom_histogram(aes(nsSumOutsideMPI / (nsSumInsideMPI + nsSumOutsideMPI)), stat = "bin", binwidth = 0.05) +
  facet_grid(vars(sr), vars(dataset), scale = "free_y", labeller = labeller(dataset = datasetLabels)) + 
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "time spent doing work / total runtime",
    y = "rank count"
  )

ggplot(data = proFileData_siteRepeats) %>%
  geom_bar(aes(x = nsSumInsideMPI + nsSumOutsideMPI))