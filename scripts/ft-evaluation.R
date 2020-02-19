library(ggplot)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

dataDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/ft-evaluation"

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
  read_csv(flnm, col_types = cols(
    rank = col_integer(),
    processor = col_character(),
    numMiniCheckpoints = col_integer(),
    numRecoveries = col_integer(),
    # uin64_t might be to large for col_integer()
    nsSumMiniCheckpoints = col_double(),
    nsSumLostWork = col_double(),
    nsSumRecalculateAssignment = col_double(),
    nsSumReloadSites = col_double(),
    nsSumRestoreModels = col_double()
    )) %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all profiling data from multiple runs and aggregate it into `data`
overallStats <-
  list.files(
    path = dataDir,  
    pattern = "*.overallStats.csv", 
    full.names = T) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset))

runtimeStats <- read_csv(
  paste(dataDir, "ft-eval-runtimes.csv", sep = "/"),
  col_types = cols(
    dataset = col_character(),
    runtime = col_double()
  )) %>%
  mutate(dataset = factor(dataset)) %>%
  rename(sRuntime = runtime)

overallStats$msPerMiniCheckpoint <- overallStats$nsSumMiniCheckpoints / overallStats$numMiniCheckpoints / 10^6

overallStatsSummary <- overallStats %>%
  group_by(dataset) %>%
  summarise(
    avg_msPerMiniCheckpoint = mean(msPerMiniCheckpoint),
    stddev_msPerMiniCheckpoint = sd(msPerMiniCheckpoint),
    numMiniCheckpoints = min(numMiniCheckpoints),
    sSumMiniCheckpoints = min(nsSumMiniCheckpoints) / 10^9)

overallStatsSummary <- inner_join(by = "dataset", overallStatsSummary, runtimeStats)
overallStatsSummary$fractionOfRuntime = overallStatsSummary$sSumMiniCheckpoints / overallStatsSummary$sRuntime
