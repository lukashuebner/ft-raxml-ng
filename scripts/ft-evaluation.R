library(ggplot)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

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
  read_csv(
    flnm,
    col_types = 
      cols(
        .default = col_double(),
        processor = col_character()
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

#overallStats$msPerMiniCheckpoint <- overallStats$nsSumMiniCheckpoints / overallStats$numMiniCheckpoints / 10^6
nsSum_long <- overallStats %>% select(- starts_with("num")) %>% gather(key = timer, value = nsSum, starts_with("nsSum"))
num_long <- overallStats %>% select(- starts_with("nsSum")) %>% gather(key = timer, value = num, starts_with("num"))
nsSum_long$timer = str_replace(nsSum_long$timer, "nsSum", "")
num_long$timer = str_replace(num_long$timer, "num", "")

stats_long <- inner_join(by = c("rank", "processor", "dataset", "timer"), num_long, nsSum_long)
stats_long$nsSum <- stats_long$nsSum / stats_long$num

stats_long <- stats_long %>% group_by(dataset, timer) %>%
  summarise(
    nsAvg = mean(nsSum),
    nsSD = sd(nsSum),
  ) %>%
  mutate(
    msAvg = nsAvg / 10^6,
    msSD = nsSD / 10^6
  )

stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "dna_ShiD9_20@1") %>% arrange(desc(msAvg)) %>% print(n=30)
stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "aa_rokasA1_10@4") %>% arrange(desc(msAvg)) %>% print(n=30)
stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "dna_rokasD1_20@20") %>% arrange(desc(msAvg)) %>% print(n=30)
