library(readr)
library(ggplot2)
library(grid)


datasetLevels =  c("dna_ShiD9_20@1", "dna_rokasD4_20@8", "dna_rokasD6_20@10", "dna_rokasD1_20@20",
                   "dna_rokasD1_20@18", "dna_rokasD1_20@4", "aa_rokasA4_20@8", "dna_rokasD2a_20@1")

### Data ;-) ###
dataDistribution <- tribble(
  ~dataset           , ~partitions, ~sites, ~weightPerThread,
  "dna_rokasD6_20@10", 1         ,  1133    ,  18128,
  "dna_rokasD2a_20@1", 1         , 57134    , 914144,
  "dna_rokasD4_20@8" , 1         ,  1037    ,  16592,
  "dna_rokasD1_20@4" , 1         ,  9331    , 149296,
  "dna_rokasD1_20@30", 1         ,  1867    ,  29872,
  "dna_rokasD1_20@18", 1         ,  2074    ,  33184,
  "dna_ShiD9_20@1"   , 1         ,   666    ,  10656,
  "aa_rokasA4_20@8"  , 1         ,  9675    , 193500
)
dataDistribution <- dataDistribution %>%
  mutate(dataset = factor(dataset2Label(dataset), levels = dataset2Label(datasetLevels)))

### Setup variables ###
csvDir <- "/home/lukas/Documents/Uni/Masterarbeit/raxml-run/collection3"

### Data loading ###
# load data from single csv file and add `dataset` column
readFile <- function(flnm) {
  read_csv(flnm) %>%
    mutate(dataset = str_extract(flnm, "(dna|aa)_.+_[:digit:]{1,2}@[:digit:]{1,2}"))
}

# Load all data from multiple runs and aggregate it into `data`
callsPerSecondData <-
  list.files(
    path = csvDir,  
    pattern = "*.callsPerSecond.csv", 
    full.names = T) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(timer = factor(timer)) %>%
  mutate(dataset = factor(dataset, levels = datasetLevels))

callsPerSecondData$dataset <-
  factor(
    map_chr(callsPerSecondData$dataset, dataset2Label),
    levels = dataset2Label(datasetLevels))

totalRuntime <- callsPerSecondData %>%
  group_by(timer, dataset) %>%
  summarise(runtime = max(secondsPassed))

callsPerSecondData <- inner_join(callsPerSecondData, totalRuntime, by = c("timer", "dataset"))

callsPerSecondData <- callsPerSecondData %>%
  mutate(fractionOfRuntime = secondsPassed / runtime)

callsPerSecondData %>%
  filter(secondsPassed > 1) %>%
  ggplot(aes(
    x = fractionOfRuntime,
    y = callsPerSecond,
    color = dataset
  )) +
  geom_point(size=0.7) +
  geom_line() + 
  theme_bw() +
  theme(legend.key.height = unit(1.5, "line")) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "fraction of runtime for this task",
    y = "number of MPI_Allreduce calls per second",
    colour = "dataset"
  )

# How well does weight per core correlate with number of MPI_Allreduce calls per second?
callsAndWeight <- callsPerSecondData %>%
  filter(fractionOfRuntime > 0.25, timer == "MPI_Allreduce") %>%
  select(-timer) %>%
  group_by(dataset) %>%
  summarise(averageCallsPerSecond = mean(callsPerSecond))

callsAndWeight <- inner_join(callsAndWeight, dataDistribution, by = c("dataset"))

ggplot(
  data = callsAndWeight,
  aes(
    x = weightPerThread,
    y = 1 / averageCallsPerSecond
  )) +
  geom_point() +
  geom_smooth(method = 'lm', se=FALSE) +
  theme_bw()

cor.test(1 / callsAndWeight$averageCallsPerSecond, callsAndWeight$weightPerThread, method = "pearson" )
