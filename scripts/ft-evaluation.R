library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

dataDirMiniCheckpointing <- "/home/lukas/Promotion/profiling/ft-evaluation/mini-checkpointing"
dataDirRecovery <- "/home/lukas/Promotion/profiling/ft-evaluation/recovery-overhead"

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
perRankStats <-
  list.files(
    path = dataDirMiniCheckpointing,  # ADJUST THIS!!
    pattern = "*.overallStats.csv", 
    full.names = TRUE
  ) %>% 
  map_df(~readFile(.)) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset))

runtimeStats <- read_csv(
  paste(dataDirMiniCheckpointing, "ft-eval-runtimes.csv", sep = "/"),
  col_types = cols(
    dataset = col_character(),
    msRuntime = col_double()
  )) %>%
  mutate(
    dataset  = factor(dataset),
    sRuntime = msRuntime / 1000
  )

# Aggregate the per rank statistics into per measurement statistics
## First, convert the measurements into a long format
perRankStats_long_time <- perRankStats %>%
  select(- starts_with("num")) %>%
  gather(
    key = timer,
    value = nsSum,
    starts_with("nsSum")
  ) %>%
  mutate (
    timer = str_replace(timer, "nsSum", "")
  )

perRankStats_long_num <- perRankStats %>%
  select(- starts_with("nsSum")) %>%
  gather(
    key = timer,
    value = num,
    starts_with("num")
  ) %>%
  mutate(
    timer = str_replace(timer, "num", "")
  )

perRankStats_long <- inner_join(
  by = c("rank", "processor", "dataset", "timer"),
  perRankStats_long_time,
  perRankStats_long_num
) %>%
  rename(
    nsRankwiseTotalTime = nsSum,
    nCallsOnThisRank = num
  ) %>%
  mutate(
    nsRankwiseTimePerCall = nsRankwiseTotalTime / nCallsOnThisRank
  )

# How many tree and models updates are there per dataset?
perRankStats_long %>%
  filter(
    rank == 0, # There should be the same number of mini-checkpoints per rank
    timer %in% c("UpdateTree", "UpdateModels")
  )

## Next, aggregate per rank measurements into per dataset summary statistics
acrossRanksStats_long <- perRankStats_long %>%
  group_by(dataset, timer) %>%
  summarise(
    .groups = "keep",
    nsAcrossRanksTimePerCallAvg = mean(nsRankwiseTimePerCall),
    nsAcrossRanksTimePerCallSD = sd(nsRankwiseTimePerCall),
    nsAcrossRanksTotalTimeAvg = mean(nsRankwiseTotalTime),
    nsAcrossRanksTotalTimeSD = sd(nsRankwiseTotalTime)
  ) %>%
  mutate(
    msAcrossRanksTimePerCallAvg = nsAcrossRanksTimePerCallAvg / 10^6,
    msAcrossRanksTimePerCallSD = nsAcrossRanksTimePerCallSD / 10^6,
    msAcrossRanksTotalTimeAvg = nsAcrossRanksTotalTimeAvg / 10^6,
    msAcrossRanksTotalTimeSD = nsAcrossRanksTotalTimeSD / 10^6,
  )

# Plot comparing the time a rank requires to load a part of the MSA {the first time, if someone else already loaded it}
ggplot() +
  geom_errorbar(
    data = filter(acrossRanksStats_long, timer == "LoadAssignmentData"),
    aes(ymin = msAvg - msSD, ymax = msAvg + msSD, x = dataset, color = "rebalancing")
  ) +
  geom_point(
    data = filter(acrossRanksStats_long, timer == "LoadAssignmentData"),
    aes(y = msAvg, x = dataset, color = "rebalancing")
  ) +
  geom_errorbar(
    data = filter(acrossRanksStats_long, timer == "LoadAssignmentDataFirstLoad"),
    aes(ymin = msAvg - msSD, ymax = msAvg + msSD, x = dataset, color = "initial load op.")
  ) +
  geom_point(
    data = filter(acrossRanksStats_long, timer == "LoadAssignmentDataFirstLoad"),
    aes(y = msAvg, x = dataset, color = "initial load op.")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.15, 0.85),
    legend.title = element_blank()
  ) +
  ylab("time [ms]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_10@4" = "AA rokasA1 @ 40 ranks\n172k sites · 60 taxa",
      "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n1806k sites · 58 taxa",
      "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n505k sites · 95 taxa",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n3011k sites · 174 taxa",
      "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n1339k sites · 37 taxa",
      "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n1240k sites · 144 taxa",
      "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n240k sites · 46 taxa",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n20k sites · 815 taxa"
    ),
    limits = (acrossRanksStats_long %>% filter(timer == "LoadAssignmentDataFirstLoad") %>% arrange(msAvg))$dataset
  )

# Quick plot showing the ratio between time spent updating the models and updating the tree
ggplot(acrossRanksStats_long) +
  geom_bar(aes(x = dataset, fill = timer, y = msAcrossRanksTotalTimeAvg), stat="identity")

# Plot the time required when restoring
bind_rows(
  acrossRanksStats_long %>%
    filter(timer %in% c(
      "BalanceLoad", "ComputePartMasters", "CreateTree", "DestroyPartitions",
      "RestoreModels", "TreeinfoDestroy", "TreeinfoResetPartitions", "UTreeGraphDestroy")
    ) %>%
    group_by(dataset) %>%
    summarize(msAcrossRanksTimePerCallAvg = sum(msAcrossRanksTimePerCallAvg), .groups = "keep") %>%
    mutate(timer = "Other"),
  acrossRanksStats_long %>%
    filter(timer %in% c("LoadAssignmentData", "RedoPartitionAssignment", "TreeinfoUpdatePartialsAndCLVs", "TreeinfoInit")) %>%
    select(dataset, timer, msAcrossRanksTimePerCallAvg)
) %>%
ggplot() +
  geom_bar(aes(x = dataset, y = msAcrossRanksTimePerCallAvg, fill = timer), stat = "identity") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.23, 0.65)
  ) +
  ylab("time [ms]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_10@4" = "AA NagyA\n 40 ranks\n172k sites · 60 taxa",
      "aa_rokasA4_20@8" = "AA ChenA4\n 160 ranks\n1806k sites · 58 taxa",
      "aa_rokasA8_20@4" = "AA YangA8\n 80 ranks\n505k sites · 95 taxa",
      "dna_PeteD8_13@20" = "DNA PeteD8\n 260 ranks\n3011k sites · 174 taxa",
      "dna_rokasD1_20@20" = "DNA SongD1\n 400 ranks\n1339k sites · 37 taxa",
      "dna_rokasD2a_20@1" = "DNA MisoD2a\n 20 ranks\n1240k sites · 144 taxa",
      "dna_rokasD4_20@8" = "DNA XiD4\n 160 ranks\n240k sites · 46 taxa",
      "dna_ShiD9_20@1" = "DNA ShiD9\n 20 ranks\n20k sites · 815 taxa"
    ),
    limits = acrossRanksStats_long %>%
        filter(timer %in% c("LoadAssignmentData", "RedoPartitionAssignment", "TreeinfoUpdatePartialsAndCLVs", "TreeinfoInit", "Other")) %>%
        group_by(dataset) %>%
        summarise(msRecovery = sum(msAcrossRanksTimePerCallAvg), .groups = "keep") %>%
        arrange(msRecovery) %>%
        pull(dataset)
  ) + 
  scale_fill_brewer(type="qual")
ggsave("~/Promotion/figures/paper/time-for-recovery.pdf", width = 15, height = 10, units = "cm")

# Plot time required for mini-checkpointing
acrossRanksStats_miniCheckpoints_long <- acrossRanksStats_long %>%
  mutate(
    # Prepare to compute the combined standard deviation via the sum of the variances
    msAcrossRanksTotalTimeVar = msAcrossRanksTotalTimeSD ** 2
  ) %>%
  select(dataset, timer, msAcrossRanksTotalTimeAvg, msAcrossRanksTotalTimeVar) %>%
  filter(timer %in% c("UpdateModels", "UpdateTree")) %>%
  pivot_wider(
    names_from = timer,
    values_from = c(msAcrossRanksTotalTimeAvg, msAcrossRanksTotalTimeVar)
  ) %>%
  mutate(
    timer = "MiniCheckpoint",
    msAcrossRanksTotalTimeAvg = (msAcrossRanksTotalTimeAvg_UpdateModels + msAcrossRanksTotalTimeAvg_UpdateTree) / 2.,
    msAcrossRanksTotalTimeSD = sqrt((msAcrossRanksTotalTimeVar_UpdateModels + msAcrossRanksTotalTimeVar_UpdateTree) / 2.)
  ) %>%
  select(dataset, timer, msAcrossRanksTotalTimeAvg, msAcrossRanksTotalTimeSD)

acrossRanksStats_miniCheckpoints_long <- inner_join(
  acrossRanksStats_miniCheckpoints_long,
  runtimeStats,
  by = c("dataset")
) %>%
  mutate(
    percAcrossRanksTotalTimeAvg = msAcrossRanksTotalTimeAvg / msRuntime * 100,
    percAcrossRanksTotalTimeSD = msAcrossRanksTotalTimeSD / msRuntime * 100
  )

### Oveall time for mini-checkpoints
acrossRanksStats_miniCheckpoints_long %>%
ggplot() +
  geom_point(
    aes(
      x = dataset,
      y = percAcrossRanksTotalTimeAvg
  )) +
  geom_errorbar(
    aes(
      ymin = percAcrossRanksTotalTimeAvg - percAcrossRanksTotalTimeSD,
      ymax = percAcrossRanksTotalTimeAvg + percAcrossRanksTotalTimeSD,
      x = dataset
  )) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = c(0.23, 0.75)
    #legend.title = element_blank()
  ) +
  ylab("time for mini-checkpointing\n[% of runtime]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_20@4" = "AA NagyA1\n80 ranks\n594 models",
      "aa_rokasA4_20@8" = "AA ChenA4\n160 ranks\n1 model",
      "aa_rokasA8_20@4" = "AA YangA8\n80 ranks\n1122 models",
      "dna_PeteD8_13@20" = "DNA PeteD8\n260 ranks\n4116 models",
      "dna_rokasD1_20@20" = "DNA SongD1\n400 ranks\n1 model",
      "dna_rokasD4_20@8" = "DNA XiD4\n160 ranks\n1 model",
      "dna_ShiD9_20@1" = "DNA ShiD9\n20 ranks\n29 models"
    ),
    limits = acrossRanksStats_miniCheckpoints_long %>% arrange(percAcrossRanksTotalTimeAvg) %>% pull(dataset)
  )
ggsave("~/Promotion/figures/paper/time-for-mini-checkpoints.pdf", width = 15, height = 8, units = "cm")

### Hard facts 
acrossRanksStats_miniCheckpoints_long %>%
  select(dataset, percAcrossRanksTotalTimeAvg, percAcrossRanksTotalTimeSD)
