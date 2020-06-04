library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(xtable)

dataDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/ft-evaluation"
plotDir <- "/home/lukas/Documents/Uni/Masterarbeit/thesis/figures"
rawDir <- "/home/lukas/Documents/Uni/Masterarbeit/thesis/rawdata"

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

jkruntimeStats <- read_csv(
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
  ) %>%
  arrange(msAvg)

D9_profile <- stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "dna_ShiD9_20@1") %>% arrange(desc(msAvg))
A8_profile <- stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "aa_rokasA8_20@4") %>% arrange(desc(msAvg))
D2a_profile <- stats_long %>% select(-starts_with("ns")) %>% filter(dataset == "dna_rokasD2a_20@1") %>% arrange(desc(msAvg))
print.xtable(xtable(D9_profile), type = "latex", file = paste(rawDir, "ft-profile-D9.tex"), sep = "/")
print.xtable(xtable(A8_profile), type = "latex", file = paste(rawDir, "ft-profile-A8.tex"), sep = "/")
print.xtable(xtable(D2a_profile), type = "latex", file = paste(rawDir, "ft-profile-D2a.tex"), sep = "/")

LoadMSA_profile <- stats_long %>% select(-starts_with("ns")) %>% filter(timer == "LoadAssignmentData")
print.xtable(xtable(LoadMSA_profile), type = "latex", file = paste(rawDir, "ft-profile-loadMSA.tex"), sep = "/")

ggplot() +
  geom_errorbar(
    data = filter(stats_long, timer == "LoadAssignmentData"),
    aes(ymin = msAvg - msSD, ymax = msAvg + msSD, x = dataset, color = "further load op.")
  ) +
  geom_point(
    data = filter(stats_long, timer == "LoadAssignmentData"),
    aes(y = msAvg, x = dataset, color = "further load op.")
  ) +
  geom_errorbar(
    data = filter(stats_long, timer == "LoadAssignmentDataFirstLoad"),
    aes(ymin = msAvg - msSD, ymax = msAvg + msSD, x = dataset, color = "initial load op.")
  ) +
  geom_point(
    data = filter(stats_long, timer == "LoadAssignmentDataFirstLoad"),
    aes(y = msAvg, x = dataset, color = "initial load op.")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 47, hjust = 1),
    legend.position = c(0.15, 0.85),
    legend.title = element_blank()
  ) +
  ylab("time [ms]") +
  scale_x_discrete(
    labels = c(
      # "aa_rokasA1_10@4" = "AA rokasA1 @ 40 ranks\n172k sites * 60 taxa",
      # "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n1806k sites * 58 taxa",
      # "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n505k sites * 95 taxa",
      # "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n3011k sites * 174 taxa",
      # "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n1339k sites * 37 taxa",
      # "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n1240k sites * 144 taxa",
      # "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n240k sites * 46 taxa",
      # "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n20k sites * 815 taxa"
      "aa_rokasA1_10@4" = "AA rokasA1 @ 40 ranks\n10 MiB",
      "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n100 MiB",
      "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n46 MiB",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n500 MiB",
      "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n48 MiB",
      "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n171 MiB",
      "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n11 MiB",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n16 MiB"
    ),
    limits = (stats_long %>% filter(timer == "LoadAssignmentDataFirstLoad") %>% arrange(msAvg))$dataset
  )
ggsave(
  filename = "time-for-loading-msa.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)

bind_rows(stats_long %>%
    filter(!(timer %in% c("LoadAssignmentData", "RedoPartitionAssignment", "LoadAssignmentDataFirstLoad", "TreeinfoUpdatePartialsAndCLVs", "MiniCheckpoint", "TreeinfoInit"))) %>%
    group_by(dataset) %>%
    summarise(msAvg = sum(msAvg)) %>%
    mutate(timer = "Other"),
  stats_long %>%
    filter(timer %in% c("LoadAssignmentData", "RedoPartitionAssignment", "TreeinfoUpdatePartialsAndCLVs", "TreeinfoInit")) %>%
    select(dataset, timer, msAvg)
) %>%
ggplot() +
  geom_bar(aes(x = dataset, y = msAvg, fill = timer), stat = "identity") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = c(0.2, 0.6)
    #legend.title = element_blank()
  ) +
  ylab("time [ms]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_10@4" = "AA rokasA1 @ 40 ranks\n172k sites * 60 taxa",
      "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n1806k sites * 58 taxa",
      "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n505k sites * 95 taxa",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n3011k sites * 174 taxa",
      "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n1339k sites * 37 taxa",
      "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n1240k sites * 144 taxa",
      "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n240k sites * 46 taxa",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n20k sites * 815 taxa"
    ),
    limits = (
      stats_long %>%
        filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "MiniCheckpoint"))) %>%
        group_by(dataset) %>%
        summarise(msAvg = sum(msAvg)) %>%
        arrange(msAvg)
    )$dataset
  ) + 
  scale_fill_brewer(type="qual")

ggsave(
  filename = "time-for-fault-recovery.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)

ggplot(filter(stats_long, timer == "MiniCheckpoint")) +
  geom_point(aes(x = dataset, y = msAvg)) +
  geom_errorbar(aes(ymin = msAvg - msSD, ymax = msAvg + msSD, x = dataset)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = c(0.23, 0.75)
    #legend.title = element_blank()
  ) +
  ylab("time for mini-checkpoint [ms]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_10@4" = "AA rokasA1 @ 40 ranks\n594 models",
      "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n1 model",
      "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n1122 models",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n4116 models",
      "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n1 model",
      "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n100 models",
      "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n1 model",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n29 models"
    ),
    limits = (
      stats_long %>%
        filter(timer == "MiniCheckpoint") %>%
        arrange(msAvg)
    )$dataset
  )

ggsave(
  filename = "time-for-mini-checkpoint.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)

  