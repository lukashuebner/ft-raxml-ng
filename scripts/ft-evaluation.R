library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(xtable)

dataDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/ft-evaluation/"
plotDir <- "/home/lukas/Documents/Uni/Masterarbeit/thesis/figures"
recoveryDataDir <- "/home/lukas/Documents/Uni/Masterarbeit/profiling/recovery"

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

recoveryStats <-
  list.files(
    path = recoveryDataDir,  
    pattern = "*.csv", 
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

# overallStats
nsSum_long <- overallStats %>% select(- starts_with("num")) %>% gather(key = timer, value = nsSum, starts_with("nsSum"))
num_long <- overallStats %>% select(- starts_with("nsSum")) %>% gather(key = timer, value = num, starts_with("num"))
nsSum_long$timer = str_replace(nsSum_long$timer, "nsSum", "")
num_long$timer = str_replace(num_long$timer, "num", "")
stats_long <- inner_join(by = c("rank", "processor", "dataset", "timer"), num_long, nsSum_long)
stats_long$nsSum <- stats_long$nsSum / stats_long$num

# recovery
nsSum_long <- recoveryStats %>% select(- starts_with("num")) %>% gather(key = timer, value = nsSum, starts_with("nsSum"))
num_long <- recoveryStats %>% select(- starts_with("nsSum")) %>% gather(key = timer, value = num, starts_with("num"))
nsSum_long$timer = str_replace(nsSum_long$timer, "nsSum", "")
num_long$timer = str_replace(num_long$timer, "num", "")
recovery_long <- inner_join(by = c("rank", "processor", "dataset", "timer"), num_long, nsSum_long)
recovery_long$nsSum <- recovery_long$nsSum / recovery_long$num


# overallStats
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

summary_stats <- stats_long %>%
  #filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "MiniCheckpoint"))) %>%
  filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "UpdateTree", "UpdateModels", "work"))) %>%
  group_by(dataset) %>%
  summarise(
    msRestoreSearchState = sum(msAvg)
  )

# recoveryStats
recovery_long <- recovery_long %>% group_by(dataset, timer) %>%
  summarise(
    nsAvg = mean(nsSum),
    nsSD = sd(nsSum),
  ) %>%
  mutate(
    msAvg = nsAvg / 10^6,
    msSD = nsSD / 10^6
  ) %>%
  arrange(msAvg)

recovery_summary_stats <- recovery_long %>%
  filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "UpdateTree", "UpdateModels", "work"))) %>%
  group_by(dataset) %>%
  summarise(
    msRestoreSearchState = sum(msAvg)
  )



stats_summary <- inner_join(
  stats_long %>%
    #filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "MiniCheckpoint"))) %>%
    filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "UpdateTree", "UpdateModels", "work"))) %>%
    group_by(dataset) %>%
    summarise(
      msRestoreSearchState = sum(msAvg),
      msSDRestoreSearchState = sum(msSD**2)
    ) %>%
    mutate(
      msSDRestoreSearchState = sqrt(msSDRestoreSearchState)
    ),
  dataDistribution <- tribble(
    ~dataset           , ~taxa, ~sites , ~ type,
    "aa_rokasA1_20@4"  ,  60  ,  172073, "AA",
    "aa_rokasA4_20@8"  , 160  , 1806035, "AA", 
    "aa_rokasA4_20@20"  , 160  , 1806035, "AA", 
    "aa_rokasA8_20@4"  ,  80  ,  504850, "AA",             
    "dna_PeteD8_13@20" , 174  , 3011099, "DNA",                  
    "dna_rokasD1_20@20",  37  , 1338678, "DNA",        
    "dna_rokasD2a_20@1", 144  , 1240377, "DNA",        
    "dna_rokasD4_20@8" , 160  ,  237363, "DNA",       
    "dna_ShiD9_20@1"   , 815  ,   20364, "DNA",
    "dna_rokasD7_20@20",  46  ,21410970, "DNA"
  ) %>% mutate(dataset = as.factor(dataset)),
  by = "dataset"
)

LoadMSA_profile <- stats_long %>% select(-starts_with("ns")) %>% filter(timer == "LoadAssignmentData")

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
    filter(!(timer %in% c("LoadAssignmentData", "RedoPartitionAssignment", "LoadAssignmentDataFirstLoad", "TreeinfoUpdatePartialsAndCLVs", "UpdateTree", "UpdateModels", "work",  "TreeinfoInit"))) %>%
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
      "aa_rokasA1_20@4" = "AA rokasA1 @ 80 ranks\n172k sites * 60 taxa",
      "aa_rokasA4_20@8" = "AA rokasA4 @ 160 ranks\n1806k sites * 58 taxa",
      "aa_rokasA4_20@20" = "AA rokasA4 @ 400 ranks\n1806k sites * 58 taxa",
      "aa_rokasA8_20@4" = "AA rokasA8 @ 80 ranks\n505k sites * 95 taxa",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n3011k sites * 174 taxa",
      "dna_rokasD1_20@20" = "DNA rokasD1 @ 400 ranks\n1339k sites * 37 taxa",
      "dna_rokasD2a_20@1" = "DNA rokasD2a @ 20 ranks\n1240k sites * 144 taxa",
      "dna_rokasD4_20@8" = "DNA rokasD4 @ 160 ranks\n240k sites * 46 taxa",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n20k sites * 815 taxa",
      "dna_rokasD7_20@20" = "DNA TarvD7 @ 400 ranks\n21M sites * 37 taxa"
    ),
    limits = (
      stats_long %>%
        #filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "MiniCheckpoint"))) %>%
        filter(!(timer %in% c("LoadAssignmentDataFirstLoad", "UpdateTree", "UpdateModels", "work"))) %>%
        
        group_by(dataset) %>%
        summarise(msAvg = sum(msAvg)) %>%
        arrange(msAvg)
    )$dataset
  ) + 
  scale_fill_brewer(type="qual") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

ggsave(
  filename = "time-for-fault-recovery.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)

ggplot(filter(stats_long, timer == "UpdateModels", dataset != "aa_KatzA10_20@2", dataset != "aa_GitzA12_20@1")) +
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
      "aa_rokasA1_20@4" = "AA NagyA1 @ 80 ranks\n594 models",
      "aa_rokasA4_20@8" = "AA ChenA4 @ 160 ranks\n1 model",
      "aa_rokasA4_20@20" = "AA ChenA4 @ 400 ranks\n1 model",
      "aa_rokasA8_20@4" = "AA YangA8 @ 80 ranks\n1122 models",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n4116 models",
      "dna_rokasD1_20@20" = "DNA SongD1 @ 400 ranks\n1 model",
      "dna_rokasD2a_20@1" = "DNA MisoD2a @ 20 ranks\n100 models",
      "dna_rokasD4_20@8" = "DNA XiD4 @ 160 ranks\n1 model",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n29 models",
      "dna_rokasD7_20@20" = "DNA TarvD7 @ 400 ranks\n1 model"
    ),
    limits = (
      stats_long %>%
        filter(timer == "UpdateModels", dataset != "aa_KatzA10_20@2", dataset != "aa_GitzA12_20@1") %>%
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

ggplot(  ) +
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
  ylab("time to backup the current best tree [ms]") +
  scale_x_discrete(
    labels = c(
      "aa_rokasA1_20@4" = "AA NagyA1 @ 80 ranks\n60 taxa",
      "aa_rokasA4_20@8" = "AA ChenA4 @ 160 ranks\n58 taxa",
      "aa_rokasA4_20@20" = "AA ChenA4 @ 400 ranks\n58 taxa",
      "aa_rokasA8_20@4" = "AA YangA8 @ 80 ranks\n95 taxa",
      "dna_PeteD8_13@20" = "DNA PeteD8 @ 260 ranks\n174 taxa",
      "dna_rokasD1_20@20" = "DNA SongD1 @ 400 ranks\n37 taxa",
      "dna_rokasD4_20@8" = "DNA XiD4 @ 160 ranks\n46 taxa",
      "dna_rokasD7_20@20" = "AA TarvD7 @ 400 ranks\n36 taxa",
      "aa_KatzA10_20@2" = "AA KatzA10 @ 40 ranks\n798 taxa",
      "dna_ShiD9_20@1" = "DNA ShiD9 @ 20 ranks\n815 taxa",
      "aa_GitzA12_20@1" = "AA GitzA12 @ 20 ranks\n1897 taxa"
    ),
    limits = (
      stats_long %>%
        filter(timer == "UpdateTree") %>%
        arrange(msAvg)
    )$dataset
  )

ggsave(
  filename = "time-for-updating-tree.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)

updateTree <- filter(stats_long, timer == "UpdateTree")
updateTree$ranks <- c(400,160,400,400,160,80,80,260)
updateTree$taxa <- c(37,46,36,58,58,60,95,174)
corr.test(updateTree$taxa, updateTree$msAvg)
