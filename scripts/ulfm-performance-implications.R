library(tibble)
library(readr)
library(ggplot2)
library(dplyr)
library(scales)

plotDir <- "/home/lukas/Documents/Uni/Masterarbeit/thesis/figures"
dataDir <- "/home/lukas/Documents/Uni/Masterarbeit/thesis/rawdata"

# Comparing runtimes of different MPI implementations and fault detection mechanisms
mpi_runtimes <- tribble(
  ~mpiVersion, ~dataset, ~runtime,
  "ULFM-no-ft", "aa_rokasA4_20@8", 5593.958,
  "ULFM-no-ft", "dna_rokasD1_20@18", 1463.537,
  "ULFM-no-ft", "dna_ShiD9_20@1", 22238.792,
  "OpenMPI-4.0", "aa_rokasA4_20@8", 5582.245,
  "OpenMPI-4.0", "dna_rokasD1_20@18", 1436.832,
  "OpenMPI-4.0", "dna_ShiD9_20@1", 20910.897,
  "ULFM-ft-thread", "aa_rokasA4_20@8", 5549.155,
  "ULFM-ft-thread", "dna_rokasD1_20@18", 1471.628,
  "ULFM-ft-thread", "dna_ShiD9_20@1", 22711.731,
)

mpi_runtimes <- inner_join(
  mpi_runtimes,
  mpi_runtimes %>% filter(mpiVersion == "OpenMPI-4.0") %>% rename(reference = runtime) %>% select(-mpiVersion),
  by = "dataset"
)

mpi_runtimes <- mpi_runtimes %>% mutate(speedup = reference / runtime)

# Plotting recovery time after failure
recoveryTime <- read.csv("~/Documents/Uni/Masterarbeit/profiling/faultingSpeed.csv") %>% as_tibble()
recoveryTimeLabels <- c(
  dt = "hbt ON, 300 ms timeout\n4 nodes á 7 ranks *",
  dt_s = "hbt ON, 1000 ms timeout\n4 nodes á 7 ranks *",
  dt_big = "hbt ON, 300 ms timeout\n20 nodes á 20 ranks",
  full_dt = "hbt ON, 300 ms timeout\n4 nodes á 20 ranks",
  full_no_dt = "hbt OFF, 300 ms timeout\n4 nodes á 20 ranks",
  no_dt = "hbt OFF, 300 ms timeout\n4 nodes á 7 ranks *"
)

ggplot(recoveryTime) +
  geom_histogram(aes(y = stat(density) * 0.1, x = time / 1000), binwidth = 0.1) +
  scale_y_continuous(label = percent) +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~method, labeller = labeller(method = recoveryTimeLabels)) +
  labs(x = "time [s]", y = "frequency")

recoveryTime$timePoint <- 1
recoveryTime <- recoveryTime %>% group_by(method) %>% mutate(timePoint = cumsum(timePoint)) %>% ungroup()

recoveryTime

ggplot(recoveryTime) +
  geom_point(aes(x = timePoint, y = time / 1000)) +
  facet_wrap(
    ~method,
    scale = "free_x",
    labeller = labeller(method = recoveryTimeLabels)
  ) +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 2 * 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 3 * 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 4 * 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 5 * 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 6 * 27), color = "darkgrey") +
  geom_vline(data = filter(recoveryTime, method == "no_dt"), aes(xintercept = 7 * 27), color = "darkgrey") +
  theme_bw() +
  labs(x = "n-th failed node", y = "recover time [s]")

ggplot(recoveryTime) +
  geom_violin(aes(x = method, y = time / 1000)) +
  scale_x_discrete(labels = recoveryTimeLabels) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.01)) +
  labs(x = "", y = "recovery time [s]")
ggsave(
  filename = "recovery-violin.pdf",
  path = plotDir,
  device = "pdf",
  width = 16.7976,
  height = 10.4332,
  units = "cm"
)
recoveryTime %>% group_by(method) %>% summarise(cnt = n())

recoveryTime %>% group_by(method) %>% summarise(median = median(time), sd = sd(time))
recoveryTime %>% filter(method == "dt_big", time > 2000)

sr_runtimes <- tribble(
  ~siteRepeats, ~dataset, ~runtime,
  "site-repeats OFF", "dna_rokasD1_20@18", 5632.729,
  "site-repeats OFF", "dna_rokasD2a_20@1", 65794.786,
  "site-repeats OFF", "dna_ShiD9_20@1", 48450.181,
  "site-repeats OFF", "aa_rokasA4_20@8", 7718.041,
  "site_repeats ON", "dna_rokasD1_20@18", 1351.651,
  "site-repeats ON", "dna_rokasD2a_20@1", 38806.859,
  "site-repeats ON", "dna_ShiD9_20@1", 20910.897,
  "site-repeats ON", "aa_rokasA4_20@8", 5582.245
)

sr_runtimes <- inner_join(
  sr_runtimes,
  sr_runtimes %>% filter(siteRepeats == "site-repeats OFF") %>% rename(reference = runtime) %>% select(-siteRepeats),
  by = "dataset"
)

sr_runtimes <- sr_runtimes %>% mutate(speedup = reference / runtime) %>% select(-reference)
