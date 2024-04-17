# test.R

# load libraries
pacman::p_load(tidyverse)

# source functions
source("func.R")

#-------------------
# graphing
#-------------------
#----
#' @title theme_custom
#' @description a custom theme for ggplots that makes publication-ready figures
theme_custom <- theme_classic(base_size = 24) +
  # adjust axis title position
  theme(axis.title.y=element_text(vjust=1.5), 
        axis.title.x=element_text(vjust=0.2)) + 
  # adjust plot margins and line element size
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), 
        line = element_line(linewidth = 1.25)) + 
  # draw x and y axes
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) + 
  #put margins around axis labels so that nothing overlaps
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) + 
  # move tickmarks inside the axes and paint black
  theme(axis.ticks.length =unit(-0.3, "cm")) + 
  #spread out facets
  theme(panel.spacing = unit(2, units = "lines")) + 
  #make tick marks black
  theme(axis.ticks = element_line(color = "black")) + 
  #remove border from facet labels
  theme(strip.background = element_blank()) 



#-------------------
# compensation
#-------------------
# set a path to the controls
path.scc <- "../controls/"
# read in the compensation controls
data.scc <- path.scc %>%
  read.it() %>%
  # arrange data for processing
  pivot_longer(-.id,
               names_to = ".ch",
               values_to = ".val") %>%
  group_by(.id, .ch) %>%
  # remove off-scale events
  filter(.val >= 1 & .val <= 250000) %>%
  arrange(.id, .ch, .val) %>%
  # identify positives and negatives as being >2sd from median
  mutate(.row = row_number(),
         .n = length(.row)) %>%
  mutate(.type = ifelse(.row < ceiling(.n * 0.0228), ".neg",
                        ifelse(.row > ceiling(.n * 0.9772), ".pos", NA))) %>%
  filter(!is.na(.type)) %>%
  # calculate the geometric mean for positives and negatives
  group_by(.id, .ch, .type) %>%
  mutate(.val = log10(.val)) %>%
  reframe(.val = geo.mean(.val)) %>%
  mutate(.val = 10^.val) %>%
  pivot_wider(names_from = ".type",
              values_from = ".val") %>%
  # for each sample, estimate the compensation value for each channel
  group_by(.id, .ch) %>%
  reframe(.comp = .pos - .neg) %>%
  unite(col = .ch.id,
        c(.ch, .id),
        sep = ".") %>%
  pivot_wider(names_from = ".ch.id",
              values_from = ".comp") %>%
  # generate the compensation matrix 
  reframe(fsc = fsc.neg - (fsc.dr + fsc.mto +
                             fsc.sox + fsc.h42),
          ssc = ssc.neg - (ssc.dr + ssc.mto +
                             ssc.sox + ssc.h42),
          af = af.neg - (af.dr + af.mto +
                           af.sox + af.h42),
          dr = dr.dr - (dr.neg + dr.mto +
                          dr.sox + dr.h42),
          mto = mto.mto - (mto.neg + mto.dr +
                             mto.sox + mto.h42),
          sox = sox.sox - (sox.neg + sox.dr +
                             sox.mto + sox.h42),
          h42 = h42.h42 - (h42.neg + h42.dr +
                             h42.mto + h42.sox)) %>%
  # pivot for joining to sample data
  pivot_longer(cols = everything(),
               names_to = ".ch",
               values_to = ".comp")

#-------------------
# data processing
#-------------------
# make path to samples
path.samp <- "../samples/"

# read in the samples
data.samp <- path.samp %>%
  read.it() %>%
  # apply channel labels
  is.mitofx() %>%
  # add cell identifier for processing
  group_by(.id) %>%
  mutate(.cell = row_number()) %>%
  ungroup() %>%
  select(.id, .cell, everything()) %>%
  # pivot to long-form for easy processing
  pivot_longer(3:ncol(.),
               names_to = ".ch",
               values_to = ".val") %>%
  # append compensation matrix
  left_join(data.scc, by = ".ch") %>%
  # compensate channel values
  group_by(.id, .ch, .cell) %>%
  reframe(.val = .val - .comp) %>%
  # pivot back to wide-form
  pivot_wider(names_from = ".ch",
              values_from = ".val",
              id_cols = c(".id", ".cell")) %>%
  select(-.cell) %>%
  # remove off-scale values
  group_by(.id) %>%
  filter(if_all(.cols = everything(),
                .fns = ~.x >= 1 & .x <= 250000)) %>%
  # calculate summary statistics for each channel in each sample
  mutate(across(.cols = everything(),
                .fns = log10)) %>%
  reframe(across(.cols = everything(),
                 .fns = list(geo.mean = geo.mean,
                             sd = r.sd),
                 .names = "{col}.{fn}")) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, ".geo.mean"))

write_csv(data.samp, "../data.csv")
#~~~~~~~~~~~~~~~~~~~~~~~

#-----------------------
# normalizing
#-----------------------

# make simple intercept models
bf.fsc <- bf(fsc | sd(fsc.sd) ~ 0 + .spec)
bf.ssc <- bf(ssc | sd(ssc.sd) ~ 0 + .spec)
bf.af <- bf(af | sd(af.sd) ~ 0 + .spec)
bf.h42 <- bf(h42 | sd(h42.sd) ~ 0 + .spec)
bf.dr <- bf(dr | sd(dr.sd) ~ 0 + .spec)
bf.mto <- bf(mto | sd(mto.sd) ~ 0 + .spec)
bf.sox <- bf(sox | sd(sox.sd) ~ 0 + .spec)

mod.mitofx <- brm(bf.fsc +
                    bf.ssc +
                    bf.af +
                    bf.h42 +
                    bf.dr +
                    bf.mto +
                    bf.sox +
                    set_rescor(rescor = T),
                  data = data.samp,
                  chains = 4, cores = 4,
                  iter = 1000, warmup = 500,
                  file_refit = "on_change")

data.norm <- posterior_samples(mod.mitofx) %>%
  select(grep("^b_", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^b_")) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, ".id")) %>%
  pivot_longer(cols = everything(),
               names_to = ".ch.id",
               values_to = ".val") %>%
  separate(col = .ch.id,
           into = c(".ch", ".id"),
           sep = "_") %>%
  separate(.id, into = c(".spec", ".samp"),
           sep = 1) %>%
  group_by(.spec, .samp, .ch) %>%
  mutate(.tracer = row_number()) %>%
  pivot_wider(names_from = ".samp",
              values_from = ".val") %>%
  group_by(.ch, .spec, .tracer) %>%
  reframe(.val = (a + b)/2) %>%
  pivot_wider(names_from = ".ch",
              values_from = ".val") %>%
  select(-.tracer)

#-----------------
# visualisation
#-----------------

data.samp %>%
  pivot_longer(4:ncol(.),
               names_to = ".ch",
               values_to = ".val") %>%
  ggplot(aes(x = .val, fill = .spec, alpha = 0.2)) +
  geom_density() +
  facet_wrap(~.ch, scales = "free_y") +
  theme_classic()

data.norm %>%
  pivot_longer(2:ncol(.),
               names_to = ".ch",
               values_to = ".val") %>%
  ggplot(aes(x = .val, fill = .id, alpha = 0.2)) +
  geom_density() +
  facet_wrap(~.ch, scales = "free_y") +
  theme_classic()

data.samp %>%
  pivot_longer(3:ncol(.),
               names_to = ".ch",
               values_to = ".val") %>%
  ggplot(aes(x = .ch, y = .val)) +
  geom_violin() +
  stat_summary(fun = geo.mean, shape = 21) +
  stat_summary(fun = median) +
  theme_classic()


data.samp %>%
  ggplot(aes(x = fsc, y = ssc, col = h42)) +
  geom_point() +
  scale_color_gradient(low = "darkblue",
                       high = "cyan") +
  theme_classic()

data.samp %>%
  ggplot(aes(x = mto, y = sox, col = dr)) +
  geom_point() +
  scale_color_gradient(low = "darkred",
                       high = "red") +
  theme_classic()

data.samp %>%
  separate(col = .id, 
           into = c(".spec", ".samp"),
           sep = 1) %>%
  group_by(.spec) %>%
  reframe(across(.cols = af:ssc,
                 .fns = geo.mean)) %>%
  glimpse()



data.norm %>%
  group_by(.spec) %>%
  reframe(across(.cols = everything(),
                 .fns = geo.mean)) %>%
  glimpse()

data.norm %>%
  ggplot(aes(x = fsc, y = ssc, 
             col = h42, shape = .spec)) +
  geom_point() + 
  scale_color_gradient(low = "cyan",
                       high = "darkblue") +
  theme_classic()

data.norm %>%
  ggplot(aes(x = mto, y = sox, 
             col = dr, shape = .spec)) +
  geom_point() + 
  scale_color_gradient(low = "red",
                       high = "darkred") +
  theme_classic()
