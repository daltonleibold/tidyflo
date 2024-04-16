#-------------------
# comp.R
# compensate flow data
#-------------------

#-------------------
# bookkeeping
#-------------------

# load libraries
pacman::p_load(tidyverse, brms, rstan)

# source functions
source("func.R")

#-------------------
# processing
#-------------------

data <- read.fcs("../data") %>%
  # standardize the dataframe
  tidy.fcs() %>%
  # run preliminary processing
  process.fcs() %>%
  # only analysing the area channels
  filter(.param == ".a") %>% 
  # separate the id column into the specimen and the sample
  separate(.id, into = c(".spec", ".samp"), sep = "_") %>%
  # apply channel labels
  mutate(.ch = re.code(.ch, mapping = c("alexa.fluor.488" = "af",
                                        "apc" = "dr",
                                        "pe" = "mto",
                                        "percp.cy5.5" = "sox",
                                        "buv.496" = "h42"))) %>%
  # summarise the dataframe - need geo.mean and r.se for modelling
  group_by(.spec, .samp, .param, .ch) %>%
  summary.fcs()


bf.norm <- bf(geo.mean | se(r.se) ~ 0 + .ch * .spec + (1|.spec:.samp)) + gaussian()

mod.norm <- brm(bf.norm,
      data = data,
      chains = 4, cores = 4,
      iter = 5000, warmup = 2500,
      file_refit = "always",
      file = "../output/norm.rds")
  
mod.norm %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(grep("^b_", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^b_")) %>%
  glimpse()

#-------------------
# compensate
#-------------------

# run the intercept model
mod.scc <- brm(bf.scc,
               data = data,
               chains = 4, cores = 4,
               iter = 5000, warmup = 2500,
               file_refit = "on_change",
               file = "../output/scc.rds")

#-------------------
# normalize
#-------------------

# extract compensated values
data.norm <- mod.scc %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(grep("^b_", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^b_")) %>%
  select(grep(":", colnames(.))) %>%
  pivot_longer(cols = everything(),
               names_to = ".id:.ch",
               values_to = ".val") %>%
  separate(`.id:.ch`, into = c(".id", ".ch"), sep = ":") %>%
  mutate(.id = str_remove(.id, ".id"),
         .ch = str_remove(.ch, ".ch")) %>%
  separate(.id, into = c(".spec", ".samp"), sep = "_") %>%
  filter(.spec != "scc") %>%
  group_by(.spec, .samp, .ch) %>%
  reframe(mean.val = mean(.val),
          se.val = se(.val)) %>%
  mutate(.spec = factor(.spec),
         .samp = factor(.samp),
         .ch = factor(.ch))

# make a formula to normalize across specimen and samples within specimen
bf.norm <- bf(mean.val | se(se.val) ~ 0 + .ch * .spec + (1|.spec:.samp)) + gaussian()

# run the normalizing model
mod.norm <- brm(bf.norm,
               data = data.norm,
               chains = 4, cores = 4,
               iter = 5000, warmup = 2500,
               file_refit = "on_change",
               file = "../output/scc.rds")

post.norm <- mod.norm %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(grep("b_", colnames(.)))

post.norm %>%
  glimpse()

#-------------------
# extracting
#-------------------

post.scc <- extract.scc(mod.scc) %>%
  mutate(.val = 10^.val) %>%
  group_by(.ch, .samp) %>%
  reframe(.val = geo.mean(.val))

model <- mod.scc
# extract posteriors function (for troubleshooting)
model %>% 
  as_draws_df() %>%
  glimpse()
  select(grep("^b_", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^b_")) %>%
  select(grep(":", colnames(.))) %>%
  glimpse()
  as_tibble() 
# true negative
# autofluorescence channel in the unstained sample
a. <- post %>%
  reframe(.ch = neg,
          .samp = neg,
          .val = Intercept)
# autofluorescence spillover / channel background
# fluorescent signal in each channel of the unstained control
b. <- post %>%
  select(Intercept, grep("^.ch", colnames(.))) %>%
  select(-grep(":", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^.ch")) %>%
  pivot_longer(-Intercept,
               names_to = ".ch",
               values_to = ".val") %>%
  mutate(.val = Intercept + .val,
         .samp = "af") %>%
  select(.ch, .samp, .val)
# dye-induce autofluorescence
# autofluorescent channel in each single color control
c. <- post %>%
  select(Intercept, grep("^.samp", colnames(.))) %>%
  select(-grep(":", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^.samp")) %>%
  pivot_longer(-Intercept,
               names_to = ".samp",
               values_to = ".val") %>%
  mutate(.val = Intercept + .val,
         .ch = neg) %>%
  select(.ch, .samp, .val)
# fluorescent spillover
# each of the fluorescent channels for each of the single color controls
d. <- post %>%
  select(Intercept, grep(":", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "^.samp")) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, ".ch")) %>%
  pivot_longer(-Intercept,
               names_to = ".samp.ch",
               values_to = ".val") %>%
  separate(col = .samp.ch,
           into = c(".samp", ".ch"),
           sep = ":") %>%
  select(.ch, .samp, .val)
# bind the extracted distributions back together
data <- rbind(a., b., c., d.)

#-------------------
# visualizing
#-------------------

# make model compensation matrix
data.mod %>%
  group_by(.samp, .ch) %>%
  # log10 transform to linearize
  mutate(.val = log10(abs(.val))) %>%
  # calculate the geometric mean
  reframe(.val = geo.mean(abs(.val))) %>%
  # exponentiate into original units
  mutate(.val = 10^.val) %>%
  # make values negative if the channel should be negative
  mutate(.val = ifelse(.samp == .ch, 
                       .val, -.val)) %>%
  # make a compensation matrix
  ggplot(aes(x = .ch, y = .samp, 
             fill = .val, group = .ch)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "coral",
                       mid = "gold",
                       high = "seagreen") +
  geom_text(aes(label = round(.val, 2)),
            color = "black", size = 4) + 
  xlab("Channel") +
  ylab("Control") +
  guides(fill = "none") +
  theme_custom

# save model compensation matrix
ggsave(filename = "../output/comp_matrix.png",
       width = 9, height = 6, dpi = 300)


