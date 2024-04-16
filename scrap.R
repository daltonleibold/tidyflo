# scrap.R
# some scrap script for testing stuff out

pacman::p_load(tidyverse)
source("func.R")


#----------------
# processing
#----------------
#----
# make path to data
path <- "./test.fcs"
# read in data
data <- read.fcs(path) %>%
  tidy.fcs() %>%
  process.fcs() %>%
  mutate(.ch = re.code(.ch, mapping = c("alexa.fluor.488" = "af",
                                        "apc" = "dr",
                                        "pe" = "mto",
                                        "percp.cy5.5" = "sox",
                                        "buv.496" = "h42")))

data %>%
  summary.fcs()

#----------------
# visualising
#----------------
#----
plot.wave(data)
plot.density(data)
plot.wave(data, batch = T)
plot.density(data, batch = T)
#----