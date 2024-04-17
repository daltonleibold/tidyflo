# test.R

#-------------------
# bookkeeping
#-------------------

# load libraries
pacman::p_load(tidyverse)

# load this library
devtools::load_all()

#-------------------
# processing
#-------------------

# make a path to the data
path <- "./raw_data"

# load in the data
data <- read.fcs(path) %>%
  tidy.fcs() %>%
  process.fcs() %>%
  separate(.id, into = c(".spec", ".samp"), sep = "_") %>%
  mutate(.ch = re.code(.ch, c("alexa.fluor.488" = "af",
                              "apc" = "dr",
                              "pe" = "mto",
                              "percp.cy5.5" = "sox",
                              "buv.496" = "h42")))

data %>%
  group_by(.spec, .samp, .ch) %>%
  summary.fcs() %>%
  filter(.spec == "scc" & .samp != "pos") %>%
  filter(.ch != "fsc" & .ch != "ssc") %>%
  group_by(.samp, .ch) %>%
  mutate(.comp = 10^mean) %>%
  reframe(.comp = ifelse(.samp == "neg" & .ch == "af", 0,
                         ifelse(.samp == .ch, 0, .comp))) %>%
  mutate(across(.cols = c(.samp, .ch),
                .fns = ~factor(.x, levels = c("neg",
                                              "af",
                                              "dr",
                                              "mto",
                                              "sox",
                                              "h42")))) %>%
  #group_by(.ch) %>%
  #reframe(.comp = sum(.comp))
  # make a compensation matrix
  ggplot(aes(x = .ch, y = .samp, 
             fill = .comp, group = .ch)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "coral",
                       mid = "gold",
                       high = "seagreen") +
  geom_text(aes(label = round(.comp, 2)),
            color = "black", size = 4) + 
  xlab("Channel") +
  ylab("Control") +
  guides(fill = "none") +
  theme_classic()

#-------------------
#
#-------------------

#-------------------
#
#-------------------

#-------------------
#
#-------------------

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


