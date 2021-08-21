## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(spatialTIME)

## ---- create------------------------------------------------------------------

x <- spatialTIME::create_mif(clinical_data = example_clinical,
                sample_data = example_summary,
                spatial_list = example_spatial,
                patient_id = "deidentified_id", 
                sample_id = "deidentified_sample")

x #prints a summary of how many patients, samples, and spatial files are present


## ---- plot_bad, fig.width=10, fig.height=6,fig.align = 'left', fig.cap='Notice that using the "bad" (B) order that one would have the impression that there is a huge number of CD3+ and would scratch there head over why there are no CD3+ FOXP3+ or CD3+ CD8+. While the "good" (G) order can  better make this distiction as well as show cells that are only positive for CD3 and not postive for FOXP3 or CD8.'----
mnames_bad <- c("CD3..CD8.","CD3..FOXP3.","CD3..Opal.570..Positive",
                "CD8..Opal.520..Positive","FOXP3..Opal.620..Positive", 
                "PDL1..Opal.540..Positive", "PD1..Opal.650..Positive")

# Used to make the legends in both plots below be in same order and use the 
# same coloring scheme for the purpose making a common legend

values = viridis::turbo(length(mnames_bad))
names(values) = mnames_bad

x<- spatialTIME::plot_immunoflo(x, plot_title = "deidentified_sample",  mnames = mnames_bad,
                   cell_type = "Classifier.Label")

bad_names <- x[["derived"]][["spatial_plots"]][[4]] + 
  ggplot2::theme(legend.position = 'bottom') + 
  ggplot2::scale_color_manual(breaks = mnames_bad, values = values)

mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
                 "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
                 "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")

x <- spatialTIME::plot_immunoflo(x, plot_title = "deidentified_sample", mnames = mnames_good, 
                    cell_type = "Classifier.Label")

good_names <- x[["derived"]][["spatial_plots"]][[4]] + 
  ggplot2::theme(legend.position = 'bottom') + 
  ggplot2::scale_color_manual(breaks = mnames_good, 
                     values = values[match(mnames_good, names(values))])

x$sample %>% dplyr::filter(deidentified_sample == 'TMA3_[9,K].tif') %>% 
  dplyr::select(c(2, 4:15)) %>%
  tidyr::pivot_longer(cols = 2:13, names_to = 'Marker', values_to = 'Count')

ggpubr::ggarrange(plotlist = list(bad_names, good_names), labels = c('B', 'G'),
                  common.legend = TRUE, legend = 'bottom')


## ---- ripleys, warning = FALSE, fig.width=10, fig.height=6, fig.align = 'center'----

x <- spatialTIME::ripleys_k(mif = x, mnames = mnames_good, num_permutations = 10,
               edge_correction = 'translation', r = seq(0,100,10),
               keep_perm_dis = FALSE, workers = 1)

# This will keeps the colors in every plot for the remainder of the vignette compatable 
values = viridis::turbo(length(unique(x$derived$univariate_Count$deidentified_sample)))
names(values) = unique(x$derived$univariate_Count$deidentified_sample)

x$derived$univariate_Count  %>%
  dplyr::filter(Marker != 'PDL1..Opal.540..Positive') %>%
  ggplot2::ggplot(ggplot2::aes(x = r, y = `Degree of Clustering Permutation`)) +
  ggplot2::geom_line(ggplot2::aes(color = deidentified_sample), show.legend = FALSE) +
  ggplot2::facet_wrap(Marker~., scales = 'free') + ggplot2::theme_bw() + 
  ggplot2::scale_color_manual(values = values)
  


## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x <- spatialTIME::bi_ripleys_k(mif = x, mnames = c("CD3..CD8.", "CD3..FOXP3."), num_permutations = 10,
               edge_correction = 'translation', r = seq(0,100,10),
               keep_perm_dis = FALSE, workers = 1, exhaustive = TRUE)

x$derived$bivariate_Count  %>%
  dplyr::filter(anchor == 'CD3..FOXP3.') %>%
  ggplot2::ggplot(ggplot2::aes(x = r, y = `Degree of Clustering Permutation`)) +
  ggplot2::geom_line(ggplot2::aes(color = deidentified_sample), show.legend = FALSE) +
  ggplot2::theme_bw() + ggplot2::scale_color_manual(values = values)

## ---- NN, warning = FALSE, fig.width=10, fig.height=6, fig.align = 'center'----

x <- spatialTIME::NN_G(mif = x, mnames = mnames_good, num_permutations = 10,
                edge_correction = 'rs', r = seq(0,100,10),
                keep_perm_dis = FALSE, workers = 1)

x$derived$univariate_NN  %>%
  dplyr::filter(Marker != 'PDL1..Opal.540..Positive') %>%
  ggplot2::ggplot(ggplot2::aes(x = r, y = `Degree of Clustering Permutation`)) +
  ggplot2::geom_line(ggplot2::aes(color = deidentified_sample)) +
  ggplot2::facet_wrap(Marker~., scales = 'free') + ggplot2::theme_bw() + 
  ggplot2::scale_color_manual(values = values)
  


## ----fig.width=10, fig.height=6, fig.align = 'center'-------------------------
x <- spatialTIME::bi_NN_G(mif = x, mnames = c("CD3..CD8.", "CD3..FOXP3."), num_permutations = 10,
               edge_correction = 'rs', r = seq(0,100,10),
               keep_perm_dis = FALSE, workers = 1)

x$derived$bivariate_NN  %>%
  dplyr::filter(anchor == 'CD3..FOXP3.') %>%
  ggplot2::ggplot(ggplot2::aes(x = r, y = `Degree of Clustering Permutation`)) +
  ggplot2::geom_line(ggplot2::aes(color = deidentified_sample), show.legend = FALSE) +
  ggplot2::theme_bw() +  ggplot2::scale_color_manual(values = values)

