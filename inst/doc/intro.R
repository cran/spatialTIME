## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


library(spatialTIME)

## ---- create------------------------------------------------------------------

x <- create_mif(clinical_data = example_clinical,
                sample_data = example_summary,
                spatial_list = example_spatial,
                patient_id = "deidentified_id", sample_id = "deidentified_sample",
                clean_columns = TRUE)
x


## ---- plot, fig.width=10, fig.height=10,fig.align = 'center'------------------
mnames <- c("foxp3_opal_620_positive", "cd3_opal_570_positive", "cd8_opal_520_positive",
            "pd1_opal_650_positive", "pdl1_opal_540_positive")

mlabels <- c("FOXP3", "CD3", "CD8", "PD1", "PDL1")

x[["derived"]][["spatial_plots"]] <- plot_immunoflo(x, plot_title = "deidentified_sample", 
                                                    mnames = mnames, mlabels = mlabels, 
                                                    cell_type = "classifier_label")

x[["derived"]][["spatial_plots"]][[4]]

## ---- ripleys, warning = FALSE, fig.width=13, fig.height=20, fig.align = 'center', fig.cap="Locations of FOXP3+ cells (left) and trajectory of Observed Ripley's K as well as the theoretical and permuted estimates of CSR (right)."----

x[["derived"]][["ripleys"]] <- ripleys_k(mif = x, id = "deidentified_sample", 
                                         mnames = mnames[1], num_permutations = 1,
                                         edge_correction = 'translation', r = seq(0,100,10),
                                         kestimation = TRUE, keep_perm_dis = FALSE, 
                                         mlabels = mlabels[1])

x[["derived"]][["spatial_plotsfoxp3"]] <- plot_immunoflo(x, plot_title = "deidentified_sample", 
                                                    mnames = mnames[1], mlabels = mlabels[1], 
                                                    cell_type = "classifier_label", 
                                                    mcolors = '#522D80')


#For additional plots, uncomment out this section
#library(rlist)
#library(tidyr)
# library(ggplot2)
# library(dplyr)
# library(grid)

# plot_data = x[["derived"]][["ripleys"]] %>%
#   select(deidentified_sample, marker, r_value,
#          `Observed K`, `Theoretical CSR`, `Permuted CSR`) %>%
#   pivot_longer(cols = 4:6, names_to = 'Estimate', values_to = 'K')
# 
# plotlist = lapply(setNames(as.character(unique(x[["derived"]][["ripleys"]]$deidentified_sample)),
#                            as.character(unique(x[["derived"]][["ripleys"]]$deidentified_sample))),
#                   function(a){
#                     left = x[["derived"]][["spatial_plotsfoxp3"]][[a]] +
#                       theme(legend.position = 'none',
#                             title = element_blank(),
#                             axis.text = element_blank(),
#                             axis.ticks = element_blank(),
#                             aspect.ratio = 1)
# 
#                     right = ggplot(data = plot_data %>% filter(deidentified_sample == a),
#                                    aes(x = r_value, K, color = Estimate)) +
#                       geom_line() +
#                       theme_bw()
#                     return(ggarrange(plotlist = list(left, right), align = 'h'))
# })
# 
# 
# 
# ggarrange(plotlist = plotlist, ncol = 1, common.legend = TRUE, legend = 'bottom',
#           labels = names(plotlist), hjust = 0.0)


## ---- bi_ripleys, warning=FALSE, fig.width=13, fig.height=20,fig.align = 'center', fig.cap="Locations of CD8+ and FOXP3+ cells (left) and trajectory of Observed Ripley's K as well as the theoretical and permuted estimates of CSR using FOXP3+ cells as the anchor (right)."----
#Define the makerers on interest
mnames_pairs <- list(list("cd8_opal_520_positive", "foxp3_opal_620_positive"))

mlabels_pairs <- list(list("CD8+", "FOXP3"))

x[["derived"]][["bi-ripleys"]] <- bi_ripleys_k(mif = x, id = "deidentified_sample", 
                                               mnames = mnames_pairs, num_permutations = 1,
                                               r_range = seq(0, 50, 5), edge_correction = "translation", 
                                               kestimation = TRUE, keep_perm_dis = FALSE,
                                               mlabels = mlabels_pairs)


x[["derived"]][["spatial_bivariate"]] <- plot_immunoflo(x, plot_title = "deidentified_sample", 
                                                    mnames = unlist(mnames_pairs[[1]]), mlabels = unlist(mlabels_pairs[[1]]), 
                                                    cell_type = "classifier_label", 
                                                    mcolors = c('#522D80','#F56600'))

#For additonal plots uncomment this section

#library(rlist)
#library(tidyr)
# library(ggplot2)
# library(dplyr)
# library(grid)
# plot_data = x[["derived"]][["bi-ripleys"]] %>% 
#   select(deidentified_sample, anchor_marker, comparison_marker, r_value,
#          `Observed K`, `Theoretical CSR`, `Permuted CSR`) %>%
#   pivot_longer(cols = 5:7, names_to = 'Estimate', values_to = 'K')
# 
# plotlist = lapply(setNames(as.character(unique(x[["derived"]][["bi-ripleys"]]$deidentified_sample)),
#                            as.character(unique(x[["derived"]][["bi-ripleys"]]$deidentified_sample))),
#                   function(a){
#                     left = x[["derived"]][["spatial_bivariate"]][[a]] +
#                       theme(legend.position = 'none',
#                             title = element_blank(),
#                             axis.text = element_blank(),
#                             axis.ticks = element_blank(),
#                             aspect.ratio = 1)
# 
#                     right = ggplot(data = plot_data %>% filter(deidentified_sample == a),
#                                    aes(x = r_value, K, color = Estimate)) +
#                       geom_line() +
#                       theme_bw()
#                     return(ggarrange(plotlist = list(left, right),align = 'h'))
# })
# 
# 
# 
# ggarrange(plotlist = plotlist, ncol = 1, common.legend = TRUE, legend = 'bottom',
#           labels = names(plotlist), hjust = 0)



## ---- hist, warning=FALSE, fig.width=13, fig.height=10,fig.align = 'center', fig.cap="Histogram of individual permutation distribution of observed values and the vertical line corresponds to the theoretical value. Note that for permutation estimate of CSR 100 permutations is adequate for estimation and 1000 may be required for hypothesis testing."----

x[["derived"]][["ripleys"]] <- ripleys_k(mif = x, id = "deidentified_sample", 
                                         mnames = mnames[1], num_permutations = 1,
                                         edge_correction = 'translation', r = seq(0,100,10),
                                         kestimation = TRUE, keep_perm_dis = TRUE, 
                                         mlabels = mlabels[1])

#library(rlist)
#library(tidyr)
# library(ggplot2)
# library(dplyr)
# library(grid)

# hist_data = x[["derived"]][["ripleys"]] %>% 
#   select(deidentified_sample, marker, r_value,`Permuted CSR`) %>%
#   group_by(deidentified_sample, r_value) %>% 
#   summarize(`Permuted Mean` = mean(`Permuted CSR`)) %>%
#   mutate(`Theoretical CSR` = pi*r_value^2) %>%
#   right_join(x[["derived"]][["ripleys"]]) %>%
#   pivot_longer(cols = c(`Observed K`, `Permuted Mean`, `Theoretical CSR`),
#                values_to = 'vert_line', names_to = 'Estimate') %>%
#   filter(r_value >= 50) %>%
#   mutate(r_value = factor(r_value))
#   
# 
# ggplot(data = hist_data, aes(x = `Permuted CSR`, y=..count../3)) + 
#   geom_histogram(col = 'white') + facet_grid(r_value~deidentified_sample, 
#                                              scale = 'fixed') + 
#   geom_vline(aes(xintercept = vert_line, color = Estimate)) +
#   ylab('Count')


