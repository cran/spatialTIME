% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bi_nn_g2.R
\name{bi_NN_G}
\alias{bi_NN_G}
\title{Bivariate Nearest Neighbor G(r)}
\usage{
bi_NN_G(
  mif,
  mnames,
  r_range = 0:100,
  num_permutations = 50,
  edge_correction = "rs",
  keep_perm_dis = FALSE,
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL
)
}
\arguments{
\item{mif}{object of class `mif` created by function `create_mif()`}

\item{mnames}{character vector of column names within the spatial files, indicating whether a cell row is positive for a phenotype}

\item{r_range}{numeric vector of radii around marker positive cells which to use for G(r)}

\item{num_permutations}{integer number of permutations to use for estimating core specific complete spatial randomness (CSR)}

\item{edge_correction}{character vector of edge correction methods to use: "rs", "km" or "han"}

\item{keep_perm_dis}{boolean for whether to summarise permutations to a single value or maintain each permutations result}

\item{workers}{integer number for the number of CPU cores to use in parallel to calculate all samples/markers}

\item{overwrite}{boolean whether to overwrite previous run of NN G(r) or increment "RUN" and maintain  previous measurements}

\item{xloc, yloc}{the x and y location columns in the spatial files that indicate the center of the respective cells}
}
\value{
object of class `mif` containing a new slot under `derived` got nearest neighbor distances
}
\description{
Bivariate Nearest Neighbor G(r)
}
\examples{
x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical \%>\% 
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  sample_data = spatialTIME::example_summary \%>\% 
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  spatial_list = spatialTIME::example_spatial[1:2],
  patient_id = "deidentified_id", 
  sample_id = "deidentified_sample")
    
mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
  "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
  "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
\dontrun{
x2 = bi_NN_G(mif = x, mnames = mnames_good[1:2], 
      r_range = 0:100, num_permutations = 10, 
      edge_correction = "rs", keep_perm_dis = FALSE, 
      workers = 1, overwrite = TRUE)
}

}
