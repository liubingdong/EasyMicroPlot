\name{hello}
\alias{hello}
\title{beta_plot}
\usage{
hello()
}
\description{
Prints 'A easy R script plot for Microbiome analysis'.
}
\arguments{
  \item{dir}{An directory contain microbial relative abundance at different levels.}
  \item{group_level}{Set the group order in Mapping File [Not necessary]}
  \item{min_relative}{Set the min relative abundance for filtering}
  \item{min_ratio}{Set the min ratio presence for taxonomy}
  \item{design}{Mapping File}
  \item{adjust}{Set addtional dummy tax in empty data when filtering [Default:False]}
  \item{pattern}{Set key string character for relative abundance files}
  \item{palette}{Colour palette}
  \item{method}{Set method in multiple comparison [Tukey(Default),Dunnett,Sequen,AVE,Changepoint,Williams,Marcus,McDermott,UmbrellaWilliams,GrandMean]}
  \item{distance}{bray(Default),manhattan, euclidean, canberra, clark, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao or mahalanobis}
}
\examples{
beta_result=beta_plot(group_level=c('Ctrl','CtrlD', 'HFD', 'HFDDH'),
                      dir = '.',min_relative = 0.001,min_ratio = 0.85,
                      design = 'mapping.txt',adjust = F,pattern = 'L')
}
