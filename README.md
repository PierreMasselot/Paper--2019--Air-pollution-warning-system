# Paper--2019--Air-pollution-warning-system
***
R code implementing the air pollution-health warning system of the province of Quebec proposed in the paper:

> Masselot, P., Chebana, F., Lavigne, É., Campagna, C., Gosselin, P., Ouarda, T.B.M.J., 2019. Toward an improved air pollution warning system in Quebec. *International Journal of Environmental Research and Public Health.*

This repository contains the following files:
* Main.R: Main code executing all four steps of the methodology as well as the sensitivity analysis. The code also produce the figures found in the article;
* Other functions.R: Supplementary custom functions used to produce figures;
* CMMdata.csv: Data of the metropolitan area of Montreal used to produce the paper's results;
* CMQdata.csv: Data of the metropolitan area of Quebec City used to produce the paper's results.

In addition, this code uses functions from the package `hhws` which can be installed from R using the following command (necessitates the `devtools` package):
```R 
devtools::install_github("PierreMasselot/hhws")
```
