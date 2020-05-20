# omixr #
### Multivariate randomization and reproducible generation of sample layouts to optimally combat batch effects in -omics data ###
***
![alt text](http://www.molepi.nl/images/logo.png)

Batch effects can have a major impact on the results of omics studies [(Leek et al, 2010)](https://www.nature.com/articles/nrg2825). Randomization is the first, and arguably most influential, step in handling them. However, its implementation suffers from a few key issues:

* A single random draw can inadvertently result in high correlation between technical covariates and biological factors. Particularly in studies with high numbers of batches and outcomes of interest, it is not trivial to minimize these correlations.
* Unless the lab is fully automated, long, randomized sample lists translate poorly into the wet lab. This monotonous presentation of sample layouts is unintuitive and this can result in errors and mixups.
* The randomization process is inherently unclear in many publications, and may even not be reproducible depending on the code used. Some studies may have more or less effective methods of randomization, and this variability is rarely described in papers.

To combat these problems, we developed **omixr** - an R package for multivate randomization and reproducible generation of sample layouts.

## The Workflow ##

Using a sample list input, omixr will generate the specified number of randomized sample lists (default: 10,000). It can handle paired samples, keeping these adjacent but shuffling their order.

These lists are then combined with plate layouts, which can be selected from commonly used setups or custom-made. Explicitly masking wells is possible if, for example, you are sharing plates with another study.

After calculating correlations between defined technical covariates and selected randomization factors, a layout is chosen that using the following criteria:
* No test provided sufficient evidence to suggest correlation between the variables (all p-values over 0.05)
* From the remaining layouts, return one where the absolute sum of correlations is minimized

The resulting sample layout can then be printed as a list for automated setups, or processed by `omixr_sheet` which returns easy-to-read visual plate layouts for the wet lab.

## Installation ##

The **omixr** package can be installed using [**devtools**](https://github.com/hadley/devtools) in R.

```{r devtools, eval=FALSE}
library(devtools)
install_github("molepi/omixr")
```    
