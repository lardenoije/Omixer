# omixr #
### Multivariate and reproducible randomization of omics samples with intuitive visualization of resulting layouts ###
***
![](http://www.molepi.nl/images/logo.png)

Batch effects can have a major impact on the results of omics studies [(Leek et al, 2010)](https://www.nature.com/articles/nrg2825). Randomization is the first, and arguably most influential, step in handling them. However, its implementation suffers from a few key issues:

* A single random draw can inadvertently result in high correlation between technical covariates and biological factors. Particularly in studies with large numbers of batches and outcomes of interest, minimizing these correlations is crucial.
* Long, randomized sample lists are unintuitive and translate poorly into any wet lab that is not fully automated. This can result in errors and sample mixups.
* The randomization process is inherently unclear in many publications, rarely described despite the varying efficacy of methods.
* Randomized layouts are not always reproducible, resulting in inconsistent results.

To combat these problems, we developed **Omixer** - an R package for multivariate randomization and reproducible generation of intuitive sample layouts.

## Workflow ##

**Omixer** randomizes input sample lists multiple times (default: 1,000) and then combines these randomized lists with plate layouts, which can be selected from commonly used setups or custom-made. It can handle paired samples, keeping these adjacent but shuffling their order, and allows explicit masking of specific wells if, for example, plates are shared between different studies.

After performing robust tests of correlation between technical covariates and selected randomization factors, a layout is chosen using these criteria:

* No test provided sufficient evidence to suggest correlation between the variables (all p-values over 0.05).
* From the remaining layouts, return one where the absolute sum of correlations is minimized.

The optimal randomized list can then be processed by `omixerSheet`, returning intuitive sample layouts for the wet lab.
