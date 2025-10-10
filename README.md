
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

This repository contains all the scripts required to reproduce results
in the manuscript **Species overlap and phylogenetic relatedness result
in community statistical non-independence (and what to do about it)**.
If you are interested in using EcoCoMix, please visit the
[webiste](https://tpaknok.github.io/EcoCoMix/articles/Empirical_single.html)

| Folder name | Description                                                                                        |
|-------------|----------------------------------------------------------------------------------------------------|
| Data        | Data needed for the empirical analyses and obtaining summary statistics of the simulation analyses |
| Code        | Code needed for the empirical analyses and simulation.                                             |
| Figure      | A folder containing all result figures.                                                            |
| Table       | A folder containing TableS2 (pre-formatting) in the manuscript.                                    |

# Code

Two scripts are provided, one for the simulation and another for the
empirical analyses.

| Scripts             | Description                                                                                     |
|---------------------|-------------------------------------------------------------------------------------------------|
| updated_sim_spaMM.R | The script for the simulation analyses in the main text.                                        |
| sim_supp.R          | The script for the simulation analyses in the supporting information (High diversity scenario). |
| KSR_empirical.R     | The script for the empirical analyses.                                                          |

# Data - Overview

Please install the [EcoCoMix
package](https://github.com/tpaknok/EcoCoMix), which contains the
necessary functions and also data for the empirical analyses. The data
have also been provided in the Data folder.

| Data          | Description                                                                                      |
|---------------|--------------------------------------------------------------------------------------------------|
| KSR.csv       | Species compositional data across the plots.                                                     |
| KSR_EF.csv    | Ecosystem functions measured in each plot.                                                       |
| KSR_MLtree    | Phylogeny tree for the 14 species.                                                               |
| sim500.Rdata  | The results of the simulation described in the main text.                                        |
| simSupp.Rdata | The results of the simulation described in the supporting information (High diversity scenario). |

## Data - Species composition

KSR.csv (or run data(KSR) in R) describes the species composition of
each plot (Row: 88 sites, column: 14 species).

| Species                         | Abbreviation |
|---------------------------------|--------------|
| <i>Andropogon gerardii</i>      | ANGE         |
| <i>Asclepias tuberosa</i>       | ASTU         |
| <i>Elymus canadensis</i>        | ELCA         |
| <i>Elymus trachycaulus</i>      | ELTR         |
| <i>Desmodium canadense</i>      | DECA         |
| <i>Lespedeza capitata</i>       | LECA         |
| <i>Monarda fistulosa</i>        | MOFI         |
| <i>Penstemon hirsutus</i>       | PEHI         |
| <i>Pycnanthemum tenuifolium</i> | PYTE         |
| <i>Pycnanthemum virginianum</i> | PYVI         |
| <i>Rudbeckia hirta</i>          | RUHI         |
| <i>Schizachyrium scoparium</i>  | SCSC         |
| <i>Solidago altissima</i>       | SOAL         |
| <i>Solidago nemoralis</i>       | SONE         |

## Data - Phylogeny

KSR_MLtree (or data(KSR_MLtree)) is a phylogeny tree containing the 14
species used in the experiment.

## Data - Ecosystem function

KSR_EF.csv (or data(KSR_EF)) describes the ten ecosystem functions
measured in the experiment.

| Variable         | Description                                                      |
|------------------|------------------------------------------------------------------|
| Plot             | Plot identity                                                    |
| Real.rich        | Number of species planted                                        |
| litter2012       | Amount of litter measured in 2012                                |
| ave.biomass      | Average biomass across 2012-2014                                 |
| LAI              | Leaf area index, a simplified dimension of structural complexity |
| mean.N.change    | delta 15N change averaged across surface and deep soil           |
| poll_total       | Total number of pollinators                                      |
| flwr_total       | Total number of flowers                                          |
| Mass.loss.2month | Decomposition after 2 months                                     |
| Damage_effect    | Damage reduction effect                                          |
| bugs             | Total number of arthropods                                       |
| bug.rich         | Species richness of arthropods                                   |

## Simulation

For the simulation, you can simply run the simulation scripts. However,
this will take a long time. It took me \>24 hours to finish the
simulation with my laptop (AMD Ryzen 7 6800H 3.20GHz, 16GB RAM).

IF you don’t want ot wait, you can load the results (sim500.RData or
simSupp.Rdata) and obtain summary statistics.

# Figure

The folder contains all result figures included in the manuscript:
p_coef(Figure S1), p_KSR (Figure 4), p_KSR_all(Figure S2), p_sim(Figure
3).

# Table

The folder contains Table S2 (before formatting). Note that this is NOT
the data.

# Contact

Toby Tsang (<tpaknok@gmail.com>)
