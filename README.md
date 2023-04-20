---
title: "ReadMe"
author: "Steffen Funk"
date: '2023-04-20'
output: html_document
---

## Background:
Here you get access to an individual-based bioenergetic growth model which is set up to model the growth of cod in ages 2 to 4 (*Gadus morhua* L.) in the Belt Sea (western part of the Western Baltic Sea) on a daily basis within one year. The model incorporates contemporary in-situ process knowledge on food intake and seasonal- and temperate-related spatial distribution of cod and allows to identify seasonal growth patterns. For a detailed description of the model (i.e., explanation of equations, submodels used etc.) we refer to the corresponding article by **Funk et al. (currently under review)**.

|![alt text](https://github.com/FunkSteffen/WBC_bioenergetic_growth_model/blob/main/figure1_proto.png?raw=true)|
|:--:| 
|**Set-up of the individual-based bioenergetic growth model.** a) Schematic representation of important model components and their influence on growth: (1) temperature conditions (SST and TDiff, see Material and methods) (2) temperature at residence depth, (3) food availability and intake, (4) physiological processes, and (5) growth in length. b) Functional pathways (arrows) within the bioenergetic model: rectangles denote temperature variables (SST and TDiff), diamonds denote estimates derived from sub-models I-III; grey ellipses indicate estimates based on published physiological functions; hexagons display cod individual length and weight. i – time step, j – size, k – diet cluster.|

## Information on the model script:
The model script contains several functions which are defined in the beginning of the script and corresponding examples which explain how to use the above defined functions.
The bioenergetic growth model relies on several submodels which are imported as rds data files. Furthermore, some data frames (i.e., cod SMALK-data and temperature predictions from the BSIOM) are needed to run the model which can be downloaded as as csv data files and which will be imported when you run the bioenergetic model code.
Note that the bioenergetic growth model for cod in the Western Baltic model was intentionally set up and tested for the modelling periods 2016-2017 and 2016-2017 (since the are the years, where the stomach data was collected). However, we here provide also temperature data for additional years to run the model for other time periods. 


## Before you are getting started:
To run the model script (**WBC_bioenergetic_growth_model_v1_0.R**), you need to download the following data files:

**temp_at_3m_depth_strata_1979_2018.csv** (contains predicted temperature values from the BSIOM for 3m depth strata for the period between 1979 and 2018 from the Belt Sea)

**SST_SBT_Tdiff_1978_2018.csv** (contains predicted sea surface temperature data, sea bottom temperature data and data for the proxy of stratification from BSIOM for the period between 1978 to 2018 from the Belt Sea)

**smalk_raw_data_males.csv** (contains a dummy data set on biological parameters for cod in the Belt Sea and Sound based on cod biological data recorded during the Baltic international trawl survey which can be downloaded from ICES database DATRAS (ICES BITS)).

**energy_dens_and_pi_table.csv** (contains a table listing prey specific gastric evacuation coefficients (pi) and energy densities for different prey types (i.e., diet clusters). Values are used for calculation of daily consumption and maintenance ration)

**observer_model.rds** (contains a fitted linear regression model which is used to calculate the residence depth of a cod. The model is based on a model which explains fishing depth of commercial gill net fishers recorded by at-sea observers by sea surface temperature, a proxy for stratification and the gill net mesh size used. For further details on the model please see also Funk et al., 2021b).

**stomach_model.rds** (contains a fitted generalized additive model which explains stomach content weight of cod with the explanatory variables cod length, temperature at catch depth and catch depth. The model is used to predict the daily stomach contents of cod for a given residence at at a given day. For further details on the model please see also Funk et al., 2021a).

**multinom.rds** (contains a multinomial logistic regression model which explains the diet composition of a cod with the explanatory variables cod length, quarter and catch depth. The model is used to predict the diet composition of cod for a given day. For further details on the model please see also Funk et al., 2021a)


## Contact information:
For further information, support or help please contact the first author of the study S. Funk (Email: steffen.funk@uni-hamburg.de).

## References:
**Funk, S., Frelat, R., Möllmann, C., Temming, A., and Krumme, U. 2021a.** The forgotten feeding ground: patterns in seasonal and depth-specific food intake of adult cod *Gadus morhua* in the western Baltic Sea. *Journal of Fish Biology*, 98(3): 707-722. https://doi.org/10.1111/jfb.14615

**Funk, S., Krumme, U., Temming, A., and Möllmann, C. 2021b.** Gillnet fishers‘ knowledge reveals seasonality in depth and habitat use of cod (*Gadus morhua*) in the Western Baltic Sea. *ICES Journal of Marine Science*, 77(5):1816-1829.

**ICES**, ICES BITS (Baltic International Trawl Survey) dataset. ICES, Copenhagen.





