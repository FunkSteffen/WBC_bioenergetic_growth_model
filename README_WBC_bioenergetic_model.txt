READ ME

Here you get access to a bioenergetic growth model which is set up to model the growth of cod (Gadus morhua L.) in the Belt Sea (western part of the Western Baltic Sea).
For details on the model (i.e., explanation of equations, submodels used etc.) we refer to the corresponding article by Funk et al.

Information on the model script

The model script contains on both functions which are defined in the beginning of the script and examples which explain how to use the above defined functions (i.e., get the bioenergetic growth model running).
The model relies on some submodels which are imported as rds data files. Furthermore, some data frames (i.e., cod SMALK-data and temperature predictions from the BSIOM) are needed and are imported as csv data files.
The bioenergetic growth model is intentionally set up to calculate growth of cod in 2016 and 2017. However, we provide also temperature for other years. Thus, the model can be also run for other years.

To run the model script (WBC_bioenergetic_growth_model_v1_0.R), you need to download the following data files:

temp_at_3_depth_strata_1979_2018.csv (contains predicted temperature values from the BSIOM for 3m depth strata for the period between 1979 and 2018 from the Belt Sea)

SST_SBT_Tdiff_1978_2018.csv (contains predicted sea surface temperature data, sea bottom temperature data and data for the proxy of stratification from BSIOM for the period between 1978 to 2018 from the Belt Sea)

smalk_raw_data_males.csv (contains a dummy data set on biological parameters for cod in the Belt Sea and Sound based on cod biological data recorded during the Baltic international trawl survey which can be downloaded from ICES database DATRAS (ICES BITS)).

energy_dens_and_pi_table.csv (contains a table listing prey specific gastric evacuation coefficients (pi) and energy densities for different prey types (i.e., diet clusters). Values are used for calculation of daily consumption and maintenance ration)

observer_model.rds (contains a fitted linear regression model which is used to calculate the residence depth of a cod. The model is based on a model which explains fishing depth of commercial gill net fishers recorded by at-sea observers by sea surface temperature, a proxy for stratification and the gill net mesh size used. For further details on the model please see also Funk et al., 2020).

stomach_model_incl_small.rds (contains a fitted generalized additive model which explains stomach content weight of cod with the explanatory variables cod length, temperature at catch depth and catch depth. The model is used to predict the daily stomach contents of cod for a given residence at at a given day. For further details on the model please see also Funk et al., 2021).

multinom.rds (contains a multinomial logistic regression model which explains the diet composition of a cod with the explanatory variables cod length, quarter and catch depth. The model is used to predict the diet composition of cod for a given day. For further details on the model please see also Funk et al., 2021)



Contact information
For further information, support or help please contact the first author (Steffen Funk, Email: steffen.funk@uni-hamburg.de).


	
References
Funk, S., Frelat, R., Möllmann, C., Temming, A., and Krumme, U. 2020a. The forgotten feeding ground: patterns in seasonal and depth-specific food intake of adult cod Gadus morhua in the western Baltic Sea. Journal of Fish Biology, 98(3): 707-722. https://doi.org/10.1111/jfb.14615

Funk, S., Krumme, U., Temming, A., and Möllmann, C. 2020b. Gillnet fishers‘ knowledge reveals seasonality in depth and habitat use of cod (Gadus morhua) in the Western Baltic Sea. ICES Journal of Marine Science, 77(5):1816-1829.

ICES BITS (Baltic International Trawl Survey) dataset. ICES, Copenhagen.
