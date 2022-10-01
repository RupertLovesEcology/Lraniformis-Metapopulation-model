# A stochastic hydroecological metapopulation model for <i>Litoria raniformis</i> which explores the role of environmental water provision in the probability of persistence over a 60 year time step
These files underpin the stochastic hydroecological metapopulation model presented in the article "Stochastic metapopulation dynamics of a threatened amphibian to guide optimal water delivery"

<strong>AUTHORS</strong>: Rupert Mathwin, Matt S. Gibbs and Corey J. A. Bradshaw

<strong>CONTACT</strong>: rupert.mathwin.ecology@gmail.com

<strong>URL</strong>: http://GlobalEcologyFlinders.com

<strong>INSTITUTION</strong>: Flinders University

<strong>INSTITUTION</strong>: Rupert.Mathwin.Ecology

<strong>RELEASE DATE</strong>: May 2022

R code accompanies article: 

<strong>Mathwin, R, Wassens, S, Gibbs, MS, Young, J, Ye, Q, Saltré, F,</strong> and <strong>Bradshaw, CJA</strong> Stochastic metapopulation dynamics of a threatened amphibian to guide optimal water delivery <i>in review</i>

<strong>AIM</strong>: This stochastic hydroecological metapopulation model examines the population viability of <i>Litoria raniformis</i> across a metapopulation of 23 wetlands between Locks 3 and 2 in southern Australia's Murray-Darling Basin. We test five competing conservation treatments which manipulate hydroperiod. We then reexamine these treatments under catastrophic drought conditions.

Repository includes the following files:
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/Lraniformis_Metapopulation_Model_V12.R">Lraniformis_Metapopulation_Model_V12.R</a>' — R code to run the metapopulation model. 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/Lraniformis_Metapopulation_globalSensitivity_V12.R">Lraniformis_Metapopulation_globalSensitivity_V12.R</a>' — R code to run the global sensitivity analysis on the metapopulation model (uses latin hypercube resampling to populate a boosted regression tree). 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/orderedAnnualWetness.csv">orderedAnnualWetness.csv</a>' — .csv containing the summed spring and winter river heights as a proxy for wetness. Determines if the wetland remains full between years (and accumulates aquatic predators).
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/orderedAnnualWetnessDrought2.csv">orderedAnnualWetnessDrought2.csv</a>' — .csv containing the summed spring and winter river heights as a proxy for wetness. Determines if the wetland remains full between years (and accumulates aquatic predators). This treatment includes a second severe drought.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/orderedHydroDrought2.csv">orderedHydroDrought2.csv</a>' — .csv containing the 10th highest daily river height which determines wetland filling throughout the reach. This treatment includes a second severe drought.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/orderedHydroForecast.csv">orderedHydroForecast.csv.csv</a>' — .csv containing the 10th highest daily river height which determines wetland filling throughout the reach. 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/PostDrought.csv">PostDrought.csv</a>' — .csv containing the 10000 starting scenarios for the reach (23 populations with age demography). 
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/matrixOperators.r">matrixOperators.R</a>' — functions to manipulate matrix models
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/movementType.csv">movementType.csv</a>' — .csv which has classified every possible single dispersal between 2 wetlands into one of five journey types. This is used to assign landscape resistance. 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/stayWet.csv">stayWet.csv</a>' — .csv containing the wetness threshold that each wetland requires to remian wet between years (verified from historic satelite data). 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetlandMetadataV2.csv">wetlandMetadataV2.csv</a>' — column 1 is the size category (which determines the local popualtion capacity) and sill height is the river height required to start filling the wetland
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetlandMovement.csv">wetlandMovement.csv</a>' — .csv which measures every possible single dispersal between 2 wetlands following several movement rules (e.g. can only cross the river once). 
