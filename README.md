# A stochastic hydroecological metapopulation model for <i>Litoria raniformis</i> which explores the role of environmental water provision in the probability of persistence over a 60 year time step
These files underpin the stochastic hydroecological metapopulation model presented in the article "Stochastic metapopulation dynamics of a threatened amphibian to guide optimal water delivery"

<strong>AUTHORS</strong>: Rupert Mathwin, Matt S. Gibbs and Corey J. A. Bradshaw

<strong>CONTACT</strong>: rupert.mathwin.ecology@gmail.com

<strong>URL</strong>: http://GlobalEcologyFlinders.com

<strong>INSTITUTION</strong>: Flinders University

<strong>INSTITUTION</strong>: Rupert.Mathwin.Ecology

<strong>RELEASE DATE</strong>: May 2022

R code accompanies article: 

<strong>Mathwin, R, Wassens, S, Gibbs, MS, Young, J, Ye, Q, Saltré, F,</strong> and <strong>Bradshaw, CJA</strong> Stochastic metapopulation dynamics of a threatened amphibian to guide optimal water delivery <i>under development</i>

<strong>AIM</strong>: This stochastic hydroecological metapopulation model examines the population viability of <i>Litoria raniformis</i> across a metapopulation of 23 wetlands between Locks 3 and 2 in southern Australia's Murray-Darling Basin. We test five competing conservation treatments which manipulate hydroperiod. We then reexamine these treatments under catastrophic drought conditions.

Repository includes the following files:
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/Lraniformis_Metapopulation_Model_V12.R">Lraniformis_Metapopulation_Model_V12.R</a>' — R code to run the metapopulation model. 
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/Lraniformis_Metapopulation_globalSensitivity_V12.R">Lraniformis_Metapopulation_globalSensitivity_V12.R</a>' — R code to run the metapopulation model. 
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/matrixOperators.r">matrixOperators.R</a>' — functions to manipulate matrix models

- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/CompHeightL3.csv">CompHeightL3.csv</a>' — wetland inundation records used to construct Markov-chains for stochastic resampling of inundation regimes
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetlandMetadataV2.csv">wetlandMetadataV2.csv</a>' — column 1 is the size category (which determines the local popualtion capacity) and sill height is the river height required to start filling the wetland
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/AnnualTransitionsofExistingSills.R">AnnualTransitionsofExistingSills.R</a>' — c
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/CalculateAnnualSillsfromRaw.R">CalculateAnnualSillsfromRaw.R</a>' — c
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/burntinNmat.csv">burntinNmat.csv</a>' — c
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/timeSeriesForecastV4.R">timeSeriesForecastV4.R</a>' — c
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetlandMovement.csv">wetlandMovement.csv</a>' — c

- TO BE UPDATED ONCE FINALISED
