# covid19-growth

This repository contains the data and the codes to generate all figures in the paper: 

Pellis L, Scarabel F, Stage HB, Overton CE, Chappell LHK, Fearon E, Bennett E, University of Manchester COVID-19 Modelling Group, Lythgoe KA, House TA and Hall I, 
"Challenges in control of COVID-19: short doubling time and long delay to effect of interventions", Philosophical Transactions of the Royal Society B (2021)
medRxiv: https://doi.org/10.1101/2020.04.12.20059972
 
Previous preprint:
31 March 2020: https://arxiv.org/abs/2004.00117
(Version 1.0 of code refers to this preprint)

## Data

Due to backfilling and retrospective updates, inconsistencies are expected when data is downloaded at different times, so we provide for convenience the data used for the analyses in this manuscript.

The data file "WHO_data_18.csv" contains time series of cumulative numbers of daily confirmed cases reported by WHO, restricted to the 18 European countries with more than 1000 cumulative cases by 27 March (the first day cumulative cases in Italy overtook those in China). The data file "Italy.csv" contains the cumulative daily numbers for multiple metrics in Italy collected from daily reports published by the Istituto Superiore di Sanità (ISS: http://www.salute.gov.it/portale/news/p3_2_1.jsp?lingua=italiano&menu=notizie&area=nuovoCoronavirus&notizie.page=0). 

For both data files, under the assumption that data reported by WHO one day are those collected the day before, the dates have been shifted backwards by one day (e.g. the cumulative number of cases published by WHO on 1 March appear in our dataset as referring to 29 February).



The data files "IncubationPeriod.xlsx", "HospitalisationDelay_HongKong.xlsx" and "HospitalisationDelay_Singapore.xlsx" were extracted from:
K.  Sun,  Spreadsheet  of  patient-level  data  until  31.01.2020. https://docs.google.com/spreadsheets/d/1Gb5cyg0fjUtsqh3hl_L-C5A23zIOXmWH5veBklfSHzg/edit?usp=sharing. From K. Sun, J. Chen, C. Viboud, "Early  epidemiological  analysis  of  the  coronavirus  disease  2019  outbreak  based  on  crowdsourced  data: a population-level observational study". The Lancet Digital Health 2, 4 (2020)

In particular:
-- "IncubationPeriod.xlsx" was obtained by filtering on patients who had fixed travel windows in Wuhan prior to symptom onset. These travel dates provided an exposure window. The columns are exposure start, exposure end and onset date. 
-- "HospitalisationDelay_HongKong.xlsx" was obtained by filtering on patients in Hong Kong with known symptom onset dates and hospitalisation dates. The columns are onset date and hospitalisation date. 
-- "HospitalisationDelay_Singapore.xlsx" was obtained by filtering on patients in Singapore with known symptom onset dates and hospitalisation dates. The columns are onset date and hospitalisation date. 

In all three files, the data are processed to remove inconsistent dates and replace "0" durations with "0.5".



The file "Literature_data.xlsx" contains the filtered list of the literature search for unconstrained (or early) doubling times and basic (or early) reproduction numbers from the first half of 2020. Publications have been extracted from Google Scholar between 31 March and 1 April 2021, with searches [(doubling time’ OR ‘growth rate) AND (covid OR SARS-CoV-2)] and [‘basic reproduction number’ AND (covid OR SARS-CoV-2)]. Google Scholar citations counts were recorded between 31 March and 1 April 2021. Publications with fewer than 100 citations are not shown, with the exception of the estimates presented in this paper "UoM_estimates.xlsx", which were given a default citation count of 150 for visibility in the plot.


## Codes

The codes are available under the MIT licence, see ‘LICENCE’ for details.
Please cite this manuscript when using the codes. 
