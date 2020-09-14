This code is associated with the paper from Shocket et al., "Transmission of West Nile and five other temperate mosquito-borne viruses peaks at temperatures between 23˚C and 26˚C". eLife, 2020. http://dx.doi.org/10.7554/eLife.58511


# Six-Viruses-Temp

This repo contains all of the data, code, and methods for the article: 
Shocket MS, Verwillow AB, Numazu MG, Slamani H, Cohen JM, El Moustaid F, Rohr J, Johnson LR, and Mordecai EA. 2020. Transmission of West Nile and five other temperate mosquito-borne viruses peaks at temperature between 23C and 26C. eLife.

Corresponding author: Marta Shocket (marta.shocket@gmail.com)

The study builds mechanistic, trait-based R0 models for transmission of 6 mosquito-borne viruses (West Nile virus [WNV], Eastern and Western Equine Encephalitis viruses [EEEV and WEEV], St. Louis Encephalitis virus [SLEV], Sindbis virus [SINV], and Rift Valley Fever virus [RVFV]) that are transmitted by many vectors (including Culex pipiens, Cx. quinquefasciatus, Cx. tarsalis, and others), and performs model validations with human case data for WNV, SLEV, and EEEV.

METADATA FOR THE DATA FILES

TRAIT DATA (IN TRAIT FITS DIRECTORY) - Figures 3, 4, 5, 6, and 7

The trait data are in 7 .csv files with file names "TraitData_x" where x is the trait name or relevant parameter in the model. In each file, the column names and definitions are: 
  Ref = series ID with the name of the person who digitized the data from the published paper
  Trait Name = trait abbreviation, which is often the model parameter (see list below for details)
  T = temperature (C)
  trait = trait value
  ErrorPosSI & ErrorBegSI = standard errors for the trait value
  Trait2.name = second variable in the experiment, if needed; for vector competence this is usually time
  host.code = four letter code for mosquito vector species (see list below for details)
  paras.code = four letter code for parasite (virus)
  Citation = first author, year of publication, and journal name for publication the data came from
  Figure = the figure or table that the data came from in the original publication
  Notes = any notes
  joint.code = eight letter code for combination of mosqutio vector species and virus (host.code + paras.code, for infection traits only)

The Trait Names are as follows: 
  a = biting rate (1/day; the inverse of gonotropic cycle duration GCD [a = 1/GCD] - the data are recorded as GCD).
  bc = vector competence (propotion; bc = full vector competence = # transmitting / # exposed; c = infection efficiency = # infected / # exposed; b = transmission efficiency = # transmitting / # infected)
  EFD = fecundity traits (model parameter = eggs per female per day; here traits are EFOC/EFGC = eggs per female per oviposition/gonotrophic cycle; ER/EPR = eggs per raft; pO = proportion of females ovipositing)
  EV = egg viability = proportion of eggs hatching
  lf = lifespan (days; the inverse of mortality rate, mu (lf = 1/mu)
  MDR = mosquito development rate (1/day)
  PDR = pathogen/parasite development rate (1/day; the inverse of extrinsic incubation period [EIP = 1/PDR])
  pLA = proportion of larvae surviving to adulthood
  
The host.codes are as follows:
  Cpip = Cx. pipiens
  Cqui = Cx. quinquefasciatus
  Ctar = Cx. tarsalis
  Cuni = Cx. univittatus
  Cthe = Cx. theileri 
  Atri = Ae. triseriatus
  Atae = Ae. taeniorhynchus
  Avex = Ae. vexans
  Cmel = Cs. melanura
  Cmol = Cx. pipiens molestus
  Cpal = Cx. pipiens pallens
  Cres = Cx. restuans
  Ador = Ae. dorsalis
  Anig = Ae. nigromaculis
  Asol = Ae. sollicitans
  Asal = Ae. salinarius
  
VALIDATION DATA (IN WNV CASE DATA DIRECTORY) - Figures 8 and 9

There are two .csv files with data for WNV mosquito vector species in each state, and the main .csv file with West Nile virus case data. The column names and defitions for each file are as follows:

WNVPosMosqbyStateYrSpp_fromPaull.csv
  Year = year
  State = two letter US state abbreviation
  PIP = # of Cx. pipiens testing positive for WNV
  QUI = # of Cx. quinquefasciatus testing positive for WNV
  TAR = # of Cx. tarsalis testing positive for WNV
  PIPPerc = proportion of WNV + mosquitoes that were Cx. pipiens
  QUIPerc = proportion of WNV + mosquitoes that were Cx. quinquefasciatus 
  TARPerc = proportion of WNV + mosquitoes that were Cx. tarsalis
  
WNVPosMosqStateAvgs.csv (avg of all years together from data above, plus values for 3 missing states interpolated from adjacent states)
  State = two letter US state abbreviation
  PIPPerc = proportion of WNV + mosquitoes that were Cx. pipiens
  QUIPerc = proportion of WNV + mosquitoes that were Cx. quinquefasciatus
  TARPerc = proportion of WNV + mosquitoes that were Cx. tarsalis
  
countiesnew_monthlytemps.csv
  County = US county name
  State = US state name
  Abbr = two letter US state abbreviation
  CountyStat = string with two letter US state abbreviation and US county name concatenated together
  StateID = unique number for each US state
  CountyID = unique number for each US county
  StateCount = unique number from joining StateID and CountyID 
  x = latitude coordinate for centroid of county
  y = longitude coordinate for centroid of county
  optional = unclear what this means, value is TRUE for all rows
  wnv2001 = number of cases of neuroinvasive WNV cases in 2001
  ...
  wnv2016 = number of cases of neuroinvasive WNV cases in 2016
  sum = sum of neuroinvasive WNV cases from 2001-2016
  pop = county population according to 2010 census
  adj2001 = incidence per 1000 people of neuroinvasive WNV in 2001
  ...
  adj2016 = incidence per 1000 people of neuroinvasive WNV in 2016
  adjsum = sum of incidence per 1000 people of neuroinvasive from 2001-2016
  Arrival = year that WNV arrived in the county
  Timeadjsum = average incidence per 1000 people of neuroinvasive WNV, calculated from year of first arrival through 2016
  Occurence = binary 0/1 = yes/no for if neuroinvasive WNV cases were observed in the county
  yearspresent = number of years that neuroinvasive WNV cases were observed in the county
  peakyear = highest incidence per 1000 people recorded during a single year for the county
  tmp = average summer (May - Sept) temperature (celcius) from 2001-2016, calculated from monthly mean temperatures
  pre = average precipitation from 2001-2016 (mm per year)
  area = county area (meters^2)
  tmp01 = monthly mean temperature for January 2001-2016
  ...
  tmp12 = monthly mean temperature for December 2001-2016
  PIPPercAvg = proportion of WNV + mosquitoes that were Cx. pipiens (from state level data)
  QUIPercAvg = proportion of WNV + mosquitoes that were Cx. quinquefasciatus (from state level data)
  TARPercAvg = proportion of WNV + mosquitoes that were Cx. tarsalis (from state level data)
  WNV.R0.01 = predicted temperature-dependent R0 for WNV in January based average January temperature, the % of each vector species, and the 3 species specific WNV R0 models from the study
  ...
  WNV.R0.12 = predicted temperature-dependent R0 for WNV in December based average December temperature, the % of each vector species, and the 3 species specific WNV R0 models from the study
  SLEV.R0.01 = predicted temperature-dependent R0 for SLEV in January based average January temperature and the SLEV R0 model from the study
  ...
  SLEV.R0.12 = predicted temperature-dependent R0 for SLEV in December based average December temperature and the SLEV R0 model from the study
  EEEV.R0.01 = predicted temperature-dependent R0 for EEEV in January based average January temperature and the SLEV R0 model from the study
  ...
  EEEV.R0.12 = predicted temperature-dependent R0 for EEEV in December based average December temperature and the SLEV R0 model from the study
  pop.weight = the proportion of total US population living in that county (calculated from pop value and the sum of all pop values in the spreadsheet)


Materials and methods

All analyses were conducted using R 3.1.3 (R Development Core Team, 2016).

Temperature-dependent trait data

Fitting thermal responses
We fit trait thermal responses with a Bayesian approach using the ‘r2jags’ package (Y-S and Yajima, 2009), an R interface for the popular JAGS program (Plummer, 2003) for the analysis of Bayesian graphical models using Gibbs sampling. It is a (near) clone of BUGS (Bayesian inference Using Gibbs Sampling) (Spiegelhalter et al., 2003). In JAGS, samples from a target distribution are obtained via Markov Chain Monte Carlo (MCMC). More specifically, JAGS uses a Metropolis-within-Gibbs approach, with an Adaptive Rejection Metropolis sampler used at each Gibbs step (for more infor- mation on MCMC algorithms see Gilks et al., 1998).

For each thermal response being fit to trait data, we visually identified the most appropriate functional form (quadratic, Briére, or linear; Equations 3–5) for that specific trait–species combination (Mordecai et al., 2019). For traits with ambiguous functional responses, we fit the quadratic and Briere and used the deviance information criterion (DIC) (Spiegelhalter et al., 2002) to pick the best fit. We assumed normal likelihood distributions with temperature-dependent mean values described by the appropriate function (Equations 3–5) and a constant standard deviation (s) described by an additional fitted parameter (t = 1/s2). The 95% credible intervals in Figures 3–6 estimate the uncer- tainty in the mean thermal response.

We set all thermal response functions to zero when T < Tmin and T > Tmax (for Equation 3 and 4) or when T > -z/m (Equation 5) to prevent trait values from becoming negative. For traits that were proportions or probabilities, we also limited the thermal response functions at 1. For the linear ther- mal responses, we calculated the predicted thermal response in a piecewise manner in order to be conservative: for temperatures at or above the coldest observed data point, we used the trait values predicted by the fitted thermal response (i.e. the typical method); for temperatures below the cold- est observed data point, we substituted the trait estimate at the coldest observed data point (i.e. forcing the thermal response to plateau, rather than continue increasing beyond the range of observed data).

For the fitting process, we ran three concurrent MCMC chains for 25,000 iterations each, discard- ing the first 5000 iterations for burn-in (convergence was checked visually). We thinned the resultant chains, saving every eighth step. These settings resulted in 7500 samples in the full posterior distri- bution that we kept for further analysis.

Generating priors
We used data-informed priors to decrease the uncertainty in our estimated thermal responses and constrain the fitted thermal responses to be biologically plausible, particularly when data were sparse. These priors used our total dataset, which contained temperature-dependent trait data for all of the main species in the analysis (but with the focal species removed, see below), as well as from additional temperate Aedes and Culex species (Buth et al., 1990; Ciota et al., 2014; Kiarie- Makara et al., 2015; Madder et al., 1983; McHaffey, 1972b; McHaffey and Harwood, 1970; Mogi, 1992; Muturi et al., 2011; Oda et al., 1999; Oda et al., 1980; Olejnı ́cek and Gelbic, 2000; Parker, 1982).

We fit each thermal response with a sequential two-step process, where both steps employed the same general fitting method (described above in Fitting Thermal Responses) but used different pri- ors and data. In step 1, we generated high-information priors by fitting a thermal response to data from all species except the focal species of interest (i.e. a ‘leave-one-out’ approach). For example, for the prior for biting rate for Cx. pipiens, we used the biting rate data for all species except Cx. pipiens. For this step, we set general, low-information priors that represented minimal biological constrains on these functions (e.g. typically mosquitoes die if temperatures exceed 45 ̊C, so all bio- logical processes are expected to cease; Tmin must be less than Tmax). The bounds of these uniformly distributed priors were: 0 < Tmin < 24, 26 < Tmax < 45 (quadratic) or 28 < Tmax < 45 (Brie ́re), 0 < q < 1,–10 < m < 10, and 0 < b < 250. Then in step 2, we fit a thermal response to data from the focal species using the high-information priors from step 1.

Because we cannot directly pass posterior samples from JAGS as a prior, we modified the results from step 1 to use them in step 2. We used the ‘MASS’ package (Venables and Ripley, 2002) to fit a gamma probability distribution to the posterior distributions for each thermal response parameter (Tmin, Tmax, and q [Equation 3 and 4]; or m and z [Equation 5]) obtained in step 1. The resulting gamma distribution parameters can be used directly to specify the priors in the JAGS model. Because the prior datasets were often very large, in many cases the priors were too strong and over- determined the fit to the focal data. In a few other cases, we had philosophical reasons to strongly constrain the fit to the focal data even when they were sparse (e.g. to constrain Tmax to very high temperatures so that other traits with more information determine the upper thermal limit for R0). Thus, we deflated or inflated the variance as needed (i.e., we fixed the gamma distribution mean but altered the variance by adjusting the parameters that describe the distribution accordingly). See Appendix 1 for more details and specific variance modifications for each thermal response.

Constructing R0 models
When data were missing for a vector–virus pair, we used two criteria to decide which thermal response to use as a substitute: 1) the ecological similarly (i.e. geographic range overlap) of species with available thermal responses, and 2) how restrictive the upper and lower bounds of the available thermal responses were. All else being equal, we chose the more conservative (i.e. least restrictive) option so that R0 would be less likely to be determined by trait thermal responses that did not origi- nate from the focal species. See Appendix 1 for more information about specific models.

When there was more than one option for how to parameterize a model (e.g. vector competence data for WEEV in Cx. tarsalis were available in two forms: separately as b and c, and combined as bc), we calculated R0 both ways. The results were very similar, except for the model for RVFV with lifespan data from Cx. pipiens lifespan in place of Ae. taeniorhynchus (Appendix 1-figure 22). See Appendix 1 for sensitivity and uncertainty methods and Appendix 1—figures 11–20 for results.

Model validation: spatial analysis
We obtained county-level neuroinvasive WNV disease data from 2001 to 2016 for the contiguous US (n = 3109) through the CDC’s county-level disease monitoring program (Centers for Disease Con- trol and Prevention, 2018c). Data were available as total human cases per year, which we adjusted to average cases per 1000 people (using 2010 US county-level census data) to account for popula- tion differences. We averaged cases across years beginning with the first year that had reported cases in a given county to account for the initial spread of WNV and the strong impact of immunity on interannual variation (Paull et al., 2017). Ninety-eight percent of human cases of WNV in the US occur between June and October (data described below), and cases of mosquito-borne disease often lag behind temperature by 1–2 months (Shocket et al., 2018; Stewart Ibarra et al., 2013). Thus, we extracted monthly mean temperature data between the months of May–September for all years between 2001 and 2016 and averaged the data to estimate typical summer conditions for each county. Specifically, we took the centroid geographic coordinate for every county in the contig- uous US with the ‘rgeos’ package Bivand and Rundel, 2012 and extracted corresponding historic climate data for monthly mean temperatures (Climate Research Unit 3.1 rasters) (Harris et al., 2014) from 0.5 ̊2 cells (approximately 2500–3000 km2) using the ‘raster’ package (Hijmans, 2020). The monthly mean temperatures in this climate product are calculated by averaging daily mean tempera- tures at the station level (based on 4–8 observations per day at regular intervals) and interpolating these over a grid (World Meteorological Organization, 2009).

We fit a generalized additive model (GAM) for average incidence as a function of average sum- mer temperature using the ‘mgcv’ package (Wood, 2006). We used a gamma distribution with a log-link function to restrict incidence to positive values and capture heteroskedasticity in the data (i.e. higher variance with higher predicted means), adding a small, near-zero constant (0.0001) to all incidence values to allow the log-transformation for counties with zero incidence. GAMs use additive functions of smooth predictor effects to fit responses that are extremely flexible in the shape of the response. We restricted the number of knots to minimize overfitting (k = 7; see Appendix 1-figure 24 for results across varying values of k). For comparison, we also used the ‘loess’ function in base R ‘stats’ package (R Development Core Team, 2016) to fit locally estimated scatterplot smoothing (LOESS) regressions of the same data. LOESS regression is a simpler but similarly flexible method for estimating the central tendency of data. See Appendix 1-figure 25 for LOESS model results. See Appendix 1-figure 26 for non-binned county-level data.

Model validation: seasonality analysis
We calculated monthly temperature-dependent relative R0 to compare with month-of-onset data for neuroinvasive WNV, EEEV, and SLEV disease aggregated nationwide (the only spatial scale available) from 2001 to 2016 (Centers for Disease Control and Prevention, 2018c; Curren et al., 2018; Lindsey et al., 2018), using the same county-level monthly mean temperature data as above. For WNV, we used the subset of counties with reported cases (68% of counties). For SLEV and EEEV we used all counties from states with reported cases (16 and 20 states, respectively). We calculated a monthly R0(T) for each county, and then weighted each county R0(T) by its population size to calculate a national monthly estimate of R0(T). For WNV, the county-level estimates of R0(T) used models for three Culex species (Cx. pipiens, Cx. quinquefasciatus, and Cx. tarsalis) weighted according to the proportion of WNV-positive mosquitoes reported at the state level, reported in Paull et al., 2017. SLEV and EEEV both only had one R0 model. The estimated monthly temperature-dependent relative R0 values and month-of-onset data were compared visually.

Appendix 1: Model specifications
The equation for R0 (Equation 2 in main text) as a function of temperature (T) that was used in previous analyses (Johnson et al., 2015; Mordecai et al., 2017; Mordecai et al., 2013; Parham and Michael, 2010; Shocket et al., 2018; Tesla et al., 2018) has fecundity measured as eggs per female per day (EFD). Fecundity data were not available directly as eggs per female per day, so we had to transform the available data to obtain the quantities needed for these models. The data for Cx. pipiens were reported as eggs per female per gonotrophic cycle (EFGC). To obtain EFD, we needed to divide EFGC by the length of the gonotrophic cycle. In general, the gonotrophic cycle is assumed to be approximately the inverse of the biting rate. In fact, our ‘biting rate’ (a) data were observations of gonotrophic cycle duration. Accordingly, EFD = EFGC * a, resulting in the equation for R0 below. All but two of the vector–virus parameterizations used this form (Equation A1) of the R0 model (see Appendix 1—table 1, exceptions described below):

R0(T) = ((a^3*bc*exp(mu/PDR)*EFGC*EV*pLA*MDR)/mu^3)^1/2 (eq. A1)

The fecundity data for Cx. quinquefasciatus were reported as eggs per raft (ER). Females lay rafts once per gonotrophic cycle. Thus, in order to obtain an approximation to EFD (eggs per female per day), we again divide by the number of days per gonotrophic cycle and, further, we multiply by the proportion of females ovipositing (pO), since not every female lays an egg raft. These changes result in the following equation for R0:

R0(T) = ((a^3*bc*exp(mu/PDR)*ER*pO*EV*pLA*MDR)/mu^3)^1/2 (eq. A2)

The Cx. quinquefasciatus–WNV model used Equation A2. The Ae. triseriatus–EEEV model also used Equation A2 (i.e., included pO) but substituted the Cx. pipiens thermal response for EFGC in place of the Cx. quinquefasciatus thermal response for ER for the following reasons. There were no fecundity trait data available for Ae. triseriatus. (Ae. triseratus was chosen as the focal species for the EEEV model because it is the only species with temperature-dependent vector competence data available, and it is a possible bridge vector for EEEV transmis- sion to humans). Cs. melanura is the primary vector for maintaining enzootic cycles of EEEV in birds (Mahmood and Crans, 1998), more often cited in the literature in association with EEEV (e.g. [Weaver and Barrett, 2004]), and had data for pO (proportion ovipositing) available. Thus, we chose to include this thermal response in model because it contained information that could affect the upper and lower bounds of transmission (even though most models did not include pO [proportion ovipositing], because they use the Cx. pipiens EFGC [eggs per female per gonotrophic cycle] thermal response that includes pO implicitly). Then we needed to choose which egg production metric to include. We chose the Cx. pipiens EFGC thermal response over the Cx. quinquefasciatus ER thermal response because the former was the better choice according to both criteria: Cx. pipiens has a more similar species range to Ae. triseriatus and Cs. melanura and its thermal response was slightly more conservative (less restrictive = cooler lower thermal limit and warmer upper thermal limit). Although technically the units are not correct (see above), the thermal responses for Cx. pipiens EFGC and Cx. quinquefasciatus ER are so similar despite having different units (Figure 4B), we decided that the other two criteria were more important than being strict with regard to the units, as it is feasible to have an ER thermal response that is quite similar to the EFGC thermal response. Ultimately, because the thermal responses for EFGC and ER are so similar, this decision only has a small impact on the R0 results (see Appendix 1-figure 22 comparing four alternative model specifications/parameterizations for the Ae. triseriatus-EEEV model).
