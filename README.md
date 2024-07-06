#  Exploring Associations of Pediatric Medical Complexity and Neighborhood Opportunity with Developmental Outcomes in High-Risk Infants Living in Southern California

This code assessed key determinants of infant development among those from a HRIF clinic located in an urban children’s hospital in Southern California. We investigated associations between infant medical complexity, as defined by the Pediatric Medical Complexity Algorithm (PMCA), and infant neighbourhoods using the Child Opportunity Index (COI), with infant developmental domains, as measured by the Bayley Scales of Infant and Toddler Development (BSID), which is considered the gold standard for evaluating infant cognitive, language and motor development from birth to 42 months. 

The input data utilized was deidentified medical record data. This may be available upon request with corresponding authors. The scripts included in this repo, along with their functionality, are outlined below. 

- `1.PMCA_classification.R` : This script runs the pediatric medical complexity algorithm on first encounter diagnostic codes of infants. 

- `2.cleaning_data.Rmd` : A script to read in raw data and clean and process the dataframes for downstream use. This includes things such as changing variable names, dealing with NA or out of range values, and assessing baseline summary stats.

- `3.stats_figures.Rmd` : This script runs the main models, statistics, and figures included in the paper.

Should you have any further questions, please reach out to emily.yeo@colorado.edu. Details on where to find the paper will be included upon publication. 

![05D2CCEF-80A8-46AD-B70D-AD2C5F6A1C0F_1_201_a](https://github.com/emilynyeo/infant_developement_PMCA_COI/assets/104112036/95ab97f6-3dca-464f-828f-c9bb0d3cab11)

# Reproducible Code Artifacts 

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/emilynyeo/infant_developement_PMCA_COI">Code Artifact for "Exploring Associations of Pediatric Medical Complexity and Neighborhood Opportunity with Developmental Outcomes in High-Risk Infants Living in Southern California"</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://www.linkedin.com/in/emily-nadine-yeo/">Emily Nadine Yeo, University of Colorado Boulder.</a> is licensed under <a href="https://creativecommons.org/licenses/by-sa/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY-SA 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/sa.svg?ref=chooser-v1" alt=""></a></p>
