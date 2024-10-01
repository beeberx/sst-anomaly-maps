# sst-anomaly-maps

## Name
Sea Surface Temperature Anomaly Maps of different SST data products.  
 

## Description
Create Sea Surface Temperature Anomaly Maps of different SST data products.  Anomalies (difference from reference period mean) and Normalised Anomalies (difference from reference period mean, scaled by standard deviation) are calculated and plots can be generated for specific year, for range of years, or for entire period.  Currently able to create monthly averages of OISST V2 High Res (from NOAA, USA) and of HADISST (from Met Office Hadley Centre, UK)


## How to use these scripts

- All necessary functions and files, beyond those included in Matlab as standard, are part of the repository. 
- The repository includes reproduced copies of cbrewer functions and of the coastline shared under m_map. 
- Open script in Matlab, make necessary edits to point to locations of your data files, time period, latitudinal and longitudinal region of interest. 
- Run scripts.  Figures are saved in sub-folders of main repository location.  


## Data Dependencies

- Hadley Centre HADISST data can be downloaded from https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
- NOAA OISST V2 High Res data can be downloaded from ftp2.psl.noaa.gov/Projects/Datasets/noaa.oisst.v2.highres


## Authors and acknowledgment
Project created by Dr Bee Berx. 

## License
Published under Open Government License. 

## Project status
Under development. 
