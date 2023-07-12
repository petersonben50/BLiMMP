#### code/cleaning_scripts/DO_data_for_heatmap.R ####
# Benjamin D. Peterson

# This script prepares a set of DO data to generate
# a heatmap of DO over 2020 and 2021 summers.
# This data was collected by Mark Gahler for the 

# Links:
# NTL-LTER dataset entry: https://lter.limnology.wisc.edu/lake_mendota_multiparameter_sonde_profiles
# Project on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.400.4
# Code generated from EDI site

#### This is your mess, you need to clean it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(tidyverse)



#### Read in and combine data ####
# Package ID: knb-lter-ntl.400.4 Cataloging System:https://pasta.edirepository.org.
# Data set title: Lake Mendota Multiparameter Sonde Profiles: 2017 - current.
# Data set creator:  John Magnuson - University of Wisconsin-Madison 
# Data set creator:  Stephen Carpenter - University of Wisconsin-Madison 
# Data set creator:  Emily Stanley - University of Wisconsin-Madison 
# Metadata Provider:  NTL Information Manager - University of Wisconsin-Madison 
# Contact:    -  NTL LTER  - ntl.infomgr@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/400/4/8fd57df54b19ddd592dbb45321bbd548" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lakeid",     
                 "sampledate",     
                 "sampletime",     
                 "depth",     
                 "wtemp",     
                 "do_raw",     
                 "do_sat",     
                 "ph",     
                 "spec_cond",     
                 "chlor_rfu",     
                 "phyco_rfu",     
                 "turbidity",     
                 "fdom",     
                 "flag_wtemp",     
                 "flag_do_raw",     
                 "flag_do_sat",     
                 "flag_ph",     
                 "flag_spec_cond",     
                 "flag_chlor_rfu",     
                 "flag_phyco_rfu",     
                 "flag_turbidity",     
                 "flag_fdom"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)                                   
# attempting to convert dt1$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1sampledate<-as.Date(dt1$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1sampledate) == length(tmp1sampledate[!is.na(tmp1sampledate)])){dt1$sampledate <- tmp1sampledate } else {print("Date conversion failed for dt1$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1sampledate) 
if (class(dt1$depth)=="factor") dt1$depth <-as.numeric(levels(dt1$depth))[as.integer(dt1$depth) ]               
if (class(dt1$depth)=="character") dt1$depth <-as.numeric(dt1$depth)
if (class(dt1$wtemp)=="factor") dt1$wtemp <-as.numeric(levels(dt1$wtemp))[as.integer(dt1$wtemp) ]               
if (class(dt1$wtemp)=="character") dt1$wtemp <-as.numeric(dt1$wtemp)
if (class(dt1$do_raw)=="factor") dt1$do_raw <-as.numeric(levels(dt1$do_raw))[as.integer(dt1$do_raw) ]               
if (class(dt1$do_raw)=="character") dt1$do_raw <-as.numeric(dt1$do_raw)
if (class(dt1$do_sat)=="factor") dt1$do_sat <-as.numeric(levels(dt1$do_sat))[as.integer(dt1$do_sat) ]               
if (class(dt1$do_sat)=="character") dt1$do_sat <-as.numeric(dt1$do_sat)
if (class(dt1$ph)=="factor") dt1$ph <-as.numeric(levels(dt1$ph))[as.integer(dt1$ph) ]               
if (class(dt1$ph)=="character") dt1$ph <-as.numeric(dt1$ph)
if (class(dt1$spec_cond)=="factor") dt1$spec_cond <-as.numeric(levels(dt1$spec_cond))[as.integer(dt1$spec_cond) ]               
if (class(dt1$spec_cond)=="character") dt1$spec_cond <-as.numeric(dt1$spec_cond)
if (class(dt1$chlor_rfu)=="factor") dt1$chlor_rfu <-as.numeric(levels(dt1$chlor_rfu))[as.integer(dt1$chlor_rfu) ]               
if (class(dt1$chlor_rfu)=="character") dt1$chlor_rfu <-as.numeric(dt1$chlor_rfu)
if (class(dt1$phyco_rfu)=="factor") dt1$phyco_rfu <-as.numeric(levels(dt1$phyco_rfu))[as.integer(dt1$phyco_rfu) ]               
if (class(dt1$phyco_rfu)=="character") dt1$phyco_rfu <-as.numeric(dt1$phyco_rfu)
if (class(dt1$turbidity)=="factor") dt1$turbidity <-as.numeric(levels(dt1$turbidity))[as.integer(dt1$turbidity) ]               
if (class(dt1$turbidity)=="character") dt1$turbidity <-as.numeric(dt1$turbidity)
if (class(dt1$fdom)=="factor") dt1$fdom <-as.numeric(levels(dt1$fdom))[as.integer(dt1$fdom) ]               
if (class(dt1$fdom)=="character") dt1$fdom <-as.numeric(dt1$fdom)
if (class(dt1$flag_wtemp)!="factor") dt1$flag_wtemp<- as.factor(dt1$flag_wtemp)
if (class(dt1$flag_do_raw)!="factor") dt1$flag_do_raw<- as.factor(dt1$flag_do_raw)
if (class(dt1$flag_do_sat)!="factor") dt1$flag_do_sat<- as.factor(dt1$flag_do_sat)
if (class(dt1$flag_ph)!="factor") dt1$flag_ph<- as.factor(dt1$flag_ph)
if (class(dt1$flag_spec_cond)!="factor") dt1$flag_spec_cond<- as.factor(dt1$flag_spec_cond)
if (class(dt1$flag_chlor_rfu)!="factor") dt1$flag_chlor_rfu<- as.factor(dt1$flag_chlor_rfu)
if (class(dt1$flag_phyco_rfu)!="factor") dt1$flag_phyco_rfu<- as.factor(dt1$flag_phyco_rfu)
if (class(dt1$flag_turbidity)!="factor") dt1$flag_turbidity<- as.factor(dt1$flag_turbidity)
if (class(dt1$flag_fdom)!="factor") dt1$flag_fdom<- as.factor(dt1$flag_fdom)



#### Select needed data and read it out ####
DO_data <- dt1 %>%
  select(sampledate, depth, do_raw) %>%
  filter(year(sampledate) %in% c(2020, 2021)) %>%
  rename(sampleDate = sampledate)
write.csv(DO_data,
          "dataFinal/DO_profiles_data.csv",
          row.names = FALSE)
