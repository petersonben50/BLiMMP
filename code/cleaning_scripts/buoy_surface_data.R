#### code/cleaning_scripts/buoy_surface_data.R ####
# Benjamin D. Peterson

# This script prepares a set of data to generate
# a time course of the chlorophyll

# Links:
# NTL-LTER dataset entry: https://lter.limnology.wisc.edu/lake_mendota_multiparameter_sonde_profiles
# Project on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.400.4
# Code generated from EDI site

#### This is your mess, you need to clean it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(tidyverse)



#### Read in and combine data ####

inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/129/35/72494d432fe1e977f5326100a733cece" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "year4",     
                 "sampledate",     
                 "hour",     
                 "avg_air_temp",     
                 "flag_avg_air_temp",     
                 "avg_rel_hum",     
                 "flag_avg_rel_hum",     
                 "avg_wind_speed",     
                 "flag_avg_wind_speed",     
                 "avg_wind_dir",     
                 "flag_avg_wind_dir",     
                 "avg_chlor_rfu",     
                 "flag_avg_chlor_rfu",     
                 "avg_phyco_rfu",     
                 "flag_avg_phyco_rfu",     
                 "avg_par",     
                 "flag_avg_par",     
                 "avg_par_below",     
                 "flag_avg_par_below",     
                 "avg_do_wtemp",     
                 "flag_avg_do_wtemp",     
                 "avg_do_sat",     
                 "flag_avg_do_sat",     
                 "avg_do_raw",     
                 "flag_avg_do_raw",     
                 "avg_pco2_ppm",     
                 "flag_avg_pco2_ppm",     
                 "avg_ph",     
                 "flag_avg_ph",     
                 "avg_spec_cond",     
                 "flag_avg_spec_cond",     
                 "avg_fdom",     
                 "flag_avg_fdom",     
                 "avg_turbidity",     
                 "flag_avg_turbidity"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$year4)=="factor") dt2$year4 <-as.numeric(levels(dt2$year4))[as.integer(dt2$year4) ]               
if (class(dt2$year4)=="character") dt2$year4 <-as.numeric(dt2$year4)                                   
# attempting to convert dt2$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2sampledate<-as.Date(dt2$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2sampledate) == length(tmp2sampledate[!is.na(tmp2sampledate)])){dt2$sampledate <- tmp2sampledate } else {print("Date conversion failed for dt2$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2sampledate) 
if (class(dt2$hour)=="factor") dt2$hour <-as.numeric(levels(dt2$hour))[as.integer(dt2$hour) ]               
if (class(dt2$hour)=="character") dt2$hour <-as.numeric(dt2$hour)
if (class(dt2$avg_air_temp)=="factor") dt2$avg_air_temp <-as.numeric(levels(dt2$avg_air_temp))[as.integer(dt2$avg_air_temp) ]               
if (class(dt2$avg_air_temp)=="character") dt2$avg_air_temp <-as.numeric(dt2$avg_air_temp)
if (class(dt2$flag_avg_air_temp)!="factor") dt2$flag_avg_air_temp<- as.factor(dt2$flag_avg_air_temp)
if (class(dt2$avg_rel_hum)=="factor") dt2$avg_rel_hum <-as.numeric(levels(dt2$avg_rel_hum))[as.integer(dt2$avg_rel_hum) ]               
if (class(dt2$avg_rel_hum)=="character") dt2$avg_rel_hum <-as.numeric(dt2$avg_rel_hum)
if (class(dt2$flag_avg_rel_hum)!="factor") dt2$flag_avg_rel_hum<- as.factor(dt2$flag_avg_rel_hum)
if (class(dt2$avg_wind_speed)=="factor") dt2$avg_wind_speed <-as.numeric(levels(dt2$avg_wind_speed))[as.integer(dt2$avg_wind_speed) ]               
if (class(dt2$avg_wind_speed)=="character") dt2$avg_wind_speed <-as.numeric(dt2$avg_wind_speed)
if (class(dt2$flag_avg_wind_speed)!="factor") dt2$flag_avg_wind_speed<- as.factor(dt2$flag_avg_wind_speed)
if (class(dt2$avg_wind_dir)=="factor") dt2$avg_wind_dir <-as.numeric(levels(dt2$avg_wind_dir))[as.integer(dt2$avg_wind_dir) ]               
if (class(dt2$avg_wind_dir)=="character") dt2$avg_wind_dir <-as.numeric(dt2$avg_wind_dir)
if (class(dt2$flag_avg_wind_dir)!="factor") dt2$flag_avg_wind_dir<- as.factor(dt2$flag_avg_wind_dir)
if (class(dt2$avg_chlor_rfu)=="factor") dt2$avg_chlor_rfu <-as.numeric(levels(dt2$avg_chlor_rfu))[as.integer(dt2$avg_chlor_rfu) ]               
if (class(dt2$avg_chlor_rfu)=="character") dt2$avg_chlor_rfu <-as.numeric(dt2$avg_chlor_rfu)
if (class(dt2$flag_avg_chlor_rfu)!="factor") dt2$flag_avg_chlor_rfu<- as.factor(dt2$flag_avg_chlor_rfu)
if (class(dt2$avg_phyco_rfu)=="factor") dt2$avg_phyco_rfu <-as.numeric(levels(dt2$avg_phyco_rfu))[as.integer(dt2$avg_phyco_rfu) ]               
if (class(dt2$avg_phyco_rfu)=="character") dt2$avg_phyco_rfu <-as.numeric(dt2$avg_phyco_rfu)
if (class(dt2$flag_avg_phyco_rfu)!="factor") dt2$flag_avg_phyco_rfu<- as.factor(dt2$flag_avg_phyco_rfu)
if (class(dt2$avg_par)=="factor") dt2$avg_par <-as.numeric(levels(dt2$avg_par))[as.integer(dt2$avg_par) ]               
if (class(dt2$avg_par)=="character") dt2$avg_par <-as.numeric(dt2$avg_par)
if (class(dt2$flag_avg_par)!="factor") dt2$flag_avg_par<- as.factor(dt2$flag_avg_par)
if (class(dt2$avg_par_below)=="factor") dt2$avg_par_below <-as.numeric(levels(dt2$avg_par_below))[as.integer(dt2$avg_par_below) ]               
if (class(dt2$avg_par_below)=="character") dt2$avg_par_below <-as.numeric(dt2$avg_par_below)
if (class(dt2$flag_avg_par_below)!="factor") dt2$flag_avg_par_below<- as.factor(dt2$flag_avg_par_below)
if (class(dt2$avg_do_wtemp)=="factor") dt2$avg_do_wtemp <-as.numeric(levels(dt2$avg_do_wtemp))[as.integer(dt2$avg_do_wtemp) ]               
if (class(dt2$avg_do_wtemp)=="character") dt2$avg_do_wtemp <-as.numeric(dt2$avg_do_wtemp)
if (class(dt2$flag_avg_do_wtemp)!="factor") dt2$flag_avg_do_wtemp<- as.factor(dt2$flag_avg_do_wtemp)
if (class(dt2$avg_do_sat)=="factor") dt2$avg_do_sat <-as.numeric(levels(dt2$avg_do_sat))[as.integer(dt2$avg_do_sat) ]               
if (class(dt2$avg_do_sat)=="character") dt2$avg_do_sat <-as.numeric(dt2$avg_do_sat)
if (class(dt2$flag_avg_do_sat)!="factor") dt2$flag_avg_do_sat<- as.factor(dt2$flag_avg_do_sat)
if (class(dt2$avg_do_raw)=="factor") dt2$avg_do_raw <-as.numeric(levels(dt2$avg_do_raw))[as.integer(dt2$avg_do_raw) ]               
if (class(dt2$avg_do_raw)=="character") dt2$avg_do_raw <-as.numeric(dt2$avg_do_raw)
if (class(dt2$flag_avg_do_raw)!="factor") dt2$flag_avg_do_raw<- as.factor(dt2$flag_avg_do_raw)
if (class(dt2$avg_pco2_ppm)=="factor") dt2$avg_pco2_ppm <-as.numeric(levels(dt2$avg_pco2_ppm))[as.integer(dt2$avg_pco2_ppm) ]               
if (class(dt2$avg_pco2_ppm)=="character") dt2$avg_pco2_ppm <-as.numeric(dt2$avg_pco2_ppm)
if (class(dt2$flag_avg_pco2_ppm)!="factor") dt2$flag_avg_pco2_ppm<- as.factor(dt2$flag_avg_pco2_ppm)
if (class(dt2$avg_ph)=="factor") dt2$avg_ph <-as.numeric(levels(dt2$avg_ph))[as.integer(dt2$avg_ph) ]               
if (class(dt2$avg_ph)=="character") dt2$avg_ph <-as.numeric(dt2$avg_ph)
if (class(dt2$flag_avg_ph)!="factor") dt2$flag_avg_ph<- as.factor(dt2$flag_avg_ph)
if (class(dt2$avg_spec_cond)=="factor") dt2$avg_spec_cond <-as.numeric(levels(dt2$avg_spec_cond))[as.integer(dt2$avg_spec_cond) ]               
if (class(dt2$avg_spec_cond)=="character") dt2$avg_spec_cond <-as.numeric(dt2$avg_spec_cond)
if (class(dt2$flag_avg_spec_cond)!="factor") dt2$flag_avg_spec_cond<- as.factor(dt2$flag_avg_spec_cond)
if (class(dt2$avg_fdom)=="factor") dt2$avg_fdom <-as.numeric(levels(dt2$avg_fdom))[as.integer(dt2$avg_fdom) ]               
if (class(dt2$avg_fdom)=="character") dt2$avg_fdom <-as.numeric(dt2$avg_fdom)
if (class(dt2$flag_avg_fdom)!="factor") dt2$flag_avg_fdom<- as.factor(dt2$flag_avg_fdom)
if (class(dt2$avg_turbidity)=="factor") dt2$avg_turbidity <-as.numeric(levels(dt2$avg_turbidity))[as.integer(dt2$avg_turbidity) ]               
if (class(dt2$avg_turbidity)=="character") dt2$avg_turbidity <-as.numeric(dt2$avg_turbidity)
if (class(dt2$flag_avg_turbidity)!="factor") dt2$flag_avg_turbidity<- as.factor(dt2$flag_avg_turbidity)


#### Select needed data and read it out ####
pigment_data <- dt2 %>%
  filter(year(sampledate) %in% c(2020, 2021)) %>%
  select(sampledate, avg_chlor_rfu, avg_phyco_rfu, avg_turbidity) %>%
  rename(sampleDate = sampledate) %>%
  filter(!is.na(avg_chlor_rfu)) %>%
  group_by(sampleDate) %>%
  summarise(mean_chlor_rfu = mean(avg_chlor_rfu),
            sd_chlor_rfu = sd(avg_chlor_rfu),
            mean_phyco_rfu = mean(avg_phyco_rfu),
            sd_phyco_rfu = sd(avg_phyco_rfu),
            mean_turb = mean(avg_turbidity),
            sd_turb = sd(avg_turbidity))
write.csv(pigment_data,
          "dataFinal/pigment_profiles_data.csv",
          row.names = FALSE)
