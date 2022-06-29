
library(logspline)
setOldClass("oldlogspline")
library(KernSmooth)
library(scales) #for alpha function
library(zoo) #to fill in NAs from the dpt1 estimates
library(survival) #for the PDF of age of vaccination
library(countrycode) #to transfer b/t uncode and iso3codes and country names
library(readxl)
library(dplyr)
library(wpp2017)

dyn.load("./source/MRModel-funcs.so") #run "R CMD SHLIB source/MRModel-funcs.c" in terminal to compile
source("./source/build.R")
source("./source/base.R")
source("./source/user_interface.R")
source("./source/vimc_inputs.R")

setup <- setupCountry_202110gavi_v3(country="ZMB")
year <- 1980 
t.max <- length(year:2100) 

#Function to set up experiment and run the transients out
EXt0 <- EX <- EX.Country.part1(uncode=setup$uncode,
                               generation.time = 0.5, 
                               age.classes = c(1:240, seq(252,1212,12)),
                               maternal.decay.rt = 0.45, #based on Leuridan + others wanes b/w 3-9 months
                               exponent = 0.97,
                               frequency.dep=T,
                               yr.births.per.1000.acrossyears=rep(setup$cbr.1950.2100[year-1950+1],20),
                               targeted.intro=F,
                               intro.rate=1/24/320,
                               tot.pop=NULL, #it will pull a pop automatically from pop distribution if NULL
                               R0 = 5,
                               t.max = 20,
                               get.births=setup$get.births.here,
                               seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                               flat.WAIFW=F,
                               country.specific.WAIFW=T,
                               vynnycky.waifw=F,
                               vynnycky.waifw.betas = c(2,1),
                               asdr.object=setup$asdr.object,
                               year=year,
                               use_montagu_demog=T,
                               routine.vac=0,
                               routine.vac.age.index=12)
print("transients done")


#NO VACCINATION SCENARIO
tmpn <- EX.Country.part2(uncode = setup$uncode,
                         generation.time = 0.5, 
                         pop.rescale=setup$pop.total.1950.2100[seq(1990,(year+t.max-10),10)-1950+1], 
                         pop.time=seq(1990,(year+t.max-10),10)-year, 
                         is.stochastic=TRUE,
                         get.births=setup$get.births.here,
                         t.max = length(year:2100),
                         rescale.WAIFW=F,
                         yr.births.per.1000.acrossyears=setup$cbr.1950.2100[(year-1950+1):((year-1950+1)+t.max)], 
                         asdr.object=setup$asdr.object, 
                         year=1980,
                         EXt0=EXt0,
                         time.specific.MR1cov =  rep(0, length(1:t.max)), 
                         time.specific.MR2cov = rep(0, length(1:t.max)),
                         time.specific.SIAcov =   rep(0, length(1:t.max)), 
                         time.specific.min.age.MR1 = rep(0, length(1:t.max)),
                         time.specific.max.age.MR1 = rep(0, length(1:t.max)),
                         time.specific.min.age.MR2 = rep(0, length(1:t.max)),
                         time.specific.max.age.MR2 = rep(0, length(1:t.max)),
                         time.specific.min.age.SIA = rep(0, length(1:t.max)),
                         time.specific.max.age.SIA = rep(0, length(1:t.max)),
                         obj.vcdf.MR1 = get.MR1cdf.survival(uncode=24, 1, 24), #get.vcdf.uniform(12, 23)
                         obj.vcdf.MR2 = get.vcdf.normal(25, 36), #get.vcdf.uniform(24, 35)
                         obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                         sia.timing.in.year = (3/12),
                         MR1MR2correlation = FALSE,
                         MR1SIAcorrelation = FALSE,
                         MR2SIAcorrelation = FALSE,
                         SIAinacc = FALSE,
                         prop.inacc = rep(0, length(1:t.max)),
                         SIAinefficient = FALSE)
print("no vaccine done")
plot(tmpn@result)
tmpr <- tmpn
#Get output and storing as _n for no vaccine scenario
output_n <- getOutput_201910gavi_v4(sim.res=tmpr, t.max=t.max, setup=setup, year=year, 
                                    crs_rate_gestational_age_nbiweeks=c(6,8), 
                                    crs_rate_agespecific=c(0.4748415, 0.1664943), 
                                    crs_rate_overall=(0.4748415*0.75)+(0.1664943*0.25), 
                                    fetal_death_rate=0.08802422, 
                                    infant_death_rate=0.24344971,
                                    scenario="no vaccine", run.number=s, R0=R0)



#RCV2 BESTCASE SCENARIO
cov.scenario <- setup$coverage$rcv2_bestcase
tmp_rcv2_best <- EX.Country.part2(uncode = setup$uncode,
                                  generation.time = 0.5, 
                                  pop.rescale=setup$pop.total.1950.2100[seq(1990,(year+t.max-10),10)-1950+1], 
                                  pop.time=seq(1990,(year+t.max-10),10)-year, 
                                  is.stochastic=TRUE,
                                  get.births=setup$get.births.here,
                                  t.max = length(year:2100),
                                  rescale.WAIFW=F,
                                  yr.births.per.1000.acrossyears=setup$cbr.1950.2100[(year-1950+1):((year-1950+1)+t.max)],  
                                  asdr.object=setup$asdr.object, 
                                  year=1980,
                                  EXt0=EXt0,
                                  time.specific.MR1cov =  cov.scenario@MR1.cov, 
                                  time.specific.MR2cov = cov.scenario@MR2.cov,  
                                  time.specific.SIAcov =   cov.scenario@sia.cov, 
                                  time.specific.min.age.MR1 = rep(1, length(cov.scenario@MR1.cov)),
                                  time.specific.max.age.MR1 = rep(24, length(cov.scenario@MR1.cov)),
                                  time.specific.min.age.MR2 = rep(25, length(cov.scenario@MR1.cov)),
                                  time.specific.max.age.MR2 = rep(36, length(cov.scenario@MR1.cov)),
                                  time.specific.min.age.SIA = cov.scenario@sia.min.age,
                                  time.specific.max.age.SIA = cov.scenario@sia.max.age,
                                  obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24), #get.vcdf.uniform(12, 23)
                                  obj.vcdf.MR2 = get.vcdf.normal(25, 36), #get.vcdf.uniform(24, 35)
                                  obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                                  sia.timing.in.year = (1/12),
                                  MR1MR2correlation = TRUE,
                                  MR1SIAcorrelation = FALSE,
                                  MR2SIAcorrelation = FALSE,
                                  SIAinacc = TRUE,
                                  prop.inacc = cov.scenario@inaccessible.prop,
                                  SIAinefficient = TRUE) 
print("RCV2 BESTCASE done")
plot(tmp_rcv2_best@result)
tmpr <- tmp_rcv2_best
#Get output and storing as _rcv2_best for RCV2 Bestcase scenario
output_rcv2_best <- getOutput_201910gavi_v4(sim.res=tmpr, t.max=t.max, setup=setup, year=year, 
                                            crs_rate_gestational_age_nbiweeks=c(6,8), 
                                            crs_rate_agespecific=c(0.4748415, 0.1664943), 
                                            crs_rate_overall=(0.4748415*0.75)+(0.1664943*0.25), 
                                            fetal_death_rate=0.08802422, 
                                            infant_death_rate=0.24344971,
                                            scenario="rcv2 bestcase", run.number=s, R0=R0)


