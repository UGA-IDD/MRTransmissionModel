## Nov 8, 2023 - I am working from this file to modify the R package - note there are lots of changes that have not been pushed to Git
##             - My goal right now is to make sure the model is working for no school or non-spatial models.

rm(list = ls())
gc()

library(MRTransmissionModel)

setup <- setupCountry.Nov2023(country="Montenegro")
year <- 2022
t.max <- length(year:2100)

#Function to set up experiment and run the transients out
EXt0 <- EX <- EX.Country.part1(uncode=setup$uncode,
                               generation.time = 0.5,
                               age.classes = c(1:240, seq(252,1212,12)),
                               maternal.decay.rt = 0.45, #based on Leuridan + others wanes b/w 3-9 months
                               exponent = 0.97,
                               frequency.dep=TRUE,
                               yr.births.per.1000.acrossyears=rep(setup$cbr.1950.2100[year-1950+1],20),
                               intro.rate=1/24/320,
                               tot.pop=NULL, #it will pull a pop automatically from pop distribution if NULL
                               targeted.intro=FALSE,
                               R0 = 5,
                               t.max = 20,
                               get.births=setup$get.births.here,
                               seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                               flat.WAIFW=FALSE,
                               country.specific.WAIFW=TRUE,
                               vynnycky.waifw=FALSE,
                               vynnycky.waifw.betas = c(2,1),
                               asdr.object=setup$asdr.object,
                               year=year,
                               use_montagu_demog=FALSE,
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




## VACCINATION SCENARIO
tmp1 <- EX.Country.part2(uncode = setup$uncode,
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
                         time.specific.MR1cov =  setup$MCV1.coverage.1980to2100[(year-1980+1):t.max],
                         time.specific.MR2cov = setup$MCV2.coverage.1980to2100[(year-1980+1):t.max],
                         time.specific.SIAcov =  setup$measlesSIA.coverage.1980to2100[(year-1980+1):t.max],
                         time.specific.min.age.MR1 = rep(9, length((year-1980+1):t.max)),
                         time.specific.max.age.MR1 = rep(24, length((year-1980+1):t.max)),
                         time.specific.min.age.MR2 = rep(25, length((year-1980+1):t.max)),
                         time.specific.max.age.MR2 = rep(36, length((year-1980+1):t.max)),
                         time.specific.min.age.SIA = setup$age.min.sia.measles[(year-1980+1):t.max],
                         time.specific.max.age.SIA = setup$age.max.sia.measles[(year-1980+1):t.max],
                         obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24), #get.vcdf.uniform(12, 23)
                         obj.vcdf.MR2 = get.vcdf.normal(25, 36), #get.vcdf.uniform(24, 35)
                         obj.prob.vsucc = pvacsuccess(c(1:240, seq(252,1212,12)), get.boulianne.vsucc()),
                         sia.timing.in.year = (3/12),
                         MR1MR2correlation = FALSE,
                         MR1SIAcorrelation = FALSE,
                         MR2SIAcorrelation = FALSE,
                         SIAinacc = FALSE,
                         prop.inacc = rep(0, length(1:t.max)),
                         SIAinefficient = FALSE)
plot(tmp1@result)


## Susceptible population
#Function to set up experiment and run the transients out
EXt0 <- EX.Country.part1.setSUS(uncode=setup$uncode,
                               generation.time = 0.5,
                               age.classes = c(1:240, seq(252,1212,12)),
                               maternal.decay.rt = 0.45, #based on Leuridan + others wanes b/w 3-9 months
                               exponent = 0.97,
                               frequency.dep=TRUE,
                               yr.births.per.1000.acrossyears=rep(setup$cbr.1950.2100[year-1950+1],20),
                               intro.rate=1/24/320,
                               tot.pop=NULL, #it will pull a pop automatically from pop distribution if NULL
                               targeted.intro=FALSE,
                               R0 = 5,
                               t.max = 20,
                               get.births=setup$get.births.here,
                               seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                               flat.WAIFW=FALSE,
                               country.specific.WAIFW=TRUE,
                               vynnycky.waifw=FALSE,
                               vynnycky.waifw.betas = c(2,1),
                               asdr.object=setup$asdr.object,
                               year=year,
                               use_montagu_demog=FALSE,
                               routine.vac=0,
                               routine.vac.age.index=12,
                               starting.prop.immune = c(0.95, 0.95, 0.95, 0.95),
                               starting.prop.immune.ages.in.months = c(9, 12, 24, 12*40))
print("transients done")



