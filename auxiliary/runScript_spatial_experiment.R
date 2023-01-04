rm(list = ls())
gc()

library(MRTransmissionModel)

setup <- setupCountry.Dec2021(country="Zambia")
setup$asdr.object
year <- 1980
generation.time <- 0.5
age.classes <- c(1:240, seq(252,1212,12))
t.max <- length(year:2100)

## Setting up demography spatial example - this is temporary b/c still need to create two new functions
# 1. space.setupCountry()
# 2. space.getDemography()
setup$pop.total.1950.2100 <- rbind(setup$pop.total.1950.2100, setup$pop.total.1950.2100*0.5)
setup$pop.age.byageclasses.1950.2100 <- rbind(setup$pop.age.byageclasses.1950.2100, setup$pop.age.byageclasses.1950.2100*0.5)
setup$tfr.1950.2100 <- rbind(setup$tfr.1950.2100, setup$tfr.1950.2100)
setup$e0.1950.2100 <- rbind(setup$e0.1950.2100, setup$e0.1950.2100)
setup$asfr.1950.2100 <- rbind(setup$asfr.1950.2100, setup$asfr.1950.2100)
setup$repro.age.sex.dist.1950.2100 <- rbind(setup$repro.age.sex.dist.1950.2100, setup$repro.age.sex.dist.1950.2100)
setup$births.1950.2100 <- rbind(setup$births.1950.2100, setup$births.1950.2100)
setup$cbr.1950.2100 <- rbind(setup$cbr.1950.2100, setup$cbr.1950.2100*0.95)
setup$cbr.1950.2100 <- cbind(setup$cbr.1950.2100, setup$cbr.1950.2100[,ncol(setup$cbr.1950.2100)])
setup$yr.agespecificbirths.per.1000 <- rbind(setup$yr.agespecificbirths.per.1000, setup$yr.agespecificbirths.per.1000)
setup$asdr.1950.2100.by5 <- rbind(setup$asdr.1950.2100.by5, setup$asdr.1950.2100.by5)
setup$asdr.object <- space.nMx <- new("space.nMx",
                                      rate.years = setup$asdr.object@rate.years,
                                      rates = rbind(setup$asdr.object@rates, setup$asdr.object@rates*0.9),
                                      mid.age = setup$asdr.object@mid.age,
                                      n.subpops = 2)

# Run the transients out

# uncode=setup$uncode
# generation.time = 0.5
# age.classes = c(1:240, seq(252,1212,12))
# maternal.decay.rt = 0.45 #based on Leuridan + others wanes b/w 3-9 months
# exponent = 0.97
# frequency.dep=TRUE
# yr.births.per.1000.bysubpop <- setup$cbr.1950.2100[,year-1950+1]
# targeted.intro=FALSE
# intro.rate=1/24/320
# tot.subpop=NULL #it will pull a pop automatically from pop distribution if NULL
# R0 = 5
# t.max = 20
# get.births=setup$get.births.here
# seasonal.amp = 0.15 #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
# flat.WAIFW=FALSE
# country.specific.WAIFW=TRUE
# vynnycky.waifw=FALSE
# vynnycky.waifw.betas = c(2,1)
# space.asdr.object=setup$asdr.object
# year=year
# routine.vac=0
# routine.vac.age.index=12
# n.subpops=2
# coupling = matrix(1, nrow=n.subpops, ncol=n.subpops)

EXt0 <- EX <- Spatial.EX.Country.part1(uncode=setup$uncode,
                                       generation.time = 0.5,
                                       age.classes = c(1:240, seq(252,1212,12)),
                                       maternal.decay.rt = 0.45, #based on Leuridan + others wanes b/w 3-9 months
                                       exponent = 0.97,
                                       frequency.dep=TRUE,
                                       yr.births.per.1000.bysubpop <- setup$cbr.1950.2100[,year-1950+1],
                                       targeted.intro=FALSE,
                                       intro.rate=1/24/320,
                                       tot.subpop=NULL, #it will pull a pop automatically from pop distribution if NULL
                                       R0 = 5,
                                       t.max = 20,
                                       get.births=setup$get.births.here,
                                       seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                                       flat.WAIFW=FALSE,
                                       country.specific.WAIFW=TRUE,
                                       vynnycky.waifw=FALSE,
                                       vynnycky.waifw.betas = c(2,1),
                                       space.asdr.object=setup$asdr.object,
                                       year=year,
                                       routine.vac=0,
                                       routine.vac.age.index=12,
                                       n.subpops=2,
                                       coupling = matrix(1, nrow=n.subpops, ncol=n.subpops))


print("transients done")

# Run Vaccination Experiment
rcv_exp <- Spatial.EX.Country.part2(uncode = setup$uncode,
                                    generation.time = 0.5,
                                    t.max = t.max,
                                    pop.rescale=setup$pop.total.1950.2100[,seq(1990,(year+t.max-10),10)-1950+1] ,
                                    pop.time=seq(1990,(year+t.max-10),10)-year,
                                    is.stochastic=TRUE,
                                    get.births=setup$get.births.here,
                                    rescale.WAIFW=F,
                                    yr.births.per.1000.acrossyears.bysupop=setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                                    space.asdr.object=setup$asdr.object,
                                    year=1980,
                                    EXt0=EXt0,
                                    time.specific.MR1cov =  rbind(setup$RCV1.coverage.1980to2100, setup$RCV1.coverage.1980to2100 *0.9),
                                    time.specific.MR2cov = rbind(rep(0, length(setup$RCV1.coverage.1980to2100)), rep(0, length(setup$RCV1.coverage.1980to2100))),
                                    time.specific.SIAcov =  rbind(setup$rubellaSIA.coverage.1980to2100, setup$rubellaSIA.coverage.1980to2100*0.5),
                                    time.specific.min.age.MR1 = rep(1, length(setup$RCV1.coverage.1980to2100)),
                                    time.specific.max.age.MR1 = rep(24, length(setup$RCV1.coverage.1980to2100)),
                                    time.specific.min.age.MR2 = rep(25, length(setup$RCV1.coverage.1980to2100)),
                                    time.specific.max.age.MR2 = rep(36, length(setup$RCV1.coverage.1980to2100)),
                                    time.specific.min.age.SIA = setup$age.min.sia.rubella,
                                    time.specific.max.age.SIA = setup$age.max.sia.rubella,
                                    obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24), #get.vcdf.uniform(12, 23)
                                    obj.vcdf.MR2 = get.vcdf.normal(25, 36), #get.vcdf.uniform(24, 35)
                                    obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                                    sia.timing.in.year = (1/12),
                                    MR1MR2correlation = NULL,
                                    MR1SIAcorrelation = FALSE,
                                    MR2SIAcorrelation = FALSE,
                                    SIAinacc = TRUE,
                                    prop.inacc = setup$inaccessible.prop.1980.to.2100,
                                    SIAinefficient = TRUE)
print("vaccine experiment done")
plot(rcv_exp@result)
