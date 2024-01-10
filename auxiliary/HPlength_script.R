

setup <- setupCountry.Nov2023(country="Samoa")

year <- 2022
t.max <- 40
year.index.1980 <- which(1980:2100==year)
flat.WAIFW=F
country.specific.WAIFW=T
R0=18

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
                                R0 = R0,
                                t.max = 20,
                                get.births=setup$get.births.here,
                                seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                                flat.WAIFW= flat.WAIFW,
                                country.specific.WAIFW= country.specific.WAIFW,
                                vynnycky.waifw=FALSE,
                                vynnycky.waifw.betas = c(2,1),
                                asdr.object=setup$asdr.object,
                                year=year,
                                use_montagu_demog=FALSE,
                                routine.vac=0,
                                routine.vac.age.index=12,
                                starting.prop.immune = c(rep(0.95, length(9:260)), rep(1, length(261:321))),
                                starting.prop.immune.ages.in.months = c(9:240, seq(252,1212,12)),
                                max.immunity = 1)
#print("part 1 done")

#Remove any infected individuals
EXt0@state.t0[EXt0@trans@i.inds] <- 0
EXt0@trans@introduction.rate <- rep(0, length(EXt0@trans@introduction.rate))
#Check starting immunity correct
#EXt0@state.t0[EXt0@trans@r.inds] /  (EXt0@state.t0[EXt0@trans@r.inds] + EXt0@state.t0[EXt0@trans@s.inds])

# I am revising MR1 and MR2 so that all future years (i.e., >=2023 is the max of 2018:2022, given the pandemic)
setup$MCV1.coverage.1980to2100 <- c(setup$MCV1.coverage.1980to2100[1:year.index.1980],
                                    rep(max(setup$MCV1.coverage.1980to2100[(year.index.1980-4):year.index.1980]), length((year.index.1980+1):length(setup$MCV1.coverage.1980to2100))))
setup$MCV2.coverage.1980to2100 <- c(setup$MCV2.coverage.1980to2100[1:year.index.1980],
                                    rep(max(setup$MCV2.coverage.1980to2100[(year.index.1980-4):year.index.1980]), length((year.index.1980+1):length(setup$MCV2.coverage.1980to2100))))
#setup$MCV1.coverage.1980to2100[setup$MCV1.coverage.1980to2100==1] <- 0.99
#setup$MCV2.coverage.1980to2100[setup$MCV2.coverage.1980to2100==1] <- 0.99

# Updating vaccination coverage is an improvement scenario
if (s==2) {
  index <- 2022-1980+1
  setup$MCV2.coverage.1980to2100[index:(index+t.max)] <- min(setup$MCV1.coverage.1980to2100[index],
                                                             max(0.80, setup$MCV2.coverage.1980to2100[index]))
}
if (s==3) {
  index <- 2022-1980+1
  setup$MCV1.coverage.1980to2100[index:(index+t.max)] <- max(0.9,setup$MCV1.coverage.1980to2100[index])
}
if (s==4) {
  index <- 2022-1980+1
  setup$MCV1.coverage.1980to2100[index:(index+t.max)] <- max(0.9,setup$MCV1.coverage.1980to2100[index])
  setup$MCV2.coverage.1980to2100[index:(index+t.max)] <- min(setup$MCV1.coverage.1980to2100[index],
                                                             max(0.80, setup$MCV2.coverage.1980to2100[index]))
}

tmp <- EX.Country.part2(uncode = setup$uncode,
                        generation.time = 0.5,
                        pop.rescale=setup$pop.total.1950.2100[seq(2030,(year+t.max-10),10)-1950+1],
                        pop.time=seq(2030,(year+t.max-10),10)-year,
                        is.stochastic=FALSE,
                        get.births=setup$get.births.here,
                        t.max = t.max,
                        rescale.WAIFW=F,
                        yr.births.per.1000.acrossyears=setup$cbr.1950.2100[(year-1950+1):((year-1950+1)+t.max)],
                        asdr.object=setup$asdr.object,
                        year=year,
                        EXt0=EXt0,
                        time.specific.MR1cov =  setup$MCV1.coverage.1980to2100[year.index.1980:(year.index.1980+t.max)],
                        time.specific.MR2cov = setup$MCV2.coverage.1980to2100[year.index.1980:(year.index.1980+t.max)],
                        time.specific.SIAcov =  rep(0, length(year:(year+t.max))),
                        time.specific.min.age.MR1 = rep(9, length(year.index.1980:(year.index.1980+t.max))),
                        time.specific.max.age.MR1 = rep(24, length(year.index.1980:(year.index.1980+t.max))),
                        time.specific.min.age.MR2 = rep(25, length(year.index.1980:(year.index.1980+t.max))),
                        time.specific.max.age.MR2 = rep(36, length(year.index.1980:(year.index.1980+t.max))),
                        time.specific.min.age.SIA = rep(9, length(year:(year+t.max))),
                        time.specific.max.age.SIA = rep(59, length(year:(year+t.max))),
                        obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24), #get.vcdf.uniform(12, 23)
                        obj.vcdf.MR2 = get.vcdf.normal(25, 36), #get.vcdf.uniform(24, 35)
                        obj.prob.vsucc = pvacsuccess(EXt0@trans@age.class, get.boulianne.vsucc()),
                        sia.timing.in.year = (3/12),
                        MR1MR2correlation = TRUE,
                        MR1SIAcorrelation = FALSE,
                        MR2SIAcorrelation = FALSE,
                        SIAinacc = FALSE,
                        prop.inacc = rep(0, length(1:t.max)),
                        SIAinefficient = FALSE)



uncode=setup$uncode
generation.time = 0.5
age.classes = c(1:240, seq(252,1212,12))
maternal.decay.rt = 0.45 #based on Leuridan + others wanes b/w 3-9 months
exponent = 0.97
frequency.dep=TRUE
yr.births.per.1000.acrossyears=rep(setup$cbr.1950.2100[year-1950+1],20)
intro.rate=1/24/320
tot.pop=NULL #it will pull a pop automatically from pop distribution if NULL
targeted.intro=FALSE
t.max = 20
get.births=setup$get.births.here
seasonal.amp = 0.15 #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
flat.WAIFW= flat.WAIFW
country.specific.WAIFW= country.specific.WAIFW
vynnycky.waifw=FALSE
vynnycky.waifw.betas = c(2,1)
asdr.object=setup$asdr.object
year=year
use_montagu_demog=FALSE
routine.vac=0
routine.vac.age.index=12
starting.prop.immune = c(rep(0.95, length(9:260)), rep(1, length(261:321)))
starting.prop.immune.ages.in.months = c(9:240, seq(252,1212,12))
max.immunity = 1









uncode = setup$uncode
generation.time = 0.5
pop.rescale=setup$pop.total.1950.2100[seq(2030,(year+t.max-10),10)-1950+1]
pop.time=seq(2030,(year+t.max-10),10)-year
is.stochastic=FALSE
get.births=setup$get.births.here
t.max = t.max
rescale.WAIFW=F
yr.births.per.1000.acrossyears=setup$cbr.1950.2100[(year-1950+1):((year-1950+1)+t.max)]
asdr.object=setup$asdr.object
EXt0=EXt0
time.specific.MR1cov =  setup$MCV1.coverage.1980to2100[year.index.1980:(year.index.1980+t.max)]
time.specific.MR2cov = setup$MCV2.coverage.1980to2100[year.index.1980:(year.index.1980+t.max)]
time.specific.SIAcov =  rep(0, length(year:(year+t.max)))
time.specific.min.age.MR1 = rep(9, length(year.index.1980:(year.index.1980+t.max)))
time.specific.max.age.MR1 = rep(24, length(year.index.1980:(year.index.1980+t.max)))
time.specific.min.age.MR2 = rep(25, length(year.index.1980:(year.index.1980+t.max)))
time.specific.max.age.MR2 = rep(36, length(year.index.1980:(year.index.1980+t.max)))
time.specific.min.age.SIA = rep(9, length(year:(year+t.max)))
time.specific.max.age.SIA = rep(59, length(year:(year+t.max)))
obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24) #get.vcdf.uniform(12, 23)
obj.vcdf.MR2 = get.vcdf.normal(25, 36) #get.vcdf.uniform(24, 35)
obj.prob.vsucc = pvacsuccess(EXt0@trans@age.class, get.boulianne.vsucc())
sia.timing.in.year = (3/12)
MR1MR2correlation = TRUE
MR1SIAcorrelation = FALSE
MR2SIAcorrelation = FALSE
SIAinacc = FALSE
prop.inacc = rep(0, length(1:t.max))
SIAinefficient = FALSE

