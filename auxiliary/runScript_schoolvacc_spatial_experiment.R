## Amy Winter
## 24 Janauary 2024
## This is the code used to conduct deterministic runs to Gavi Preliminary Report
## Note: I couldn't get parallel processing to work on my Mac Studio b/c couldn't download Rmpi


## Part 1 - Get Initial Population ####
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.RData")
year=2016
R0=12 #18
EXt0 <- Spatial.EX.Country.part1(year = 2016,
                         uncode=setup$uncode,
                         iso3code="ZMB",
                         generation.time = 0.5,
                         age.classes = c(1:240, seq(252,960,12)),
                         maternal.decay.rt = 0.45, #based on Leuridan + others wanes b/w 3-9 months
                         exponent = 0.97,
                         frequency.dep=TRUE,
                         yr.births.per.1000.bysubpop <- setup$cbr.1950.2100[,year-1950+1],
                         targeted.intro=FALSE,
                         intro.rate=1/24/320,
                         tot.subpop=NULL, #it will pull a pop automatically from pop distribution if NULL
                         R0 = R0,
                         t.max = 20,
                         get.births=setup$get.births.here,
                         seasonal.amp = 0.15, #metcalf et al 2012 used 0.35, age-structured metcalf paper used 0.2, ferrari 2007 nature 0.2 to 0.6 for measles
                         flat.WAIFW=FALSE,
                         country.specific.WAIFW=TRUE,
                         vynnycky.waifw=FALSE,
                         vynnycky.waifw.betas = c(2,1),
                         space.asdr.object=setup$asdr.object,
                         routine.vac=0,
                         routine.vac.age.index=12,
                         n.subpops=length(setup$district.order),
                         coupling = as.matrix(setup$coupling),
                         starting.prop.immune = setup$seroprev.byageclasses.subpop %>%
                           select(age, mean.prob, district) %>%
                           pivot_wider(names_from=district, values_from=mean.prob) %>%
                           select(-age) %>% as.data.frame(),
                         starting.prop.immune.ages.in.months = c(unique(setup$seroprev.byageclasses.subpop$age)),
                         max.immunity = 0.98)

# Saving it as hacked b/c hacked survival matrices
# Note: assumed MR1 and MR2 to be dependent (hard coded in the spatial experiments in run.R)
#save(EXt0,  file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_hacked.RData")
#save(EXt0,  file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_12_hacked.RData")


## Part 2 - Run the 10 scnearios ####
#library(foreach)
#library(doMPI)
#library(doRNG)
#cl <- startMPIcluster()
#print(clusterSize(cl)) # this just tells you how many you've got
#registerDoMPI(cl)
#foreach(s=c(1:10))  %dopar%  {

#load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_hacked.RData")
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_12_hacked.RData")
t.max = 30
year = 2016

#Routine Only
EXt0@state.t0[EXt0@trans@i.inds] <- 0  #remove all I
EXt0@trans@introduction.rate <- rep(0, length(EXt0@trans@introduction.rate))
prob.vsucc <- new("prob.vsucc.byage",
                  prob.vsucc=rep(1, length(1:(20*12))),
                  ages=1:(20*12))
#if (s==1) {
mr1 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = matrix(rep(0,116*30), nrow=116, ncol=30),
                               time.specific.SIAcov = matrix(rep(0, 116*30), nrow=116, ncol=30),
                               time.specific.schoolvacc.cov=NULL,
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = rep(0, 30),
                               time.specific.max.age.SIA = rep(0, 30),
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = prob.vsucc,
                               sia.timing.in.year = (10/12),
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
save(mr1, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/mr1_R012.RData")
mr2 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                                generation.time = 0.5,
                                pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                                pop.time = seq(2016,(year+t.max-10),10)-year,
                                is.stochastic = FALSE,
                                get.births = setup$get.births.here,
                                t.max = t.max,
                                rescale.WAIFW = F,
                                yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                                space.asdr.object = setup$asdr.object,
                                year = year,
                                EXt0 = EXt0,
                                time.specific.MR1cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                                time.specific.MR2cov = matrix(rep(0,116*30), nrow=116, ncol=30),
                                time.specific.SIAcov = matrix(rep(0, 116*30), nrow=116, ncol=30),
                                time.specific.schoolvacc.cov=NULL,
                                time.specific.min.age.MR1 = rep(25,30),
                                time.specific.max.age.MR1 = rep(36,30),
                                time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                                time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                                time.specific.min.age.SIA = rep(0, 30),
                                time.specific.max.age.SIA = rep(0, 30),
                                obj.vcdf.MR1 = get.vcdf.normal(25, 36),
                                obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                                obj.prob.vsucc = prob.vsucc,
                                sia.timing.in.year = (10/12),
                                schoolvacc.timing.in.year = (9/12),
                                MR1MR2correlation=NULL)
save(mr2, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/mr2_R012.RData")
#}

#Baseline
#if (s==2) {
bs <- Spatial.EX.Country.part2(uncode = setup$uncode,
                                generation.time = 0.5,
                                pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                                pop.time = seq(2016,(year+t.max-10),10)-year,
                                is.stochastic = FALSE,
                                get.births = setup$get.births.here,
                                t.max = t.max,
                                rescale.WAIFW = F,
                                yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                                space.asdr.object = setup$asdr.object,
                                year = year,
                                EXt0 = EXt0,
                                time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                                time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                                time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                                time.specific.schoolvacc.cov=NULL,
                                time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                                time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                                time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                                time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                                time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                                time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                                obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                                obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                                obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                                sia.timing.in.year = (10/12),
                                schoolvacc.timing.in.year = (9/12),
                                MR1MR2correlation=NULL)
save(bs, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/bs_12.RData")
#plot(bs@result)
#plot(colSums(bs@result@.Data[bs@result@i.inds,])/colSums(bs@result@.Data))
#plot(colSums(bs@result@.Data[bs@result@v.inds,])/colSums(bs@result@.Data))
#plot(colSums(bs@result@.Data[bs@result@r.inds,])/colSums(bs@result@.Data))
#}

#if (s==3) {
# Scenario 1
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (4*12)+11
s1 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.schoolvacc.cov=NULL,
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
save(s1, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s1_12.RData")
#plot(s1@result)
#}

#if (s==4) {
# Scenario 2
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (9*12)+11
s2 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.schoolvacc.cov=NULL,
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
save(s2, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s2.RData")
#plot(s2@result)
#colSums(s2@result@.Data[s2@result@i.inds,])/colSums(s2@result@.Data)
#colSums(s2@result@.Data[s2@result@i.inds,])/colSums(s1@result@.Data[s1@result@i.inds,])
#colSums(s2@result@.Data[s2@result@v.inds,])/colSums(s1@result@.Data[s1@result@v.inds,]) #this is not good, no differences!
#colSums(s2@result@.Data[s2@result@v.inds,])/colSums(bs@result@.Data[bs@result@v.inds,]) #that is good, only see changes by year 8.83!
#}

#if (s==5) {
# Scenario 3
#load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.RData")
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.ForcedSchoolRate.RData")
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(0.8, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
s3 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
#save(s3, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s3.RData")
#save(s3, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s3_12.RData")
save(s3, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s3_12_forced.RData")

#plot(s3@result)
#plot(colSums(s3@result@.Data[bs@result@i.inds,])/colSums(s3@result@.Data))
#plot(colSums(s3@result@.Data[bs@result@v.inds,])/colSums(s3@result@.Data))
#plot(colSums(s3@result@.Data[bs@result@r.inds,])/colSums(s3@result@.Data))
#}

#if (s==6) {
# Scenario 4
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(0.8, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (4*12)+11
s4 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
#save(s4, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s4_12.RData")
save(s4, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s4_12_forced.RData")

#}

#if (s==7) {
# Scenario 5
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(0.8, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (9*12)+11
s5 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
save(s5, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s5.RData")
#}

#if (s==8) {
# Scenario 6
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.RData")
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(1, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
s6 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
#save(s6, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s6_12.RData")
save(s6, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s6_12_forced.RData")
plot(s6@result)
#}

#if (s==9) {
# Scenario 7
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(1, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (4*12)+11
s7 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
#save(s7, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s7_12.RData")
save(s7, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s7_12_forced.RData")

#}

#if (s==10) {
# Scenario 8
time.specific.schoolvacc.cov = matrix(c(rep(0,length(year:2023)*116), rep(1, length(2024:(year+t.max))*116)),
                                      nrow=116, ncol=length(year:(year+t.max)))
time.specific.min.age.schoolvacc = rep(3*12, length(year:(year+t.max)))
time.specific.max.age.schoolvacc = rep(15*12, length(year:(year+t.max)))
setup$measlesSIA.coverage.1980to2100[,which((1980:2100)%in%seq(2024,2100,4))] <- 0.90
setup$time.specific.min.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- 9
setup$time.specific.max.age.SIA[which((1980:2100)%in%seq(2024,2100,4))] <- (9*12)+11
s8 <- Spatial.EX.Country.part2(uncode = setup$uncode,
                               generation.time = 0.5,
                               pop.rescale = setup$pop.total.1950.2100[,seq(2016,(year+t.max-10),10)-1950+1],
                               pop.time = seq(2016,(year+t.max-10),10)-year,
                               is.stochastic = FALSE,
                               get.births = setup$get.births.here,
                               t.max = t.max,
                               rescale.WAIFW = F,
                               yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)],
                               space.asdr.object = setup$asdr.object,
                               year = year,
                               EXt0 = EXt0,
                               time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)],
                               time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)],
                               time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)],
                               obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24),
                               obj.vcdf.MR2 = get.vcdf.normal(25, 36),
                               obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                               sia.timing.in.year = (10/12),
                               time.specific.schoolvacc.cov = time.specific.schoolvacc.cov,
                               time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc,
                               time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc,
                               list.obj.vcdf.schoolenroll = setup$list.obj.vcdf.schoolvacc,
                               schoolvacc.timing.in.year = (9/12),
                               MR1MR2correlation=NULL)
save(s8, file="/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/s8.RData")
#}

#}
