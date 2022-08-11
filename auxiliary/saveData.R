save(contact_all, file = "data/contact_all.rda",
     compression_level = 9)


save(country_codes, file = "data/country_codes.rda",
     compression_level = 9)


save(DHS_ageMCV1_Data, file = "data/DHS_ageMCV1_Data.rda",
     compression_level = 9)


save(DHS_ageMCV1_Data_AFRO, file = "data/DHS_ageMCV1_Data_AFRO.rda",
     compression_level = 9)


save(DHS_ageMCV1_Data_EMRO, file = "data/DHS_ageMCV1_Data_EMRO.rda",
     compression_level = 9)


save(DHS_ageMCV1_Data_SEARO, file = "data/DHS_ageMCV1_Data_SEARO.rda",
     compression_level = 9)


tools::resaveRdaFiles(paths = "data/DHS_ageMCV1_Data.rda", compress = "bzip2")
tools::resaveRdaFiles(paths = "data/DHS_ageMCV1_Data_AFRO.rda", compress = "bzip2")
tools::resaveRdaFiles(paths = "data/DHS_ageMCV1_Data_EMRO.rda", compress = "bzip2")
tools::resaveRdaFiles(paths = "data/DHS_ageMCV1_Data_SEARO.rda", compress = "bzip2")

tools::resaveRdaFiles(paths = "data/contact_all.rda", compress = "xz")




dtp1_estimates <- read.csv(file = "data/dtp1_estimates.csv")
save(dtp1_estimates, file = "data/dtp1_estimates.rda")

Encoding(dtp1_estimates$Cname) <- "latin1"
save(dtp1_estimates, file = "data/dtp1_estimates.rda")

#tools::resaveRdaFiles(paths = "data/dtp1_estimates.rda", compress = "")

polymodRaw <- read.table(file = "data/polymodRaw.txt", sep = ",", header = TRUE)
save(polymodRaw, file = "data/polymodRaw.rda")



