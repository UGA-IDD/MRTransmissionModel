## Amy Winter
## Issue - space.create.ID.transition.MSIRV() takes way too long
## This is an attempt to vectorize the nested loops. I was successful at getting diagnoal values and indeces, when in line 158 it wasn't allowing me to input
## due to the type of object.  Also, I never did it for the off (i.e, bottom) diagonals


#  Vectorized calculation for diagonal elements
tmp.template <- template.mtrx
tmp.template[2,1] <- 0

#aging rate as a vector - each age-specific rate listed 25 times for each district
aging.rate.vec <- rep(rep((1-aging.rate),each=length(tmp.template)),nrow(survival.rate))
#survival rate as a vector - each age-specific rate listed 25 times for each district
survival.rate.vec <- rep(as.vector(t(survival.rate)), each=length(tmp.template))
#tmp.template as vector - each matrix listed for each age and subpop
tmp.template.vec <- rep(as.vector(tmp.template),length(survival.rate))

#multiply all for diagonal elements - order is 5x5 matrix for district 1 all age groups, then district 2 all age groups, etc.
diag_values <- aging.rate.vec*survival.rate.vec*tmp.template.vec

#vectorized index calculation for diagonal elements
vec_start <- rep(seq(1, 1496, 5), each=5)
vec_end <- rep(seq(5, 1500, 5), each=5)
inds_diag_d1 <-rep(c(1500*0,(174000+(1500*0))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end)) #district 1
inds_diag_d2 <-rep(c(1500*1,(174000+(1500*1))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end)) #district 2
inds_diag_d3 <-rep(c(1500*2,(174000+(1500*2))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end)) #district 3
inds_diag_d4 <-rep(c(1500*3,(174000+(1500*3))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d5 <-rep(c(1500*4,(174000+(1500*4))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d6 <-rep(c(1500*5,(174000+(1500*5))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d7 <-rep(c(1500*6,(174000+(1500*6))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d8 <-rep(c(1500*7,(174000+(1500*7))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d9 <-rep(c(1500*8,(174000+(1500*8))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d10 <-rep(c(1500*9,(174000+(1500*9))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d11 <-rep(c(1500*10,(174000+(1500*10))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d12 <-rep(c(1500*11,(174000+(1500*11))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d13 <-rep(c(1500*12,(174000+(1500*12))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d14 <-rep(c(1500*13,(174000+(1500*13))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d15 <-rep(c(1500*14,(174000+(1500*14))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d16 <-rep(c(1500*15,(174000+(1500*15))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d17 <-rep(c(1500*16,(174000+(1500*16))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d18 <-rep(c(1500*17,(174000+(1500*17))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d19 <-rep(c(1500*18,(174000+(1500*18))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d20 <-rep(c(1500*19,(174000+(1500*19))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d21 <-rep(c(1500*20,(174000+(1500*20))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d22 <-rep(c(1500*21,(174000+(1500*21))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d23 <-rep(c(1500*22,(174000+(1500*22))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d24 <-rep(c(1500*23,(174000+(1500*23))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d25 <-rep(c(1500*24,(174000+(1500*24))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d26 <-rep(c(1500*25,(174000+(1500*25))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d27 <-rep(c(1500*26,(174000+(1500*26))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d28 <-rep(c(1500*27,(174000+(1500*27))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d29 <-rep(c(1500*28,(174000+(1500*28))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d30 <-rep(c(1500*29,(174000+(1500*29))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d31 <-rep(c(1500*30,(174000+(1500*30))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d32 <-rep(c(1500*31,(174000+(1500*31))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d33 <-rep(c(1500*32,(174000+(1500*32))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d34 <-rep(c(1500*33,(174000+(1500*33))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d35 <-rep(c(1500*34,(174000+(1500*34))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d36 <-rep(c(1500*35,(174000+(1500*35))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d37 <-rep(c(1500*36,(174000+(1500*36))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d38 <-rep(c(1500*37,(174000+(1500*37))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d39 <-rep(c(1500*38,(174000+(1500*38))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d40 <-rep(c(1500*39,(174000+(1500*39))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d41 <-rep(c(1500*40,(174000+(1500*40))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d42 <-rep(c(1500*41,(174000+(1500*41))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d43 <-rep(c(1500*42,(174000+(1500*42))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d44 <-rep(c(1500*43,(174000+(1500*43))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d45 <-rep(c(1500*44,(174000+(1500*44))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d46 <-rep(c(1500*45,(174000+(1500*45))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d47 <-rep(c(1500*46,(174000+(1500*46))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d48 <-rep(c(1500*47,(174000+(1500*47))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d49 <-rep(c(1500*48,(174000+(1500*48))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d50 <-rep(c(1500*49,(174000+(1500*49))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d51 <-rep(c(1500*50,(174000+(1500*50))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d52 <-rep(c(1500*51,(174000+(1500*51))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d53 <-rep(c(1500*52,(174000+(1500*52))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d54 <-rep(c(1500*53,(174000+(1500*53))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d55 <-rep(c(1500*54,(174000+(1500*54))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d56 <-rep(c(1500*55,(174000+(1500*55))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d57 <-rep(c(1500*56,(174000+(1500*56))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d58 <-rep(c(1500*57,(174000+(1500*57))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d59 <-rep(c(1500*58,(174000+(1500*58))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d60 <-rep(c(1500*59,(174000+(1500*59))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d61 <-rep(c(1500*60,(174000+(1500*60))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d62 <-rep(c(1500*61,(174000+(1500*61))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d63 <-rep(c(1500*62,(174000+(1500*62))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d64 <-rep(c(1500*63,(174000+(1500*63))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d65 <-rep(c(1500*64,(174000+(1500*64))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d66 <-rep(c(1500*65,(174000+(1500*65))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d67 <-rep(c(1500*66,(174000+(1500*66))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d68 <-rep(c(1500*67,(174000+(1500*67))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d69 <-rep(c(1500*68,(174000+(1500*68))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d70 <-rep(c(1500*69,(174000+(1500*69))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d71 <-rep(c(1500*70,(174000+(1500*70))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d72 <-rep(c(1500*71,(174000+(1500*71))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d73 <-rep(c(1500*72,(174000+(1500*72))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d74 <-rep(c(1500*73,(174000+(1500*73))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d75 <-rep(c(1500*74,(174000+(1500*74))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d76 <-rep(c(1500*75,(174000+(1500*75))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d77 <-rep(c(1500*76,(174000+(1500*76))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d78 <-rep(c(1500*77,(174000+(1500*77))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d79 <-rep(c(1500*78,(174000+(1500*78))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d80 <-rep(c(1500*79,(174000+(1500*79))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d81 <-rep(c(1500*80,(174000+(1500*80))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d82 <-rep(c(1500*81,(174000+(1500*81))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d83 <-rep(c(1500*82,(174000+(1500*82))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d84 <-rep(c(1500*83,(174000+(1500*83))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d85 <-rep(c(1500*84,(174000+(1500*84))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d86 <-rep(c(1500*85,(174000+(1500*85))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d87 <-rep(c(1500*86,(174000+(1500*86))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d88 <-rep(c(1500*87,(174000+(1500*87))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d89 <-rep(c(1500*88,(174000+(1500*88))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d90 <-rep(c(1500*89,(174000+(1500*89))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d91 <-rep(c(1500*90,(174000+(1500*90))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d92 <-rep(c(1500*91,(174000+(1500*91))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d93 <-rep(c(1500*91,(174000+(1500*92))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d94 <-rep(c(1500*93,(174000+(1500*93))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d95 <-rep(c(1500*94,(174000+(1500*94))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d96 <-rep(c(1500*95,(174000+(1500*95))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d97 <-rep(c(1500*96,(174000+(1500*96))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d98 <-rep(c(1500*97,(174000+(1500*97))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d99 <-rep(c(1500*98,(174000+(1500*98))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d100 <-rep(c(1500*99,(174000+(1500*99))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d101 <-rep(c(1500*100,(174000+(1500*100))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d102 <-rep(c(1500*101,(174000+(1500*101))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d103 <-rep(c(1500*102,(174000+(1500*102))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d104 <-rep(c(1500*103,(174000+(1500*103))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d105 <-rep(c(1500*104,(174000+(1500*104))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d106 <-rep(c(1500*105,(174000+(1500*105))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d107 <-rep(c(1500*106,(174000+(1500*106))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d108 <-rep(c(1500*107,(174000+(1500*107))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d109 <-rep(c(1500*108,(174000+(1500*108))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d110 <-rep(c(1500*109,(174000+(1500*109))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d111 <-rep(c(1500*110,(174000+(1500*110))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d112 <-rep(c(1500*111,(174000+(1500*111))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d113 <-rep(c(1500*112,(174000+(1500*112))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d114 <-rep(c(1500*113,(174000+(1500*113))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d115 <-rep(c(1500*114,(174000+(1500*114))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag_d116 <-rep(c(1500*115,(174000+(1500*115))*1:1499),each=5) + as.vector(mapply(seq, vec_start, vec_end))
inds_diag <- c(inds_diag_d1, inds_diag_d2, inds_diag_d3, inds_diag_d4, inds_diag_d5, inds_diag_d6, inds_diag_d7,
               inds_diag_d8, inds_diag_d9, inds_diag_d10, inds_diag_d11, inds_diag_d12, inds_diag_d13, inds_diag_d14,
               inds_diag_d15, inds_diag_d16, inds_diag_d17, inds_diag_d18, inds_diag_d19, inds_diag_d20, inds_diag_d21,
               inds_diag_d22, inds_diag_d23, inds_diag_d24, inds_diag_d25, inds_diag_d26, inds_diag_d27, inds_diag_d28,
               inds_diag_d29, inds_diag_d30, inds_diag_d31, inds_diag_d32, inds_diag_d33, inds_diag_d34, inds_diag_d35,
               inds_diag_d36, inds_diag_d37, inds_diag_d38, inds_diag_d39, inds_diag_d40, inds_diag_d41, inds_diag_d42,
               inds_diag_d43, inds_diag_d44, inds_diag_d45, inds_diag_d46, inds_diag_d47, inds_diag_d48, inds_diag_d49,
               inds_diag_d50, inds_diag_d51, inds_diag_d52, inds_diag_d53, inds_diag_d54, inds_diag_d55, inds_diag_d56,
               inds_diag_d57, inds_diag_d58, inds_diag_d59, inds_diag_d60, inds_diag_d61, inds_diag_d62, inds_diag_d63,
               inds_diag_d64, inds_diag_d65, inds_diag_d66, inds_diag_d67, inds_diag_d68, inds_diag_d69, inds_diag_d70,
               inds_diag_d71, inds_diag_d72, inds_diag_d73, inds_diag_d74, inds_diag_d75, inds_diag_d76, inds_diag_d77,
               inds_diag_d78, inds_diag_d79, inds_diag_d80, inds_diag_d81, inds_diag_d82, inds_diag_d83, inds_diag_d84,
               inds_diag_d85, inds_diag_d86, inds_diag_d87, inds_diag_d88, inds_diag_d89, inds_diag_d90, inds_diag_d91,
               inds_diag_d92, inds_diag_d93, inds_diag_d94, inds_diag_d95, inds_diag_d96, inds_diag_d97, inds_diag_d98,
               inds_diag_d99, inds_diag_d100, inds_diag_d101, inds_diag_d102, inds_diag_d103, inds_diag_d104, inds_diag_d105,
               inds_diag_d106, inds_diag_d107, inds_diag_d108, inds_diag_d109, inds_diag_d110, inds_diag_d111, inds_diag_d112,
               inds_diag_d113, inds_diag_d114, inds_diag_d115, inds_diag_d116)

#fill it in
rc@age.surv.matrix[inds_diag]  <- diag_values




## Here is a way to get diag_values by age group and then by district ####
#aging rate as a vector - each district-specific rate listed 25 times for each age
aging.rate.vec <- rep((1-aging.rate),each=length(tmp.template)*nrow(survival.rate))
#survival rate as a vector - each district-specific rate listed 25 times for each age
survival.rate.vec <- rep(as.vector(survival.rate), each=length(tmp.template))
#tmp.template as vector - each matrix listed for each age and subpop
tmp.template.vec <- rep(as.vector(tmp.template),length(survival.rate))

#multiply all for diagonal elements - order is 5x5 matrix for age group1 for all districts, then agegroup2 for all districts, etc.
diag_values <- aging.rate.vec*survival.rate.vec*tmp.template.vec






