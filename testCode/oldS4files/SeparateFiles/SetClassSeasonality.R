#generic virtual class for seasonal forcing
setClass("seasonal.force")

setClass("seasonal.cosine",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                        offset = "numeric", #offset from t=0 for nadir
                        period = "numeric"), #number of peaks per year
         prototype = list(amplitude=.2, offset = 0, period = 1),
         contains="seasonal.force"
)

setClass("seasonal.periodic.bspline",
         slots = list(parameters = "numeric",
                        nbasis = "numeric",
                        degree = "numeric",
                        period = "numeric"
         ),
         prototype = list(parameters= c(0.5,1,1), nbasis=3, degree=3, period=1),
         contains="seasonal.force",
)

setClass("seasonal.age.specific",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                        offset = "numeric", #offset from t=0 for nadir
                        period = "numeric", #number of peaks per year
                        age.class = "numeric",#age classes in the WAIFW
                        chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(amplitude=.2, offset = 0, period = 1, age.class=c(1:60,60*12),chosen.age=c(50,60)),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.piecewise.constant",
         slots = list(parameters = "numeric",
                        time.in.year = "numeric"),
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.age.specific.piecewise.constant",
         slots = list(parameters = "numeric",
                        time.in.year = "numeric",
                        age.class = "numeric",#age classes in the WAIFW
                        chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)
