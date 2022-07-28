#### Maternal Immunity ####

#assumes that protection is the same as exponetial decay
setClass("maternal.exp.decay",
         slots = list(decay.rt = "numeric")#the exponential decay rate
)


#assume that there is some threshold age at which
#you are 100% protected if this age or younger
setClass("maternal.thresh",
         slots = list(thresh="numeric"))
