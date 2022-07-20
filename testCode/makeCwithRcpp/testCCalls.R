# test C calls
getwd()
setwd("C:/Users/lhoskovec/Documents/Repos/MRTransmissionModel/makeCwithRcpp")

#######################################
### Three ways to call C/C++ from R ### ----------------------------------------
#######################################

# 1) the .C function
## limitations: all arguments must be pointers, must return void

# R CMD SHLIB pointersOnly.c
dyn.load("pointersOnly.dll") # worked
# use .C to call the C function


# 2) the .Call function
## limitations: must return void or SEXP
## requires header files

# R CMD SHLIB cFunctions.c
dyn.load("cFunctions.dll")

