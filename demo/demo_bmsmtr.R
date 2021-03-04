#load packages
library(rstan)

#load model and utility libraries
require("bmsr")

#set random seed
set.seed(101)

#run demo for bmsmtr
demo_bmsmtr(dY = 3)
