#impotring libraries and code####

setwd("C:/Users/konkons/Desktop/ENF/ENFdata/Konstantinos") # cHANGE THE PATH
source('./exscript.R') #load data,and libaries
source('./Helpfunctions.R') #load functions
source('./Angles_Lengths.R')
source('./Angular_Distributions.R')
library(fitdistrplus)
require(extraDistr)

#FUNCTION THAT FITS DISTRIBUTIONS TO DATA #####
subject_param_estimates = function(list_element_FB, list_element_EP){
  L = list_element_FB$Length
  L1 = list_element_EP$Length
  fit_L= fitdist(L, 'gamma')
  shape = fit_L$estimate[1]
  rate = fit_L$estimate[2]
  fit_L1 = fitdist(L1, 'gamma')
  shape_ep = fit_L1$estimate[1]
  rate_ep = fit_L1$estimate[2]
  shpptsep<- spherepoints(list_element_FB$Colatidute ,list_element_FB$XYangle)
  R = sqrt(sum(shpptsep$X)^2+sum(shpptsep$Y)^2)/length(shpptsep$X)
  KAPPA= (R*2-R^3)/(1-R^2)
  beta = betamle(100,list_element_EP$Colatidute,1)
  res= data.frame(shape = shape , rate=rate , shape_ep = shape_ep , rate_ep = rate_ep,
                  kappa= KAPPA , beta = beta , SID=list_element_FB$SID[1])
}



#NORMAL GROUP PARAMTERS ##############
Group = 'Normal'
FB = Get_Firstbranch(Group)
EP = Get_endbranch(Group)
EP1 = subset(EP, EP$Zangle>0)
SID_Normal = unique(FB$SID)
Subject_information_FB = lapply(seq_along(SID_Normal), function(x) subset(FB,FB$SID==SID_Normal[x]))
Subject_information_EP = lapply(seq_along(SID_Normal), function(x) subset(EP1,EP1$SID==SID_Normal[x]))

PARAM_ESTIMATES = lapply(seq_along(Subject_information_FB),
        function(x) subject_param_estimates(Subject_information_FB[[x]],Subject_information_EP[[x]]))
shape = unlist(lapply(PARAM_ESTIMATES, function(x) x$shape))
shape_ep = unlist(lapply(PARAM_ESTIMATES, function(x) x$shape_ep))
rate = unlist(lapply(PARAM_ESTIMATES, function(x) x$rate))
rate_ep = unlist(lapply(PARAM_ESTIMATES, function(x) x$rate_ep))
beta = unlist(lapply(PARAM_ESTIMATES, function(x) x$beta))
kappa = unlist(lapply(PARAM_ESTIMATES, function(x) x$kappa))
params = data.frame(shape = shape , rate=rate , shape_ep = shape_ep , rate_ep = rate_ep,
                    kappa= kappa , beta = beta , SID = unique(FB$SID) )



##MILD GROUP PARAMETERS  ############
Group = 'Mild'
FB_m = Get_Firstbranch(Group)
EP_m = Get_endbranch(Group)
EP1_m = subset(EP_m, EP_m$Zangle>0)
SID_Mild = unique(FB_m$SID)
Subject_information_FB_m = lapply(seq_along(SID_Mild), function(x) subset(FB_m,FB_m$SID==SID_Mild[x]))
Subject_information_EP_m = lapply(seq_along(SID_Mild), function(x) subset(EP1_m,EP1_m$SID==SID_Mild[x]))

PARAM_ESTIMATES = lapply(seq_along(Subject_information_FB_m),
                         function(x) subject_param_estimates(Subject_information_FB_m[[x]],
                                                             Subject_information_EP_m[[x]]))
shape_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$shape))
shape_ep_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$shape_ep))
rate_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$rate))
rate_ep_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$rate_ep))
beta_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$beta))
kappa_m = unlist(lapply(PARAM_ESTIMATES, function(x) x$kappa))
params_m = data.frame(shape = shape_m , rate=rate_m , shape_ep = shape_ep_m , rate_ep = rate_ep_m,
                    kappa= kappa_m , beta = beta_m , SID = unique(FB_m$SID) )







subject_param_estimates = function(list_element_FB, list_element_EP){
  L = list_element_FB$Length
  L1 = list_element_EP$Length
  fit_L= fitdist(L, 'gamma')
  shape = fit_L$estimate[1]
  rate = fit_L$estimate[2]
  fit_L1 = fitdist(L1, 'gamma')
  shape_ep = fit_L1$estimate[1]
  rate_ep = fit_L1$estimate[2]
  shpptsep<- spherepoints(list_element_EP$Colatidute ,list_element_EP$XYangle)
  R = sqrt(sum(shpptsep$X)^2+sum(shpptsep$Y)^2+sum(shpptsep$Z)^2)/length(shpptsep$X)
  KAPPA= (R*3-R^3)/(1-R^2)
  beta = betamle(100,list_element_EP$Colatidute,1)
  res= data.frame(shape = shape , rate=rate , shape_ep = shape_ep , rate_ep = rate_ep,
                  kappa= KAPPA , beta = beta , SID=list_element_FB$SID[1])
}
