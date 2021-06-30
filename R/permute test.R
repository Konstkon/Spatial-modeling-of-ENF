require(spatstat)
require(spptest)

# Function that performs a permutation test for a marked point proccess
# Data= marked point pattern I.e realbpm that stands for real/base points/marked
# mark_type =  Name of the mark

permutation_test<- function(data,mark_type){ 
  
  rand_permutation = rlabel(data, labels=marks(data), permute=TRUE, nsim=250, drop=TRUE)
  
  if (mark_type=="Max_ord"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate",method="density",
                                                                  kernel="epanechnikov", f=  function(m1,m2) {(m1*m2)})$Max_ord)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $Max_ord$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$Max_ord$r
  }
  if (mark_type=="AvgDiameter_microm"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate", f=  function(m1,m2) {(m1*m2)})$AvgDiameter_microm)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $AvgDiameter_microm$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$AvgDiameter_microm$r
  }
  if (mark_type=="SurfaceArea_microm2"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate", f=  function(m1,m2) {(m1*m2)})$SurfaceArea_microm2)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $SurfaceArea_microm2$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$SurfaceArea_microm2$r
  }
  if (mark_type=="Volume_microm"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate", f=  function(m1,m2) {(m1*m2)})$Volume_microm)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $Volume_microm$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$Volume_microm$r
  }
  
  if (mark_type =="BaseDiameter_microm"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate", f=  function(m1,m2) {(m1*m2)})$BaseDiameter_microm)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $BaseDiameter_microm$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$BaseDiameter_microm$r
  }
  if (mark_type =="Volume"){
    markcorr_list = lapply(rand_permutation, function(x) markcorr(x, correction ="translate", f=  function(m1,m2) {(m1*m2)})$Volume)
    observed = markcorr(data , f=  function(m1,m2) {(m1*m2)} ) $Volume$trans
    r=markcorr(data, f=  function(m1,m2) {(m1*m2)} )$Volume$r
  }
  
  s = zeros(length(r), length(markcorr_list))
  for (k in 1:length(markcorr_list)){
    s[,k] = markcorr_list[[k]]$trans
  }
  
  curve_list <-  list( r= r , obs =  as.vector(observed), sim_m = s)
  curve_set <- create_curve_set(curve_list)
  res <- rank_envelope(curve_set, alpha=0.05)
  plot(res,use_ggplot2=TRUE)
}

#FIGURE 7 CAN BE CREATED WITH THE ABOVE FUNCTION