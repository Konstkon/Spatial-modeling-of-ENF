library(latex2exp)
source('K_overall_help.R')
direction_list=readRDS('./datarelated/model_params/normal_directions_noc.RDA')
direction_list_mild = readRDS('./datarelated/model_params/mild_directions_noc.RDA')

#GETS NOC DIRECTION
get_direction_to_open_space <-function(Basepoints){
  Basep = Basepoints[,c('X','Y')]
  basedist = as.matrix(dist(Basep))
  
  direction = data.frame()
  for ( i in 1:length(basedist[1,])){
    min_dist = sort(basedist[i,])[2]
    pos = which(basedist[i,] == min_dist)
    direction_temp = data.frame(X = Basep[i,]$X - Basep[pos,]$X , Y = Basep[i,]$Y-Basep[pos,]$Y)    
    direction = rbind( direction, direction_temp)
  }
  as.matrix(direction)
  
}   
#HOW TO USE
#direction_list = lapply(seq_along(realep),function(x) cart2pol(get_direction_to_open_space(realbpdf[[x]]))[,1])

Group ='Normal'
if (Group=='Normal'){
  #THIS MATRIX SAVES THE CURVES(OVERALL K) CREATED AT EVERY ITERATION
  Curve_sim = zeros(300,2500)
  for (nsim in 1:2500){
    tic()
    print(paste(nsim))
    #SIMULATE
    Noclike_3d = lapply(seq_along(realbp), function(x) Noc_model(realbpdf[[x]],params,direction_list[[x]],SID[x]))
    #GET ENDPOINTS AND MAKE PPP OBJECTS
    endsim_2d_noc =lapply(Noclike_3d , function (x) subset(x, Endpoint=='TRUE'))
    endsimNOC =lapply(seq_along(endsim_2d_noc) ,  function(x) pp3(endsim_2d_noc[[x]]$X,endsim_2d_noc[[x]]$Y,endsim_2d_noc[[x]]$Z,
                                                                  window = realep3D[[x]]$domain))
    toc()
    Curve_sim[,nsim] = K_overall_function(endsimNOC,samplesid_per_subj,3,TRUE)$K
  }
  #COMPUTE THE L FUNCTION
  curve_sim_L=zeros(300,2500)
  for (nsim in 1:2500){
    curve_sim_L[,nsim] = (Curve_sim[,nsim]/(4*pi/3))^0.33-r
  }
  
  #OBSERVED CURVE
  observed = (K_overall_function(realep3D, samplesid_per_subj,3,TRUE)$K/(4*pi/3))^0.33-K_overall_function(realep3D,samplesid_per_subj,3,TRUE)$r
  #RADIUS
  r = K_overall_function(realep3D,samplesid_per_subj,3,TRUE)$r
  #CREATE A CURVE LIST
  curve_list <-  list( r= r , obs =  as.vector(observed), sim_m = curve_sim_L)
  curve_set <- create_curve_set(curve_list)
  curve_set <- crop_curves(curve_set, r_min = 0, r_max = 30)
  res <- rank_envelope(curve_set,alpha=0.05)
  #Plot overall K
  #FIGURE 11 CAN BE CREATED LIKE THIS
  plot(res, use_ggplot2 = TRUE ,ylab = TeX('L(r)-r') )
}

#MILD CASE(CODE IS SIMILAR)
if (Group=='Mild'){
    Curve_sim = zeros(300,2500)
  for (nsim in 1:2500){
    tic()
    print(paste(nsim))
    
    Noclike_3d_mild = apply(seq_along(realbpMILD), function(x) Noc_model(realbpdfmild[[x]],params_m,direction_list_mild[[x]],SIDM[x]))
    endsim_2d_noc_mild =lapply(Noclike_3d_mild , function (x) subset(x, Endpoint=='TRUE'))
    
    endsimNOC_mild =lapply(seq_along(endsim_2d_noc_mild) ,  function(x) pp3(endsim_2d_noc_mild[[x]]$X,endsim_2d_noc_mild[[x]]$Y,
                                                                            endsim_2d_noc_mild[[x]]$Z,window = realep3Dmild[[x]]$domain))
    toc()
    Curve_sim[,nsim] = K_overall_function(endsimNOC_mild,samplesid_per_subj,3,TRUE)$K
  }
  
  curve_sim_L=zeros(300,2500)
  r = K_overall_function(realep3Dmild,samplesid_per_subj,3,TRUE)$r
  for (nsim in 1:2500){
    curve_sim_L[,nsim] = (Curve_sim[,nsim]/(4*pi/3))^0.33 - r
  }
  
  observed = (K_overall_function(realep3Dmild, samplesid_per_subj,3,TRUE)$K/(4*pi/3))^0.33 - r
  #r = K_overall_function(realep3D,samplesid_per_subj,3,TRUE)$r
  curve_list <-  list( r= r , obs =  as.vector(observed), sim_m = curve_sim_L)
  
  curve_set <- create_curve_set(curve_list)
  curve_set <- crop_curves(curve_set, r_min = 0, r_max = 30)
  res <- rank_envelope(curve_set,alpha=0.05)
  plot(res, use_ggplot2 = TRUE ,ylab = TeX('L(r)-r') )
}


#FIGURES 11 AND 13 CAN BE CREATED WITH ABOVE CODE BY CHANGING THE SUMMARY FUNCTION OF INTEREST.