# function K_overall to compute an overall K function
# ppp_list: a list with ppp
# sample_list: a list with the sample ids of every subject
# dim: dimension
# weighted : if TRUE it returns the weighted mean function if not the median
# returns the overal K_functions


K_overall_function <- function(ppp_list,sample_list, dim ,weighted){
  list_size = length(ppp_list)
  subject_size = length(sample_list)
  
  #dim =2
  if (dim ==2){
    K = lapply(ppp_list, function(x) Kest(x,nrval =300,rmax=80))
    #weighted = true
    if(weighted==TRUE){
      numset = c(1:subject_size)
      l = c()
      w = c()
      w1 = c()
      K_subject_wise= c()
      #calculate weights
      for (num in numset){
        l = c()
        for (i in 1:length(sample_list[[num]])){    #########CHANGE ACCORDINGLY
          l=c(l, ppp_list[[i]]$n)
        }
        w1 = c(w1,sum(l))
        w = c(w,l /sum(l))
      }
      #K_subject_wise
      countw=0
      K_subject_wise= c()
      for (num in numset){
        K_sub_temp = 0
        for (i in  1:length(sample_list[[num]])){
          countw = countw+1
          K_sub_temp = K[[i]]$trans*w[countw] + K_sub_temp
        }
        K_subject_wise= cbind(K_subject_wise,K_sub_temp)
      }
      #calculate weights
      w2 = w1^2
      w2 = w2/sum(w2)
      K_overal = 0
      countw=0
      
      for (i in  1:length(numset)){
        countw = countw+1
        K_overal = K_subject_wise[,i]*w2[countw] + K_overal}
      radius = K[[1]]$r
      df = data.frame(K = K_overal, r=radius)
      df
      
    }
    else{
       K_r_sub = c()
      for (k in 1:length(K)){
        K_r_sub = cbind(K_r_sub, K[[k]]$trans)
      }
      K_med_real= c()
      for (j in 1:length(K_r_sub[,1])){
        K_med_real = c(K_med_real,median(K_r_sub[j,]))
      }
      radius = K[[1]]$r
      df = data.frame(K = K_med_real, r=radius)
      df
    }
    
  }
  else if(dim == 3){
    K = lapply(ppp_list, function(x) K3est(x,nrval =128,rmax=80))
    if(weighted==TRUE){
      numset = c(1:subject_size)
      l = c()
      w = c()
      w1 = c()
      K_subject_wise= c()
      #calculate weights
      for (num in numset){
        l = c()
        for (i in 1:length(sample_list[[num]])){
          l=c(l, length(ppp_list[[i]]$data$x))
        }
        w1 = c(w1,sum(l))
        w = c(w,l /sum(l))
      }
      #K_subject_wise
      countw=0
      K_subject_wise= c()
      for (num in numset){
        K_sub_temp = 0
        for (i in  1:length(sample_list[[num]])){
          countw = countw+1
          K_sub_temp = K[[i]]$trans*w[countw] + K_sub_temp
        }
        K_subject_wise= cbind(K_subject_wise,K_sub_temp)
      }
      #calculate weights
      w2 = w1^2
      w2 = w2/sum(w2)
      K_overal = 0
      countw=0
      
      for (i in  1:length(numset)){
        countw = countw+1
        K_overal = K_subject_wise[,i]*w2[countw] + K_overal}
      radius = K[[1]]$r
      df = data.frame(K = K_overal, r=radius)
      df
      
    }
    else{
      K_r_sub = c()
      for (k in 1:length(K)){
        K_r_sub = cbind(K_r_sub, K[[k]]$trans)
      }
      K_med_real= c()
      for (j in 1:length(K_r_sub[,1])){
        K_med_real = c(K_med_real,median(K_r_sub[j,]))
      }
      radius = K[[1]]$r
      df = data.frame(K = K_med_real, r=radius)
      df
    }
  }
  else {
    print('Please provide a valid dimension(2 or 3)')
    }
}
K_overall_function_group <- function(ppp_list,sample_list, statistic = "K",dim ,weighted){
  list_size = length(ppp_list)
  subject_size = length(sample_list)
  
  #dim =2
  if (dim ==2){
    if(statistic=="K"){
    K = lapply(ppp_list, function(x) Kest(x,nrval =300,rmax=80))
    }   
    if(statistic=="G"){
      K = lapply(ppp_list, function(x) Gest(x))
    }
    if(statistic=="F"){
      K = lapply(ppp_list, function(x) Fest(x))
    }
    if(statistic=="J"){
      K = lapply(ppp_list, function(x) Fest(x,))
    }
    #weighted = true
    if(weighted==TRUE){
      numset = c(1:subject_size)
      l = c()
      w = c()
      w1 = c()
      K_subject_wise= c()
      #calculate weights
      for (num in numset){
        l = c()
        for (i in sample_list[[num]]){    #########CHANGE ACCORDINGLY
          l=c(l, ppp_list[[i]]$n)
        }
        w1 = c(w1,sum(l))
        w = c(w,l /sum(l))
      }
      #K_subject_wise
      countw=0
      K_subject_wise= c()
      for (num in numset){
        K_sub_temp = 0
        for (i in  sample_list[[num]]){
          countw = countw+1
          if(statistic=="K"){
          K_sub_temp = K[[i]]$trans*w[countw] + K_sub_temp
          }
          if(statistic!="K"){
            K_sub_temp = K[[i]]$rs*w[countw] + K_sub_temp
            
          }
        }
        K_subject_wise= cbind(K_subject_wise,K_sub_temp)
      }
      #calculate weights
      w2 = w1^2
      w2 = w2/sum(w2)
      K_overal = 0
      countw=0
      
      for (i in  1:length(numset)){
        countw = countw+1
        K_overal = K_subject_wise[,i]*w2[countw] + K_overal}
      radius = K[[1]]$r
      df = data.frame(K = K_overal, r=radius)
      df
      
    }
    else{
      K_r_sub = c()
      for (k in 1:length(K)){
        K_r_sub = cbind(K_r_sub, K[[k]]$trans)
      }
      K_med_real= c()
      for (j in 1:length(K_r_sub[,1])){
        K_med_real = c(K_med_real,median(K_r_sub[j,]))
      }
      radius = K[[1]]$r
      df = data.frame(K = K_med_real, r=radius)
      df
    }
    
  }
  else if(dim == 3){
    K = lapply(ppp_list, function(x) K3est(x,nrval =128,rmax=80))
    if(weighted==TRUE){
      numset = c(1:subject_size)
      l = c()
      w = c()
      w1 = c()
      K_subject_wise= c()
      #calculate weights
      for (num in numset){
        l = c()
        for (i in sample_list[[num]]){
          l=c(l, length(ppp_list[[i]]$data$x))
        }
        w1 = c(w1,sum(l))
        w = c(w,l /sum(l))
      }
      #K_subject_wise
      countw=0
      K_subject_wise= c()
      for (num in numset){
        K_sub_temp = 0
        for (i in  sample_list[[num]]){
          countw = countw+1
          K_sub_temp = K[[i]]$trans*w[countw] + K_sub_temp
        }
        K_subject_wise= cbind(K_subject_wise,K_sub_temp)
      }
      #calculate weights
      w2 = w1^2
      w2 = w2/sum(w2)
      K_overal = 0
      countw=0
      
      for (i in  1:length(numset)){
        countw = countw+1
        K_overal = K_subject_wise[,i]*w2[countw] + K_overal}
      radius = K[[1]]$r
      df = data.frame(K = K_overal, r=radius)
      df
      
    }
    else{
      K_r_sub = c()
      for (k in 1:length(K)){
        K_r_sub = cbind(K_r_sub, K[[k]]$trans)
      }
      K_med_real= c()
      for (j in 1:length(K_r_sub[,1])){
        K_med_real = c(K_med_real,median(K_r_sub[j,]))
      }
      radius = K[[1]]$r
      df = data.frame(K = K_med_real, r=radius)
      df
    }
  }
  else {
    print('Please provide a valid dimension(2 or 3)')
  }
}


##EXAMPLES HOW TO RUN
#a1 = K_overall_function(realep, samplesid_per_subj,2,FALSE)
#a2 = K_overall_function(realep, samplesid_per_subj,2,TRUE)

#a3 = K_overall_function(realep3D, samplesid_per_subj,3,TRUE)
#a4 = K_overall_function(realep3D, samplesid_per_subj,3,FALSE)
#plot(a1$r, sqrt(a1$K/pi)-a1$r,'l')
#lines(a2$r, sqrt(a2$K/pi)-a2$r,'l')

#plot(a3$r, (a3$K/(4*pi/3))^0.33-a3$r ,'l')
#lines(a4$r, (a4$K/(4*pi/3))^0.33-a4$r ,'l')
