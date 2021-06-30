#k weighted 

K_weighted <- function(K_samples , weights, transl = TRUE ){
  K_weighted_overall= 0
 for (num in 1:length(K_samples)){
    if(transl==TRUE){
    K_weighted_overall = K_weighted_overall + K_samples[[num]]$trans*weights[num]
    }
    if(transl==FALSE){
      K_weighted_overall = K_weighted_overall + K_samples[[num]]*weights[num]
    }
    }
  K_weighted_overall
}



# Digle test using the Overal K function
Kreal3D_end = lapply(realep3D, function(x) K3est(x,nrvals =200,rmax=80))
R = Kreal3D_end[[1]]$r

n_ij = unlist(lapply(realep, function(x) x$n))
n_ij_sqr = n_ij^2
total_points = sum(n_ij)
total_points_sqr = sum(n_ij_sqr)


weights_sqr = n_ij_sqr/total_points_sqr
weights = n_ij/total_points

K_group_= K_weighted(Kreal3D_end, weights) # Overal K function from the replicated data

#MILD
Kreal3D_end_mild = lapply(realep3Dmild, function(x) K3est(x,nrvals =200,rmax=80))
R = Kreal3D_end[[1]]$r

n_ij_mild = unlist(lapply(realepmild, function(x) x$n))
n_ij_sqr_mild = n_ij_mild^2
total_points_sqr_mild = sum(n_ij_sqr_mild)
total_points_mild = sum(n_ij_mild)

weights_mild = n_ij_mild/total_points_mild

K_group_mild= K_weighted(Kreal3D_end_mild, weights_mild) # Overal K function from the replicated data



####

K_ij = lapply(Kreal3D_end , function(x) x$trans)           # K_ij for sample i of subkect j





############ Calculate R_ij

R_ij=c()
for (num in 1:32){
  for (samp in samplesid_per_subj[[num]]){
    R_ij = cbind(R_ij, sqrt(n_ij[samp]) *(K_ij[[samp]]   - K_group_))
  }
}
R_ij

Kstar_group = c()
for (p in 1:10000){
  
  ################## Permute r_ij
  Rstar_ij = R_ij[,sample(ncol(R_ij),replace=TRUE)]
  Kstar_ij = c()
  idx = sample.int(ncol(R_ij),112,replace=TRUE)
  for (num in idx){
    Kstar_ij =  cbind(Kstar_ij,K_group_ + Rstar_ij[,num]/ sqrt(n_ij[num]))
  }
 K = 0
 for (n in 1:112){
  K = K + weights[n]*Kstar_ij[,n]
    }

  Kstar_group = cbind(Kstar_group,K)
#lines(R,(K/(4*pi/3))^0.33 - R)
  
}


Kstar_lowq = c()
Kstar_hiq = c()
for (j in 1:length(R)){
  Kstar_lowq= c(Kstar_lowq , quantile(Kstar_group[j,],prob=0.0275))
  Kstar_hiq= c(Kstar_hiq , quantile(Kstar_group[j,],prob=0.975))
}

#mild
R_ij_mild=c()
for (num in 1:8){
  for (samp in samplesid_per_subj_mild[[num]]){
    R_ij_mild = cbind(R_ij_mild, sqrt(n_ij_mild[samp]) *(K_ij_mild[[samp]]   - K_group_mild))
  }
}
R_ij_mild

Kstar_group_mild = c()
for (p in 1:10000){
  
  ################## Permute r_ij
  Rstar_ij_mild = R_ij_mild[,sample(ncol(R_ij_mild),replace=TRUE)]
  Kstar_ij_mild = c()
  idx = sample.int(ncol(R_ij_mild),28,replace=TRUE)
  for (num in idx){
    Kstar_ij_mild =  cbind(Kstar_ij_mild,K_group_mild + Rstar_ij_mild[,num]/ sqrt(n_ij_mild[num]))
  }
  K = 0
  for (n in 1:28){
    K = K + weights_mild[n]*Kstar_ij_mild[,n]
  }
  
  Kstar_group_mild = cbind(Kstar_group_mild,K)
  #lines(R,(K/(4*pi/3))^0.33 - R)
  
}


Kstar_lowq_mild = c()
Kstar_hiq_mild = c()
for (j in 1:length(R)){
  Kstar_lowq_mild= c(Kstar_lowq_mild , quantile(Kstar_group_mild[j,],prob=0.0275))
  Kstar_hiq_mild= c(Kstar_hiq_mild , quantile(Kstar_group_mild[j,],prob=0.975))
}
  
  
  
id = which(K_group_>=0)
par(mfrow=c(1,1))

#Figure 10
plot(R, (K_group_/(4*pi/3))^0.33 - R,'l',xlim=c(0,60), 
     lwd=2,col=2,ylim=c(-5,35),ylab='L(r)-r' ,xlab='r')

lines(R, (K_group_mild[id]/(4*pi/3))^0.33 - R,'l',xlim=c(0,60), 
     lwd=2,col=1,ylim=c(-5,30),ylab='L(r)-r' ,xlab='r')





lines(R , (Kstar_hiq/(4*pi/3))^0.33 - R,lty=3,col=2 )
lines(R , (Kstar_lowq/(4*pi/3))^0.33 - R,lty=3,col=2)
lines(R , (Kstar_hiq_mild/(4*pi/3))^0.33 - R,lty=3,col=1 )
lines(R , (Kstar_lowq_mild/(4*pi/3))^0.33 - R,lty=3,col=1)
legend("topright",lty=c(1,1),col=c(1,2),c("Mild","Healthy"),lwd=c(2,2))
