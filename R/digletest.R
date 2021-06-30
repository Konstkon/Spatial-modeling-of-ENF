#envelope file
source('./helpool.R')
source('./Data_load.R')

count_sub <- function(data, fun , SID , Group){
Nsamp = samplerpersubj( SID, unique(SID))                  # samples per subject
Change_Subj = changeofsubj( Nsamp ) 
count = unlist(lapply(data, function(x) length(x$data$y))) # Number of endpoints per sample
subject_count = totc( count , Change_Subj )
subject_count
}

count_sam <- function(data, fun , SID , Group){
  Nsamp = samplerpersubj( SID, unique(SID))                  # samples per subject
  Change_Subj = changeofsubj( Nsamp ) 
  count = unlist(lapply(data, function(x) length(x$data$y))) # Number of endpoints per sample
 count
}

K_subject_wise<- function(data, fun , SID , Group){
Nsamp = samplerpersubj( SID, unique(SID))                  # samples per subject
Change_Subj = changeofsubj( Nsamp )                        # INDEX WHERE WE HAVE CHANGE IN SUBJECT
count = unlist(lapply(data, function(x) length(x$data$y))) # Number of endpoints per sample
subject_count = totc( count , Change_Subj )    # Number of endpoints per subjec
weights = calculateweights( subject_count , count , Change_Subj  ) # WEIGHTS FOR THE SUBJECT SPECIFIC SUMMARY FUNCTION
Nsubjectsqr = subject_count^2 # number of endpoints squared per subject
Ng = sum ( Nsubjectsqr )      #sum of square number of endpoints per subject
Ksubjectspecific = Subjectspecificnormal( fun , Change_Subj , weights ,Group) #K_i for subject i
Ksubjectspecific
}



K_group <- function ( data, fun , SID , K_sub ){
  Nsamp = samplerpersubj( SID, unique(SID))                  # samples per subject
  Change_Subj = changeofsubj( Nsamp )                        # INDEX WHERE WE HAVE CHANGE IN SUBJECT
  count = unlist(lapply(data, function(x) length(x$data$y))) # Number of endpoints per sample
  subject_count = totc( count , Change_Subj )    # Number of endpoints per subjec
  weights = calculateweights( subject_count , count , Change_Subj  ) # WEIGHTS FOR THE SUBJECT SPECIFIC SUMMARY FUNCTION
  Nsubjectsqr = subject_count^2 # number of endpoints squared per subject
  Ng = sum ( Nsubjectsqr )      #sum of square number of endpoints per subject
  sum = zeros( length( K_sub[[ 1 ]]$pooltrans) , 1 )
  for (k in 1:length( K_sub ) ){
    sum =  sum + Nsubjectsqr[ k ] * K_sub[[ k ]]$pooltrans
  }
  Kg = sum / Ng
}

create_fv <- function (K_sub, K_group){
df = data.frame (r = K_sub[[1]]$r , trans = K_group)
K_g =fv(df,valu ="trans")
}

get_group_K <- function(data,fun,SID,Group){
  K_i = K_subject_wise (data,fun,SID,Group)
  K_n = K_group(data, fun , SID , K_i)
  K = create_fv (K_i,K_n)
  K
  }

get_overall_K <- function (K_normal, K_mild , K_sub){
  N_mild = 1439
  N_norm = 8326
  N_tot = N_mild + N_norm
  K = N_norm/N_tot * K_normal$trans +N_mild/N_tot * K_mild$trans
  K_over = create_fv(K_sub,K)
  K_over
}
############
data= realep3D
fun = Kreal3D_end
SID = SID
Group = 'Normal'
K_sub_normal = K_subject_wise (data,fun,SID,'Normal')
K_normal = get_group_K(data,fun,SID,Group)
count_normal = count_sub(data,fun,SID,Group)
sqrt_count_norm = sqrt(count_normal)#(K_normal$trans) 

R_sub_norm=c()
for (k in 1:31){
  R_sub_norm = cbind(R_sub_norm, sqrt_count_norm[k] * (K_sub_normal[[k]]$pooltrans-K_normal$trans))
}
R_sub_norm


##
data= realep3Dmild
fun = Kreal3D_end_mild
SIDM = SIDM
Group = 'Mild'
K_sub_mild = K_subject_wise (data,fun,SIDM,'Mild')
K_mild = get_group_K(data,fun,SIDM,Group)
count_mild = count_sub(data,fun,SIDM,Group)
sqrt_count_mild = sqrt(count_mild)

R_sub_mild=c()
for (k in 1:8){
  R_sub_mild = cbind(R_sub_mild, sqrt_count_mild[k] * (K_sub_mild[[k]]$pooltrans-K_mild$trans))
}
R_sub_mild


##
K_overall = get_overall_K(K_normal , K_mild, K_sub_normal)
tot_count_norm = sum(count_normal)
tot_count_mild = sum (count_mild)

#APROXIMATE D
D = tot_count_norm*(K_normal$trans-K_overall$trans)^2+tot_count_mild*(K_mild$trans-K_overall$trans)^2
D = D/K_normal$r^2 
D = 0.334482 * D
D = sum(D[4:length(D)])

#Residual K function
count_all = c(count_mild,count_normal)
R_all = cbind(R_sub_norm,R_sub_mild)

#

get_K_star <- function(R_all, count_all,K_overall){
R_permuted = R_all[,sample(ncol(R_all))]
K_star_sub = c()
for (j in 1:39){
K_star_sub = cbind(K_star_sub, K_overall$trans + R_permuted[,j]/ sqrt(count_all[j]))
}
K_star_sub
}

#
res = c()
for (p in 1:9999){
K_star_subj = get_K_star(R_all,count_all,K_overall)
K_star_subj_mild = K_star_subj[,1:8]
K_star_subj_norm = K_star_subj[,9:39]
#

N_mild =sum ( count_mild^2 )
N_norm = sum ( count_normal^2 )
################

sum = zeros( length( K_star_subj_mild[,1]) , 1 )
for (k in 1:8 ){
  sum =  sum + count_mild[k]^2 * K_star_subj[,k]
}
Kstr_mild = sum / N_mild
sum = zeros( length( K_star_subj_mild[,1]) , 1 )
for (k in 9:39 ){
  sum =  sum + count_all[k]^2 * K_star_subj[,k]
}
Kstr_norm = sum / N_norm

N_mild = 1439
N_norm = 8326
N_tot = N_mild + N_norm
K = N_norm/N_tot * Kstr_norm +N_mild/N_tot * Kstr_mild
Kstr_overall = create_fv(K_sub_normal,K)
#Kstr_overall = get_overall_K(Kstr_norm , Kstr_mild, K_sub_normal)
D1 = tot_count_norm*(Kstr_norm - Kstr_overall$trans)^2+tot_count_mild*(Kstr_mild - Kstr_overall$trans)^2
D1 = D1/K_normal$r^2 
D1 = 0.334482 * D1
D1 = sum(D1[is.na(D1)==FALSE])

res = c (res, D1)
}
res1 = c(res,D)
res_sort = sort(res1)
which(res_sort == D)



