# Functions used for pooling the summary functions

#samples per subject
samplerpersubj <- function ( l , unl ){
  z = as.vector( ones( length( unl ) , 1 ) )
  i = 1
  for ( k in 2:length( l ) ){
    if (l[ k ] == l[ k - 1 ] ){
      z[ i ] = z[ i ] + 1
    }
    else{ i = i + 1 }
  }
  z
}

# index where we have a change of subject
changeofsubj <- function ( z ) {
  sums = zeros ( length( z ) , 1 )
  sums[ 1 ]  = z[ 1 ]
  for ( k in 2:length( z ) ){
    sums[ k ] = sums[ k - 1 ] + z[ k ]
  }
  sums = as.vector( sums )
}


#calculate weights
totc = function ( count , sums  ){
  totalcount = sum(count[1:sums[1]])
  for (k in 2:length(sums)){
    tc <-  sum(count[sums[k-1]+1:(sums[k]-sums[k-1])] )
    totalcount <- c(totalcount, tc)
  }
  totalcount
}

#weights
calculateweights <- function ( totalcount , count, sums ){
  weights = totalcount[ 1 ] * ones( sums[ 1 ] , 1 )
  for ( k in 2:length( sums ) ){
    w = totalcount[ k ] * ones( sums[ k ] - sums[ k - 1 ] , 1 )
    weights = c( weights , w )
  }
  weights = count/weights
  weights
}

#Compute subject specific weighted K function
Subjectspecificnormal <- function( K , sums ,weights, Group  ){
  
    Ksubj = list( anylist( K[[ 1 ]] , K[[ 2 ]] , K[[ 3 ]] , K[[ 4 ]]) ) # first subject
    
   # if (Group=='Mild'){
   #   Ksubj = list(anylist(K[[1]], K[[ 2 ]] , K[[ 3 ]]))
   # }
  for (k in 2: length( sums ) ){
    Ktemp = list( K[[ sums[ k - 1 ] + 1 ]])
    if (sums[ k ] - sums[ k - 1 ] > 1 ){
      for (j in 2:(sums[ k ] - sums[ k - 1 ])){
        Ktemp = append( Ktemp , K[ sums[ k - 1 ] + j ] )
      }}
      Ktemp = list( as.anylist( Ktemp ) )
      Ksubj = append( Ksubj , Ktemp )

  }
  Ksubjectspecific = list(pool ( Ksubj[[1]] , weights = weights[1:sums[1]]))
  for (k in 2:length(Ksubj)){
    if (sums[ k ] - sums[ k - 1 ] > 1 ){
    Ksubjectspecific = append(Ksubjectspecific, list(pool ( Ksubj[[ k ]] , weights = weights[(sums[ k-1 ] + 1 ):sums[k]])))
    }}
  Ksubjectspecific
}

###
createmarkpppsim <- function( sim ){
  bp = ppp ( sim$X[sim$Basepoint==TRUE] ,sim$Y[sim$Basepoint==TRUE] ,window = owin(xrange=c(0,432), yrange=c(-330,0)))
  bp = setmarks(bp,as.factor('BP'))
  ep = ppp ( sim$X[sim$Endpoint==TRUE] ,sim$Y[sim$Endpoint==TRUE] ,window = owin(xrange=c(0,432), yrange=c(-330,0)))
  ep = setmarks(ep,as.factor('EP'))
  pp =  superimpose(bp,ep)
  pp
}
# Function that Returns the poin process of the end points given a simulation
geteppp<- function(sim,bp){
  win =bp$window
  ep = ppp ( sim$X[sim$Endpoint==TRUE] ,sim$Y[sim$Endpoint==TRUE],
            window = win,2 )
   ep
}
# Function that Returns the poin process of the end points given a simulation 3D
geteppp3D<- function(sim,bp){
  win = bp$domain
  ep = pp3 ( sim$X[sim$Endpoint==TRUE] ,sim$Y[sim$Endpoint==TRUE], sim$Z[sim$Endpoint==TRUE],
             win )
  ep
}
# Function that Returns the poin process of the first branching points given a simulation
getfbppp <- function (sim){
  fb =ppp ( sim$X[sim$Branchpoint==TRUE] ,sim$Y[sim$Branchpoint==TRUE] ,window = owin(xrange=c(0,432), yrange=c(-330,0.5)))
    fb
}
# Function that Returns the poin process of the base points given a simulation
getbppp <- function (sim){
  bp = ppp( sim$X[sim$Basepoint==TRUE] ,sim$Y[sim$Basepoint==TRUE] ,
window = owin(xrange=c(min(sim$X[sim$Basepoint==TRUE])-10,min(sim$X[sim$Basepoint==TRUE])+422),
              yrange=c(min(sim$Y[sim$Basepoint==TRUE])-10,min(sim$Y[sim$Basepoint==TRUE]) +320)))
  bp
}

# Function that creates a marked point process for base and end points for the data
# Translates the process and data to the correct observation window

createmarkpppdat <- function( bp , ep ){

  #bp = setmarks( bp, as.factor( 'BP' ) )
  #ep = setmarks( ep, as.factor( 'EP' ) )
  pp =  superimpose( bp , ep )
  pp1 = shift(pp,c(-window(pp)[[1]][[2]][[1]] ,-window(pp)[[1]][[3]][[2]]))
  pp1

}

Confidenseplot <- function ( observed , simulation ){
  lbobs = lohboot( observed , Ldot , nsim = 100  , confidence = 0.98)
  lbsim = lohboot( simulation , Ldot , nsim = 100 , confidence = 0.98)
  plot( lbsim$r , lbsim$lo - lbsim$r , 'l' , col = 2  , xlim = c(0,80) , ylim = c(-5,30) , xlab= 'r', ylab = 'Lbp,ep(r)-r')
  points( lbsim$r , lbsim$hi - lbsim$r , 'l' , col = 2 )
  points( lbsim$r , lbsim$hi - lbsim$r , 'l' , col = 2 )
  points( lbobs$r , lbobs$hi - lbobs$r , 'l' , col = 1 )
  points( lbobs$r , lbobs$lo - lbobs$r , 'l' , col = 1 )
  points( lbobs$r , lbobs$iso - lbobs$r , 'l' , col = 3 )
  legend('bottomright',legend = c('conf band for observed','conf band for simulation','Isotropic correction for observed'),col=c(1,2,3),lty=c(1,1,1),cex=0.7)
  ###
}

getnum_ep<- function( pp ){
  N=pp$n
  N
}

Distoff <- function(df){
  res=c()
  for (i in 1:max(df$Tree)){
    if (length(df$X[df$Tree==i])>0){
    res=c(res,length(df$X[df$Tree==i]))
  }}
  res
}


Simulate_pattern<- function( Group){
  if (Group=='Normal'){
  simulation = lapply(footlist, function(x) SimulateENFdata(get.basepoints(x),Group))
  for (x in 1 : length( simulation ) ){
    simulation[[ x ]]$SID = footlist[[x]]$Covariates$SubjectID
  }}
  if (Group== 'Mild'){
    simulation = lapply(footlistmild, function(x) SimulateENFdata(get.basepoints(x),Group))
    for (x in 1 : length( simulation ) ){
      simulation[[ x ]]$SID = footlistmild[[x]]$Covariates$SubjectID
    }}
  
  simulation


}


Simulate_pattern_SUBJECT<- function( Group, params){
  if (Group=='Normal'){
    simulation = lapply(footlist, function(x) SimulateENFdata_subject(get.basepoints(x),Group, params , as.numeric(Get_subject_from_str(x$name))))
    for (x in 1 : length( simulation ) ){
      simulation[[ x ]]$SID = footlist[[x]]$Covariates$SubjectID
    }}
  if (Group== 'Mild'){
    simulation = lapply(footlistmild, function(x) SimulateENFdata_subject(get.basepoints(x),Group,  params , as.numeric(Get_subject_from_str(x$name))))
    for (x in 1 : length( simulation ) ){
      simulation[[ x ]]$SID = footlistmild[[x]]$Covariates$SubjectID
    }}
  
  simulation
  
  
}

GetsubjectID <- function (list ){

  list$Covariates$SubjectID
}


Simulate_pattern1<- function( Group){
  if (Group=='Mild'){
  simulation = lapply(footlistmild, function(x) SimulateENFdata(get.basepoints(x),Group))
  for (x in 1 : length( simulation ) ){
    simulation[[ x ]]$SID = footlistmild[[x]]$Covariates$SubjectID
  }}
  if (Group=='Normal'){
    simulation = lapply(footlist, function(x) SimulateENFdata(get.basepoints(x),Group))
    for (x in 1 : length( simulation ) ){
      simulation[[ x ]]$SID = footlist[[x]]$Covariates$SubjectID
    }
  }
  simulation
}


plot_summary<- function(data, fun , SID , Group, cond){
  Nsamp = samplerpersubj( SID, unique(SID))
  Change_Subj = changeofsubj( Nsamp ) # INDEX WHERE WE HAVE CHANGE IN SUBJECT
  count = unlist(lapply(data, function(x) x$n)) # Number of endpoints per sample
  totalcount = totc( count , Change_Subj )    # Number of endpoints per subjec
  weights = calculateweights( totalcount , count , Change_Subj  ) # WEIGHTS FOR THE SUBJECT SPECIFIC SUMMARY FUNCTION
  Nsubjectsqr = totalcount^2 # number of endpoints squared per subject
  Ng = sum ( Nsubjectsqr )
  Ksubjectspecific = Subjectspecificnormal( fun , Change_Subj , weights ,Group)
  sum1 = zeros( length( Ksubjectspecific[[ 1 ]]$pooliso) , 1 )
  for (k in 1:length( Ksubjectspecific ) ){
    sum1=  sum1 + Nsubjectsqr[ k ] * Ksubjectspecific[[ k ]]$pooliso
  }
  Kg = sum1 / Ng
  Lg = sqrt( Kg / pi ) - Ksubjectspecific[[1]]$r #overal centred group L function
  if (cond==TRUE){
    plot(Ksubjectspecific[[1]]$r,Lg,'l',xlab='r',ylab=' Lbp,ep(r)-r ',col=3, ylim=c(0,25),cex=1.2)
  }
  else {
    points(Ksubjectspecific[[1]]$r,Lg,'l',xlab='r',ylab=' Lbp,ep(r)-r ',col=1, ylim=c(0,25),cex=1.2)
  }
  Lg
}



plot_summary_3D<- function(data, fun , SID , Group, cond){
  Nsamp = samplerpersubj( SID, unique(SID))
  Change_Subj = changeofsubj( Nsamp ) # INDEX WHERE WE HAVE CHANGE IN SUBJECT
  count = unlist(lapply(data, function(x) length(x$data$y))) # Number of endpoints per sample
  totalcount = totc( count , Change_Subj )    # Number of endpoints per subjec
  weights = calculateweights( totalcount , count , Change_Subj  ) # WEIGHTS FOR THE SUBJECT SPECIFIC SUMMARY FUNCTION
  Nsubjectsqr = totalcount^2 # number of endpoints squared per subject
  Ng = sum ( Nsubjectsqr )
  Ksubjectspecific = Subjectspecificnormal( fun , Change_Subj , weights ,Group)
  sum1 = zeros( length( Ksubjectspecific[[ 1 ]]$pooltrans) , 1 )
  for (k in 1:length( Ksubjectspecific ) ){
    sum1=  sum1 + Nsubjectsqr[ k ] * Ksubjectspecific[[ k ]]$pooltrans
  }
  Kg = sum1 / Ng
  Lg = ( Kg / (4*pi/3 ))^(1/3) - Ksubjectspecific[[1]]$r #overal centred group L function
  if (cond==TRUE){
    plot(Ksubjectspecific[[1]]$r,Lg,'l',xlab='r',ylab=' L(r)-r ',xlim=c(0,100),ylim=c(0,30),col=2,cex=1.2)
  }
  else {
    points(Ksubjectspecific[[1]]$r,Lg,'l',xlab='r',ylab=' L(r)-r ',col=1,cex=1.2)
    
  }
  Lg
}

