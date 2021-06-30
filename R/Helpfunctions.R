#load libraries
require( ENFdata )
require(circular)
require(MASS)
require(pracma)


source( './tree.structs.R')  #Code to get information of the tree structure from a pattern
source( './get.nonorf.ep.R') #Code to get the non-orphan points from the pattern
source( './point.indices.R') #Gets point indices


#GET NUMBER OF ENPOINTS
getnum_ep<- function( pp ){
  N=pp$n
  N
}
#GET SUBJECT ID
GetsubjectID <- function (list ){
  
  list$Covariates$SubjectID
}

#CHANGE NAME(MAINLY USED TO LOAD DATA)
namechange <- function(x){
  y <- get( x )
  y$name <- x
  y
  
}
#RETURNS A MARKED PP3 
# Data : a ppp or pp3
# Dimension is 2 or 3
# df contating info about bp
# df_ord contains info about order


create_mark_bp <- function(data , df,df_ord, dimension){
  if (dimension == 2){
    res= lapply(seq_along(data),function(x) setmarks(data[[x]],
                                                     data.frame(Tree=df[[x]]$Tree,type = ones(data[[x]]$n,1),Volume_microm=df[[x]]$Volume_microm,
                                                                AvgDiameter_microm=df[[x]]$AvgDiameter_microm, SurfaceArea_microm2=df[[x]]$SurfaceArea_microm2,
                                                                BaseDiameter_microm=df[[x]]$BaseDiameter_microm, Order=df[[x]]$Order,Max_ord=df_ord[[x]]$Max_ord)))
    
  }
  if (dimension==3){
    res= lapply(seq_along(data),function(x) setmarks(data[[x]],
                                                     data.frame(Tree=df[[x]]$Tree,type = ones(length(data[[x]]$data$z),1),Volume_microm=df[[x]]$Volume_microm,
                                                                AvgDiameter_microm=df[[x]]$AvgDiameter_microm, SurfaceArea_microm2=df[[x]]$SurfaceArea_microm2,
                                                                BaseDiameter_microm=df[[x]]$BaseDiameter_microm, Order=df[[x]]$Order,Max_ord=df_ord[[x]]$Max_ord)))
  }
  res
}


# Return a pp3 with marks Tree and type
# df= dataframe with X,Y,Z coords and marks
# Base point type == 1
# End point type == 0
create_pp3_from_df <- function( df, type){
  
  res = lapply( df, function(x) pp3(x$X,x$Y,x$Z,box3(c(0,430),c(-330,0),c( min(x$Z)-5, max(x$Z)+35))))
  res = lapply(seq_along(df),function(x) setmarks(res[[x]],
                                                  data.frame(Tree=df[[x]]$Tree,type = type*ones(length(res[[x]]$data$z),1))))
  res
}

#GET BASEPOINTS MARKS
get.basepoints.marks <- function(x, names=c("Tree", "X", "Y", "Z","Order","BaseDiameter_microm","AvgDiameter_microm","SurfaceArea_microm2","Volume_microm")) {
  # x = an ENFdata object	
  if(class(x)!="ENFdata") stop("The function get.basepoints is only for ENFdata objects.\n")
  #cat("Basepoints chosen by: TerminalType=c(\"Branch\",\"Normal\"), Order=1, ENFtype=\"Axon\"\n")
  x$df[(x$df$TerminalType=="Branch" | x$df$TerminalType=="Normal") & x$df$Order==1 & x$df$ENFtype=="Axon", names]
}
#GET NON ORPHAN EP MARKS
get.nonorf.ep.df.marks <- function (x, X = "X", Y = "Y", unique = F) {
  
  
  Endpoints <- get.endpoints(x, names = c(X, Y,'Z','Tree','BranchOrder','Distance_microm'))
  Basepoints <- get.basepoints( x )
  Endpoints <- Endpoints[ which(Endpoints$Tree %in% Basepoints$Tree), ]
  
  
  Endpoints
  
}
#GET TREE ORDER
get_tree_order<- function(x, names=c("Tree","Max_order") ){
  treeID= unique(x$Tree)
  df= data.frame(Tree=c(), Max_ord=c())
  for (j in 1:length(treeID)){
    ind = which(x$Tree %in% treeID[j])
    Num_end = length(ind)
    Max_ord = max(x$BranchOrder[ind])
    df1 = data.frame(Tree = treeID[j], Max_ord= Max_ord)
    df=rbind(df,df1)
  }
  df
}
#GET ENDPOINT ORDER
get_end_order <-function(x, names=c("Tree","Order") ){
  treeID= unique(x$Tree)
  df= data.frame(Tree=c(), Order=c())
  for (j in 1:length(treeID)){
    ind = which(x$Tree %in% treeID[j])
    Num_end = length(ind)
    Order = x$BranchOrder[ind]
    df1 = data.frame(Tree = treeID[j]*ones(Num_end,1), Order= Order)
    df=rbind(df,df1)
  }
  df
}


#superimpose two pp3 patterns with marks

superimpose_pp3 <- function( pattern1, pattern2){
  marks1 = data.frame(Tree=marks(pattern1)$Tree, type=marks(pattern1)$type)
  marks2= data.frame(Tree=marks(pattern2)$Tree, type=marks(pattern2)$type)
  marks = rbind(marks1,marks2)
  
  w1 = pattern1$domain
  w2 = pattern2$domain
  zlow = min(c(w1$zrange[1],w2$zrange[1]))
  zhigh = max(c(w1$zrange[2],w2$zrange[2]))
  X = c( pattern1$data$x , pattern2$data$x )
  Y = c( pattern1$data$y , pattern2$data$y )
  Z = c( pattern1$data$z , pattern2$data$z )
  res = pp3(X,Y,Z,marks=marks, box3(w1$xrange,w2$yrange,c(zlow,zhigh)))
  res
}
#gETS DAUGHTERS
Getdaughep <-function(Segment,group,treestats){
  res=c()
  SID =c()
  if ( length( Segment ) > 0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length( Segment )){
        res<-c(res , length( Segment[[ i ]]$ep [ Segment[[i]]$ep == TRUE ]))
      SID= c(SID, Segment[[i]]$SubjectID[ 1 ])
        }

    }}
    df=data.frame(offspring=res,SID=SID)
df
  }

#function for number of daughters at 1st braching point

Getdaugh <- function(Segment,group,treestats){
  res=c()
  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        res<-c(res , Segment[[i]]$Daughters [ Segment[[i]]$Index == 2 ])
      }
    }
  }
  res
}

#returns distribution of number of branches p=(p1,p2,p3,p4) etc
getdaugt_distr<- function(dau){
  p <- c()
  for (k in 1:max(dau)){
    p <- c(p , length(dau[dau==k])/length(dau))
  }
  p
}

#Base points vs first branch
# returns xy angle, zangle and length between base points and 1st branching point
bp2fb<- function(Segment,group,treestats){
  res = data.frame(xyangle = c( ) ,
                   zangle = c( ),
                   length = c( ) )
  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        X2   <- Segment[[i]]$X [ Segment[[i]]$Index == 2 ]
        Y2   <- Segment[[i]]$Y [ Segment[[i]]$Index == 2 ]
        Z2   <- Segment[[i]]$Z [ Segment[[i]]$Index == 2 ]
        X1   <- Segment[[i]]$X [ Segment[[i]]$Index == 1 ]
        Y1   <- Segment[[i]]$Y [ Segment[[i]]$Index == 1 ]
        Z1   <- Segment[[i]]$Z [ Segment[[i]]$Index == 1 ]
        X <- X2-X1
        Y <- Y2-Y1
        Z <- Z2-Z1
        sph <- cart2sph(t(matrix(c(X,Y,Z))))
        d <- list(xyangle =sph[1], zangle =sph[2],length =sph[3])
        res <- rbind(res,d)
      }
    }
  }
  res
}

#Base points vs first branch
# returns xy angle, SID,zangle and length between base points and 1st branching point
bp2fb1<- function(Segment,group,treestats){
  res = data.frame(xyangle = c( ) ,
                   zangle = c( ),
                   length = c( ),
                   SID =  c())
                    
  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        X2   <- Segment[[i]]$X [ Segment[[i]]$Index == 2 ]
        Y2   <- Segment[[i]]$Y [ Segment[[i]]$Index == 2 ]
        Z2   <- Segment[[i]]$Z [ Segment[[i]]$Index == 2 ]
        X1   <- Segment[[i]]$X [ Segment[[i]]$Index == 1 ]
        Y1   <- Segment[[i]]$Y [ Segment[[i]]$Index == 1 ]
        Z1   <- Segment[[i]]$Z [ Segment[[i]]$Index == 1 ]
        X <- X2-X1
        Y <- Y2-Y1
        Z <- Z2-Z1
        sph <- cart2sph(t(matrix(c(X,Y,Z))))
        d <- list(xyangle =sph[1], zangle =sph[2],length =sph[3],SID=Segment[[1]]$SubjectID[1])
        res <- rbind(res,d)
      }
    }
  }
  res
}
#fUNCTION THAT RETURNS THE FIRST BRANCHING POINTS
Getbranchpoints <- function( Segment , group , treestats ){
  res = data.frame(X = c( ) ,
                   Y = c( ),
                   Z = c( ),
                   Tree = c() )
  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        X   <- Segment[[i]]$X [ Segment[[i]]$Index == 2 ]
        Y   <- Segment[[i]]$Y [ Segment[[i]]$Index == 2 ]
        Z   <- Segment[[i]]$Z [ Segment[[i]]$Index == 2 ]
        Tree   <- Segment[[i]]$Tree [ Segment[[i]]$Index == 2 ]
        SID <-Segment[[i]]$SubjectID[ Segment[[i]]$Index == 1 ]
        d <- list(X=X,Y=Y,Z=Z,Tree=Tree, SID=SID)
        res <- rbind(res,d)
      }
    }
  }
  res
}



#First branch vs end points
# returns xy angle, zangle and length between end points and 1st branching point
fb2end1 <- function(Segment,group,treestats){
  res = data.frame(xyangle = c( ) ,
                   zangle = c( ),
                   length = c( ) ,SID =  c())

  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        X1   <- Segment[[i]]$X [ Segment[[i]]$Index == 2 ]
        Y1   <- Segment[[i]]$Y [ Segment[[i]]$Index == 2 ]
        Z1   <- Segment[[i]]$Z [ Segment[[i]]$Index == 2 ]
        X2   <- Segment[[i]]$X[Segment[[i]]$ep==TRUE ]
        Y2   <- Segment[[i]]$Y[Segment[[i]]$ep==TRUE ]
        Z2   <- Segment[[i]]$Z[Segment[[i]]$ep==TRUE ]
        for (k in 1:length(Z2)){
          X= X2[k]-X1
          Y= Y2[k]-Y1
          Z= Z2[k]-Z1

            sph <- cart2sph(t(matrix(c(X,Y,Z))))
            d <- list(xyangle =sph[1], zangle =sph[2],length =sph[3],SID=Segment[[1]]$SubjectID[1])
            res <- rbind(res,d)}
      }


    }
  }
  res
}
#First branch vs end points
# returns xy angle, zangle and length between end points and 1st branching point
bp2end1 <- function(Segment,group,treestats){
  res = data.frame(xyangle = c( ) ,
                   zangle = c( ),
                   length = c( ) ,SID =  c())
  
  if (length(Segment)>0 ){
    if (treestats$diseasestatus[treestats$SID==Segment[[1]]$SubjectID[1]]==group){
      for (i in 1:length(Segment)){
        X1   <- Segment[[i]]$X [ Segment[[i]]$Index == 1 ]
        Y1   <- Segment[[i]]$Y [ Segment[[i]]$Index == 1 ]
        Z1   <- Segment[[i]]$Z [ Segment[[i]]$Index == 1 ]
        X2   <- Segment[[i]]$X[Segment[[i]]$ep==TRUE ]
        Y2   <- Segment[[i]]$Y[Segment[[i]]$ep==TRUE ]
        Z2   <- Segment[[i]]$Z[Segment[[i]]$ep==TRUE ]
        for (k in 1:length(Z2)){
          X= X2[k]-X1
          Y= Y2[k]-Y1
          Z= Z2[k]-Z1
          
          sph <- cart2sph(t(matrix(c(X,Y,Z))))
          d <- list(xyangle =sph[1], zangle =sph[2],length =sph[3],SID=Segment[[1]]$SubjectID[1])
          res <- rbind(res,d)}
      }
      
      
    }
  }
  res
}

#get functions
#PLANAR ANGLE
getxyangle <-function(x){
  res=c()
  for (i in x$xyangle){
    if (i<0) i =i+2*pi
    res<-c(res,i)
  }


  res
}
#LENGTH
getlen <-function(x){
  res=c()
  for (i in x$length){
    res<- c(res,i)}

  res
}
#ZANGLE
getzangle<-function(x){
  res=c()
  for (i in x$zangle){
    res<- c(res,i)
  }
  res
}


#recieves latitude and returns colatidute
lat2colat <- function(x){
  if ((x<=pi/2) & (x>=-pi/2))
  {colat<- pi/2-x
  colat
  }
  else print('Error: The values for the latitude x is not in the interval [-pi/2,pi/2]')
}


#SPHERICAL COORDINATES TO POINTS
spherepoints <- function(theta,phi){
  res = data.frame(X = c( ) ,
                   Y = c( ),
                   Z = c( ) )
  line <- list(X =sin(theta)*cos(phi), Y = sin(theta)*sin(phi),Z = cos(theta))
  res <- rbind(res,line)
  res
}


#SIMULATE FROM VONMISES
simulatevonmises <- function(N , k , mu ){
  res <-data.frame( thetap=c(),
                    phip=c()
                    )
  labda <- exp(-2*k)
  a <- mu[1]
  b <- mu[2]
  theta <- c()
  phi   <- c()
  thetaprime <- c()
  phiprime   <- c()
  for (i in 1:N){
    R1 = runif( 1 , 0 , 1 )
    R2 = runif( 1 , 0 , 1 )
    theta <-c( theta , 2*asin( sqrt( -log( ( R1*( 1-labda )+labda))/( 2*k )  ) ))
    phi <- c( phi , 2*pi*R2 )
  }
  dat <- c(cos (a)*cos(b) , cos(a)*sin(b), -sin(a) , -sin(b) , cos(b) , 0, sin(a)*cos(b) , sin(a)*sin(b), cos(a) )
  A <- matrix(dat ,nrow=3 , ncol=3 )
  A <- t(A)
  for (i in 1:length(theta)){
  coords = t(as.matrix(spherepoints(theta[i],phi[i]))) #give x,y,z coordinates on sphere
  B <- A %*% coords
  th <- acos( B[3] )
  ph <- acos( B[1]/sin(th))
  thetaprime <- c( thetaprime , th )
  phiprime <- c( phiprime , ph )
  }
  line <- list(thetap =thetaprime , phip = phiprime)
  res <- rbind(res,line)
  res
}

#Create a mark point process of endpoints and base points (for pp from simulation)

createmarkpppsim <- function( sim ){
  bp = ppp ( Sim$X[Sim$Basepoint==TRUE] ,Sim$Y[Sim$Basepoint==TRUE] ,window = owin(xrange=c(0,432), yrange=c(-330,0.5)))
  bp = setmarks(bp,as.factor('BP'))
  ep = ppp ( Sim$X[Sim$Endpoint==TRUE] ,Sim$Y[Sim$Endpoint==TRUE] ,window = owin(xrange=c(0,432), yrange=c(-330,0.5)))
  ep = setmarks(ep,as.factor('EP'))
  pp =  superimpose(bp,ep)
  pp
}

#Create a mark point process of endpoints and base points (for pp from data)

createmarkpppdat <- function( bp , ep ){
  bp = setmarks( bp, as.factor( 'BP' ) )
  ep = setmarks( ep, as.factor( 'EP' ) )
  pp =  superimpose( bp , ep )
  pp
}
#COUNTS ENDPOINTS
Countendpoints <- function (df){
  length(df$Endpoint[df$Endpoint==TRUE])
}


#GETS SUBJECT NAME
Get_subject_from_str <- function(name){
  if(substr(name ,9, 9)=='.'){
    ID= substr(name , 6 , 8 )
  }
  else if (substr(name,6,7)=='99'){ID='99'}else{ ID = substr(name , 6 , 9 )}
}
