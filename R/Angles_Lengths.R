## Load data and get length/angles
# RETURN ANGLES AND LENGTHS BETWEEN BASE AND FIRST BRANCHING POINTS GIVEN A DISEASEGROUP

Get_Firstbranch <- function(Group){

#get angles/lengths between base and first branch points

fb <- lapply(foot.trees, function(x) bp2fb1(x,Group,tree.stats))
xyangfb <- unlist(lapply(fb, function(x) getxyangle(x)))
zangfb <- unlist(lapply(fb, function(x) getzangle(x)))
lenfb <- unlist(lapply(fb, function(x) getlen(x)))
SID<- unlist(lapply(seq_along(fb), function(x) fb[[x]]$SID))
xyangfb <-xyangfb[zangfb>0] #remove possible wrong measurements
lenfb   <- lenfb[zangfb>0]
SID<-SID[zangfb>0]
zangfb  <- zangfb[zangfb>0]

colatfb <- lat2colat(zangfb)

First_branch = data.frame( Length = lenfb , XYangle = xyangfb ,  Zangle = zangfb, Colatidute = colatfb, SID=SID)
First_branch
}


## Load data and get length/angles
# RETURN ANGLES AND LENGTHS BETWEEN END AND FIRST BRANCHING POINTS GIVEN A DISEASEGROUP

Get_endbranch<- function(Group){


# get angles/lengths betwwen first branch and end points
fbend <- lapply(foot.trees, function(x) fb2end1(x,Group,tree.stats))
xyangep <- unlist(lapply(fbend, function(x) getxyangle(x)))
zangep <- unlist(lapply(fbend, function(x) getzangle(x)))
lenep <-  unlist(lapply(fbend, function(x) getlen(x)))
SID<- unlist(lapply(seq_along(fbend), function(x) fbend[[x]]$SID))
xyangep <-xyangep[lenep!=0] #remove info from 1st branching points that are also end points
zangep  <- zangep[lenep!=0]
SID<-SID[lenep!=0]
lenep   <- lenep[lenep!=0]
colatep <- lat2colat(zangep)
Endpoint_Branch = data.frame(Length = lenep, XYangle = xyangep , Zangle = zangep , Colatidute = colatep,SID=SID )
Endpoint_Branch
}



Get_endbranch_bp<- function(Group){
  
  
  # get angles/lengths betwwen first branch and end points
  fbend <- lapply(foot.trees, function(x) bp2end1(x,Group,tree.stats))
  xyangep <- unlist(lapply(fbend, function(x) getxyangle(x)))
  zangep <- unlist(lapply(fbend, function(x) getzangle(x)))
  lenep <-  unlist(lapply(fbend, function(x) getlen(x)))
  SID<- unlist(lapply(seq_along(fbend), function(x) fbend[[x]]$SID))
  xyangep <-xyangep[lenep!=0] #remove info from 1st branching points that are also end points
  zangep  <- zangep[lenep!=0]
  SID<-SID[lenep!=0]
  lenep   <- lenep[lenep!=0]
  colatep <- lat2colat(zangep)
  Endpoint_Branch = data.frame(Length = lenep, XYangle = xyangep , Zangle = zangep , Colatidute = colatep,SID=SID )
  Endpoint_Branch
}

fixZ <- function(Segment){
  Z=c()
  K=list()
  if (length(Segment)>0 ){       #check if there are trees
    
    for (i in 1:length(Segment)){ #loop through trees
      Z2   <- Segment[[i]]$Z [ Segment[[i]]$Index == 2 ]
    random = runif(1,0,1)          #uniform number in 0,1
    Z2<- Z2 + random
    Z=c(Segment[[i]]$Z [ Segment[[i]]$Index == 1 ],Z2)
    
    if (length(Segment[[i]]$X)>2){
  for (j in 3:length(Segment[[i]]$X)){
  random2 = runif(1,0,1-random)  #uniform number in 0,1-random
  Z <- c(Z,Segment[[i]]$Z [ Segment[[i]]$Index == j ]+random+random2)
  }
    }
    Segment[[i]]$Znew <-Z
    K [[i]]= Segment[[i]]
    }

  }
K
}


fixZ(foot.trees[[23]])
foot.trees.newZ <-lapply(foot.trees, function(x) fixZ(x))

