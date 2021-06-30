#Function that returns the tree structure of every pattern
#The following information about every point in each tree is saved

# Index: Index of the point
# Connect: Index of its parent points (0 = basepoint)
# X,Y,Z coordinates of the point
# Order: The order of the points
# ep: boolean variable that is TRUE if the point is an endpoint
# Daughters : Number of daughters.
# SegLength : Segmen length connecting the point with its parent.
# Angle information of that segment
# Tree: Tree index
# Subjet ID, Age, Gender, Height, Weight,BMI, disease Group
# Blister and Sample IDs.

tree.structs <- function( sample ){
  
  #Which covariates to store
  covkeep <- c('SubjectID', 'Age', 'Gender', 'Height_cm', 
               'Weight_kg', 'BMI', 'Group')
  
  # Extract the trees that are fully in the image
  T0 <- get.basepoints( sample )$Tree
  epp <- get.nonorf.ep( sample )
  
  # Extract blister and sample 
  spl <- strsplit( sample$name, split = '[.]')
  blister <- as.numeric( spl[[1]][3] )
  sampnum <- as.numeric( gsub( 'pgp','',spl[[1]][4]) )
  
  # A distance function
  dfun <- function( x1, x2 ){
    
    if( length(x2)==3){
      
      sqrt( sum( (x1-x2)^2 ))
      
    }
    else{
      sqrt(rowSums( (x1[rep(1,nrow(x2)),]-as.matrix(x2) )^2))
    }
  }
  
  # Extract the information in the df that we need
  toUse <- which( sample$df$Tree %in% T0 )
  df0 <- sample$df[toUse,]
  
  # Same thing for the end point df
  toUse <- which( sample$endpoints$Tree %in% T0 )
  epdf <- sample$endpoints[toUse,]
  
  # Get the (3D) end point coordinates
  ep <- get.endpoints( sample, names = c('X','Y','Z'))
  ep <- ep[toUse,]
  trees <- unique( epdf$Tree )
  totd <- vector()
  
  # All trees in sample
  tree.list <- list()
  
  # Function to shorten distance computations
  distfun <- function( dfrow1, dfrow2 ){
    
    sqrt(sum((dfrow1[,c("X","Y","Z")]-dfrow2[,c("X","Y","Z")])^2))
    
  }
  
  count.idx <- 1
  for( t in trees ){
    
    dftmp <- df0[df0$Tree == t,]
    epdftmp <- epdf[epdf$Tree==t,]
    ep3 <- ep[epdf$Tree ==t,]
    
    epseg <- which( dftmp$Centripete == 0 )
    maxord <- max( dftmp$Order )
    tort <- dftmp$Tortuosity
    
    pt.idx <- point.indices( dftmp )
    ep.idx <- 1:length(ep3$X)+max(pt.idx)
    all.idx <- c(pt.idx, ep.idx)
    seg.length <- vector(length=length(all.idx))
    PlanarAngle <- XYAngle <- ZAngle <- rep(NA, length(all.idx))
    bpts <- cbind( dftmp$X, dftmp$Y, dftmp$Z)
    conn.to <- vector(length=length(all.idx))
    conn.to[all.idx==1] <- 0
    seg.length[all.idx==1] <- NA
    
    npts <- max( all.idx )
    
    #Connect end points to branch points
    
    if( npts == 2){
      conn.to[2] <- 1
      
      bpts <- data.frame( Index =pt.idx, 
                          Connect=conn.to[1:length(pt.idx)], 
                          X = bpts[,1],Y=bpts[,2],Z=bpts[,3],
                          Order=dftmp$Order,ep = rep(F,length(pt.idx)),
                          Daughters = 1)
      epts <- cbind( Index=ep.idx, 
                     Connect=conn.to[(length(pt.idx)+1):length(all.idx)],ep3,
                     Order=epdftmp$BranchOrder+1,
                     ep=rep(T,length(ep.idx)), Daughters = 0)
      
      pts <- rbind( bpts, epts)
      seg.length[2] <- distfun( pts[1,],pts[2,] )
      pts$SegLength <- seg.length
      pts <- unique( pts )
      pts$Tree <- rep(t,2)
      pts <- cbind( pts, sample$Covariates[rep(1,2),covkeep])
      pts$Blister <- rep(blister,2)
      pts$Sample <- rep( sampnum,2)
      PlanarAngle[2] <- dftmp$PlanarAngle
      pts$PlanarAngle <- PlanarAngle
      XYAngle[2] <- dftmp$XYAngle
      pts$XYAngle <- XYAngle
      ZAngle[2] <- dftmp$ZAngle
      pts$ZAngle <- ZAngle
      rownames(pts) <- NULL
      tree.list[[count.idx]] <- pts
      count.idx <- count.idx + 1
      
    }
    else{
      for( i in 1:length(ep.idx) ){
        
        possible <- which( dftmp$Order == epdftmp$BranchOrder[i] )
        xe <- cbind( ep3[i,1], ep3[i,2], ep3[i,3])
        xb <- bpts[ possible, ] 
        suggd <- dfun( xe, xb)
        
        if( length( possible )>1){  
          obsd <- dftmp$Length_microm[possible]/tort[possible]          
          best.match <- which.min( abs(suggd-obsd ) )
          possible <- possible[best.match]
          suggd <- suggd[best.match]
        }
        conn.to[ which(all.idx==ep.idx[i]) ] <- 
          pt.idx[possible]
        PlanarAngle[ which(all.idx==ep.idx[i])] <- dftmp$PlanarAngle[possible]
        XYAngle[ which(all.idx==ep.idx[i])] <- dftmp$XYAngle[possible]
        ZAngle[ which(all.idx==ep.idx[i])] <- dftmp$ZAngle[possible]
        seg.length[which(all.idx==ep.idx[i])] <- suggd
        
        
      }
      
      #Connect branch points
      
      for( i in 2:max(pt.idx) ){
        
        ind <- which( pt.idx == i)[1]
        ordtmp <- unique( dftmp$Order[pt.idx==i] )
        
        possible <- pt.idx[which( dftmp$Order == (ordtmp-1))] 
        poss.idx <- which(pt.idx%in%possible)
        xe <- cbind( bpts[ind,1], bpts[ind,2], bpts[ind,3])
        xb <- bpts[poss.idx , ] 
        suggd <- dfun( xe, xb)
        best.match <- 1
        if( length( possible)>1){
          
          obsd <- dftmp$Length_microm[poss.idx]/tort[poss.idx]
          best.match <- which.min( abs(suggd-obsd ) )
          possible <- pt.idx[poss.idx[best.match]]
          suggd <- suggd[best.match]
        }
        
        PlanarAngle[ which(all.idx==i)] <- dftmp$PlanarAngle[poss.idx[best.match]]
        XYAngle[ which(all.idx==i)] <- dftmp$XYAngle[poss.idx[best.match]]
        ZAngle[ which(all.idx==i)] <- dftmp$ZAngle[poss.idx[best.match]]
        conn.to[ which(all.idx==i) ] <- possible
        seg.length[which(all.idx==i)] <- unique( suggd )
        
        
      }
      
      bpts <- data.frame( Index =pt.idx, 
                          Connect=conn.to[1:length(pt.idx)], 
                          X = bpts[,1],Y=bpts[,2],Z=bpts[,3],
                          Order=dftmp$Order,ep = rep(F,length(pt.idx)),
                          Daughters = rep(0,length(pt.idx)))
      
      epts <- cbind( Index=ep.idx, 
                     Connect=conn.to[(length(pt.idx)+1):length(all.idx)],
                     ep3,
                     Order=epdftmp$BranchOrder+1,
                     ep=rep(T,length(ep.idx)), 
                     Daughters = rep(0,length(ep.idx)))  
      
      pts <- rbind( bpts, epts)
      pts$SegLength <- seg.length
      pts$PlanarAngle <- PlanarAngle
      pts$XYAngle <- XYAngle
      pts$ZAngle <- ZAngle
      pts <- unique( pts )
      rownames(pts) <- NULL
      Daughters <- vapply( pts$Index, function(x) 
        sum( c(pts$Connect) == x ), 1)
      pts$Daughters <- Daughters
      pts$Tree <- rep( t, length(Daughters))
      
      pts <- cbind( pts, sample$Covariates[rep(1,length(Daughters)),covkeep])
      pts$Blister <- rep(blister,length(Daughters))
      pts$Sample <- rep( sampnum, length( Daughters ))
      tree.list[[count.idx]] <- pts
      count.idx <- count.idx+1
      
      
    }
  }  
  if( !is.empty( tree.list )){
  return(tree.list)
  }
  else{
    
    data.frame( )
  }
}
