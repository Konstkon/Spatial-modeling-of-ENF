point.indices <- function( df ){
  
  coords <- cbind( df$X, df$Y, df$Z )
  
  pt.idx <- vector( length=nrow( df ))
  tosort <- 1:nrow(df)
  
  pt.idx[df$Order==1] <- 1
  tosort <- tosort[-which(df$Order==1)]
  
  ordtmp <- 2
  idx <- 2
  while( !is.empty(tosort)){
    
    idtmp <- which( df$Order == ordtmp )
    
    while( !is.empty( idtmp )){
      
      if( length(idtmp)<=2){
        pt.idx[idtmp] <- idx
        idx <- idx+1
        tosort <- tosort[-which(tosort%in%idtmp)]
       
        idtmp <- vector()
      }
      else{
        idtmp1 <- idtmp[1]
        idtmp<-idtmp[-1]
        tosort <- tosort[-which(tosort==idtmp1)]
        pt.idx[idtmp1] <- idx
        dupl <- which(apply( coords[idtmp,], 1, 
                             function(x) identical(x,coords[idtmp1,])) )
        if( !is.empty( dupl )){
        pt.idx[idtmp[dupl]] <- idx
        tosort <- tosort[-which(tosort%in%idtmp[dupl])]
        idtmp <- idtmp[-dupl]
        }
        idx <- idx+1
      }
    }
    
    
    ordtmp <- ordtmp+1
  }  
  
  pt.idx
}
