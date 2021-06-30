#Help functions to get the non orphan points from the patterns
# i.e endpoints where their basepoint is in the window



#Returns a ppp object of the non orphan enpoints in the pattern

get.nonorf.ep <- function (x, X = "X", Y = "Y", unique = FALSE) {
  
  
  Endpoints <- get.endpoints(x, names = c("X", "Y",'Tree'))
  Basepoints <- get.basepoints( x )
  Endpoints <- Endpoints[ which(Endpoints$Tree %in% Basepoints$Tree), ]
  
  
  if( unique==TRUE ){
    unique( ppp(x = Endpoints[, X], y = Endpoints[, Y], window = x$window) )
  }
 else{ppp(x = Endpoints[, "X"], y = Endpoints[, "Y"], window = x$window) }
  
}

#Function that returns a dataframe of the non-orphan enpoints(2D) in the patterns

get.nonorf.ep.df <- function (x, X = "X", Y = "Y", unique = F) {
  
  
  Endpoints <- get.endpoints(x, names = c(X, Y,'Tree'))
  Basepoints <- get.basepoints( x )
  Endpoints <- Endpoints[ which(Endpoints$Tree %in% Basepoints$Tree), ]
  
  
  Endpoints
  
}

#Function that returns a dataframe of the non-orphan enpoints(3D) in the patterns

get.nonorf.ep.df3d <- function (x, X = "X", Y = "Y",Z="Z", unique = F) {
  
  
  Endpoints <- get.endpoints(x)
  Basepoints <- get.basepoints( x )
  Endpoints <- Endpoints[ which(Endpoints$Tree %in% Basepoints$Tree), ]
  
  
  Endpoints
  
}

#Function that returns a dataframe of the non-orphan basepoints(3D) in the patterns

get.nonorf.bp<-function (x, X = "X", Y = "Y",Z="Z", unique = F){
  Endpoints <- get.endpoints(x)
  Basepoints <- get.basepoints( x )
  Trees <- Endpoints$Tree
  Basepoints <- Basepoints[ which(Basepoints$Tree %in% Trees), ]
  Basepoints
}

#Returns a ppp object of the non orphanan basepoints(2D) in the pattern

get.nonorf.bp.ppp<-function (x, X = "X", Y = "Y",Z="Z", unique = F){
  Endpoints <- get.endpoints(x)
  Basepoints <- get.basepoints( x )
  Trees <- Endpoints$Tree
  Basepoints <- Basepoints[ which(Basepoints$Tree %in% Trees), ]
  if( unique ){
    unique( ppp(x = Basepoints[, X], y = Basepoints[, Y], window = x$window) )
  }
  else{ppp(x = Basepoints[, X], y = Basepoints[, Y], window = x$window) }
}