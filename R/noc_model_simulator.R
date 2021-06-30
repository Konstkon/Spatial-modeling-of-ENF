#Lengths first branch
Simulate_length_first_branch <- function( Num_bp , params){
  shape= params$shape
  rate= params$rate
  Lengths = rgamma( Num_bp , shape , rate=rate )
  
  Lengths
  }

#Angles first branch
Simulate_Direction_first_branch <- function ( Num_bp,phi ,par ){
  radius = sqrt ( runif ( Num_bp , 0 , 1 ) )
  kappa = par$kappa
  phi_ = c()
  for (l in 1:Num_bp){
  phi_ = c(phi_, as.numeric(rvonmises(1, phi[l], kappa)))
  }
  x = radius * cos( phi_ )
  y = radius * sin( phi_ )
  z = sqrt( 1 - x * x - y * y )
  theta = acos( z )
  df = data.frame( theta = theta , phi = phi_ )
  df
}

#Simulatepoints
Simulate_points <- function( Length , angles, Motherpoints ) {
  Z = Length * cos (angles$theta) + Motherpoints$Z
  X = Length * sin( angles$theta ) * cos( angles$phi ) + Motherpoints$X
  Y = Length * sin( angles$theta ) * sin( angles$phi ) + Motherpoints$Y
  df=data.frame( X = X, Y = Y , Z = Z , Tree = Motherpoints$Tree, Basepoint = FALSE )
  df
}
#offsprings
Simulate_offspring_subject <- function( prob,size ){
  offspring = rnbinom(1 , size , prob )
  offspring
}


Simulate_Length_endpoints <- function( Num_ep  ,par){
  shape= par$shape_ep
  rate = par$rate_ep
 Length = rgamma(Num_ep , shape , rate=rate  )
  while(max(Length>60)){
    Length = rgamma(Num_ep , shape , rate=rate  )
  }
  Length
}

#andlge ends
Simulate_Direction_endpoints <- function ( Num_ep , par ){
  beta=par$beta
  Num_sim = Num_ep 
  ang=data.frame(theta=c(),phi=c())
   angs= SimulateShladitz(Num_sim,beta,  c(0,0))
   angs
}
# most recent simulate noc function


Noc_model <- function ( Basepoints , parameters,directions, ID )
{
  end_p = data.frame()
  Basepoints= Basepoints[,1:4]
  # Simulating first branch point
  Basepoints$Basepoint = TRUE
  Basepoints$Endpoint = FALSE
  Basepoints$Branchpoint = FALSE
   par = subset(parameters, parameters$SID == ID)
  Num_bp =    length(Basepoints$X)                                       # Number of basepoints
  Length_fb = Simulate_length_first_branch( Num_bp , par)                     # Simulates N lengths
  angles_fb = Simulate_Direction_first_branch( Num_bp,directions,par  )                        # Simulates N directions
  fb_points = Simulate_points( Length_fb , angles_fb , Basepoints)     # Simulates points for the first branch
  fb_points$Endpoint = FALSE
  fb_points$Branchpoint = TRUE
  
  # Simulating end points for every first branch point
  
  for (k in 1:Num_bp) {
    Num_ep =     Simulate_offspring_subject( par$prob,par$size )  + 1                           # Number of end points for branch point k
    if (Num_ep > 1 ){
      Length_ep =  Simulate_Length_endpoints ( Num_ep ,par)                    # Simulates Lengths
      angles_ep =  Simulate_Direction_endpoints ( Num_ep , par )                 # Simulates Directions
      end_points = Simulate_points( Length_ep , angles_ep , fb_points [k,]  ) # Simulates end points
      end_points$Endpoint = TRUE
      end_points$Branchpoint = FALSE
      end_p = rbind ( end_p , end_points )
    }
    else {
      fb_points[k,]$Endpoint = TRUE
      
    }
  }
  DF = rbind(end_p,fb_points, Basepoints)
  DF
}


#HOW TO SIMULATE
#
#  simulation = lapply(seq_along(realbp), function(x) Noc_model(realbpdf[[x]],params,direction_list[[x]],SID[x]))
#

