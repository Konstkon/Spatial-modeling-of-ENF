# Function that simulates Lengths of the first branches given group
# Lengths are simulated from a weibull distribution
# Num_bp =  Number of basepoints (  Number of lengths to be simulated)
# Group =  patien Group (i.e 'Normal', 'Mild')
# Returns lengths

Simulate_Length_fb_subject <- function( Num_bp , params){
  shape= params$shape
  rate= params$rate
  Lengths = rgamma( Num_bp , shape , rate=rate )
  Lengths
}


# Function that simulates angles of colatidute and longitude for the first branch
# It initially simulates points uniformly on disc and then compute z on a sphere.
# Then  Colatidute and Longitude are calculated
# Num_bp =  Number of basepoints (  Number of angles to be simulated)
# Returns a dataframe containing angles of colatidute and longitude

Simulate_Direction_fb_subject <- function ( Num_bp, par ){
  radius = sqrt ( runif ( Num_bp , 0 , 1 ) )
  k = par$kappa
  
  mean_dir= runif( 1 , 0 , 2 * pi )
  
  phi = rvonmises( Num_bp , mean_dir,k  )
  x = radius * cos( phi )
  y = radius * sin( phi )
  z = sqrt( 1 - x * x - y * y )
  theta = acos( z )
  df = data.frame( theta = theta , phi = phi )
  
}

# Function that calculates points given directions and lengths and mother point
# Length =  a vector containing Lengths
# angles =  a dataframe containing colatitude (theta) and longityde(phi ) angles
# Motherpoints = A data frame containing the cartesian cooridinates of the mother points
# Returns a  dataframe with the points of interest

Simulate_points <- function( Length , angles, Motherpoints ) {
  Z = Length * cos (angles$theta) + Motherpoints$Z
  X = Length * sin( angles$theta ) * cos( angles$phi ) + Motherpoints$X
  Y = Length * sin( angles$theta ) * sin( angles$phi ) + Motherpoints$Y
  df=data.frame( X = X, Y = Y , Z = Z , Tree = Motherpoints$Tree, Basepoint = FALSE )
  df
}

## offspring by subject
Simulate_offspring_subject<- function( prob,size ){
  offspring = rnbinom(1 , size , prob )
  offspring
}


# Function that simulates Lengths between first branch and end points given group
# Lengths are simulated from gamma distribution
# Num_ep =  Number of endpoints (  Number of lengths to be simulated)
# Group =  patien Group (i.e 'Normal', 'Mild')
# Returns lengths

Simulate_Length_ep_subject <- function( Num_ep , Group ,par){
  shape= par$shape_ep
  rate = par$rate_ep
  
  Length = rgamma(Num_ep , shape , rate=rate  )
  while(max(Length>60)){
    Length = rgamma(Num_ep , shape , rate=rate  )
  }
  Length
}


# Function that simulates angles between first branch and end points given group
# Angles are simulated from a Schladitz distribution
# Num_ep =  Number of endpoints (  Number of lengths to be simulated)
# Group =  patien Group (i.e 'Normal', 'Mild')
# Returns angles
Simulate_Direction_ep_subject <- function ( Num_ep , Group, par ){
  beta=par$beta
  Num_sim = Num_ep #- Num_pi2
  # Main simulation
  ang = data.frame(theta=c(),phi=c())
  angs = SimulateShladitz(Num_sim,beta,  c(0,0))
   angs
}