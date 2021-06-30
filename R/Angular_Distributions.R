#3d distribtuions
Smatrix <- function(df){
  res = matrix( 0 , 3 , 3 )
  Size = length( df[ , 1 ] )
  for ( k in 1:Size ){
    incr = t( df[ 1 , ] ) %*% as.matrix( df[ 1 , ] )
    res = res + incr
    
  }
  res = res / Size
  res
}

#Simulate from Schladitz distribution
SimulateShladitz <- function( N ,  beta , direction ){
  
  theta  = c( )
  phii = c( )
  eps = 1e-8
  # MAIN SIMULATION STEP
  for (i in 1: N){
    R1 = runif( 1 , 0 , 1 )
    R2 = runif( 1 , 0 , 1 )
    phi = 2 * pi * R1
    ksi = 2 * R2 - 1
    eta = sqrt ( 1 - ksi * ksi )
    if ( abs( beta - 1) > eps ){
      norm = sqrt ( ksi * ksi - ksi * ksi * beta * beta + beta * beta )
      ksi = ksi / norm
      eta = eta * beta  / norm
    }
    theta = c( theta , acos( ksi ))
    phii = c( phii , phi )
  }
   res =data.frame(theta = theta, phi = phii)
}


#Estimate beta using mle
betamle <- function ( N , theta, betainit ){
  beta = betainit
  for (i in 1:N){
    Xib = beta * beta* cos(theta)*cos(theta)/(1+(beta*beta-1)*cos(theta)*cos(theta))
    gb = 3* sum( Xib )-length(Xib)
    Yib = 3/2*beta*sin(2*theta)*sin(2*theta)/(1+(beta*beta-1)*cos(theta)*cos(theta))^2
    gbprime = sum ( Yib)
    beta = beta - gb / gbprime
    
  }
  
  beta
}

