source('./Angular_Distributions.R') # code to simulate from Schladitz
source('./simulatorhelp.R')

#load parametrs
params_m = readRDS("./datarelated/model_params/params_subject_mild.RDA")
params = readRDS("./datarelated/model_params/subject_params_normal.RDA")
# Function that simulates endpoints of ENF given Base points
# Basepoints = Base point of an ENF object 
# Group = Patient Group ('Normal' or 'Mild')
# Returns a Dataframe.

SimulateENFdata_subject <- function ( Basepoints , Group , parameters, ID )
{
  end_p = data.frame()
  Basepoints= Basepoints[,1:4]
  # Simulating first branch point
  Basepoints$Basepoint = TRUE
  Basepoints$Endpoint = FALSE
  Basepoints$Branchpoint = FALSE
  par = subset(parameters, parameters$SID == ID)
  Num_bp =    length(Basepoints$X)                                       # Number of basepoints
  Length_fb = Simulate_Length_fb_subject( Num_bp , par)                     # Simulates N lengths
  angles_fb = Simulate_Direction_fb_subject ( Num_bp, par  )                        # Simulates N directions
  fb_points = Simulate_points( Length_fb , angles_fb , Basepoints)     # Simulates points for the first branch
  fb_points$Endpoint = FALSE
  fb_points$Branchpoint = TRUE
  
  # Simulating end points for every first branch point
  
  for (k in 1:Num_bp) {
    Num_ep =     Simulate_offspring_subject( par$prob,par$size )  + 1                           # Number of end points for branch point k
    if (Num_ep > 1 ){
      Length_ep =  Simulate_Length_ep_subject ( Num_ep , Group ,par)                    # Simulates Lengths
      angles_ep =  Simulate_Direction_ep_subject ( Num_ep , Group, par )                 # Simulates Directions
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

