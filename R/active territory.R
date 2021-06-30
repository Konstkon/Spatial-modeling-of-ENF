source('./conv.R') #ACTIVE TERRITORY HELP FUNCTIONS
source('./noc_model_simulator.R') # NOC MODEL
source('./helpool.R') # HELP FUNCTION FOR POOLING

#unfortunatly the data and model params cannot be shared
#direction_list_mild = readRDS("./datarelated/model_params/mild_directions_noc.RDA")
#direction_list = readRDS("./datarelated/model_params/normal_directions_noc.RDA")

Reactive_normals = lapply ( seq_along(realepdf), function(x) reactive_territory_3D(realepdf[[x]],realbpdf[[x]]) )
Reactive_mild = lapply ( seq_along(realepdfmild), function(x) reactive_territory_3D(realepdfmild[[x]],realbpdfmild[[x]]) )

#SIMULATIONS FROM MODEL####
simulation = lapply(seq_along(realbp), 
                                   function(x) Noc_model(realbpdf[[x]],
                                                         params,
                                                         direction_list[[x]],
                                                         SID[x]))

simulation_mild = lapply(seq_along(realbpMILD), 
                         function(x) Noc_model(realbpdfmild[[x]],
                                               params_m,
                                               direction_list_mild[[x]],
                                               SIDM[x]))

#COMPUTE REACTIVE TERRITORIES####
Reactive_simulation_normals = lapply(seq_along(simulation), function(x)
     reactive_territory_3D(subset(simulation[[x]],simulation[[x]]$Endpoint==TRUE),
                           subset(simulation[[x]],simulation[[x]]$Basepoint==TRUE)))

Reactive_simulation_mild = lapply(seq_along(simulation_mild), function(x)
  reactive_territory_3D(subset(simulation_mild[[x]],simulation_mild[[x]]$Endpoint==TRUE),
                        subset(simulation_mild[[x]],simulation_mild[[x]]$Basepoint==TRUE)))


total_volume_healthy = unlist(lapply(Reactive_normals, function(x) sum(x)))
total_volume_healthy_sim = unlist(lapply(Reactive_simulation_normals, function(x) sum(x)))

total_volume_sim_mild = unlist(lapply(Reactive_simulation_mild, function(x) sum(x)))
total_volume_mild= unlist(lapply(Reactive_mild, function(x) sum(x)))


#### Subject wise 

Subject_EAT <- function(Volume,SID){
Nsamp = samplerpersubj( SID, unique(SID))
Change_Subj = changeofsubj( Nsamp ) # INDEX WHERE WE HAVE CHANGE IN SUBJECT
Subject_coverage = mean(Volume[1:Change_Subj[1]])
for (j in 2:length(unique(SID))){
Subject_coverage = c(Subject_coverage, mean(Volume[(Change_Subj[j-1]+1):Change_Subj[j]]))
}
Subject_coverage
}

Mild_subject_EAT = Subject_EAT(total_volume_mild, SIDM)
Healthy_subject_EAT = Subject_EAT(total_volume_healthy, SID)
Healthy_subject_EAT_sim = Subject_EAT(total_volume_healthy_sim, SID)



#FIGURE 6 ####
par(mfrow=c(1,2))
boxplot(list(mild=total_volume_mild/1000,Healthy=total_volume_healthy/1000),main='Volume of EATs',ylim=c(0,60))

shiftplot(total_volume_healthy/1000, total_volume_mild/1000, pch=20, 
          xlab='Healthy quantiles', ylab='Mild quantiles - Healthy quantiles')
abline(h=0)

#fIGURE 12####
shiftplot(total_volume_healthy/1000, total_volume_healthy_sim/1000, pch=20,
          xlab='Healthy quantiles', ylab='Model quantiles - Healthy quantiles')
abline(h=0)
shiftplot(total_volume_mild/1000,total_volume_sim_mild/1000, pch=20,
          xlab='Mild quantiles', ylab='Model quantiles - Mild quantiles')
abline(h=0)
