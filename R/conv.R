setwd("C:/Users/konkons/Desktop/ENF/ENFdata/Konstantinos") # Chnage path
source('./Data_load.R')  #Load the data
require(geometry)


# Computes the 3D reactive territory
# df: dataframe containing information about the end points
# df2: dataframe with information about the base points
# returns: volume of the 3D reactive territory.

reactive_territory_3D <- function( df, df2 ){
  
treeunique = unique(df$Tree)
names=c('X','Y','Z')
volume_react = c()

for (j in treeunique){  # loop over trees
ep_same_tree = as.matrix(df[df$Tree==j,names])  #get ep
bp_of_tree = as.matrix(df2[df2$Tree==j,names])  #get bp
Tree = rbind(ep_same_tree , bp_of_tree)         #get points in tree
# IF there are more than two end points in a tree
# compute the volume of the delaunayn triangulation of those points
if (length(Tree[,1])>3){
del = delaunayn( Tree , output.options = 'Fa')
volume_react = c( volume_react , sum(del$areas) )  
}
# if there are two end points in a tree
# compute the volume of the cone with diameter distance between ep 
# and heigth distance between ep and x1+x2/2 xi ep.

if (length(Tree[,1])==3){
  coords = c(mean(ep_same_tree[,1]),mean(ep_same_tree[,2]),mean(ep_same_tree[,3])) + c(dist(ep_same_tree)/4,0,0)
  coords2 = c(mean(ep_same_tree[,1]),mean(ep_same_tree[,2]),mean(ep_same_tree[,3]))-  c(dist(ep_same_tree)/4,0,0)
  Tree1=rbind(Tree,coords,coords2)
  del = delaunayn( Tree1 , output.options = 'Fa')
  volume_react = c( volume_react , sum(del$areas) ) 
}

# if there is only one ep in the tree
# compute the distance between the ep and its mother

if (length(Tree[,1]) == 2){
  volume_react  = c( volume_react , dist(Tree))
}

}
volume_react
}

#Get reactive territories for every sample
Reactive_normals = lapply ( seq_along(realepdf), function(x) reactive_territory_3D(realepdf[[x]],realbpdf[[x]]) )
Reactive_mild = lapply ( seq_along(realepdfmild), function(x) reactive_territory_3D(realepdfmild[[x]],realbpdfmild[[x]]) )

#Compute the total volume of the EATs per sample
total_volume_healthy = unlist(lapply(Reactive_normals, function(x) sum(x)))
total_volume_mild= unlist(lapply(Reactive_mild, function(x) sum(x)))

#Add volumes of EATs as marks to the basepoint patterns(both mild and healthy patterns)
realbpm =lapply(seq_along(realbpm), function(x) setmarks(realbpm[[x]], 
                                                         cbind(marks(realbpm[[x]]),data.frame(Volume=Reactive_normals[[x]]))))
realbpmmild =lapply(seq_along(realbpmmild), function(x) setmarks(realbpmmild[[x]], 
                                                                 cbind(marks(realbpmmild[[x]]),data.frame(Volume=Reactive_mild[[x]]))))
