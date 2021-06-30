require( ENFdata )
require(circular)
require(MASS)
source( './tree.structs.R')
source( './get.nonorf.ep.R')
source( './point.indices.R')

#script to save statistics from data
data( Foot.all)
data(  covariatetable.all)
data(Calf.all)

namechange <- function(x){
  y <- get( x )
  y$name <- x
  y
  
}

footlist <- lapply(   Foot.all.names, 
                     function(x) namechange(x) )
calflist <- lapply(  Calf.all.names, 
             function(x) namechange(x) )

#number of base points
nbp <- vapply( footlist, function(x) get.basepoints.ppp(x)$n, 1 )
nbp_calf <-vapply( calflist, function(x) get.basepoints.ppp(x)$n, 1 )
#footlist = footlist[ nbp>10 ]
#footlist=footlist[-c(16,34:37)]

sidtmp <- vapply( footlist, function(x) x$Covariates$SubjectID, 1)
sidtmp_calf <- vapply( calflist, function(x) x$Covariates$SubjectID, 1)

disstmp <- vapply( sidtmp, function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'normal')
disstmp_calf <- vapply( sidtmp_calf, function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'normal')

dissub <- vapply( unique( sidtmp ), function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'normal')
dissub_calf<- vapply( unique( sidtmp_calf ), function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'normal')

sid <- unique( sidtmp ) #subject idFoot
sid_calf<- unique( sidtmp_calf )

# Tree stats contains the info on tree level (in a df)
# foot.trees the info on segment level (per sample, in a list)
foot.trees <- lapply( footlist, tree.structs )
calf.trees <- lapply(calflist, tree.structs)

tree.stats <- lapply( foot.trees, function(y) t(vapply( y, function(x) c(x$SubjectID[1],max(x$Order),sum(x$SegLength,na.rm=T),nrow(x)),c(1,1,1,1) )))
tree.stats_calf<- lapply( calf.trees, function(y) t(vapply( y, function(x) c(x$SubjectID[1],max(x$Order),sum(x$SegLength,na.rm=T),nrow(x)),c(1,1,1,1) )))

tree.stats <- do.call(  rbind, tree.stats )
tree.stats_calf <- do.call(  rbind, tree.stats_calf )

tree.stats <- data.frame(tree.stats)
tree.stats_calf<- data.frame(tree.stats_calf)

names( tree.stats ) <- c( 'SID', 'MaxOrd', 'SumSeg','NSeg')
names( tree.stats_calf ) <- c( 'SID', 'MaxOrd', 'SumSeg','NSeg')

tree.dis <- vapply( tree.stats[,1], function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'Normal' )
tree.dis_calf<- vapply( tree.stats_calf[,1], function(x) covariatetable.all$diseasestatus[covariatetable.all$SubjectID==x],'Normal' )

tree.dis <- factor( tree.dis )
tree.dis_calf <- factor( tree.dis_calf )

tree.dis <- factor( tree.dis, levels( tree.dis)[c(3,1,2,4)] )
tree.dis_calf<- factor( tree.dis_calf, levels( tree.dis_calf)[c(3,1,2,4)] )

tree.stats$diseasestatus <- tree.dis
tree.stats_calf$diseasestatus <- tree.dis_calf
#saveRDS(tree.stats,'./tree.stats.rds')
#saveRDS(tree.stats_calf,'./tree.stats_calf.rds')

