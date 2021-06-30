#Load some code
source('./Helpfunctions.R')  # CONTAINS HELP FUNCTIONS CREATED TO ANALYSE THE DATA
#Load libraries
library(extraDistr)
library(spptest)

#########################################################################################

# This loads data

'
Base points:

realbpdf : Dataframe with information about base points
realbpdfmild : Dataframe with information about base points mild group

nbp : number of basepoints per sample
nbpmild :number of bp per sample mild group

realbp : marked ppp of base points with marks type and Tree
realbpMILD :ppp for mild group

realbp3D : marked pp3 of base points with marks type and Tree in 3d
realbp3Dmild : pp3 for bp of the mild group with marks tree and tupe

realbpm: marked ppp of base points with more marks of interest I.e volume,avgdiameter,max order ...
realbpmmild :marked ppp of base points with more marks of interest I.e volume,avgdiameter,max order ... for mild
realbpm3D: similar for pp3 healthy group( marked in 3D)
realbpm3Dmild : similar for pp3 mild group




End points:

realepdf: Dataframe with information about end points
realepdfmild:Dataframe with information about end points for mild

nep : number of endpoints per sample
nepmild : number of endpoints per sample mild group

realep : marked ppp of end points with marks type and Tree
realepmild: ppp for mild group

realep3D: marked pp3 of end points with marks type and Tree
realep3Dmild :marked pp3 of end points with marks type and Tree mild group

realepm: marked ppp of end points with more marks of interest I.e distance and order ...



Point pattern:

realpp: marked ppp, basically a superposition of realbp and realep
realppmild :ppp, basically a superposition of realbpMILD and realepmild

realpp3D:marked pp3, basically a superposition of realbp3D and realep3D
realppmild3D:marked pp3, basically a superposition of realbp3Dmild and realep3Dmild

dATA:
footlistmild: list that contains the data for mild group
footlist: list that contains the data for normal patients
SID : SUBJECT ID
SIDM: SUBJECT ID MILD
'



############################



#load data
data(Foot.allnormals)
dat = Foot.allnormals.names

footlist <- lapply(   dat,
                      function(x) namechange(x) )

#consider only samples with basepoint in (4,46)
nbp <- vapply( footlist, function(x) get.basepoints.ppp(x)$n, 1 )
footlist = footlist[ nbp>10 ]
footlist=footlist[-c(16,34:37)]
#footlist = footlist[ nbp>10 ]
realep = lapply (footlist , function ( x) get.nonorf.ep(x) )
nep = vapply(realep, function(x) getnum_ep(x),1)
nbp <- vapply( footlist, function(x) get.basepoints.ppp(x)$n, 1 )

# get real basepoints + non orphan endpoints = observed pattern
# We create marked point processes with two marks
# mark Tree corresponds to the index of the tree
# mark type correspond to the type of point (base or end)
# i.e type=1 implies that this point is a base point

realbpdf = lapply ( footlist, function(x) get.basepoints.marks(x) ) #DATA FRAME FOR BP
realbp3D = create_pp3_from_df(realbpdf,1)
realbp = lapply ( footlist, function(x) get.basepoints.ppp(x) )     #PPP
realbp = lapply(seq_along(realbp),function(x) setmarks(realbp[[x]],
                                                       data.frame(Tree=realbpdf[[x]]$Tree,type = ones(realbp[[x]]$n,1))))


#branchpoints
realbranchdf = lapply ( footlist, function(x) subset(x$df,x$df$Order==2) ) #DATA FRAME FOR BP
realbranch = lapply(seq_along(footlist), function(x) ppp(realbranchdf[[x]]$X,realbranchdf[[x]]$Y,window = footlist[[x]]$window))



remove_dub = lapply(realbranch , function(x) which(duplicated(x)==FALSE))
realbranchdf = lapply(seq_along(realbranch) , function(x) realbranchdf[[x]][remove_dub[[x]],])
realbranch = lapply(seq_along(footlist), function(x) ppp(realbranchdf[[x]]$X,realbranchdf[[x]]$Y,window = footlist[[x]]$window))
realbranch3D = create_pp3_from_df(realbranchdf,2)

realbranch_new = lapply(seq_along(realbranch),
                        function(x)superimpose(realbranch[[x]],attr(realbranch[[x]],"rejects")))
realbranch_new = lapply(seq_along(1:112),function(x) setmarks(realbranch_new[[x]],
                                                              data.frame(Tree=realbranchdf[[x]]$Tree,type = 2*ones(realbranch_new[[x]]$n,1))))


realepdf =lapply (footlist , function ( x) get.nonorf.ep.df.marks(x) ) #DATA FRAME FOR EP
realep3D = create_pp3_from_df(realepdf,0)
realbpord = lapply(realepdf, function(x) get_tree_order(x))            #TREE ORDER

#MARKED BASE POINTS

realbpm = create_mark_bp(realbp, realbpdf, realbpord,2)
realbpm3D = create_mark_bp(realbp3D, realbpdf, realbpord,3)

realepord = lapply(realepdf, function(x) get_end_order(x))   #EP ORDER
realep = lapply (footlist , function ( x) get.nonorf.ep(x) ) #PPP
realepdf =lapply (footlist , function ( x) get.nonorf.ep.df.marks(x) ) #DF FOR EP
realep = lapply(seq_along ( realep) , function(x) setmarks(realep[[x]],
                                                           data.frame(Tree= realepdf[[x]]$Tree, type =zeros(realep[[x]]$n,1))))
realepm = lapply(seq_along ( realep) , function(x) setmarks(realep[[x]], 
                                                            data.frame(Tree= realepdf[[x]]$Tree, type =zeros(realep[[x]]$n,1), Distance_microm=realepdf[[x]]$Distance_microm, Order=realepord[[x]]$Order)))
realpp = lapply (seq_along ( realep) , function (x) createmarkpppdat(realbp[[x]],realep[[x]]))

realpp3D = lapply(seq_along (realep3D), function(x) superimpose_pp3(realbp3D[[x]],realep3D[[x]]))


SID = unlist(lapply ( footlist, function(x) GetsubjectID(x))) #Subject ID

################################### MILD GROUP
data(Foot.all)
foot.allmild = list(Foot.151.1.pgp1,Foot.151.1.pgp2, Foot.151.2.pgp1,Foot.151.2.pgp2,Foot.249.1.pgp1,Foot.249.1.pgp2,Foot.255.1.pgp1,Foot.255.1.pgp2,Foot.255.2.pgp1,Foot.255.2.pgp2,
                    Foot.261.1.pgp1,Foot.261.1.pgp2,Foot.264.1.pgp1,Foot.264.2.pgp1,Foot.264.2.pgp2,Foot.264.1.pgp2,Foot.273.1.pgp1,Foot.273.1.pgp2,Foot.273.1.pgp3,Foot.273.1.pgp4,
                    Foot.275.1.pgp1,Foot.275.2.pgp1,Foot.275.1.pgp2,Foot.275.2.pgp2,Foot.286.1.pgp1,Foot.286.2.pgp1, Foot.286.1.pgp2,Foot.286.2.pgp2)

nbpmild = vapply(foot.allmild, function(x) get.basepoints.ppp(x)$n,1)
footlistmild = foot.allmild
realepmild =lapply(footlistmild, function(x) get.nonorf.ep(x))
realepdfmild = lapply ( footlistmild, function(x) get.nonorf.ep.df.marks(x) )

realepmild = lapply(seq_along ( realepmild) , function(x) setmarks(realepmild[[x]],
                                                           data.frame(Tree= realepdfmild[[x]]$Tree,
                                                                      type =zeros(realepmild[[x]]$n,1))))

nepmild = vapply(realepmild, function(x) getnum_ep(x),1)
realbpdfmild = lapply ( footlistmild, function(x) get.basepoints.marks(x) )

realbpMILD = lapply ( footlistmild, function(x) get.basepoints.ppp(x) )
realbpMILD = lapply(seq_along(realbpMILD),function(x) setmarks(realbpMILD[[x]],
                                                       data.frame(Tree=realbpdfmild[[x]]$Tree,
                                                                  type = ones(realbpMILD[[x]]$n,1))))

realppmild = lapply (footlistmild , function (x) createmarkpppdat(get.basepoints.ppp(x),get.nonorf.ep(x)))
#MILD
realbpdfmild = lapply ( footlistmild, function(x) get.basepoints.marks(x) )
SIDM = unlist(lapply ( footlistmild, function(x) GetsubjectID(x)))
realbp3Dmild = create_pp3_from_df(realbpdfmild, 1)
realepdfmild = lapply ( footlistmild, function(x) get.nonorf.ep.df.marks(x) )
realep3Dmild = create_pp3_from_df(realepdfmild, 0)

realbpordMILD = lapply(realepdfmild, function(x) get_tree_order(x))
realepordMILD = lapply(realepdfmild, function(x) get_end_order(x))

realbpmmild = create_mark_bp(realbpMILD, realbpdfmild, realbpordMILD,2)
realbpm3Dmild = create_mark_bp(realbp3Dmild , realbpdfmild, realbpordMILD,3)


realpp3Dmild = lapply(seq_along (realep3Dmild), function(x) superimpose_pp3(realbp3Dmild[[x]],realep3Dmild[[x]]))

#mild branch points
realbranchdfmild = lapply ( footlistmild, function(x) subset(x$df,x$df$Order==2) ) #DATA FRAME FOR BP
realbranchmild = lapply(seq_along(footlistmild), function(x) ppp(realbranchdfmild[[x]]$X,realbranchdfmild[[x]]$Y,window = footlistmild[[x]]$window))
remove_dub_mild = lapply(realbranchmild , function(x) which(duplicated(x)==FALSE))
realbranchdfmild = lapply(seq_along(realbranchmild) , function(x) realbranchdfmild[[x]][remove_dub_mild[[x]],])
realbranchmild = lapply(seq_along(footlistmild), function(x) ppp(realbranchdfmild[[x]]$X,realbranchdfmild[[x]]$Y,window = footlistmild[[x]]$window))
realbranch3Dmild = create_pp3_from_df(realbranchdfmild,2)

#mild branch points




##moderate cases and severe cases
data("covariatetable.all")


footlistmoderate = list(Foot.126.1.pgp1, Foot.126.1.pgp2,
                        Foot.127.1.pgp1 , Foot.127.1.pgp2, Foot.127.2.pgp1, Foot.127.2.pgp2,
                        Foot.145.1.pgp1, Foot.145.1.pgp2, Foot.145.2.pgp1, Foot.145.2.pgp2, 
                        Foot.148.1.pgp1, Foot.148.1.pgp2, Foot.148.2.pgp1, Foot.148.2.pgp2,
                        Foot.223.1.pgp1, Foot.223.1.pgp2,
                        Foot.258.1.pgp1, Foot.258.1.pgp2, Foot.258.2.pgp2,
                        Foot.259.1.pgp1, Foot.259.1.pgp2, Foot.259.2.pgp1, Foot.259.2.pgp2,
                        Foot.276.1.pgp1, Foot.276.1.pgp2, Foot.276.2.pgp1, Foot.276.2.pgp2,
                        Foot.278.1.pgp1, Foot.278.1.pgp2, Foot.278.2.pgp1, Foot.278.2.pgp2,
                        Foot.282.1.pgp1, Foot.282.1.pgp2, Foot.282.2.pgp1, Foot.282.2.pgp2)


nbpmod = vapply(footlistmoderate, function(x) get.basepoints.ppp(x)$n,1)
footlistmoderate = footlistmoderate[ nbpmod>7 ]


realepmoderate =lapply(footlistmoderate, function(x) get.nonorf.ep(x))
realepdfmoderate = lapply ( footlistmoderate, function(x) get.nonorf.ep.df.marks(x) )

realepmoderate = lapply(seq_along ( realepmoderate) , function(x) setmarks(realepmoderate[[x]],
                                                                   data.frame(Tree= realepdfmoderate[[x]]$Tree,
                                                                              type =zeros(realepmoderate[[x]]$n,1))))

nepmod = vapply(realepmoderate, function(x) getnum_ep(x),1)
realbpdfmoderate = lapply ( footlistmoderate, function(x) get.basepoints.marks(x) )

realbpmod = lapply ( footlistmoderate, function(x) get.basepoints.ppp(x) )
realbpmod = lapply(seq_along(realbpmod),function(x) setmarks(realbpmod[[x]],
                                                               data.frame(Tree=realbpdfmoderate[[x]]$Tree,
                                                                          type = ones(realbpmod[[x]]$n,1))))

realppmod = lapply (footlistmoderate , function (x) createmarkpppdat(get.basepoints.ppp(x),get.nonorf.ep(x)))
#Moderate
SIDMod = unlist(lapply ( footlistmoderate, function(x) GetsubjectID(x)))
realbp3Dmod = create_pp3_from_df(realbpdfmoderate, 1)
realepdfmoderate = lapply ( footlistmoderate, function(x) get.nonorf.ep.df.marks(x) )
realep3Dmod = create_pp3_from_df(realepdfmoderate, 0)

realbpordmod = lapply(realepdfmoderate, function(x) get_tree_order(x))
realepordmod = lapply(realepdfmoderate, function(x) get_end_order(x))

realbpmmod = create_mark_bp(realbpmod, realbpdfmoderate, realbpordmod,2)
realbpm3Dmod = create_mark_bp(realbp3Dmod , realbpdfmoderate, realbpordmod,3)


realpp3Dmod = lapply(seq_along (realep3Dmod), function(x) superimpose_pp3(realbp3Dmod[[x]],realep3Dmod[[x]]))

#moderate branch points
realbranchdfmod = lapply ( footlistmoderate, function(x) subset(x$df,x$df$Order==2) ) #DATA FRAME FOR BP
realbranchmod = lapply(seq_along(footlistmoderate), function(x) ppp(realbranchdfmod[[x]]$X,realbranchdfmod[[x]]$Y,
                                                                    window = footlistmoderate[[x]]$window))
remove_dub_mod = lapply(realbranchmod , function(x) which(duplicated(x)==FALSE))
realbranchdfmod = lapply(seq_along(realbranchmod) ,
                         function(x) realbranchdfmod[[x]][remove_dub_mod[[x]],])
realbranchmod = lapply(seq_along(footlistmoderate),
                       function(x) ppp(realbranchdfmod[[x]]$X,realbranchdfmod[[x]]$Y,
                                       window = footlistmoderate[[x]]$window))
realbranch3Dmod = create_pp3_from_df(realbranchdfmod,2)



####
#severe
footlistsev = list( Foot.251.1.pgp1, Foot.251.1.pgp2, Foot.251.2.pgp1, Foot.251.2.pgp2,
                       Foot.257.1.pgp1, Foot.257.1.pgp2, Foot.257.2.pgp1, Foot.257.2.pgp2,
                       Foot.263.1.pgp1, Foot.263.1.pgp2, Foot.263.2.pgp1, Foot.263.2.pgp2,
                       Foot.265.1.pgp1, Foot.265.1.pgp2, Foot.265.2.pgp1, Foot.265.2.pgp2,
                       Foot.269.1.pgp1, Foot.269.1.pgp2, Foot.269.2.pgp1, Foot.269.2.pgp2,
                       Foot.271.1.pgp1, Foot.271.1.pgp2, Foot.271.2.pgp1, Foot.271.2.pgp2
                       )



nbpsev = vapply(footlistsev, function(x) get.basepoints.ppp(x)$n,1)
footlistsev = footlistsev[ nbpsev>7 ]


realepsev =lapply(footlistsev, function(x) get.nonorf.ep(x))
realepdfsev = lapply ( footlistsev, function(x) get.nonorf.ep.df.marks(x) )

realepsev = lapply(seq_along ( realepsev) , function(x) setmarks(realepsev[[x]],
                                                                           data.frame(Tree= realepdfsev[[x]]$Tree,
                                                                                      type =zeros(realepsev[[x]]$n,1))))

nepsev = vapply(realepsev, function(x) getnum_ep(x),1)
realbpdfsev = lapply ( footlistsev, function(x) get.basepoints.marks(x) )

realbpsev = lapply ( footlistsev, function(x) get.basepoints.ppp(x) )
realbpsev = lapply(seq_along(realbpsev),function(x) setmarks(realbpsev[[x]],
                                                             data.frame(Tree=realbpdfsev[[x]]$Tree,
                                                                        type = ones(realbpsev[[x]]$n,1))))

realppsev = lapply (footlistsev , function (x) createmarkpppdat(get.basepoints.ppp(x),get.nonorf.ep(x)))
#Moderate
SIDsev = unlist(lapply ( footlistsev, function(x) GetsubjectID(x)))
realbp3Dsev = create_pp3_from_df(realbpdfsev, 1)
realepdfsev = lapply ( footlistsev, function(x) get.nonorf.ep.df.marks(x) )
realep3Dsev = create_pp3_from_df(realepdfsev, 0)

realbpordsev = lapply(realepdfsev, function(x) get_tree_order(x))
realepordsev = lapply(realepdfsev, function(x) get_end_order(x))

realbpmsev = create_mark_bp(realbpsev, realbpdfsev, realbpordsev,2)
realbpm3Dsev = create_mark_bp(realbp3Dsev , realbpdfsev, realbpordsev,3)


realpp3Dsev = lapply(seq_along (realep3Dsev), function(x) superimpose_pp3(realbp3Dsev[[x]],realep3Dsev[[x]]))

#moderate branch points
realbranchdfsev = lapply ( footlistsev, function(x) subset(x$df,x$df$Order==2) ) #DATA FRAME FOR BP
realbranchsev = lapply(seq_along(footlistsev), function(x) ppp(realbranchdfsev[[x]]$X,realbranchdfsev[[x]]$Y,
                                                                    window = footlistsev[[x]]$window))
remove_dub_sev = lapply(realbranchsev , function(x) which(duplicated(x)==FALSE))
realbranchdfsev = lapply(seq_along(realbranchsev) ,
                         function(x) realbranchdfsev[[x]][remove_dub_sev[[x]],])
realbranchsev = lapply(seq_along(footlistsev),
                       function(x) ppp(realbranchdfsev[[x]]$X,realbranchdfsev[[x]]$Y,
                                       window = footlistsev[[x]]$window))
realbranch3Dsev = create_pp3_from_df(realbranchdfsev,2)

A=" The following variables are created

Base points:
  
realbpdf : Dataframe with information about base points
realbpdfmild : Dataframe with information about base points mild group

nbp : number of basepoints per sample
nbpmild :number of bp per sample mild group

realbp : marked ppp of base points with marks type and Tree
realbpMILD :ppp for mild group

realbp3D : marked pp3 of base points with marks type and Tree in 3d
realbp3Dmild : pp3 for bp of the mild group with marks tree and tupe

realbpm: marked ppp of base points with more marks of interest I.e volume,avgdiameter,max order ...
realbpmmild :marked ppp of base points with more marks of interest I.e volume,avgdiameter,max order ... for mild
realbpm3D: similar for pp3 healthy group( marked in 3D)
realbpm3Dmild : similar for pp3 mild group

End points:
  
realepdf: Dataframe with information about end points
realepdfmild:Dataframe with information about end points for mild

nep : number of endpoints per sample
nepmild : number of endpoints per sample mild group

realep : marked ppp of end points with marks type and Tree
realepmild: ppp for mild group

realep3D: marked pp3 of end points with marks type and Tree
realep3Dmild :marked pp3 of end points with marks type and Tree mild group

realepm: marked ppp of end points with more marks of interest I.e distance and order ...

Point pattern:
  
realpp: marked ppp, basically a superposition of realbp and realep
realppmild :ppp, basically a superposition of realbpMILD and realepmild

realpp3D:marked pp3, basically a superposition of realbp3D and realep3D
realppmild3D:marked pp3, basically a superposition of realbp3Dmild and realep3Dmild

dATA:
footlistmild: list that contains the data for mild group
footlist: list that contains the data for normal patients
SID : SUBJECT ID
SIDM: SUBJECT ID MILD"
cat(A, sep="\n")


