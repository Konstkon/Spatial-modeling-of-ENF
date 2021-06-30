source('./SimulateENF_subject.R')  #UCC like model
source('./noc_model_simulator.R')         #NOC model
source('./Data_load.R')                  #load datas

samplesid_per_subj=readRDS('./DATAFOLDER/DATASTORED/samplesid_per_subj')

Group = "Normal"
#simulate from UCC-like model
simulation <-lapply(seq_along(realbpdf), 
                    function(x) SimulateENFdata_subject(realbpdf[[x]],Group ,params,SID[x])) 
#Get branching points and create ppp objects
FB_simul = lapply(simulation, function(x) subset(x,Branchpoint==TRUE))                 
fb_simul_2d = lapply(FB_simul, function(x) x[c('X','Y')])
fb_simul_2d_pattern = lapply(seq_along(fb_simul_2d),
                      function(x) ppp(fb_simul_2d[[x]]$X,fb_simul_2d[[x]]$Y,
                                      window=realep[[x]]$window))

#simulate from NOC-like model
simulation_noc = lapply(seq_along(realbp),
          function(x) Noc_model(realbpdf[[x]],params,direction_list[[x]],SID[x]))
#Get branching points and create ppp objects
FB_simul_noc = lapply(simulation_noc, function(x) subset(x,Branchpoint==TRUE))                 
fb_simul_2d_noc = lapply(FB_simul_noc, function(x) x[c('X','Y')])
fb_simul_2d_pattern_noc = lapply(seq_along(fb_simul_2d_noc), 
                    function(x) ppp(fb_simul_2d_noc[[x]]$X,fb_simul_2d_noc[[x]]$Y,
                                    window=realep[[x]]$window))



#Figure B1
#OVERALL K
#K_bp = K_overall_function(realbp,samplesid_per_subj,2,TRUE)
K_branchpoint =  K_overall_function(realbranch,samplesid_per_subj,2,TRUE)
K_bp_UCC = K_overall_function(fb_simul_2d_pattern,samplesid_per_subj,2,TRUE)
K_bp_noc = K_overall_function(fb_simul_2d_pattern_noc,samplesid_per_subj,2,TRUE)
radius = Kest(realbranch[[1]])$r


par(mfrow=c(1,1))
plot(radius,sqrt(envelope(rpoispp(0.00026,win = realep[[1]]$window),Kest,nsim=39)$lo/pi)-radius,'l',xlab='r',ylab='L(r)-r',lty=3,ylim=c(-25,25))
lines(radius,sqrt(envelope(rpoispp(0.00026,win = realep[[1]]$window),Kest,nsim=39)$hi/pi)-radius,'l',lty=3)
lines(radius,sqrt(K_bp_UCC$K/pi)-radius, col=1,lty=1,lwd=2)
lines(radius,sqrt(K_branchpoint$K/pi)-radius,lty=1,lwd=2,col=3)
lines(radius, zeros(length(radius),1),col=2,lty=2,lwd=1)
lines(radius, sqrt(K_bp_noc$K/pi )- radius, col=4,lty=1,lwd=2)
legend("topright", c('Branching points(NOC)','Branching points(UCC)','Branching points(observed)',TeX('L_{theo}(r)-r'),'Poisson Envelopes'),lty=c(1,1,1,2,3),col=c(4,1,3,2,1),lwd=c(2,2,2,1,1))




#Figure 9
par(mfrow=c(1,2))
K_bp =  K_overall_function(realbp,samplesid_per_subj,2,TRUE)#basepoints
radius = K_bp$r
plot(radius,sqrt(envelope(rpoispp(0.00026,win = realep[[1]]$window),Kest,nsim=39)$lo/pi)-radius,'l',xlab='r',ylab='L(r)-r',lty=3,ylim=c(-25,25))
lines(radius,sqrt(envelope(rpoispp(0.00026,win = realep[[1]]$window),Kest,nsim=39)$hi/pi)-radius,'l',lty=3)
lines(radius,sqrt(K_branchpoint$K/pi)-radius,lty=1,lwd=2,col=3)
lines(radius,sqrt(K_bp$K/pi)-radius,lty=1,lwd=2,col=1)
lines(radius, zeros(length(radius),1),col=2,lty=2,lwd=1)
legend("topright", c('Basepoints(observed)','Branchingpoints(observed)',TeX('L_{theo}(r)-r'),'Poisson Envelopes'),lty=c(1,1,2,3),col=c(1,3,1),lwd=c(2,2,1,1))

K_bp3 =  K_overall_function(realbp3D,samplesid_per_subj,3,TRUE)#basepoints
K_branch3 =  K_overall_function(realbp3D,samplesid_per_subj,3,TRUE)#basepoints
KBP <- lapply(realbp3D, K3est, ratio=TRUE,r=seq(0,100,50))
KBP_pooled =pool(as.anylist(KBP))
KBP_pooled$r
DATA <-   hyperframe(
  g = factor(rep.int(c('BASE', 'BRANCH'),
                     c(length(realep), length(realep)))),
  ppp = c(realbp3D, realbranch3D)
)

DATA$K <- with(DATA, K3est(ppp, ratio=TRUE, rmax = 80))
KTOTAL <- anylapply(split(DATA$K, DATA$g), pool)
radius = KTOTAL$BASE$r
A = envelope(rpoispp(0.00026,win = realep[[1]]$window),Kest,nsim=39)
plot(A$r, sqrt(A$lo/pi)-A$r,'l',xlab='r',ylab='L(r)-r',lty=3,ylim=c(-25,25))
lines(A$r,sqrt(A$hi/pi)-A$r,'l',lty=3)
lines(radius,(KTOTAL$BRANCH$pooltrans/(4*pi/3))^0.33-radius,lty=1,lwd=2,col=3)
lines(radius,(KTOTAL$BASE$pooltrans/(4*pi/3))^0.33-radius,lty=1,lwd=2,col=1,'l')
lines(radius, zeros(length(radius),1),col=2,lty=2,lwd=1)
legend("topright", c('Basepoints(observed)','Branchingpoints(observed)',TeX('L_{theo}(r)-r'),'Poisson Envelopes'),lty=c(1,1,2,3),col=c(1,3,2,1),lwd=c(2,2,1,1))

