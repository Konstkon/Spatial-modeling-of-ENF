#IMPORTING LIBRARIES AND CODE ####
source('./exscript.R') #load data,and libaries
source('./Helpfunctions.R') #load functions
source('./Angles_Lengths.R')
library(fitdistrplus)
require(extraDistr)
library(latex2exp)
library(ggplot2)
library(gridExtra)

#Extract relevant information####

Group ='Normal'
FB = Get_Firstbranch(Group)  
EP = Get_endbranch(Group)

Group = 'Mild'
FB_mild = Get_Firstbranch(Group)
EP_mild=Get_endbranch(Group)


###Create Figure 4#####
data <- data.frame(
  Group = c( rep("Mild", length(FB_mild$XYangle)), rep("Healthy", length(FB$XYangle)) ),
  value = c( FB_mild$XYangle,FB$XYangle )
)

p1 = ggplot(data) + 
  stat_density(aes(x=value, color=Group),
               geom="line",position="identity")+
  labs(title="First segment",x=TeX('$\\theta'), y = "Density")+ylim(c(0,0.2))


p1new <- data %>%
  ggplot( aes(x=value, fill=Group)) +
  geom_histogram( aes(y=..density..),color="#e9ecef", alpha=0.6,binwidth =1, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="",title="First segment",x=TeX('$\\theta'), y = "Frequency")+ylim(c(0,0.2))



data1 <- data.frame(
  Group = c( rep("Mild", length(EP_mild$XYangle)), rep("Healthy", length(EP$XYangle)) ),
  value = c( EP_mild$XYangle,EP$XYangle )
)

p2 = ggplot(data1) + 
  stat_density(aes(x=value, color=Group),
               geom="line",position="identity")+
  labs(title="Endpoint segments",x=TeX('$\\phi'), y = "Density")+ylim(c(0,0.2))

p2new <- data1 %>%
  ggplot( aes(x=value, fill=Group)) +
  geom_histogram( aes(y=..density..),color="#e9ecef", alpha=0.6,binwidth =1, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="",title="Endpoint segments",x=TeX('$\\phi'), y = "Frequency")+ylim(c(0,0.2))
p2
grid.arrange(p1new, p2new,nrow=1)  #Figure 4




#Figures 3 and 5####
EP1 = subset(EP, EP$Zangle>0)
EP1_m = subset(EP_mild, EP_mild$Zangle>0)
EP_m =EP_mild
shpptsep<- spherepoints(EP1$Colatidute ,EP1$XYangle)

shpptsep<- spherepoints(FB$Colatidute ,FB$XYangle)
R = sqrt(sum(shpptsep$X)^2+sum(shpptsep$Y)^2+sum(shpptsep$Z)^2)/length(shpptsep$X)
library(extRemes)
KAPPA= (R*3-R^3)/(1-R^2)
par(mfrow=c(1,3))
par(mar=c(3,3,2,2))

#FIG 5
shiftplot(EP1$Zangle,FB$Zangle,main='Z angle',xlab = 'Endpoint quantiles',ylab='First segment quantiles - Endpoint quantiles')
abline(h=0)
shiftplot(EP1_m$Zangle,FB_mild$Zangle,main='Z angle',xlab = 'Endpoint quantiles',ylab='First segment quantiles - Endpoint quantiles')
abline(h=0)
shiftplot(EP1$Zangle,EP1_m$Zangle,main='Z angle',xlab = 'Healthy quantiles',ylab='Mild quantiles - Healthy quantiles')
abline(h=0)

#Fig 3

shiftplot(EP$Length[EP$Length<60],FB$Length[FB$Length<60],main='Length',xlab = 'Endpoint quantiles',ylab='First segment quantiles - Endpoint quantiles')
abline(h=0)
shiftplot(EP_mild$Length[EP_mild$Length<60],FB_mild$Length[FB_mild$Length<60],main='Length',xlab = 'Healthy quantiles',ylab='Mild quantiles - Healthy quantiles')
abline(h=0)
shiftplot(EP$Length[EP$Length<60],EP_m$Length[EP_m$Length<60],main='Length',xlab = 'Healthy quantiles',ylab='Mild quantiles - Healthy quantiles')
abline(h=0)



#Figure 8####

ep = realepdf[[32]]
bp = realbpdf[[32]]
branchp = realbranchdf[[32]]

completed_trees =branchp$Tree
par(mfrow=c(1,1))
plot(bp$X[bp$Tree==completed_trees[1]],bp$Y[bp$Tree==completed_trees[1]],col=1, xlim=c(15,125),
     ylim=c(-150,-60),pch=19,cex=1,xlab='X',ylab='Y')

#points(branchp$X[1],branchp$Y[1],pch=4, cex=1)
numset=c(1:length(completed_trees))

for (ntree in numset){
  segments(bp$X[bp$Tree==completed_trees[ntree]], bp$Y[bp$Tree==completed_trees[ntree]] , 
           branchp$X[ntree],branchp$Y[ntree],col=8,lwd=1
  )
  segments(branchp$X[ntree],branchp$Y[ntree],
           ep$X[ep$Tree==completed_trees[ntree]],ep$Y[ep$Tree==completed_trees[ntree]],col=8,lwd=1
  )
  points(bp$X[bp$Tree==completed_trees[ntree]],bp$Y[bp$Tree==completed_trees[ntree]],col=ntree+1,pch=19,cex=2)
  points(branchp$X[ntree],branchp$Y[ntree],pch=4, cex=2, col=ntree+1)
  points(ep$X[ep$Tree==completed_trees[ntree]],ep$Y[ep$Tree==completed_trees[ntree]], col=ntree+1, pch=24, cex=2)
  points(mean(ep$X[ep$Tree==completed_trees[ntree]]),mean(ep$Y[ep$Tree==completed_trees[ntree]]), col=ntree+1, pch=3, cex=2)
  
}
legend('topleft',c('Basepoints','Branching point','Endpoints','mean(endpoints)'),pch=c(19,4,24,3))

ntree=1
