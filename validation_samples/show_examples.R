#read output files traces

ex1 <- read.table(file="output_example__newfeat_01.txt",skip=29,comment.char="")
ex2 <- read.table(file="output_example__newfeat_02.txt",skip=31,comment.char="")
ex3 <- read.table(file="output_example__newfeat_03.txt",skip=35,comment.char="")
head(ex1)
head(ex2)
head(ex3)

#calculate the frequency of each mutation in the populations for each example
#example1
p1g1.ex1 <- ex1[(ex1[,4]=="p1" & ex1[,6]=="1")==T,c(2,11)]
p2g2.ex1 <- ex1[(ex1[,4]=="p2" & ex1[,6]=="2")==T,c(2,11)]

#example2
p1g1.ex2 <- ex2[(ex2[,4]=="p1" & ex2[,6]=="1")==T,c(2,11)]
p2g2.ex2 <- ex2[(ex2[,4]=="p2" & ex2[,6]=="2")==T,c(2,11)]

#example3
p1g1.ex3 <- ex3[(ex3[,4]=="p1" & ex3[,6]=="1")==T,c(2,11)]
p1g2.ex3 <- ex3[(ex3[,4]=="p1" & ex3[,6]=="2")==T,c(2,11)]
p2g1.ex3 <- ex3[(ex3[,4]=="p2" & ex3[,6]=="1")==T,c(2,11)]
p2g2.ex3 <- ex3[(ex3[,4]=="p2" & ex3[,6]=="2")==T,c(2,11)]

#plot the trajectories of the neutral and the selective mutations
pdf(file="Examples_01_02.pdf",height=10,width=10)
par(mfrow=c(2,1))

plot(p1g1.ex1,type="l",main=sprintf("example01: \nblack -> neutral mutation at p1g1\n red -> selective mutation at p2g2"),col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p2g2.ex1,col="red",type="l")

plot(p1g1.ex2,type="l",main=sprintf("example02: \nblack -> neutral mutation at p1g1\n red -> selective mutation at p2g2"),col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p2g2.ex2,col="red",type="l")

par(mfrow=c(1,1))
plot(p1g1.ex3[-1,],type="l",main=sprintf("example03: \nblack,blue -> neutral mutations p1g1 and p1g2\n red,violet -> selective mutations p2g1 and p2g2"),col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p1g2.ex3[-1,],col="blue",type="l")
lines(p2g1.ex3[-1,],col="red",type="l")
lines(p2g2.ex3[-1,],col="violet",type="l")

dev.off()

