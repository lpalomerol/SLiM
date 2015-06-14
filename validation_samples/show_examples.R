#read output files traces

doPlot <- function(p1g1,p1g2,p2g1,p2g2,main,labels) {
  plot(p1g1[-1,],type="l",
       main=sprintf(main),
       
       col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
  lines(p1g2[-1,],col="blue",type="l")
  lines(p2g1[-1,],col="red",type="l")
  lines(p2g2[-1,],col="violet",type="l")
  legend('bottomright', lty = c(1,1,1,1),
         labels,
         col=c('black','blue','red','violet'))  
}


ex1 <- read.table(file="luis_output1.txt",skip=29,comment.char="")
ex2 <- read.table(file="luis_output2.txt",skip=31,comment.char="")
ex3 <- read.table(file="luis_output3.txt",skip=35,comment.char="")
ex4 <- read.table(file="luis_output4.txt",skip=39,comment.char="")
ex5 <- read.table(file="luis_output5.txt",skip=39,comment.char="")
ex6 <- read.table(file="luis_output6.txt",skip=39,comment.char="")
ex7 <- read.table(file="luis_output7.txt",skip=41,comment.char="")
head(ex1)
head(ex2)
head(ex3)
head(ex4)
head(ex5)
head(ex6)
head(ex7)
#calculate the frequency of each mutation in the populations for each example
#example1
p1g1.ex1 <- ex1[(ex1[,4]=="p1" & ex1[,6]=="1")==T,c(2,11)]
p2g2.ex1 <- ex1[(ex1[,4]=="p2" & ex1[,6]=="2")==T,c(2,11)]

#example2
p1g1.ex2 <- ex2[(ex2[,4]=="p1" & ex2[,6]=="1")==T,c(2,11)]
p2g2.ex2 <- ex2[(ex2[,4]=="p2" & ex2[,6]=="2")==T,c(2,11)]

#example3
p1g1.ex3 <- ex3[(ex3[,4]=="p1" & ex3[,6]=="1")==T,c(2,11)] # black ->  neutral p1,m1
p1g2.ex3 <- ex3[(ex3[,4]=="p1" & ex3[,6]=="2")==T,c(2,11)] # blue ->  selectiva p1,g2
p2g1.ex3 <- ex3[(ex3[,4]=="p2" & ex3[,6]=="1")==T,c(2,11)] # red -> neutral, p2,g1
p2g2.ex3 <- ex3[(ex3[,4]=="p2" & ex3[,6]=="2")==T,c(2,11)] # violet ->selectiva, p2,g2


#example4
p1g1.ex4 <- ex4[(ex4[,4]=="p1" & ex4[,6]=="1")==T,c(2,11)] # black ->  neutral p1,m1
p1g2.ex4 <- ex4[(ex4[,4]=="p1" & ex4[,6]=="2")==T,c(2,11)] # blue ->  selectiva p1,g2
p2g1.ex4 <- ex4[(ex4[,4]=="p2" & ex4[,6]=="1")==T,c(2,11)] # red -> neutral, p2,g1
p2g2.ex4 <- ex4[(ex4[,4]=="p2" & ex4[,6]=="2")==T,c(2,11)] # violet ->selectiva, p2,g2

#example5
p1g1.ex5 <- ex5[(ex5[,4]=="p1" & ex5[,6]=="1")==T,c(2,11)] # black ->  neutral p1,m1
p1g2.ex5 <- ex5[(ex5[,4]=="p1" & ex5[,6]=="2")==T,c(2,11)] # blue ->  selectiva p1,g2
p2g1.ex5 <- ex5[(ex5[,4]=="p2" & ex5[,6]=="1")==T,c(2,11)] # red -> neutral, p2,g1
p2g2.ex5 <- ex5[(ex5[,4]=="p2" & ex5[,6]=="2")==T,c(2,11)] # violet ->selectiva, p2,g2


#example6
p1g1.ex6 <- ex6[(ex6[,4]=="p1" & ex6[,6]=="1")==T,c(2,11)] # black ->  neutral p1,m1
p1g2.ex6 <- ex6[(ex6[,4]=="p1" & ex6[,6]=="2")==T,c(2,11)] # blue ->  selectiva p1,g2
p2g1.ex6 <- ex6[(ex6[,4]=="p2" & ex6[,6]=="1")==T,c(2,11)] # red -> neutral, p2,g1
p2g2.ex6 <- ex6[(ex6[,4]=="p2" & ex6[,6]=="2")==T,c(2,11)] # violet ->selectiva, p2,g2

#example7
p1g1.ex7 <- ex7[(ex7[,4]=="p1" & ex7[,6]=="1")==T,c(2,11)] # black ->  neutral p1,m1
p1g2.ex7 <- ex7[(ex7[,4]=="p1" & ex7[,6]=="2")==T,c(2,11)] # blue ->  selectiva p1,g2
p2g1.ex7 <- ex7[(ex7[,4]=="p2" & ex7[,6]=="1")==T,c(2,11)] # red -> neutral, p2,g1
p2g2.ex7 <- ex7[(ex7[,4]=="p2" & ex7[,6]=="2")==T,c(2,11)] # violet ->selectiva, p2,g2


par(mfrow=c(1,1))
plot(p1g1.ex3[-1,],type="l",
     main=sprintf("example03: 
                  black -> neutral mutations p1g1
                  blue -> selective mutation p1g2
                  red -> neutral mutation p2g1
                  violet -> selective mutation p2ge"),
                 
     col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p1g2.ex3[-1,],col="blue",type="l")
lines(p2g1.ex3[-1,],col="red",type="l")
lines(p2g2.ex3[-1,],col="violet",type="l")
legend('bottomright', lty = c(1,1,1,1),
       c('p1g1(neut)',
         'p1g2(sel)',
         'p2g1(neut)',
         'p2ge(sel)'),
       col=c('black','blue','red','violet'))



#plot the trajectories of the neutral and the selective mutations
pdf(file="Examples_01_02.pdf",height=10,width=10)
par(mfrow=c(2,1))

plot(p1g1.ex1,type="l",
     main=sprintf("example01: \nblack -> neutral mutation at p1g1\n red -> selective mutation at p2g2"),
     col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p2g2.ex1,col="red",type="l")

plot(p1g1.ex2,type="l",
     main=sprintf("example02: \nblack -> neutral mutation at p1g1\n red -> selective mutation at p2g2"),
     col="black",xlab="generations",ylab="frequency",ylim=c(0,1000),xlim=c(1,1000))
lines(p2g2.ex2,col="red",type="l")


doPlot(p1g1.ex3,
       p1g2.ex3,
       p2g1.ex3,
       p2g2.ex3,
       "Example 3",
       c('p1g1(neut)',
         'p1g2(sel)',
         'p2g1(neut)',
         'p2ge(sel)'))

doPlot(p1g1.ex4,
       p1g2.ex4,
       p2g1.ex4,
       p2g2.ex4,
       "Example 4 (hermaphrodites)",
       c('p1g1(neut) Threshold',
         'p1g2(sel) Threshold',
         'p2g1(neut) No Threshold',
         'p2ge(sel) No threshold'))

doPlot(p1g1.ex5,
       p1g2.ex5,
       p2g1.ex5,
       p2g2.ex5,
       "Example 5 (sexed)",
       c('p1g1(neut) Threshold',
         'p1g2(sel) Threshold',
         'p2g1(neut) No Threshold',
         'p2ge(sel) No threshold'))


doPlot(p1g1.ex6,
       p1g2.ex6,
       p2g1.ex6,
       p2g2.ex6,
       "Example 6 (ratio)",
       c('p1g1(neut) 1/10',
         'p1g2(sel) 1/10',
         'p2g1(neut) 1/1',
         'p2ge(sel) 1/1'))


doPlot(p1g1.ex7,
       p1g2.ex7,
       p2g1.ex7,
       p2g2.ex7,
       "Example 7 (ratio + threshold)",
       c('p1g1(neut) 1/10 Thr',
         'p1g2(sel) 1/10 Thr',
         'p2g1(neut) 1/1 NThr',
         'p2ge(sel) 1/1 NThr'))


dev.off()

