rm(list = ls())
source("function new.R")
data.dir <- "~/Box Sync/Strauss, Steven/E7/"
file.name <- "E7_master_in.rds"
treatment.name <- expression(paste("Melatonin (", mu, "M) - Serotonin (", mu, "M)", sep = ""))
rotation <- TRUE

data <- readRDS(paste(data.dir, file.name, sep = ""))

######################################
############ Data Cleaning ###########
######################################

source("data.cleaning.R")

##############################################
############ Finish Data Cleaning ############
##############################################

length(unique(data$Genotype))
length(unique(data$Treatment))
length(unique(data$Block))
data$Treatment <- factor(paste(data$`Melatonin_(uM)`, " - ", data$`Serotonin_(uM)`, sep = ""), 
                         levels = c("0 - 0", "50 - 0", "100 - 0", "200 - 0", "0 - 50", "0 - 100", "0 - 200", "100 - 100"))

########################################
############ Generate Plots ############
########################################

source("plots.R")

###############Table X E9############
###heritablity table 
data <- E7_data

list.herit <- heritability(data= E7_data,trait,comp= NULL,formulah =as.formula(paste(trait,"~(Genotype)",sep = "")))
outputd1 <- spread(list.herit$herit[c(1,2,3)],Explant,heritability)
pvaluelist <- pvalue_herit(data=E7_data,trait,comp= NULL,formulah =as.formula(paste(trait,"~(1|Genotype)")))

outputd2 <- spread(pvaluelist,Explant,pvalue)
outputd2$Leaf[outputd2$Leaf<0.05]<-"*"
outputd2$Petiole[outputd2$Petiole<0.05]<-"*"
outputd2$Stem[outputd2$Stem<0.05]<-"*"
outputd2$Leaf[outputd2$Leaf>0.05]<-" "
outputd2$Petiole[outputd2$Petiole>0.05]<-" "
outputd2$Stem[outputd2$Stem>0.05]<-" "

tr <- unique(outputd1$Treatment)
le <- rep(1,length(tr))
pe <- rep(1,length(tr))
se <- rep(1,length(tr))
outdata <- data.frame(tr,le,pe,se)
for(i in 1:length(unique(outputd2$Treatment))){
  outdata[i,2] <- paste(outputd1[i,2],outputd2[i,2])
  outdata[i,3] <- paste(outputd1[i,3],outputd2[i,3])
  outdata[i,4] <- paste(outputd1[i,4],outputd2[i,4])
}

outputd <- outdata
trait
write.csv(outputd, file="/Users/yz/Desktop/output.csv")

############table mean###############
data.new <- average_geno(data=E7_data,y=c("Treatment","Explant"),trait=c("prop.callus","size.callus.mean","prop.shoot","number.shoot.mean","prop.root","number.root.mean"),
                         cont=NULL)
data.new$Treatment <- as.character(data.new$Treatment)

outputd1<- spread(data.new[c(1,2,3)],key=Explant,value=prop.callus)
outputd2<- spread(data.new[c(1,2,4)],key=Explant,value=size.callus.mean)
outputd3 <- spread(data.new[c(1,2,5)],key=Explant,value=prop.shoot)
outputd4 <-spread(data.new[c(1,2,6)],key=Explant,value=number.shoot.mean)
outputd5 <-spread(data.new[c(1,2,7)],key=Explant,value=prop.root)
outputd6 <-spread(data.new[c(1,2,8)],key=Explant,value=number.root.mean)
outputd <- data.frame(outputd1$Treatment,outputd1$Leaf,outputd1$Petiole,outputd1$Stem,
                      outputd2$Leaf,outputd2$Petiole,outputd2$Stem,
                      outputd3$Leaf,outputd3$Petiole,outputd3$Stem,
                      outputd4$Leaf,outputd4$Petiole,outputd4$Stem,
                      outputd5$Leaf,outputd5$Petiole,outputd5$Stem,
                      outputd6$Leaf,outputd6$Petiole,outputd6$Stem)
write.csv(outputd, file="/Users/yz/Desktop/output.csv")

trait <- "prop.callus"
trait <- "size.callus.mean"
trait <- "prop.root"
trait <- "number.root.mean"

trait <- "prop.shoot"
trait <- "number.shoot.mean"
trait <- "length.shoot.mean"
formula <- "~Genotype*Treatment*Explant"
fit <- lm(as.formula(paste(trait, formula , sep = "")), data = E7_data )
summary(aov(fit))
d <- summary(aov(fit))[[1]]$`Pr`
d[d<0.05] <- "*"
d[d>0.05] <- ""
a <- round(summary(aov(fit))[[1]]$`Sum Sq`/sum(summary(aov(fit))[[1]]$`Sum Sq`),2)
round(summary(aov(fit))[[1]]$`Sum Sq`)
b <- c("clone ","treatment","tissue  ","clone:treatment","clone:tissue","treatment:tissue ","clone:treatment:tissue","Residuals ")
c <- data.frame(b,a)
names(c) <- c("factor","SS propotion")
outputd <- data.frame(paste(a[1]*100,"%",d[1]),
                      paste(a[2]*100,"%",d[2]),
                      paste(a[3]*100,"%",d[3]),
                      paste(a[4]*100,"%",d[4]),
                      paste(a[5]*100,"%",d[5]),
                      paste(a[6]*100,"%",d[6]),
                      paste(a[7]*100,"%",d[7]),
                      paste(sum(a[4],a[5],a[6],a[7])*100,"%"))

write.csv(outputd, file="/Users/yz/Desktop/output.csv")

