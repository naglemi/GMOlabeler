library(readr)
library(tidyverse)
library(data.table)
library(lattice)
library(lme4)
library(VCA)
library(olsrr)
library(grid)
library(GGally)


#####A FUNCTION TO AVERAGE THE GENOTYPE AND GENERATE A NEW DATA FRAME###############################

average_geno <- function(data=data,y=c("Treatment","Block"),trait=c("prop.callus","size.callus.mean"),
                         cont=c("NAA","Two_4_D"))
{
  label <- apply(unique(data[, y]), 2, as.character)
  data.new <- NULL
  for(i in 1 : nrow(label)){
    data.sub <- data[apply(data[, y] , 1, function(data){all(data == label[i, ])}),]
    data.ti <- data.frame(t(colMeans(data.sub[, trait], na.rm = TRUE)))
    label.ti <- data.frame(t(label[i, ]))
    if(!is.null(cont)) {
      continuous <- unique(data.sub[, cont])
      onerow <- data.frame(label.ti, data.ti, continuous)
    }
    else {
      onerow <- data.frame(label.ti, data.ti)
    }
    data.new <- rbind(data.new, onerow)
  }
  data.new
}

###############Generate the themes corresponding to the different traits##########################

themenew <- function(trait, capitalization = TRUE){
  if (trait == "portion_transgenic_Shoot"){
    ifelse(capitalization, paste("Proportion of explants with transgenic shoot"), paste("proportion of explants with transgenic shoot"))
  }
  else if(trait == "portion_Shoot"){
    ifelse(capitalization, paste("Proportion of explants with shoot"), paste("proportion of explants with shoot"))
  }
  if (trait == "portion_transgenic_Callus"){
    ifelse(capitalization, paste("Proportion of explants with transgenic callus"), paste("proportion of explants with transgenic callus"))
  }
  else if(trait == "portion_Callus"){
    ifelse(capitalization, paste("Proportion of explants with callus"), paste("proportion of explants with callus"))
  }
  if (trait == "portion_transgenic_Stem"){
    ifelse(capitalization, paste("Proportion of explants with transgenic stem"), paste("proportion of explants with transgenic stem"))
  }
  else if(trait == "portion_transgenic_All_regenerated_tissue"){
    ifelse(capitalization, paste("Proportion of explants with transgenic callus or shoot"), paste("proportion of explants with transgenic callus or shoot"))
  }
  if(trait == "portion_All_regenerated_tissue"){
    ifelse(capitalization, paste("Proportion of explants with callus or shoot"), paste("proportion of explants with callus or shoot"))
  }
  else if(trait == "portion_transgenic_All_tissue"){
    ifelse(capitalization, paste("Proportion of explants with transgenic tissue"), paste("proportion of explants with transgenic tissue"))
  }
}

####################### Boxplot for the optimal Treatment ########################################
box_opt_treat <- function(trait=trait,data=data,treatment.name=treatment.name,rotation=TRUE){
  if(!rotation) {
    bwplot(as.formula(paste(trait, " ~ Treatment", sep = "")), data = data,
           xlab=list(treatment.name, cex=1.6) ,ylab=list(themenew(trait), cex=1.6), horizontal =FALSE,
           scales = list(x=list(cex=1.5), y=list(cex(1.5))),
           par.strip.text=list(cex=1.5)) 
  }
  else {
    bwplot(as.formula(paste(trait, " ~ Treatment", sep = "")), data = data,
           xlab=list(treatment.name, cex=1.6), ylab=list(themenew(trait), cex=1.6), horizontal =FALSE, 
           scales=list(x=list(rot=90, cex=1.5), y=list(cex=1.5)),
           par.strip.text=list(cex=1.5))
  }
}

#####################################################################################################
# complete ###################### Lineplot for the optimal heritability #########################
#####################################################################################################
data_summary1 <- function(x) {
  if(length(x)>=2){
    m <- mean(x)
    ymin <- min(x)-0.2
    ymax <- max(x)+0.2
    return(c(y=m,ymin=ymin,ymax=ymax))}
  else{
    return(c(y=x,ymin=x,ymax=x))
  }
}

data_summary2 <- function(x) {
  if(length(x)>=2){
    m <- mean(x)
    ymin <- min(x)-0.05
    ymax <- max(x)+0.05
    return(c(y=m,ymin=ymin,ymax=ymax))}
  else{
    return(c(y=x,ymin=x,ymax=x))
  }
}

change_clone <- function(data){
  for(i in 1:length(unique(data$Genotype))){
    t <- data$Genotype
    
    data$Genotype2[t == unique(t)[i]] <- i
  }
  
  data$Genotype2 <- as.factor(data$Genotype2)
  data
  
}


lineplot_herit<- function(data_summary,data,trait,list.herit){
  
#  a <- list.herit$herit$Explant
  b <- list.herit$herit$Treatment
  c <- list.herit$herit$heritability
  try.labs <-  paste(#a, 
    b, "heritability = ", round(c, digits = 2))
  #names(try.labs) <- a
  ggplot(data, aes(x=Genotype, y=data[,trait]))+geom_point()+
    stat_summary(fun.data=data_summary,size=0.3,color="blue")+
    facet_wrap(~ Treatment)+ 
    xlab("Genotype") 
}

lineplot_herit_new<- function(data_summary,data,trait,list.herit){
  #data$tis_treat <- paste(data$Explant,data$Treatment)
  #data$tis_treat <- factor(data$tis_treat,levels=c(paste("Leaf",levels(data$Treatment)),paste("Petiole",levels(data$Treatment)),paste("Stem",levels(data$Treatment))))
  #a <- list.herit$herit$Explant
  b <- list.herit$herit$Treatment
  c <- list.herit$herit$heritability
  try.labs <- paste(b, ", H2 = ", round(c, digits = 2), sep = "")
  names(try.labs) <- paste(a,b)
  ggplot(data, aes(x=Genotype, y=data[,trait]))+geom_point()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    stat_summary(fun.data=data_summary,size=0.3,color="blue")+
    facet_wrap(~ tis_treat, labeller = as_labeller(try.labs), nrow = 3)+
    xlab("Genotype")+ylab(themenew(trait)) }
#########################################################################################################
############################# Calculation of Heritability ###############################################
#########################################################################################################
sdh <- function(fitr){
  a <- fitr$aov.tab$VC[2]
  p <- fitr$aov.tab$VC[1]
  va <- fitr$VarCov[1,1]
  ve <- fitr$VarCov[2,2]
  covae <- fitr$VarCov[1,2]
  vh <- (a/p)^2*(va/(a^2)+(va+ve+2*covae)/(p^2)-2*(covae+va)/a/p)
  sqrt(vh)
}

heritability <- function(data=data,trait="prop.callus",comp= NULL,formula =as.formula(paste("prop.callus~(Genotype)"))){
  datah<- data[c(trait,
                 #"Explant",
                 "Treatment",comp,"Genotype")]
  print(head(datah))
  datah$Genotype<- as.factor(datah$Genotype)
  #stop()
  datah$Treatment <- factor(datah$Treatment)
  datah$Treatment <- as.character(datah$Treatment)
  #stop()
  #datah <- datah[!is.na(datah[trait]),]
  datah <- na.omit(datah)
  ntr <- length(unique(datah$Treatment))
  #nti <- length(levels(as.factor(datah$Explant)))
  Treatment <- rep(0,ntr)
  #Explant <- rep(0,ntr)
  heritability <- rep(0,ntr)
  sd_heritability <- rep(0,ntr)
  herit<- data.frame(Treatment,
                     #Explant,
                     heritability,sd_heritability)
  print("Unique treatments are")
  print(unique(datah$Treatment))
  k <- 0
  for(i in 1 : ntr)
  {

    tr <- unique(datah$Treatment)[i]
    print(paste0("Calculating heritability with treatment: ", tr))
    k <- k +1
    data.ti.tr <- subset(datah, Treatment == tr)
    herit$Treatment[k] <- tr
    #herit$Explant[k] <- ti
    #stop()
    tryrm <- try(remlMM(formula, data.ti.tr, cov=TRUE), silent=TRUE) 
    if('try-error' %in% class(tryrm)== TRUE){
      herit$heritability[k]<- NA
      herit$sd_heritability[k] <- NA
    }
    else{
      fitr <- remlMM(formula, data.ti.tr, cov=TRUE)
      print(fitr)
      herit$heritability[k] <- round(fitr$aov.tab$VC[2]/fitr$aov.tab$VC[1],2)
      herit$sd_heritability[k] <-round(sdh(fitr),2)
    }

  }
  herit$Treatment <- as.factor(herit$Treatment)
  list(herit = herit, trait = trait)
}


########################################################################################################
################################# The plot of Heritability #############################################
########################################################################################################


view_herit <- function(list.herit, rotation = TRUE){
  if(!rotation){
    ggplot(list.herit$herit,aes(x=Treatment,y=heritability)) + 
      geom_point()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=14),
            axis.text.y=element_text(size=14),
            axis.title=element_text(size=16),
            legend.title=element_text(size=16),
            legend.text=element_text(size=14),
            strip.text.x = element_text(size=14)) +
      geom_errorbar(aes(ymin=heritability-1.96*sd_heritability,ymax=heritability+1.96*sd_heritability),width=0.2,colour="blue")+
      #facet_wrap(~ Explant)+
      xlab(treatment.name) + ylab(paste("Heritability:", themenew(trait, TRUE)))
  }
  else {
    ggplot(list.herit$herit,aes(x=Treatment,y=heritability)) + 
      geom_point()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=14),
            axis.text.y=element_text(size=14),
            axis.title=element_text(size=16),
            legend.title=element_text(size=16),
            legend.text=element_text(size=14),
            strip.text.x = element_text(size=14)) +
      geom_errorbar(aes(ymin=heritability-1.96*sd_heritability,ymax=heritability+1.96*sd_heritability),width=0.2,colour="blue")+
      #facet_wrap(~ Explant)+
      xlab(treatment.name) + ylab(paste("Heritability:", themenew(trait, TRUE)))
  }
    
}


order_view_herit <- function(list.herit){
  orderlist.herit <- list.herit
  orderform <- orderlist.herit$herit
  orderform  <- transform(orderform , Treatment=reorder(Treatment, -heritability)) 
  orderlist.herit$herit <- orderform
  view_herit(orderlist.herit)
}
# order_view_herit(list.herit)


########################################################################################################
############################# The calculation of P-value################################################
########################################################################################################


pvalue_herit <- function(data=data,trait="prop.callus",comp= NULL,formula =as.formula(paste("prop.callus~(1|Genotype)"))){
  datah<- data[c(trait,
                 #"Explant",
                 "Treatment",comp,"Genotype")]
  datah <- na.omit(datah)
  #datah <- datah[!is.na(datah[trait]),]
  ntr <- length(unique(as.factor(datah$Treatment)))
  #nti <- length(unique(as.factor(datah$Explant)))
  Treatment <- rep(0,ntr)
  #Explant <- rep(0,ntr)
  pvalue <- rep(0,ntr)
  pvaluetable<- data.frame(Treatment,
                           #Explant,
                           pvalue)
  
  k <- 0
  for(i in 1 : ntr)
  {
    tr <- unique(datah$Treatment)[i]

    k <- k +1
    print(k)
    #ti <- unique(datah$Explant)[j]
    data.ti.tr <- subset(datah, Treatment == tr)
    pvaluetable$Treatment[k] <- tr
    #pvaluetable$Explant[k] <- ti
    tryrm <- try(lmer(formula, data = data.ti.tr), silent=TRUE) 
    if('try-error' %in% class(tryrm)== TRUE){
      pvaluetable$pvalue[k]<- NA
    }
    else{
      fit.3 <- lmer(formula, data = data.ti.tr)
      formula <- as.formula(paste(trait,"~1",sep = ""))
      fit.4 <- lm(formula, data = data.ti.tr)
      s34 <- anova(fit.3, fit.4,test="LRT")
      pvaluetable$pvalue[k]<- round(s34$ Pr[2],3)
    }
  }
  pvaluetable
}


#########################################################################################################
############################# The plot of P-value #######################################################
#########################################################################################################

view_pvalue <- function(data=pvaluetable,trait= trait){
  ggplot(data, aes(x=Treatment, y=pvalue)) +
    geom_bar(stat="identity")+
    #facet_wrap(~ Explant)+
    geom_hline(yintercept=0.05,colour="red")+ ggtitle(paste("The histogram of the pvalue for ",themenew(trait),sep = "")) +
    xlab("Treatment") + ylab("pvalue")
}



############################################################################################################
#############scatterplot####################################################################################
############################################################################################################


#################callus vs shoot############################################################################
cor_cal_sho <- function(cal_sho,trait1,trait2){
  corr <- round(cor(as.numeric(cal_sho[,trait1]),as.numeric(cal_sho[,trait2]),use="complete.obs"),3)
  ggplot(data = cal_sho) + 
    geom_point(mapping = aes(x=cal_sho[,trait1], y =cal_sho[,trait2],
                             #color = Treatment,
                             #shape=Explant,
                             size=1.0)) +
    xlab(themenew(trait1)) + ylab(themenew(trait2))+
    annotate("text",x=-Inf,y=Inf,label=paste("Correlation =",corr),hjust=-0.2,vjust=2,
             size = 6) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          legend.title=element_text(size=16),
          legend.text=element_text(size=14))
}






















