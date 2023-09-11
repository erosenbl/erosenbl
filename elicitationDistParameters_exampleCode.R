#Script to calculate distributions from expert elicitation responses
#Deer Ecology Panel
#Elias Rosenblatt 
#With code adapted from Molly Bletz
#13 Dec 22

## load up the packages we will need
library(fitdistrplus)
library(ggplot2)
library(logitnorm)
library(stats) 
library(rmutil) 
library(tidyverse)
library(readxl)
library(formatR)
library(gridExtra)
library(RColorBrewer)

#Read in Elicitation Results
main_data <- read.csv("DeerEcology_EE_Estimates_Raw_proofed.csv")

#Rename Fields
colnames(main_data) <- c("Name", "Response",
      "lo.durationDeerDeer", "hi.durationDeerDeer", "best.durationDeerDeer", "ci.durationDeerDeer","notes.durationDeerDeer",
      "lo.dcDeer", "hi.dcDeer", "best.dcDeer", "ci.dcDeer", "notes.dcDeer",
      "lo.baiting", "hi.baiting", "best.baiting", "ci.baiting", "notes.baiting",
      "lo.proximityRural", "hi.proximityRural", "best.proximityRural", "ci.proximityRural", "notes.proximityRural",
      "lo.durationHumanRural", "hi.durationHumanRural", "best.durationHumanRural", "ci.durationHumanRural", "notes.durationHumanRural",
      "lo.proximitySuburban", "hi.proximitySuburban", "best.proximitySuburban", "ci.proximitySuburban", "notes.proximitySuburban",
      "lo.durationHumanSuburban", "hi.durationHumanSuburban", "best.durationHumanSuburban", "ci.durationHumanSuburban", "notes.durationHumanSuburban",
      "lo.proximityCaptive", "hi.proximityCaptive", "best.proximityCaptive", "ci.proximityCaptive", "notes.proximityCaptive",
      "lo.durationHumanCaptive", "hi.durationHumanCaptive", "best.durationHumanCaptive", "ci.durationHumanCaptive", "notes.durationHumanCaptive",
      "lo.proximityDeerCaptive", "hi.proximityDeerCaptive", "best.proximityDeerCaptive", "ci.proximityDeerCaptive", "notes.proximityDeerCaptive")

dat.summary <- main_data %>% 
  mutate(ci.durationDeerDeer.D = ci.durationDeerDeer/100, ci.dcDeer.D = ci.dcDeer/100, ci.baiting.D = ci.baiting/100,
         ci.proximityRural.D = ci.proximityRural/100, ci.durationHumanRural.D = ci.durationHumanRural/100, ci.proximitySuburban.D = ci.proximitySuburban/100,
         ci.durationHumanSuburban.D = ci.durationHumanSuburban/100, ci.proximityCaptive.D = ci.proximityCaptive/100, ci.durationHumanCaptive.D = ci.durationHumanCaptive/100,
         ci.proximityDeerCaptive.D = ci.proximityDeerCaptive/100) %>%  # turn CI into probability
  filter(Response==2) ## deals will any missing responses

# making sure the variables are reading as numeric
str(dat.summary)

#Split data by questions (will just focus on Question 1)
#Q1: Given that two individual deer are in proximity (within 1.5 m of each other), how long do you expect these individuals to stay in proximity on average (minutes)?

Q1 <- select(dat.summary, "Name", 
             "lo.durationDeerDeer", "hi.durationDeerDeer", "best.durationDeerDeer", "ci.durationDeerDeer.D") %>% 
  dplyr::rename(lo = lo.durationDeerDeer, hi = hi.durationDeerDeer, best = best.durationDeerDeer, CID = ci.durationDeerDeer.D) %>% 
  mutate(.,"Expert" = seq(1,nrow(.)))

#Determine how many experts participated
n.experts <- nrow(dat.summary)

#combine datasets into list (important for use of code with multiple questions)
resp.list <- list("Q1" = Q1)

#Summarize distributions for each expert. In this set of questions, all distributions are log-normal.
IndividQuantiles<-function(x) {
  
  # gets the data out of tibble format because for loop code doesn't like tibble apparently
  param_data <- as.data.frame(x)
  
  # define number of experts in group that responded (needed later in for loops)
  n.experts <- nrow(param_data)
  
  # make empty vectors
  parms <- data.frame(Name = rep(NA, n.experts),Expert = rep(NA, n.experts),LCL = rep(NA, n.experts), Med = rep(NA, n.experts), UCL = rep(NA, n.experts), mean = rep(NA, n.experts),sd = rep(NA, n.experts))
  
  #Record name of expert (not to be used in circulated materials)
  parms[,1] <- x$Name
  
  #Record expert identifier
  parms[,2] <- x$Expert
  
  #Add expert estimates by question
  parms <- left_join(parms, x, by = "Name") %>% 
    select(., -Expert.y)
  
  #Cycle through each expert's response...
  for(j in 1:n.experts){
    
    #Isolate low, best, and high estimates.
    w <- c(param_data[j,"lo"],param_data[j,"best"],param_data[j,"hi"])
    
    #Calculate lower and upper percentiles
    lower <- (1 - param_data[j,"CID"])/2
    upper <- param_data[j,"CID"] + (1 - param_data[j,"CID"])/2

    #Draw quantiles for each expert
    if(deparse(substitute(x))!="Q2") { #Q2 is a logit normal, and all other questions were lognormal
      
      logw <- log(w)
      IPDF <- qnorm(p = c(lower,0.5,upper))
      int <- as.vector(lm(IPDF~logw)$coefficients["(Intercept)"])
      slope <- as.vector(lm(IPDF~logw)$coefficients["logw"])
      draws <- rlnorm(1000, mean = (-int/slope), sd = (1/slope))
      est <- quantile(draws,c(.025,.5,.975))
      est2 <- qmedist(data = w, distr = "lnorm",
                      probs = c(lower, upper))
      parms[j,3:5] <- as.numeric(est)
      parms[j,6] <- est2$estimate[1]
      parms[j,7] <- est2$estimate[2]
      
    } #End if statement for all questions other than Q2 (all questions that used lognormal distributions.
    
    if(deparse(substitute(x))=="Q2") { #else logit distribution for Q2
      logitw <- logit(w)
      IPDF <- qnorm(p = c(lower,0.5,upper))
      int <- as.vector(lm(IPDF~logitw)$coefficients["(Intercept)"])
      slope <- as.vector(lm(IPDF~logitw)$coefficients["logitw"])
      draws <- rlogitnorm(1000, mu = (-int/slope), sigma = (1/slope))
      est <- quantile(draws,c(.025,.5,.975))
      est2 <- qmedist(data = w, distr = "logitnorm",
                      probs = c(lower, upper), start = list(mu = -int/slope,sigma = 1/slope))
      parms[j,3:5] <- as.numeric(est)
      parms[j,6] <- est2$estimate[1]
      parms[j,7] <- est2$estimate[2]
    } #End if statement for Question 2
  } #End expert loop
  parms
}

#Visualize data

##Q1
#Summarize data
Q1_prms <- IndividQuantiles(Q1)

#Draw individual expert distribution
Q1_indiv_densities <- purrr::map(1:nrow(Q1_prms),
                                 function(y) stat_function(fun = dlnorm,
                                                           args = list(meanlog = Q1_prms$mean[y],
                                                                       sdlog = Q1_prms$sd[y]),
                                                           color = brewer.pal(n = nrow(Q1), name = "Dark2")[y],
                                                           linetype = "twodash", size = 1))

#Summarize median distribution characteristics
Q1_prms %>%
  summarise(across(c(lo, best, hi,mean,sd), median)) %>%
  as.numeric() -> Q1_agg_prms

#Estimate mean and sd for distribution
Q1_est.agg <- qmedist(Q1_agg_prms[1:3], "lnorm" ,probs = c(0.025, 0.975),start = list(meanlog = Q1_agg_prms[4],sdlog = Q1_agg_prms[5]))

# Build the density curve for aggregate distribution
aggr_densityQ1 <- stat_function(fun = dlnorm,
                                args = list(meanlog = Q1_est.agg$estimate[1],
                                            sdlog = Q1_est.agg$estimate[2]), color = "black", size = 1.5)

#Summary stats for reporting (can change percentile if you would like)
Q1_agg_prms[6] <- qlnorm(c(0.1), meanlog = Q1_est.agg$estimate[1],sdlog = Q1_est.agg$estimate[2]) #Add 80%LCL

Q1_agg_prms[7] <- qlnorm(c(0.5), meanlog = Q1_est.agg$estimate[1],sdlog = Q1_est.agg$estimate[2]) #median

Q1_agg_prms[8] <- qlnorm(c(0.9), meanlog = Q1_est.agg$estimate[1],sdlog = Q1_est.agg$estimate[2]) #Add 80%UCL

#Plotting Code
#P1: Expert responses
P1Q1 <- ggplot(data = Q1)+
  geom_linerange(mapping = aes(x = Expert, ymin = lo,ymax = hi),
                 size = 1.5,color = "grey") +
  geom_point(aes(x = Expert,y = best,color = as.factor(Expert)),
             size = 4,shape = 1,stroke = 2) +
  geom_text(aes(label = paste0(CID*100,"%"),x = Expert,y = hi),
            hjust = -0.1, size = 6) +
  theme_classic() +
  scale_color_manual(values = brewer.pal(n = nrow(Q1), name = "Dark2"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.1))) +
  xlab('Expert Identifier') +
  ylab("Duration of deer-deer proximity event (minutes)") +
  scale_x_continuous(limits = c(1, n.experts),breaks = seq(1, n.experts, by = 1))+
  scale_y_continuous(limits = c(0,max(Q1$hi)+100))+
  coord_flip() +
  guides(color = "none", size = "none")


#Plot 2, show distributions of experts and aggrigated distribution
P2Q1 <- ggplot() +
  Q1_indiv_densities +
  aggr_densityQ1 +
  geom_linerange(aes(xmin = Q1_agg_prms[6], xmax = Q1_agg_prms[8], y = 0), color = "grey", size = 3) +
  geom_point(aes(x = Q1_agg_prms[7], y = 0)) +
  scale_x_continuous(limits = c(0, 50)) +
  theme_classic() +
  xlab("Duration of deer-deer proximity event (minutes)")+
  ylab("Density") +
  guides(color = "none", size = "none")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.1)))

grid.arrange(P1Q1,P2Q1, ncol = 2) #Expert responses with agg dist

ggpubr::ggarrange(ggpubr::ggarrange(P1Q1,P2Q1, ncol = 2)+ #Plot expert responses with agg curve
                    annotate("label", x = 0.1, y = .9, label = "A", fill = "white", size = 6) +
                    annotate("label", x = 0.9, y = .9, label = "B", fill = "white", size = 6))
