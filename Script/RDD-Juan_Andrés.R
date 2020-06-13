"Causal Inference, 2020-1
Juan Andrés Rincón Galvis
juana.rincon@urosario.edu.co
ja.rincong@uniandes.edu.co"

"Assignment 4: Regression Discontinuity Design"


"Preliminary Organization"

  #Cleaning Working space
  rm(list = ls())
  cat("\f")
  
  #Setting Working Directory
  Juancho <- "/Users/ja.rincong/Git/RDD"
  setwd(Juancho)
  
  #Packages
  packs <- c("tidyverse","doBy","gdata","ggforce","ggpubr","haven","Hmisc","lubridate","rdd","readxl",
             "sandwich","stargazer","dagitty","rdrobust","lmtest")
  sapply(packs,require,character=TRUE)

"Data"
  dt <- read.csv("Data/hansen_dwi.csv")
  dt <- dt %>% mutate(Date=as.Date(Date, format = "%d%b%Y"))
  str(dt)


"------------------------------Data Manipulation------------------------------"
  #Treatment variable
  dt <- dt %>% mutate(DUI=ifelse(bac1>=0.08,1,0))
  #Descriptive Statistics of treatment var
  stargazer(dt,
            title = "Treatment variable descriptive statistics.", 
            summary = T, 
            out = "Tables/TreatmentStats.tex",
            keep = "DUI")

"------------------------McCrary Test for manipulation------------------------"
  #Test running
  McCraryTest <- rdd::DCdensity(runvar = dt$bac1,cutpoint = 0.08,
                                ext.out = T, 
                                htest = T, 
                                plot = T,
                                ) #With bin-width = standard
  #Results Parameters
  MC_results1 <- McCraryTest[1:7]
  MC_results1
  #Saving the data to graph with ggplot2 (prettier)
  MC <- McCraryTest$data
  MC <- MC %>% mutate(cellval=ifelse(cellval==0,NA,cellval))
  #Pretty ggplot2 graph
  ggplot(data = MC, mapping = aes(x = cellmp, y = cellval)) +
    geom_point(alpha=.5)+
    geom_col()+
    labs(title = "McCrary Density for Running Variable", x = "Running Variable: Blood Alcohol Content", y = "Density")+
    geom_smooth(se = F, size = .5)+
    geom_vline(xintercept = 0.08, linetype = "dashed")+
    theme_classic()
  ggsave("Figures/McCraryDensity.pdf", height = 5, width = 9)

  #Test Running (bin-width=0.001 as in the paper)
  McCraryTest2 <- rdd::DCdensity(runvar = dt$bac1,cutpoint = 0.08,
                                ext.out = T, 
                                htest = T, 
                                plot = T,
                                bin = 0.001
  ) #With bin-width = 0.001 Ho is not rejected
  MC_results2 <- McCraryTest2[1:7]
  #Both results
  MC_results <- list(MC_results1,MC_results2)
  
  #Remove unnecesary objects
  rm(MC_results1,MC_results2)
"------------------------------Covariate Balance------------------------------"
  
  # Model (1): Y = DUI + bac1 + bac1*DUI + X

  #Bandwidth
  h <- 0.05
  #Centered Running Variable, Rectangular Kernel and Triangular Kernel
  dt <- dt %>% mutate(bac1_cent=bac1-0.08) %>% mutate(kernel_r=ifelse(abs(bac1_cent)<h,1/(2*h),0)) %>% mutate(kernel_t=ifelse(abs(bac1_cent)<h,1-(abs(bac1_cent)/h),0))
  
  #Estimating local linear regression discontinuity design for 4 variables
    #Rectangular Kernel
    #Bandwidth of 0.05
  dep_vars <- c("white","male","aged","acc")
  results <- list()
  for (i in 1:4) {
    print(dep_vars[i])
    mdl <- lm(formula = dt[,dep_vars[i]] ~ dt[,"DUI"] + dt[,"bac1_cent"] + dt[,"bac1_cent"]*dt[,"DUI"], weights = dt[,"kernel_r"])
    mdl_r <- coeftest(mdl, vcov = vcovHC(mdl,"HC1"))
    results[[i]] <- mdl_r
  }
  
  #Output Table
  stargazer(results[1],results[2],results[3],results[4],
            title = "Regression Discontinuity Estimates for the effect of exceeding BAC threshold o predetermined characteristics.",
            out = "Tables/RDD_out.tex",
            out.header = T,
            header = F,
            column.labels = c("White","Male","Age","Accident"),
            covariate.labels = c("DUI","BAC","BAC$*$DUI"),
            keep = "DUI", nobs = T, mean.sd = T, label = "RDDout"
            )
  #Mean at 0.079
    print(mean(dt$white[dt$bac1<=0.079]))  
    print(mean(dt$male[dt$bac1<=0.079]))  
    print(mean(dt$aged[dt$bac1<=0.079]))  
    print(mean(dt$acc[dt$bac1<=0.079]))  
  
    
  plot(RDestimate(formula = male ~ bac1, data = dt, cutpoint = 0.08, bw = 0.05, kernel = "rectangular", se.type = "HC1"))
  abline(v=0.08, lty = 2)
  
"------------------------------Covariate Balance Graphs------------------------------"

  #Pretty Plots
  bin <- cut(dt$bac1, seq(min(dt$bac1),max(dt$bac1),0.002)) #bins = 0.002 as in the paper
  data_bin <- aggregate(dt,list(bin),function(x) { return(c(mean(x),length(x)))})
  
  #Linear Fit
    #Accident
    ggacc <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=acc[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = acc, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.05,0.25)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Fit","Treatment Fit"))+
      labs(title = "Panel A. Accident at scene")+
      theme_classic()+
      theme(legend.position = "bottom")
    #Male
    ggmale <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=male[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = male, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.74,0.82)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Fit","Treatment Fit"))+
      labs(title = "Panel B. Male")+
      theme_classic()+
      theme(legend.position = "bottom")
    #Age
    ggage <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=aged[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = aged, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(33,39)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Fit","Treatment Fit"))+
      labs(title = "Panel C. Age")+
      theme_classic()+
      theme(legend.position = "bottom")
    #White
    ggwhite <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=white[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = white, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.8,0.9)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Fit","Treatment Fit"))+
      labs(title = "Panel D. White")+
      theme_classic()+
      theme(legend.position = "bottom")
  #Linear fit group
    ggarrange(ggacc,ggmale,ggage,ggwhite,ncol = 2, nrow = 2)
    ggsave("Figures/CovBalanceLinear.pdf", height = 6, width = 9)
    
  #Quadratic Fit
    #Accident
    ggacc <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=acc[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = acc, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm, formula = y ~ poly(x,2))+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.05,0.25)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control QFit","Treatment QFit"))+
      labs(title = "Panel A. Accident at scene")+
      theme_classic()+
      theme(legend.position = "bottom")
    #Male
    ggmale <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=male[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = male, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm, formula = y ~ poly(x,2))+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.74,0.82)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control QFit","Treatment QFit"))+
      labs(title = "Panel B. Male")+
      theme_classic()+
      theme(legend.position = "bottom")
    #Age
    ggage <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=aged[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = aged, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm, formula = y ~ poly(x,2))+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(33,39)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control QFit","Treatment QFit"))+
      labs(title = "Panel C. Age")+
      theme_classic()+
      theme(legend.position = "bottom")
    #White
    ggwhite <- ggplot() + 
      geom_point(data = data_bin, aes(x=bac1[,1],y=white[,1], alpha=.5)) + 
      stat_smooth(data = dt, aes(x = bac1, y = white, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm, formula = y ~ poly(x,2))+
      geom_vline(xintercept = 0.08, linetype = "longdash")+
      scale_x_continuous(name = "BAC", limits = c(0,0.2))+
      scale_y_continuous(name = "")+
      coord_cartesian(ylim=c(0.8,0.9)) +
      scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
      scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control QFit","Treatment QFit"))+
      labs(title = "Panel D. White")+
      theme_classic()+
      theme(legend.position = "bottom")
    #Linear fit group
    ggarrange(ggacc,ggmale,ggage,ggwhite,ncol = 2, nrow = 2)
    ggsave("Figures/CovBalanceQuad.pdf", height = 6, width = 9)
    
"------------------------------RDD Estimation on Recidivism------------------------------"

# Model (1): Y = DUI + bac1 + bac1*DUI + X
dt <- dt %>% mutate(bac1_centDUI2=(bac1_cent*DUI)^2)
#Models
  mdl1 <- recidivism ~ DUI + bac1_cent + white + male + aged + acc
  mdl2 <- recidivism ~ DUI + bac1_cent + bac1_cent*DUI + white + male + aged + acc
  mdl3 <- recidivism ~ DUI + bac1_cent + bac1_cent*DUI + bac1_centDUI2 + white + male + aged + acc

#Panel A: Bandwidth h = 0.05
  mdl_a1 <- lm(formula = mdl1, data = dt, weights = kernel_r)
  mdl_a1_r <- coeftest(mdl_a1, vcov = vcovHC(mdl_a1,"HC1"))
  
  mdl_a2 <- lm(formula = mdl2, data = dt, weights = kernel_r)
  mdl_a2_r <- coeftest(mdl_a2, vcov = vcovHC(mdl_a2,"HC1"))

  mdl_a3 <- lm(formula = mdl3, data = dt, weights = kernel_r)
  mdl_a3_r <- coeftest(mdl_a3, vcov = vcovHC(mdl_a3,"HC1"))
  
  #Generate table A
  stargazer(mdl_a1_r,mdl_a2_r,mdl_a3_r, title = "", 
            out = "Tables/RDresultsA.tex", 
            out.header = T,
            header = F,
            column.labels = c("Linear Control","Interaction","Squared Interaction"),
            covariate.labels = "DUI",
            keep = "DUI", nobs = T, mean.sd = T, label = "RDresultsA")
  
#Panel B: Bandwidth h = 0.025
  #New Rectangular Kernel
  h = 0.025
  dt <- dt %>% mutate(kernel_r2=ifelse(abs(bac1_cent)<h,1/(2*h),0))
  
  mdl_b1 <- lm(formula = mdl1, data = dt, weights = kernel_r2)
  mdl_b1_r <- coeftest(mdl_b1, vcov = vcovHC(mdl_b1,"HC1"))
  
  mdl_b2 <- lm(formula = mdl2, data = dt, weights = kernel_r2)
  mdl_b2_r <- coeftest(mdl_b2, vcov = vcovHC(mdl_b2,"HC1"))
  
  mdl_b3 <- lm(formula = mdl3, data = dt, weights = kernel_r2)
  mdl_b3_r <- coeftest(mdl_b3, vcov = vcovHC(mdl_b3,"HC1"))  
  
  #Generate table B
  stargazer(mdl_b1_r,mdl_b2_r,mdl_b3_r, title = "", 
            out = "Tables/RDresultsB.tex", 
            out.header = T,
            header = F,
            column.labels = c("Linear Control","Interaction","Squared Interaction"),
            covariate.labels = "DUI",
            keep = "DUI", nobs = T, mean.sd = T, label = "RDresultsB")
  
  "------------------------------RDD Estimation on Recidivism Graphs------------------------------"
  
  #Filter Data on BAC<0.15
  dt <- dt %>% filter(bac1<0.15)
  #Update bins
  bin <- cut(dt$bac1, seq(min(dt$bac1),max(dt$bac1),0.002)) #bins = 0.002 as in the paper
  data_bin <- aggregate(dt,list(bin),function(x) { return(c(mean(x),length(x)))})
  
  #Graphs: BAC on Recidivism 
  #Linear Fit
  ggrecidl <- ggplot() + 
    geom_point(data = data_bin, aes(x=bac1[,1],y=recidivism[,1], alpha=.5)) + 
    stat_smooth(data = dt, aes(x = bac1, y = recidivism, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
    geom_vline(xintercept = 0.08, linetype = "longdash")+
    scale_x_continuous(name = "BAC")+
    scale_y_continuous(name = "")+
    scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
    scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control L. Fit","Treatment L. Fit"))+
    labs(title = "BAC on Recidivism", subtitle = "Linear Fit")+
    theme_classic()+
    theme(legend.position = "bottom")
  ggrecidl  
  
  #Quadratic Fit
  ggrecidq <- ggplot() + 
    geom_point(data = data_bin, aes(x=bac1[,1],y=recidivism[,1], alpha=.5)) + 
    stat_smooth(data = dt, aes(x = bac1, y = recidivism, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm, formula = y ~ poly(x,2))+
    geom_vline(xintercept = 0.08, linetype = "longdash")+
    scale_x_continuous(name = "BAC")+
    scale_y_continuous(name = "")+
    scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
    scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Q. Fit","Treatment Q. Fit"))+
    labs(title = "BAC on Recidivism", subtitle = "Quadratic Fit")+
    theme_classic()+
    theme(legend.position = "bottom")
  ggrecidq 
  
  #Arrange the plots
  ggarrange(ggrecidl,ggrecidq,ncol = 2, nrow = 1)
  ggsave("Figures/RDDrecid.pdf", height = 4, width = 9)
    