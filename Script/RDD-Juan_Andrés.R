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
  packs <- c("tidyverse","doBy","gdata","ggforce","haven","Hmisc","lubridate","rdd","readxl","sandwich","stargazer","dagitty","rdrobust","lmtest")
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
  
  ggplot(data = dt, mapping = aes(x = bac1, y = aged, color = DUI, group = DUI)) +
    geom_point() + stat_smooth(size = 0.5, method = loess, formula = y ~ poly(x,2)) + stat_smooth(size = 0.5, method = lm)
    geom_vline(xintercept = 0.08, linetype = "longdash") +
    theme_classic()
  
  
  
  
  
  ##creating a fake dataset (N=1000, 500 at treated, 500 at control group)
  #outcome variable
  outcome <- c(rnorm(500, mean = 50, sd = 10),  rnorm(500, mean = 70, sd = 10))
  
  #running variable
  running.var <- seq(0, 1, by = .0001)
  running.var <- sample(running.var, size = 1000, replace = T)
  
  ##Put negative values for the running variable in the control group
  running.var[1:500] <- -running.var[1:500]
  
  #treatment indicator (just a binary variable indicating treated and control groups)
  treat.ind <- c(rep(0,500), rep(1,500))
  data <- data.frame(cbind(outcome, running.var))
                           
   ggplot(data, aes(running.var, outcome, color = treat.ind, )) +
     geom_point() + geom_smooth(method=loess) +
     geom_vline(xintercept=0, linetype="longdash") +
     xlab("Running variable") +
     ylab("Outcome variable") +
     theme_classic()










  