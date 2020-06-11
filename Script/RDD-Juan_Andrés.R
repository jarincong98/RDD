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
packs <- c("tidyverse","doBy","gdata","ggforce","haven","Hmisc","lubridate","rdd","readxl","sandwich","stargazer","dagitty")
sapply(packs,require,character=TRUE)

"Data"
dt <- read.csv("Data/hansen_dwi.csv")
dt <- dt %>% mutate(Date=as.Date(Date, format = "%d%b%Y"))
str(dt)

