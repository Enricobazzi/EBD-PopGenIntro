#Packages needed for the practical lesson
install.packages("vcfR")
install.packages("rehh")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("XML", repos= "http://www.omegahat.net/R"
if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("topGO")
    
