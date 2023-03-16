
# ---
# title: "Soil USDA Classification"
# author: "Moazzam Rind"
# date: "3/16/2023"
# output: "Soil_classify.txt"
# ---


#setwd (........)                   #Set working directory

# soil Texture triangle.

#Install Package

#install.packages( pkgs = "soiltexture" )

#Load Library

library(soiltexture)

#Empty Triangle
TT.plot( class.sys = "none" )

#USDA Soil Triangle (Page 13)
TT.plot( class.sys = "USDA.TT" )


#For Classification Refer Page 68

#install Package to import raster

#install.packages("raster")
#install.packages ("rgdal")
#install.packages("writexl")

library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(writexl)

#Load Clay raster (change location as needed )
Clay <- raster("Percent Clay_0 to 30cm.tif")


#verify if the data is uploaded
Clay

a <- data.frame(levels(Clay))

data_clay <- data.frame(a$ID,a$CLAY_DCP)

#Renames column names
colnames(data_clay)<- c("ID","CLAY")


#Load Silt raster (change location as needed for you)
Silt <- raster("Percent Silt_0 to 30cm.tif")

#verify if the data is uploaded
Silt

b <- data.frame(levels(Silt))

data_silt <- data.frame(b$ID,b$SILT_DCP)

#Renames column names
colnames(data_silt)<- c("ID","SILT")


#Load Sand raster (change location as needed for you)
Sand <- raster("Percent Sand_0 to 30cm.tif")

#verify if the data is uploaded
Sand

c <- data.frame(levels(Sand))

data_sand <- data.frame(c$ID,c$SAND_DCP)

#Renames column names
colnames(data_sand)<- c("ID","SAND")

#Merge three datasets into one

Soil_Data1 <- merge(data_clay,data_sand, by=c("ID"))

Soil_Data <- merge(Soil_Data1,data_silt, by=c("ID"))

#Rearrange order of the columns

Soil_Data <- Soil_Data [,c(1,2,4,3)]

write_xlsx(Soil_Data,"Soil_Data.xlsx")     #writing Soil Data as xlsx. Just for check



#classification of data

classify <- TT.points.in.classes(
  tri.data = Soil_Data [,c(2,3,4)],
  tri.sum.tst = FALSE,
  class.sys = "USDA.TT"
  )


#Export the file as text.file

write.table(classify, file = "Soil_classify.txt", sep = "\t",
             col.names = TRUE)











