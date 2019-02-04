require("here")
require("keras")
require("textreuse")
require("dplyr")
require("keras")
require("EBImage")
require("stringr")
require("pbapply")
require("OpenImageR")
require("jpeg")

diag <- read.csv("MRI_Diagnosis.csv", header = T, na.strings = "")
imgpath <- list.files("~/Documents/TDI/Challenge/jpg/")
imgpath <- c(imgpath)
imgnames <- filenames(imgpath)
imgsub <- factor(unlist(lapply(strsplit(imgnames, split = "[_]"), function(x){x[1]}))) # 535 unique subjects
imgdiag <- diag[diag$Subject%in%imgsub, c("Subject", "Gender", "NORMCOG", "DEMENTED")]
imgdiag[is.na(imgdiag)] <- 0
imgdiag <- imgdiag[-which(imgdiag$NORMCOG==imgdiag$DEMENTED),] # Remove rows where diagnosis is the same for normal and demented. This could be due to missing data or mislabelling
imgdiag <- imgdiag[imgsub%in%imgdiag$Subject,]
imgdiag <- imgdiag[!duplicated(imgdiag$Subject),]
pdata <- data.frame("Path" = imgpath, "Subject" = imgsub)
pdata <- left_join(pdata, imgdiag, by="Subject")
pdata <- pdata[complete.cases(pdata),]


train <- sample_n(pdata, 0.75*nrow(pdata))
test <- setdiff(pdata, train)
train$Path <- as.character(gsub(train$Path, pattern = "OAS", replacement = "./jpg/OAS"))
test$Path <- as.character(gsub(test$Path, pattern = "OAS", replacement = "./jpg/OAS"))
train_img <- sapply(train$Path, imager::load.image)

width <- 50
height <- 50



