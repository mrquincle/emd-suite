#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Make sure the data is properly formatted.
# It should be separated by tabs. It should start with non-white space.
# It should have three columns, x, y and class
# The last column will be used for coloring

separator <- ','

library(ggplot2)

if (length(args) < 4) {
	stop("Usage: Rscript --vanilla visualize_input.R <cloud 1> <cloud 2> <match> <output folder>.\n", call.=FALSE)
} 

source("read.octave.R")
	
cloud1 <- args[1]
cloud2 <- args[2]
match  <- args[3]
dir_out <- args[4]

bname <- basename(match)

datapoints1 <- read.table(cloud1, sep=separator)
datapoints2 <- read.table(cloud2, sep=separator)

datamatch <- read.table(match, sep=separator)

dataframe1 <- data.frame(
	datapoints1, 1
)
names(dataframe1) <- c("x","y","index")
dataframe2 <- data.frame(
	datapoints2, 2
)
names(dataframe2) <- c("x","y","index")

dataframe <- rbind(dataframe1, dataframe2)

matchframe <- data.frame(
	datamatch
)

#print(matchframe)

pdf(NULL)

p <- ggplot(data=dataframe, aes(x=x, y=y, color=index)) + 
	geom_point() + xlab("") + ylab("") +
	theme(legend.position="none") 

bold.text <- element_text(face = "bold")

dir.create(file.path(dir_out), showWarnings = FALSE)
ggsave(file.path(dir_out, paste(bname, ".png", sep="")))
