#!/usr/bin/Rscript

# Workaround to allow calling from R (interactively) as well as from command line
if(!exists("vargs")) {
	vargs = commandArgs(trailingOnly=TRUE)
} 

# Make sure the data is properly formatted.
# It should be separated by tabs. It should start with non-white space.
# It should have three columns, x, y and class
# The last column will be used for coloring

# Regretfully, ggplot cannot be used to visualize 3D scatter plots
separator <- '\t'

library(plot3D)
library(rgl)
library(plot3Drgl)

if (length(vargs) < 2) {
	stop("\n\nUsage:\t./display_cube.R <input file> <output folder>\n\nor\n\tR> vargs=c(<input file>,<output folder>)\n\tR> source('display_cube.R')", call.=FALSE)
} 

source("read.octave.R")
	
file_in <- vargs[1]
dir_out <- vargs[2]

bname <- basename(file_in)

datapoints <- read.table(file_in, sep=separator)

dataframe <- data.frame(
	datapoints
)

#filename <- paste(dir_out, "/", "original_", bname, ".png", sep="")
filename <- paste(dir_out, "/", "original_", bname, ".pdf", sep="")
dir.create(file.path(dir_out), showWarnings = FALSE)
pdf(filename)
#png(filename, width = 400, height = 400, units = "px", res = 300)
#png(filename, width = 2048, height = 2048)
#png(filename, width = 1024, height = 1024)
#png(filename, width = 512, height = 512)

x <- dataframe$V1
y <- dataframe$V2
z <- dataframe$V3
color <- dataframe$V4

# byt="g" gives gray background with white grid lines
# pch is point shape, cex is point size
# colkey = FALSE removes legend
# theta and phi are viewing angles (default is 40 for both, theta rotates from right to left, phi from bottom to top)
# with theta = 60 we look more at the right side, with phi = 25 we look less from the top
# main = "Title"
#scatter3D(x, y, z, colvar=color, colkey = FALSE, theta = 60, phi = 25, bty = "g", pch = 20, cex = 2)
scatter3D(x, y, z, colvar=color, colkey = FALSE, theta = 60, phi = 0, bty = "g", pch = 20, cex = 1)

#addlines = FALSE, length = 20, width = 20)

invisible(dev.off())

cat("Written to \"", filename, "\"\n", sep="")
