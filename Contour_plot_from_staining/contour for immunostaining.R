# This script uses kernel density estimation to figure out the position of 
# in situ hybridization signal

# Spinal cord is first normalized:
# - Distance from midline to lateral border = 300 um
# - Distance from ventral border to dorsal border = 650 um

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Open the raw file (.lsm or .czi)
# Image > Property > Pixel width
# Use this number as the size factor
size_factor <- 0.454
# The distance after normalization from midline to edge
x_dist <- 300
# The distance after normalization from ventral-most to dorsal-most
y_dist <- 650

# # Inferring file path from the file you choose
path <- dirname(file.choose())
setwd(path)

# # Load every .txt file in the path
file_list <- list.files(path = path)
file_list <- grep(pattern = "*.txt", x = file_list,
                  ignore.case = TRUE, value = TRUE)



raw <- lapply(file_list, function(x) {
  raw_temp <- read.table(x, header = FALSE)
  
  # Prepare y coordinate for gather() (Convert to long form)
  raw_temp$y <- seq(nrow(raw_temp))
  raw_temp <- gather(raw_temp, key = "x", value = "z", -y)
  
  # Because R would name the column number with "V"
  # I am going to remove it here or our x coordinate would be
  # "V1" "V2" instead of 1 2
  raw_temp$x <- gsub(x = raw_temp$x, pattern = "^V", replacement = "")
  raw_temp$x <- as.numeric(raw_temp$x)
  
  # We assume the sample size is x_dist horizontally / y_dist vertically
  # So we scale the coordinate to achieve this effect
  norm_factor_x <- x_dist/max(raw_temp$x)
  norm_factor_y <- y_dist/max(raw_temp$y)
  raw_temp$x <- raw_temp$x * norm_factor_x
  raw_temp$y <- raw_temp$y * norm_factor_y
  
  # The input is already binarized, and kernel density requires a input
  # of the coordinates of the positive points, so I am going to leave only
  # the positive (a.k.a, not 0) points here
  raw_temp <- filter(raw_temp, z != 0)
  
  raw_temp$sample <- x
  return(raw_temp)
})

raw <- do.call("rbind", raw)

# Interactively reformat sample labeling
cat("Current labeling for your sample is:\n")
print(unique(raw$sample))
cat("It is advisable to change them into gene names or replicate number.\n")
labels <- readline(prompt = 
                     "Rename the samples, and separate with comma (,):")

# Remove space in the input
labels <- gsub(pattern = " ", replacement = "", x = labels)
labels <- unlist(strsplit(labels, ","))

# Show the processed label in console for confirmation
cat("Your labels would be:\n")
for (i in seq(length(file_list))) {
  cat(paste(labels[i], "=", file_list[i],"\n"))
}

# Confirmation. If not correct, stop running.
go <- readline(prompt = "Is this label correct? (y/n):")
while (!go %in% c("y", "n")) {
  go <- readline(prompt = "Please answer with only y/n:")
}

if (go == "n") {
  stop("Let's start over")
}


raw$sample <- factor(x = raw$sample, levels = unique(raw$sample),
                     labels = labels)

kde <- ggplot(data = raw, aes(x = x, y = y, color = sample)) +
  geom_density_2d() +
  labs(x = "Distance from midline (um)",
       y = "Distance from upper border",
       color = "Gene") +
  scale_y_reverse()
plot(kde)
