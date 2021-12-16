library(data.table)

rm(list = ls())
setwd("/Users/NilsGubela/Documents/GitHub/Softwareprojekt_Nils_Gubela/code")
data <- read.csv("result.csv", header = F)
colnames(data) <- c("sequence","structure","old occupancy","new occupancy","pauses","closest distance","closest structure","closest pause","closest occupancy")
data$length <- nchar(data$sequence)

# how many structures per length
table(data$length)

data$binary_pause <- ifelse(data$pauses == "no improvement found with pausing", 0 ,1)

# pauses found per length
table(data$length, data$binary_pause)

data$occur_only_with_pause <- ifelse((data$`old occupancy` == 0 & data$binary_pause == 1), 1, 0)

# new structures found after pausing
table(data$length, data$occur_only_with_pause)

data$not_observable <- ifelse((data$`old occupancy` == 0 & data$binary_pause == 0), 1, 0)

# not observable without and with pausing
table(data$length, data$not_observable)

hist( data$`new occupancy`-data$`old occupancy`)

seqs <- unique(data$sequence)

for(seq in seqs){
  mfe = data[data$sequence == seq,][1,]
  nmfe = data[data$sequence == seq,][2,]
}

dt_data <- as.data.table(data)
dt_data[, .N, by = sequence][N>2,]

data[data$occur_only_with_pause == 1,]

