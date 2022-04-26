library(stringr)
library(progress)
library(mclust)

##### 1. Download Station Data #####
# Based on station list from R
# Download data files into new folder
# https://statistics.berkeley.edu/computing/faqs/reading-web-pages-r


# Specify station
stations = read.csv("stations20.txt", header=F)

download_data <- function(s, obsCode){
  url = paste0('https://reg.bom.gov.au/jsp/ncc/cdio/weatherData/av?p_nccObsCode=',
               obsCode,'&p_display_type=dataFile&p_startYear=&p_c=&p_stn_num=',s)
  
  x = readLines(url)
  line = grep("All years of data", x)
  
  pattern ='p_c=-[0-9]+'
  phrase = str_match(x[line], pattern)
  
  dat_url = paste0('https://reg.bom.gov.au/jsp/ncc/cdio/weatherData/av?p_display_type=monthlyZippedDataFile&p_stn_num=', s,
                   '&', phrase, '&p_nccObsCode=',obsCode, '&p_startYear=')
  
  folder = paste0('~/Documents/UNSW/Math5271/Portfolio/Math5271 Portfolio Project/Station Data/',obsCode)
  
  dest = paste0(folder, '/', station, '.zip')
  
  download.file(dat_url, dest, quiet = TRUE)
}

rainCode = 139
maxCode = 36
minCode = 38

pb <- progress_bar$new(format = "[:bar] :percent || Estimated time remaining: :eta",
                       total = length(stations),
                       complete = "-",   # Completion bar character
                       incomplete = " ", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar

# Download for all stations
for (station in stations){
  pb$tick()
  download_data(station, rainCode)
  download_data(station, maxCode)
  download_data(station, minCode)
}
##### End #####

##### 2. Exploratory Analysis #####
# Ensure you're in the directory with station data 

# Specify stations
stations = read.csv("stations20.txt", header=F)

# Check that Data exist between 2000 and 2020 for each station & type
ymin = 2000
ymax = ymin + 20
rem <- c()
for (i in 1:length(stations)){
  for (type in c(1,2,4)){
    dat = read.dat1(type, stations[1,i])
    # remove stations which don't have data for the specified period
    if(min(dat$Year) >= ymin || max(dat$Year) <= ymax){
      # add to list of stations to remove
      rem <- c(rem, stations[1,i])
      # don't check other types if one is missing
      break
    } 
    # Include below to remove poor rainfall QC stations
  }
}

stations.sub = stations[!(stations %in% rem)]
length(stations.sub)


# Exploratory Analysis of Variables
# Plot over time
create_time_plots(1, 'Rain')
create_time_plots(2, "Max Temp")
create_time_plots(4, "Min Temp")
create_averaged_time_plots(1, "Average Rain")
create_averaged_time_plots(2, "Average Max")
create_averaged_time_plots(3, "Average Min")

# Exploratory Analysis
# Create dataframe of Mean Temperature and Mean Rainfall for each station
par(mfrow = c(1,1))
pdf(file = "Mean.pdf")
mean.df <- as.data.frame(t(sapply(stations.sub, mean.read.station)))
colnames(mean.df) <- c("Rainfall", "Temperature")
rownames(mean.df) <- stations.sub
saveRDS(mean.df, file = "M.rds")

hist(mean.df$Temp,
     main = "Temp",
     xlab = "Amount",
     ylab = "Frequency")
hist(mean.df$Rain,
     main = "Rain",
     xlab = "Amount",
     ylab = "Frequency")
dev.off()

seasonal.mean.df <- as.data.frame(read.all.stations(stations.sub, seasonal.means.station))
rownames(seasonal.mean.df) <- stations.sub
# Save data
saveRDS(seasonal.mean.df, file = "SDF.rds")

pdf(file = "Seasonal Mean Hist.pdf")
for(i in 1:24){
  hist(seasonal.mean.df[,i],
       main = colnames(seasonal.mean.df)[i],
       xlab = "")
}
dev.off()

pdf(file = "Seasonal Mean Stations.pdf")
par(mfrow = c(2,2))
for(i in 1:length(stations.sub)){
  plot(1:12, c(seasonal.mean.df[i, 1:12]),
       main = paste0("Station ", stations.sub[i], " Rainfall"),
       type = "l",
       xlab = "Month",
       ylab = "log(Precipitation)")
  plot(1:12, c(seasonal.mean.df[i, 13:24]),
       main = paste0("Station ", stations.sub[i], " Temp"),
       type = "l",
       xlab = "Month",
       ylab = "Temp")
}
dev.off()

SD <- as.data.frame(read.all.stations(stations.sub,
                                      seasonal.means.diff.station))
saveRDS(SD, "SDDF.rds")

# Example figure for distributions of temperature and rainfall
dat <- read.all.stations(stations.sub, read.station.plain)
colnames(dat) <- c("Rainfall", "Temperature")

# Create histogram for report
par(mfrow = c(1,2))
hist(dat$Rainfall, xlab = "Rainfall (mm)", 
     main = "Histogram of Rainfall")
hist(dat$Temperature, xlab = expression("Temperature ("*~degree*C*")"),
     main = "Histogram of Temperature")

##### End #####

##### 3. Clustering #####
##### Setup #####

# Read in Data
S <- readRDS("SDF.rds")
SD <- readRDS("SDDF.rds")
M <- readRDS("M.rds")
stations.sub = readRDS("stations_sub.rds")

# Split sample into training & test set
set.seed(1234) # for reproducibility
sample <- sample(c(rep(1, 0.7*dim(M)[1]), rep(0, 0.3*dim(M)[1])))

S.train <- S[sample == 1,]
S.test <- S[sample == 0,]

M.train <- M[sample == 1,]
M.test <- M[sample == 0,]

SD.train <- SD[sample == 1,]
SD.test <- SD[sample == 0,]
##### End #####

##### FMM New #####
# create initialisations based on PCR w/ set seed for replication
hcm <- hc(M.train, use = "PCR", seed = 5271)
hcs <- hc(S.train, use = "PCR", seed = 5271)
hcsd <- hc(SD.train, use = "PCR", seed = 5271)

# find optimal model & number of clusters
# limit itmax to reduce computation time
(m.icl <- mclustICL(M.train, G=2:11, control = emControl(itmax = 1000),
                    initialization = list(hcPairs = hcm))) #VVE 2
(s.icl <- mclustICL(S.train, G=2:11, control = emControl(itmax = 1000),
                    initialization = list(hcPairs = hcs))) # VVE 5
(sd.icl <- mclustICL(SD.train, G=2:11, control = emControl(itmax = 1000),
                     initialization = list(hcPairs = hcsd))) # VVE 6

# plot ICL
par(mfrow = c(2,2))
plot(m.icl)
plot(s.icl)
plot(sd.icl)

# Fit the optimal models
fmm.M <- Mclust(M.train, G = 2, model = "EVV")
fmm.S <- Mclust(S.train, G = 6, model = "VVE")
fmm.SD <- Mclust(SD.train, G = 6, model = "VVE")

# Predict the test set
pred.fmm.M <- predict.Mclust(fmm.M, M.test)
pred.fmm.S <- predict.Mclust(fmm.S, S.test)
pred.fmm.SD <- predict.Mclust(fmm.SD, SD.test)

# Save results for reproducibility
saveRDS(fmm.M, file = "fmmM.rds")
saveRDS(fmm.S, file = "fmmS.rds")
saveRDS(fmm.SD, file = "fmmSD.rds")

saveRDS(pred.fmm.M, file = "fmmMp.rds")
saveRDS(pred.fmm.S, file = "fmmSp.rds")
saveRDS(pred.fmm.SD, file = "fmmSDp.rds")
##### End #####

##### Plotting #####

# Plot clusters individually
plot.fmm.test("S clusters.pdf", S.test, fmm.S, pred.fmm.S)
plot.fmm.test("SD clusters.pdf", SD.test, fmm.SD, pred.fmm.SD)

# Plot clusters on one plot for comparison
plot.fmm.one("Sall.pdf", S.train, fmm.S, pred.fmm.S, S.test)
plot.fmm.one("SDall.pdf", SD.train, fmm.SD, pred.fmm.SD, SD.test)
##### End #####

##### Determine CZ #####
summer <- c(1,2,12)
winter <- c(6,7,8)

S.cluster.summary <- data.frame(row.names = c("Summer Rain",
                                              "Winter Rain",
                                              "Summer Temp",
                                              "Winter Temp"))
for (i in 1:fmm.S$G){
  summer.rain <- exp(mean(fmm.S$parameters$mean[summer, i]))
  summer.temp <- mean(fmm.S$parameters$mean[summer + 12, i])
  winter.rain <- exp(mean(fmm.S$parameters$mean[winter, i]))
  winter.temp <- mean(fmm.S$parameters$mean[winter + 12, i])
  
  
  S.cluster.summary <- cbind(S.cluster.summary,
                             c(summer.rain,winter.rain, summer.temp,  winter.temp)) 
}
colnames(S.cluster.summary) <- seq(1, fmm.S$G)

SD.cluster.summary <- data.frame(row.names = c("Summer Rain",
                                               "Winter Rain",
                                               "Summer Temp",
                                               "Winter Temp",
                                               "Summer Diff",
                                               "Winter Diff"))
for (i in 1:fmm.SD$G){
  summer.rain <- exp(mean(fmm.SD$parameters$mean[summer, i]))
  summer.temp <- mean(fmm.SD$parameters$mean[summer + 12, i])
  summer.diff <- mean(fmm.SD$parameters$mean[summer + 24, i])
  winter.rain <- exp(mean(fmm.SD$parameters$mean[winter, i]))
  winter.temp <- mean(fmm.SD$parameters$mean[winter + 12, i])
  winter.diff <- mean(fmm.SD$parameters$mean[winter + 24, i])
  
  
  SD.cluster.summary <- cbind(SD.cluster.summary,
                              c(summer.rain, winter.rain,summer.temp, 
                                winter.temp, summer.diff, winter.diff)) 
}

colnames(SD.cluster.summary) <- seq(1, fmm.SD$G)
##### End #####

##### Log Likelihood #####
ll <- list(fmm.M$loglik, fmm.S$loglik, fmm.SD$loglik)
which.max(ll)
##### End #####

##### Output results for plotting in Python #####
res.train <- cbind(fmm.M$classification,
                   fmm.S$classification,
                   fmm.SD$classification,
                   rep(1, length(fmm.M$classification)))

res.test <- cbind(pred.fmm.M$classification,
                  pred.fmm.S$classification,
                  pred.fmm.SD$classification,
                  rep(0, length(pred.fmm.M$classification)))

res.full <- rbind(res.train, res.test)
res.full <- cbind(c(stations.sub[sample == 1], stations.sub[sample == 0]), res.full)

colnames(res.full) <- c("Station", "M", "S", "SD", "Train")
write.csv2(res.full, file = "clusterMax.csv", row.names=FALSE)
##### End #####