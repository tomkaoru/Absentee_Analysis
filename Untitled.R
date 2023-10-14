
setwd('/DIRECTORY')

# libraries
library("tidyverse")
library("sf") # sf object
library("tmap") # map
library("readxl") #read excel
library(stringr) #change strings
library("spdep") #convert sf to sp
library("rstan") #stan
library("geostan") #adjacency matrix
library("SpatialEpi") #compute expected number
library("tidybayes") #clean model

# use multiple cores
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# #make sure stan is latest version
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# read shapefile
jpn <- read_sf('shapefile/jpn_admbnda_adm1_2019.shp')
jpn <- st_make_valid(jpn)

tm_shape(jpn) + tm_polygons() # map shapefile
# remove unnecessary columns
jpn <- jpn[-c(1:3,5:6)]
# convert accented prefectures (e.g. Ōita, Kōchi) to nonaccented
jpn[13, 1] = "Hyogo"
jpn[20, 1] = "Kochi"
jpn[30, 1] = "Oita"
jpn[10, 1] = "Gumma" # change misspelling
colnames(jpn)[1]  <- "area" 
jpn$area <- str_trim(jpn$area, "left")

# read absentee data
absentee <- read_excel("data/absentees.xls",)
absentee <- absentee[ , c(10, 99,166)] #only extract neessary columns
# column 10 contains name of prefecture
# column 99 is # of elementary school students for 2015
# column 66 is 	Number of long-term absentees from elementary school for 2015
absentee <- tail(absentee, -11) # remove unnecessary row
absentee <- head(absentee, - 2)  # remove unnecessary row
names(absentee) <- c('area','total','absentee') #change col nname

absentee$area <- gsub("\\-ken", "", absentee$area) #remove Japanese word 'ken' from prefecture name e.g. Hiroshima-ken -> Hiroshima
absentee$area <- gsub("\\-to", "", absentee$area) #remove Japanese word 'to' from prefecture name e.g. Tokyo-to -> Tokyo
absentee$area <- gsub("\\-fu", "", absentee$area) #remove Japanese word 'fu' from prefecture name e.g. Osaka-fu -> Osaka

# read expenditure data
expenditure <- read.csv(file = 'data/expenditure.csv', skip = 11)
expenditure <- expenditure[expenditure[,2] == 2015,] # only keep data from 2015
expenditure <- tail(expenditure, - 1) # remove unnecessary row
expenditure <- expenditure[ , c(4,49)] # keep necessary columns
# column 4 prefecture name'
# column 49 elementary school expenditure per student
names(expenditure) <- c('area','expenditure')
expenditure$area <- gsub("\\-ken", "", expenditure$area)
expenditure$area <- gsub("\\-to", "", expenditure$area)
expenditure$area <- gsub("\\-fu", "", expenditure$area)
expenditure$expenditure <- gsub(",", "", expenditure$expenditure) #remove "," e.g. 1,030 -> 1030

# read income data
income <- read.csv(file = 'data/income.csv', skip = 11)
income <- income[income[,2] == 2015,] # only keep data from 2015
income <- tail(income, - 1) # remove unnecessary row
income <- income[ , c(4,9)] # keep necessary columns
# column 4 prefecture name'
# column 9 average household income
names(income) <- c('area','income')
income$area <- gsub("\\-ken", "", income$area)
income$area <- gsub("\\-to", "", income$area)
income$area <- gsub("\\-fu", "", income$area)

#merged all data at prefecture level including the shapefile
merged <- merge(jpn, absentee, by.x = "area", by.y = "area",all.x=TRUE)
merged <- merge(merged, expenditure, by.x = "area", by.y = "area",all.x=TRUE)
merged <- merge(merged, income, by.x = "area", by.y = "area",all.x=TRUE)

# convert columns to numeric
merged[, c(2,3,4,5)] <- lapply(c(2,3,4,5), function(x) as.numeric(merged[[x]]))

# calculate expected number of absentees for each prefecture
merged$expectedabsentee <- round(expected(population = merged$total, cases = merged$absentee, n.strata = 1), 0)

# only keep prefectures in Honshu island by removing prefectures not in Honshu
# [Fukuoka,Hokkaido,Okinawa,Saga,Oita,Mizayaki,Kumamoto,Nagasaki,Kagoshima,Ehime, Kochi, Tokushima, Kagawa]
merged_cleaned = merged[-c(5,7,12,17,18,20,21,25,27,30,32,34,40),]
rownames(merged_cleaned) <- NULL 
hist(merged_cleaned$absentee, breaks = 200) #histogram of dependent variable
# Data cleaning is done

tokyo <- merged_cleaned[28,]
osaka <- merged_cleaned[22,]
tm_shape(merged_cleaned) + tm_polygons() + tm_shape(tokyo) + tm_fill('red') + tm_shape(osaka) + tm_fill('green') + 
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


# CREATE SPATIAL INFORMATION
# reorder columns so geometry come last
merged_cleaned <- merged_cleaned[, c(1,2,3,4,5,7,6)]
# change scale of independent variables
merged_cleaned$expenditure <- merged_cleaned$expenditure / 100
merged_cleaned$income <- merged_cleaned$income / 100

sp_merged <- as(merged_cleaned, "Spatial") # convert to spatial object
adjacencyMatrix <- shape2mat(sp_merged)
# get adjacency matrix and prepare data for icar model
extractComponents <- prep_icar_data(adjacencyMatrix)

n <- as.numeric(extractComponents$group_size) # number of observations
nod1 <- extractComponents$node1 
nod2 <- extractComponents$node2
n_edges <- as.numeric(extractComponents$n_edges)
y <- merged_cleaned$absentee # dependent variable
x <- st_drop_geometry(merged_cleaned[c(4,5)])# independent variables
e <- merged_cleaned$expectedabsentee #offset

stan_data <- list(N=n, N_edges=n_edges, node1=nod1, node2=nod2, Y=y, X=x, E=e, k=2)
# k is the # of independent variable

icar_poisson_fit = stan("Untitled_1.stan", data=stan_data, iter=20000, chains=6, verbose = FALSE)

options(scipen = 999) # remove scientific notation 
summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))$summary # estimates not exponentiated yet

# area
# summary(icar_poisson_fit, pars=c("phi"), probs=c(0.025, 0.975))$summary

# ANALYSIS
print(icar_poisson_fit, pars=c("alpha", "beta", "sigma", "phi"), probs=c(0.025, 0.975))

# RELATIVE RISK  - - - - - - - - - - - - - -
# get exponentiated risks
rr_result <- as.data.frame(summary(icar_poisson_fit, pars=c("mu"), probs=c(0.025, 0.975))$summary)
row.names(rr_result) <- 1:nrow(rr_result)
rr_result <- rr_result[, c(1,4,5,7)]
colnames(rr_result)[1] <- "rr"
colnames(rr_result)[2] <- "rrlower"
colnames(rr_result)[3] <- "rrupper"
colnames(rr_result)[4] <- "rHAT"

# view clean table 
head(rr_result)

# now, we proceed to generate our risk maps
# align the results to the areas in shapefile
sp_merged$rr <- rr_result[, "rr"]
sp_merged$rrlower <- rr_result[, "rrlower"]
sp_merged$rrupper <- rr_result[, "rrupper"]

# create categories to define if an area has significant increase or decrease in risk, or nothing all (with 95% cconfidence)
sp_merged$Significance <- NA
sp_merged$Significance[sp_merged$rrlower<1 & sp_merged$rrupper>1] <- 0    # NOT SIGNIFICANT
sp_merged$Significance[sp_merged$rrlower==1 | sp_merged$rrupper==1] <- 0  # NOT SIGNIFICANT
sp_merged$Significance[sp_merged$rrlower>1 & sp_merged$rrupper>1] <- 1    # SIGNIFICANT INCREASE
sp_merged$Significance[sp_merged$rrlower<1 & sp_merged$rrupper<1] <- -1   # SIGNIFICANT DECREASE

summary(sp_merged$rr)
hist(sp_merged$rr)

RiskCategorylist <- c(">0.0 to 0.25", "0.26 to 0.50", "0.51 to 0.75", "0.76 to 0.99", "1.00 & <1.01",
                      "1.01 to 1.10", "1.11 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.00", "2.01 to 3.00")

RRPalette <- c("#65bafe","#98cffe","#cbe6fe","#dfeffe","white","#fed5d5","#fcbba1","#fc9272","#fb6a4a","#de2d26","#a50f15")

sp_merged$RelativeRiskCat <- NA
sp_merged$RelativeRiskCat[sp_merged$rr>= 0 & sp_merged$rr <= 0.25] <- -4
sp_merged$RelativeRiskCat[sp_merged$rr> 0.25 & sp_merged$rr <= 0.50] <- -3
sp_merged$RelativeRiskCat[sp_merged$rr> 0.50 & sp_merged$rr <= 0.75] <- -2
sp_merged$RelativeRiskCat[sp_merged$rr> 0.75 & sp_merged$rr < 1] <- -1
sp_merged$RelativeRiskCat[sp_merged$rr>= 1.00 & sp_merged$rr < 1.01] <- 0
sp_merged$RelativeRiskCat[sp_merged$rr>= 1.01 & sp_merged$rr <= 1.10] <- 1
sp_merged$RelativeRiskCat[sp_merged$rr> 1.10 & sp_merged$rr <= 1.25] <- 2
sp_merged$RelativeRiskCat[sp_merged$rr> 1.25 & sp_merged$rr <= 1.50] <- 3
sp_merged$RelativeRiskCat[sp_merged$rr> 1.50 & sp_merged$rr <= 1.75] <- 4
sp_merged$RelativeRiskCat[sp_merged$rr> 1.75 & sp_merged$rr <= 2.00] <- 5
sp_merged$RelativeRiskCat[sp_merged$rr> 2.00 & sp_merged$rr <= 10] <- 6
table(sp_merged$RelativeRiskCat) # check the distribution of relative risk category

# map of relative risk
tm_shape(sp_merged) + tm_fill("RelativeRiskCat", style = "cat", title = "Relative Risk", palette = RRPalette, labels = RiskCategorylist) + 
  tm_shape(merged_cleaned) + tm_polygons(alpha = 0.05) + tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))

# map of significance regions
tm_shape(sp_merged) + 
  tm_fill("Significance", style = "cat", title = "Significance Categories", palette = c("#33a6fe", "white", "#fe0000"), labels = c("Significantly low", "Not Significant", "Significantly high")) +
  tm_shape(merged_cleaned) + tm_polygons(alpha = 0.05) + tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


# EXCEEDANCE PROBABILITY 
# extract the exceedence probabilities from the icar_possion_fit object
# calculate the probability that an area has a relative risk ratio greater than 1
threshold <- function(x){mean(x > 1.00)}
exceedance <- icar_poisson_fit %>% spread_draws(mu[i]) %>% 
  group_by(i) %>% summarise(mu=threshold(mu)) %>%
  pull(mu)

# insert the exceedance values into the spatial data frame
sp_merged$exceedance <- exceedance 

ProbCategorylist <- c("<0.01", "0.01-0.33", "0.34-0.66", "0.67-0.99", "1.00")

# categorising the probabilities in bands of 10s
sp_merged$ProbCat <- NA
sp_merged$ProbCat[sp_merged$exceedance>=0 & sp_merged$exceedance< 0.01] <- 1
sp_merged$ProbCat[sp_merged$exceedance>=0.01 & sp_merged$exceedance< 0.34] <- 2
sp_merged$ProbCat[sp_merged$exceedance>=0.34 & sp_merged$exceedance< 0.67] <- 3
sp_merged$ProbCat[sp_merged$exceedance>=0.67 & sp_merged$exceedance< 1.00] <- 4
# sp_merged$ProbCat[sp_merged$exceedance>=0.50 & sp_merged$exceedance< 0.60] <- 7
# sp_merged$ProbCat[sp_merged$eexceedance>=0.60 & sp_merged$exceedance< 0.70] <- 8
# sp_merged$ProbCat[sp_merged$exceedance>=0.70 & sp_merged$exceedance< 0.80] <- 9
# sp_merged$ProbCat[sp_merged$exceedance>=0.80 & sp_merged$exceedance< 0.90] <- 10
# sp_merged$ProbCat[sp_merged$exceedance>=0.90 & sp_merged$exceedansce< 1.00] <- 11
sp_merged$ProbCat[sp_merged$exceedance == 1.00] <- 5
table(sp_merged$ProbCat) # check if categories contain at least 1 prefecture

# map of exceedance probabilities
tm_shape(sp_merged) + 
  tm_fill("ProbCat", style = "cat", title = "Probability", palette = "GnBu", labels = ProbCategorylist) +
  tm_shape(merged_cleaned) + tm_polygons(alpha = 0.05) + tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom")) + tm_facets(drop.empty.facets = FALSE)
