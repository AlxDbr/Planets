# Goal of the project: we have discovered a lot of exoplanets, but it is
# difficult to estimate whether a discovered planet is telluric, like
# Mercury, Venus, Earth and Mars, is a gas giant like Jupiter and Saturn,
# or is an icy giant like Uranus and Neptune.
# The dataset of the csv that we will import doesn't contain this target column
# that can answer to our question
# --> We will use UNSUPERVISED LEARNING METHOD and create 3 clusters which
# will divide the exoplanets into these mentioned categories.
rm(list=ls(all=TRUE))
library(dplyr)
library(GGally)
library(factoextra)
# We must be sure to have a clean R environment, without other R objects
# We must also be sure to have the packages dplyr, GGally and factoextra
# (if ggplot2 is not installed, it's necessary to add it too)

# 1 - DATA EXPLORATION PHASE
planets=read.table(file="C:/Users/u515192/R_folders/data/exoplaneteu_catalog-1.csv",header=TRUE,dec=".",sep=";",quote="\"")
# Import of the dataset containing the exoplanets and their characteristics as a dataframe named planets
head(planets)
dim(planets)
# Our dataset contains 5159 lines and 98 columns
sapply(planets,class)
table(sapply(planets,class))
# Among the 98 columns, 1 contains integer variables (Discovery), 83 contains numeric variables and 14 contains characters
colSums(is.na(planets))
head(planets$Molecules)
head(planets$Log.G..g.)
# Seeing the first lines of "Molecules", we see that there's no NA for this reason:
# when we don't know molecules that compose the the planet, an empty string of
# characters is returned
planets[planets==""]=NA
number_na=colSums(is.na(planets))
perc_na = sort(number_na)/nrow(planets)
View(perc_na)
# allows us to see directly, during the development of the script, all the
# variables, sorted according to the number of NA



# 2 - FEATURES SELECTION
# We create planets_b to allow us to find data that are in columns that we have
# deleted, they can always be found in the dataset planets.
planets_b = planets %>% select(-contains("error"),-contains("star"),-contains("Mag"))
names(planets_b)
# We have removed all rows that concern error on measures because this will not
# speak to us when we will create our clusters, 
# We also do not need the data concerning the stars around which each of the
# planets revolve, because we already have the variables that will be affected,
# mainly the temperature and the orbital period
# We can also remove the 5 columns Mag. that give the magnitude of this star,
# for the same reason
# We there are 34 columns left out
# of the 98 concerned
planets_b = planets_b %>% mutate (Mass=pmax(Mass..M.Jup.,Mass.sini..M.Jup...M.Earth., na.rm=TRUE))
planets_b = planets_b %>% select(-Mass..M.Jup.,-Mass.sini..M.Jup...M.Earth.)
# We have defined Mass as the max between its estimated value and Mass.sini...
# that represents the minimal mass of the planet due to inclinaison effect.
# Habitually, the first mass is greater and represent the best estimation,
# but if the minimal mass is greater, we will prefer keep this one
head(planets_b)
planets_b[!is.na(planets$Albedo),c("Name","Albedo")]
# Albedo is the reflective power of a surface; it was calculated in particular
# on the Kepler planets, we are going to delete this data which is not 
# sufficiently represented in the dataset but it could be interesting to see
# the clusters in which we will find these planets if we want to get an idea
# of the influence of the Albedo on the category of the planet
planets_b = planets_b %>% select(-Albedo)
head(planets_b)
# The clustering should also not be influenced by the year of discovery, the
# update, the time of passage at the periapse (T.peri..JD.), time of star-planet
# upper conjunction (T.conj..JD.), other observations times (T0..JD., T0.sec..JD.,
# TVR..JD.), the planet status, the publication, the alternative names
# and the detection, mass detection and radius detection types.
#
# We will also remove the right ascention and declination angles, which are
# angles of observation when seeing the planet from Earth, and the longitude
# of the periapsis from which it's calculated, as well as the 
# lambda angle (Rossiter-McLaughlin anomaly)
#
# Let's see what we have with this :
planets_b = planets_b %>% select(-c("Discovery","Update","T.peri..JD.","T.conj..JD.","T0..JD.","T0.sec..JD.","TVR..JD.","Planet.Status","publication","Alternate.names","Detection.type","Mass.detection.type","Radius.detection.type","RA.α..deg.","Dec.δ..deg.","ω..deg.","λ.angle..deg."))
# We will create a variable Temperature that will be primarily given by the
# measured temperature, and if this value is unknown, the calculated one
planets_b = planets_b %>% mutate(Temperature=case_when(!is.na(Temp.measured..K.)~Temp.measured..K.,!is.na(Temp.calculated..K.)~Temp.calculated..K.))
planets_b = planets_b %>% select(-Temp.measured..K.,-Temp.calculated..K.)
head(planets_b)
perc_na2 = sort(colSums(is.na(planets_b)))/nrow(planets_b)
View(perc_na2)
# We see that there are still 4 columns with more than 80% of NA; we will not
# get a lot of information out of it and we will delete it, but il will
# be interesting to see in the old dataset (planets) the molecular compositions
# and logarithm of gravity that we often have for each cluster
planets_b = planets_b[,colSums(is.na(planets_b))/nrow(planets_b) < 0.8]
head(planets_b)
# We will redefine the names of the variables (we will not need to precise the
# units of the values because we will do scaling to give to the features the
# same importance before doing clustering
names(planets_b)=c("Name","Radius","Period","SemiMajor_axis","Eccentricity","Inclination","Impact_param","Radial_velocity","Mass","Temperature")
summary(planets_b)

# 3 - PREPROCESSING PHASE:
planets_b[is.na(planets_b$Period..day.),]
planets_datas = planets_b %>% select(-Name)
head(planets_datas)
ggpairs(planets_datas)
# At first, we will remove the 3 exoplanets that we can see as outliers in the
# pairplot:
planets_b = planets_b %>%
filter(Period < 7000000 | is.na(Period)) %>%
filter(Impact_param < 50 | is.na(Impact_param)) %>%
filter(Radial_velocity > -150000 | is.na(Radial_velocity))


# Correlation between Period and SemiMajor_axis:
planets_b1 = planets_b %>% filter(Period < 50 | is.na(Period))
dim(planets_b1)
ggplot(planets_b1) + aes(x=Period, y=SemiMajor_axis) + geom_point() + coord_cartesian(ylim=c(0,0.5))
# We can see that, for the exoplanets with a period < 50, that represent almost
# 80% of then, SemiMajor_axis seems to have a dependance according to a Period
# put to a power <1.
modele_1 = nls(SemiMajor_axis ~ a*Period^b, data=planets_b1, start=list(a=0.02,b=0.65))
modele_1
# We can see that the best power regression model estimate SemiMajor_axis
# as 0.01915*Period^0.66014.
# The mean squared error is +/- 1.5, which is strogly preferable to what we
# would have by taking the high-period planets
# If we have a value only for 1 of the 2 features, we will 
# estimate the other one by using this estimation:
a = 0.01915
b = 0.66014
planets_b = planets_b %>% mutate(SemiMajor_axis = case_when(is.na(SemiMajor_axis) ~ a*Period^b, TRUE ~ SemiMajor_axis))
# if y=a*x^b, x=(y/a)^(1/b)
# Some platets with unknown Period have a high SM_axis; to avoid non-realistic
# values, we have capped the period at a value of 200000
planets_b = planets_b %>% mutate(Period = case_when(is.na(Period) ~ pmin((SemiMajor_axis/a)^(1/b),200000), TRUE ~ Period))
ggplot(planets_b) + aes(x=Period, y=SemiMajor_axis) + geom_point() + coord_cartesian(xlim=c(0,200000),ylim=c(0,65))


# Anti-correlation between Inclination and Impact_param:
ggplot(planets_b) + aes(x=Inclination, y=Impact_param) + geom_point()
ggplot(planets_b) + aes(x=Inclination, y=Impact_param) + geom_point() + coord_cartesian(xlim=c(75,95))
# We can see at the first plot that the majority of planet orbits have an
# inclination between 75 and 90°, and an impact parameter from a to 1, with a
# that varies linearly with the inclination.
# We can also estimate that, for a given inclination, the mean impact param
# is the "third of the way" between a and 1.
# So, we will estimate Impact_param as: 1 if Inclination < 75
#						    (97.5-Inc.)/(97.5-75) if 75 < Inclination < 97.5
#						    0 if Inclination > 97.5
# Conversely, if we only know Impact_param, we will estimate Inclination as
# 90 - 5*Imp.param
planets_b = planets_b %>% mutate(Impact_param = case_when(is.na(Impact_param) ~ pmax(pmin((97.5-Inclination)/(97.5-75),1),0), TRUE ~ Impact_param))
planets_b = planets_b %>% mutate(Inclination = case_when(is.na(Inclination) ~ 90-5*Impact_param, TRUE ~ Inclination))
ggplot(planets_b) + aes(x=Inclination, y=Impact_param) + geom_point() + coord_cartesian(xlim=c(75,95))
# We clearly can see at this plot that our 2 linear predictions to replace the
# NA values have been realised


# Correlation between Radial_velocity and Mass:
ggplot(planets_b) + aes(x=Radial_velocity, y=Mass) + geom_point()+ coord_cartesian(xlim=c(0,1000), ylim=c(0,25))
# We can clearly observe that the variation of the Mass with the velocity can
# take 2 different behaviours, and we will try to find out which parameter can
# help us to estimate in which of the 2 cases we are; the threshold will make
# this estimate as precise as possible.
planets_b1 = planets_b %>% 
filter(!is.na(Mass) & !is.na(Radial_velocity)) %>%
filter(Mass <= 0.02*Radial_velocity)
dim(planets_b1)
dim(planets_b1 %>% filter(Period <= 50))
dim(planets_b1 %>% filter(Period > 50))
planets_b2 = planets_b %>% 
filter(!is.na(Mass) & !is.na(Radial_velocity)) %>%
filter(Mass > 0.02*Radial_velocity)
dim(planets_b2)
dim(planets_b2 %>% filter(Period <= 50))
dim(planets_b2 %>% filter(Period > 50))
# We see that when Period is <=50, in 93% of cases, we are in the sample
# planets_b1, and when period is >50, in 93% of cases, we are in the sample
# planets_b2. So, 50 for period is a very good threshold to separate our dataset
# into 2 datasets, and find the best regression for each of then
ggplot(planets_b) + aes(x=Radial_velocity, y=Mass, colour=Period<=50) + geom_point()+ coord_cartesian(xlim=c(0,1000), ylim=c(0,25))
test_1 = planets_b %>% filter(is.na(Mass) & !is.na(Radial_velocity))
dim(test_1)
test_1
test_2 = planets_b %>% filter(!is.na(Mass) & is.na(Radial_velocity))
dim(test_2)
# We can see that there are only 2 exoplanets that have an unknown mass but a
# known radial velocity. So, we will not need to find a model that estimate
# the mass, but simply estimate the 1412 unknown radial velocity's for which
# we know the mass.
planets_b1 = planets_b %>% filter(Period <= 50 & !is.na(Mass) & !is.na(Radial_velocity))
dim(planets_b1)
planets_b1 = planets_b %>% 
filter(Period <= 50) %>% 
filter(Mass < 0.02*Radial_velocity) %>%
filter(Radial_velocity < 900)
dim(planets_b1)
# We can see that, for the exoplanets with a period < 50, that represent almost
# 95% of then, and we see on the plot that Mass seems to have a linear dependance
# with Radial_velocity
planets_b2 = planets_b %>% filter(Period > 50 & !is.na(Mass) & !is.na(Radial_velocity))
dim(planets_b2)
planets_b2 = planets_b %>% 
filter(Period > 50) %>%
filter(Mass < Radial_velocity) %>%
filter(Radial_velocity < 1000)
dim(planets_b2)
# This filter select more than 95% of exoplanets with a period > 50, and we also
# see a linear dependance between Mass and Radial_velocity, with a higher coefficient
ggplot(planets_b1) + aes(x=Mass, y=Radial_velocity) + geom_point()
ggplot(planets_b2) + aes(x=Mass, y=Radial_velocity) + geom_point()
modele_b1 = lm(Radial_velocity~Mass, data=planets_b1)
modele_b1
a = 108.444
b = 9.321
modele_b2 = lm(Radial_velocity~Mass, data=planets_b2)
modele_b2
cc = 14.49 # c should define a column
d = 22.24
# We can see that the best linear regression model estimate Radial_velocity
# as 108.444*Mass + 9.321 for exoplanets with period < 50,
# and 14.49*Mass + 22.24 for exoplanets with period > 50
planets_b = planets_b %>% mutate(Radial_velocity = case_when(is.na(Radial_velocity)&Period <= 50 ~ a*Mass+b, is.na(Radial_velocity)&Period > 50 ~ cc*Mass+d, TRUE ~ Radial_velocity))
ggplot(planets_b) + aes(x=Radial_velocity, y=Mass, colour=Period<=50) + geom_point()+ coord_cartesian(xlim=c(0,1000), ylim=c(0,25))
# We clearly can see at this plot that our 2 linear predictions to replace the
# NA values for radial velocity have been realised


# Correlation between Radius and Temperature:
ggplot(planets_b) + aes(x=Radius, y=Temperature) + geom_point()
# The variation of the temperature with the radius can also take 2 different 
# behaviours, and we can see at the next figure that we can separate these 2 
# cases by defining a threshold of 0.1 for the mass:
ggplot(planets_b) + aes(x=Radius, y=Temperature, colour=Mass<=0.1) + geom_point()
planets_bbis = planets_b %>% filter(Mass>0.1)
ggplot(planets_bbis) + aes(x=Radius, y=Temperature) + geom_point()
# If the mass is lower than 0.1, we can't find any coherent relation, but we
# can at least make an estimation by linear regression if it is > 0.1:
modele_c = lm(Temperature~Radius, data=planets_bbis)
modele_c
# Best linear estimation calculated for Temperature: 874*Radius + 358.3
a = 874
b = 358.3
planets_b = planets_b %>% mutate(Temperature = case_when(is.na(Temperature)& Mass>0.1 ~ a*Radius+b, TRUE ~ Temperature))
# if y=a*x+b, x=(y/a)-(b/a)
planets_b = planets_b %>% mutate(Radius = case_when(is.na(Radius)& Mass>0.1 ~ (1/a)*Temperature-(b/a), TRUE ~ Radius))
planets_bbis = planets_b %>% filter(Mass>0.1)
ggplot(planets_bbis) + aes(x=Radius, y=Temperature) + geom_point()


# Feature scaling:
perc_na2 = sort(colSums(is.na(planets_b)))/nrow(planets_b)
View(perc_na2)
# The proportion of NA has passed from 0.08 to 0.02 for Period, from 0.33 to
# 0.02 for SemiMajor_axis, from 0.75 to 0.5 for Radial_velocity from 0.72 to
# 0.55 for inclination, from 0.78 to 0.7 for Temperature, from 0.66 to 0.55
# for Impact_param and from 0.27 to 0.26 for Radius.
# The last NA's will be replaced by 0 for the following reason:
# We must do scaling because we must be aware that the distance of metrics
# is unweighted according to the feature
# All these scaled features will follow a normal distribution, with a mean of
# 0 and a standard deviation of 1. So, it will be very simple to estimate the
# last NA's as we can by calculating their average values!
planets_datas = planets_b %>% select(-Name)
planets_names = planets_b %>% select(Name)
planets_datas_scale = as.data.frame(scale(planets_datas))
ggpairs(planets_datas_scale)
planets_datas_scale[is.na(planets_datas_scale)] = 0
head(planets_datas_scale)
head(planets_names)
# It could be interesting not to crush the dataset planets_b if we are looking
# for a value that is not standardized, so we will create a new dataset named 
# planets_c and keep planets_b
planets_c = bind_cols(planets_names,planets_datas_scale)

head(planets_c)

# 4 - CLUSTERING:
# We will at first calculate the optimal number of clusters with the Elbow method:
fviz_nbclust(planets_datas_scale, kmeans, method="wss") +
  labs(subtitle = "Elbow method")
# We see that the best compromise to minimize distance between points and the
# number of categories is to use 4 clusters. Let's use the KMeans algorithm now:
datas_clust = kmeans(planets_datas_scale,centers=4,nstart=100)
class(datas_clust)
class(datas_clust$cluster)
rownames(planets_datas_scale) = planets_names$Name
fviz_cluster(list(data=planets_datas_scale, cluster=datas_clust$cluster))
planets_with_cluster = bind_cols(planets_c,as.data.frame(datas_clust$cluster))
head(planets_with_cluster)
planets_with_cluster = planets_with_cluster %>% rename("Cluster"="datas_clust$cluster")
head(planets_with_cluster)
# Test: we will see if the cluster of the planets Kepler-539 c, HD 5388 b, 
# Kepler-1698 b and Ross 19 b is different, and we will then be able to guess
# their nature by observing the clusters of certain planets whose we know the
# molecular composition
(planets_with_cluster %>% filter(Name=="Kepler-539 c"))$Cluster
(planets_with_cluster %>% filter(Name=="HD 5388 b"))$Cluster
(planets_with_cluster %>% filter(Name=="Kepler-1698 b"))$Cluster
(planets_with_cluster %>% filter(Name=="Ross 19 b"))$Cluster
planets %>% filter(!is.na(Molecules)) %>% select("Name","Molecules")

# The telluric planets are composed in particular of iron, magnesium and silicon,
# Among the 102 planets for which we have identified the molecules, we can
# suppose that KELT-9 b, HAT-P-70 b, WASP-12 b, WASP-189 b, HD 189733 b,
# HD 209458 b, WASP-76 b, MASCARA-2 b/KELT-20 b, WASP-121 b and WASP-33 b are
# telluric.
(planets_with_cluster %>% filter(Name=="KELT-9 b"))$Cluster
(planets_with_cluster %>% filter(Name=="HAT-P-70 b"))$Cluster
(planets_with_cluster %>% filter(Name=="WASP-12 b"))$Cluster
(planets_with_cluster %>% filter(Name=="WASP-189 b"))$Cluster
(planets_with_cluster %>% filter(Name=="HD 189733 b"))$Cluster
(planets_with_cluster %>% filter(Name=="HD 209458 b"))$Cluster
(planets_with_cluster %>% filter(Name=="WASP-76 b"))$Cluster
(planets_with_cluster %>% filter(Name=="MASCARA-2 b/KELT-20 b"))$Cluster
(planets_with_cluster %>% filter(Name=="WASP-121 b"))$Cluster
(planets_with_cluster %>% filter(Name=="WASP-33 b"))$Cluster

# The icy giants planets are composed of volatile elements such as methane,
# ammonia and water, we can suppose that HR 8799 b, Ross 458 (AB) c and
# WISE J0458+6434 A are icy giants.
(planets_with_cluster %>% filter(Name=="HR 8799 b"))$Cluster
(planets_with_cluster %>% filter(Name=="Ross 458 (AB) c"))$Cluster
(planets_with_cluster %>% filter(Name=="WISE J0458+6434 A"))$Cluster

# The gas giants planets are mainly composed of hydrogen and helium, we can
# suppose that GJ 3470 b, HAT-P-11 b, HAT-P-18 b and GJ 1214 b are gas giants.
(planets_with_cluster %>% filter(Name=="GJ 3470 b"))$Cluster
(planets_with_cluster %>% filter(Name=="HAT-P-11 b"))$Cluster
(planets_with_cluster %>% filter(Name=="HAT-P-18 b"))$Cluster
(planets_with_cluster %>% filter(Name=="GJ 1214 b"))$Cluster

# We can observe that the cluster at the bottom right of the figure contains
# mostly telluric planets, the one at the bottom left the gas giants, and the
# one at the top left (and we can also include the last one which is closer),
# the icy giants.
# To answer the 4 questions of our earlier test, we can assume that Kepler-539
# c is a gas giant, HD 5388 b and Ross 19 b are icy giants and Kepler-1698 b
# is telluric !

