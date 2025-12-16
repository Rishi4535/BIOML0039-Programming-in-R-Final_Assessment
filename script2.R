#Setting up our working directory
setwd("C:/Rproject")

#loading all the required libraries
library(vroom)
library(dplyr)
library(tidyr)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(readxl)
library(rstatix)
library(dplyr)
library(performance)
library(fitdistrplus)
library(effectsize)
library(car)
library(lme4)
library(ggeffects)
library(DHARMa)
library(lattice)
library(forcats)

#listing the files in our working directory
files <- fs::dir_ls(path="C:/Rproject/mass_lv25974")
#listing files in the folder
files
#loading the sample data for analysis
analysis_data <- vroom(files, id="C:/Rproject/mass_lv25974/mass_lv25974.tsv")
#checking object type and contents
class(analysis_data)
#summarising the dataframe
describe(analysis_data)
summary(analysis_data)

#cleaning our data by getting rid of unnecessary columns
analysis_data$`C:/Rproject/mass_lv25974/mass_lv25974.tsv`<-NULL
analysis_data$`student_id`<-NULL
#Modifying data entries to keep data homogenous
analysis_data$lion_presence[analysis_data$lion_presence=="Asbent"]<-"Absent"
analysis_data$lion_presence[analysis_data$lion_presence=="Preesnt"]<-"Present"
analysis_data$lion_presence[analysis_data$lion_presence=="Preseet"]<-"Present"


#As per the given instructions, sequential triplicates of data are considered to be of the same individual
#Removing triplicate value and reducing them to one observation point
#Retaining other column values and try to make 2nd last data column entries as columns and assign data to it

transformed_data <- analysis_data %>%
  pivot_wider(
    id_cols     = 1:7,        # retain first 7 columns exactly
    names_from  = 8,         # column with 3 specific entries
    values_from = 9          # column with values to assign
  )

#dataframe summary
class(transformed_data)
describe(transformed_data)
summary(transformed_data)
#visualising our trasnformed data
hist(transformed_data$body_mass)


#creating a new dataset for model fit andanalysis
#re-assigning NA as unknowns so model will still consider these values and not ignore them wile modelling
df <- transformed_data %>%
  mutate(
    landscape = fct_explicit_na(landscape, na_level = "Unknown"),
    vegetation = fct_explicit_na(vegetation, na_level = "Unknown"),
    lion_presence = fct_explicit_na(lion_presence, na_level = "Unknown")
  )
#indicating categorical variables so they are not taken as numerical values for modelling
df <- df %>%
  mutate(
    park = factor(park),
    observer = factor(observer),
    landscape = factor(landscape),
    vegetation = factor(vegetation),
    lion_presence = factor(lion_presence)
  )
#missing data in predictor variables should not be deleted, rather imputed with median of the data column to preserve data frame
df <- df %>%
  mutate(
    across(
      c(distance_to_water_km, average_annual_rainfall_mm, predator_density_per_sq_km),
      ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
    )
  )
#removing outliers, body mass > 2000kg as maximum recordedweight of giraffe is around 2000kg
df <- df %>%
  filter(body_mass <= 2000)
#log of body mass is taken as extra data point for modelling analysis
df <- df %>%
  mutate(log_body_mass = log(body_mass))

#dataframe summary
class(df)
describe(df)
summary(df)
#visualising distrubution of body mass
hist(df$body_mass)


#visualisation of body mass across various categorical variables or factors
factor_vars <- df %>%
  select_if(where(is.factor)) %>%
  names()
#plotting
for (var in factor_vars) {
  
  p <- ggplot(df, aes_string(x = var, y = "body_mass")) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", size = 4, colour = "red") +
    labs(
      x = var,
      y = "Body mass (kg)",
      title = paste("Giraffe body mass across", var)
    ) +
    theme_bw()
  
  print(p)
}


#visualisation of body mass over environmental/ predictor variables
ggplot(df, aes(x = distance_to_water_km, y = body_mass)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Distance to water (km)",
    y = "Body mass (kg)",
    title = "Relationship between body mass and distance to water"
  )

ggplot(df, aes(x = average_annual_rainfall_mm, y = body_mass)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Average annual rainfall (mm)",
    y = "Body mass (kg)",
    title = "Relationship between body mass and average annual rainfall"
  )

ggplot(df, aes(x = predator_density_per_sq_km, y = body_mass)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Predator Density (sq.km)",
    y = "Body mass (kg)",
    title = "Relationship between body mass and Predator Density"
  )

#fitting a distribution over our data
f_norm <- fitdist(df$body_mass, "norm")
gofstat(f_norm)
plot(f_norm)

#fitting a model on our data to get predictions
m1 <- lmer(
  body_mass ~ landscape + vegetation +
    distance_to_water_km + average_annual_rainfall_mm +
    predator_density_per_sq_km + lion_presence +
    (1 | park) + (1 | observer),
  data = df,
  REML = FALSE
)

#model summary and fit evaluation
summary(m1)
standardize_parameters(m1)
plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))
check_predictions(m1)
# Check assumptions with DHARMa
sim_res <- simulateResiduals(m1, plot = FALSE)
plot(sim_res)  # residual diagnostics 
testDispersion(m1)


#predictions
pred <- ggpredict(m1, terms = "distance_to_water_km")
plot(pred) +
  labs(
    x = "Distance to Water (Km)",
    y = "Predicted body mass",
    title = "Distance to water sources effect on giraffe body mass"
  ) +
  theme_bw()

#predictions
pred <- ggpredict(m1, terms = "average_annual_rainfall_mm")
plot(pred) +
  labs(
    x = "average_annual_rainfall_mm",
    y = "Predicted body mass",
    title = "average_annual_rainfall_mm effect on giraffe body mass"
  ) +
  theme_bw()

#predictions
pred <- ggpredict(m1, terms = "predator_density_per_sq_km")
plot(pred) +
  labs(
    x = "predator_per_sq_km",
    y = "Predicted body mass",
    title = "predator_per_sq_km effect on giraffe body mass"
  ) +
  theme_bw()

glimpse(df)
