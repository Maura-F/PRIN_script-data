library(dplyr)
library(lme4)
library(lmerTest)
library(car)
library(performance)
library(ggplot2)
library(tidyr)
library(betapart)
library(reshape2)
library(glmmTMB)

## Carico Dati grezzi:
# File senza variabili ambeintali
df_boreal       <-read.csv("dati\\df_boreal_subm.csv")
df_mediterranean<-read.csv("dati\\df_mediterranean_subm.csv")
df_oak          <-read.csv("dati\\df_oak_subm.csv")
df_beech        <-read.csv("dati\\df_beech_subm.csv")

# File con TUTTE le variabili
Predict_matrix<- read.csv("dati\\Env_var_matrix.csv")


### APLHA: 
# Preparo df
Boreal_matrix           <-df_boreal        %>%  distinct(Site, Year, .keep_all = TRUE)%>%  select(Site, Year, Richness, N.Sus, Forest.type)
Mediterranean_matrix    <-df_mediterranean %>%  distinct(Site, Year, .keep_all = TRUE)%>%  select(Site, Year, Richness, N.Sus, Forest.type)
Oak_matrix              <-df_oak           %>%  distinct(Site, Year, .keep_all = TRUE)%>%  select(Site, Year, Richness, N.Sus, Forest.type)
Beech_matrix            <-df_beech         %>%  distinct(Site, Year, .keep_all = TRUE)%>%  select(Site, Year, Richness, N.Sus, Forest.type)

#### Sp.RICHNESS~YEAR
LMM_Boreal          <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Boreal_matrix       )
LMM_Mediterranean   <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Mediterranean_matrix)
LMM_Oak             <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Oak_matrix          )
LMM_Beech           <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Beech_matrix        )

##Sp. Richness~Year FINO AL 2014
Boreal_matrix_2014          <- subset(Boreal_matrix, Year <= 2014)
Mediterranean_matrix_2014   <- subset(Mediterranean_matrix, Year <= 2014)
Oak_matrix_2014     <- subset(Oak_matrix, Year <= 2014)
Beech_matrix_2014   <- subset(Beech_matrix, Year <= 2014)

LMM_Boreal_14          <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Boreal_matrix_2014)
LMM_Mediterranean_14   <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Mediterranean_matrix_2014)
LMM_Oak_14             <- lmer(Richness ~ Year + (1 | Site), data = Oak_matrix_2014)
LMM_Beech_14           <- lmer(Richness ~ Year + (1 | Site) + (1 | N.Sus), data = Beech_matrix_2014)


#### 1. Standardizzo tutte le variabili indip.
exclude_cols <- c("X", "Site", "Year", "Forest_type", "N.Sus", "Richness") # Colonne da escludere dalla standardizzazione
cols_to_standardize <- setdiff(names(Predict_matrix), exclude_cols) # Seleziono le colonne da standardizzare
# Applico standardizzazione
Predict_matrix_std <- Predict_matrix
Predict_matrix_std[cols_to_standardize] <- lapply(
  Predict_matrix_std[cols_to_standardize],
  function(x) {
    x <- as.numeric(x)  
    scale(x)            # Standardizza (media = 0, sd = 1)
  })


# df dei singoli Biomi con variabili ambientali standardizzate, da utilizzare per i Modelli:
Boreal_matrix_std        <- subset(Predict_matrix_std, Forest_type == "Boreal forest")
Mediterranean_matrix_std <- subset(Predict_matrix_std, Forest_type == "Mediterranean forest")
Nemoral_oak_matrix_std   <- subset(Predict_matrix_std, Forest_type == "Nemoral oak forest")
Nemoral_beech_matrix_std <- subset(Predict_matrix_std, Forest_type == "Nemoral beech forest")


## 2. Sp.RICHNESS~ Env.Var.:
variabili_clima<- c("Mean_temp", "Annual_prec", "cv_prec", "tot_prec_growing_season", "mean_temp_growing_season", "CDD_annuale", "CDD_growing_season", "TX90p_annual", "Temp_seasonality")
formula_clima_lm <- as.formula(paste("Richness ~", paste(variabili_clima, collapse = " + ")))

variabili_struttura<- c("Tree_cover", "DefM", "Shrub_cover")
formula_struttura_lm <- as.formula(paste("Richness ~", paste(variabili_struttura, collapse = " + ")))
formula_struttura <- as.formula(paste("Richness ~", paste(variabili_struttura, collapse = " + "), "+ (1|Site)", "+ (1|N.Sus)"))

variabili_suolo <- c("K", "NH4", "SO4", "NO3", "pH")
formula_suolo_lm <- as.formula(paste("Richness ~", paste(variabili_suolo, collapse = " + ")))
formula_suolo <- as.formula(paste("Richness ~", paste(variabili_suolo, collapse = " + "), "+ (1|Site)"))

# Per ogni gruppo di variabili (suolo, clima, struttura), per ogni bioma faccio:
# -> A. Verifico la MULTICOLLINEARITA' (con VIF) delle variabili indipendenti
# -> B. Utilizzo lo STEPWHISE selection per selezionare le variabili predittive da mantenere nel modello 
# -> C. Verifico il rispetto delle assunzioni dei modelli

# 2.1. Clima:
lm_full_model <- lm(formula_clima_lm, data =Boreal_matrix_std)
full_model<-lmer(Richness~ cv_prec + tot_prec_growing_season + mean_temp_growing_season + CDD_annuale + CDD_growing_season + TX90p_annual + Temp_seasonality  + (1|Site) + (1|N.Sus), data = Boreal_matrix_std)

lm_full_model <- lm(formula_clima_lm, data =Mediterranean_matrix_std)
full_model<-lmer(Richness~ Annual_prec + cv_prec + tot_prec_growing_season + mean_temp_growing_season + CDD_growing_season + TX90p_annual + Temp_seasonality + (1|Site) + (1|N.Sus), data = Mediterranean_matrix_std)

lm_full_model <- lm(formula_clima_lm, data =Nemoral_beech_matrix_std )
full_model<-lmer(Richness~ Annual_prec + cv_prec + tot_prec_growing_season + mean_temp_growing_season + CDD_growing_season + CDD_annuale + TX90p_annual + Temp_seasonality + (1|Site) + (1|N.Sus), data = Nemoral_beech_matrix_std)

lm_full_model <- lm(formula_clima_lm, data =Nemoral_oak_matrix_std)
full_model<-lmer(Richness~  Annual_prec + cv_prec + tot_prec_growing_season + mean_temp_growing_season + CDD_growing_season + CDD_annuale + TX90p_annual + Temp_seasonality + (1|Site) + (1|N.Sus), data = Nemoral_oak_matrix_std)

#VIF
vif(lm_full_model)

#STEWISE SELECTION
step_model<- step(full_model, direction="both", trace=FALSE)
final_model <- get_model(step_model)
summary(final_model)


# 2.2. Struttura:
lm_full_model <- lm(formula_struttura_lm, data = subset(Boreal_matrix_std,  !is.na(DefM)))          
full_model<- lmer(formula_struttura, data =  subset(Boreal_matrix_std,          !is.na(DefM)))          

lm_full_model <- lm(formula_struttura_lm, data = subset(Mediterranean_matrix_std, !is.na(DefM)))    
full_model<- lmer(formula_struttura, data =  subset(Mediterranean_matrix_std,   !is.na(DefM)))    

lm_full_model <- lm(formula_struttura_lm, data = subset(Nemoral_beech_matrix_std, !is.na(DefM)))    
full_model<- lmer(formula_struttura, data =  subset(Nemoral_beech_matrix_std,   !is.na(DefM)))    

lm_full_model <- lm(formula_struttura_lm, data = subset(Nemoral_oak_matrix_std, !is.na(DefM)))      
full_model<- lmer(formula_struttura, data =  subset(Nemoral_oak_matrix_std,     !is.na(DefM)))      

#VIF
vif(lm_full_model)
#STEWISE SELECTION
step_model<- step(full_model, direction="both", trace=FALSE)
final_model <- get_model(step_model)
summary(final_model)


# 2.1. Suolo:
lm_full_model <- lm(formula_suolo_lm , data =subset(Boreal_matrix_std,       !is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(K) & !is.na(NO3)))
full_model<- lmer(formula_suolo, data = subset(Boreal_matrix_std, !is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(NO3)))

lm_full_model <- lm(formula_suolo_lm , data =subset(Nemoral_beech_matrix_std,!is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(K) & !is.na(NO3)))
full_model<- lmer(sqrt(Richness)~ NO3 + NH4 + SO4 + pH + (1|Site), data = subset(Nemoral_beech_matrix_std, !is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(NO3)))

lm_full_model <- lm(formula_suolo_lm , data =subset(Nemoral_oak_matrix_std,  !is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(K) & !is.na(NO3)))
full_model<- lmer(sqrt(Richness)~ K + NO3 + NH4 + SO4 + pH + (1|Site), data = subset(Nemoral_oak_matrix_std, !is.na(pH) & !is.na(NH4) & !is.na(SO4) & !is.na(NO3)))

#VIF
vif(lm_full_model)
#STEPWISE SELECTION
step_model<- step(full_model, direction="both", trace=FALSE)
print(step_model)
final_model <- get_model(step_model)
summary(final_model)


#%VARIANCE EXPLAINED by de model
r2_results <- r2(full_model) 
print(r2_results)


### VERIFICA ASSUNZIONI MODELLI
modello<-final_model
residui<- residuals(modello)
valori_predetti <- fitted(modello)

#Grafico dei residui vs valori predetti
plot(fitted(modello), residuals(modello))
#
ggplot(data = NULL, aes(x = valori_predetti, y = residui)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # Linea di tendenza
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()
#
plot(residuals(modello))
qqnorm(residuals(modello),main = "QQ Plot of Residuals")
qqline(residuals(modello), col = "red")
#
hist(residui, probability = TRUE,
     main = "Histogram of Residuals with Normal Curve",
     xlab = "Residuals", ylab = "Density")
curve(dnorm(x, mean = mean(residui), sd = sd(residui)),
      add = TRUE, col = "red", lwd = 2)
densita_empirica <- density(residui)
lines(densita_empirica, col = "blue", lwd = 2)





### ELLENBERG
# APRI
df_boreal_Ell       <-read.csv( "dati\\boreal_Ellemberg.csv"         )
df_mediterranean_Ell<-read.csv( "dati\\mediterranean_Ellemberg.csv"  )
df_oak_Ell          <-read.csv( "dati\\oak_Ellemberg.csv"            )
df_beech_Ell        <-read.csv( "dati\\beech_Ellemberg.csv"          )

## CALCOLO la MEDIA dei valori Ell (separatamente per sito e anno di ogni bioma)
df_boreal_Ell_mean <- df_boreal_Ell %>%
  group_by(Site_Year) %>%
  summarise(
    mean_Ell_Temp = mean(Ell_Temp, na.rm = TRUE),
    mean_Ell_Moist = mean(Ell_Moisture, na.rm = TRUE),
    mean_Ell_Light = mean(Ell_Light, na.rm = TRUE)
  )

df_mediterranean_Ell_mean <- df_mediterranean_Ell %>%
  group_by(Site_Year) %>%
  summarise(
    mean_Ell_Temp = mean(Ell_Temp, na.rm = TRUE),
    mean_Ell_Moist = mean(Ell_Moisture, na.rm = TRUE),
    mean_Ell_Light = mean(Ell_Light, na.rm = TRUE)
  )

df_oak_Ell_mean <- df_oak_Ell %>%
  group_by(Site_Year) %>%
  summarise(
    mean_Ell_Temp = mean(Ell_Temp, na.rm = TRUE),
    mean_Ell_Moist = mean(Ell_Moisture, na.rm = TRUE),
    mean_Ell_Light = mean(Ell_Light, na.rm = TRUE)
  )

df_beech_Ell_mean <- df_beech_Ell %>%
  group_by(Site_Year) %>%
  summarise(
    mean_Ell_Temp = mean(Ell_Temp, na.rm = TRUE),
    mean_Ell_Moist = mean(Ell_Moisture, na.rm = TRUE),
    mean_Ell_Light = mean(Ell_Light, na.rm = TRUE)
  )
# Genero le colonne "Site" e "Year"
df_boreal_Ell_mean <- df_boreal_Ell_mean %>%
  separate(Site_Year, into = c("Site", "Year"), sep = "_") %>%
  mutate(Year = as.integer(Year))

df_mediterranean_Ell_mean <- df_mediterranean_Ell_mean %>%
  separate(Site_Year, into = c("Site", "Year"), sep = "_") %>%
  mutate(Year = as.integer(Year))

df_oak_Ell_mean <- df_oak_Ell_mean %>%
  separate(Site_Year, into = c("Site", "Year"), sep = "_") %>%
  mutate(Year = as.integer(Year))

df_beech_Ell_mean <- df_beech_Ell_mean %>%
  separate(Site_Year, into = c("Site", "Year"), sep = "_") %>%
  mutate(Year = as.integer(Year))

# LMM Ellenberg nel tempo
# LIGHT ☀️
LMM_Boreal_Ell_LIGHT <- lmer(log(mean_Ell_Light) ~ Year + (1 | Site), data = df_boreal_Ell_mean)
summary(LMM_Boreal_Ell_LIGHT)
LMM_Mediterranean_Ell_LIGHT <- lmer(mean_Ell_Light ~ Year + (1 | Site), data = df_mediterranean_Ell_mean)
summary(LMM_Mediterranean_Ell_LIGHT) 
LMM_oak_Ell_LIGHT <- lmer(mean_Ell_Light ~ Year + (1 | Site), data = df_oak_Ell_mean)
summary(LMM_oak_Ell_LIGHT) 
LMM_beech_Ell_LIGHT <- lmer(sqrt(mean_Ell_Light) ~ Year + (1 | Site), data = df_beech_Ell_mean)
summary(LMM_beech_Ell_LIGHT) 

# MOISTURE 💧
LMM_Boreal_Ell_MOIST          <- lmer(mean_Ell_Moist ~ Year + (1 | Site), data = df_boreal_Ell_mean)
summary(LMM_Boreal_Ell_MOIST) 
LMM_Mediterranean_Ell_MOIST          <- lmer(mean_Ell_Moist ~ Year + (1 | Site), data = df_mediterranean_Ell_mean)
summary(LMM_Mediterranean_Ell_MOIST) 
LMM_oak_Ell_MOIST          <- lmer(mean_Ell_Moist ~ Year + (1 | Site), data = df_oak_Ell_mean)
summary(LMM_oak_Ell_MOIST) 
LMM_beech_Ell_MOIST          <- lmer(mean_Ell_Moist ~ Year + (1 | Site), data = df_beech_Ell_mean)
summary(LMM_beech_Ell_MOIST)

# TEMPERATURE 🌡️
LMM_Boreal_Ell_TEMPERATURE          <- lmer(mean_Ell_Temp ~ Year + (1 | Site), data = df_boreal_Ell_mean)
summary(LMM_Boreal_Ell_TEMPERATURE )  
LMM_Mediterranean_Ell_TEMPERATURE           <- lmer(mean_Ell_Temp ~ Year + (1 | Site), data = df_mediterranean_Ell_mean)
summary(LMM_Mediterranean_Ell_TEMPERATURE )
LMM_oak_Ell_TEMPERATURE           <- lmer(mean_Ell_Temp ~ Year + (1 | Site), data = df_oak_Ell_mean)
summary(LMM_oak_Ell_TEMPERATURE )
LMM_beech_Ell_TEMPERATURE           <- lmer(mean_Ell_Temp ~ Year + (1 | Site), data = df_beech_Ell_mean)
summary(LMM_beech_Ell_TEMPERATURE ) 



### BETA: 
# Creo matrici di presenza-assenza
presence_absence_matrix_boreal          <- dcast(df_boreal       ,  Site + Year ~ Species, value.var = "Presence", fun.aggregate = max, fill = 0)
presence_absence_matrix_mediterranean   <- dcast(df_mediterranean,  Site + Year ~ Species, value.var = "Presence", fun.aggregate = max, fill = 0)
presence_absence_matrix_oak             <- dcast(df_oak          ,  Site + Year ~ Species, value.var = "Presence", fun.aggregate = max, fill = 0)
presence_absence_matrix_beech           <- dcast(df_beech        ,  Site + Year ~ Species, value.var = "Presence", fun.aggregate = max, fill = 0)

## FUNZIONE per calcolare Turnover e Nestedness CUMULATIVI
calculate_turnover_nestedness_cumulative <- function(presence_absence_data) {
  sites <- unique(presence_absence_data$Site)
  beta_values <- data.frame(Site = character(), Year1 = integer(), Year2 = integer(), 
                            N.Sus = numeric(), Turnover = numeric(), Nestedness = numeric())
  for (site in sites) {  # Loop per ogni sito
    site_subset <- presence_absence_data[presence_absence_data$Site == site, ]
    years <- sort(unique(site_subset$Year))  # Ordina gli anni in modo crescente
    presence_absence_matrix <- site_subset[, -(1:2)]  
    if (length(years) > 1) {
      # Calcola la beta diversity tra l'anno meno recente e tutti gli altri anni
      for (j in 2:length(years)) {
        # Seleziona i dati dell'anno base e dell'anno di confronto
        base_year_data <- presence_absence_matrix[site_subset$Year == years[1], ]
        comparison_year_data <- presence_absence_matrix[site_subset$Year == years[j], ]
        two_years_data <- rbind(base_year_data[, -(1:3)], comparison_year_data[, -(1:3)])
        # Calcola le componenti della beta-diversity: turnover e nestedness
        beta_result <- beta.pair(two_years_data, index.family = "sorensen")
        # Estrai le componenti turnover e nestedness
       if (is.matrix(beta_result$beta.sim) && is.matrix(beta_result$beta.sne)) {
          turnover_value <- as.numeric(beta_result$beta.sim[1, 2])
          nestedness_value <- as.numeric(beta_result$beta.sne[1, 2])
        } else {
          turnover_value <- as.numeric(beta_result$beta.sim)
          nestedness_value <- as.numeric(beta_result$beta.sne)
        }
        beta_values <- rbind(beta_values, data.frame(Site = site, Year1 = years[1], Year2 = years[j], 
                                                     Turnover = turnover_value, Nestedness = nestedness_value))
      }
    }
  }
  return(beta_values)
}

## FUNZIONE per calcolare Turnover e Nestedness IMMEDIATE
calculate_turnover_nestedness_immediate <- function(presence_absence_data) {
  sites <- unique(presence_absence_data$Site)
  beta_values <- data.frame(Site = character(), Year1 = integer(), Year2 = integer(), 
                            Turnover = numeric(), Nestedness = numeric())
  for (site in sites) {  # Loop per ogni sito
    site_subset <- presence_absence_data[presence_absence_data$Site == site, ]
    years <- sort(unique(site_subset$Year))
    presence_absence_matrix <- site_subset[, -(1:2)] 
    if (length(years) > 1) {
      for (j in 1:(length(years) - 1)) {
        current_year_data <- presence_absence_matrix[site_subset$Year == years[j], ]
        next_year_data <- presence_absence_matrix[site_subset$Year == years[j + 1], ]
        two_years_data <- rbind(current_year_data, next_year_data)
        # Calcola le componenti della beta-diversity: turnover e nestedness
        beta_result <- beta.pair(two_years_data, index.family = "sorensen")
        if (is.matrix(beta_result$beta.sim) && is.matrix(beta_result$beta.sne)) {
          turnover_value <- as.numeric(beta_result$beta.sim[1, 2])
          nestedness_value <- as.numeric(beta_result$beta.sne[1, 2])
        } else {
          turnover_value <- as.numeric(beta_result$beta.sim)
          nestedness_value <- as.numeric(beta_result$beta.sne)
        }
        beta_values <- rbind(beta_values, data.frame(Site = site, Year1 = years[j], Year2 = years[j + 1], 
                                                     Turnover = turnover_value, Nestedness = nestedness_value))
      }
    }
  }
  return(beta_values)
}

## FUNZIONE per calcolare la somma di N.Sus per ogni sito e coppia di anni (Year1, Year2)
calculate_n_sus_sum <- function(presence_absence_data, beta_values) {
  n_sus_years <- presence_absence_data %>%
    select(Site, Year, N.Sus)
  n_sus_sums <- n_sus_years %>%
    inner_join(n_sus_years, by = "Site", suffix = c(".Year1", ".Year2")) %>%
    filter(Year.Year1 != Year.Year2) %>%
    group_by(Site, Year1 = Year.Year1, Year2 = Year.Year2) %>%
    summarise(N.Sus_sum = sum(N.Sus.Year1, N.Sus.Year2, na.rm = TRUE), .groups = "drop")
  final_result <- beta_values %>%
    left_join(n_sus_sums, by = c("Site", "Year1", "Year2"))
  return(final_result)
}

#File
Boreal_nsus           <-read.csv        ("dati\\Boreal_NSus.csv"         )
Beech_nsus           <-read.csv         ("dati\\Nemoral_beech_NSus.csv"  )
Oak_nsus           <-read.csv           ("dati\\Nemoral_oak_NSus.csv"    )
Mediterranean_nsus           <-read.csv ("dati\\Mediterranean_NSus.csv"  )

## Applico funzioni per BETA-DIVERSITY CUMULATIVA
risultati_beta_boreal_cumulative        <- calculate_n_sus_sum(Boreal_nsus, calculate_turnover_nestedness_cumulative(presence_absence_matrix_boreal))
risultati_beta_beech_cumulative         <- calculate_n_sus_sum(Beech_nsus,calculate_turnover_nestedness_cumulative(presence_absence_matrix_beech))
risultati_beta_oak_cumulative           <- calculate_n_sus_sum(Oak_nsus,calculate_turnover_nestedness_cumulative(presence_absence_matrix_oak))
risultati_beta_mediterranean_cumulative <- calculate_n_sus_sum(Mediterranean_nsus,calculate_turnover_nestedness_cumulative(presence_absence_matrix_mediterranean))

# Calcola la Distanza Temporale -CUMULATIVE Turnover-
risultati_beta_boreal_cumulative       $temporal_distance <- risultati_beta_boreal_cumulative       $Year2 - risultati_beta_boreal_cumulative       $Year1 
risultati_beta_beech_cumulative        $temporal_distance <- risultati_beta_beech_cumulative        $Year2 - risultati_beta_beech_cumulative        $Year1 
risultati_beta_oak_cumulative          $temporal_distance <- risultati_beta_oak_cumulative          $Year2 - risultati_beta_oak_cumulative          $Year1 
risultati_beta_mediterranean_cumulative$temporal_distance <- risultati_beta_mediterranean_cumulative$Year2 - risultati_beta_mediterranean_cumulative$Year1 

# LMM per vedere se ci sono TREND nel TURNOVER 
Cumulative_TURNOVER_boreal                  <-lmer(Turnover~temporal_distance + (1|Site) + (1|N.Sus_sum), data = risultati_beta_boreal_cumulative          )
Cumulative_TURNOVER_boreal_snus             <-lmer(Turnover~temporal_distance + (1|Site), data = risultati_beta_boreal_cumulative                          )
Cumulative_TURNOVER_beech                   <-lmer(Turnover~temporal_distance + (1|Site), data = risultati_beta_beech_cumulative                           )
Cumulative_TURNOVER_oak                     <-lmer(Turnover~temporal_distance + (1|Site) , data = risultati_beta_oak_cumulative                            )
Cumulative_TURNOVER_mediterranean           <-lmer(Turnover~temporal_distance + (1|Site) + (1|N.Sus_sum), data = risultati_beta_mediterranean_cumulative   )
Cumulative_TURNOVER_mediterranean_snus      <-lmer(Turnover~temporal_distance + (1|Site), data = risultati_beta_mediterranean_cumulative                   )

# LMM per vedere se ci sono TREND nella NESTEDNESS
Cumulative_NESTEDNESS_beech            <-lmer(Nestedness~temporal_distance + (1|Site) + (1|N.Sus_sum), data = risultati_beta_beech_cumulative   )
Cumulative_NESTEDNESS_beech_nsus       <-lmer(Nestedness~temporal_distance + (1|Site), data = risultati_beta_beech_cumulative                   )

Cumulative_NESTEDNESS_oak              <-lmer(Nestedness~temporal_distance + (1|Site) + (1|N.Sus_sum), data = risultati_beta_oak_cumulative     )
Cumulative_NESTEDNESS_oak_nsus         <-lmer(Nestedness~temporal_distance + (1|Site), data = risultati_beta_oak_cumulative                     )

Cumulative_NESTEDNESS_mediterranean         <-lmer(Nestedness~temporal_distance + (1|Site) + (1|N.Sus_sum), data = risultati_beta_mediterranean_cumulative  )
Cumulative_NESTEDNESS_mediterranean_nsus    <-lmer(Nestedness~temporal_distance + (1|Site), data = risultati_beta_mediterranean_cumulative                  )

# Cumulative_NESTEDNESS_boreal LMM non rispetta le assunzioni, quindi uso GLMM:
epsilon <- 1e-5  # Aggiungo piccola costante
risultati_beta_boreal_cumulative$Nestedness_adjusted <- 
    ifelse(risultati_beta_boreal_cumulative$Nestedness == 0, 
           epsilon, 
           risultati_beta_boreal_cumulative$Nestedness)
risultati_beta_boreal_cumulative$Nestedness_adjusted <- 
    ifelse(risultati_beta_boreal_cumulative$Nestedness_adjusted == 1, 
           1 - epsilon, 
           risultati_beta_boreal_cumulative$Nestedness_adjusted)

Cumulative_NESTEDNESS_boreal_modello_beta_sqrt <- glmmTMB(sqrt(Nestedness_adjusted) ~ temporal_distance + (1 | Site),
                         family = beta_family(link = "logit"), 
                         data = risultati_beta_boreal_cumulative)
summary(Cumulative_NESTEDNESS_boreal_modello_beta_sqrt)



# Applico funzioni per BETA-DIVERSITY IMMEDIATA
risultati_beta_boreal_immediate         <- calculate_n_sus_sum(Boreal_nsus, calculate_turnover_nestedness_immediate(presence_absence_matrix_boreal))
risultati_beta_beech_immediate         <- calculate_n_sus_sum(Beech_nsus,calculate_turnover_nestedness_immediate(presence_absence_matrix_beech))
risultati_beta_oak_immediate           <- calculate_n_sus_sum(Oak_nsus,calculate_turnover_nestedness_immediate(presence_absence_matrix_oak))
risultati_beta_mediterranean_immediate <- calculate_n_sus_sum(Mediterranean_nsus,calculate_turnover_nestedness_immediate(presence_absence_matrix_mediterranean))

# Calcola la Distanza Temporale -IMMEDIATE Turnover-
risultati_beta_boreal_immediate       $temporal_distance <- risultati_beta_boreal_immediate       $Year2 - risultati_beta_boreal_immediate       $Year1 
risultati_beta_beech_immediate        $temporal_distance <- risultati_beta_beech_immediate        $Year2 - risultati_beta_beech_immediate        $Year1 
risultati_beta_oak_immediate          $temporal_distance <- risultati_beta_oak_immediate          $Year2 - risultati_beta_oak_immediate          $Year1 
risultati_beta_mediterranean_immediate$temporal_distance <- risultati_beta_mediterranean_immediate$Year2 - risultati_beta_mediterranean_immediate$Year1 

# LMM per vedere se ci sono TREND nel TURNOVER 
Immediate_TURNOVER_boreal                <-lmer(Turnover~Year2 + (1|Site) + (1|temporal_distance)+ (1|N.Sus_sum), data = risultati_beta_boreal_immediate    )
Immediate_TURNOVER_boreal_nsus_sqrt      <-lmer(sqrt(Turnover)~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_boreal_immediate             )

Immediate_TURNOVER_beech                  <-lmer(Turnover~Year2 + (1|Site) + (1|temporal_distance) + (1|N.Sus_sum), data = risultati_beta_beech_immediate   )
Immediate_TURNOVER_beech_nsus             <-lmer(Turnover~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_beech_immediate                   )

Immediate_TURNOVER_oak_nsus              <-lmer(Turnover~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_oak_immediate)

Immediate_TURNOVER_mediterranean_nsus    <-lmer(Turnover~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_mediterranean_immediate)

# LMM per vedere se ci sono TREND nella NESTEDNESS - Immediate-
Immediate_NESTEDNESS_boreal             <-lmer(Nestedness~Year2 + (1|Site) + (1|temporal_distance) + (1|N.Sus_sum), data = risultati_beta_boreal_immediate  )
Immediate_NESTEDNESS_boreal_nsus_sqrt   <-lmer(sqrt(Nestedness)~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_boreal_immediate            )

Immediate_NESTEDNESS_beech            <-lmer(Nestedness~Year2 + (1|Site) + (1|temporal_distance) + (1|N.Sus_sum), data = risultati_beta_beech_immediate )
Immediate_NESTEDNESS_beech_nsus       <-lmer(Nestedness~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_beech_immediate                 )

Immediate_NESTEDNESS_oak_nsus         <-lmer(Nestedness~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_oak_immediate          )
Immediate_NESTEDNESS_oak_nsus_sqrt    <-lmer(sqrt(Nestedness)~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_oak_immediate    )

Immediate_NESTEDNESS_mediterranean_nsus    <-lmer(Nestedness~Year2 + (1|Site) + (1|temporal_distance), data = risultati_beta_mediterranean_immediate)


# CONTROLLO rispetto delle ASSUNZIONI dei modelli
modello<-
residui<- residuals(modello)
valori_predetti <- fitted(modello)
# 2. Grafico dei residui vs valori predetti
plot(fitted(modello), residuals(modello))
# 2.
ggplot(data = NULL, aes(x = valori_predetti, y = residui)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # Linea di tendenza
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()

# 3. Verifica dell'indipendenza degli errori:
plot(residuals(modello))
# 4. Verifica della normalità degli errori:
qqnorm(residuals(modello),main = "QQ Plot of Residuals") # QQ plot
qqline(residuals(modello), col = "red")
# 4.
hist(residui, probability = TRUE,
     main = "Histogram of Residuals with Normal Curve",
     xlab = "Residuals", ylab = "Density")
curve(dnorm(x, mean = mean(residui), sd = sd(residui)),
      add = TRUE, col = "red", lwd = 2)
densita_empirica <- density(residui)
lines(densita_empirica, col = "blue", lwd = 2)
