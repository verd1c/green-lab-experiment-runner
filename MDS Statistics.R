library(dplyr)
library(ggplot2)
library(e1071)

data <- read.csv("~/run_table_full.csv")
data$force_field <- as.factor(data$force_field)
data$integrator <- as.factor(data$integrator)
data$nsteps <- as.factor(data$nsteps)


get_normalised_inverse <- function(a) {
  return(1 - ((a - min(a)) / (max(a) - min(a))))
}

data$Accuracy <- get_normalised_inverse(data$avg_rmsd)
data$Efficiency <- get_normalised_inverse(data$energy_usage)


get_harmonic_mean <- function(a, e) { 
  return(2 * a * e / (a + e))
}

data$Energy_Accuracy <- mapply(get_harmonic_mean, data$Accuracy, data$Efficiency)


############### SUMMARY STATISTICS 
get_density_at <- function(x, m) {
  den <- density(x)
  approx(den$x, den$y, xout = m)$y
}

data_force_field <- data %>%
group_by(force_field) %>%
summarise(
  Mean = mean(Energy_Accuracy),
  Median = median(Energy_Accuracy),
  Skewness = skewness(Energy_Accuracy),
  Density_at_Mean = get_density_at(Energy_Accuracy, Mean),
  Density_at_Median = get_density_at(Energy_Accuracy, Median),
  Standard_Deviation = sd(Energy_Accuracy),
  Variance = var(Energy_Accuracy),
  Range = diff(range(Energy_Accuracy)),
  Minimum = min(Energy_Accuracy),
  Maximum = max(Energy_Accuracy),
  Coefficient_of_Variation = 100 * (Standard_Deviation / Mean)
)

data_integrator <- data %>%
group_by(integrator) %>%
summarise(
  Mean = mean(Energy_Accuracy),
  Median = median(Energy_Accuracy),
  Density_at_Mean = get_density_at(Energy_Accuracy, Mean),
  Density_at_Median = get_density_at(Energy_Accuracy, Median),
  Skewness = skewness(Energy_Accuracy),
  Standard_Deviation = sd(Energy_Accuracy),
  Variance = var(Energy_Accuracy),
  Range = diff(range(Energy_Accuracy)),
  Minimum = min(Energy_Accuracy),
  Maximum = max(Energy_Accuracy),
  Coefficient_of_Variation = 100 * (Standard_Deviation / Mean)
)

data_nsteps <- data %>%
group_by(nsteps) %>%
summarise(
  Mean = mean(Energy_Accuracy),
  Median = median(Energy_Accuracy),
  Skewness = skewness(Energy_Accuracy),
  Density_at_Mean = get_density_at(Energy_Accuracy, Mean),
  Density_at_Median = get_density_at(Energy_Accuracy, Median),
  Standard_Deviation = sd(Energy_Accuracy),
  Variance = var(Energy_Accuracy),
  Range = diff(range(Energy_Accuracy)),
  Minimum = min(Energy_Accuracy),
  Maximum = max(Energy_Accuracy),
  Coefficient_of_Variation = 100 * (Standard_Deviation / Mean)
)

data <- data %>%
mutate(Configuration = interaction(integrator, force_field, nsteps))
data$Configuration <- as.factor(data$Configuration)
  
data_configuration <- data %>%
group_by(Configuration) %>%
summarise(
  Mean = mean(Energy_Accuracy),
  Median = median(Energy_Accuracy),
  Skewness = skewness(Energy_Accuracy),
  Density_at_Mean = get_density_at(Energy_Accuracy, Mean),
  Density_at_Median = get_density_at(Energy_Accuracy, Median),
  Standard_Deviation = sd(Energy_Accuracy),
  Variance = var(Energy_Accuracy),
  Range = diff(range(Energy_Accuracy)),
  Minimum = min(Energy_Accuracy),
  Maximum = max(Energy_Accuracy),
  Coefficient_of_Variation = 100 * (Standard_Deviation / Mean)
)

############### BOX PLOTS 
ggplot(data, aes(x = Energy_Accuracy, y = force_field)) +
geom_boxplot(aes(fill = force_field), alpha = 0.5) +
scale_y_discrete(labels = c("charmm" = "CHARMM27", "opls" = "OPLS-AA/L")) +
stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5, stroke = 1, fill = "yellow") +
labs(title = "", x = "Energy-Accuracy", y = "") +
theme_minimal() +
theme(legend.position = "none")

ggplot(data, aes(x = Energy_Accuracy, y = integrator)) +
geom_boxplot(aes(fill = integrator), alpha = 0.5) +
stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5, stroke = 1, fill = "yellow") +
labs(title = "", x = "Energy-Accuracy", y = "") +
theme_minimal()+
theme(legend.position = "none")

ggplot(data, aes(x = Energy_Accuracy, y = nsteps)) +
geom_boxplot(aes(fill = nsteps), alpha = 0.5) +
stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5, stroke = 1, fill = "yellow") +
labs(title = "", x = "Energy-Accuracy", y = "") +
theme_minimal()+
theme(legend.position = "none")


############### SCATTER PLOTS
ggplot(data, aes(x = Efficiency, y = Accuracy)) +
geom_point(aes(color = force_field), alpha = 0.5, size = 2) +
scale_color_discrete(labels = c("charmm" = "CHARMM27", "opls" = "OPLS-AA/L")) +
geom_smooth(method = "lm", color = "black", se = FALSE) +
labs(title = "", x = "Energy Efficiency", y = "Accuracy", color = "Force Field:") +
theme_minimal() + 
theme(legend.position = "bottom")

ggplot(data, aes(x = Efficiency, y = Accuracy)) +
geom_point(aes(color = integrator), alpha = 0.5, size = 2) +
geom_smooth(method = "lm", color = "black", se = FALSE) +
labs(title = "", x = "Energy Efficiency", y = "Accuracy", color = "Integrator:") +
theme_minimal() + 
theme(legend.position = "bottom")

ggplot(data, aes(x = Efficiency, y = Accuracy)) +
geom_point(aes(color = nsteps), alpha = 0.5, size = 2) +
geom_smooth(method = "lm", color = "black", se = FALSE) +
labs(title = "", x = "Energy Efficiency", y = "Accuracy", color = "Step Count:") +
theme_minimal() + 
theme(legend.position = "bottom")


############### DENSITY PLOTS
config_labels = c("md.charmm.50000" = "CHARMM27, md, 50000", 
                  "sd.charmm.50000" = "CHARMM27, sd, 50000", 
                  "md.opls.50000" = "OPLS-AA/L, md, 50000", 
                  "sd.opls.50000" = "OPLS-AA/L, sd, 50000", 
                  "md.charmm.250000" = "CHARMM27, md, 250000", 
                  "sd.charmm.250000" = "CHARMM27, sd, 250000", 
                  "md.opls.250000" = "OPLS-AA/L, md, 250000", 
                  "sd.opls.250000" = "OPLS-AA/L, sd, 250000", 
                  "md.charmm.500000" = "CHARMM27, md, 500000", 
                  "sd.charmm.500000" = "CHARMM27, sd, 500000", 
                  "md.opls.500000" = "OPLS-AA/L, md, 500000", 
                  "sd.opls.500000" = "OPLS-AA/L, sd, 500000")

ggplot(data, aes(x = Energy_Accuracy)) +
geom_density(aes(fill = Configuration), alpha = 0.5, color = "black", bw = "nrd0") +
scale_fill_discrete(labels = config_labels) +
facet_wrap(~ Configuration, labeller = labeller(Configuration = config_labels)) +
labs(title = "",x = expression(bold("Energy-Accuracy")), y = expression(bold("Density")), fill = "Configuration") +
theme_minimal() +
theme(legend.position = "none")

ggplot(data, aes(x = Energy_Accuracy, fill = force_field)) +
geom_density(alpha = 0.5) +
scale_color_manual(values = c("Mean" = "black", "Median" = "black")) +
geom_point(data = data_force_field, aes(x = Mean, y = Density_at_Mean, color = "Mean"), size = 3.5, stroke = 1, fill = "yellow", shape = 21) +
geom_point(data = data_force_field, aes(x = Median, y = Density_at_Median, color = "Median"), size = 3.5, stroke = 1, fill = "orange", shape = 22) +
scale_fill_discrete(labels = c("charmm" = "CHARMM27", "opls" = "OPLS-AA/L")) +
labs(title = "", x = "Energy-Accuracy", y = "Density", fill = "Force Field:", color = "Properties:") +
theme_minimal() + 
theme(legend.position = "bottom")

ggplot(data, aes(x = Energy_Accuracy, fill = integrator)) +
geom_density(alpha = 0.5) +
scale_color_manual(values = c("Mean" = "black", "Median" = "black")) +
geom_point(data = data_integrator, aes(x = Mean, y = Density_at_Mean, color = "Mean"), size = 3.5, stroke = 1, fill = "yellow", shape = 21) +
geom_point(data = data_integrator, aes(x = Median, y = Density_at_Median, color = "Median"), size = 3.5, stroke = 1, fill = "orange", shape = 22) +
guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
labs(title = "", x = "Energy-Accuracy", y = "Density", fill = "Integrator:", color = "Properties:") +
theme_minimal() + 
theme(legend.position = "bottom")

ggplot(data, aes(x = Energy_Accuracy, fill = nsteps)) +
geom_density(alpha = 0.5) +
scale_color_manual(values = c("Mean" = "black", "Median" = "black")) +
geom_point(data = data_nsteps, aes(x = Mean, y = Density_at_Mean, color = "Mean"), size = 3.5, stroke = 1, fill = "yellow", shape = 21) +
geom_point(data = data_nsteps, aes(x = Median, y = Density_at_Median, color = "Median"), size = 3.5, stroke = 1, fill = "orange", shape = 22) +
labs(title = "", x = "Energy-Accuracy", y = "Density", fill = "Step Count:", color = "Properties:") +
theme_minimal() + 
theme(legend.position = "bottom")


############### HYPOTHESIS TESTING
ANOVA_SQ <- aov(Energy_Accuracy ~ force_field + integrator + nsteps, data = data)
summary(ANOVA_SQ)
Shapiro_SQ <- shapiro.test(residuals(ANOVA_SQ))
Tukey_SQ <- TukeyHSD(ANOVA_SQ)

ANOVA_RQ <- aov(Energy_Accuracy ~ force_field * integrator * nsteps, data = data)
summary(ANOVA_RQ)
Shapiro_RQ <- shapiro.test(residuals(ANOVA_RQ))
Tukey_RQ <- TukeyHSD(ANOVA_RQ)
