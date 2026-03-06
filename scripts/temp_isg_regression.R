library(dplyr)
library(car)
library(ResourceSelection)
library(ggplot2)
library(effects)
library(pROC)
library(ggtext)


#select and clean metadata
metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")
metadata <- 
    metadata %>% select(subid1, age, sex, year, HR_4x_max_gmfr, SheddingD2Sum, contains("v0_gmt"))

metadata$n_prior_exposure <- rowSums(metadata[, c("h1_v0_gmt", "h3_v0_gmt", "b_v0_gmt")] > 5)
metadata <- metadata %>% select(!contains("v0_gmt"))

#merge and filter MX1 fold change
nasal_2017 <- read.table("./output/processed_data/Nasal_2017_FiltredToSharedGenes_log2FC_Ranks.tsv", sep = "\t")
nasal_2018 <- read.table("./output/processed_data/Nasal_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
blood_2017 <- read.table("./output/processed_data/Blood_2017FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
blood_2018 <- read.table("./output/processed_data/Blood_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")

nasal_2017 <- nasal_2017["MX1",]
nasal_2018 <- nasal_2018["MX1",]
nasal <- as.data.frame(t(bind_cols(nasal_2017, nasal_2018)))
nasal$subid1 <- rownames(nasal)
rownames(nasal) <- NULL
nasal <- nasal %>% dplyr::rename(mx1_nose = MX1)

blood_2017 <- blood_2017["MX1",]
blood_2018 <- blood_2018["MX1",]
blood <- as.data.frame(t(bind_cols(blood_2017, blood_2018)))
blood$subid1 <- rownames(blood)
rownames(blood) <- NULL
blood <- blood %>% dplyr::rename(mx1_blood = MX1)

#merge datasets
metadata <- left_join(metadata, nasal)
metadata <- left_join(metadata, blood)

#logistic regression
metadata$HR_4x_max_gmfr <- as.factor(metadata$HR_4x_max_gmfr)
metadata$mx1_nose_z  <- scale(metadata$mx1_nose,  center = TRUE, scale = TRUE)
metadata$mx1_blood_z <- scale(metadata$mx1_blood, center = TRUE, scale = TRUE)


# Blood MX1 model
# Fit logistic regression
blood_metadata <- filter(metadata, !is.na(metadata$mx1_blood))
model <- glm(HR_4x_max_gmfr ~ age + sex + year + SheddingD2Sum + 
               n_prior_exposure + mx1_blood_z,
             data = blood_metadata, family = binomial)
summary(model)

#check colinarity
vif(model) #Less than 2 is good

# Cook's distance
plot(model, which = 4)
blood_metadata[c(17, 31, 37),]
#rerun without these three observations
model_excluded <- glm(HR_4x_max_gmfr ~ age + sex + year + SheddingD2Sum + 
               n_prior_exposure + mx1_blood_z,
             data = blood_metadata[-c(17, 31, 37),], family = binomial)
summary(model_excluded)
#P values get more significant when excluding influential observations
# This doesn't change the interpretation so keep original model

#Hosmer-Lemeshow goodness-of-fit test
hoslem.test(as.numeric(blood_metadata$HR_4x_max_gmfr) - 1, fitted(model))

# Predicted probabilities
blood_metadata$pred_prob <- predict(model, type = "response")

# Plot predicted probabilities vs actual outcomes
ggplot(blood_metadata, aes(x = pred_prob, fill = HR_4x_max_gmfr)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  labs(title = "Predicted probabilities by actual outcome",
       x = "Predicted Probability", y = "Count")                       

#marginal effects
plot(allEffects(model))

#roc curve
roc_obj <- roc(blood_metadata$HR_4x_max_gmfr, blood_metadata$pred_prob)
plot(roc_obj, print.auc = TRUE)

# Make forrest plot
# Extract coefficients and CIs
coefs <- summary(model)$coefficients
conf <- confint(model)  # 95% CI

# Create dataframe for plotting
plot_df <- data.frame(
  variable = rownames(coefs),
  OR = exp(coefs[, "Estimate"]),
  lower = exp(conf[, 1]),
  upper = exp(conf[, 2]),
  p = coefs[, "Pr(>|z|)"]
)

plot_df$label <- dplyr::recode(plot_df$variable,
                        "sexM" = "Sex (Male)",
                        "mx1_blood_z" = "*MX1* Fold change",  # <- italic markdown
                        "year" = "Year",
                        "SheddingD2Sum" = "Number Strains Shed",
                        "age" = "Age",
                        "n_prior_exposure" = "Number Prior Strains"
)

plot_df$significant <- plot_df$p < 0.05

ggplot(plot_df[-1, ], aes(x = reorder(label, OR), y = OR, colour = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  scale_y_log10() +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Blood",
       x = "", y = "Odds Ratio (log scale)") +
  theme_minimal() +
  theme(axis.text.y = element_markdown(),
        legend.position = "none") 

ggsave("./output/figs/fig5a.svg")

# Nose MX1 model
# Fit logistic regression
nose_metadata <- filter(metadata, !is.na(metadata$mx1_nose))
model <- glm(HR_4x_max_gmfr ~ age + sex + year + SheddingD2Sum + 
               n_prior_exposure + mx1_nose_z,
             data = nose_metadata, family = binomial)
summary(model)

#check colinarity
vif(model) #Less than 2 is good

# Cook's distance
plot(model, which = 4)
nose_metadata[c(19, 67, 76),]
#rerun without these three observations
model_excluded <- glm(HR_4x_max_gmfr ~ age + sex + year + SheddingD2Sum + 
                        n_prior_exposure + mx1_nose_z,
                      data = nose_metadata[-c( 19, 67, 76),], family = binomial)
summary(model_excluded)
#P values get more significant when excluding influential observations
# This doesn't change the interpretation so keep original model

#Hosmer-Lemeshow goodness-of-fit test
hoslem.test(as.numeric(nose_metadata$HR_4x_max_gmfr) - 1, fitted(model))

# Predicted probabilities
nose_metadata$pred_prob <- predict(model, type = "response")

# Plot predicted probabilities vs actual outcomes
ggplot(nose_metadata, aes(x = pred_prob, fill = HR_4x_max_gmfr)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  labs(title = "Predicted probabilities by actual outcome",
       x = "Predicted Probability", y = "Count")                       

#marginal effects
plot(allEffects(model))

#roc curve
roc_obj <- roc(nose_metadata$HR_4x_max_gmfr, nose_metadata$pred_prob)
plot(roc_obj, print.auc = TRUE)

# Make forrest plot
# Extract coefficients and CIs
coefs <- summary(model)$coefficients
conf <- confint(model)  # 95% CI

# Create dataframe for plotting
plot_df <- data.frame(
  variable = rownames(coefs),
  OR = exp(coefs[, "Estimate"]),
  lower = exp(conf[, 1]),
  upper = exp(conf[, 2]),
  p = coefs[, "Pr(>|z|)"]
)

plot_df$significant <- plot_df$p < 0.05

plot_df$label <- dplyr::recode(plot_df$variable,
                               "sexM" = "Sex (Male)",
                               "mx1_nose_z" = "*MX1* Fold change",  # <- italic markdown
                               "year" = "Year",
                               "SheddingD2Sum" = "Number Strains Shed",
                               "age" = "Age",
                               "n_prior_exposure" = "Number Prior Strains"
)

ggplot(plot_df[-1, ], aes(x = reorder(label, OR), y = OR, colour = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  scale_y_log10() +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Nasopharyngeal ",
       x = "", y = "Odds Ratio (log scale)") +
  theme_minimal() +
  theme(axis.text.y = element_markdown(),
        legend.position = "none") 

ggsave("./output/figs/fig5b.svg")
