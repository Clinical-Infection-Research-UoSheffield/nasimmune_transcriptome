library(BayesPrism)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)

# load and extract data

# Compare default signature matrix generation with signature gene deconvolution
blood_deconvolve <- readRDS("./output/processed_data/blood_bayesprism_sigmat.rds")
fractions <- get.fraction (bp=blood_deconvolve,
                       which.theta="final",
                       state.or.type="type")

cell_type_fraction_sigmat <- as.data.frame(colMeans(fractions))
names(cell_type_fraction_sigmat) <- c("mean")

cv <- blood_deconvolve@posterior.theta_f@theta.cv
cell_type_cv_sigmat <- as.data.frame(colMeans(cv))
names(cell_type_cv_sigmat) <- c("mean")

as.data.frame(fractions) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "proportion"       
  ) %>% 
  mutate(visit0 = str_detect(subid, "V0")) %>%
  ggplot(aes(x = proportion, colour = visit0)) +
  geom_density() +
  facet_wrap(cell_type~., scales = "free",)
  
ggsave("./output/figs/bayesprism_fractions_sigmat.jpeg", width = 7, height = 4)

as.data.frame(cv) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "cv"       
  ) %>% 
  mutate(visit0 = str_detect(subid, "V0")) %>%
  ggplot(aes(x = cv, colour = visit0)) +
  geom_density() +
  facet_wrap(cell_type~., scales = "free",) +
  theme(axis.text.x = element_text(size = 5, angle = 20))

ggsave("./output/figs/bayesprism_cv_sigmat.jpeg", width = 7, height = 4)

as.data.frame(fractions) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "proportion"       
  ) %>%
  left_join(., 
            as.data.frame(cv) %>%
              rownames_to_column(var = "subid") %>%
              pivot_longer(
                cols = -subid,                 
                names_to = "cell_type",        
                values_to = "cv"       
              )) %>%
  ggplot(aes(x = proportion, y = log(cv))) +
  geom_point( size = 0.2) +
  facet_wrap(cell_type~., scales = "free",) 

ggsave("./output/figs/bayesprism_cv_prop_sigmat.jpeg", width = 7, height = 4)

# Gene signature
blood_deconvolve <- readRDS("./output/processed_data/blood_bayesprism_genesig.rds")
fractions <- get.fraction (bp=blood_deconvolve,
                           which.theta="final",
                           state.or.type="type")

cell_type_fraction_genesig <- as.data.frame(colMeans(fractions))
names(cell_type_fraction_genesig) <- c("mean")

cv <- blood_deconvolve@posterior.theta_f@theta.cv
cell_type_cv_genesig <- as.data.frame(colMeans(cv))
names(cell_type_cv_genesig) <- c("mean")

as.data.frame(fractions) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "proportion"       
  ) %>% 
  mutate(visit0 = str_detect(subid, "V0")) %>%
  ggplot(aes(x = proportion, colour = visit0)) +
  geom_density() +
  facet_wrap(cell_type~., scales = "free") +
  theme(axis.text.x = element_text(size = 5, angle = 20))

ggsave("./output/figs/bayesprism_fractions_genesig.jpeg", width = 7, height = 4)

as.data.frame(cv) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "cv"       
  ) %>% 
  mutate(visit0 = str_detect(subid, "V0")) %>%
  ggplot(aes(x = cv, colour = visit0)) +
  geom_density() +
  facet_wrap(cell_type~., scales = "free") +
  theme(axis.text.x = element_text(size = 5, angle = 20))

ggsave("./output/figs/bayesprism_cv_genesig.jpeg", width = 7, height = 4)

as.data.frame(fractions) %>%
  rownames_to_column(var = "subid") %>%
  pivot_longer(
    cols = -subid,                 
    names_to = "cell_type",        
    values_to = "proportion"       
  ) %>%
  left_join(., 
            as.data.frame(cv) %>%
              rownames_to_column(var = "subid") %>%
              pivot_longer(
                cols = -subid,                 
                names_to = "cell_type",        
                values_to = "cv"       
              )) %>%
  ggplot(aes(x = proportion, y = log(cv))) +
  geom_point( size = 0.2) +
  facet_wrap(cell_type~., scales = "free",) 

ggsave("./output/figs/bayesprism_cv_prop_genesig.jpeg", width = 7, height = 4)

# 1. Filter out cell types with low fraction and high cv


# 2. Describe average composition at day 0 and day 2.

# 3. Look at cell type specific gene expression. Cluster vst() normalised by day or immune response.