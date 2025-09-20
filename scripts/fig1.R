library(dplyr)
library(tidyr)
library(ggimage)
library(ape)
library(ggtree)
library(treeio)
library(tidytree)
library(phylogram)
library(ggnewscale)
library(stringr)
library(svglite)

# import data
metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")

####
# Fig 1A - Study design made in biorender
####

####
# Figure 1B - Circle plot
####

# make tree structure
responses <- 
    metadata %>%
    select(HR_4x_max_gmfr, HR_2x_max_iga_fc, HR_2x_max_mnp_cd4_fc, HR_2x_max_mnp_cd8_fc) %>%
    drop_na() 

names(responses) <- c("HAI (max)", "IgA (max)", "CD4 T-cell", "CD8 T-cell")

# make arbritary variable to force order of single responess
responses$n_responses <- rowSums(responses)
responses$grouping <- ifelse(responses$n_responses == 0, 0, responses$n_responses) 
responses$grouping <- ifelse(responses$n_responses == 1 & responses$`HAI (max)` == 1 , 0.4, responses$n_responses)
responses$grouping <- ifelse(responses$n_responses == 1 & responses$`IgA (max)` == 1 , 0.5, responses$grouping)
responses$grouping <- ifelse(responses$n_responses == 1 & responses$`CD4 T-cell` == 1 , 0.6, responses$grouping)
responses$grouping <- ifelse(responses$n_responses == 1 & responses$`CD8 T-cell` == 1 , 0.7, responses$grouping)

# make tree
hc <- hclust(dist(responses))
tree <- as.phylo(hc)

# reorder clades
circ <- ggtree(tree, layout = "fan", open.angle = 30, colour = "white")
circ <- rotate(circ, node = 190)
circ <- rotate(circ, node = 193)
circ <- rotate(circ, node = 199)
circ <- rotate(circ, node = 200)
circ <- rotate(circ, node = 202)
circ <- rotate_tree(circ, angle = 85)


gheatmap(circ, responses[,c(1,2,3,4)], offset = 0, width = 0.5, 
         colnames_angle = 0, colnames_offset_y = 0, hjust = 0, font.size = 2.5,) +
    scale_fill_gradient(low = "lightgoldenrodyellow", high = "gray30", name = "Response") +
    theme(legend.position = "none", 
          axis.line = element_blank(),  
          axis.ticks = element_blank(),  
          axis.text = element_blank())  +
    
    geom_cladelab(node = 199, label = "0", offset = 1.75,  offset.text = 0.4, align = TRUE, barsize = 1, textangle = 0, barcolour = "lightblue") +
    geom_cladelab(node = 200, label = "1", offset = 1.75 , offset.text = 0.45, align = TRUE, barsize = 1, textangle = 0, barcolour = "peachpuff") +
    geom_cladelab(node = 194, label = "2", offset = 1.75 , offset.text = 0.35, align = TRUE, barsize = 1, textangle = 0, barcolour = "lightpink") +
    geom_cladelab(node = 198, label = "3", offset = 1.75 , offset.text = 0.2, align = TRUE, barsize = 1, textangle = 0, barcolour = "lightgreen") +
    geom_cladelab(node = 197, label = "4", offset = 1.75 , offset.text = 0.25, align = TRUE, barsize = 1, textangle = 0, barcolour = "lavender") 

ggsave("./output/figs/fig1a.svg", height = 4)

ggplot() + 
    geom_tile(aes(x = 5.1, y = 2, width = 0.3, height = 0.6), fill = "gray30") +
    geom_tile(aes(x = 5.1, y = 1, width = 0.3, height = 0.6), fill = "lightgoldenrodyellow") +
    annotate("text", x = 4.3, y = 2, label = "Responder", size = 8) +
    annotate("text", x = 4.44, y = 1, label = "Non-responder", size = 8) +
    coord_cartesian(clip = "off") +
    theme_nothing() +
    theme(plot.margin = unit(c(0.5, 1, 0.5, 2.5), "cm"))

ggsave("./output/figs/fig1b_legend.svg", width = 4)

####
## Figure 1C - Response barchart
####
## reload data to bring back missing data
metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")

responses <- 
  metadata %>%
  select(HR_4x_max_gmfr, HR_2x_max_iga_fc, HR_2x_max_mnp_cd4_fc, HR_2x_max_mnp_cd8_fc)

names(responses) <- c("HAI (max)", "IgA (max)", "CD4 T-cell", "CD8 T-cell")


response_proportions <- 
  responses %>%
  select(
    "HAI (max)",   
    "IgA (max)",   
    "CD4 T-cell",  
    "CD8 T-cell"
  ) %>%
  pivot_longer(c("HAI (max)",   
                 "IgA (max)",   
                 "CD4 T-cell",  
                 "CD8 T-cell")) %>%
  filter(!is.na(value)) %>%
  group_by(name) %>%
  summarise(
    n = sum(value),
    total = n()
  ) %>%
  rowwise() %>%
  mutate(
    percent = n*100/total,
    lower_ci = binom.test(n, total)$conf.int[1]*100,
    upper_ci = binom.test(n, total)$conf.int[2]*100
  ) 

response_proportions$name <- factor(response_proportions$name, levels = c("CD4 T-cell", "HAI (max)", "IgA (max)", "CD8 T-cell"))

response_proportions %>%
  ggplot(aes(x = name, y = percent)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, width = 0.05)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  xlab("Immune response") +
  ylab("% response") +
  theme_bw() +
  theme(plot.margin = margin(t = 20, r = 20, b = 10, l = 10),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=15, face = "bold"))

ggsave("./output/figs/fig1c.svg", width = 5.25, height = 3.75)  


####
## Fig 1D - Pairwise chi2 forest plot
####

pairs <- 
  expand.grid(response1 = names(responses)[1:4], response2 = names(responses)[1:4]) %>%
  mutate(response1 = as.character(response1),
         response2 = as.character(response2))

pairs <- pairs[c(2,3,4,7,8,12),]


pairwise_chisq <- 
  lapply(1:nrow(pairs), 
         function(row_n){
           two_by_two <- table(responses[[pairs$response1[row_n]]], responses[[pairs$response2[row_n]]])
           odds_ratio <- fisher.test(two_by_two)$estimate
           or_ci <- fisher.test(two_by_two)$conf.int
           chi_sq <- chisq.test(two_by_two)
           tibble(
             response1 = pairs$response1[row_n],
             response2 = pairs$response2[row_n],
             odds_ratio = odds_ratio,
             lower_95ci = or_ci[1],
             upper_95ci = or_ci[2],
             p_value = chi_sq$p.value
           ) 
         }) %>% 
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "holm")) %>%
  arrange(desc(p_adj))

pairwise_chisq$response1 <- str_replace(pairwise_chisq$response1, "\\(max\\)", "")
pairwise_chisq$response2 <- str_replace(pairwise_chisq$response2, "\\(max\\)", "")
pairwise_chisq$response1 <- str_replace(pairwise_chisq$response1, " T-cell", "")
pairwise_chisq$response2 <- str_replace(pairwise_chisq$response2, " T-cell", "")

pairwise_chisq$comparison <- paste0(pairwise_chisq$response1, " & ", pairwise_chisq$response2)
pairwise_chisq <- pairwise_chisq[order(pairwise_chisq$odds_ratio, decreasing = T),]
pairwise_chisq$comparison <- factor(pairwise_chisq$comparison, levels = pairwise_chisq$comparison)

pairwise_chisq %>%
  ggplot(aes(x = odds_ratio, y = comparison)) +
  geom_tile(width = 0.045, height = 0.1) +
  geom_errorbarh(aes(xmin = lower_95ci, xmax = upper_95ci, height = 0.1)) +
  geom_vline(xintercept = 1, linetype = 2) +
  annotate(geom="text", x = 3, y=1.25, label="Adjusted P value = 0.01", size = 3) +
  xlab("Odds ratio") + 
  ylab("Paired response") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10, face = "bold"),
        axis.text.x = element_text(size=12, face = "bold"),
        axis.title = element_text(size=12, face = "bold")
  )

ggsave("./output/figs/fig1d.svg")
