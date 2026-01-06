library(dplyr)
library(metafor)

####
## Functions
####
meta_analysis <- 
  function(file2017, file2018, p_value, effect_size, comparison){
    
    deg2017 <- read.table(file2017, sep = "\t", header = TRUE)
    deg2018 <- read.table(file2018, sep = "\t", header = TRUE)
    
    deg2017 <- dplyr::filter(deg2017, geneID %in% deg2018$geneID)
    deg2018 <- dplyr::filter(deg2018, geneID %in% deg2017$geneID)
    
    deg2017$year <- 2017
    deg2018$year <- 2018
    
    meta <- bind_rows(deg2017, deg2018) %>%
      group_by(geneID) %>%
      group_modify(~{
        dat <- .x
        k   <- nrow(dat)
        
        # Study-specific effects
        eff_2017 <- if (any(dat$year == 2017)) dat$log2FoldChange[dat$year == 2017][1] else NA_real_
        se_2017  <- if (any(dat$year == 2017)) dat$lfcSE[dat$year == 2017][1]          else NA_real_
        eff_2018 <- if (any(dat$year == 2018)) dat$log2FoldChange[dat$year == 2018][1] else NA_real_
        se_2018  <- if (any(dat$year == 2018)) dat$lfcSE[dat$year == 2018][1]          else NA_real_
        
        # If we don't have at least 2 studies or SEs are missing, return NA for meta outputs
        if (k < 2 || all(is.na(dat$log2FoldChange)) || all(is.na(dat$lfcSE))) {
          return(tibble(
            k = k,
            eff_2017 = eff_2017,
            se_2017  = se_2017,
            eff_2018 = eff_2018,
            se_2018  = se_2018,
            eff_meta_fixed   = NA_real_,
            se_meta_fixed    = NA_real_,
            z_fixed          = NA_real_,
            p_fixed          = NA_real_,
            eff_meta_random  = NA_real_,
            se_meta_random   = NA_real_,
            z_random         = NA_real_,
            p_random         = NA_real_,
            tau2             = NA_real_,
            I2               = NA_real_
          ))
        }
        
        # Fixed-effect meta-analysis
        fe <- tryCatch(
          metafor::rma.uni(yi = dat$log2FoldChange,
                           sei = dat$lfcSE,
                           method = "FE"),
          error = function(e) NULL
        )
        
        # Random-effects meta-analysis (REML)
        re <- tryCatch(
          metafor::rma.uni(yi = dat$log2FoldChange,
                           sei = dat$lfcSE,
                           method = "REML"),
          error = function(e) NULL
        )
        
        tibble(
          k = k,
          eff_2017 = eff_2017,
          se_2017  = se_2017,
          eff_2018 = eff_2018,
          se_2018  = se_2018,
          
          # Fixed effect
          eff_meta_fixed  = if (!is.null(fe)) as.numeric(fe$b[1]) else NA_real_,
          se_meta_fixed   = if (!is.null(fe)) fe$se[1]           else NA_real_,
          z_fixed         = if (!is.null(fe)) fe$zval[1]         else NA_real_,
          p_fixed         = if (!is.null(fe)) fe$pval[1]         else NA_real_,
          
          # Random effects (REML)
          eff_meta_random = if (!is.null(re)) as.numeric(re$b[1]) else NA_real_,
          se_meta_random  = if (!is.null(re)) re$se[1]            else NA_real_,
          z_random        = if (!is.null(re)) re$zval[1]          else NA_real_,
          p_random        = if (!is.null(re)) re$pval[1]          else NA_real_,
          tau2            = if (!is.null(re)) re$tau2[1]          else NA_real_,
          I2              = if (!is.null(re)) re$I2[1]            else NA_real_  # metafor gives % I²
        )
      }) %>%
      ungroup() %>%
      mutate(
        p_adj_fixed  = p.adjust(p_fixed,  method = "BH"),
        p_adj_random = p.adjust(p_random, method = "BH")
      ) %>%
      arrange(p_random)
    
    write.csv(meta, paste0("./output/processed_data/degs/meta/", comparison, "all_genes.csv"))
    
    signif_genes <- meta %>% 
      filter(p_adj_random < p_value, 
             (eff_meta_random > effect_size | eff_meta_random < -effect_size) )
    
    signif_genes$up_down <- ifelse(signif_genes$eff_meta_random > effect_size, "Up", "Down")
    
    write.csv(signif_genes, paste0("./output/processed_data/degs/meta/", comparison, "signif_genes.csv"))
    
    return(signif_genes)
  }

get_unique_degs <-
  function(nonresponder_degs, responder_degs){
    up_degs_nonresponder <- dplyr::filter(nonresponder_degs, up_down == "Up")
    up_degs_responder <- dplyr::filter(responder_degs, up_down == "Up")  
    up_degs_responder_unique <- dplyr::filter(up_degs_responder, !(geneID %in% up_degs_nonresponder$geneID))
    up_degs_non_responder_unique <- dplyr::filter(up_degs_nonresponder, !(geneID %in% up_degs_responder$geneID))
    
    down_degs_nonresponder <- dplyr::filter(nonresponder_degs, up_down == "Down")
    down_degs_responder <- dplyr::filter(responder_degs, up_down == "Down")  
    down_degs_responder_unique <- dplyr::filter(down_degs_responder, !(geneID %in% down_degs_nonresponder$geneID))
    down_degs_non_responder_unique <- dplyr::filter(down_degs_nonresponder, !(geneID %in% down_degs_responder$geneID))
    
    print(paste0("Number of unique non-responder Up-DEGS: ", nrow(up_degs_non_responder_unique)))
    print(paste0("Number of unique non-responder Down-DEGS: ", nrow(down_degs_non_responder_unique)))
    print(paste0("Number of unique responder Up-DEGS: ", nrow(up_degs_responder_unique)))
    print(paste0("Number of unique responder Down-DEGS: ", nrow(down_degs_responder_unique)))
    
    list(
      unique_nonresponders = bind_rows(up_degs_non_responder_unique, down_degs_non_responder_unique),
      unique_responders = bind_rows(up_degs_responder_unique, down_degs_responder_unique)
    )
  }


####
## Combine study years with REML random-effect meta-analysis
####

##Blood HAI
blood_hai <-
  meta_analysis("./output/processed_data/degs/Blood2017_HIGH-HR_4x_max_gmfr_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_HIGH-HR_4x_max_gmfr_DESEq2.tsv",
                0.1, 
                0.322,
                "blood_hai")

blood_nohai <-
  meta_analysis("./output/processed_data/degs/Blood2017_LOW-HR_4x_max_gmfr_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_LOW-HR_4x_max_gmfr_DESEq2.tsv",
                0.1,
                0.322,
                "blood_nohai")

##Blood IgA
blood_iga <-
  meta_analysis("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_iga_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_iga_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_iga")

blood_noiga <-
  meta_analysis("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_iga_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_iga_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_noiga")

##CD4
blood_cd4 <-
  meta_analysis("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_cd4")

blood_nocd4 <-
  meta_analysis("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_nocd4")

## CD8
blood_cd8 <-
  meta_analysis("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_cd8")

blood_nocd8 <-
  meta_analysis("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                0.1,
                0.322,
                "blood_nocd8")

## shedding
blood_shed <-
  meta_analysis("./output/processed_data/degs/Blood2017_HIGH-sheddD2_2orMore_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_HIGH-sheddD2_2orMore_DESEq2.tsv",
                0.1,
                0.322,
                "blood_shed")

blood_noshed <-
  meta_analysis("./output/processed_data/degs/Blood2017_LOW-sheddD2_2orMore_DESEq2.tsv",
                "./output/processed_data/degs/Blood2018_LOW-sheddD2_2orMore_DESEq2.tsv",
                0.1,
                0.322,
                "blood_noshed")

####
## Unique DEGS
####

##Blood
blood_hai_degs <- get_unique_degs(blood_nohai, blood_hai)
saveRDS(blood_hai_degs, file = "./output/processed_data/degs/unique/blood_hai_degs.rds")
blood_iga_degs <- get_unique_degs(blood_noiga, blood_iga)
saveRDS(blood_iga_degs, file = "./output/processed_data/degs/unique/blood_iga_degs.rds")
blood_cd4_degs <- get_unique_degs(blood_nocd4, blood_cd4)
saveRDS(blood_cd4_degs, file = "./output/processed_data/degs/unique/blood_cd4_degs.rds")
blood_cd8_degs <- get_unique_degs(blood_nocd8, blood_cd8)
saveRDS(blood_cd8_degs, file = "./output/processed_data/degs/unique/blood_cd8_degs.rds")
blood_shed_degs <- get_unique_degs(blood_noshed, blood_shed)
saveRDS(blood_shed_degs, file = "./output/processed_data/degs/unique/blood_shed_degs.rds")
