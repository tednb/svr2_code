# marginal conditional model
# regl <- function(est,frac.c,covs,age){ # all vectors except age
#   betav <-sum(est[1:5]*frac.c)+sum(frac.c*age*est[6:10])+sum(est[11:12]*covs)
#   return(betav)
# }

regl <- function(est,age){ # all vectors except age
  betav <- as.numeric(age*est[2]+est[1])
  return(betav)
}
# regl <- function(est,age){ # all vectors except age
#   betav <- as.numeric(age*est[2]+est[1]+est[3])
#   return(betav)
# }
regl_bulk <- function(est,frac.c,covs,age){
  betav <- c()
  for (k in 1:nrow(frac.c)){
    betav[k] <-sum(est[1:5]*frac.c[k,])+sum(frac.c[k,]*age*est[6:10])+sum(est[11:12]*covs[,k])
  }
  return(betav)
}


# binary transform
bt <- function(cell.m,...){
  ## idx
  sc.idx <- apply(cell.m,2,function(x) {which(x>=0.9 | x<=0.1)})
  sc.lv <- list()
  # extract and transform
  for (i in 1:length(sc.idx)){
    cell <- cell.m[sc.idx[[i]],i]
    print(length(cell))
    cell[cell>=0.9]<-1
    cell[cell<=0.1]<-0
    sc.lv[[i]]<-cell
  }
  names(sc.lv) <- colnames(cell.m)
  return(sc.lv)
}

# Scage
compute_probabilities <- function(cell, # vector: one cell
                                  cell.c, # vector effect size
                                  cell.e, # vector: formula coefs
                                  selection_mode = "percentile",
                                  CpG_parameter=0.9,
                                  zero_met_replacement = 0.001,
                                  one_met_replacement = 0.999,
                                  min_age = 0,
                                  max_age = 100,
                                  age_step = 0.1,
                                  uncertainty = 1){
  library(progress)
  # get arguments
  cell <- cell[names(cell.c)]
  single_cell_name <- names(cell)
  single_cell_met <- cell
  start <- Sys.time()
  # judge <- TRUE
  cuts <- 1
  # profiling mode selection
  # effect size sort
  if (selection_mode == "percentile") { # ex: top 1% age-associated CpGs per cell
    
    quantile <- quantile(abs(cell.c), probs = c(cuts,CpG_parameter))
    abs_top <- subset(cell, abs(cell.c) >= quantile[2] & abs(cell.c) <= quantile[1])
    #print(quantile[2])
  } else if (selection_mode == "numCpGs") { # ex: top 500 age-associated CpGs per cell
    abs_top <- cell[order(-abs(cell.c))][1:CpG_parameter]
  }
  print(length(abs_top))
  cell_subset <- cell[names(abs_top)]
  # isolate selected sites
  selected_sites <- names(cell_subset)
  
  # get age steps from min_age to max_age (inclusive of both)
  age_steps <- seq(min_age, max_age + age_step, age_step)
  
  # create list to store probability profiles
  list_of_profile_probabilities_per_age <- list()
  
  # loop through each age step
  pb <- progress_bar$new(total = length(age_steps)) # progress bar
  for (age in age_steps) {
    # create list to hold probability for all chosen CpGs for a given age
    probability_list_one_age <- list()
    # loop through each site
    for (site in selected_sites) {
      # compute methylation probability
      methylation_probability <- regl(cell.e[[site]],age)
      if (methylation_probability >= 1) {
        methylation_probability <- one_met_replacement
      } else if (methylation_probability <= 0) {
        methylation_probability <- zero_met_replacement
      } else {
        methylation_probability <- methylation_probability
      }
      # get single cell binary methylation level
      site_methylation_sc <- cell[site]
      # if the CpG is methylated, append log(methylation_probability)
      if (site_methylation_sc == 1) {
        probability_list_one_age <- append(probability_list_one_age, log(methylation_probability))
      } else if (site_methylation_sc == 0) {
        probability_list_one_age <- append(probability_list_one_age, log(1 - methylation_probability))
      } else {
        stop("Encountereda non-binary methylation value. Methylation values must be either 0 (unmethylated) or 1 (methylated).")
      }
      
      
    }
    sum_age_step <- sum(unlist(probability_list_one_age))
    # append sum_age_step to list_of_profile_probabilities_per_age
    list_of_profile_probabilities_per_age <- append(list_of_profile_probabilities_per_age, sum_age_step)
    pb$tick()
  }
  # transform list_of_profile_probabilities_per_age into a data frame
  profile_probabilities_df <- data.frame(Age = age_steps, Probability = unlist(list_of_profile_probabilities_per_age))
  
  # scale probabilities within 0-1 range
  profile_probabilities_df$Probability <- (profile_probabilities_df$Probability - min(profile_probabilities_df$Probability)) / (max(profile_probabilities_df$Probability) - min(profile_probabilities_df$Probability))
  # 
  # # add uncertainty factor
  #profile_probabilities_df$Probability <- profile_probabilities_df$Probability * (1 - uncertainty) + uncertainty / length(age_steps)
  # 
  # # normalize probabilities
  profile_probabilities_df$Probability <- profile_probabilities_df$Probability / sum(profile_probabilities_df$Probability)
  
  age_pre <- profile_probabilities_df$Age[which.max(profile_probabilities_df$Probability)]
  
  end <- Sys.time()
  print(end - start)
  
  return(age_pre)
}

# Scage bulk
compute_probabilities_bulk <- function(m_obe, # matrix: bulk samples
                                  m.c, # pcc
                                  m.e, # formula
                                  selection_mode = "percentile",
                                  CpG_parameter=0.9,
                                  zero_met_replacement = 0.001,
                                  one_met_replacement = 0.999,
                                  min_age = 0,
                                  max_age = 100,
                                  age_step = 0.1,
                                  uncertainty = 1){
  library(progress)
  start <- Sys.time()
  cuts <- 1
  # profiling mode selection
  # effect size sort
  if (selection_mode == "percentile") { # ex: top 1% age-associated CpGs
    quantile <- quantile(abs(m.c), probs = c(cuts,CpG_parameter))
    abs_top <- m_obe[abs(m.c) >= quantile[2] & abs(m.c) <= quantile[1],]
    #print(quantile[2])
  } else if (selection_mode == "numCpGs") { # ex: top 500 age-associated CpGs
    abs_top <- m_obe[order(-abs(m.c)),][1:CpG_parameter,]
  }
  m_subset <- abs_top
  # isolate selected sites
  selected_sites <- rownames(m_subset)
  # get age steps from min_age to max_age (inclusive of both)
  age_steps <- seq(min_age, max_age + age_step, age_step)
  # create list to store probability profiles
  profile_probabilities_per_age <- matrix(NA,nrow = length(age_steps),ncol = ncol(m_obe))
  rownames(profile_probabilities_per_age) <- age_steps
  # loop through each age step
compute_methylation <- function(site, m.e, age, one_met_replacement, zero_met_replacement, m_subset) {
      #methylation_probability <- regl_bulk(m.e[[site]],frac.c,covs,age)
      methylation_probability <- regl(m.e[[site]],age)
      if (any(methylation_probability >= 1)) {
        methylation_probability[methylation_probability >= 1] <- one_met_replacement
      } else if (any(methylation_probability <= 0)) {
        methylation_probability[methylation_probability <= 0] <- zero_met_replacement
      } else {
        methylation_probability <- methylation_probability
      }
      return(log(1-abs(methylation_probability - m_subset[site,])))
}
pb <- progress_bar$new(total = length(age_steps)) # progress bar
  for (age in age_steps) {
    
    probability_one_age <- matrix(NA,nrow = length(selected_sites),ncol = ncol(m_obe))
   
    rownames(probability_one_age) <- selected_sites
    # loop through each site
    probability_one_age <- lapply(selected_sites, compute_methylation, m.e = m.e, age = age, one_met_replacement = one_met_replacement, zero_met_replacement = zero_met_replacement, m_subset = m_subset)
    probability_one_age <- do.call(rbind, probability_one_age)
    
    # sum per-site probabilities for this age-step
    sum_age_step <- colSums(probability_one_age)
    profile_probabilities_per_age[as.character(age),] <- sum_age_step
    pb$tick()
  }
  
  for ( g in 1:ncol(profile_probabilities_per_age)){
    # scale probabilities within 0-1 range
    profile_probabilities_per_age[,g] <- (profile_probabilities_per_age[,g] - min(profile_probabilities_per_age[,g])) / (max(profile_probabilities_per_age[,g]) - min(profile_probabilities_per_age[,g]))
    # normalize probabilities
    profile_probabilities_per_age[,g] <- profile_probabilities_per_age[,g] / sum(profile_probabilities_per_age[,g])
    }
  # # add uncertainty factor
  #profile_probabilities_df$Probability <- profile_probabilities_df$Probability * (1 - uncertainty) + uncertainty / length(age_steps)

  age_pre <- as.numeric(unlist(apply(profile_probabilities_per_age, 2,function(x) rownames(profile_probabilities_per_age)[which.max(x)])))
  end <- Sys.time()
  print(end - start)
  
  return(age_pre)
}



# Scage parallel

compute_probabilities_parallel <- function(cell, # vector: one cell
                                  cell.c, # vector effect size
                                  cell.e, # vector: formula coefs
                                  selection_mode = "percentile",
                                  CpG_parameter=0.9,
                                  zero_met_replacement = 0.001,
                                  one_met_replacement = 0.999,
                                  min_age = 0,
                                  max_age = 100,
                                  age_step = 0.1,
                                  uncertainty = 1){

  # get arguments
  cell <- cell[names(cell.e)]
  single_cell_name <- names(cell)
  single_cell_met <- cell
  start <- Sys.time()
  # judge <- TRUE
  cuts <- 1
  # profiling mode selection
  # effect size sort
  if (selection_mode == "percentile") { # ex: top 1% age-associated CpGs per cell
    
    quantile <- quantile(abs(cell.c), probs = c(cuts,CpG_parameter))
    abs_top <- subset(cell, abs(cell.c) >= quantile[2] & abs(cell.c) <= quantile[1])
    #print(quantile[2])
  } else if (selection_mode == "numCpGs") { # ex: top 500 age-associated CpGs per cell
    abs_top <- cell[order(-abs(cell.c))][1:CpG_parameter]
  }
  cell_subset <- cell[names(abs_top)]
  # isolate selected sites
  selected_sites <- names(cell_subset)
  
  # get age steps from min_age to max_age (inclusive of both)
  age_steps <- seq(min_age, max_age + age_step, age_step)
  
  # create list to store probability profiles
  list_of_profile_probabilities_per_age <- list()
  
  # loop through each age step
  # Load required libraries
  library(doParallel)
  library(foreach)
  
  # Register the parallel backend
  registerDoParallel(cores=50)
  # Parallelized loop
  list_of_profile_probabilities_per_age <- foreach(age=age_steps, .combine='c', .packages='doParallel') %dopar% {
    
    # create list to hold probability for all chosen CpGs for a given age
    probability_list_one_age <- list()
    
    # loop through each site
    for (site in selected_sites) {
      # compute methylation probability
      methylation_probability <- regl(cell.e[[site]],age)
      if (methylation_probability >= 1) {
        methylation_probability <- one_met_replacement
      } else if (methylation_probability <= 0) {
        methylation_probability <- zero_met_replacement
      } else {
        methylation_probability <- methylation_probability
      }
      
      # get single cell binary methylation level
      site_methylation_sc <- cell[site]
      
      # if the CpG is methylated, append log(methylation_probability)
      if (site_methylation_sc == 1) {
        probability_list_one_age <- append(probability_list_one_age, log(methylation_probability))
      } else if (site_methylation_sc == 0) {
        probability_list_one_age <- append(probability_list_one_age, log(1 - methylation_probability))
      } else {
        stop("Encountered a non-binary methylation value. Methylation values must be either 0 (unmethylated) or 1 (methylated).")
      }
    }
    
    sum_age_step <- sum(unlist(probability_list_one_age))
    return(sum_age_step)
  }
  
  # Stop parallel backend
  stopImplicitCluster()
  
  # transform list_of_profile_probabilities_per_age into a data frame
  profile_probabilities_df <- data.frame(Age = age_steps, Probability = unlist(list_of_profile_probabilities_per_age))
  
  # scale probabilities within 0-1 range
  profile_probabilities_df$Probability <- (profile_probabilities_df$Probability - min(profile_probabilities_df$Probability)) / (max(profile_probabilities_df$Probability) - min(profile_probabilities_df$Probability))
  # 
  # # add uncertainty factor
  #profile_probabilities_df$Probability <- profile_probabilities_df$Probability * (1 - uncertainty) + uncertainty / length(age_steps)
  # 
  # # normalize probabilities
  profile_probabilities_df$Probability <- profile_probabilities_df$Probability / sum(profile_probabilities_df$Probability)
  
  age_pre <- profile_probabilities_df$Age[which.max(profile_probabilities_df$Probability)]
  
  end <- Sys.time()
  print(end - start)
  
  return(age_pre)
}
