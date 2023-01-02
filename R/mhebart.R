#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Multiple Hierarchical Embedded Bayesian Additive Regression Trees
#' @description This function runs a mhebart model and returns the tree and
#' other results obtained in the last iteration of the MCMC
#' @param formula The model formula
#' @param data The data to be used in the modeling
#' @param group_variable The name of the grouping variable
#' @param num_trees The number of trees (P)
#' @param control A list with control settings
#' @param priors A list with prior hyperparameters as defined by the model
#' @param inits A list with initial values for parameters
#' @param MCMC A list with MCMC parameters
#' @return A list containing:
#'  Everything
#' @details
#' Priors used ----------------------------------------------------------
#' y_{ij} ~ Normal(m_j, tau^-1)
#' tau    ~ Gamma(nu/2, nu*lambda/2)
#' mu     ~ Normal(0, tau_mu^-1)
#' phi    ~ Normal(mu, sigma_phi^2 / T)
#' ----------------------------------------------------------------------

mhebart <- function(formula,
                   data,
                   group_variables,
                   
                   # X is the feature matrix, y is the target,
                   # groups, # groups is the group number of each obs
                   num_trees = 10, # Number of trees
                   control = list(node_min_size = 5), # Size of smallest nodes
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     nu = 2,
                     lambda = 0.1,
                     tau_mu = 16 * num_trees,
                     shape_sigma_phi = 0.5,
                     scale_sigma_phi = 1,
                     sample_sigma_phi = TRUE
                   ),
                   inits = list(
                     tau = 1,
                     tau_phi = 1, 
                     sigma_phi = 1
                   ), # Initial values list
                   MCMC = list(
                     iter = 250, # Number of iterations
                     burn = 50, # Size of burn in
                     thin = 1,
                     sigma_phi_sd = 2
                   )) {

  # Handling formula interface
  formula_int <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula_int)[1]
  names_x <- all.vars(formula_int[[3]])
  n_grouping_variables <- length(group_variables)
  grouping_variables_names <- paste0("group_", 1:n_grouping_variables)

  # Used in create_S to avoid error with stumps
  mod.mat <<- function(f) {
    if(nlevels(f) == 1) {
      m <- matrix(1, nrow = length(f), ncol = 1)
    } else {
      m <- stats::model.matrix(~ f - 1)  
    }
    return(m)
  }
  
  #-------------------------------------------------------
  data <- dplyr::select(data, c(!!response_name, !!names_x, !!group_variables))
  # data    <- dplyr::select(data, c(!!response_name, !!names_x, "group"))
  names(data)[names(data) %in% group_variables] <- grouping_variables_names
  
  data <- data |> 
    dplyr::mutate_if(is.character, stringr::str_squish) |> 
    dplyr::mutate_if(is.character, stringr::str_to_lower) |>
    dplyr::mutate_if(is.character, abjutils::rm_accent) |> 
    dplyr::mutate_if(is.character, ~stringr::str_replace(.x, " ", "_"))
  
  mf <- stats::model.frame(formula_int, data = data)
  X <- as.matrix(stats::model.matrix(formula_int, mf))
  y <- stats::model.extract(mf, "response")

  # Extract control parameters
  node_min_size <- control$node_min_size

  # Extract hyper-parameters
  alpha  <- priors$alpha # Tree shape parameter 1
  beta   <- priors$beta # Tree shape parameter 2
  tau_mu <- priors$tau_mu # Overall mean precision
  shape_sigma_phi  <- priors$shape_sigma_phi # Weibull prior parameters 
  scale_sigma_phi  <- priors$scale_sigma_phi # the same will be used for all the groups
  sample_sigma_phi <- priors$sample_sigma_phi

  # Extract initial values
  tau       <- inits$tau
  sigma     <- 1 / sqrt(tau)
  #sigma_phi <- inits$sigma_phi
  #tau_phi   <- 1 / (sigma_phi^2)
  log_lik   <- 0
  
  # Tree numbers 
  initial <- round(num_trees/n_grouping_variables)
  n_trees <- seq(initial, num_trees, by = initial)
  
  # Extract MCMC details
  iter <- MCMC$iter # Number of iterations
  burn <- MCMC$burn # Size of burn in
  thin <- MCMC$thin # Amount of thinning
  sigma_phi_sd <- MCMC$sigma_phi_sd # SD parameter for sigma_phi MH update
  
  # Storage containers
  store_size      <- (iter - burn) / thin
  tree_store      <- vector("list", store_size)
  sigma_store     <- rep(NA, store_size)
  y_hat_store     <- matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store   <- rep(NA, store_size)
  sigma_phi_store <- list()
  tau_store       <- rep(NA, store_size)
  mse_store       <- rep(NA, store_size)

  # Scale the response target variable
  y_min <- min(y)
  y_max <- max(y)
  y_scale <- (y - y_min)/(y_max - y_min) - 0.5
  n <- length(y_scale)

  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to tau
  #---------------------------------------------------------------------
  str_groups <- paste(paste0("+ (1|", grouping_variables_names, ")"), collapse = ' ')
  lme_form <- paste(paste(c(formula)), str_groups, collapse = '')
  
  lme_form <- stats::as.formula(lme_form)
  data_lme <- dplyr::mutate(data, y = y_scale)
  my_lme <- lme4::lmer(lme_form, data_lme)
  # my_lme   <- brms::brm(lme_form, data_lme,
  #                      silent = TRUE)
  #res      <- stats::sd(stats::fitted(my_lme)[, "Estimate"] - y_scale)
  res      <- lme4:::sigma.merMod(my_lme)
  
  nu     <- priors$nu         # Parameter 1 for precision
  lambda <- priors$lambda # Parameter 2 for precision
  p_inv  <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
  
  # Putting high probabilities of the BCART improving a linear model
  while(p_inv < 0.95){
    p_inv <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
    if(p_inv < 0.95){
      nu = abs(nu + stats::rnorm(1))
      lambda = abs(lambda + stats::rnorm(1))
    }
  }
  
  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to the sigma_phi
  # In this case we need two priors 
  #---------------------------------------------------------------------
  shape_sigma_phi <- vector(length = n_grouping_variables)
  scale_sigma_phi <- vector(length = n_grouping_variables)
  sigma_phi <- rep(1, length = n_grouping_variables)
  tau_phi <-  1 / (sigma_phi^2)
  
  for(i in 1:n_grouping_variables){
    str_groups <- paste0("+ (1|", grouping_variables_names[i], ")")
    lme_form <- paste(paste(c(formula)), str_groups, collapse = '')
    lme_form <- stats::as.formula(lme_form)
    data_lme <- dplyr::mutate(data, y = y_scale)
    my_lme <- lme4::lmer(lme_form, data_lme)
    random_effect     <- sqrt(as.data.frame(lme4::VarCorr(my_lme))$vcov[1])
    pr <- parameters::model_parameters(my_lme, effects = "random",
                                       ci_random = TRUE,
                                       verbose = FALSE)
    se <- pr$SE[1]
    random_effect_var <- se^2
    
    if(!is.na(se)){
      shape_sigma_phi[i]  <-  (random_effect^2)/random_effect_var
      scale_sigma_phi[i]  <-  random_effect_var/random_effect
    }
  }
    

  #---------------------------------------------------------------------
  # Get the group objects
  
  M_all <- list()
  group_sizes_all <- list()
  num_groups_all <- vector(length = n_grouping_variables)
  curr_trees  <- list()
  predictions <- list()
  
  for(n_g in 1:n_grouping_variables){
    
    groups <- data[, grouping_variables_names[n_g]]
    if(!is.vector(groups)){
      groups <- dplyr::pull(groups, !!grouping_variables_names[n_g])
    }
    
    M_all[[n_g]] <- stats::model.matrix(~ factor(groups) - 1)
    group_sizes_all[[n_g]] <- table(groups)
    num_groups_all[n_g]  <- length(group_sizes_all[[n_g]])
    num_trees <- n_trees[n_g]
  
  # Create a list of trees for the initial stump
    curr_trees[[n_g]] <- create_stump(
      num_trees = initial,
      groups = groups,
      y = y_scale,
      X = X
    )
    predictions[[n_g]] <- get_group_predictions(
      trees = curr_trees[[n_g]], X, groups, 
      single_tree = num_trees == 1, 
      old_groups = groups)
    
  }
  
  
  # Set up a progress bar
  pb <- utils::txtProgressBar(
    min = 1, max = iter,
    style = 3, width = 60,
    title = "Running rBART..."
  )

  # Start the iterations loop
  for (i in 1:iter) {
    utils::setTxtProgressBar(pb, i)
    
      
      # If at the right place store everything
      if ((i > burn) & ((i %% thin) == 0)) {
        curr <- (i - burn) / thin
        tree_store[[curr]] <- curr_trees
        sigma_store[curr] <- sigma
        y_hat_store[curr, ] <- predictions_all
        #log_lik_store[curr] <- log_lik
        sigma_phi_store[[curr]] <- sigma_phi
        tau_store[curr] <- tau
        mse_store[curr] <- mse 
      }
    
    
    for(n_g in 1:n_grouping_variables){  
      #groups <- data[, grouping_variables_names[n_g]]
      
      # Setting this groups'parameters
      num_trees <-  initial
      M <- M_all[[n_g]]
      group_sizes <- group_sizes_all[[n_g]]
      num_groups  <- num_groups_all[n_g]
      predictions_all <- vector(length = length(y))
        
      # Start looping through trees
      for (j in 1:num_trees) {
        # Calculate partial residuals for current tree
        if (num_trees > 1) {
          if(n_g == 1){
            groups <- data[, grouping_variables_names[n_g]]
            if(!is.vector(groups)){
              groups <- dplyr::pull(groups, !!grouping_variables_names[n_g])
            }
            partial_trees <- curr_trees[[n_g]]
            partial_trees[[j]] <- NULL # Blank out that element of the list
            current_partial_residuals <- y_scale -
              get_group_predictions(
                trees = partial_trees, X, groups,
                single_tree = num_trees == 2, 
                old_groups = groups
              )
          } else {
            
            current_predictions_all <- vector(length = length(y_scale))
            for(n_g_i in 1:n_g){
              partial_trees <- curr_trees[[n_g_i]]
              groups <- data[, grouping_variables_names[n_g_i]]
              if(!is.vector(groups)){
                groups <- dplyr::pull(groups, !!grouping_variables_names[n_g_i])
              }
              
              partial_trees[[j]] <- NULL # Blank out that element of the list
              res_predictions    <- get_group_predictions(
                trees = partial_trees, X, groups,
                single_tree = num_trees == 2, 
                old_groups = groups
              )  
              current_predictions_all <-  current_predictions_all + res_predictions
              current_partial_residuals <- y_scale - current_predictions_all
            }
          }
        } else {
          current_partial_residuals <- y_scale
        }
        
        # Propose a new tree via grow/change/prune/swap
        type <- sample(c("grow", "prune", "change", "swap"), 1)
        if (i < max(floor(0.1 * burn), 10)) type <- "grow" # Grow for the first few iterations
        
        # Get a new tree!
        new_trees <- curr_trees[[n_g]]
        
        new_trees[[j]] <- update_tree(
          y = y_scale,
          X = X,
          groups = groups,
          type = type,
          curr_tree = curr_trees[[n_g]][[j]],
          node_min_size = node_min_size
        )
        # Calculate the complete conditional and acceptance probability
        l_new <- full_conditional_hebart(
          tree = new_trees[[j]],
          R = current_partial_residuals,
          num_trees,
          tau,
          tau_phi[n_g],
          tau_mu,
          M
        ) +
          get_tree_prior(new_trees[[j]], alpha, beta)
        
        
        l_old <- full_conditional_hebart(
          curr_trees[[n_g]][[j]],
          current_partial_residuals,
          num_trees,
          tau,
          tau_phi[n_g],
          tau_mu,
          M
        ) +
          get_tree_prior(curr_trees[[n_g]][[j]], alpha, beta)
        
        # If accepting a new tree update all relevant parts
        log.alpha <- (l_new - l_old)
        accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
        if (accept) curr_trees[[n_g]] <- new_trees
        
        # Update mu whether tree accepted or not
        curr_trees[[n_g]][[j]] <- simulate_mu_hebart(
          tree = curr_trees[[n_g]][[j]],
          R = current_partial_residuals,
          tau,
          tau_phi[n_g],
          tau_mu,
          M,
          num_trees
        )
        
        # Update phi as well
        curr_trees[[n_g]][[j]] <- simulate_phi_hebart(
          tree = curr_trees[[n_g]][[j]],
          R = current_partial_residuals,
          groups,
          tau,
          tau_phi = tau_phi[n_g],
          M,
          num_trees
        )
        
      # Check the trees
      if (any(curr_trees[[n_g]][[j]]$tree_matrix[, "node_size"] < node_min_size)) browser()
    } # End loop through trees
      
      preds <- get_group_predictions(
        trees = curr_trees[[n_g]], 
        X, 
        groups, 
        single_tree = num_trees == 1,
        old_groups = groups
      )
      
      predictions_all <- predictions_all + preds 
    }
    
    # Calculate full set of predictions
    

    mse <- mean((y_scale - predictions_all)^2)
    # Update tau
    
    tau <- update_tau(
      y_scale,
      predictions_all,
      nu,
      lambda
    )
    sigma <- 1 / sqrt(tau)

    
    for(n_g in 1:n_grouping_variables){
      
      groups <- data[, grouping_variables_names[n_g]]
      if(!is.vector(groups)){
        groups <- dplyr::pull(groups, !!grouping_variables_names[n_g])
      }
      
      curr_trees_ng <- curr_trees[[n_g]]
      
      # Update tau_phi
      S1 <- create_S(curr_trees_ng, groups)
      S2 <- create_S(curr_trees_ng)
      
      if(sample_sigma_phi){
        sigma_phi[n_g] <- update_sigma_phi(
          y_scale, S1, S2, sigma_phi[n_g], tau_mu, tau,
          shape_sigma_phi[n_g], scale_sigma_phi[n_g], 
          num_trees = initial, sigma_phi_sd
        )
      }
      tau_phi[n_g] <- 1 / (sigma_phi[n_g]^2)
      
    }
    sigma_phi_store[[i]] <- sigma_phi

    # Get the overall log likelihood
    #Omega_y <- diag(n)/tau + tcrossprod(S1)/(num_trees * tau_phi) + 
    #  tcrossprod(S2)/tau_mu
    #log_lik <- mvnfast::dmvn(
    #  y, rep(0, n), Omega_y, log = TRUE)
      
  } # End iterations loop
  cat("\n") # Make sure progress bar ends on a new line
  final_groups <- list()
  for(n_g in 1:n_grouping_variables){
    groups_final <- unique(data[, grouping_variables_names[n_g]])
    if(!is.vector(groups_final)){
      groups_final <- dplyr::pull(groups_final, 1)
    }
    final_groups[[n_g]] <-  unique(groups_final)
  }
  
  result <- list(
    trees = tree_store,
    sigma = sigma_store,
    y_hat = (y_hat_store + 0.5) * (max(y) - min(y)) + min(y),
    #log_lik = log_lik_store,
    sigma_phi = sigma_phi_store,
    tau = tau_store, 
    mse = mse_store, 
    y = y,
    X = X,
    groups = final_groups, 
    iter = iter,
    burn = burn,
    thin = thin,
    store_size = store_size,
    num_trees = num_trees,
    formula = formula,
    y_min = y_min,
    y_max = y_max
  )

  # RMSE calculation
  pred <- predict_mhebart(
    newX = data, group_variables = group_variables, 
    hebart_posterior = result, type = "mean")
  
  mse <- mean((pred - y)^2)
  rmse <- sqrt(mse)
  r.squared <- 1 - mse / stats::var(y)

  result$rmse <- rmse
  result$r.squared <- r.squared
  result$num_variables <- length(names_x)

  class(result) <- "hebart"

  return(result = result)
} # End main function