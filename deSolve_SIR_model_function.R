# SIR with ME model 
sim_SIR_ME_BG1 <- function(inits, # must be named vector with S, I, R
                           popsize = 10000, 
                           params, #vector of parameters
                           model_data, # data frame to calculate transmission
                           start_date = as.Date("2020-03-01"), end_date = as.Date("2020-06-30")){
  
  
  trans_df <- model_data %>%
    filter(between(date, start_date, end_date)) %>%
    mutate(tm = exp(log(b1) + beta_ld1*lockdown1 + beta_BG1*BG1))
  
  
  SIR_model <- function(t, Di, ME_Ia, ME_Ib, trans){
    time = seq(1, t, 1)
    tm <- approxfun(x = trans$day, y = trans$tm, method = "linear", rule = 2)
    
    initState <- inits
    names(initState) <- c("S", "I", "R")
    
    N <- popsize
    parameters <- c(gamma = 1/Di, ME_I_a = ME_Ia, ME_I_b = ME_Ib)
    names(parameters) <- c("gamma", "ME_I_a", "ME_I_b")
    
    eqn <- function(time, initState, parameters, trans_input){
      with(as.list(c(initState, parameters)), {
        dS = -tm(time)*S*I/N
        dI = tm(time)*S*I/N - gamma*I
        dR = gamma*I 
        return(list(c(S = dS, I = dI, R = dR), transmission = tm(time), 
                    IncI_scaled = tm(time)*S*I/N,  
                    IncI_scaled_ME = rnorm(1, tm(time)*S*I/N, (ME_I_a + ME_I_b*(tm(time)*S*I/N))/1.96)))
      })
    }
    
    out <- ode(y = initState, times = time, func = eqn, parms = parameters)
    out.df <- as.data.frame(out)
    
    return(out.df)
  }
  
  sim <- SIR_model(t = nrow(trans_df), trans = trans_df, 
                   Di = params["Di"], ME_Ia = params["ME_Ia"], ME_Ib = params["ME_Ib"])
  return(sim)
}


Rt_comp_fun <- function(Inc_df, Inc_col_name, model_name, tstart = NULL, tend = NULL, Rt_ref_df, 
                        Rt_prior = 4, Rt_sd_prior = 3, meansi = 5.1, stdsi = 5.2){
  # select input data according to desired input for Rt calculation
  Inc_data <- Inc_df[c("dept_id", "time", Inc_col_name)]
  colnames(Inc_data) <- c("id", "time", "Inc")
  
  # Rt calculation
  Rt_list <- vector(mode = "list")
  
  if(is_null(tstart) | (is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(mean_si = meansi, 
                                                         std_si = stdsi, 
                                                         mean_prior = Rt_prior, 
                                                         std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = ceiling((t_end+t_start)/2)-1,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  } 
  
  if(!is_null(tstart) & !(is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(
                                 mean_si = meansi, std_si = stdsi, 
                                 t_start = tstart, t_end = tend, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = ceiling((t_end+t_start)/2),
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  Rt_comp_df <- Rt_df %>%
    left_join(., Rt_ref_df, by = c("dept_id", "day")) %>%
    mutate(Rt_residual = Rt - Rt_real, 
           model = model_name)
  
  RMSE_df <- Rt_comp_df %>%
    filter(!is.na(Rt)) %>%
    group_by(dept_id, model) %>%
    summarize(RMSE = sqrt(sum(Rt_residual^2)/n())) %>%
    ungroup() %>%
    mutate(RMSE_total = mean(RMSE))
  
  return(list(Rt_comp = Rt_comp_df, RMSE = RMSE_df))
  
}


Rt_comp_fun_uncertain <- function(Inc_df, Inc_col_name, model_name, tstart = NULL, tend = NULL, Rt_ref_df, 
                                  meansi = 5, stdsi = 2, Rt_prior = 4, Rt_sd_prior = 3){
  # select input data according to desired input for Rt calculation
  Inc_data <- Inc_df[c("dept_id", "time", Inc_col_name)]
  colnames(Inc_data) <- c("id", "time", "Inc")
  
  # Rt calculation
  Rt_list <- vector(mode = "list")
  
  if(is_null(tstart) | (is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "uncertain_si",
                               config = make_config(list(
                                 mean_si = meansi, 
                                 min_mean_si = meansi - 2,
                                 max_mean_si = meansi + 2,
                                 std_mean_si = 1.5,
                                 std_std_si = 1.5,
                                 std_si = stdsi,
                                 min_std_si = stdsi*.8,
                                 max_std_si = stdsi*1.2,
                                 n1 = 50,
                                 n2 = 100, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = ceiling((t_end+t_start)/2)-1,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  } 
  
  if(!is_null(tstart) & !(is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "uncertain_si",
                               config = make_config(list(
                                 mean_si = meansi, 
                                 min_mean_si = meansi - 2,
                                 max_mean_si = meansi + 2,
                                 std_mean_si = 1.5,
                                 std_std_si = 1.5,
                                 std_si = stdsi,
                                 min_std_si = stdsi*.8,
                                 max_std_si = stdsi*1.2,
                                 n1 = 50,
                                 n2 = 100, 
                                 t_start = tstart, 
                                 t_end = tend, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = ceiling((t_end+t_start)/2),
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  Rt_comp_df <- Rt_df %>%
    left_join(., Rt_ref_df, by = c("dept_id", "day")) %>%
    mutate(Rt_residual = Rt - Rt_real, 
           model = model_name)
  
  RMSE_df <- Rt_comp_df %>%
    filter(!is.na(Rt)) %>%
    group_by(dept_id, model) %>%
    summarize(RMSE = sqrt(sum(Rt_residual^2)/n())) %>%
    ungroup() %>%
    mutate(RMSE_total = mean(RMSE))
  
  return(list(Rt_comp = Rt_comp_df, RMSE = RMSE_df))
  
}
