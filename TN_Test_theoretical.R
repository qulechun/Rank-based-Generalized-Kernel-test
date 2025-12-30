source("Test_power")
Run_Theoretical_Test <- function(n_list, d_list, scenario_list, theta_list, n_cores=14) {
  crit_val_theo <- 2 * qchisq(0.95, df=1)
  cat(sprintf(">>> Using Theoretical Critical Value: %.4f (Based on 2*Chi^2(1)) <<<\n", crit_val_theo))
  
  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, 12345) 
  clusterExport(cl, varlist = ls(globalenv()))
  clusterEvalQ(cl, { library(mvtnorm); library(randtoolbox); library(clue); library(MASS) })
  
  results <- data.frame()
  
  for(scen in scenario_list) {
    for(n in n_list) {
      for(d in d_list) {
        cat(sprintf("\nRunning: %s | n=%d | d=%d\n", scen, n, d))
        
        for(theta in theta_list) {
          res_vec <- parSapply(cl, 1:1000, function(i) {
            dat <- Gen_Data_Deb(scen, n, n, d, theta)
            sig <- med_sigma(dat$X, dat$Y)
            stat <- test_Tn_standard(dat$X, dat$Y, sig)
            return(ifelse(stat > crit_val_theo, 1, 0))
          })
          
          power <- mean(res_vec)
          cat(sprintf("  Theta=%.2f -> Power=%.3f\n", theta, power))
          
          results <- rbind(results, data.frame(
            Scenario = scen, Theta = theta, n = n, d = d, 
            Method = "Tn_Theoretical", 
            Power = power
          ))
        }
      }
    }
  }
  
  stopCluster(cl)
  return(results)
}

# parameter
n_list <- c(100, 200)
d_list <- c(2, 4)
scenarios <- c("S1", "S2", "S3", "S4")
thetas <- seq(0.3, 0.7, 0.1)
my_cores <- parallel::detectCores() - 14
if(my_cores < 1) my_cores <- 1

res_final <- Run_Theoretical_Test(n_list, d_list, scenarios, thetas, n_cores=my_cores)

write.csv(res_final, "Sim_Results_Theoretical.csv", row.names = FALSE)
