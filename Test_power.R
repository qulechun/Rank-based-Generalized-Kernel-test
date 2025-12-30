library(mvtnorm)
library(randtoolbox)
library(clue)
library(MASS)
library(energy)
library(parallel)
library(pbapply)

# --- 1.1 Functions ---
# Bandwidth Selection (Median Heuristic)
med_sigma = function(X, Y) {
  Z <- rbind(X, Y)
  Z_rank <- compurank(Z) 
  return(median(dist(Z_rank))) 
}

# Rank Transform
compurank=function(x)
{
  m=dim(x)[1]
  dim=dim(x)[2]
  gridch=halton(m,dim)
  if(dim==1)
    gridch=matrix((1:m)/m)
  
  distmat=matrix(0,nrow=m,ncol=m)
  for(i in 1:(m))
    distmat[i,]=apply((x[i,]-t(gridch)),2,norm, type="2")^2 # Fixed norm call for stability
  assignmentFUN=solve_LSAP(distmat)
  assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
  return(gridch[assignmentSOL[,2],])
}

# MMD Statistic Calculation
mmd = function(X, Y, sigma) {
  m = nrow(X)
  n = nrow(Y)
  
  Z = rbind(X,Y)
  N = nrow(Z)
  
  K = exp(-as.matrix(dist(Z)^2)/sigma^2/2) - diag(N)
  
  Kx = sum(K[1:m,1:m])/m/(m-1)
  Ky = sum(K[(m+1):N,(m+1):N])/n/(n-1)
  Kxy = sum(K[1:m,(m+1):N])/m/n
  
  MMD = Kx + Ky - 2*Kxy
  return(MMD)
}

# --- 1.2 RGK Statistic (T_N) ---
test_Tn_standard <- function(X, Y, sigma) {
  m <- nrow(X)
  n <- nrow(Y)
  N <- m + n
  Z <- rbind(X, Y)
  Z_rank <- compurank(Z)
  
  K <- exp(-as.matrix(dist(Z_rank)^2) / (2 * sigma^2))
  diag(K) <- 0 
  
  S1 <- sum(K[1:m, 1:m]) 
  S2 <- sum(K[(m + 1):N, (m + 1):N]) 
  mu_K <- sum(K) 
  S12 <- (mu_K - S1 - S2)/2
  
  A <- sum(K^2)
  B <- sum(rowSums(K)^2) - A
  C <- sum(K)^2 - 2 * A - 4 * B
  
  # Tn1 Component
  Tn1 <- S1 / (m*(m-1)) - S12 / (m*n)
  Var_Tn1 <- A/(N*(N-1)) * (m+2*n-1)/(m*n*(m-1)) + 
    B/(N*(N-1)*(N-2)) * (m^2+m*n-3*m-5*n+2)/(m*n*(m-1)) - 
    C/(N*(N-1)*(N-2)*(N-3)) * (m^2+m*n-2*m-3*n+1)/(m*n*(m-1))
  Tn1_sd <- Tn1 / sqrt(Var_Tn1)
  
  # Tn2 Component
  Tn2 <- S2 / (n*(n-1)) - S12 / (m*n)
  Var_Tn2 <- A/(N*(N-1)) * (n+2*m-1)/(m*n*(n-1)) + 
    B/(N*(N-1)*(N-2)) * (n^2+m*n-3*n-5*m+2)/(m*n*(n-1)) - 
    C/(N*(N-1)*(N-2)*(N-3)) * (n^2+m*n-2*n-3*m+1)/(m*n*(n-1))
  Tn2_sd <- Tn2 / sqrt(Var_Tn2)
  
  Tn_final <- Tn1_sd^2 + Tn2_sd^2
  
  return(Tn_final) 
}

# --- 1.3 T_CMZ Statistic (from Tool.R - Chen's Code)---
Perm.homo <- function(X, Y, dX, dY, dXY, index, alpha, R){
  n = dim(dX)[1]
  l = length(index)
  Perm = matrix(NA, nrow=1, ncol=l)
  Z = as.matrix(rbind(X,Y))
  for (i in 1:l){
    Stat = n * ( mean( exp(-dX^index[i]) + exp(-dY^index[i]) - 2*exp(-dXY^index[i]) ) )
    Test.Stat = sapply(1:R, function(bb){
      pn = sample(2*n)
      Dis.pn = as.matrix(dist(Z[pn,]))
      dpX = Dis.pn[1:n,1:n]
      dpY = Dis.pn[(n+1):(2*n),(n+1):(2*n)]
      dpXY = Dis.pn[1:n,(n+1):(2*n)]
      n * ( mean( exp(-dpX^index[i]) + exp(-dpY^index[i]) - 2*exp(-dpXY^index[i]) ) )
    })
    Perm[i] = Stat > sort(Test.Stat)[floor((1-alpha)*R)]
  }
  return(Perm)
}

# --- 1.4 GPK Method (from kerTests.R - Song's Code) ---
kertests = function(X, Y, sigma, r1=1.2, r2=0.8, perm=0) {
  m = nrow(X)
  n = nrow(Y)
  
  Z = rbind(X,Y)
  N = nrow(Z)
  
  K = exp(-as.matrix(dist(Z)^2)/sigma^2/2) - diag(N) # kernel matrix
  
  Kx = sum(K[1:m,1:m])/m/(m-1)
  Ky = sum(K[(m+1):N,(m+1):N])/n/(n-1)
  Kxy = sum(K[1:m,(m+1):N])/m/n
  
  mu_Kx = sum(K)/N/(N-1)
  mu_Ky = sum(K)/N/(N-1)
  mu_Kxy = sum(K)/N/(N-1)
  
  A = sum(K^2)
  B = sum(rowSums(K)^2) - A
  C = sum(K)^2 - 2*A - 4*B
  
  p1 = m*(m-1)/N/(N-1)
  p2 = p1*(m-2)/(N-2)
  p3 = p2*(m-3)/(N-3)
  
  q1 = n*(n-1)/N/(N-1)
  q2 = q1*(n-2)/(N-2)
  q3 = q2*(n-3)/(N-3)
  
  var_Kx = (2*A*p1 + 4*B*p2 + C*p3)/m/m/(m-1)/(m-1) - mu_Kx^2
  var_Ky = (2*A*q1 + 4*B*q2 + C*q3)/n/n/(n-1)/(n-1) - mu_Ky^2
  cov_Kx_Ky = C/N/(N-1)/(N-2)/(N-3) - mu_Kx*mu_Ky
  
  COV = matrix(c(var_Kx,cov_Kx_Ky,cov_Kx_Ky,var_Ky), nrow=2)
  Sinv = solve(COV)
  kmv = c(Kx-mu_Kx, Ky-mu_Ky)
  GPK = as.numeric(kmv%*%Sinv%*%kmv)
  
  result = list()
  result$teststat$GPK = GPK
  
  if (perm>0) {
    temp1 = rep(0, perm)
    for (i in 1:perm) {
      id = sample(N, replace = FALSE)
      K_i = K[id,id]
      
      Kx_i = sum(K_i[1:m,1:m])/m/(m-1)
      Ky_i = sum(K_i[(m+1):N,(m+1):N])/n/(n-1)
      Kxy_i = sum(K_i[1:m,(m+1):N])/m/n
      
      kmv_i = c(Kx_i-mu_Kx, Ky_i-mu_Ky)
      GPK_i = as.numeric(kmv_i%*%Sinv%*%kmv_i)
      
      temp1[i] = GPK_i
    }
    GPK_perm = length(which(temp1>=GPK))/perm
    result$pval$GPK_perm = min(1,GPK_perm)
  }
  
  return(result)
}

# --- 2.Data Generation ---
Gen_Data <- function(scenario, m, n, d, theta) {
  mu0 <- rep(0, d)
  
   if(scenario == "S1") {
    Sigma_X <- outer(1:d, 1:d, function(i,j) 0.3^abs(i-j))
    Sigma_Y <- outer(1:d, 1:d, function(i,j) theta^abs(i-j))
    X <- rmvnorm(m, mean = mu0, sigma = Sigma_X)
    Y <- rmvnorm(n, mean = mu0, sigma = Sigma_Y)
    
  } else if(scenario == "S3") {
    Sigma_X <- matrix(0.3, d, d); diag(Sigma_X) <- 1
    Sigma_Y <- matrix(theta, d, d); diag(Sigma_Y) <- 1
    X <- rmvnorm(m, mean = mu0, sigma = Sigma_X)
    Y <- rmvnorm(n, mean = mu0, sigma = Sigma_Y)
    
  } else if(scenario == "S2") {
    Sigma_X <- outer(1:d, 1:d, function(i,j) 0.3^abs(i-j))
    Sigma_Y <- outer(1:d, 1:d, function(i,j) theta^abs(i-j))
    X <- exp(rmvnorm(m, mean = mu0, sigma = Sigma_X))
    Y <- exp(rmvnorm(n, mean = mu0, sigma = Sigma_Y))
  }else if(scenario == "S4") {
    Sigma_X <- matrix(0.3, d, d); diag(Sigma_X) <- 1
    Sigma_Y <- matrix(theta, d, d); diag(Sigma_Y) <- 1
    
    X_raw <- rmvnorm(m, mean = mu0, sigma = Sigma_X)
    Y_raw <- rmvnorm(n, mean = mu0, sigma = Sigma_Y)
    
    X <- exp(X_raw)
    Y <- exp(Y_raw)
  }
  return(list(X=X, Y=Y))
}

# --- 3.Pre-computation of Critical Values (Tn & RMMD) ---
Get_Crit_Values <- function(n_list, d_list, B=10000, alpha=0.05, n_cores=14) {
  file_path <- "Critical_Values_Tn_RMMD.RData"
  
  if(file.exists(file_path)) {
    load(file_path)
  } else {
    Crit_List <- list()
  }
  
  tasks <- list()
  for(n in n_list) {
    for(d in d_list) {
      key <- paste0("n", n, "_d", d)
      if(is.null(Crit_List[[key]])) {
        tasks[[key]] <- list(n=n, d=d)
      }
    }
  }
  
  if(length(tasks) == 0) {
    return(Crit_List)
  }
  
  
  cl <- makeCluster(n_cores)
  
  clusterExport(cl, varlist = c("test_Tn_standard", "med_sigma", "compurank", "mmd", 
                                "halton", "solve_LSAP", "dist"))
  clusterEvalQ(cl, { 
    library(mvtnorm)
    library(randtoolbox)
    library(clue)
    library(MASS)
  })
  
  pboptions(type = "txt", char = "=")
  
  for(task_name in names(tasks)) {
    n <- tasks[[task_name]]$n
    d <- tasks[[task_name]]$d
    
    res_matrix <- pbsapply(1:B, function(b) {
      X0 <- matrix(runif(n*d), n, d)
      Y0 <- matrix(runif(n*d), n, d)
      
      sig <- med_sigma(X0, Y0)
      
      tn_val <- test_Tn_standard(X0, Y0, sig)
      
      Z_rank <- compurank(rbind(X0, Y0))
      rmmd_val <- mmd(Z_rank[1:n,], Z_rank[(n+1):(2*n),], sig)
      
      return(c(tn_val, rmmd_val))
    }, cl = cl) 
    
    tn_vec <- res_matrix[1, ]
    rmmd_vec <- res_matrix[2, ]
    
    Crit_List[[task_name]] <- list(
      Tn = quantile(tn_vec, 1-alpha), 
      RMMD = quantile(rmmd_vec, 1-alpha)
    )
    save(Crit_List, file = file_path)
  }
  
  stopCluster(cl)
  return(Crit_List)
}

# --- 4.Main Simulation ---
Run_Simulation_Batch <- function(n, d, scenario_list, theta_list, Crit_Table, n_cores, R_perm=500) {
  
  cl <- makeCluster(n_cores)
  
  clusterExport(cl, varlist = c("test_Tn_standard", "med_sigma", "compurank", "mmd", 
                                "Perm.homo", "kertests", "Gen_Data_Deb", 
                                "halton", "solve_LSAP", "Crit_Table"))
  
  clusterEvalQ(cl, {
    library(mvtnorm)
    library(randtoolbox)
    library(clue)
    library(MASS)
    library(energy)
  })
  
  results <- data.frame()
  
  for(scen in scenario_list) {
    for(theta in theta_list) {
      cat(sprintf("Scen: %s | Theta: %.2f | n: %d | d: %d\n", scen, theta, n, d))
      
      res_matrix <- parSapply(cl, 1:500, function(i) {
        dat <- Gen_Data(scen, n, n, d, theta)
        X <- dat$X; Y <- dat$Y
        sig <- med_sigma(X, Y)
        
        # 1. TN
        stat_tn <- test_Tn_standard(X, Y, sig)
        crit_tn <- Crit_Table[[paste0("n",n,"_d",d)]]$Tn
        rej_tn <- stat_tn > crit_tn
        
        # 2. RMMD
        z_r <- compurank(rbind(X, Y))
        stat_rmmd <- mmd(z_r[1:n,], z_r[(n+1):(2*n),], sig)
        crit_rmmd <- Crit_Table[[paste0("n",n,"_d",d)]]$RMMD
        rej_rmmd <- stat_rmmd > crit_rmmd
        
        # 3. T_CMZ (Perm)
        Dis <- as.matrix(dist(rbind(X, Y)))
        s_res <- Perm.homo(X, Y, Dis[1:n,1:n], Dis[(n+1):(2*n),(n+1):(2*n)], Dis[1:n,(n+1):(2*n)], c(1), 0.05, R_perm)
        rej_simos <- s_res[1]
        
        # 4. GPK (Perm)
        gpk_res <- kertests(X, Y, sig, perm=R_perm)
        rej_gpk <- gpk_res$pval$GPK_perm < 0.05
        
        return(c(rej_tn, rej_rmmd, rej_simos, rej_gpk))
      })
      
      pow <- rowMeans(res_matrix)
      results <- rbind(results, data.frame(Scenario=scen, Theta=theta, n=n, d=d, Method=c("Tn","RMMD","Simos","GPK"), Power=pow))
    }
  }
  stopCluster(cl)
  return(results)

}

# --- 5.Execution Control ---
Crit_Table <- Get_Crit_Values(n_list=c(100), d_list=c(2, 4), B=10000, n_cores=14)
d_list <- c(2, 4)
run_n <- c(100,200)
n_cores <- 12    # Adjust based on your CPU

all_results <- data.frame()
for (d in d_list) {

  res_S1 <- Run_Simulation_Batch(run_n, d, c("S1"), 
                                     seq(0.3, 0.7, 0.1), 
                                     Crit_Table, n_cores)
  res_S2 <- Run_Simulation_Batch(run_n, d, c("S2"), 
                                     seq(0.3, 0.7, 0.1), 
                                     Crit_Table, n_cores)
  res_S3 <- Run_Simulation_Batch(run_n, d, c("S3"), 
                                     seq(0.3, 0.7, 0.1), 
                                     Crit_Table, n_cores)
  res_S4 <- Run_Simulation_Batch(run_n, d, c("S4"), 
                                     seq(0.3, 0.7, 0.1), 
                                     Crit_Table, n_cores)
  
  all_results <- rbind(all_results, res_S1, res_S2,res_S3,res_S4)
  
  write.csv(all_results, paste0("Sim_Results_n", run_n, ".csv"), row.names = F)
}