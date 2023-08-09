library(MASS)
library(statmod)
library(Matrix)
library(matrixcalc)
library(corpcor)

library(nloptr)
library(glmnet)
library(BayesLogit)


##################### working graph ###########################
## x : 0 refers  working graph = no graph
##     1 refers  working graph = star-like
##     2 refers  working graph = adding edges between pathways  
##                               (based on star-like pathway)
##                               p = 0.3
##     3 refers  working graph = removing edges 
##                              (based on star-like pathway)
##                               p = 0.3
##     if ne=='T then add noise p=0.2 based on the corresponding graph

working_graph <- function(x,pathway_list,H,data_dim,ind_s,ind_e) {
  graph = matrix(0,p,p)
  if(x==1){
    graph = matrix(0,p,p)
    return(graph)
  }else if(x==2){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }else if(x==3){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        # within pathway
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
        # across pathway
        if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
          graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.2)
        }
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }
}


###########################  empty chain ##########################
empty_chain <- function(n,p,T,L,rho_ini,phi_ini,w_ini,z_ini,alpha_ini,tauz_ini,xiz_ini){
  chain_m   = matrix(0,ncol = p, nrow = T+2)
  
  chain_rho   = matrix(0,ncol = p*n, nrow = T+2)
  chain_rho[1,]=as.vector(t(rho_ini))
  
  chain_w   = matrix(0,ncol = p*L, nrow = T+2)
  chain_w[1,]=as.vector(t(w_ini))
  
  chain_phi   = matrix(0,ncol = H*L, nrow = T+2)
  chain_phi[1,]=as.vector(t(phi_ini))
  
  chain_tau   = matrix(0,ncol = p*L, nrow = T+2)
  chain_tau[1,]=as.vector(t(tau_ini))
  
  chain_tauz   = matrix(0,ncol = L*n, nrow = T+2)
  chain_tauz[1,]=as.vector(t(tauz_ini))
  
  chain_z   = matrix(0,ncol = L*n, nrow = T+2)
  chain_z[1,]=as.vector(t(z_ini))
  
  chain_xiz   = matrix(0,ncol = L*n, nrow = T+2)
  chain_xiz[1,]=as.vector(t(xiz_ini))
  
  chain_alpha =matrix(0,ncol = p*L, nrow = T+2)
  chain_alpha[1,]=as.vector(t(alpha_ini))
  
  list(chain_tau=chain_tau,chain_alpha=chain_alpha,chain_phi=chain_phi,chain_rho=chain_rho,chain_m=chain_m,chain_w=chain_w,chain_z=chain_z,chain_tauz=chain_tauz,chain_xiz=chain_xiz)
}



## propto density of alpha_l
# alpha_l : candidate value
# alpha_h : current value
# alpha   : p by 1 alpha vector 
# omega   : the precison matrix 
# tau_l   : tau[j,l]^2
# j : the j-th entry of alpha 

p_density_mmh <- function(alpha_l,alpha_h,tau_l,alpha,omega,j) {
  lambda_l = exp(alpha_l)
  lambda_h = exp(alpha_h)
  alp_l_aug = replace(alpha,j,alpha_l)
  omega_j = omega[j,]
  ans = (lambda_l^2)/(lambda_h^2)*exp(-(lambda_l^2-lambda_h^2)*tau_l/2)*exp((-1/(2*nu_2))*(alpha_l-alpha_h)*(omega_j%*%(alpha+alp_l_aug-2*nu_1)))
  return(ans)
}

sample_quantile = function(x){
  probs = c(0.025,0.975)
  ans = quantile(x,probs = probs)
  return(ans)
} 


#### MCMC_algorithm #### 

## arguments of the function

# dt : {0,1,2,3} = {gaussian, binary, binomial,mix}
# T : the number of iterations of MCMC
# L : the number of factors  
# p : the number of features
# n : the number of subjects 
# X : data 


# nu_1  :  mean 
# nu_2  :  variance 
# Sigma :  variance for each entry of m
# Q : covariance matrix for the proposal density
# eta : hyper parameter specifying the prior of omega
# eps : hyper parameter specifying the prior of omega 

# rho_temp : initial value for rho 
# tau_temp : initial value for tau 
# omega_temp : initial value for omega
# inv_omega_temp : initial value for inverse omega
# w_temp : initial value for w
# z_temp : initial value for z
# alpha_temp : initial value for 
# m_temp : initial value for

# th  : threshold to determine the degree of freedom (of W)

##########################################

## values of the function
# w_est : estimation of W
# z_est : estimation of Z
# m_est : estimation of m 
# df    : degree of freedom 
# BIC   : BIC criterion 
# chain_alpha : markov chain of alpha 
# chain_m     : markov chain of m 
# chain_w     : markov chain of W
# chain_z     : markov chain of Z 

## new 
BFGA_MCMC <- function(dt,T,L,p,n,X,trials,nu_1,nu_2,nu_3,nu_4,Sigma,Q,eta,eps,rho_temp,tau_temp,omega_temp,inv_omega_temp,w_temp,z_temp,alpha_temp,m_temp,tauz_temp,xiz_temp,start,end,mu){
  
  # preparation based on the data type 
  
  if(dt==0){
    psi = X
    phi =  as.matrix(rep(1,p))
    kappa = matrix(0,nrow = p,ncol = n)
    rowsum_kappa = rep(0,p)
  }else if(dt==1){
    psi = matrix(0,nrow = p,ncol = n)
    # b = matrix(1,ncol = n,nrow = p) # binary data 
    b = trials # binomial  data 
    kappa = X-b/2
    rowsum_kappa= rowSums(kappa)
  }else if(dt==2){
    psi = matrix(0,nrow = p,ncol = n)
    b = X + trials
    kappa = X-b/2
    rowsum_kappa= rowSums(kappa)
  }else if(dt==3){
    psi =  matrix(0,nrow = p,ncol = n)
    psi[1:n_nz,] = X[1:n_nz,]
    phi = as.matrix(rep(1,n_nz))
    b  = matrix(1,ncol = n,nrow = p)
    b[(2*n_nz+1):p,] = X[(2*n_nz+1):p,] + trials[(2*n_nz+1):p,] 
    kappa = X - b/2
    kappa[1:n_nz,] = 0
    rowsum_kappa= rowSums(kappa)
  }
  
  for (t in 1:T) {
    
    if(dt==0){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rowSums(rho_temp*psi)
    }else if(dt==1 |dt==2){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rep(0,p)
    }else if(dt==3){
      rowsum_rho = rowSums(rho_temp)
      row_rho_psi = rowSums(rho_temp*psi)
    }
    
    ## update m 
    for (j in 1:p) {
      covar_m = (rowsum_rho[j]+Sigma[j]^-1)^-1
      mean_m  = covar_m*(rowsum_kappa[j]+row_rho_psi[j]- (w_temp[j,]%*%z_temp)%*%rho_temp[j,])
      m_temp[j] = rnorm(1,mean_m,sqrt(covar_m))
    }
    
    chain_m[t+1,]=as.vector(m_temp)
    
    
    if(dt==0){
      for (j in 1:p) {
        # rho_j
        r = phi[j,1]+sum((X[j,]-(m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])))^2)
        rho_j = rgamma(1,shape=(phi[j,1]+n)/2 ,rate = r/2 )*rep(1,n)
        rho_temp[j,] = rho_j
      }
    }else if(dt==1 | dt==2){
      for (j in 1:p) {
        mu_j = (m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])) 
        for (i in 1:n) {
          rho_temp[j,i] = rpg.devroye(1,b[j,i],mu_j[i])  
        }
      }
    }else if(dt==3){
      if(j<=n_nz){
        r = phi[j,1]+sum((X[j,]-(m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])))^2)
        rho_j = rgamma(1,shape=(phi[j,1]+n)/2 ,rate = r/2 )*rep(1,n)
        rho_temp[j,] = rho_j
      }else{
        mu_j = (m_temp[j]*rep(1,n)+as.vector(t(z_temp)%*%w_temp[j,])) 
        rho_temp[j,] = rpg(n,b[j,],mu_j)
      }
    }
    
    
    
    chain_rho[t+1,]=as.vector(t(rho_temp))
    
    # generate alpha_l 
    
    for (l in 1:L) {
      for (j in 1:p) {
        candi = rnorm(1,alpha_temp[j,l],1)
        comp =  as.numeric(p_density_mmh(candi,alpha_temp[j,l],tau_temp[j,l],alpha_temp[,l],omega_temp,j))   
        if(is.na(comp)){
          comp=0
        }
        # print(comp)
        prob = min(comp,1)
        u = runif(1,0,1)
        #  print(u)
        if(u<prob){
          alpha_temp[j,l]=candi
        }
      }
    }
    
    chain_alpha[t+1,]=as.vector(t(alpha_temp))
    
    
    ## calculate A  
    
    A = eta*(diag(eps,p) + 1) + (alpha_temp-nu_1)%*%t(alpha_temp-nu_1)/nu_2
    
    ## Begin to update Omega
    
    for (j in 1:p) {
      ## j-th diagonal element and j-th column excluding the diag 
      A_jj = A[j,j]     
      a_j = A[-j,j]   
      
      # construct the inverse of Omega_11 via inverse of Omega 
      inv_omega_11 = inv_omega_temp[-j,-j]
      inv_omega_12 = inv_omega_temp[-j,j]
      inv_omega_22 = inv_omega_temp[j,j]
      
      O_11 = inv_omega_11 - inv_omega_12%*%t(inv_omega_12)/inv_omega_22
      
      ## find the index of non-zeros of j-th col 
      nz_ind = which(omega_temp[,j]!=0) # include diagonal
      
      ## remove the diag index: index==j
      nz_ind = nz_ind[nz_ind!=j]        # exclude diagonal
      
      
      ## see if there is any non-zero entry except for the diag
      ## if length(nz_ind)=0: only generate the diag
      if(length(nz_ind)!=0){
        
        ## define the block matrices (mean part)
        ## i) select nz_ind rows and remove the j-th col 
        ## ii) remove the (nz_ind,j) rows and the j-th col 
        ## construct the precision matrix for \omega_12^{(1)}
        
        if(length(nz_ind)>1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          prec_nz_sqrt = chol(prec_nz)  # upper tri matrix 
          mid = rnorm(length(nz_ind),0,1/sqrt(A_jj))
          pt_mean_o1 = (-1/A_jj)*A[nz_ind,j]
          pt_mean_o2 = forwardsolve(t(prec_nz_sqrt),pt_mean_o1)
          mean_nz = backsolve(prec_nz_sqrt,pt_mean_o2)
          non_zero_omega = backsolve(prec_nz_sqrt,mid) + mean_nz 
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          partial = prec_nz_sqrt%*%non_zero_omega
          diag_omega = xi + t(partial)%*%partial
          
        }else if(length(nz_ind)==1){
          prec_nz = inv_omega_temp[nz_ind, nz_ind]-inv_omega_temp[nz_ind,j]%*%t(inv_omega_temp[nz_ind,j])/inv_omega_temp[j,j]
          mean_nz = (-1/(A_jj*prec_nz))*A[nz_ind,j]
          non_zero_omega = rnorm(1,mean_nz, 1/sqrt(prec_nz*A_jj))
          
          ## generate the diag entry
          
          xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
          diag_omega = xi + (non_zero_omega^2)*prec_nz
        }
        
        ## fill in the non-zero entries with the new sample
        ## both col and row d
        omega_temp[nz_ind,j] = non_zero_omega   # col 
        omega_temp[j,nz_ind] = non_zero_omega   # row
        omega_temp[j,j]  = diag_omega
      }else {
        
        ## when each entry = 0 except for the diag  
        ## we only generate the xi with Gamma
        xi = rgamma(1,shape = 1+(eta*(1+eps)+L)/2,rate = A_jj/2)
        diag_omega = xi
        omega_temp[j,j]  = diag_omega
      }
      
      # update inv_omega_temp 
      inv_omega_temp[j,j] = 1/xi
      pp = O_11%*%omega_temp[-j,j]
      inv_omega_11 = O_11 + inv_omega_temp[j,j]*pp%*%t(pp)
      inv_omega_12 = -inv_omega_11%*%omega_temp[-j,j]/omega_temp[j,j]
      
      inv_omega_temp[-j,-j] = inv_omega_11
      inv_omega_temp[-j,j] = inv_omega_12
      inv_omega_temp[j,-j] = inv_omega_12 
      
    }
    
    # generate z_i & tau_z
    
    for (i in 1:n) {
      pt_cov = t(w_temp)%*%(w_temp*rho_temp[,i])
      diag(pt_cov) = diag(pt_cov)+ (1/tauz_temp[,i])*xiz_temp[,i]
      chol_pt_cov = chol(pt_cov) # upper matrix 
      
      pt_mean_z1 = t(w_temp)%*%(rho_temp[,i]*(psi[,i]-m_temp)+kappa[,i])
      pt_mean_z2 = forwardsolve(t(chol_pt_cov),pt_mean_z1)
      mean_z_i  = backsolve(chol_pt_cov,pt_mean_z2)
      z_temp[,i] = backsolve(chol_pt_cov,rnorm(L,0,1)) + mean_z_i
      tauz_temp[,i]=1/rinvgauss(L,abs(1/z_temp[,i]),xiz_temp[,i])
      #print(xiz_temp[,i])
    }
    
    
    chain_z[t+1,]=as.vector(t(z_temp))
    chain_tauz[t+1,]=as.vector(t(tauz_temp))
    
    # generate xi_z
    for (i in 1:n){
      for (l in 1:L){
        xiz_shape = nu_4 + 3/2
        xiz_rate  = nu_3 + (tauz_temp[l,i])/2 + (z_temp[l,i])^2/(2*tauz_temp[l,i]) 
        xiz_temp[l,i] = rgamma(1,shape=xiz_shape,rate=xiz_rate) 
      }
    }
    
    chain_xiz[t+1,]=as.vector(t(xiz_temp))
    
    
    # update w_j & tau
    
    for (j in 1:p) {
      cov_pt = z_temp%*%(t(z_temp)*rho_temp[j,])
      diag(cov_pt)= diag(cov_pt)+ 1/tau_temp[j,]
      chol_cov_pt = chol(cov_pt)
      
      pt_mean_w1 = t(t(z_temp)*rho_temp[j,])%*%(psi[j,]-m_temp[j]*rep(1,n)+(rho_temp[j,]^-1)*kappa[j,])
      pt_mean_w2 = forwardsolve(t(chol_cov_pt),pt_mean_w1)
      mean_w_j = backsolve(chol_cov_pt,pt_mean_w2)
      w_temp[j,] = backsolve(chol_cov_pt,rnorm(L,0,1)) + mean_w_j
      tau_temp[j,]=1/rinvgauss(L,abs(exp(alpha_temp[j,])/w_temp[j,]),exp(2*alpha_temp[j,]))
    }
    chain_w[t+1,]=as.vector(t(w_temp))
    chain_tau[t+1,]=as.vector(t(tau_temp))
    #print(paste0("Task Progress: ", t/T ))
  }
  
  
  # degree of freedom + estimation of W and Z
  
  interval_z = apply(chain_z[start:end,],2,sample_quantile)
  df_z = sum(interval_z[1,]*interval_z[2,]>0)
  z_est = apply(chain_z[start:end,],2,mean)
  zero_index_z = which(interval_z[1,]*interval_z[2,]<0,arr.ind = TRUE)
  z_est[zero_index_z] = 0
  z_est=matrix(z_est,nrow=L ,ncol=n,byrow=TRUE)
  
  
  interval_w = apply(chain_w[start:end,],2,sample_quantile)
  df_w = sum(interval_w[1,]*interval_w[2,]>0)
  w_est= apply(chain_w[start:end,],2,mean)
  zero_index_w = which(interval_w[1,]*interval_w[2,]<0,arr.ind = TRUE)
  w_est[zero_index_w] = 0
  w_est=matrix(w_est,nrow=p,ncol=L,byrow=TRUE)
  
  m_est  = apply(chain_m[start:end,],2,mean)
  interval_m = apply(chain_m[start:end,],2,sample_quantile)  
  zero_index_m = which(interval_m[1,]*interval_m[2,]<0,arr.ind = TRUE)
  m_est[zero_index_m] = 0
  m_est = matrix(m_est,nrow=p ,ncol=1,byrow=T)
  
  mu_est = w_est%*%z_est + as.vector(m_est)
  
  error = norm(mu_est-mu,type = "F")/norm(mu,type = "F")
  
  # sparsity structures of W, Z and m  
  mask_m_w = rep(1,p*L)
  mask_m_w[zero_index_w] = 0
  mask_m_w = matrix(mask_m_w,nrow=p,ncol=L,byrow = TRUE)
  
  mask_m_z = rep(1,L*n)
  mask_m_z[zero_index_z]=0 
  mask_m_z = matrix(mask_m_z,nrow=L,ncol=n,byrow=TRUE)
  
  mask_m_m = rep(1,p)
  mask_m_m[zero_index_m] = 0
  mask_m_m = matrix(mask_m_m,nrow=p,ncol=1,byrow = TRUE)
  
  prec_mat = matrix(apply(chain_rho[start:end,], 2, mean),nrow=p,ncol=n,byrow=T)
  
  # likelihood chain over MC 
  likelihood_chain = c() 
  likelihood_sparse_chain = c()
  
  nllikelihood = matrix(0,nrow=p,ncol=n) # not log
  nllikelihood_sparse = matrix(0,nrow=p,ncol=n) # not log
  
  muchain_est=matrix(0,nrow=p,ncol=n)
  muchain_est_sparse=matrix(0,nrow=p,ncol=n)
  for (t in start:end) {
    prec_mat_t = matrix(chain_rho[t,],nrow=p,ncol=n,byrow=T)
    w_est_t = matrix(chain_w[t,],nrow=p,ncol=L,byrow=TRUE)
    w_est_t_sparse = w_est_t * mask_m_w
    z_est_t = matrix(chain_z[t,],nrow=L ,ncol=n,byrow=TRUE)
    z_est_t_sparse = z_est_t*mask_m_z
    m_est_t = matrix(chain_m[t,],nrow=p ,ncol=1,byrow=T)
    m_est_t_sparse = m_est_t*mask_m_m
    
    mu_est_t = w_est_t%*%z_est_t + as.vector(m_est_t)
    mu_est_t_sparse = w_est_t_sparse %*% z_est_t_sparse + as.vector(m_est_t_sparse)
    
    muchain_est = muchain_est + mu_est_t
    muchain_est_sparse = muchain_est_sparse +mu_est_t_sparse
    if(dt==0){
      # matrix form 
      likelihood_g  = 0.5*log(prec_mat_t)-0.5*log(2*pi)-0.5*prec_mat_t*(X-mu_est_t)^2
      likelihood_g_sparse  = 0.5*log(prec_mat_t)-0.5*log(2*pi)-0.5*prec_mat_t*(X-mu_est_t_sparse)^2
      
      nllikelihood = nllikelihood+sqrt(prec_mat_t/(2*pi))*exp(-(prec_mat_t*(X-mu_est_t)^2)/2)
      nllikelihood_sparse=nllikelihood_sparse+sqrt(prec_mat_t/(2*pi))*exp(-(prec_mat_t*(X-mu_est_t_sparse)^2)/2)
      
    }else if(dt==1|dt==2){
      likelihood_g = log(choose(trials,X))+mu_est_t*X-log(1+exp(mu_est_t))*trials 
      likelihood_g_sparse = log(choose(trials,X))+mu_est_t_sparse*X-log(1+exp(mu_est_t_sparse))*trials
      
      nllikelihood =  nllikelihood+choose(trials,X)*exp(mu_est_t*X)/((1+exp(mu_est_t))^trials)
      nllikelihood_sparse =nllikelihood_sparse+choose(trials,X)*exp(mu_est_t_sparse*X)/((1+exp(mu_est_t_sparse))^trials)
    }
    # sum form 
    likelihood_chain = c(likelihood_chain,sum(likelihood_g))
    likelihood_sparse_chain = c(likelihood_sparse_chain,sum(likelihood_g_sparse))
  }
  
  # posterior est of mu over the MC 
  muchain_est = muchain_est/(end-start+1)
  muchain_est_sparse = muchain_est_sparse/(end-start+1)
  
  # lppd
  lppd = sum(log(nllikelihood/(end-start+1)))
  lppd_sparse =sum(log(nllikelihood_sparse/(end-start+1)))
  
  # scalar 
  mean_like = mean(likelihood_chain)  
  mean_like_sparse= mean(likelihood_sparse_chain)
  
  # BIC 1 and 2
  bic_11 = -2*(mean_like) + log(n)*(df_w+df_z)
  bic_12 = -2*(mean_like_sparse) + log(n)*(df_w+df_z)
  
  # DIC 1 and 2
  likelihood_pointest =sum(0.5*log(prec_mat)-0.5*log(2*pi)-0.5*prec_mat*(X-muchain_est)^2)
  likelihood_pointest_sparse = sum(0.5*log(prec_mat)-0.5*log(2*pi)-0.5*prec_mat*(X-muchain_est_sparse)^2)
  dic_11 = -2*likelihood_pointest + 2*2*(likelihood_pointest-mean_like)
  dic_12 = -2*likelihood_pointest_sparse + 2*2*(likelihood_pointest_sparse-mean_like_sparse)
  
  # WAIC 1 and 2
  waic_11 = -2*lppd + 2*2*(lppd- mean_like)
  waic_12 = -2*lppd_sparse + 2*2*(lppd_sparse-mean_like_sparse)
  
  #list(lppd=lppd,lppd_sparse=lppd_sparse,likelihood_sparse_chain = likelihood_sparse_chain,bic_11=bic_11,bic_12=bic_12,dic_11=dic_11,dic_12=dic_12,waic_11=waic_11,waic_12=waic_12,likelihood_chain=likelihood_chain,error=error,w_est=w_est,z_est=z_est,m_est=m_est,interval_w=interval_w,df_w=df_w,interval_z=interval_z,df_z=df_z,interval_m=interval_m,chain_m=chain_m[start:end,],chain_w=chain_w[start:end,],chain_z=chain_z[start:end,],mu_est=mu_est)
  # save memory 
  list(lppd=lppd,lppd_sparse=lppd_sparse,likelihood_sparse_chain = likelihood_sparse_chain,bic_11=bic_11,bic_12=bic_12,dic_11=dic_11,dic_12=dic_12,waic_11=waic_11,waic_12=waic_12,likelihood_chain=likelihood_chain,error=error,w_est=w_est,z_est=z_est,m_est=m_est,interval_w=interval_w,df_w=df_w,interval_z=interval_z,df_z=df_z,interval_m=interval_m,mu_est=mu_est)
}




