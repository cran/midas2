



#' @title An Information borrowing Bayesian platform design with subgroup efficacy exploration(midas-2)
#' @description
#' The rapid screening of effective and optimal therapies from large numbers of candidate combinations, as well as exploring subgroup efficacy, remains challenging, which necessitates innovative, integrated, and efficient trial designs.
#' MIDAS-2 package enables quick and continuous screening of promising combination strategies and exploration of their subgroup effects within a unified platform design framework. We used a regression model to characterize the efficacy pattern in subgroups. Information borrowing was applied through Bayesian hierarchical model to improve trial efficiency considering the limited sample size in subgroups. MIDAS-2 provides an adaptive drug screening and subgroup exploring framework to accelerate immunotherapy development in an efficient, accurate, and integrated fashion.
#' @param seed set a random seed to maintain the repeatability of the simulation results.
#' @param p a matrix indicating the efficacy. Row number represents the number of candidate drugs.
#' @param p_tox a vector indicating the toxicity.
#' @param n.burnin the number of iterations in burn-in period, which means the first \emph{n.burnin} iterations are discarded.
#' @param n.iter the number of posterior iterations.
#' @param n.thin every \emph{n.thin} iterations after the burn-in period is retained for analysis.
#' @param C_T early toxicity stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to unacceptable levels of toxicity or adverse events in the study participants. This threshold is established to ensure the safety and well-being of the trial participants and to prevent further harm.
#' @param C_E1 early futility stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to lack of efficacy or futility. It is established to prevent the continuation of a trial that is unlikely to demonstrate a significant treatment effect, thus saving time, resources, and participant exposure to ineffective treatments.
#' @param C_E2 early efficacy stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to the demonstration of significant efficacy or positive treatment effects. This threshold is established to allow for timely decision-making and saves sample size.
#'
#' @return term.tox the indicator of whether early stopping for toxicity
#' @return term.fut the indicator of whether early stopping for futility
#' @return term.eff the indicator of whether early stopping for efficacy
#' @return final.eff a vector of final decision, either efficacy or inefficacy
#' @return post.subg subgroup analysis for treatments
#' @return post.sign signature analysis for treatments
#' @return best selection of best treatment for each subgroup
#' @return N sample size, which refers to the number of participants included in a study or experiment.
#' @export
#' @details The sample size for a particular subgroup is limited. Therefore, it is difficult to accurately estimate treatment-related effects in each experiment arm separately, and the corresponding subgroup power may also be insufficient. Given that all candidate combination therapies in each arm consist of investigational drugs, it is justifiable to leverage information from specific subgroups across arms. Specifically, we adopt the Bayesian hierarchical model and assign half-cauchy prior distributions to the the standard deviation parameters.
#'
#' @examples
#'
#'
#'   # Example 1
#'   p0 <- c(0.1,0.1,0.1,0.1)
#'   p1 <- c(0.1,0.1,0.1,0.1)
#'
#'   p <- rbind(p0,p1)
#'   p_tox <- c(0.1,0.4)
#'
#'   # consider 1 candidate drugs with 4 subgroups
#'   result <- platform_midas2(seed=20,p,p_tox,n.burnin=1000,
#'             n.iter=2000,n.thin=2,C_T=0.85,C_E1=0.15,C_E2=0.999)
#'   result
#'
#'
#'   \donttest{
#'   # Example 2
#'   p0 <- c(0.05,0.10,0.05,0.10)
#'   p1 <- c(0.24,0.40,0.12,0.22)
#'   p2 <- c(0.24,0.40,0.12,0.22)
#'   p3 <- c(0.12,0.22,0.05,0.10)
#'   p4 <- c(0.24,0.40,0.12,0.22)
#'   p5 <- c(0.28,0.45,0.12,0.22)
#'   p6 <- c(0.24,0.40,0.12,0.22)
#'   p7 <- c(0.12,0.22,0.05,0.10)
#'
#'   p <- rbind(p0, p1, p2, p3, p4, p5, p6, p7)
#'   p_tox <- c(0.10,0.10,0.10,0.10,0.10,0.10,0.15,0.20)
#'
#'   # consider 7 candidate drugs with 4 subgroups
#'   result <- platform_midas2(seed=12,p,p_tox,n.burnin=10000,
#'             n.iter=20000,n.thin=2,C_T=0.85,C_E1=0.15,C_E2=0.999)
#'   result
#'   }
#'



platform_midas2 <- function(seed,p,p_tox,n.burnin=10000,n.iter=20000,n.thin=2,C_T=0.85,C_E1=0.15,C_E2=0.999){

  random2 <- function(alpha, n_arm, cohortsize){
    Nmin <- 15
    Nmax <- 75
    pi_arm <- rep(0, length(n_arm));
    c <- ifelse(n_arm[1:length(n_arm)] >= Nmin, 1, 0);  ##从第2个臂开始，人数是否大于N_min，是=1，否=0
    s <- sum(c)
    if(s == 0 & n_arm[1] < Nmin){
      pi_arm <- rep(1/length(n_arm), length(n_arm));
    }
    if(s >= 1){
      pi_arm[c==0] <- pmax(alpha[c==0]/sum(alpha), rep(1/length(n_arm),length(alpha[c==0])));   ## 向量中小于1/k的都用1/k代替 《追赶规则》
      if(n_arm[1] < (0.35*Nmax)){    #beta=0.35,Nmax=40    ##  对照组被分配的患者数是否大于beta*N_max
        pi_arm[1] <- (1-sum(pi_arm[c==0])) * (alpha[1] / (alpha[1] + sum(alpha[c==1])));  ##  未大于beta*N_max，按比例随机化
        pi_arm[c==1] <- (1-sum(pi_arm[c==0])-pi_arm[1]) * (alpha[c==1] / sum(alpha[c==1]));
      }else{                                                                  ##  已大于beta*N_max，按比例削减规则
        pi_arm[1] <- (1-sum(pi_arm[c==0])) * (alpha[1] / (alpha[1]+sum(alpha[c==1]))) * (0.35*Nmax/n_arm[1])^2;
        pi_arm[2:length(n_arm)][(c==1)[2:length(n_arm)]] <- (1-sum(pi_arm[c==0])-pi_arm[1]) * (alpha[2:length(n_arm)][(c==1)[2:length(n_arm)]] / sum(alpha[2:length(n_arm)][(c==1)[2:length(n_arm)]]))
      }
    }
    pi_arm[pi_arm < 0.1] <- 0.1;
    pi_arm[pi_arm > 0.9] <- 0.9;
    pi_arm <- pi_arm/sum(pi_arm);
    n_arm <-  round(cohortsize * pi_arm);
    return(list(n_arm = n_arm, pi_arm = pi_arm));
  }
  getprior <- function(data,index,seed,n.burnin,n.iter,n.thin){

    jags.data <- list("y" = data$y, "x1" = data$x1, "x2" = data$x2, "Trt" = data$Trt, "N" = nrow(data))
    ###  jags model
    #index <- seq(1:5)
    indx <- index[-1]-1
    drug_vec <- rep(FALSE,7)
    drug_vec[indx] <- TRUE


    regre <- "
       model {
       for (i in 1:N) {
       y[i] ~ dbern(p[i])
       logit(p[i]) <- beta0 + beta1*x1[i] +  beta2*x2[i] + "
    for(i in indx){
      if (i == indx[length(indx)]){
        regre <- paste(regre,paste("equals(Trt[i],",i,")*theta",i," + equals(Trt[i],",i,")*gamma1_",i,"*x1[i] + equals(Trt[i],",i,")*gamma2_",i,"*x2[i]}", sep=""))
      }else{
        regre <- paste(regre,paste("equals(Trt[i],",i,")*theta",i," + equals(Trt[i],",i,")*gamma1_",i,"*x1[i] + equals(Trt[i],",i,")*gamma2_",i,"*x2[i] + ", sep=""))
      }
    }

    prior_expr <- "
       beta0 ~ dnorm(-2.2, 1)

       prob1 ~ dunif(0,1)
       spike1 ~ dbern(prob1)
       slab1_1 ~ dnorm(0, 0.01)
       beta1  <- spike1 * slab1_1

       prob2 ~ dunif(0,1)
       spike2 ~ dbern(prob2)
       slab2_1 ~ dnorm(0, 0.01)
       beta2  <- spike2 * slab2_1
     "

    for(i in indx){
      prior_expr <- paste(prior_expr, paste("theta",i," ~ dnorm(0, 0.001)","  \n  ",sep= ""))
    }

    for(i in indx){
      prior_expr <- paste(prior_expr, paste("gamma1_",i," ~ dnorm(0, 0.001)","  \n  ",sep= ""))
      prior_expr <- paste(prior_expr, paste("gamma2_",i," ~ dnorm(0, 0.001)","  \n  ",sep= ""))
    }


    model_text <- paste(regre,prior_expr,"}", sep = " \n ")
    modstr = textConnection(model_text)


    para <- c("beta0", "beta1", "beta2",  paste("theta",indx,sep= ""), paste("gamma1_",indx,sep= ""), paste("gamma2_",indx,sep= ""))

    output <- R2jags::jags(data = jags.data, model.file=modstr, parameters.to.save=para, n.burnin = n.burnin,
                   n.chains = 1, n.iter =  n.iter, n.thin = n.thin, set.seed(seed), jags.module = c("dic"))

    beta1 <-  as.numeric( coda::as.mcmc(output)[,"beta1"][[1]])
    beta2 <-  as.numeric( coda::as.mcmc(output)[,"beta2"][[1]])

    p.b <- c(stats::median(beta1),stats::median(beta2))
    return(p.b)
  }
  # estimation
  estimate <- function(data,index,seed, p.beta,n.burnin,n.iter,n.thin){

    jags.data <- list("y" = data$y, "x1" = data$x1, "x2" = data$x2, "Trt" = data$Trt, "N" = nrow(data),
                      "beta1"=p.beta[1],"beta2"=p.beta[2])
    ###  jags model
    #index <- seq(1:5)
    indx <- index[-1]-1
    drug_vec <- rep(FALSE,7)
    drug_vec[indx] <- TRUE


    regre <- "
       model {
       for (i in 1:N) {
       y[i] ~ dbern(p[i])
       logit(p[i]) <- beta0 + beta1*x1[i] +  beta2*x2[i] + "
    for(i in indx){
      if (i == indx[length(indx)]){
        regre <- paste(regre,paste("equals(Trt[i],",i,")*theta",i," + equals(Trt[i],",i,")*gamma1_",i,"*x1[i] + equals(Trt[i],",i,")*gamma2_",i,"*x2[i]}", sep=""))
      }else{
        regre <- paste(regre,paste("equals(Trt[i],",i,")*theta",i," + equals(Trt[i],",i,")*gamma1_",i,"*x1[i] + equals(Trt[i],",i,")*gamma2_",i,"*x2[i] + ", sep=""))
      }
    }

    prior_expr <- "
        beta0 ~ dnorm(-2.2, 1)
      "

    for(i in indx){
      prior_expr <- paste(prior_expr, paste("gamma1_",i,"  ~ dnorm(gamma1, vtau1)   \n  ",sep= ""))
    }

    for(i in indx){
      prior_expr <- paste(prior_expr, paste("gamma2_",i,"  ~ dnorm(gamma2, vtau2)   \n  ",sep= ""))
    }

    for(i in indx){
      prior_expr <- paste(prior_expr, paste("theta",i," ~ dnorm(mu, vta)","  \n  ",sep= ""))
    }


    prior_expr <- paste(prior_expr,"
      gamma1 ~ dnorm(0, 0.001)
      gamma2 ~ dnorm(0, 0.001)
      mu ~ dnorm(0, 0.001)
      vtau1 <- 1/(tau1*tau1)
      vtau2 <- 1/(tau2*tau2)
      vta <- 1/(ta*ta)
      tau1 ~ dt(0,0.0016,1)I(0,)
      tau2 ~ dt(0,0.0016,1)I(0,)
      ta ~ dt(0,0.0016,1)I(0,)  "," \n  ",sep= "")

    model_text <- paste(regre,prior_expr,"}", sep = " \n ")
    modstr = textConnection(model_text)


    para <- c("beta0", "beta1", "beta2",  paste("theta",indx,sep= ""), paste("gamma1_",indx,sep= ""), paste("gamma2_",indx,sep= ""))

    output <- R2jags::jags(data = jags.data, model.file=modstr, parameters.to.save=para, n.burnin = n.burnin,
                   n.chains = 1, n.iter =  n.iter, n.thin = n.thin, set.seed(seed), jags.module = c("dic"))


    beta0 <- as.numeric( coda::as.mcmc(output)[,"beta0"][[1]])
    beta1 <-  as.numeric( coda::as.mcmc(output)[,"beta1"][[1]])
    beta2 <-  as.numeric( coda::as.mcmc(output)[,"beta2"][[1]])

    theta <- matrix(NA,round((n.iter-n.burnin)/n.thin),length(indx))
    for(i in 1:length(indx)){
      theta[,i] <- as.numeric( coda::as.mcmc(output)[,paste("theta",indx[i],sep="")][[1]])
    }
    theta <- cbind(rep(0,round((n.iter-n.burnin)/n.thin)),theta)

    gamma1 <- matrix(NA,round((n.iter-n.burnin)/n.thin),length(indx))
    gamma2 <- matrix(NA,round((n.iter-n.burnin)/n.thin),length(indx))
    for(i in 1:length(indx)){
      gamma1[,i] <- as.numeric( coda::as.mcmc(output)[,paste("gamma1_",indx[i],sep="")][[1]])
      gamma2[,i] <- as.numeric( coda::as.mcmc(output)[,paste("gamma2_",indx[i],sep="")][[1]])
    }
    gamma1 <- cbind(rep(0,round((n.iter-n.burnin)/n.thin)),gamma1)
    gamma2 <- cbind(rep(0,round((n.iter-n.burnin)/n.thin)),gamma2)


    p <- array(NA,dim=c(length(indx)+1,4,round((n.iter-n.burnin)/n.thin)))
    for(i in 1:(length(indx)+1)){
      p[i,1,] <-  exp(beta0 + beta1 + beta2 + theta[,i] + gamma1[,i] + gamma2[,i])/(1+exp(beta0 + beta1 + beta2 + theta[,i] + gamma1[,i] + gamma2[,i]))   ###  type 1
      p[i,2,] <-  exp(beta0 + beta1 + theta[,i] + gamma1[,i])/(1+exp(beta0 + beta1 + theta[,i] + gamma1[,i]))   ###  type 2
      p[i,3,] <-  exp(beta0 + beta2 + theta[,i] + gamma2[,i])/(1+exp(beta0 + beta2 + theta[,i] + gamma2[,i]))   ###  type 3
      p[i,4,] <-  exp(beta0 + theta[,i])/(1+exp(beta0 + theta[,i]))   ###  type 4
    }
    p[is.na(p)] <- stats::runif(1,0.1,0.9)


    phi1 = phi2 = phi3 = phi4 = 0.1  #Prior
    phi <- c(sum(data$Type==1) + phi1, sum(data$Type==2) + phi2, sum(data$Type==3) + phi3, sum(data$Type==4) + phi4)  #Posterior
    set.seed(seed)
    di <- MCMCpack::rdirichlet(round((n.iter-n.burnin)/n.thin), phi)

    p_d <- matrix(NA,nrow=round((n.iter-n.burnin)/n.thin),ncol=length(indx)+1)
    for(i in 1:(length(indx)+1)){
      p_d[,i] <-  di[,1] *  p[i,1,] + di[,2] * p[i,2,] + di[,3] *  p[i,3,] +  di[,4] * p[i,4,]
    }

    ## 每个臂的后验有效率
    p_arm <- rep(NA,length(indx))
    for(i in 1:length(indx)){
      p_arm[i] <- mean( p_d[,1+i] > p_d[,1])
    }

    p_s <-  array(NA,dim=c(length(indx)+1,4,round((n.iter-n.burnin)/n.thin)))
    for(i in 1:(length(indx)+1)){
      p_s[i,1,] <-   (di[,1] * p[i,1,] + di[,2] * p[i,2,])/(di[,1] + di[,2])     #  signature 1
      p_s[i,2,] <-   (di[,3] * p[i,3,] + di[,4] * p[i,4,])/(di[,3] + di[,4])     #  signature 2
      p_s[i,3,] <-   (di[,1] * p[i,1,] + di[,3] * p[i,3,])/(di[,1] + di[,3])     #  signature 3
      p_s[i,4,] <-   (di[,2] * p[i,2,] + di[,4] * p[i,4,])/(di[,2] + di[,4])     #  signature 4
    }


    ############   亚组分析   ##########
    eff.est.sub <- matrix(NA,length(indx)+1,4,dimnames = list(c("control",paste("arm",indx,sep="")),c("subgroup1","subgroup2","subgroup3","subgroup4")))
    for(i in 1:(length(indx)+1)){
      for(j in 1:4){
        eff.est.sub[i,j] <- mean(p[i,j,])
      }
    }

    subgroup <- matrix(NA,length(indx),4,dimnames = list(paste("arm",indx,sep=""),c("subgroup1","subgroup2","subgroup3","subgroup4")))
    for(i in 1:(length(indx))){
      for(j in 1:4){
        subgroup[i,j] <- mean(p[i+1,j,]>p[1,j,])
      }
    }
    ############   signature分析   ##########
    eff.est.sig <- matrix(NA,length(indx)+1,4,dimnames = list(c("control",paste("arm",indx,sep="")),c("signature1","signature2","signature3","signature4")))
    for(i in 1:(length(indx)+1)){
      for(j in 1:4){
        eff.est.sig[i,j] <- mean(p_s[i,j,])
      }
    }

    signature <- matrix(NA,length(indx),4,dimnames = list(paste("arm",indx,sep=""),c("signature1","signature2","signature3","signature4")))
    for(i in 1:(length(indx))){
      for(j in 1:4){
        signature[i,j] <- mean(p_s[i+1,j,]>p_s[1,j,])
      }
    }

    re <- list(arm=p_arm,  eff.est.sub=eff.est.sub,  subgroup=subgroup,
               eff.est.sig=eff.est.sig,  signature=signature)
    return(re)

  }


  set.seed(seed+100)


  tox.post <- rep(0,nrow(p))
  Nmin <- 15  ##每个臂都至少分配20人再进行AR
  Nmax <- 75  ##每个臂预设最大样本量为100
  finish <- rep(0,nrow(p))        # 记录试验最后结果，包括对照组
  term.fut <- rep(0, nrow(p)-1 )  # 记录早期无效停止，不包括对照组
  term.eff <- rep(0, nrow(p)-1 )  # 记录早期有效停止
  term.tox <- rep(0, nrow(p)-1 )  # 记录早期过毒停止
  final.eff <-  rep(0, nrow(p)-1 )  # 最后有效性判断
  post.subg <- matrix(0, nrow(p)-1, 4)
  post.sign <- matrix(0, nrow(p)-1, 4)

  #cohortsize <- 100  ##最大5个臂同时运行
  alpha <- rep(0,nrow(p))
  N <- rep(0,nrow(p))
  Trt <- Type <- y <- tox <- t <- NULL
  interim_num <- 0
  ##  主循环
  while( sum(finish) <= nrow(p)-2 ){     ##当还有药物还没结束试验时
    ind <- which(finish == 0)
    index <- if(length(ind) >= 5) ind[1:5] else ind   ##一次只取5个药物同时进行试验
    cohortsize <- length(index)*(15-round((7/length(index))^1.1))   ## 期中分析时人数 = 臂数 * 20
    rand <- random2(alpha[index], N[index], cohortsize)

    N[index]  <- N[index] + rand$n_arm
    pi[index] <- rand$pi_arm
    ##  入组患者产生 data

    Trt <- c(Trt, rep(index-1, rand$n_arm))
    Type <- c(Type, sample(1:4, sum(rand$n_arm), replace = T, prob = c(0.275, 0.135, 0.225, 0.365)))
    x1 <- as.numeric(Type == 1 | Type == 2)  ## PD-L1 Low = 0 ; PD-L1 High = 1
    x2 <- as.numeric(Type == 1 | Type == 3)  ## Non-Squamous = 0 ; Squamous = 1
    y <- rep(NA,  length(Trt))
    tox <- rep(NA,  length(Trt))
    t <- c(t, stats::runif(sum(rand$n_arm),  100 * interim_num , (interim_num + 1) * 100 ))
    data <- data.frame(Trt, Type, x1, x2, y, tox, t)

    ## 整除 %/%， e.g. 82 %/% 10 = 8
    for (i in index){
      data$tox[Trt==i-1] <- stats::rbinom(length(data$tox[Trt==i-1]),1,p_tox[i])
      for (j in 1:4){
        data$y[Trt==i-1 & Type==j] <- stats::rbinom(length(data$y[Trt==i-1 & Type==j]), 1, p[i,j])
      }
    }
    # empi_p <- matrix(NA,5,4)
    # for (i in 1:5){
    #   for (j in 1:4){
    #     empi_p[i,j] <- round(mean(data$y[Trt==(i-1) & Type==j]),2)
    #   }
    # }
    # empi_p

    data_est <- data[data$Trt  %in% (index-1),]

    p.beta <- getprior( data_est,index, seed,n.burnin,n.iter,n.thin)

    temp <- estimate(data_est, index, seed, p.beta,n.burnin,n.iter,n.thin)
    alpha[index] <- apply(temp$eff.est.sub,1,mean)/(sum(apply(temp$eff.est.sub,1,mean)))

    for(l in 1:length(index)){
      if(N[index[l]]>=Nmin){  ## 判断是否大于最小样本量
        if(l>1){   ## 只对试验药进行分析，drug 1对照药永远在试验中
          tox.post[index[l]] <- 1 - stats::pbeta(0.3,1+sum(data$tox[Trt==index[l]-1]),1+length(data$tox[Trt==index[l]-1])-sum(data$tox[Trt==index[l]-1]))
          tox_ind <- ifelse(tox.post[index[l]] > C_T, 1, 0)  ## Pr(p_T>0.3)>CT，则过毒停止
          if(tox_ind==1){  ## 判断是否因毒性停止
            finish[index[l]] <- 1
            term.tox[index[l]-1] <- 1  #有效性指标咋办？

          }
          else{    ## 如果不是，进入有效性判断
            futility <- ifelse(temp$arm[l-1] < C_E1, 1, 0)
            effgrad <- ifelse(temp$arm[l-1]  > C_E2, 1, 0)
            if(futility == 1) { term.fut[index[l]-1] <- 1 }
            if(effgrad == 1) { term.eff[index[l]-1] <- 1 }
            if( (N[index[l]] >= Nmax || futility == 1 || effgrad == 1) ){        ##判断是否达到了最大样本量
              finish[index[l]] <- 1
              final.eff[index[l]-1] <- temp$arm[l-1]
              post.subg[index[l]-1,] <- temp$subgroup[l-1,]
              post.sign[index[l]-1,] <- temp$signature[l-1,]

            }
          }
        } #if(l>1)
      }##if(n[l]>=10)
    }##for(l in index)
    interim_num <- interim_num + 1
  } ## while 主循环结束，试验结束，所有药物都被test过，返回结果

  best <-  sapply(1:4,function(i) ifelse(post.subg[,i]==max(post.subg[,i]),1,0))

  re <- list(term.tox=term.tox,term.fut=term.fut,term.eff=term.eff,final.eff=final.eff,N=N,
             post.subg=post.subg,post.sign=post.sign,best=best)
  return(re)
}





