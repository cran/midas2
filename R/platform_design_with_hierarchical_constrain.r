

#' Title
#' @title An Information borrowing drug-combination Bayesian platform design(midas-2)
#' @param seed set a random seed to maintain the repeatability of the simulation results.
#' @param p a matrix indicating the efficacy. Row number represents the number of candidate drugs.
#' @param p_tox a vector indicating the toxicity.
#'
#' @return term.tox the indicator of whether early stopping for toxicity
#' @return term.fut the indicator of whether early stopping for futility
#' @return term.eff the indicator of whether early stopping for efficacy
#' @return final.eff a vector of final decision, either efficacy or inefficacy
#' @return post.subg subgroup analysis for treatments
#' @return post.sign signature analysis for treatments
#' @return post.spike posterior estimation for spike parameters
#' @return best selection of best treatment for each subgroup
#' @export
#'
#' @examples
#'


hc_platform <- function(seed,p,p_tox){
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
  # estimation
  estimate <- function(data,index,seed){

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
    slab1_1 ~ dnorm(0, 1)
    beta1  <- spike1 * slab1_1

  "
    for(i in indx){
      prior_expr <- paste(prior_expr, paste("kexi1_",i," ~ dnorm(gamma1, tau1)  \n ","gamma1_",i," <- beta1*kexi1_",i,"  \n  ",sep= ""))
    }
    prior_expr <- paste(prior_expr,"
      prob2 ~ dunif(0,1)
      spike2 ~ dbern(prob2)
      slab2_1 ~ dnorm(1, 1)
      beta2  <- spike2 * slab2_1"," \n  ",sep= "")
    for(i in indx){
      prior_expr <- paste(prior_expr, paste("kexi2_",i," ~ dnorm(gamma2, tau2)  \n ","gamma2_",i," <- beta2*kexi2_",i,"  \n  ",sep= ""))
    }
    for(i in indx){
      prior_expr <- paste(prior_expr, paste("theta",i," ~ dnorm(0, 0.001)","  \n  ",sep= ""))
    }
    prior_expr <- paste(prior_expr,"
      gamma1 ~ dnorm(0, 0.001)
      gamma2 ~ dnorm(0, 0.001)
      tau1 ~ dgamma(0.001, 0.001)
      tau2 ~ dgamma(0.001, 0.001)"," \n  ",sep= "")

    model_text <- paste(regre,prior_expr,"}", sep = " \n ")
    modstr = textConnection(model_text)


    para <- c("beta0", "beta1", "beta2",  paste("theta",indx,sep= ""), paste("gamma1_",indx,sep= ""), paste("gamma2_",indx,sep= ""),"spike1","spike2")

    output <- jags(data = jags.data, model.file=modstr, parameters.to.save=para, n.burnin = 10000,
                   n.chains = 1, n.iter = 20000, n.thin = 2, set.seed(seed), jags.module = c("dic"))


    beta0 <- as.numeric( as.mcmc(output)[,"beta0"][[1]])
    beta1 <-  as.numeric( as.mcmc(output)[,"beta1"][[1]])
    beta2 <-  as.numeric( as.mcmc(output)[,"beta2"][[1]])

    theta <- matrix(NA,5000,length(indx))
    for(i in 1:length(indx)){
      theta[,i] <- as.numeric( as.mcmc(output)[,paste("theta",indx[i],sep="")][[1]])
    }
    theta <- cbind(rep(0,5000),theta)


    spike1 <-  as.numeric( as.mcmc(output)[,"spike1"][[1]])
    spike2 <-  as.numeric( as.mcmc(output)[,"spike2"][[1]])

    gamma1 <- matrix(NA,5000,length(indx))
    gamma2 <- matrix(NA,5000,length(indx))
    for(i in 1:length(indx)){
      gamma1[,i] <- as.numeric( as.mcmc(output)[,paste("gamma1_",indx[i],sep="")][[1]])
      gamma2[,i] <- as.numeric( as.mcmc(output)[,paste("gamma2_",indx[i],sep="")][[1]])
    }
    gamma1 <- cbind(rep(0,5000),gamma1)
    gamma2 <- cbind(rep(0,5000),gamma2)


    p <- array(NA,dim=c(length(indx)+1,4,5000))
    for(i in 1:(length(indx)+1)){
      p[i,1,] <-  exp(beta0 + beta1 + beta2 + theta[,i] + gamma1[,i] + gamma2[,i])/(1+exp(beta0 + beta1 + beta2 + theta[,i] + gamma1[,i] + gamma2[,i]))   ###  type 1
      p[i,2,] <-  exp(beta0 + beta1 + theta[,i] + gamma1[,i])/(1+exp(beta0 + beta1 + theta[,i] + gamma1[,i]))   ###  type 2
      p[i,3,] <-  exp(beta0 + beta2 + theta[,i] + gamma2[,i])/(1+exp(beta0 + beta2 + theta[,i] + gamma2[,i]))   ###  type 3
      p[i,4,] <-  exp(beta0 + theta[,i])/(1+exp(beta0 + theta[,i]))   ###  type 4
    }
    p[is.na(p)] <- runif(1,0.1,0.9)


    phi1 = phi2 = phi3 = phi4 = 0.1  #Prior
    phi <- c(sum(data$Type==1) + phi1, sum(data$Type==2) + phi2, sum(data$Type==3) + phi3, sum(data$Type==4) + phi4)  #Posterior
    set.seed(seed)
    di <- rdirichlet(5000, phi)

    p_d <- matrix(NA,nrow=5000,ncol=length(indx)+1)
    for(i in 1:(length(indx)+1)){
      p_d[,i] <-  di[,1] *  p[i,1,] + di[,2] * p[i,2,] + di[,3] *  p[i,3,] +  di[,4] * p[i,4,]
    }

    ## 每个臂的后验有效率
    p_arm <- rep(NA,length(indx))
    for(i in 1:length(indx)){
      p_arm[i] <- mean( p_d[,1+i] > p_d[,1])
    }

    p_s <-  array(NA,dim=c(length(indx)+1,4,5000))
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

    ## beta1 beta2 被估计为0的概率
    prob_spike <- c(mean(spike1),mean(spike2))
    #para.est <- colMeans(matrix(unlist(as.mcmc(output)),nrow=5000,byrow=FALSE))
    re <- list(arm=p_arm,  eff.est.sub=eff.est.sub,  subgroup=subgroup,
               eff.est.sig=eff.est.sig,  signature=signature,  prob_spike=prob_spike)
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
  post.spike <- NULL
  C_T <- 0.85   ## 毒性停止阈值
  C_E1 <- 0.150   ## 无效停止阈值
  C_E2 <- 0.999   ## 早期有效停止阈值
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
    t <- c(t, runif(sum(rand$n_arm),  100 * interim_num , (interim_num + 1) * 100 ))
    data <- data.frame(Trt, Type, x1, x2, y, tox, t)

    ## 整除 %/%， e.g. 82 %/% 10 = 8
    for (i in index){
      data$tox[Trt==i-1] <- rbinom(length(data$tox[Trt==i-1]),1,p_tox[i])
      for (j in 1:4){
        data$y[Trt==i-1 & Type==j] <- rbinom(length(data$y[Trt==i-1 & Type==j]), 1, p[i,j])
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
    temp <- estimate(data_est, index, seed)
    alpha[index] <- apply(temp$eff.est.sub,1,mean)/(sum(apply(temp$eff.est.sub,1,mean)))

    for(l in 1:length(index)){
        if(N[index[l]]>=Nmin){  ## 判断是否大于最小样本量
            if(l>1){   ## 只对试验药进行分析，drug 1对照药永远在试验中
              tox.post[index[l]] <- 1 - pbeta(0.3,1+sum(data$tox[Trt==index[l]-1]),1+length(data$tox[Trt==index[l]-1])-sum(data$tox[Trt==index[l]-1]))
              tox_ind <- ifelse(tox.post[index[l]] > C_T, 1, 0)  ## Pr(p_T>0.3)>CT，则过毒停止
              if(tox_ind==1){  ## 判断是否因毒性停止
                  finish[index[l]] <- 1
                  term.tox[index[l]-1] <- 1  #有效性指标咋办？
                  post.spike <- rbind(post.spike,temp$prob_spike)
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
                    post.spike <- rbind(post.spike,temp$prob_spike)
                  }
              }
            } #if(l>1)
        }##if(n[l]>=10)
    }##for(l in index)
    interim_num <- interim_num + 1
  } ## while 主循环结束，试验结束，所有药物都被test过，返回结果

  best <-  sapply(1:4,function(i) ifelse(post.subg[,i]==max(post.subg[,i]),1,0))

  re <- list(term.tox=term.tox,term.fut=term.fut,term.eff=term.eff,final.eff=final.eff,N=N,
             post.subg=post.subg,post.sign=post.sign,post.spike=post.spike,best=best)
  return(re)
}









