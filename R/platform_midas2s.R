

#' @title An Bayesian platform design without subgroup efficacy exploration(midas-2s), which is the degenerate competing design in the simulation.
#' @description MIDAS-2s is the degenerate competing designs that do not consider subgroups. Beta-binomial model is applied for efficacy in whole population of each arm.
#' @param seed set a random seed to maintain the repeatability of the simulation results.
#' @param p a matrix indicating the efficacy. Row number represents the number of candidate drugs.
#' @param p_tox a vector indicating the toxicity.
#' @param C_T early toxicity stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to unacceptable levels of toxicity or adverse events in the study participants. This threshold is established to ensure the safety and well-being of the trial participants and to prevent further harm.
#' @param C_E1 early futility stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to lack of efficacy or futility. It is established to prevent the continuation of a trial that is unlikely to demonstrate a significant treatment effect, thus saving time, resources, and participant exposure to ineffective treatments.
#' @param C_E2 early efficacy stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to the demonstration of significant efficacy or positive treatment effects. This threshold is established to allow for timely decision-making and saves sample size.
#'
#' @return term.tox the indicator of whether early stopping for toxicity
#' @return term.fut the indicator of whether early stopping for futility
#' @return term.eff the indicator of whether early stopping for efficacy
#' @return final.eff a vector of final decision, either efficacy or inefficacy
#' @return N sample size, which refers to the number of participants included in a study or experiment.
#' @export
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
#'   result <- platform_midas2s(seed=20,p,p_tox,C_T=0.85,C_E1=0.15,C_E2=0.999)
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
#'   result <- platform_midas2s(seed=12,p,p_tox,C_T=0.85,C_E1=0.15,C_E2=0.999)
#'   result
#'   }
#'
#'




platform_midas2s <- function(seed,p,p_tox,C_T=0.85,C_E1=0.15,C_E2=0.999){
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
  set.seed(seed+100)
  prevalence <- c(0.275, 0.135, 0.225, 0.365)
  p <- p%*%prevalence
  tox.post <- rep(0,nrow(p))
  eff.post <- rep(0,nrow(p))
  Nmin <- 15  ##每个臂都至少分配20人再进行AR
  Nmax <- 75  ##每个臂预设最大样本量为100
  finish <- rep(0,nrow(p))        # 记录试验最后结果，包括对照组
  term.fut <- rep(0, nrow(p)-1 )  # 记录早期无效停止，不包括对照组
  term.eff <- rep(0, nrow(p)-1 )  # 记录早期有效停止
  term.tox <- rep(0, nrow(p)-1 )  # 记录早期过毒停止
  final.eff <-  rep(0, nrow(p)-1 )  # 最后有效性判断


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

    Trt <- rep(index-1,  N[index])
    y <- rep(NA,  length(Trt))
    tox <- rep(NA,  length(Trt))
    data <- data.frame(Trt, y, tox )

    for (i in index){
      data$y[Trt==i-1] <- stats::rbinom(length(data$y[Trt==i-1]),1,p[i])
      data$tox[Trt==i-1] <- stats::rbinom(length(data$tox[Trt==i-1]),1,p_tox[i])
    }

    posteff <- function(c1,c2,t1,t2){
      f <- function(x) { stats::pbeta(x,c1+1,c2+1)*stats::dbeta(x, t1+1, t2+1) }
      prob <- stats::integrate(f, 0, 1, rel.tol=1e-4)$value
      return(prob)
    }
    for(k in 1:length(index)){
      eff.post[index[k]] <-  (sum(data$y[Trt==index[k]-1])+1)/(length(data$y[Trt==index[k]-1])+2)
    }

    alpha[index] <- round(eff.post/(sum(eff.post)),3)[index]

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
            eff.fina <-  posteff(sum(data$y[Trt==0]),length(data$y[Trt==0])-sum(data$y[Trt==0]),
                                 sum(data$y[Trt==index[l]-1]),
                                 length(data$y[Trt==index[l]-1])-sum(data$y[Trt==index[l]-1]))
            futility <- ifelse( eff.fina < C_E1, 1, 0)
            effgrad <- ifelse(eff.fina > C_E2, 1, 0)
            if(futility == 1) { term.fut[index[l]-1] <- 1 }
            if(effgrad == 1) { term.eff[index[l]-1] <- 1 }
            if( (N[index[l]] >= Nmax || futility == 1 || effgrad == 1) ){        ##判断是否达到了最大样本量
              finish[index[l]] <- 1
              final.eff[index[l]-1] <- eff.fina
            }
          }
        } #if(l>1)
      }##if(n[l]>=10)
    }##for(l in index)
    interim_num <- interim_num + 1
  } ## while 主循环结束，试验结束，所有药物都被test过，返回结果
  re <- list(term.tox=term.tox,term.fut=term.fut,term.eff=term.eff,final.eff=final.eff,N=N)
  return(re)
}








