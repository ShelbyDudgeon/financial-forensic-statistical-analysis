# plots #
library(ggplot2)
require(gridExtra)
library(dplyr)
library(AnomDetct)
library(parallel)
library(POT)
library(rmutil)

set.seed(112)
gen_gam <- rggamma(n = 6276, s = 1.5, m = 40, f = 1)
hist(gen_gam)

v1 <- runif(100, 65, 66)
v2 <- runif(50, 130, 131)
v3 <- runif(25, 195, 196)


dat_clust <- c(gen_gam, v1, v2, v3)
df <- as.data.frame(dat_clust)
## order statistics ##
x <- sort(dat_clust,decreasing = F)
x_backup <- x


#------------------------- gamma distribution MLE -----------------------------#
sum_log_x <- sum(log(x))
sum_x <- sum(x)
gamma_den <- function(para){
  alp <- para[1]
  bet <- para[2]
  n <- length(x)
  if(alp<=0 || bet <=0)
    return(NA)
  else
    return(-(n*alp*log(bet)-n*log(gamma(alp))+(alp-1)*sum_log_x-bet*sum_x))
}

mu <- mean(x)
sd <- var(x)
# para_hat <- stats::optim(par = c(mu^2/sd,mu/sd), fn = gamma_den,
#                          method = "L-BFGS-B",lower = 0)$par
para_hat <- stats::optim(par = c(10,0.01), fn = gamma_den,
                         method = "L-BFGS-B",lower = 10^(-5))$par

print("MLE Gamma distribution parameters are")
print(round(para_hat,2))

#------------------------ transaction CDF histogram ---------------------------#
cdf_x <- AnomDetct::cdf_convert(x, "gamma", 
                                shape = para_hat[1], rate = para_hat[2])

#----------------------------- scan window length -----------------------------#
MLE_window_lth <- function(x,dist_null,..., unit = 1,
                           quant = 0.99){
  if(unit<=0 | is.infinite(unit)){
    stop("`unit` needs to be positive finite number")
  }
  
  y <- x/unit
  bds <- ceiling(range(y))
  breaks <- seq(bds[1]-1,bds[2])*unit
  x_cnt <- table(cut(x,breaks = breaks))
  
  if(dist_null == "gpd"){
    if(requireNamespace("POT", quietly = TRUE)){
      exp_cnt <- diff(POT::pgpd(breaks, ...))*length(x)
    }else{
      stop("Need package POT", call. = F)
    }
  }else{
    my_str <- paste("p",dist_null,"(breaks, ...)",sep = "")
    exp_cnt  <- diff(eval(parse(text = my_str)))*length(x)
  }
  x_diff <- x_cnt - exp_cnt
  res <- quantile(x_diff, quant)
  return(round(res))
}

window_MLE <- MLE_window_lth(x, dist_null = "gamma",
                             shape = para_hat[1], rate = para_hat[2],
                             unit = 1, quant = 0.99)

#------------------------------ scan threshold --------------------------------#
## Sample size ##
N <- length(x)

## use bandwidth n^(-1/2) ##
bandwidth <- N^(-1/2)

## window length ##
window_l <- max(30, window_MLE)

## theta ##
theta_0 <- 1

## alpha-level of test ##  
alpha_0 <- .05

## density function estimation points ##
pts <- seq(0,1,0.05)

##--------------------------- Predefine Functions ----------------------------##
## Rolling Sum, return a sequence of sums ##
rolling_sum <- function(temp,window_l){
  m <- length(temp)
  count <- sapply(1:(m-window_l+1),
                  FUN = function(x)sum(temp[x:(x+window_l-1)]))
  return(count)
}

## Probablity for sum of success via one-dependent sequence approximation ##
prob_fun <- function(k,p,window_l){
  # ## function b and F ##
  
  b_fun <- function(k,m) return(dbinom(k,m,p))
  F_fun <- function(r,s) return(pbinom(r,s,p))
  ## function q1 ##
  q1_function <- function(k, m = window_l){
    F_fun(k,m)^2 - 
      k*b_fun(k+1,m)*F_fun(k-1, m) + 
      m*p*b_fun(k+1,m)*F_fun(k-2,m-1)
  }
  
  ## function q2 ##
  
  A1 <- function(k,m){
    2*b_fun(k+1,m)*F_fun(k,m)*(k*F_fun(k-1,m) - m*p*F_fun(k-2,m-1))
  }
  
  A2 <- function(k,m){
    0.5*b_fun(k+1,m)^2*(k*(k-1)*F_fun(k-2,m) - 
                          2*(k-1)*m*F_fun(k-3,m-1) + 
                          m*(m-1)*p^2*F_fun(k-4,m-2))
  }
  
  A3 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-1:k, m = m) * mapply(F_fun, 0:(k-1), s = m)^2)
  }
  
  A4 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-2:k, m = m) * 
          mapply(b_fun, 2:k+1, m = m) *
          (2:k * mapply(F_fun, 2:k-1, s = m) - 
             m*p*mapply(F_fun, 2:k-2, s = m-1)))
  }
  
  q2_function <- function(k,m=window_l){
    F_fun(k,m)^3 - A1(k,m) + A2(k,m) + A3(k,m) - A4(k,m)
  }
  
  
  ## function qL ##
  qL_function <- function(L,k, m = window_l){
    q1 = q1_function(k,m)
    q2 = q2_function(k,m)
    (2*q1 - q2)/(1 + q1 - q2 + 2*(q1-q2)^2)^L
  }
  
  ## linear interpolation ##
  n <- (N-1)/window_l
  qL_lower <- qL_function(floor(n),k = k,m = window_l)
  qL_upper <- qL_function(ceiling(n),k = k,m = window_l)
  return(qL_upper*(n-floor(n)) + qL_lower*(1-(n-floor(n))))
}

#---------- Estimate critincal cnt on different window and alpha level --------#
crit_theta_fun <- function(window_l,alpha_0,theta = theta_0){
  ## binary search ##
  lower_crit <- 0
  upper_crit <- window_l
  
  while(lower_crit < upper_crit){
    critical_cnt <- floor(mean(c(lower_crit,upper_crit)))
    p <- prob_fun(k = critical_cnt, p = 1-exp(-theta),window_l = window_l)
    ## select the smallest p>=1-alpha_0
    if(p < 1-alpha_0){
      lower_crit <- critical_cnt + 1
    }else{
      upper_crit <- critical_cnt
    }
  }
  
  return(critical_cnt+1) ## P(x >= k)<=0.05
}

## fhat point estimation ##
## max_size to limited the number of x, machine can't process too large N
pts_est <- function(pts,x,bandwidth, max_size = 10000){
  if(length(x)>max_size) x <- sort(sample(x, max_size))
  x_ext <- c(rev(-x),x,rev(2-x)) # reflection
  # The empircal CDF on each observation is 1/n, 2/n,..., 1 #
  Yemp <- (1:length(x_ext))/length(x_ext)
  
  ## def function to estimate f(x) on each given point pt ##
  pt_est_ <- function(pt){
    ## build and solve x*beta = yhat ##
    Wp       <- diag(sqrt(stats::dnorm((x_ext-pt)/bandwidth)))
    dsg_mtx  <- cbind(rep(1,length(x_ext)), # constant term
                      (x_ext-pt), # linear term
                      (x_ext-pt)^2) # quadratic term
    X_mtx    <- Wp %*% dsg_mtx
    
    ## beta_hat[1] is F(x), beta_hat[2] is f(x) ##
    ## This is the f(x) of refected x, real f(x) should *3 ##
    beta_hat <- Matrix::solve(t(X_mtx)%*%X_mtx, t(X_mtx)%*%Wp%*%Yemp)
    
    return(beta_hat)
  }
  
  ## estimate f(x) on each pt and apply linear interp. on other points ##
  fhat_pts <- parallel::mcmapply(pt_est_, pts, SIMPLIFY = TRUE)
  
  return(stats::approxfun(pts,fhat_pts[2,]*3))
}

## main function scanning on observations ##
scan_fnc <- function(x,fhat,critical_cnt,theta,window_l){
  x_with_bd <- c(0,x,1)
  N <- length(x_with_bd)
  ## N is different as n in the paper. ##
  ## N = n+2 with boundary points 0 and 1. ##
  expd_x <- theta/((N-1)*(fhat(x_with_bd[1:(N-1)]))) + x_with_bd[1:(N-1)]
  indicator <- x_with_bd[2:N] <= expd_x
  
  cnt <- rolling_sum(temp = indicator, window_l = window_l)
  clst_loca <- which(cnt >= critical_cnt)
  ## If indicator = 1, that means the NEXT obs. is coming faster ##
  start_index <- c(clst_loca[1], 
                   ifelse(diff(clst_loca) < window_l, NA, clst_loca[-1])) + 1
  end_index <- c(ifelse(diff(clst_loca) < window_l, 
                        NA, clst_loca[-length(clst_loca)] + window_l),
                 clst_loca[length(clst_loca)] + window_l)
  start_index <- na.omit(start_index)
  end_index <- na.omit(end_index)
  
  ## all cluster are shifted to right by 1 unit becasue we added point 0 ##
  ## so here minus 1 to transfer the index back ## 
 
  clst <- matrix(c(start_index-1, end_index-1), ncol = 2, byrow = F)
  colnames(clst) <- c("start", "end")
  
  return(clst)
}

## plot function ##
clst_plot <- function(x, clst, alpha_lvl = 0.05, unit = 1,
                      plt_mgn = 0, lower_bd = 51, upper_percentile = 0.99){
  x <- c(0,x,max(x))
  org_x <- x
  x_range <- c(min(org_x),max(org_x))
  
  ## drop large alpha level
  ## transaction < 50 
  ## cluster width > 100
  clst <- lapply(clst,function(y)y[(y[,3]<=alpha_lvl)&
                                     (x[y[,1]]>lower_bd)&
                                     (y[,2] <= upper_percentile*N),
                                   ,drop = F])
  clst_x <- lapply(clst,function(foo)
    t(apply(foo[,1:2,drop = F],1,function(i)x[i])))
  
  ## cluster colors ##
  p_values <- unlist(lapply(clst,function(x)x[,3,drop = F]))
  num_col  <- (log10(alpha_lvl+10^-16) - log10(p_values+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  ## histogram breaks ##
  bds <- ceiling(range(org_x/unit))
  breaks <- seq(bds[1]-1,bds[2])*unit
  
  df_org_x <- data.frame(x = org_x)
  my_plot <- ggplot2::ggplot(df_org_x, ggplot2::aes(x)) +
    ggplot2::stat_bin(breaks = breaks,
                      closed = "right",
                      fill = "white",col = "black") +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(xlim = x_range)
  
  y_max <- ggplot2::layer_scales(my_plot)$y$range$range[2]
  
  ## cluster data.frame convert ##
  plot_data <- function(ind){
    x_lim <- clst_x[[ind]]
    y_bot <- (ind-1)*y_max/length(clst_x)
    y_up  <- ind*y_max/length(clst_x)
    
    res   <- data.frame(x = numeric(0),
                        y = numeric(0),
                        col = numeric(0),
                        group = character(0))
    
    if(ncol(x_lim)>0){
      for (i in seq.int(nrow(x_lim))){
        temp <- data.frame(x = rep(c(x_lim[i,1] - plt_mgn,
                                     x_lim[i,2] + plt_mgn), each = 2),
                           y = c(y_bot,y_up,y_up,y_bot),
                           group = paste(ind, i, sep = "-"))
        res <- rbind.data.frame(res,temp)
      }
    }
    return(res)
  }
  
  if (requireNamespace("parallel", quietly = TRUE)){
    df <- parallel::mclapply(seq.int(length(clst_x)), plot_data)
  }else{
    df <- lapply(seq.int(length(clst_x)), plot_data)
  }
  
  df_clst <- do.call("rbind.data.frame",df)
  df_clst$col <- rep(num_col,each = 4)
  
  ## cluster plot ##
  p_round <- round(alpha_lvl,3)
  a <- nchar(p_round)-2
  b <- p_round*10^a
  p_labels <- b*10^-(seq(a,16,by = 2))
  p_breaks  <- (log10(alpha_lvl+10^-16) - log10(p_labels+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  p_labels <- formatC(p_labels,format = "e",digits = 0)
  
  my_plot <- my_plot +
    ggplot2::geom_polygon(data = df_clst,
                          ggplot2::aes(x = x, y = y,
                                       group = group, fill = col)) +
    ggplot2::scale_fill_gradient(
      low = grDevices::rgb(red = 0.1,green = 0.1,blue = 1,alpha = 0.6),
      high = grDevices::rgb(red = 1,green = 0.1,blue = 0.1,alpha = 0.6),
      breaks = p_breaks, name = expression( ~ alpha ~level),
      labels = p_labels, guide = "colourbar",limits = c(0,1)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.1,"npc"),
                   plot.title = ggplot2::element_text(hjust = 0.5,size = 40),
                   axis.title = ggplot2::element_text(size = 25),
                   axis.text = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size = 25)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(reverse=T)) +
    ggplot2::xlab("")
  
  return(my_plot)
}


critical_cnt <- crit_theta_fun(window_l = window_l,
                               alpha_0,
                               theta = theta_0)

#------------------------------- histogram CDF --------------------------------#
df <- data.frame(x = cdf_x)
p <- ggplot(df,aes(x)) + 
  stat_bin(breaks = seq(0,1,by = 0.05),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Transaction Amount Percentile") +
  ggtitle("Transaction Percentile Histogram") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

plot(p)

#------------------------------- estimate fhat --------------------------------#
## The first time run this code

fhat_fun <- stats::approxfun(c(0,1), c(1,1))
fhat <- fhat_fun


#-------------- detect clusters on different theta level ----------------------#
seq_theta <- seq(0.5, 1, by = 0.05)*theta_0
clst <- parallel::mcmapply(seq_theta, FUN = scan_fnc,
                           MoreArgs = list(x = cdf_x,
                                           fhat = fhat,
                                           critical_cnt = critical_cnt,
                                           window_l = window_l),
                           SIMPLIFY = F)

q_values <- mapply(FUN = prob_fun, 1-exp(-seq_theta),
                   k = critical_cnt, window_l = window_l)

## q_values may exceed 1 as the approx. error
clst_p_values <- mapply(function(x,y)if(nrow(x)>0){cbind(x,y)}
                        else{cbind(x,numeric(0L))},
                        clst, mapply(max, 1 - q_values, 0),
                        SIMPLIFY = F)

#------------------------------- clusters plot --------------------------------#
my_plot <- clst_plot(x = x, clst = clst_p_values, unit = 1, 
                     alpha_lvl = 0.05, plt_mgn = 0.5)
plot(my_plot)
#----------------------------- cluster ranking  -------------------------------#
anom_score <- function(bds, alpha_lvl = 0.05, 
                       x, amt_bd = 50, upper_percentile = 0.99){
  ## check clusters from top to bottom
  ## is clusters overlap in x axis, two clusters are should group to a same one
  ## clusters bds = (min(pre_bds, cluster), max(pre_bds, cluster))
  ## loop all to the lowest level of cluster
  x <- c(0, x, max(x)) ## add pseudo bounds
  res <- data.frame(lower_bd = numeric(0L),
                    upper_bd = numeric(0L),
                    average_amount = numeric(0L),
                    score = numeric(0L))
  for(clsts in rev(clst_p_values)){
    for(row in seq.int(nrow(clsts))){
      lb <- x[clsts[row,1]] ## transaction amount
      ub <- x[clsts[row,2]] ## transaction amount
      if(lb<=amt_bd || ub > quantile(x, upper_percentile)) next
      ## observations
      x_obs <- sum(x[x>=lb&x<=ub])
      loc_logic <- pmax(lb, res$lower_bd) <= pmin(ub, res$upper_bd)
      ## new cluster, can only happen on the first clsts
      if(sum(loc_logic, na.rm = TRUE) == 0){
        AA <- mean(x[x>=lb&x<=ub])
        res[nrow(res)+1,] <- c(lb, ub, AA, x_obs)
      } 
      ## merge to old clusters, only one loc_logic can be TRUE
      else{
        res[which(loc_logic),4] <- res[which(loc_logic),4] + x_obs
      }
    }
  }
  return(res)
}

clst_score <- anom_score(bds = clst_p_values, alpha_lvl = 0.05, 
                         x = x_backup)

clst_score %>%
  arrange(desc(score))



#--------------------------Decluster-----------------------------

pos1l <- which(x == clst_score[[1]][[1]])
pos1u <- which(x == clst_score[[2]][[1]])
pos2l <- which(x == clst_score[[1]][[2]])
pos2u <- which(x == clst_score[[2]][[2]])


list_cl <- c(seq(pos1l, pos1u), seq(pos2l, pos2u))

list_out <- list()

positions <- c(list_cl)
A = x[-positions]

#------------------------------refill------------------------------------
length(list_cl)
Tao <- positions[1]-1
Tao_i <- A[Tao]

kappa <- positions[1]
kappa_i <- A[kappa]
d <- A[kappa] - A[Tao]
d.2 <- d/2
t <- A[Tao] - d.2
k <- A[kappa] + d.2

c1 <- A[A < Tao_i & A > t]
c2 <- A[A < k & A > kappa_i]

length(c1)
length(c2)

lambda1 <- length(c1) + length(c2)
set.seed(25)
Y1 <- rpois(1, lambda1)
set.seed(22)
y1 <- runif(Y1, min = Tao_i, max = kappa_i)
print(y1)
hist(y1, breaks = 5)



Tao2c <- positions[151]-1
Tao2c_i <- x[Tao2c]

kappa2c <- positions[151]
kappa2c_i <- x[kappa2c]
d2c <- x[kappa2c] - x[Tao2c]
d.2c <- d2c/2

t2c <- x[Tao2c] - d.2c
k2c <- x[kappa2c] + d.2c

c2c1 <- x[x < Tao2c_i & x > t2c]
c2c2 <- x[x < k2c & x > kappa2c_i]
lambda2 <- length(c2c1) + length(c2c2)

set.seed(15)
Y <- rpois(1, lambda2)
set.seed(25)
y <- runif(Y, min = Tao2c_i, max = kappa2c_i)
print(y)
hist(y, breaks = 5)

B <- c(A, y1, y)
hist(B, breaks = seq(0, 300, 1))
#-------------------Threshold selection----------------------
library(selectiveInference)
library(lmom)
if(!require("eva")) {
  install.packages("eva")
  library(eva)
}


setwd("E:/R research/NorthropColeman2014")
# Reads in functions and nidd.thresh data
source("NorthropColeman2014.fns")


quan <- c(quantile(B, probs=seq(0.10, 0.99, 0.01)))

set.seed(3)
quan.dat <- score.fitrange(B,quan)


sorted <- sort(quan.dat$e.p.values)

khat <- forwardStop(sorted, alpha=.05)
print(khat)

pvalue <- sorted[khat]
print(pvalue)

matched <- which(quan.dat$e.p.values == pvalue)
print(matched)
print(quan)
threshold <- quan[matched]
print(threshold)

print(length(B[B>threshold]))


excedences <- (B[B> threshold]-threshold)

mle_fit <- gpdFit(excedences, threshold = 0, method = "mle")

cdf <- cdfgpa(excedences, para = c(0, mle_fit$par.ests[[1]], mle_fit$par.ests[[2]]))
hist(cdf, freq = 
       FALSE, main = expression(hat(G)*""[xi]*""[beta]*(Y*""[i])), ylim = c(0, 1.2), yaxp=c(0, 1, 2))


#------------------method on excedences--------------------------

x <- sort(excedences,decreasing = F)
e_backup <- x

df <- data.frame(x = x)
#----------------------------- scan window length -----------------------------#
y = x/1
ceiling(range(y))

bds = ceiling(range(y))
breaks = seq(bds[1]-1, bds[2])*1
cnt = table(cut(x, breaks = breaks))

N_excedences <- length(excedences)

expect_cnt = diff(pgpd(breaks))*N_excedences
WL = cnt - expect_cnt
window_l <- round(max(WL), 1)

#------------------------------ scan threshold --------------------------------#
## Sample size ##
N <- length(x)

## use bandwidth n^(-1/2) ##
bandwidth <- N^(-1/2)

## window length ##
#window_l <- max(15, window_MLE)

## theta ##
theta_0 <- 1

## alpha-level of test ##  
alpha_0 <- .05

## density function estimation points ##
pts <- seq(0,1,0.05)

## Probablity for sum of success via one-dependent sequence approximation ##
prob_fun <- function(k,p,window_l){
  # ## function b and F ##

  
  b_fun <- function(k,m) return(dbinom(k,m,p))
  F_fun <- function(r,s) return(pbinom(r,s,p))
  ## function q1 ##
  q1_function <- function(k, m = window_l){
    F_fun(k,m)^2 - 
      k*b_fun(k+1,m)*F_fun(k-1, m) + 
      m*p*b_fun(k+1,m)*F_fun(k-2,m-1)
  }
  
  ## function q2 ##
  
  A1 <- function(k,m){
    2*b_fun(k+1,m)*F_fun(k,m)*(k*F_fun(k-1,m) - m*p*F_fun(k-2,m-1))
  }
  
  A2 <- function(k,m){
    0.5*b_fun(k+1,m)^2*(k*(k-1)*F_fun(k-2,m) - 
                          2*(k-1)*m*F_fun(k-3,m-1) + 
                          m*(m-1)*p^2*F_fun(k-4,m-2))
  }
  
  A3 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-1:k, m = m) * mapply(F_fun, 0:(k-1), s = m)^2)
  }
  
  A4 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-2:k, m = m) * 
          mapply(b_fun, 2:k+1, m = m) *
          (2:k * mapply(F_fun, 2:k-1, s = m) - 
             m*p*mapply(F_fun, 2:k-2, s = m-1)))
  }
  
  q2_function <- function(k,m=window_l){
    F_fun(k,m)^3 - A1(k,m) + A2(k,m) + A3(k,m) - A4(k,m)
  }
  
  
  ## function qL ##
  qL_function <- function(L,k, m = window_l){
    q1 = q1_function(k,m)
    q2 = q2_function(k,m)
    (2*q1 - q2)/(1 + q1 - q2 + 2*(q1-q2)^2)^L
  }
  
  ## linear interpolation ##
  n <- (N-1)/window_l
  qL_lower <- qL_function(floor(n),k = k,m = window_l)
  qL_upper <- qL_function(ceiling(n),k = k,m = window_l)
  return(qL_upper*(n-floor(n)) + qL_lower*(1-(n-floor(n))))
}

#---------- Estimate critincal cnt on different window and alpha level --------#
crit_theta_fun <- function(window_l,alpha_0,theta = theta_0){
  ## binary search ##
  lower_crit <- 0
  upper_crit <- window_l
  
  while(lower_crit < upper_crit){
    critical_cnt <- floor(mean(c(lower_crit,upper_crit)))
    p <- prob_fun(k = critical_cnt, p = 1-exp(-theta),window_l = window_l)
    ## select the smallest p>=1-alpha_0
    if(p < 1-alpha_0){
      lower_crit <- critical_cnt + 1
    }else{
      upper_crit <- critical_cnt
    }
  }
  
  return(critical_cnt+1) ## P(x >= k)<=0.05
}

## fhat point estimation ##
## max_size to limited the number of x, machine can't process too large N
pts_est <- function(pts,x,bandwidth, max_size = 10000){
  if(length(x)>max_size) x <- sort(sample(x, max_size))
  x_ext <- c(rev(-x),x,rev(2-x)) # reflection
  # The empircal CDF on each observation is 1/n, 2/n,..., 1 #
  Yemp <- (1:length(x_ext))/length(x_ext)
  
  ## def function to estimate f(x) on each given point pt ##
  pt_est_ <- function(pt){
    ## build and solve x*beta = yhat ##
    Wp       <- diag(sqrt(stats::dnorm((x_ext-pt)/bandwidth)))
    dsg_mtx  <- cbind(rep(1,length(x_ext)), # constant term
                      (x_ext-pt), # linear term
                      (x_ext-pt)^2) # quadratic term
    X_mtx    <- Wp %*% dsg_mtx
    
    ## beta_hat[1] is F(x), beta_hat[2] is f(x) ##
    ## This is the f(x) of refected x, real f(x) should *3 ##
    beta_hat <- Matrix::solve(t(X_mtx)%*%X_mtx, t(X_mtx)%*%Wp%*%Yemp)
    
    return(beta_hat)
  }
  
  ## estimate f(x) on each pt and apply linear interp. on other points ##
  fhat_pts <- parallel::mcmapply(pt_est_, pts, SIMPLIFY = TRUE)
  
  return(stats::approxfun(pts,fhat_pts[2,]*3))
}

## main function scanning on observations ##
scan_fnc <- function(x,fhat,critical_cnt,theta,window_l){
  x_with_bd <- c(0,x,1)
  N <- length(x_with_bd)
  ## N is different as n in the paper. ##
  ## N = n+2 with boundary points 0 and 1. ##
  expd_x <- theta/((N-1)*(fhat(x_with_bd[1:(N-1)]))) + x_with_bd[1:(N-1)]
  indicator <- x_with_bd[2:N] <= expd_x
  
  cnt <- rolling_sum(temp = indicator, window_l = window_l)
  clst_loca <- which(cnt >= critical_cnt)
  ## If indicator = 1, that means the NEXT obs. is coming faster ##
  start_index <- c(clst_loca[1], 
                   ifelse(diff(clst_loca) < window_l, NA, clst_loca[-1])) + 1
  end_index <- c(ifelse(diff(clst_loca) < window_l, 
                        NA, clst_loca[-length(clst_loca)] + window_l),
                 clst_loca[length(clst_loca)] + window_l)
  start_index <- na.omit(start_index)
  end_index <- na.omit(end_index)
  
  ## all cluster are shifted to right by 1 unit becasue we added point 0 ##
  ## so here minus 1 to transfer the index back ## 
 
  clst <- matrix(c(start_index-1, end_index-1), ncol = 2, byrow = F)
  colnames(clst) <- c("start", "end")
  
  return(clst)
}

## plot function ##
clst_plot <- function(x, clst, alpha_lvl = 0.05, unit = 1,
                      plt_mgn = 0, lower_bd = 51, upper_percentile = 0.99){
  x <- c(0,x,max(x))
  org_x <- x
  x_range <- c(min(org_x),max(org_x))
  
  ## drop large alpha level
  ## transaction < 50 
  ## cluster width > 100
  clst <- lapply(clst,function(y)y[(y[,3]<=alpha_lvl)&
                                     (x[y[,1]]>lower_bd)&
                                     (y[,2] <= upper_percentile*N),
                                   ,drop = F])
  clst_x <- lapply(clst,function(foo)
    t(apply(foo[,1:2,drop = F],1,function(i)x[i])))
  
  ## cluster colors ##
  p_values <- unlist(lapply(clst,function(x)x[,3,drop = F]))
  num_col  <- (log10(alpha_lvl+10^-16) - log10(p_values+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  ## histogram breaks ##
  bds <- ceiling(range(org_x/unit))
  breaks <- seq(bds[1]-1,bds[2])*unit
  
  df_org_x <- data.frame(x = org_x)
  my_plot <- ggplot2::ggplot(df_org_x, ggplot2::aes(x)) +
    ggplot2::stat_bin(breaks = breaks,
                      closed = "right",
                      fill = "white",col = "black") +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(xlim = x_range)
  
  y_max <- ggplot2::layer_scales(my_plot)$y$range$range[2]
  
  ## cluster data.frame convert ##
  plot_data <- function(ind){
    x_lim <- clst_x[[ind]]
    y_bot <- (ind-1)*y_max/length(clst_x)
    y_up  <- ind*y_max/length(clst_x)
    
    res   <- data.frame(x = numeric(0),
                        y = numeric(0),
                        col = numeric(0),
                        group = character(0))
    
    if(ncol(x_lim)>0){
      for (i in seq.int(nrow(x_lim))){
        temp <- data.frame(x = rep(c(x_lim[i,1] - plt_mgn,
                                     x_lim[i,2] + plt_mgn), each = 2),
                           y = c(y_bot,y_up,y_up,y_bot),
                           group = paste(ind, i, sep = "-"))
        res <- rbind.data.frame(res,temp)
      }
    }
    return(res)
  }
  
  if (requireNamespace("parallel", quietly = TRUE)){
    df <- parallel::mclapply(seq.int(length(clst_x)), plot_data)
  }else{
    df <- lapply(seq.int(length(clst_x)), plot_data)
  }
  
  df_clst <- do.call("rbind.data.frame",df)
  df_clst$col <- rep(num_col,each = 4)
  
  ## cluster plot ##
  p_round <- round(alpha_lvl,3)
  a <- nchar(p_round)-2
  b <- p_round*10^a
  p_labels <- b*10^-(seq(a,16,by = 2))
  p_breaks  <- (log10(alpha_lvl+10^-16) - log10(p_labels+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  p_labels <- formatC(p_labels,format = "e",digits = 0)
  
  my_plot <- my_plot +
    ggplot2::geom_polygon(data = df_clst,
                          ggplot2::aes(x = x, y = y,
                                       group = group, fill = col)) +
    ggplot2::scale_fill_gradient(
      low = grDevices::rgb(red = 0.1,green = 0.1,blue = 1,alpha = 0.6),
      high = grDevices::rgb(red = 1,green = 0.1,blue = 0.1,alpha = 0.6),
      breaks = p_breaks, name = expression( ~ alpha ~level),
      labels = p_labels, guide = "colourbar",limits = c(0,1)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.1,"npc"),
                   plot.title = ggplot2::element_text(hjust = 0.5,size = 40),
                   axis.title = ggplot2::element_text(size = 25),
                   axis.text = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size = 25)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(reverse=T)) +
    ggplot2::xlab("")
  
  return(my_plot)
}

critical_cnt <- crit_theta_fun(window_l = window_l,
                               alpha_0,
                               theta = theta_0)

#------------------------------- histogram CDF --------------------------------#
df <- data.frame(x = cdf)
p <- ggplot(df,aes(x)) + 
  stat_bin(breaks = seq(0,1,by = 0.05),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Transaction Amount Percentile") +
  ggtitle("Transaction Percentile Histogram") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
plot(p)
#------------------------------- estimate fhat --------------------------------#


fhat_fun <- stats::approxfun(c(0,1), c(1,1))

fhat <- fhat_fun
#-------------- detect clusters on different theta level ----------------------#
seq_theta <- seq(0.5, 1, by = 0.05)*theta_0
clst <- parallel::mcmapply(seq_theta, FUN = scan_fnc,
                           MoreArgs = list(x = cdf,
                                           fhat = fhat,
                                           critical_cnt = critical_cnt,
                                           window_l = window_l),
                           SIMPLIFY = F)



q_values <- mapply(FUN = prob_fun, 1-exp(-seq_theta),
                   k = critical_cnt, window_l = window_l)

## q_values may exceed 1 as the approx. error
clst_p_values <- mapply(function(x,y)if(nrow(x)>0){cbind(x,y)}
                        else{cbind(x,numeric(0L))},
                        clst, mapply(max, 1 - q_values, 0),
                        SIMPLIFY = F)

#------------------------------- clusters plot --------------------------------#
my_plot <- clst_plot(x = x, clst = clst_p_values, unit = 1, 
                     alpha_lvl = 0.05, plt_mgn = 0.5)
plot(my_plot)

#----------------------------- cluster ranking  -------------------------------#
anom_score <- function(bds, alpha_lvl = 0.05, 
                       x, amt_bd = 50, upper_percentile = 0.99){
  ## check clusters from top to bottom
  ## is clusters overlap in x axis, two clusters are should group to a same one
  ## clusters bds = (min(pre_bds, cluster), max(pre_bds, cluster))
  ## loop all to the lowest level of cluster
  x <- c(0, x, max(x)) ## add pseudo bounds
  res <- data.frame(lower_bd = numeric(0L),
                    upper_bd = numeric(0L),
                    average_amount = numeric(0L),
                    score = numeric(0L))
  for(clsts in rev(clst_p_values)){
    for(row in seq.int(nrow(clsts))){
      lb <- x[clsts[row,1]] ## transaction amount
      ub <- x[clsts[row,2]] ## transaction amount
      if(lb<=amt_bd || ub > quantile(x, upper_percentile)) next
      ## observations
      x_obs <- sum(x[x>=lb&x<=ub])
      loc_logic <- pmax(lb, res$lower_bd) <= pmin(ub, res$upper_bd)
      ## new cluster, can only happen on the first clsts
      if(sum(loc_logic, na.rm = TRUE) == 0){
        AA <- mean(x[x>=lb&x<=ub])
        res[nrow(res)+1,] <- c(lb, ub, AA, x_obs)
      } 
      ## merge to old clusters, only one loc_logic can be TRUE
      else{
        res[which(loc_logic),4] <- res[which(loc_logic),4] + x_obs
      }
    }
  }
  return(res)
}

clst_score <- anom_score(bds = clst_p_values, alpha_lvl = 0.05, 
                         x = e_backup)

clst_score %>%
  arrange(desc(score))
#--------------------------Decluster-----------------------------

pos3l <- which(B == (clst_score[[1]][[1]] + threshold)) 
pos3u <- which(B == (clst_score[[2]][[1]] + threshold))

list_cl_small <- c(seq(pos3l, pos3u))
list_out <- list()

positions2 <- c(list_cl_small)

G <- B[-positions2]


#----------------------refill_data---------------------
length(list_cl_small)
Tao_2 <- positions2[1]-1
Tao_2i <- G[Tao_2]

kappa_2 <- positions2[1]
kappa_2i <- G[kappa_2]
d3 <- G[kappa_2] - G[Tao_2]
d3.3 <- d3/2
t_3 <- G[Tao_2] - d3.3
k_3 <- G[kappa_2] + d3.3
Avg_Amt_1Thresh <- (Tao_2i + kappa_2i)/2
c1_3 <- G[G < Tao_2i & G > t_3]
c2_3 <- G[G < k_3 & G > kappa_2i]

length(c1_3)
length(c2_3)

lambda_3 <- length(c1_3) + length(c2_3)

set.seed(1)
Y_2 <- rpois(1, lambda_3)
set.seed(3)
y_2 <- runif(Y_2, min = Tao_2i, max = kappa_2i)
print(y_2)
hist(y_2, breaks = 5)

L <- c(G, y_2)
hist(L)
hist(L, breaks = seq(from = -.50, to =300, by = 1))

L_df <- as.data.frame(L)

#-----------------------excedences refilled histogram----------------------------------
excedence_refilled <- L[L > threshold] - threshold
excedence_refilled_df <- as.data.frame(excedence_refilled)
Excedence_Refilled_df <- ggplot(excedence_refilled_df, aes(excedence_refilled)) +
  stat_bin(breaks = seq(-.50,150,by = 1),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Trasaction Amount $") +
  ggtitle("Excedences Declustered and Refilled") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylim(0, 20)
plot(Excedence_Refilled_df)

#-------------------Threshold selection----------------------



setwd("E:/R research/NorthropColeman2014")
# Reads in functions and nidd.thresh data
source("NorthropColeman2014.fns")


quan <- c(quantile(erf, probs=seq(0.10, 0.99, 0.01)))

set.seed(14)
quan.dat <- score.fitrange(erf,quan)


sorted <- sort(quan.dat$e.p.values)

khat <- forwardStop(sorted, alpha=.05)
print(khat)

pvalue <- sorted[khat]
print(pvalue)

matched <- which(quan.dat$e.p.values == pvalue)
print(matched)
print(quan)
threshold <- quan[matched]

print(threshold)
u2 <- erf + 135.65
print(length(u2[u2>threshold + 135.65]))


excedences <- (u2[u2>threshold + 135.65])

print(length(excedences))


hist(excedences, breaks = 10)
hist(excedences, breaks = 100)


mle_fit <- gpdFit(excedences, threshold = 0, method = "mle")

cdf <- cdfgpa(excedences, para = c(0, mle_fit$par.ests[[1]], mle_fit$par.ests[[2]]))
hist(cdf, freq = 
       FALSE, main = expression(hat(G)*""[xi]*""[beta]*(Y*""[i])), ylim = c(0, 1.2), yaxp=c(0, 1, 2))


#------------------method on excedences--------------------------

x <- sort(excedences,decreasing = F)
e_backup <- x

df <- data.frame(x = x)
p <- ggplot(df,aes(x)) + 
  stat_bin(breaks = seq(0,max(x),by = 1),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Transaction Amount $") +
  ggtitle("Transaction Histogram") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
plot(p)



#------------------------ transaction CDF histogram ---------------------------#

p <- ggplot(data.frame(x = cdf),aes(x)) + 
  stat_bin(breaks = seq(0,1,by = 0.05),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Quantile") +
  ggtitle("CDF Histogram") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

plot(p)

#----------------------------- scan window length -----------------------------#
y = x/1
ceiling(range(y))

bds = ceiling(range(y))
breaks = seq(bds[1]-1, bds[2])*1
cnt = table(cut(x, breaks = breaks))

N_excedences <- length(excedences)

expect_cnt = diff(pgpd(breaks))*N_excedences
WL = cnt - expect_cnt
window_l <- round(max(WL), 1)

#------------------------------Bernoulli construction--------------------------#

i <- seq(1, length(excedences))

K1 <- c(cdf[i])
J1 <-c(0, cdf[i-1]) 

num <- K1 - J1
denom <- (1-J1)/(length(excedences)-i+1)

V_i <- num/denom
#V_i <- (K1 - J1)/((1-J1)/length(excedences)-i+1)
hist(V_i, breaks = 100)
Y_i <- ifelse(V_i <= 1, 1, 0)


print(Y_i)
#------------------------------ scan threshold --------------------------------#
## Sample size ##
N <- length(x)

## use bandwidth n^(-1/2) ##
bandwidth <- N^(-1/2)

## window length ##
#window_l <- max(15, window_MLE)

## theta ##
theta_0 <- 1

## alpha-level of test ##  
alpha_0 <- .05

## density function estimation points ##
pts <- seq(0,1,0.05)

## Probablity for sum of success via one-dependent sequence approximation ##
prob_fun <- function(k,p,window_l){
  # ## function b and F ##

  
  b_fun <- function(k,m) return(dbinom(k,m,p))
  F_fun <- function(r,s) return(pbinom(r,s,p))
  ## function q1 ##
  q1_function <- function(k, m = window_l){
    F_fun(k,m)^2 - 
      k*b_fun(k+1,m)*F_fun(k-1, m) + 
      m*p*b_fun(k+1,m)*F_fun(k-2,m-1)
  }
  
  ## function q2 ##
  
  A1 <- function(k,m){
    2*b_fun(k+1,m)*F_fun(k,m)*(k*F_fun(k-1,m) - m*p*F_fun(k-2,m-1))
  }
  
  A2 <- function(k,m){
    0.5*b_fun(k+1,m)^2*(k*(k-1)*F_fun(k-2,m) - 
                          2*(k-1)*m*F_fun(k-3,m-1) + 
                          m*(m-1)*p^2*F_fun(k-4,m-2))
  }
  
  A3 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-1:k, m = m) * mapply(F_fun, 0:(k-1), s = m)^2)
  }
  
  A4 <- function(k,m){
    sum(mapply(b_fun, 2*(k+1)-2:k, m = m) * 
          mapply(b_fun, 2:k+1, m = m) *
          (2:k * mapply(F_fun, 2:k-1, s = m) - 
             m*p*mapply(F_fun, 2:k-2, s = m-1)))
  }
  
  q2_function <- function(k,m=window_l){
    F_fun(k,m)^3 - A1(k,m) + A2(k,m) + A3(k,m) - A4(k,m)
  }
  
  
  ## function qL ##
  qL_function <- function(L,k, m = window_l){
    q1 = q1_function(k,m)
    q2 = q2_function(k,m)
    (2*q1 - q2)/(1 + q1 - q2 + 2*(q1-q2)^2)^L
  }
  
  ## linear interpolation ##
  n <- (N-1)/window_l
  qL_lower <- qL_function(floor(n),k = k,m = window_l)
  qL_upper <- qL_function(ceiling(n),k = k,m = window_l)
  return(qL_upper*(n-floor(n)) + qL_lower*(1-(n-floor(n))))
}

#---------- Estimate critincal cnt on different window and alpha level --------#
crit_theta_fun <- function(window_l,alpha_0,theta = theta_0){
  ## binary search ##
  lower_crit <- 0
  upper_crit <- window_l
  
  while(lower_crit < upper_crit){
    critical_cnt <- floor(mean(c(lower_crit,upper_crit)))
    p <- prob_fun(k = critical_cnt, p = 1-exp(-theta),window_l = window_l)
    ## select the smallest p>=1-alpha_0
    if(p < 1-alpha_0){
      lower_crit <- critical_cnt + 1
    }else{
      upper_crit <- critical_cnt
    }
  }
  
  return(critical_cnt+1) ## P(x >= k)<=0.05
}

## fhat point estimation ##
## max_size to limited the number of x, machine can't process too large N
pts_est <- function(pts,x,bandwidth, max_size = 10000){
  if(length(x)>max_size) x <- sort(sample(x, max_size))
  x_ext <- c(rev(-x),x,rev(2-x)) # reflection
  # The empircal CDF on each observation is 1/n, 2/n,..., 1 #
  Yemp <- (1:length(x_ext))/length(x_ext)
  
  ## def function to estimate f(x) on each given point pt ##
  pt_est_ <- function(pt){
    ## build and solve x*beta = yhat ##
    Wp       <- diag(sqrt(stats::dnorm((x_ext-pt)/bandwidth)))
    dsg_mtx  <- cbind(rep(1,length(x_ext)), # constant term
                      (x_ext-pt), # linear term
                      (x_ext-pt)^2) # quadratic term
    X_mtx    <- Wp %*% dsg_mtx
    
    ## beta_hat[1] is F(x), beta_hat[2] is f(x) ##
    ## This is the f(x) of refected x, real f(x) should *3 ##
    beta_hat <- Matrix::solve(t(X_mtx)%*%X_mtx, t(X_mtx)%*%Wp%*%Yemp)
    
    return(beta_hat)
  }
  
  ## estimate f(x) on each pt and apply linear interp. on other points ##
  fhat_pts <- parallel::mcmapply(pt_est_, pts, SIMPLIFY = TRUE)
  
  return(stats::approxfun(pts,fhat_pts[2,]*3))
}

## main function scanning on observations ##
scan_fnc <- function(x,fhat,critical_cnt,theta,window_l){
  x_with_bd <- c(0,x,1)
  N <- length(x_with_bd)
  ## N is different as n in the paper. ##
  ## N = n+2 with boundary points 0 and 1. ##
  expd_x <- theta/((N-1)*(fhat(x_with_bd[1:(N-1)]))) + x_with_bd[1:(N-1)]
  indicator <- x_with_bd[2:N] <= expd_x
  
  cnt <- rolling_sum(temp = indicator, window_l = window_l)
  clst_loca <- which(cnt >= critical_cnt)
  ## If indicator = 1, that means the NEXT obs. is coming faster ##
  start_index <- c(clst_loca[1], 
                   ifelse(diff(clst_loca) < window_l, NA, clst_loca[-1])) + 1
  end_index <- c(ifelse(diff(clst_loca) < window_l, 
                        NA, clst_loca[-length(clst_loca)] + window_l),
                 clst_loca[length(clst_loca)] + window_l)
  start_index <- na.omit(start_index)
  end_index <- na.omit(end_index)
  
  ## all cluster are shifted to right by 1 unit becasue we added point 0 ##
  ## so here minus 1 to transfer the index back ## 
  clst <- matrix(c(start_index-1, end_index-1), ncol = 2, byrow = F)
  colnames(clst) <- c("start", "end")
  
  return(clst)
}

## plot function ##
clst_plot <- function(x, clst, alpha_lvl = 0.05, unit = 1,
                      plt_mgn = 0, lower_bd = 51, upper_percentile = 0.99){
  x <- c(0,x,max(x))
  org_x <- x
  x_range <- c(min(org_x),max(org_x))
  
  ## drop large alpha level
  ## transaction < 50 
  ## cluster width > 100
  clst <- lapply(clst,function(y)y[(y[,3]<=alpha_lvl)&
                                     (x[y[,1]]>lower_bd)&
                                     (y[,2] <= upper_percentile*N),
                                   ,drop = F])
  clst_x <- lapply(clst,function(foo)
    t(apply(foo[,1:2,drop = F],1,function(i)x[i])))
  
  ## cluster colors ##
  p_values <- unlist(lapply(clst,function(x)x[,3,drop = F]))
  num_col  <- (log10(alpha_lvl+10^-16) - log10(p_values+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  ## histogram breaks ##
  bds <- ceiling(range(org_x/unit))
  breaks <- seq(bds[1]-1,bds[2])*unit
  
  df_org_x <- data.frame(x = org_x)
  my_plot <- ggplot2::ggplot(df_org_x, ggplot2::aes(x)) +
    ggplot2::stat_bin(breaks = breaks,
                      closed = "right",
                      fill = "white",col = "black") +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(xlim = x_range)
  
  y_max <- ggplot2::layer_scales(my_plot)$y$range$range[2]
  
  ## cluster data.frame convert ##
  plot_data <- function(ind){
    x_lim <- clst_x[[ind]]
    y_bot <- (ind-1)*y_max/length(clst_x)
    y_up  <- ind*y_max/length(clst_x)
    
    res   <- data.frame(x = numeric(0),
                        y = numeric(0),
                        col = numeric(0),
                        group = character(0))
    
    if(ncol(x_lim)>0){
      for (i in seq.int(nrow(x_lim))){
        temp <- data.frame(x = rep(c(x_lim[i,1] - plt_mgn,
                                     x_lim[i,2] + plt_mgn), each = 2),
                           y = c(y_bot,y_up,y_up,y_bot),
                           group = paste(ind, i, sep = "-"))
        res <- rbind.data.frame(res,temp)
      }
    }
    return(res)
  }
  
  if (requireNamespace("parallel", quietly = TRUE)){
    df <- parallel::mclapply(seq.int(length(clst_x)), plot_data)
  }else{
    df <- lapply(seq.int(length(clst_x)), plot_data)
  }
  
  df_clst <- do.call("rbind.data.frame",df)
  df_clst$col <- rep(num_col,each = 4)
  
  ## cluster plot ##
  p_round <- round(alpha_lvl,3)
  a <- nchar(p_round)-2
  b <- p_round*10^a
  p_labels <- b*10^-(seq(a,16,by = 2))
  p_breaks  <- (log10(alpha_lvl+10^-16) - log10(p_labels+10^-16))/
    (log10(alpha_lvl+10^-16)-min(log10(10^-16)))
  
  p_labels <- formatC(p_labels,format = "e",digits = 0)
  
  my_plot <- my_plot +
    ggplot2::geom_polygon(data = df_clst,
                          ggplot2::aes(x = x, y = y,
                                       group = group, fill = col)) +
    ggplot2::scale_fill_gradient(
      low = grDevices::rgb(red = 0.1,green = 0.1,blue = 1,alpha = 0.6),
      high = grDevices::rgb(red = 1,green = 0.1,blue = 0.1,alpha = 0.6),
      breaks = p_breaks, name = expression( ~ alpha ~level),
      labels = p_labels, guide = "colourbar",limits = c(0,1)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.1,"npc"),
                   plot.title = ggplot2::element_text(hjust = 0.5,size = 40),
                   axis.title = ggplot2::element_text(size = 25),
                   axis.text = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size = 25)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(reverse=T)) +
    ggplot2::xlab("")
  
  return(my_plot)
}

critical_cnt <- crit_theta_fun(window_l = window_l,
                               alpha_0,
                               theta = theta_0)

#------------------------------- histogram CDF --------------------------------#
df <- data.frame(x = cdf)
p <- ggplot(df,aes(x)) + 
  stat_bin(breaks = seq(0,1,by = 0.05),
           closed = "right",
           fill = "white",col = "black") +
  theme_bw() +
  xlab("Transaction Amount Percentile") +
  ggtitle("Transaction Percentile Histogram") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

#------------------------------- estimate fhat --------------------------------#


fhat_fun <- stats::approxfun(c(0,1), c(1,1))

fhat <- fhat_fun
#-------------- detect clusters on different theta level ----------------------#
seq_theta <- seq(0.5, 1, by = 0.05)*theta_0
clst <- parallel::mcmapply(seq_theta, FUN = scan_fnc,
                           MoreArgs = list(x = cdf,
                                           fhat = fhat,
                                           critical_cnt = critical_cnt,
                                           window_l = window_l),
                           SIMPLIFY = F)



q_values <- mapply(FUN = prob_fun, 1-exp(-seq_theta),
                   k = critical_cnt, window_l = window_l)

## q_values may exceed 1 as the approx. error
clst_p_values <- mapply(function(x,y)if(nrow(x)>0){cbind(x,y)}
                        else{cbind(x,numeric(0L))},
                        clst, mapply(max, 1 - q_values, 0),
                        SIMPLIFY = F)

#------------------------------- clusters plot --------------------------------#
my_plot <- clst_plot(x = x, clst = clst_p_values, unit = 1, 
                     alpha_lvl = 0.05, plt_mgn = 0.5)
plot(my_plot)





