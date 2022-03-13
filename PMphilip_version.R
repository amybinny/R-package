#' Liptak test
#'
#' Finding the p-value of the modified t-statistic of Kim et al. for partially matched samples problem where the data follows a Gaussian distribution.
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values
#' @param alternative the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return p-value for the test
#' @examples
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:6]<-NA
#' y[14:16]<-NA
#' my.Liptak(x,y,alternative="two.sided", a=0.05)
#' @export
my.Liptak<-function(x,y,alternative="two.sided", a=0.01){
  if (length(x)!=length(y)){
    stop('Number of values in x is different from number of values in y')
  }
  else if (length(is.na(x))==0 & length(is.na(y))==0){
    stop('There are no missing values. Then the test has no purpose here.')
  }
  #taking tumor paired samples
  x_both_samples<-x[!is.na(x) & !is.na(y)]
  #taking normal paired samples
  y_both_samples<-y[!is.na(x) & !is.na(y)]
  #taking  tumor  unpaired samples 
  x_samples_part<-x[!is.na(x) & is.na(y)]
  #taking  normal unpaired samples
  y_samples_part<-y[is.na(x) & !is.na(y)]
  #number of paired sample
  n1<-length(x_both_samples)
  #number of tumor unpaired samples 
  n2<-length(x_samples_part)
  #number of normal unpaired samples
  n3<-length(y_samples_part)
  #Considering different scenarios
  #when we do not have enough paired and unpaired samples
  if (n1<2 & (n2<2 | n3<2)){
    stop('Do not have enough data to perform the test')
  }
  var.eq=TRUE
  if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value<a){
    var.eq = FALSE
  }
  else if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value>=a){
    var.eq = TRUE
  }
  #when we do not have enough paired samples  
  if (n1<2 & (n2>=2 & n3>=2)){
    warning('do not have enough paired samples, so applying unpaired t.test()')
    p_c<-t.test(x_samples_part,y_samples_part,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  else if (n1>=2 & (n2<2 | n3<2)){
    warning('do not have enough unpaired samples, so applying paired t.test()')
    p_c<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  w1<-sqrt(2*n1)
  w2<-sqrt(n2+n3)
  if (alternative=="greater" |alternative=="two.sided"){
    p1<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative="greater", var.equal = var.eq)$p.value
    p2<-t.test(x_samples_part,y_samples_part,paired=FALSE,alternative="greater", var.equal = var.eq)$p.value
    z1<-qnorm(1-p1)
    z2<-qnorm(1-p2)
    p_c<-1-pnorm((w1*z1+w2*z2)/sqrt(w1^2+w2^2))
    if (alternative=="two.sided"){
      if (p_c<1/2){
        p_comb<-2*p_c
      }
      else {
        p_comb<-2*(1-p_c)
      }
    }
    else if (alternative=="greater"){
      p_comb<-p_c
    }
  }
  else if (alternative=="less"){
    p1<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative="less", var.equal = var.eq)$p.value
    p2<-t.test(x_samples_part,y_samples_part,paired=FALSE,alternative="less", var.equal = var.eq)$p.value
    z1<-qnorm(1-p1)
    z2<-qnorm(1-p2)
    p_comb<-1-pnorm((w1*z1+w2*z2)/sqrt(w1^2+w2^2), lower.tail=TRUE)
  }
  return(pvalue=p_comb)
}

#' Kim et al. t-statistic test
#'
#' Finding the p-value of the modified t-statistic of Kim et al. for partially matched samples problem where the data follows a Gaussian distribution.
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values
#' @param alternative the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return p-value for the test
#' @examples 
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:7]<-NA
#' y[14:20]<-NA
#' my.Kim(x,y,alternative="greater", a=0.01)
#' @export
#Kim t-statistic
my.Kim<-function(x,y,alternative="two.sided", a=0.01){
  if (length(x)!=length(y)){
    stop('Number of values in x is different from number of values in y')
  }
  else if (length(is.na(x))==0 & length(is.na(y))==0){
    stop('There are no missing values. Then the test has no purpose here.')
  }
  #taking tumor paired samples
  x_both_samples<-x[!is.na(x) & !is.na(y)]
  #taking normal paired samples
  y_both_samples<-y[!is.na(x) & !is.na(y)]
  #taking  tumor  unpaired samples 
  x_samples_part<-x[!is.na(x) & is.na(y)]
  #taking  normal unpaired samples
  y_samples_part<-y[is.na(x) & !is.na(y)]
  #number of paired sample
  n1<-length(x_both_samples)
  #number of tumor unpaired samples 
  n2<-length(x_samples_part)
  #number of normal unpaired samples
  n3<-length(y_samples_part)
  #Considering different scenarios
  #when we do not have enough paired and unpaired samples
  if (n1<2 & (n2<2 | n3<2)){
    stop('Do not have enough data to perform the test')
  }
  var.eq=TRUE
  if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value<a){
    var.eq = FALSE
  }
  else if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value>=a){
    var.eq = TRUE
  }
  #when we do not have enough paired samples  
  if (n1<2 & (n2>=2 & n3>=2)){
    warning('do not have enough paired samples, so applying unpaired t.test()')
    p_c<-t.test(x_samples_part,y_samples_part,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  else if (n1>=2 & (n2<2 | n3<2)){
    warning('do not have enough unpaired samples, so applying paired t.test()')
    p_c<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  #mean of the difference of paired samples
  D<-mean(x_both_samples-y_both_samples)
  #mean of tumor paired samples
  T<-mean(x_samples_part)
  #mean of normal paired samples
  N<-mean(y_samples_part)
  #standard deviation of the difference of paired samples
  SD_D<-sd(x_both_samples-y_both_samples)
  #standard deviation of tumor paired samples
  SD_T<-sd(x_samples_part)
  #standard deviation of normal paired samples
  SD_N<-sd(y_samples_part)
  #harmonic mean of n2 and n3 
  nH<-2/(1/n2+1/n3)
  #modified t-statistic of Kim et al
  t<-(n1*D+nH*(T-N))/sqrt(n1*SD_D^2+nH^2*(SD_N^2/n3+SD_T^2/n2))
  if (alternative=="two.sided"){
    pvalue_t<-pnorm(t,lower.tail=FALSE)
    if (pvalue_t<1/2){
      pvalue_t3<-2*pvalue_t
    }
    else {
      pvalue_t3<-2*(1-pvalue_t)
    }
  }
  else if (alternative=="greater"){
    pvalue_t3<-pnorm(t,lower.tail=FALSE)
  }
  else if (alternative=="less"){
    pvalue_t3<-pnorm(t,lower.tail=TRUE)
  }
  return(pvalue_t3) 
}

#' Looney and Jones Z-statistic test
#'
#' The corrected Z-test of Looney and Jones is based on a modified variance estimation of the standard Z-test by accounting for the correlation among the matched pairs where the data follows a Gaussian distribution.
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values
#' @param alternative the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return p-value for the test
#' @examples 
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:13]<-NA
#' y[14:26]<-NA
#' my.Looney.Jones(x,y)
#' @export
#Looney and Jones Zcorr: accounts for correlation among paired samples
my.Looney.Jones<-function(x,y,alternative="two.sided", a=0.01){
  if (length(x)!=length(y)){
    stop('Number of values in x is different from number of values in y')
  }
  else if (length(is.na(x))==0 & length(is.na(y))==0){
    stop('There are no missing values. Then the test has no purpose here.')
  }
  #taking tumor paired samples
  x_both_samples<-x[!is.na(x) & !is.na(y)]
  #taking normal paired samples
  y_both_samples<-y[!is.na(x) & !is.na(y)]
  #taking  tumor  unpaired samples 
  x_samples_part<-x[!is.na(x) & is.na(y)]
  #taking  normal unpaired samples
  y_samples_part<-y[is.na(x) & !is.na(y)]
  #number of paired sample
  n1<-length(x_both_samples)
  #number of tumor unpaired samples 
  n2<-length(x_samples_part)
  #number of normal unpaired samples
  n3<-length(y_samples_part)
  #Considering different scenarios
  #when we do not have enough paired and unpaired samples
  if (n1<2 & (n2<2 | n3<2)){
    stop('Do not have enough data to perform the test')
  }
  var.eq=TRUE
  if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value<a){
    var.eq = FALSE
  }
  else if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value>=a){
    var.eq = TRUE
  }
  #when we do not have enough paired samples  
  if (n1<2 & (n2>=2 & n3>=2)){
    warning('do not have enough paired samples, so applying unpaired t.test()')
    p_c<-t.test(x_samples_part,y_samples_part,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  else if (n1>=2 & (n2<2 | n3<2)){
    warning('do not have enough unpaired samples, so applying paired t.test()')
    p_c<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  #mean of tumor paired samples
  T1<-mean(x_both_samples)
  #mean of normal paired samples
  N1<-mean(y_both_samples)
  #mean of tumor paired samples and tumor unpaired samples
  T_<-mean(c(x_both_samples,x_samples_part))
  #mean of normal paired samples and normal unpaired samples
  N_<-mean(c(y_both_samples,y_samples_part))
  #standard deviation of tumor paired samples and tumor unpaired samples
  SD_T_<-sd(c(x_both_samples,x_samples_part))
  #standard deviation of normal paired samples normal unpaired samples
  SD_N_<-sd(c(y_both_samples,y_samples_part))
  #covariance of tumor and normal paired samples
  SD_TN1<-sum((x_both_samples-T1)*(y_both_samples-N1))/(n1-1)
  #SD_TN1<-cov(x_both_samples,y_both_samples)
  #modified variance estimation of the standard Z-test
  Z_corr<-(T_-N_)/sqrt(SD_T_^2/(n1+n2)+SD_N_^2/(n1+n3)-2*n1*SD_TN1/((n1+n2)*(n1+n3)))
  if (alternative=="two.sided"){
    pvalue_Z_corr<-pnorm(Z_corr,lower.tail=FALSE)
    if (pvalue_Z_corr<1/2){
      pvalue_Zcorr<-2*pvalue_Z_corr
    }
    else {
      pvalue_Zcorr<-2*(1-pvalue_Z_corr)
    }
  }
  else if (alternative=="greater"){
    pvalue_Zcorr<-pnorm(Z_corr,lower.tail=FALSE)
  }
  else if (alternative=="less"){
    pvalue_Zcorr<-pnorm(Z_corr,lower.tail=TRUE)
  }
  return(pvalue_Zcorr) 
}
#' Lin and Stivers MLE test
#'
#' Finding the p-value of the MLE based test statistic under heteroscedasticity by Lin and Stivers for partially matched samples problem where the data follows a Gaussian distribution.
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values
#' @param alternative the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param var.equal a logical which says whether we need to treat the variances as equal. If true then pooled variance is calculated otherwise alternate approximation is done
#' @return p-value for the test
#' @examples 
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:13]<-NA
#' y[14:26]<-NA
#' my.Lin.Stivers(x,y)
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:7]<-NA; y[14:20]<-NA
#' my.Lin.Stivers(x,y)
#' @export
#Lin and Stivers MLE test for heteroscedasticity
my.Lin.Stivers<-function(x,y,alternative="two.sided", a=0.01){
  if (length(x)!=length(y)){
    stop('Number of values in x is different from number of values in y')
  }
  else if (length(is.na(x))==0 & length(is.na(y))==0){
    stop('There are no missing values. Then the test has no purpose here.')
  }
  #taking tumor paired samples
  x_both_samples<-x[!is.na(x) & !is.na(y)]
  #taking normal paired samples
  y_both_samples<-y[!is.na(x) & !is.na(y)]
  #taking  tumor  unpaired samples 
  x_samples_part<-x[!is.na(x) & is.na(y)]
  #taking  normal unpaired samples
  y_samples_part<-y[is.na(x) & !is.na(y)]
  #number of paired sample
  n1<-length(x_both_samples)
  #number of tumor unpaired samples 
  n2<-length(x_samples_part)
  #number of normal unpaired samples
  n3<-length(y_samples_part)
  #Considering different scenarios
  #when we do not have enough paired and unpaired samples
  if (n1<2 & (n2<2 | n3<2)){
    stop('Do not have enough data to perform the test')
  }
  var.eq=TRUE
  if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value<a){
    var.eq = FALSE
  }
  else if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value>=a){
    var.eq = TRUE
  }
  #when we do not have enough paired samples  
  if (n1<2 & (n2>=2 & n3>=2)){
    warning('do not have enough paired samples, so applying unpaired t.test()')
    p_c<-t.test(x_samples_part,y_samples_part,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  else if (n1>=2 & (n2<2 | n3<2)){
    warning('do not have enough unpaired samples, so applying paired t.test()')
    p_c<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  #mean of tumor unpaired samples
  T<-mean(x_samples_part)
  #mean of normal unpaired samples
  N<-mean(y_samples_part)
  #mean of tumor paired samples
  T1<-mean(x_both_samples)
  #mean of normal paired samples
  N1<-mean(y_both_samples)
  #standard deviation of tumor paired samples
  SD_T1<-sd(x_both_samples)
  #standard deviation of normal paired samples
  SD_N1<-sd(y_both_samples)
  #covariance of tumor and normal paired samples
  SD_TN1<-sum((x_both_samples-T1)*(y_both_samples-N1))/(n1-1)
  #SD_TN1<-cov(x_both_samples,y_both_samples) 
  r<-SD_TN1/(SD_T1*SD_N1)
  f<-(n1*(n1+n3+n2*SD_TN1/SD_T1^2))/((n1+n2)*(n1+n3)-n2*n3*r^2)
  g<-(n1*(n1+n2+n3*SD_TN1/SD_N1^2))/((n1+n2)*(n1+n3)-n2*n3*r^2)
  V1<-(f^2/n1+((1-f)^2)/n2)*SD_T1^2+(g^2/n1+((1-g)^2)/n3)*SD_N1^2-2*f*g*SD_TN1/n1
  #MLE based test statistic under heteroscedasticity by Lin and Stivers
  Z_ls<-(f*(T1-T)-g*(N1-N)+T-N)/sqrt(V1)
  if (alternative=="two.sided"){
    pvalue_Z_ls<-pt(Z_ls,df=n1, lower.tail=FALSE)
    if (pvalue_Z_ls<1/2){
      pvalue_Zls<-2*pvalue_Z_ls
    }
    else {
      pvalue_Zls<-2*(1-pvalue_Z_ls)
    }
  }
  else if (alternative=="greater"){
    pvalue_Zls<-pt(Z_ls,df=n1, lower.tail=FALSE)
  }
  else if (alternative=="less"){
    pvalue_Zls<-pt(Z_ls,df=n1, lower.tail=TRUE)
  }
  if (var.eq==TRUE){
    warning('Lin and Stivers test works best for the heteroscedastic case, that is when variances are not equal. The data provided here has equal varaince')
    return(pvalue_Zls)
  }
  else {
    return(pvalue_Zls)
  }
}

#' Ekbohm's MLE test
#'
#' Finding the p-value of the MLE based test statistic under homoscedasticity by Ekbohm for partially matched samples problem where the data follows a Gaussian distribution.
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values
#' @param alternative the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return p-value for the test
#' @examples 
#' x<-subset(ChickWeight, Time==8)$weight
#' y<-subset(ChickWeight, Time==10)$weight
#' x[2:6]<-NA
#' y[14:16]<-NA
#' my.Ekbohm(x,y, a=0.05)
#' if(!require(mvtnorm)){
#' install.packages("mvtnorm")}
#' library(mvtnorm)
#' set.seed(123); 
#' s<-matrix(c(1,0.2,0.2,4), ncol=2)
#' d<-rmvnorm(100,c(0,0),s)
#' x<-d[,1]; y<-d[,2]
#' x[2:13]<-NA
#' y[14:26]<-NA
#' my.Ekbohm(x,y, a=0.5)
#' @export
#Ekbohm's MLE test when variances of normal and tumor samples are equal ()
my.Ekbohm<-function(x,y,alternative="two.sided", a=0.01){
  if (length(x)!=length(y)){
    stop('Number of values in x is different from number of values in y')
  }
  else if (length(is.na(x))==0 & length(is.na(y))==0){
    stop('There are no missing values. Then the test has no purpose here.')
  }
  #taking tumor paired samples
  x_both_samples<-x[!is.na(x) & !is.na(y)]
  #taking normal paired samples
  y_both_samples<-y[!is.na(x) & !is.na(y)]
  #taking  tumor  unpaired samples 
  x_samples_part<-x[!is.na(x) & is.na(y)]
  #taking  normal unpaired samples
  y_samples_part<-y[is.na(x) & !is.na(y)]
  #number of paired sample
  n1<-length(x_both_samples)
  #number of tumor unpaired samples 
  n2<-length(x_samples_part)
  #number of normal unpaired samples
  n3<-length(y_samples_part)
  #Considering different scenarios
  #when we do not have enough paired and unpaired samples
  if (n1<2 & (n2<2 | n3<2)){
    stop('Do not have enough data to perform the test')
  }
  var.eq=TRUE
  if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value<a){
    var.eq = FALSE
  }
  else if (var.test(c(x_both_samples,x_samples_part),c(y_both_samples,y_samples_part))$p.value>=a){
    var.eq = TRUE
  }
  #when we do not have enough paired samples  
  if (n1<2 & (n2>=2 & n3>=2)){
    warning('do not have enough paired samples, so applying unpaired t.test()')
    p_c<-t.test(x_samples_part,y_samples_part,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  else if (n1>=2 & (n2<2 | n3<2)){
    warning('do not have enough unpaired samples, so applying paired t.test()')
    p_c<-t.test(x_both_samples,y_both_samples,paired=TRUE,alternative=alternative, var.equal = var.eq)$p.value
    return(pvalue=p_c)
  }
  #mean of tumor unpaired samples
  T<-mean(x_samples_part)
  #mean of unpaired normal samples
  N<-mean(y_samples_part)
  #mean of tumor paired samples
  T1<-mean(x_both_samples)
  #mean of normal paired samples
  N1<-mean(y_both_samples)
  #standard deviation of tumor unpaired samples
  SD_T<-sd(x_samples_part)
  #standard deviation of normal unpaired samples
  SD_N<-sd(y_samples_part)
  #standard deviation of tumor paired samples
  SD_T1<-sd(x_both_samples)
  #standard deviation of normal paired samples
  SD_N1<-sd(y_both_samples)
  #covariance of tumor and normal paired samples
  SD_TN1<-sum((x_both_samples-T1)*(y_both_samples-N1))/(n1-1)
  #SD_TN1<-cov(x_both_samples,y_both_samples)
  r<-SD_TN1/(SD_T1*SD_N1)
  f_<-(n1*(n1+n3+n2*r))/((n1+n2)*(n1+n3)-n2*n3*r^2)
  g_<-(n1*(n1+n2+n3*r))/((n1+n2)*(n1+n3)-n2*n3*r^2)
  s<-(SD_T1^2*(n1-1)+SD_N1^2*(n1-1)+(1+r^2)*(SD_T^2*(n2-1)+SD_N^2*(n3-1)))/(2*(n1-1)+(1+r^2)*(n2+n3-2))
  V1_<-s*(2*n1*(1-r)+(n2+n3)*(1-r^2))/((n1+n2)*(n1+n3)-n2*n3*r^2)
  #MLE based test statistics for homoscedasticity by Ekbohm
  Z_E<-(f_*(T1-T)-g_*(N1-N)+T-N)/sqrt(V1_)
  if (alternative=="two.sided"){
    pvalue_Z_E<-pt(Z_E,df=n1,lower.tail = FALSE)
    if (pvalue_Z_E<1/2){
      pvalue_ZE<-2*pvalue_Z_E
    }
    else {
      pvalue_ZE<-2*(1-pvalue_Z_E)
    }
  }
  else if (alternative=="greater"){
    pvalue_ZE<-pt(Z_E,df=n1,lower.tail = FALSE)
  }
  else if (alternative=="less"){
    pvalue_ZE<-pt(Z_E,df=n1,lower.tail = TRUE)
  }
  if (var.eq==FALSE){
    warning('Ekbohm test works best for the homoscedastic case, that is when variances are equal. The data provided here has unequal varaince')
    return(pvalue_ZE)
  }
  else {
    return(pvalue_ZE)
  }
}
