######################################################## 0 Setup
library(tidyverse)
library(forestploter)
require(grid)
library(ggplot2)
library(dplyr)
library(readr)

# 设置工作目录和疾病版本
directory <- "D:/studio/dietary/"  # 工作目录的基础路径
setwd(directory)

# 设置疾病名称(全因死亡all-cause设为"all")
version <- "all"  # 疾病名称，这里为“CRC”，#疾病名称(包括:EC、SC、CRC、LC、PC、GCBTC、LCOCC) #all

# 创建各种数据和图表所需的子目录
# 注释掉的部分是用于创建目录的代码，当前没有被执行
dir.create(paste0("data/GDD/",version,"/"))
dir.create(paste0("data/RR/",version,"/"))
dir.create(paste0("data/gamma/",version,"/"))
dir.create(paste0("data/PAF/",version,"/"))
dir.create(paste0("figure/",version,"/"))
dir.create(paste0("table/",version,"/"))

code_directory     <- paste0(directory, "code/")
GDD_directory      <- paste0(directory, "data/GDD/")
RR_directory       <- paste0(directory, "data/RR/",version,"/")
gamma_directory    <- paste0(directory, "data/gamma/",version,"/")
tmred_directory    <- paste0(directory, "data/TMRED/")
pop_directory      <- paste0(directory, "data/Pop/")
paf_directory      <- paste0(directory, "data/PAF/",version,"/")
figure_directory   <- paste0(directory, "figure/",version,"/")
table_directory    <- paste0(directory, "table/",version,"/")

tm <- forest_theme(base_size = 10,
                   refline_lty = "dashed",   #参考线类型
                   refline_col = "#FF9796",
                   ci_pch = c(18),
                   ci_col = c("black"),
                   ci_Theight = 0.15)

ageRR_alloca=function(XXX_RR,XXX_RRlower,XXX_RRupper,XXX_age,XXX_id){
  age_pattern=tribble(~age_group,~median_age,~RR_95CI,"25-34",30,"0.92(0.87,0.97)","35-44",40,"0.92(0.87,0.97)",
                      "45-54",50,"0.93(0.89,0.97)","55-64",60,"0.94(0.91,0.98)","65-74",70,"0.95(0.92,0.98)","75+",80,"0.97(0.96,0.99)")
  age_pattern=age_pattern %>% mutate(lowerCI=as.numeric(str_extract(RR_95CI,"(?<=\\().*(?=,)")),
                                     upperCI=as.numeric(str_extract(RR_95CI,"(?<=,).*(?=\\))")),
                                     RR=as.numeric(str_extract(RR_95CI,"^.*(?=\\()")),
                                     logRR=log(RR),loglowerCI=log(lowerCI),logupperCI=log(upperCI))
  age_group=tibble(age_group=c(paste(seq(20,90,by=5),seq(24,94,by=5),sep="-"),"95+"),
                   median_age=seq(22.5,97.5,by=5),RR=c(rep(NA,16)),lowerCI=c(rep(NA,16)),upperCI=c(rep(NA,16)))
  age_group[str_which(age_group$median_age,XXX_age),c(3)]=XXX_RR
  age_group[str_which(age_group$median_age,XXX_age),c(4)]=XXX_RRlower
  age_group[str_which(age_group$median_age,XXX_age),c(5)]=XXX_RRupper
  age_group=age_group %>% mutate(logRR=log(RR),loglowerCI=log(lowerCI),logupperCI=log(upperCI))
  ggplot(age_pattern)+
    geom_line(aes(x=median_age,y=logRR))+geom_point(aes(x=median_age,y=logRR))+
    geom_line(aes(x=median_age,y=loglowerCI),color="red")+geom_point(aes(x=median_age,y=loglowerCI),color="red")+
    geom_line(aes(x=median_age,y=logupperCI),color="red")+geom_point(aes(x=median_age,y=logupperCI),color="red")
  
  age_group$logRR_assign=NA
  age_group$logRR_assign[str_which(age_group$median_age,"32.5")]=age_pattern$logRR[1]
  age_group$logRR_assign[str_which(age_group$median_age,"37.5")]=age_pattern$logRR[1]
  age_group$logRR_assign[5]=(age_pattern$logRR[3]-age_pattern$logRR[2])*
    (age_group$median_age[5]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$logRR[2]
  age_group$logRR_assign[6]=(age_pattern$logRR[3]-age_pattern$logRR[2])*
    (age_group$median_age[6]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$logRR[2]
  age_group$logRR_assign[7]=(age_pattern$logRR[4]-age_pattern$logRR[3])*
    (age_group$median_age[7]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$logRR[3]
  age_group$logRR_assign[8]=(age_pattern$logRR[4]-age_pattern$logRR[3])*
    (age_group$median_age[8]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$logRR[3]
  age_group$logRR_assign[9]=(age_pattern$logRR[5]-age_pattern$logRR[4])*
    (age_group$median_age[9]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$logRR[4]
  age_group$logRR_assign[10]=(age_pattern$logRR[5]-age_pattern$logRR[4])*
    (age_group$median_age[10]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$logRR[4]
  age_group$logRR_assign[11]=(age_pattern$logRR[6]-age_pattern$logRR[5])*
    (age_group$median_age[11]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$logRR[5]
  age_group$logRR_assign[12]=(age_pattern$logRR[6]-age_pattern$logRR[5])*
    (age_group$median_age[12]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$logRR[5]
  
  
  age_group$loglowerCI_assign=NA
  age_group$loglowerCI_assign[str_which(age_group$median_age,"32.5")]=age_pattern$loglowerCI[1]
  age_group$loglowerCI_assign[str_which(age_group$median_age,"37.5")]=age_pattern$loglowerCI[1]
  age_group$loglowerCI_assign[5]=(age_pattern$loglowerCI[3]-age_pattern$loglowerCI[2])*
    (age_group$median_age[5]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$loglowerCI[2]
  age_group$loglowerCI_assign[6]=(age_pattern$loglowerCI[3]-age_pattern$loglowerCI[2])*
    (age_group$median_age[6]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$loglowerCI[2]
  age_group$loglowerCI_assign[7]=(age_pattern$loglowerCI[4]-age_pattern$loglowerCI[3])*
    (age_group$median_age[7]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$loglowerCI[3]
  age_group$loglowerCI_assign[8]=(age_pattern$loglowerCI[4]-age_pattern$loglowerCI[3])*
    (age_group$median_age[8]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$loglowerCI[3]
  age_group$loglowerCI_assign[9]=(age_pattern$loglowerCI[5]-age_pattern$loglowerCI[4])*
    (age_group$median_age[9]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$loglowerCI[4]
  age_group$loglowerCI_assign[10]=(age_pattern$loglowerCI[5]-age_pattern$loglowerCI[4])*
    (age_group$median_age[10]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$loglowerCI[4]
  age_group$loglowerCI_assign[11]=(age_pattern$loglowerCI[6]-age_pattern$loglowerCI[5])*
    (age_group$median_age[11]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$loglowerCI[5]
  age_group$loglowerCI_assign[12]=(age_pattern$loglowerCI[6]-age_pattern$loglowerCI[5])*
    (age_group$median_age[12]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$loglowerCI[5]
  
  
  age_group$logupperCI_assign=NA
  age_group$logupperCI_assign[str_which(age_group$median_age,"32.5")]=age_pattern$logupperCI[1]
  age_group$logupperCI_assign[str_which(age_group$median_age,"37.5")]=age_pattern$logupperCI[1]
  age_group$logupperCI_assign[5]=(age_pattern$logupperCI[3]-age_pattern$logupperCI[2])*
    (age_group$median_age[5]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$logupperCI[2]
  age_group$logupperCI_assign[6]=(age_pattern$logupperCI[3]-age_pattern$logupperCI[2])*
    (age_group$median_age[6]-age_pattern$median_age[2])/
    (age_pattern$median_age[3]-age_pattern$median_age[2])+age_pattern$logupperCI[2]
  age_group$logupperCI_assign[7]=(age_pattern$logupperCI[4]-age_pattern$logupperCI[3])*
    (age_group$median_age[7]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$logupperCI[3]
  age_group$logupperCI_assign[8]=(age_pattern$logupperCI[4]-age_pattern$logupperCI[3])*
    (age_group$median_age[8]-age_pattern$median_age[3])/
    (age_pattern$median_age[4]-age_pattern$median_age[3])+age_pattern$logupperCI[3]
  age_group$logupperCI_assign[9]=(age_pattern$logupperCI[5]-age_pattern$logupperCI[4])*
    (age_group$median_age[9]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$logupperCI[4]
  age_group$logupperCI_assign[10]=(age_pattern$logupperCI[5]-age_pattern$logupperCI[4])*
    (age_group$median_age[10]-age_pattern$median_age[4])/
    (age_pattern$median_age[5]-age_pattern$median_age[4])+age_pattern$logupperCI[4]
  age_group$logupperCI_assign[11]=(age_pattern$logupperCI[6]-age_pattern$logupperCI[5])*
    (age_group$median_age[11]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$logupperCI[5]
  age_group$logupperCI_assign[12]=(age_pattern$logupperCI[6]-age_pattern$logupperCI[5])*
    (age_group$median_age[12]-age_pattern$median_age[5])/
    (age_pattern$median_age[6]-age_pattern$median_age[5])+age_pattern$logupperCI[5]
  
  
  line=lm(logRR~median_age,age_pattern)
  age_group$extrapolation=predict(line,age_group)
  age_group=age_group %>% mutate(logRR_assign=ifelse(is.na(logRR_assign),extrapolation,logRR_assign))
  
  line=lm(loglowerCI~median_age,age_pattern)
  age_group$extrapolation_lower=predict(line,age_group)
  age_group=age_group %>% mutate(loglowerCI_assign=ifelse(is.na(loglowerCI_assign),extrapolation_lower,loglowerCI_assign))
  
  line=lm(logupperCI~median_age,age_pattern)
  age_group$extrapolation_upper=predict(line,age_group)
  age_group=age_group %>% mutate(logupperCI_assign=ifelse(is.na(logupperCI_assign),extrapolation_upper,logupperCI_assign))
  
  colnames(age_group)
  for (i in seq_along(age_group$age_group)) {
    age_group$logRR[i]=age_group$logRR[XXX_id]*age_group$logRR_assign[i]/age_group$logRR_assign[XXX_id]
  }
  for (i in seq_along(age_group$age_group)) {
    age_group$loglowerCI[i]=age_group$loglowerCI[XXX_id]*age_group$loglowerCI_assign[i]/age_group$loglowerCI_assign[XXX_id]
  }
  for (i in seq_along(age_group$age_group)) {
    age_group$logupperCI[i]=age_group$logupperCI[XXX_id]*age_group$logupperCI_assign[i]/age_group$logupperCI_assign[XXX_id]
  }
  ## 抽样
  outcome=list()
  outcome_final=tibble(lowerUI=rep(NA,16),upperUI=rep(NA,16))
  set.seed(20230826)
  for (i in 1:16) {
    outcome[[i]]=runif(1000, min = age_group$loglowerCI[i], max = age_group$logupperCI[i])
    # outcome[[i]]=rnorm(1000, mean=age_group$logRR[i], sd=log(exp(age_group$logupperCI[i])/exp(age_group$logupperCI[i]))/(2*1.96))
    outcome_final$lowerUI[i]=quantile(outcome[[i]], 0.025)
    outcome_final$upperUI[i]=quantile(outcome[[i]], 0.975)
  }
  outcome_final=outcome_final %>% mutate(median_age=seq(22.5,97.5,by=5))
  
  age_group=age_group %>% left_join(outcome_final,by="median_age")
  age_group$RR_final=round(exp(age_group$logRR),2)
  age_group$lowerUI_final=round(exp(age_group$lowerUI),2)
  age_group$upperUI_final=round(exp(age_group$upperUI),2)
  age_group=age_group %>% select(age_group,median_age,RR_final,lowerUI_final,upperUI_final) %>% 
    mutate(RR_95UI=paste0(RR_final," ","(",lowerUI_final,","," ",upperUI_final,")"))
  return(age_group)
}

gammaGDD1000=function(GDD_ditery=potatoes,YEAR=YEAR,iters=iters){
  # 注意此时有两个年份
  potatoes_wave=potatoes
  potatoes=potatoes_wave %>% filter(year==YEAR)
  
  # 1 调整能量已经完成
  # 2 拟合lm
  regression_result=lm(log(sd) ~ log(median), data = potatoes)
  beta0_hat=regression_result$coefficients[1];beta1_hat=regression_result$coefficients[2];sigma_hat=summary(regression_result)$sigma
  
  # 3 求未加误差的1000 ln(sd)
  n = length(potatoes$age_average)
  # x_hat=tibble(xjk=sample(min(potatoes$median):max(potatoes$median),iters,replace = T))
  # normal distribution for the log mean
  set.seed(20230826)
  x_hat=tibble(xjk=exp(rnorm(n=iters,mean=mean(log(potatoes$median)), sd=mean(log(potatoes$upperci_95/potatoes$lowerci_95)/(2*1.96)))))
  x_hat %>% mutate(max(xjk),min(xjk),mean(xjk))
  potatoes %>% mutate(max(median),min(median),mean(median))
  x_hat$y_x=beta1_hat*log(x_hat$xjk)+beta0_hat
  
  # 4 1000 ln(sd)
  t_sample = rt(n=iters,df=n-1);ln_sd_hat <- x_hat$`y_x`+ sigma_hat*sqrt(1+(1/n)*t_sample)
  
  # 5 1000 gamma-distribution
  sd_result=exp(ln_sd_hat)
  shape=(x_hat$xjk/sd_result)^2;rate=(x_hat$xjk/sd_result^2)
  lowerpoint=function(x){quantile(x,0.025)};upperpoint=function(x){quantile(x,0.975)}
  gamma_distribution=tibble(shape=mean(shape),rate=mean(rate))
  gamma_distribution$shape_upper=quantile(shape,0.975);gamma_distribution$shape_lower=quantile(shape,0.025)
  gamma_distribution$rate_upper=quantile(rate,0.975);gamma_distribution$rate_lower=quantile(rate,0.025)
  gamma_distribution_YEAR=gamma_distribution %>% mutate(year=YEAR)
  gamma1000shape=tibble(shape=shape) %>% arrange(shape) %>% mutate(year=YEAR)
  gamma1000rate=tibble(rate=rate) %>% arrange(rate) %>% mutate(year=YEAR)
  gamma1000_YEAR=cbind(gamma1000shape,gamma1000rate)[-4]
  return(list(gamma1000_YEAR,gamma_distribution_YEAR))
}

gammaGDD1000NA=function(GDD_ditery=potatoes,YEAR=YEAR,iters=iters){
  # 注意此时有两个年份
  potatoes_wave=na.omit(potatoes)
  potatoes=potatoes_wave %>% filter(year==YEAR)
  
  # 1 调整能量已经完成
  # 2 拟合lm
  regression_result=lm(log(sd) ~ log(median), data = potatoes)
  beta0_hat=regression_result$coefficients[1];beta1_hat=regression_result$coefficients[2];sigma_hat=summary(regression_result)$sigma
  
  # 3 求未加误差的1000 ln(sd)
  n = length(potatoes$age_average)
  # x_hat=tibble(xjk=sample(min(potatoes$median):max(potatoes$median),iters,replace = T))
  # normal distribution for the log mean
  set.seed(20230826)
  x_hat=tibble(xjk=exp(rnorm(n=iters,mean=mean(log(potatoes$median)), sd=mean(log(potatoes$upperci_95/potatoes$lowerci_95)/(2*1.96)))))
  x_hat %>% mutate(max(xjk),min(xjk),mean(xjk))
  potatoes %>% mutate(max(median),min(median),mean(median))
  x_hat$y_x=beta1_hat*log(x_hat$xjk)+beta0_hat
  
  # 4 1000 ln(sd)
  t_sample = rt(n=iters,df=n-1);ln_sd_hat <- x_hat$`y_x`+ sigma_hat*sqrt(1+(1/n)*t_sample)
  
  # 5 1000 gamma-distribution
  sd_result=exp(ln_sd_hat)
  shape=(x_hat$xjk/sd_result)^2;rate=(x_hat$xjk/sd_result^2)
  lowerpoint=function(x){quantile(x,0.025)};upperpoint=function(x){quantile(x,0.975)}
  gamma_distribution=tibble(shape=mean(shape),rate=mean(rate))
  gamma_distribution$shape_upper=quantile(shape,0.975,na.rm = T);gamma_distribution$shape_lower=quantile(shape,0.025,na.rm = T)
  gamma_distribution$rate_upper=quantile(rate,0.975,na.rm = T);gamma_distribution$rate_lower=quantile(rate,0.025,na.rm = T)
  gamma_distribution_YEAR=gamma_distribution %>% mutate(year=YEAR)
  gamma1000shape=tibble(shape=shape) %>% arrange(shape) %>% mutate(year=YEAR)
  gamma1000rate=tibble(rate=rate) %>% arrange(rate) %>% mutate(year=YEAR)
  gamma1000_YEAR=cbind(gamma1000shape,gamma1000rate)[-4]
  return(list(gamma1000_YEAR,gamma_distribution_YEAR))
}

PAFcalc=function(XXX_year,XXX_iters,XXX_tmred_mean,XXX_unit,dietary_data,gamma_distribution1000,XXX_riskfactor){
  gamma_distribution1000 = gamma_distribution1000 %>% filter(year==XXX_year)
  gamma_distribution1000$year=as.numeric(gamma_distribution1000$year)
  
  # # 2. RR(X)
  # 2. calculate n categories intervals
  n_intervals=tibble(from=c(0,pnorm((seq(2,121)-2)*0.1-6)),to=c(pnorm((seq(1,121)-1)*0.1-6)))
  gamma_intervals_from=list();gamma_intervals_to=list()
  iters=XXX_iters
  for (i in 1:(iters)){gamma_intervals_from[[i]]=c(qgamma(n_intervals$from,shape = gamma_distribution1000$shape[i],rate = gamma_distribution1000$rate[i]),NA)
  gamma_intervals_to[[i]]=c(qgamma(n_intervals$to,shape = gamma_distribution1000$shape[i],rate = gamma_distribution1000$rate[i]),NA)}
  gamma_intervals=list();gamma_intervals_diff=list()
  
  for (i in 1:(iters)) {gamma_intervals_diff[[i]]=gamma_intervals_from[[i]]-lag(gamma_intervals_from[[i]])
  gamma_intervals[[i]]=tibble(gamma_intervals_from=gamma_intervals_from[[i]],gamma_intervals_diff=lead(gamma_intervals_diff[[i]]),gamma_intervals_to=gamma_intervals_to[[i]])
  gamma_intervals[[i]]=gamma_intervals[[i]] %>% slice(1:121);gamma_intervals[[i]]=gamma_intervals[[i]] %>% mutate(year_id=i)}
  gamma_intervals=do.call(rbind, gamma_intervals)
  rm(gamma_intervals_from,gamma_intervals_to,gamma_intervals_diff)
  gamma_intervals=gamma_intervals %>% mutate(gamma_intervals_diff=ifelse(is.na(gamma_intervals_diff),(gamma_intervals_to-gamma_intervals_from),gamma_intervals_diff))
  options(scipen = 999)
  
  # 3.calculate P(x)
  gamma_distribution1000=gamma_distribution1000 %>% rownames_to_column("year_id")
  gamma_distribution1000$year_id=as.numeric(gamma_distribution1000$year_id)
  gamma_intervals=gamma_intervals %>% left_join(gamma_distribution1000,by=c("year_id"))
  
  gamma_intervals$p=dgamma((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2,shape=gamma_intervals$shape,rate=gamma_intervals$rate)
  gamma_intervals$p=gamma_intervals$p*gamma_intervals$gamma_intervals_diff
  
  # 4.calculate RR(x)=exp(beta(x-y(x))),y(x)=pgamma(qgamma(x))
  ##################### sample
  lnRR_mean=log(dietary_data$RR_final)
  lnRR_sd=log(dietary_data$upperUI_final/dietary_data$lowerUI_final)/(2*1.96)
  RR_result=list()
  set.seed(20230826)
  for (i in 1:16) {RR_result[[i]]=tibble(RR=rnorm(n=iters,mean=lnRR_mean[i], sd=lnRR_sd[i]),agegroup_id=i,agegroup_ind=1:XXX_iters)
  RR_result[[i]]$RR=exp(RR_result[[i]]$RR)}
  RR_result = do.call(rbind,RR_result)
  
  XXX_tmred_mean=XXX_tmred_mean
  gamma_intervals$tmred_mean=XXX_tmred_mean
  gamma_intervals$tmred_sd=0.1*gamma_intervals$tmred_mean
  gamma_intervals$tmred_shape=(gamma_intervals$tmred_mean/gamma_intervals$tmred_sd)^2
  gamma_intervals$tmred_rate=(gamma_intervals$tmred_mean/gamma_intervals$tmred_sd^2)
  gamma_intervals$tmred_point=qgamma(pgamma((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2,shape=gamma_intervals$shape,rate=gamma_intervals$rate),
                                     shape=gamma_intervals$tmred_shape,rate=gamma_intervals$tmred_rate)
  if (XXX_tmred_mean==0){
    gamma_intervals$difference=((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2)
  }else{
    gamma_intervals$difference=((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2)-gamma_intervals$tmred_point
  }
  gamma_intervals_list=list()
  for (i in seq_along(RR_result$RR)) {gamma_intervals_list[[i]]=gamma_intervals %>% 
    # filter(median_age==(ceiling(i/300)-1)*5+22.5) %>% 
    mutate(RR_sample=RR_result$RR[i],agegroup_id=RR_result$agegroup_id[i],agegroup_ind=RR_result$agegroup_ind[i]) %>% select(year_id,p,difference,RR_sample,agegroup_id,agegroup_ind)}
  
  if (XXX_riskfactor>1){
    for (i in seq_along(gamma_intervals_list)){gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% mutate(rr=ifelse((difference)>=0,exp(log(RR_sample)/XXX_unit*(difference)),1),
                                                                                                              rr_p=(rr-1)*p) %>% select(rr_p,year_id,agegroup_id,agegroup_ind)}
  }else{
    for (i in seq_along(gamma_intervals_list)){gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% mutate(rr=ifelse((difference)<=0,exp(log(RR_sample)/XXX_unit*(difference)),1),
                                                                                                              rr_p=(rr-1)*p) %>% select(rr_p,year_id,agegroup_id,agegroup_ind)}
  }
  for (i in seq_along(gamma_intervals_list)) {gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% group_by(year_id) %>% 
    mutate(sum = ifelse(row_number() == 1, 0, cumsum(rr_p))) %>%
    filter(row_number() == n()) %>%
    select(sum,agegroup_id,agegroup_ind) %>% mutate(PAF=sum/(sum+1)*100)}
  
  for (i in seq_along(gamma_intervals_list)){
    gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% select(-sum)
  }
  gamma_intervals_list=do.call(rbind,gamma_intervals_list)
  
  gamma_intervals_list=gamma_intervals_list %>% left_join(age_percentage)
  gamma_intervals_list=gamma_intervals_list %>% mutate(PAF_pop=pop*PAF) %>% select(-pop,-PAF)
  PAF=gamma_intervals_list %>% group_by(year_id,agegroup_ind) %>% summarise(sum(PAF_pop))
  
  PAF_final=tibble(quantile(PAF$`sum(PAF_pop)`,0.025),mean(PAF$`sum(PAF_pop)`),quantile(PAF$`sum(PAF_pop)`,0.975),quantile(PAF$`sum(PAF_pop)`,0.5))
  colnames(PAF_final)=c("PAF_lower","PAF_mean","PAF_upper","PAF_50")
  PAF_final=PAF_final %>% mutate(PAFUI_m=paste0(round(PAF_mean,2)," ","(",round(PAF_lower,2),","," ",round(PAF_upper,2),")"),
                                 PAFUI50=paste0(round(PAF_50,2)," ","(",round(PAF_lower,2),","," ",round(PAF_upper,2),")"))
  return(PAF_final)
}

PAFcalc1000=function(XXX_year,XXX_iters,XXX_tmred_mean,XXX_unit,dietary_data,gamma_distribution1000,XXX_riskfactor){
  gamma_distribution1000 = gamma_distribution1000 %>% filter(year==XXX_year)
  gamma_distribution1000$year=as.numeric(gamma_distribution1000$year)
  
  # # 2. RR(X)
  # 2. calculate n categories intervals
  n_intervals=tibble(from=c(0,pnorm((seq(2,121)-2)*0.1-6)),to=c(pnorm((seq(1,121)-1)*0.1-6)))
  gamma_intervals_from=list();gamma_intervals_to=list()
  iters=XXX_iters
  for (i in 1:(iters)){gamma_intervals_from[[i]]=c(qgamma(n_intervals$from,shape = gamma_distribution1000$shape[i],rate = gamma_distribution1000$rate[i]),NA)
  gamma_intervals_to[[i]]=c(qgamma(n_intervals$to,shape = gamma_distribution1000$shape[i],rate = gamma_distribution1000$rate[i]),NA)}
  gamma_intervals=list();gamma_intervals_diff=list()
  
  for (i in 1:(iters)) {gamma_intervals_diff[[i]]=gamma_intervals_from[[i]]-lag(gamma_intervals_from[[i]])
  gamma_intervals[[i]]=tibble(gamma_intervals_from=gamma_intervals_from[[i]],gamma_intervals_diff=lead(gamma_intervals_diff[[i]]),gamma_intervals_to=gamma_intervals_to[[i]])
  gamma_intervals[[i]]=gamma_intervals[[i]] %>% slice(1:121);gamma_intervals[[i]]=gamma_intervals[[i]] %>% mutate(year_id=i)}
  gamma_intervals=do.call(rbind, gamma_intervals)
  rm(gamma_intervals_from,gamma_intervals_to,gamma_intervals_diff)
  gamma_intervals=gamma_intervals %>% mutate(gamma_intervals_diff=ifelse(is.na(gamma_intervals_diff),(gamma_intervals_to-gamma_intervals_from),gamma_intervals_diff))
  options(scipen = 999)
  
  # 3.calculate P(x)
  gamma_distribution1000=gamma_distribution1000 %>% rownames_to_column("year_id")
  gamma_distribution1000$year_id=as.numeric(gamma_distribution1000$year_id)
  gamma_intervals=gamma_intervals %>% left_join(gamma_distribution1000,by=c("year_id"))
  
  gamma_intervals$p=dgamma((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2,shape=gamma_intervals$shape,rate=gamma_intervals$rate)
  gamma_intervals$p=gamma_intervals$p*gamma_intervals$gamma_intervals_diff
  
  # 4.calculate RR(x)=exp(beta(x-y(x))),y(x)=pgamma(qgamma(x))
  ##################### sample
  lnRR_mean=log(dietary_data$RR_final)
  lnRR_sd=log(dietary_data$upperUI_final/dietary_data$lowerUI_final)/(2*1.96)
  RR_result=list()
  set.seed(20230826)
  for (i in 1:16) {RR_result[[i]]=tibble(RR=rnorm(n=iters,mean=lnRR_mean[i], sd=lnRR_sd[i]),agegroup_id=i,agegroup_ind=1:XXX_iters)
  RR_result[[i]]$RR=exp(RR_result[[i]]$RR)}
  RR_result = do.call(rbind,RR_result)
  
  XXX_tmred_mean=XXX_tmred_mean
  gamma_intervals$tmred_mean=XXX_tmred_mean
  gamma_intervals$tmred_sd=0.1*gamma_intervals$tmred_mean
  gamma_intervals$tmred_shape=(gamma_intervals$tmred_mean/gamma_intervals$tmred_sd)^2
  gamma_intervals$tmred_rate=(gamma_intervals$tmred_mean/gamma_intervals$tmred_sd^2)
  gamma_intervals$tmred_point=qgamma(pgamma((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2,shape=gamma_intervals$shape,rate=gamma_intervals$rate),
                                     shape=gamma_intervals$tmred_shape,rate=gamma_intervals$tmred_rate)
  if (XXX_tmred_mean==0){
    gamma_intervals$difference=((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2)
  }else{
    gamma_intervals$difference=((gamma_intervals$gamma_intervals_from+gamma_intervals$gamma_intervals_to)/2)-gamma_intervals$tmred_point
  }
  gamma_intervals_list=list()
  for (i in seq_along(RR_result$RR)) {gamma_intervals_list[[i]]=gamma_intervals %>% 
    # filter(median_age==(ceiling(i/300)-1)*5+22.5) %>% 
    mutate(RR_sample=RR_result$RR[i],agegroup_id=RR_result$agegroup_id[i],agegroup_ind=RR_result$agegroup_ind[i]) %>% select(year_id,p,difference,RR_sample,agegroup_id,agegroup_ind)}
  
  if (XXX_riskfactor>1){
    for (i in seq_along(gamma_intervals_list)){gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% mutate(rr=ifelse((difference)>=0,exp(log(RR_sample)/XXX_unit*(difference)),1),
                                                                                                              rr_p=(rr-1)*p) %>% select(rr_p,year_id,agegroup_id,agegroup_ind)}
  }else{
    for (i in seq_along(gamma_intervals_list)){gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% mutate(rr=ifelse((difference)<=0,exp(log(RR_sample)/XXX_unit*(difference)),1),
                                                                                                              rr_p=(rr-1)*p) %>% select(rr_p,year_id,agegroup_id,agegroup_ind)}
  }
  for (i in seq_along(gamma_intervals_list)) {gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% group_by(year_id) %>% 
    mutate(sum = ifelse(row_number() == 1, 0, cumsum(rr_p))) %>%
    filter(row_number() == n()) %>%
    select(sum,agegroup_id,agegroup_ind) %>% mutate(PAF=sum/(sum+1)*100)}
  
  for (i in seq_along(gamma_intervals_list)){
    gamma_intervals_list[[i]]=gamma_intervals_list[[i]] %>% select(-sum)
  }
  gamma_intervals_list=do.call(rbind,gamma_intervals_list)
  
  gamma_intervals_list=gamma_intervals_list %>% left_join(age_percentage)
  gamma_intervals_list=gamma_intervals_list %>% mutate(PAF_pop=pop*PAF) %>% select(-pop,-PAF)
  PAF=gamma_intervals_list %>% group_by(year_id,agegroup_ind) %>% summarise(sum(PAF_pop))
  
  return(PAF)
}

read_PAF_csv=function(df){a=df %>% str_split_fixed("_",n=2);b=a[1] %>% str_split_fixed("/",n=4)
PAF_df=read_csv(df) %>% mutate(dietray=b[4]);return(PAF_df)}

joint_PAF=function(files){
  PAF=map_dfc(files, read_csv);PAF_factor <- colnames(PAF)[grepl("^PAF", colnames(PAF))]
  PAF=PAF %>% select(`id...3`,`year...2`,starts_with("PAF")) %>% rename(id=`id...3`,year=`year...2`)
  
  PAF = PAF %>% rowwise() %>% mutate(jointPAF = (1 - (prod(1 - c_across(starts_with("PAF")) / 100))) * 100) %>% ungroup()
  selected_columns_1990=PAF[str_which(PAF$year,"1990"),c(PAF_factor,"jointPAF")]
  selected_columns_2018=PAF[str_which(PAF$year,"2018"),c(PAF_factor,"jointPAF")]
  PAF_lower_1990 <- tibble(PAF_lower=apply(selected_columns_1990, 2, quantile, probs = 0.025),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=1990)
  PAF_upper_1990 <- tibble(PAF_upper=apply(selected_columns_1990, 2, quantile, probs = 0.975),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=1990)
  PAF_50_1990 <- tibble(PAF_50=apply(selected_columns_1990, 2, quantile, probs = 0.5),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=1990)
  PAF_mean_1990 <- tibble(PAF_50=apply(selected_columns_1990, 2, mean),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=1990)
  PAF_lower_2018 <- tibble(PAF_lower=apply(selected_columns_2018, 2, quantile, probs = 0.025),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=2018)
  PAF_upper_2018 <- tibble(PAF_upper=apply(selected_columns_2018, 2, quantile, probs = 0.975),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=2018)
  PAF_50_2018 <- tibble(PAF_50=apply(selected_columns_2018, 2, quantile, probs = 0.5),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=2018)
  PAF_mean_2018 <- tibble(PAF_50=apply(selected_columns_2018, 2, mean),dietary=c(PAF_factor,"jointPAF")) %>% mutate(year=2018)
  
  PAF_lower=rbind(PAF_lower_1990,PAF_lower_2018)
  PAF_upper=rbind(PAF_upper_1990,PAF_upper_2018)
  PAF_50=rbind(PAF_50_1990,PAF_50_2018)
  PAF_mean=rbind(PAF_mean_1990,PAF_mean_2018)
  
  PAF_outcomes=PAF_lower %>% left_join(PAF_upper) %>% left_join(PAF_50) %>% left_join(PAF_mean)
  PAF_outcomes=PAF_outcomes[,c(2,3,1,4:max(seq_along(PAF_outcomes)))]
  return(PAF_outcomes)
}