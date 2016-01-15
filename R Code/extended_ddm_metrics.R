# extended_ddm_metrics
#
# Written by Patrick Simen, Jan 3, 2016
# This code is a transcription into R of the Matlab code written by Vaibhav Srivastava, 
# accompanying the publication:
# "Explicit moments of decision times for single- and double-threshold drift-diffusion processes"
# Submitted to Journal of Mathematical Psychology

# Import the basic ddm_metrics function and some hyperbolic functions:
source('ddm_metrics.R')

#function [error, mean_RT, var_RT, t_RT, CV, skew]= extended_ddm_metrics(a_mean, sigma_a, s, x0_mean, delta, z)
extended_ddm_metrics <- function(a_mean, sigma_a, s, x0_mean, delta, z) {
  
  # Input: 
  # a_mean =mean drift rate, sigma_a = variance in drift rate,
  # s=diffusion rate,
  # x0_mean= mean initial condition, 
  #delta = support for initial condition (x0-delta, x0+delta), 
  #z= threshold
  
  #Output
  # error = error rate
  # mean_RT = mean decision time
  # var_RT = variance of decision time
  # t_RT = third central moment of decision time 
  # CV = coefficient of variation of decision time 
  # skew = skewness of decision time 
  
  
  if (delta>z) {  
    delta=z;
  }
  
  
  if (delta==0) {
    x_set = x0_mean;
    prob_x = 1;
  } else {
    delta=pmax(0.0001,delta);
    #tmp_x = linspace(-delta,delta,100);
    tmp_x = array(seq(-delta,delta,len=100),dim=c(1,100))
    x_set = array()
    x_set = x0_mean + tmp_x;
    #prob_x = 1/2/delta*ones(size(x_set));
    prob_x = matrix(1/2/delta,nrow(x_set),ncol(x_set));
  }
  
  
  
  #mean_RT_x= zeros(size(x_set));
  #sec_moment_x= zeros(size(x_set));
  #third_moment_x= zeros(size(x_set));
  #err_x=zeros(size(x_set));
  
  mean_RT_x = matrix(0,nrow(x_set),ncol(x_set))
  sec_moment_x = matrix(0,nrow(x_set),ncol(x_set))
  third_moment_x = matrix(0,nrow(x_set),ncol(x_set))
  err_x = matrix(0,nrow(x_set),ncol(x_set))
  
  
  for (jj in 1:length(x_set)) {
    
    x0 = x_set[jj];
    
    if ( sigma_a == 0 ) {
      a_supp=a_mean;
      
      #[err,m_RT, v_RT, t_RT] = ddm_metrics(a_supp,s,z,x0);
      results = ddm_metrics(a_supp,s,z,x0);
      err = results$err
      m_RT = results$m_RT
      v_RT = results$v_RT
      t_RT = results$t_RT
      
      err_x[jj] = err;
      
      mean_RT_x[jj] = m_RT;
      
      sec_moment_x[jj] = (v_RT + m_RT^2);
      
      third_moment_x[jj] = (t_RT + m_RT^3 + 3*m_RT*v_RT);
      
    } else {
      
      sigma_a=pmax(0.0001,sigma_a);
      
      #tmp_a=linspace(-5*sigma_a, 5*sigma_a, 100);
      tmp_a_seq = seq(-5*sigma_a, 5*sigma_a, len=100);
      tmp_a = array(tmp_a_seq,dim=c(1,length(tmp_a_seq)));
      
      a_supp = matrix()
      a_supp = a_mean + tmp_a;
      
      a_supp[abs(a_supp)<0.0001] = 0.0001;
      
      #prob_a = exp(max(-tmp_a^2/2/sigma_a^2, -100))/(sqrt(2*pi)*sigma_a);
      prob_a = matrix(exp(pmax(-tmp_a^2/2/sigma_a^2, -100))/(sqrt(2*pi)*sigma_a),nrow(tmp_a),ncol(tmp_a));
      
      
      #[err,m_RT, v_RT, t_RT] = ddm_metrics(a_supp,s,z,x0);
      results = ddm_metrics(a_supp,s,z,x0);
      err = results$err
      m_RT = results$m_RT
      v_RT = results$v_RT
      t_RT = results$t_RT
      
      err_x[jj] = err%*%t(prob_a)%*%(a_supp[length(a_supp)]-a_supp[1])/99;
  
      mean_RT_x[jj] = m_RT%*%t(prob_a)%*%(a_supp[length(a_supp)]-a_supp[1])/99;
      
      sec_moment_x[jj] = (v_RT + m_RT^2)%*%t(prob_a)%*%(a_supp[length(a_supp)]-a_supp[1])/99;
  
      third_moment_x[jj] = (t_RT + m_RT^3 + 3*m_RT*v_RT)%*%t(prob_a)%*%(a_supp[length(a_supp)]-a_supp[1])/99;
    }
  }
  
  if (delta == 0) {
    
    error = err_x;
    
    mean_RT = mean_RT_x;
    
    sec_mom_RT = sec_moment_x;
    
    third_mom_RT = third_moment_x;
    
    
  } else {
    error= err_x%*%t(prob_x)%*%(x_set[length(x_set)]-x_set[1])/99;
  
    mean_RT= mean_RT_x%*%t(prob_x)%*%(x_set[length(x_set)]-x_set[1])/99;
    
    sec_mom_RT= sec_moment_x%*%t(prob_x)%*%(x_set[length(x_set)]-x_set[1])/99;
  
    third_mom_RT= third_moment_x%*%t(prob_x)%*%(x_set[length(x_set)]-x_set[1])/99;
    
  }
  
  
  var_RT = abs(sec_mom_RT - mean_RT^2);
  
  t_RT =  abs(third_mom_RT -3*var_RT*mean_RT - mean_RT^3);
  
  CV= sqrt(abs(var_RT))/mean_RT;
  
  skew = t_RT/(abs(var_RT))^1.5;
  
  
  
  return_list = list()
  #error, mean_RT, var_RT, t_RT, CV, skew
  return_list$err = error
  return_list$m_RT = mean_RT
  return_list$v_RT = var_RT
  return_list$CV = CV
  return_list$t_RT = t_RT
  return_list$skew = skew
  
  return_list
}



