# extended_ddm_metrics_cond
#
# Written by Patrick Simen, Jan 3, 2016
# This code is a transcription into R of the Matlab code written by Vaibhav Srivastava, 
# accompanying the publication:
# "Explicit moments of decision times for single- and double-threshold drift-diffusion processes"
# Submitted to Journal of Mathematical Psychology

#function [mean_RTplus, mean_RTminus, var_RTplus,var_RTminus, t_RTplus,...
#t_RTminus, CVplus, CVminus, skewplus, skewminus]= extended_ddm_metrics_cond(a_mean, sigma_a, s, x0_mean, delta, z)
extended_ddm_metrics_cond <- function(a_mean, sigma_a, s, x0_mean, delta, z) {
  
  # Input: 
  # a_mean =mean drift rate, sigma_a = variance in drift rate,
  # s=diffusion rate,
  # x0_mean= mean initial condition, 
  # delta = support for initial condition (x0-delta, x0+delta), 
  # z= threshold
  
  # Output
  # mean_RTplus = mean decision time conditioned on correct decision
  # mean_RTminus = mean decision time conditioned on error
  # var_RTplus = variance of decision time conditioned on correct decision
  # var_RTminus = variance of decision time conditioned on error
  # t_RTplus = third central moment of decision time conditioned on correct decision
  # t_RTminus = third central moment of decision time conditioned on error
  # CVplus = coefficient of variation of decision time conditioned on correct decision
  # CVminus = coefficient of variation of decision time conditioned on error
  # skewplus = skewness of decision time conditioned on correct decision
  # skewminus = skewness of decision time conditioned on error
  
  
  if (delta>=z){
    delta=z-0.0001;
  }
  
  if (delta==0) {
    x_set =x0_mean;
    prob_x=1;
  } else {
    
    delta=pmax(0.0001,delta);
    #tmp_x = linspace(-delta,delta,100);
    tmp_x = array(seq(-delta,delta,len=100),dim=c(1,100))
    x_set = array()
    x_set = x0_mean + tmp_x;
    #prob_x = 1/2/delta*ones(size(x_set));
    prob_x = matrix(1/2/delta, nrow(x_set), ncol(x_set))
  }
  
  
  
  #   mean_RT_p_x= zeros(size(x_set));
  #   mean_RT_m_x= zeros(size(x_set));
  #   sec_moment_p_x= zeros(size(x_set));
  #   sec_moment_m_x= zeros(size(x_set));
  #   third_moment_p_x= zeros(size(x_set));
  #   third_moment_m_x= zeros(size(x_set));
  
  mean_RT_p_x = matrix(0,nrow(x_set),ncol(x_set))
  mean_RT_m_x = matrix(0,nrow(x_set),ncol(x_set))
  sec_moment_p_x = matrix(0,nrow(x_set),ncol(x_set))
  sec_moment_m_x = matrix(0,nrow(x_set),ncol(x_set))
  third_moment_p_x = matrix(0,nrow(x_set),ncol(x_set))
  third_moment_m_x = matrix(0,nrow(x_set),ncol(x_set))
  
  
  
  for (jj in 1:length(x_set)) {
    
    x0=x_set[jj];
    
    if (sigma_a ==0) {
      a_supp=a_mean;
      
      #[m_RTplus, m_RTminus, v_RTplus, v_RTminus, t_RTplus, t_RTminus] =ddm_metrics_cond(a_supp,s,z,x0);
      results = ddm_metrics_cond(a_supp,s,z,x0);
      m_RTplus = results$m_RTplus
      m_RTminus = results$m_RTminus
      v_RTplus = results$v_RTplus
      v_RTminus = results$v_RTminus
      t_RTplus = results$t_RTplus
      t_RTminus = results$t_RTminus
      
      
      
      mean_RT_p_x[jj] = m_RTplus;
      
      sec_moment_p_x[jj] = (v_RTplus + m_RTplus^2);
      
      third_moment_p_x[jj] = (t_RTplus + m_RTplus^3 + 3*m_RTplus*v_RTplus);
      
      mean_RT_m_x[jj] = m_RTminus;
      
      sec_moment_m_x[jj] = (v_RTminus + m_RTminus^2);
      
      third_moment_m_x[jj] = (t_RTminus + m_RTminus^3 + 3*m_RTminus*v_RTminus);
      
    } else {
      
      sigma_a=pmax(0.0001,sigma_a);
      
      #tmp_a=linspace(-5*sigma_a, 5*sigma_a,100);
      tmp_a = array(seq(-5*sigma_a,5*sigma_a,len=100),dim=c(1,100))
      
      a_supp = matrix()
      a_supp= a_mean + tmp_a ;
      
      a_supp[abs(a_supp)<0.0001]=0.0001;
      
      prob_a = matrix( exp(pmax(-tmp_a^2/2/sigma_a^2, -100))/(sqrt(2*pi)*sigma_a), nrow(tmp_a), ncol(tmp_a) )
      
      
      #[m_RTplus, m_RTminus, v_RTplus, v_RTminus, t_RTplus, t_RTminus] =ddm_metrics_cond(a_supp,s,z,x0);
      results =ddm_metrics_cond(a_supp,s,z,x0);
      m_RTplus = results$m_RTplus
      m_RTminus = results$m_RTminus
      v_RTplus = results$v_RTplus
      v_RTminus = results$v_RTminus
      t_RTplus = results$t_RTplus
      t_RTminus = results$t_RTminus
      
      
      
      mean_RT_p_x[jj] = m_RTplus %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
      sec_moment_p_x[jj] = (v_RTplus + m_RTplus^2) %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
      third_moment_p_x[jj] = (t_RTplus + m_RTplus^3 + 3*m_RTplus*v_RTplus) %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
      mean_RT_m_x[jj] = m_RTminus %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
      sec_moment_m_x[jj] = (v_RTminus + m_RTminus^2) %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
      third_moment_m_x[jj] = (t_RTminus + m_RTminus^3 + 3*m_RTminus*v_RTminus) %*% t(prob_a) %*% (a_supp[length(a_supp)]-a_supp[1])/99;
      
    }
  }
  
  
  if (delta == 0) {
    
    mean_RTplus = mean_RT_p_x;
    
    mean_RTminus = mean_RT_m_x;
    
    sec_mom_RTplus = sec_moment_p_x;
    
    sec_mom_RTminus = sec_moment_m_x;
    
    third_mom_RTplus = third_moment_p_x;
    
    third_mom_RTminus = third_moment_m_x;
    
  } else {
    
    mean_RTplus = mean_RT_p_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
    mean_RTminus = mean_RT_m_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
    sec_mom_RTplus = sec_moment_p_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
    sec_mom_RTminus = sec_moment_m_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
    third_mom_RTplus = third_moment_p_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
    third_mom_RTminus = third_moment_m_x %*% t(prob_x) %*% (x_set[length(x_set)]-x_set[1])/99;
    
  }
  
  
  var_RTplus = abs(sec_mom_RTplus - mean_RTplus^2);
  
  t_RTplus =  abs(third_mom_RTplus -3*var_RTplus*mean_RTplus - mean_RTplus^3);
  
  CVplus= sqrt(abs(var_RTplus))/mean_RTplus;
  
  skewplus = t_RTplus/(abs(var_RTplus))^1.5;
  
  var_RTminus = abs(sec_mom_RTminus - mean_RTminus^2);
  
  t_RTminus =  abs(third_mom_RTminus -3*var_RTminus*mean_RTminus - mean_RTminus^3);
  
  CVminus= sqrt(abs(var_RTminus))/mean_RTminus;
  
  skewminus = t_RTminus/(abs(var_RTminus))^1.5;
  
  
  # [mean_RTplus, mean_RTminus, var_RTplus,var_RTminus, t_RTplus,t_RTminus, CVplus, CVminus, skewplus, skewminus]=
  return_vals = list()
  return_vals$mean_RTplus = mean_RTplus
  return_vals$mean_RTminus = mean_RTminus
  return_vals$var_RTplus = var_RTplus
  return_vals$var_RTminus = var_RTminus
  return_vals$t_RTplus = t_RTplus
  return_vals$t_RTminus = t_RTminus
  return_vals$CVplus = CVplus
  return_vals$CVminus = CVminus
  return_vals$skewplus = skewplus
  return_vals$skewminus = skewminus
  
  return_vals
  
}





