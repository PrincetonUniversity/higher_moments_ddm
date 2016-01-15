# ddm_metrics_cond
#
# Written by Patrick Simen, Jan 3, 2016
# This code is a transcription into R of the Matlab code written by Vaibhav Srivastava, 
# accompanying the publication:
# "Explicit moments of decision times for single- and double-threshold drift-diffusion processes"
# Submitted to Journal of Mathematical Psychology


# It will help to define a few trig and hyperbolic functions that are not built in to R:
csch <- function( x ){
  1/sinh(x);
}

coth <- function( x ){
  cosh(x)/sinh(x);
}



#function [m_RTplus, m_RTminus, v_RTplus, v_RTminus, t_RTplus, t_RTminus] =ddm_metrics_cond(a,s,z,x0)
ddm_metrics_cond <- function (a,s,z,x0) {
  
  
  # This function computes conditional decision time moments for the DDM
  
  # Input
  # a = drift rate;  s = diffusion rate;   z = symmetric threshold;        
  # x0 =% initial condition
  
  
  # Output 
  
  # m_RTplus=  mean decision time conditioned on correct decision
  
  # m_RTminus=  mean decision time conditioned on error
  
  # v_RTplus=  variance of decision time conditioned on correct decision
  
  # v_RTminus=  variance of decision time conditioned on error
  
  # t_RTplus = third central moment of decision time conditioned on correct decision: divide it by v_RTplus^1.5 to get skewness
  
  # t_RTplus = third central moment of decision time conditioned on error: divide it by v_RTminus^1.5 to get skewness
  
  
  
  
  
  a [abs(a) <0.01] = 0.01;
  
  
  X= a*x0/s^2; Z=a*z/s^2;   # Create normalized versions of starting point and threshold (k_x = X and k_z = Z)
  
  X=pmax(-100, pmin(100,X));
  
  Z=pmax(-100, pmin(100,Z));
  
  
  Z[abs(Z)<0.0001]=0.0001;
  
  
  m_RTplus= s^2/(a^2)*(2*Z*coth(2*Z) - (X+Z)*coth(X+Z));
  
  m_RTminus= s^2/(a^2)*(2*Z*coth(2*Z) - (-X+Z)*coth(-X+Z));
  
  
  v_RTplus= s^4/(a^4)*(4*Z^2*(csch(2*Z))^2 + 2*Z*coth(2*Z) - (Z+X)^2*(csch(Z+X))^2 - (Z+X)*coth(Z+X));
  
  v_RTminus= s^4/(a^4)*(4*Z^2*(csch(2*Z))^2 + 2*Z*coth(2*Z) - (Z-X)^2*(csch(Z-X))^2 - (Z-X)*coth(Z-X));
  
  
  t_RTplus= s^6/(a^6)*(12*Z^2*(csch(2*Z))^2 + 16*Z^3*coth(2*Z)*(csch(2*Z))^2 + 6*Z*coth(2*Z)
                       -3*(Z+X)^2*(csch(Z+X))^2 - 2*(Z+X)^3*coth(Z+X)*(csch(Z+X))^2 - 3*(Z+X)*coth(Z+X));
  
  
  t_RTminus= s^6/(a^6)*(12*Z^2*(csch(2*Z))^2 + 16*Z^3*coth(2*Z)*(csch(2*Z))^2 + 6*Z*coth(2*Z)
                        -3*(Z-X)^2*(csch(Z-X))^2 - 2*(Z-X)^3*coth(Z-X)*(csch(Z-X))^2 - 3*(Z-X)*coth(Z-X));
  
  
  
  #[m_RTplus, m_RTminus, v_RTplus, v_RTminus, t_RTplus, t_RTminus]
  return_list = list()
  return_list$m_RTplus = m_RTplus;
  return_list$m_RTminus = m_RTminus;
  return_list$v_RTplus = v_RTplus;
  return_list$v_RTminus = v_RTminus; 
  return_list$t_RTplus = t_RTplus;
  return_list$t_RTminus = t_RTminus;
  
  return_list
}