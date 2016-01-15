# ddm_metrics
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



ddm_metrics <- function( a, s, z, x0 ) {
  a [abs(a) < 0.01] = 0.01;
  
  X = a*x0/s^2; Z=a*z/s^2;  # Create normalized versions of starting point and threshold (k_x = X and k_z = Z)
  
  # Unlike Matlab, R will collapse X and Z below down to one element:
  X = pmax(-100, pmin(100,X));
  
  Z = pmax(-100, pmin(100,Z));
  
  Z[abs(Z)<0.0001]=0.0001;  # Make sure that these values are within ranges that don't cause numerical problems
  
  err=(exp(-2*X) -exp(-2*Z))/(exp(2*Z)- exp(-2*Z));
  
  m_RT= ((1-2*err)*z-x0)/a;
  
  
  v_RT=  s^4/a^4*(3*Z^2*(csch(2*Z))^2 -2*Z^2*exp(-2*X)*csch(2*Z)*coth(2*Z) -
                     4*Z*X*exp(-2*X)*csch(2*Z) - Z^2*exp(-4*X)*(csch(2*Z))^2 + Z*coth(2*Z)-
                     Z*exp(-2*X)*csch(2*Z) -X );
  
  t_RT = abs(s^6/a^6/8*( (48*X^2*Z + 48*X*Z - 32*Z^3 + 12*Z)*exp(-2*X)*(csch(2*Z))^3  - 
                           (24*X^2*Z + 24*X*Z^2 + 24*X*Z + 8*Z^3 + 12*Z^2 + 6*Z)*exp(4*Z - 2*X)*(csch(2*Z))^3 + 
                           (48*X*Z^2 + 12*Z^2 - 24*Z^3)*exp(-2*Z - 4*X)*(csch(2*Z))^3 - 
                           (48*X*Z^2 + 12*Z^2 + 24*Z^3)*exp(2*Z - 4*X)*(csch(2*Z))^3 - 
                           16*Z^3*exp(-6*X)*(csch(2*Z))^3 - 6*Z*cosh(2*Z)*(csch(2*Z))^3 + 
                           6*Z*cosh(6*Z)*(csch(2*Z))^3 + 18*X*sinh(2*Z)*(csch(2*Z))^3 - 
                           6*X*sinh(6*Z)*(csch(2*Z))^3 + 112*Z^3*cosh(2*Z)*(csch(2*Z))^3 + 72*Z^2*sinh(2*Z)*(csch(2*Z))^3 -
                           (6*Z - 12*Z^2 + 8*Z^3+ 24*X*Z - 24*X*Z^2 + 24*X^2*Z)*exp(-4*Z -2*X)*(csch(2*Z))^3));
  
  return_list = list()
  
  return_list$err = err
  return_list$m_RT = m_RT
  return_list$v_RT = v_RT
  return_list$t_RT = t_RT
  
  return_list
  
}