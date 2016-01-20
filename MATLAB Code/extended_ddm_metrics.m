
function [error, mean_RT, var_RT, t_RT, CV, skew]= extended_ddm_metrics(a_mean, sigma_a, s, x0_mean, delta, z)

% Copyright 2016, Vaibhav Srivastava

% Input: 
% a_mean =mean drift rate, sigma_a = variance in drift rate,
% s=diffusion rate,
% x0_mean= mean initial condition, 
%delta = support for initial condition (x0-delta, x0+delta), 
%z= threshold

%Output
% error = error rate
% mean_RT = mean decision time
% var_RT = variance of decision time
% t_RT = third central moment of decision time 
% CV = coefficient of variation of decision time 
% skew = skewness of decision time 


if delta>z
    delta=z;
end

if delta==0
    x_set =x0_mean;
    prob_x=1;
else
    
    delta=max(0.0001,delta);
    tmp_x = linspace(-delta,delta,100);
    x_set = x0_mean + tmp_x;
    prob_x = 1/2/delta*ones(size(x_set));
end



mean_RT_x= zeros(size(x_set));
sec_moment_x= zeros(size(x_set));
third_moment_x= zeros(size(x_set));
err_x=zeros(size(x_set));

for jj=1:length(x_set)
    
    x0=x_set(jj);
    
    if sigma_a ==0
        a_supp=a_mean;
        
        [err,m_RT, v_RT, t_RT] =ddm_metrics(a_supp,s,z,x0);
        
        err_x(jj) = err;
        
        mean_RT_x(jj) = m_RT;
        
        sec_moment_x(jj) = (v_RT + m_RT.^2);
        
        third_moment_x(jj) = (t_RT + m_RT.^3 + 3*m_RT.*v_RT);
        
    else
        
        sigma_a=max(0.0001,sigma_a);
        
        tmp_a=linspace(-5*sigma_a, 5*sigma_a,100);
        
        a_supp= a_mean + tmp_a ;
        
        a_supp(abs(a_supp)<0.0001)=0.0001;
        
        prob_a = exp(max(-tmp_a.^2/2/sigma_a^2, -100))./(sqrt(2*pi)*sigma_a);
        
        
        [err,m_RT, v_RT, t_RT] =ddm_metrics(a_supp,s,z,x0);
        
        err_x(jj) = err*prob_a'*(a_supp(end)-a_supp(1))/99;
        
        mean_RT_x(jj) = m_RT*prob_a'*(a_supp(end)-a_supp(1))/99;
        
        sec_moment_x(jj) = (v_RT + m_RT.^2)*prob_a'*(a_supp(end)-a_supp(1))/99;
        
        third_moment_x(jj) = (t_RT + m_RT.^3 + 3*m_RT.*v_RT)*prob_a'*(a_supp(end)-a_supp(1))/99;
    end
end

if delta == 0
    
    error= err_x;
    
    mean_RT= mean_RT_x;
    
    sec_mom_RT= sec_moment_x;
    
    third_mom_RT= third_moment_x;
    
    
else
    error= err_x*prob_x'*(x_set(end)-x_set(1))/99;
    
    mean_RT= mean_RT_x*prob_x'*(x_set(end)-x_set(1))/99;
    
    sec_mom_RT= sec_moment_x*prob_x'*(x_set(end)-x_set(1))/99;
    
    third_mom_RT= third_moment_x*prob_x'*(x_set(end)-x_set(1))/99;
    
end


var_RT = abs(sec_mom_RT - mean_RT.^2);

t_RT =  abs(third_mom_RT -3*var_RT*mean_RT - mean_RT^3);

CV= sqrt(abs(var_RT))/mean_RT;

skew = t_RT/(abs(var_RT))^1.5;




function [err,m_RT, v_RT, t_RT] =ddm_metrics(a,s,z,x0)

a (abs(a) <0.01) = 0.01;

X= a*x0/s^2; Z=a*z/s^2;

X=max(-200, min(200,X));


Z=max(-200, min(200,Z));

Z(abs(Z)<0.0001)=0.0001;


%k=2*a/s^2;

err=(exp(-2*X) -exp(-2*Z))./(exp(2*Z)- exp(-2*Z));

m_RT= abs(((1-2*err).*z-x0)./a);


v_RT=  abs( s^4./a.^4.*(3*Z.^2.*(csch(2*Z)).^2 -2*Z.^2.*exp(-2*X).*csch(2*Z).*coth(2*Z)...
    -4*Z.*X.*exp(-2*X).*csch(2*Z) - Z.^2.*exp(-4*X).*(csch(2*Z)).^2 + Z.*coth(2*Z)...
    -Z.*exp(-2*X).*csch(2*Z) -X ));

% 3*z.^2/a^2./(sinh(k*z))^2 + 2*z/k/a^2./tanh(k*z)-2*x0/k/a^2 ...
%      -2*z.^2.*exp(-k*x0)/a^2./sinh(k*z)./tanh(k*z)-2*z.*exp(-k*x0)/k/a^2./sinh(k*z) ...
%      - 4*z.*x0.*exp(-k*x0)/a^2./sinh(k*z) -z.^2*exp(-2*k*x0)/a^2./(sinh(k*z)).^2;

t_RT = abs(s^6./a.^6/8.*( ...
    (48*X.^2.*Z + 48*X.*Z - 32*Z.^3 + 12*Z).*exp(-2*X).*(csch(2*Z)).^3 ...
    - (24*X.^2.*Z + 24*X.*Z.^2 + 24*X.*Z + 8*Z.^3 + 12*Z.^2 + 6*Z).*exp(4*Z - 2*X).*(csch(2*Z)).^3 ...
    + (48*X.*Z.^2 + 12*Z.^2 - 24*Z.^3).*exp(-2*Z - 4*X).*(csch(2*Z)).^3 ...
    - (48*X.*Z.^2 + 12*Z.^2 + 24*Z.^3).*exp(2*Z - 4*X).*(csch(2*Z)).^3 ...
    - 16*Z.^3.*exp(-6*X).*(csch(2*Z)).^3 - 6*Z.*cosh(2*Z).*(csch(2*Z)).^3  ...
    + 6*Z.*cosh(6*Z).*(csch(2*Z)).^3 + 18*X.*sinh(2*Z).*(csch(2*Z)).^3 ...
    - 6*X.*sinh(6*Z).*(csch(2*Z)).^3 + 112*Z.^3.*cosh(2*Z).*(csch(2*Z)).^3 + 72*Z.^2.*sinh(2*Z).*(csch(2*Z)).^3  ...
     -(6*Z - 12*Z.^2 + 8*Z.^3+ 24*X.*Z - 24*X.*Z.^2 + 24*X.^2.*Z).*exp(-4*Z -2*X).*(csch(2*Z)).^3));
 


% 
% t_RT = abs( s^6./a.^6.*(csch(2*Z)).^3/8.*( ...
%     (48*X.^2.*Z + 48*X.*Z - 32*Z.^3 + 12*Z).*exp(-2*X) ...
%     - (24*X.^2.*Z + 24*X.*Z.^2 + 24*X.*Z + 8*Z.^3 + 12*Z.^2 + 6*Z).*exp(4*Z - 2*X) ...
%     + (48*X.*Z.^2 + 12*Z.^2 - 24*Z.^3).*exp(-2*Z - 4*X) ...
%     - (48*X.*Z.^2 + 12*Z.^2 + 24*Z.^3).*exp(2*Z - 4*X) ...
%     - 16*Z.^3.*exp(-6*X) - 6*Z.*cosh(2*Z) + 6*Z.*cosh(6*Z) + 18*X.*sinh(2*Z)...
%     - 6*X.*sinh(6*Z) + 112*Z.^3.*cosh(2*Z) + 72*Z.^2.*sinh(2*Z)  ...
%     -(6*Z - 12*Z.^2 + 8*Z.^3+ 24*X.*Z - 24*X.*Z.^2 + 24*X.^2.*Z).*exp(-4*Z -2*X)));





