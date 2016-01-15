function [err,m_RT, v_RT, t_RT] =ddm_metrics(a,s,z,x0)


% This function computes unconditional decision time moments for the DDM

% Input
% a = drift rate;  s = diffusion rate;   z = symmetric threshold;        
% x0 =% initial condition


% Output 
% err = error rate

% m_RT=  mean decision time 

%v_RT = variance of decision time

% t_RT = third central moment of decision time: divide it by v_RT^1.5 to get skewness


a (abs(a) <0.01) = 0.01;


X= a*x0/s^2; Z=a*z/s^2;

X=max(-100, min(100,X));

Z=max(-100, min(100,Z));


Z(abs(Z)<0.0001)=0.0001;


err=(exp(-2*X) -exp(-2*Z))./(exp(2*Z)- exp(-2*Z));

m_RT= ((1-2*err).*z-x0)./a;


v_RT=  s^4./a.^4.*(3*Z.^2.*(csch(2*Z)).^2 -2*Z.^2.*exp(-2*X).*csch(2*Z).*coth(2*Z)...
       -4*Z.*X.*exp(-2*X).*csch(2*Z) - Z.^2.*exp(-4*X).*(csch(2*Z)).^2 + Z.*coth(2*Z)...
       -Z.*exp(-2*X).*csch(2*Z) -X );

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



