function [m_RTplus, m_RTminus, v_RTplus, v_RTminus, t_RTplus, t_RTminus] =ddm_metrics_cond(a,s,z,x0)


% This function computes conditional decision time moments for the DDM

% Input
% a = drift rate;  s = diffusion rate;   z = symmetric threshold;        
% x0 =% initial condition


% Output 

% m_RTplus=  mean decision time conditioned on correct decision

% m_RTminus=  mean decision time conditioned on error

% v_RTplus=  variance of decision time conditioned on correct decision

% v_RTminus=  variance of decision time conditioned on error

% t_RTplus = third central moment of decision time conditioned on correct decision: divide it by v_RTplus^1.5 to get skewness

% t_RTplus = third central moment of decision time conditioned on error: divide it by v_RTminus^1.5 to get skewness

a (abs(a) <0.01) = 0.01;


X= a*x0/s^2; Z=a*z/s^2;

X=max(-100, min(100,X));

Z=max(-100, min(100,Z));


Z(abs(Z)<0.0001)=0.0001;


m_RTplus= s^2./(a.^2).*(2*Z.*coth(2*Z) - (X+Z).*coth(X+Z));

m_RTminus= s^2./(a.^2).*(2*Z.*coth(2*Z) - (-X+Z).*coth(-X+Z));


v_RTplus= s^4./(a.^4).*(4*Z.^2.*(csch(2*Z)).^2 + 2*Z.*coth(2*Z) - (Z+X).^2.*(csch(Z+X)).^2 - (Z+X).*coth(Z+X));

v_RTminus= s^4./(a.^4).*(4*Z.^2.*(csch(2*Z)).^2 + 2*Z.*coth(2*Z) - (Z-X).^2.*(csch(Z-X)).^2 - (Z-X).*coth(Z-X));


t_RTplus= s^6./(a.^6).*(12*Z.^2.*(csch(2*Z)).^2 + 16*Z.^3.*coth(2*Z).*(csch(2*Z)).^2 + 6*Z.*coth(2*Z)...
    -3*(Z+X).^2.*(csch(Z+X)).^2 - 2*(Z+X).^3.*coth(Z+X).*(csch(Z+X)).^2 - 3*(Z+X).*coth(Z+X));

        
t_RTminus= s^6./(a.^6).*(12*Z.^2.*(csch(2*Z)).^2 + 16*Z.^3.*coth(2*Z).*(csch(2*Z)).^2 + 6*Z.*coth(2*Z)...
    -3*(Z-X).^2.*(csch(Z-X)).^2 - 2*(Z-X).^3.*coth(Z-X).*(csch(Z-X)).^2 - 3*(Z-X).*coth(Z-X));
