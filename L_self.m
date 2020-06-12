function [L] = L_self(N,r,b,mur)
% Calculate self inductance of a solenoid
% [L]=L_self(N,r,b,mur). N-Number of turns, r-solenoid ratio, b- wire
% thikness, mur- relative permeability
%------------------------------ Verified ---------------------------------
% Parâmetros
mu0=4*pi*1e-7;
D=2*r; % Diâmetro da espira 
l=2*N*b; % Comprimento do solenoide

% fator de nagaoka
kappa=D/sqrt(D^2+l^2);
[K,E]=ellipke(kappa);
kl=4*(D/l)/(3*pi)*(((2*kappa^2-1)/kappa^3)*E+((1-kappa^2)/kappa^3)*K-1);

L=pi*mu0*mur*N^2*r^2*kl/l;
end