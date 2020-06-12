%% Trabalho de conclusão de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 14/04/2020

% Rotina teste para determinação do ponto ótimo no qual tem-se a menor
% indutância mútua entre bobinas main e drives.
%
% Parâmetros de interesse: indutância mútua e desalinhamento lateral.
%
% ------------------------------------------------------------------------

clear all;
clc;

syms a b delta d
rmax = (4*a*(b+delta)/((a+b+delta)^2+d^2))^(1/2);
[K,E] = ellipke(rmax);
G = (2/rmax-rmax)*K-2/rmax*E;
ML = a*b/(a*(b+delta))^(1/2)*G;
ML_linha = diff(ML,delta);
eqn = ML_linha == 0;
S = solve(eqn,delta);






