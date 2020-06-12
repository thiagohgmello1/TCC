%% Trabalho de conclusão de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 10/04/2020

% Rotina teste para determinação da inversa genérica da matriz de impedân-
% cia para um sistema composto por quatro bobinas, sendo uma drive trans-
% missora, uma main transmissora, uma main receptora e uma drive receptora.
% Além disso, determina os parâmetros listados acima para um circuito com
% três e seis bobinas, sendo que podem ou não serem levadas em consideração
% a indutância mutua entre todas as bobinas.
%
% Parâmetros de interesse: correntes em cada bobina e impedâncias de entra-
% da e saída.
%
% ----------------------------- Parâmetros -------------------------------

clear all;
%clc;

syms M12 M13 M14 M23 M24 M34 Z1 Z2 Z3 Z4 v real

V4 = [v; 0; 0; 0];
V3 = [v; 0; 0];

%% Sistma com quatro bobinas

%---------------------- Mútua entre todas as bobinas ----------------------

Z_todas = [Z1 1j*M12 1j*M13 1j*M14; 1j*M12 Z2 1j*M23 1j*M24; 1j*M13...
           1j*M23 Z3 1j*M34; 1j*M14 1j*M24 1j*M34 Z4];
Z_inv_todas = inv(Z_todas);
I_todas = Z_todas\V4;
Zin_todas = v/I_todas(1);
%%
%------------------- Mútua apenas com a main transmissora -----------------

Z_apenas_main = [Z1 1j*M12 1j*M13 1j*M14; 1j*M12 Z2 0 0; 1j*M13 0 Z3 0;...
                 1j*M14 0 0 Z4];
I_apenas_main = Z_apenas_main\V4;
Zin_apenas_main = 1/I_apenas_main(1)*v;

%------------------ Mútua apenas entre bobinas adjacentes -----------------

Z_apenas_adj = [Z1 1j*M12 0 0; 1j*M12 Z2 1j*M23 0; 0 1j*M23 Z3 1j*M34;...
                0 0 1j*M34 Z4];
I_apenas_adj = Z_apenas_adj\V4;
Zin_apenas_adj = 1/I_apenas_adj(1)*v;

%% Sistema com três bobinas

%---------------------- Mútua entre todas as bobinas ----------------------
Z_3b_todas = [Z1 1j*M12 1j*M13; 1j*M12 Z2 1j*M23; 1j*M13 1j*M23 Z3];
I_3b_todas = Z_3b_todas\V3;
Zin_3b_todas = 1/I_3b_todas(1)*v;

%------------------- Mútua apenas com a main transmissora -----------------

Z_3b_apenas_main = [Z1 1j*M12 1j*M13; 1j*M12 Z2 0; 1j*M13 0 Z3];
I_3b_apenas_main = Z_3b_apenas_main\V3;
Zin_3b_apenas_main = 1/I_3b_apenas_main(1)*v;

%------------------ Mútua apenas entre bobinas adjacentes -----------------

Z_3b_apenas_adj = [Z1 1j*M12 0; 1j*M12 Z2 1j*M23; 0 1j*M23 Z3];
I_3b_apenas_adj = Z_3b_apenas_adj\V3;
Zin_3b_apenas_adj = 1/I_3b_apenas_adj(1)*v;

%% Otimização da potência em função da mútua entre mains

Pout = Z4*I_todas(4)*conj(I_todas(4));
Pout_linha = diff(Pout,M23);
eqn = Pout_linha == 0;
S = solve(eqn,M23);






