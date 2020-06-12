%% Trabalho de conclus�o de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 10/04/2020

% Rotina teste para determina��o da inversa gen�rica da matriz de imped�n-
% cia para um sistema composto por quatro bobinas, sendo uma drive trans-
% missora, uma main transmissora, uma main receptora e uma drive receptora.
% Al�m disso, determina os par�metros listados acima para um circuito com
% tr�s e seis bobinas, sendo que podem ou n�o serem levadas em considera��o
% a indut�ncia mutua entre todas as bobinas.
%
% Par�metros de interesse: correntes em cada bobina e imped�ncias de entra-
% da e sa�da.
%
% ----------------------------- Par�metros -------------------------------

clear all;
clc;

syms M12 M13 M14 M15 M16 M23 M24 M25 M26 M34 M35 M36 M45 M46 M56 Z1 Z2...
    Z3 Z4 Z5 Z6 v real

V = [v; 0; 0; 0; 0; 0];

%% Sistma com seis bobinas

%---------------------- M�tua entre todas as bobinas ----------------------

Z_todas = [Z1 1j*M12 1j*M13 1j*M14 1j*M15 1j*M16;...
           1j*M12 Z2 1j*M23 1j*M24 1j*M25 1j*M26;...
           1j*M13 1j*M23 Z3 1j*M34 1j*M35 1j*M36;...
           1j*M14 1j*M24 1j*M34 Z4 1j*M45 1j*M46;...
           1j*M15 1j*M25 1j*M35 1j*M45 Z5 1j*M56;...
           1j*M16 1j*M26 1j*M36 1j*M46 1j*M56 Z6];

I_todas = Z_todas\V;
Zin_todas = v/I_todas(1);




