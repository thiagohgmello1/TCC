%% Trabalho de conclusão de curso 1

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 10/11/2019

% Rotina teste para determinação da potência transferida para uma receptora
% através de uma transmissora
% Todas as bobinas estão sem a driver
% Parâmetros variados: distância e frequência
%-------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Parâmetros---------------------------------

f_r=6.78e6; % Frequência de ressonância desejada
f_min=2; % Determinação da frequência mínima analisada, em MHz
f_max=10; % Determinação da frequência máxima analisada, em MHz
delta=0.005; % Determinação da precisão, em MHz
f=f_min:delta:f_max;
w=2*pi*f*1e6; % Frequência angular
w0=2*pi*f_r;
RL=50; % Impedância da carga
Rs=50; % Impedância da fonte
v=10; % Tensão de pico
V=[v;0];
h=0:0.001:0.5;

% --------------Características do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magnética do vácuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade elétrica do material

% -----------------------Parâmetros construtivos--------------------------

% Desalinhamento lateral
d_23=0;

% Distância vertical
h_23=0.0365+h;


% ----------------------------Transmissora--------------------------------

a2=(0.25+1.35e-3)/2; % Raio
N2=9; % Número de espiras
b2=1.35e-3; % Raio do condutor
L2=L_self(N2,a2,b2,mur);
C2=1/(w0^2*L2); % Capacitor para ressonância

% -----------------------------Receptora----------------------------------

a3=(0.104+0.5e-3)/2; % Raio
N3=10; % Número de espiras
b3=0.5e-3; % Raio do condutor
L3=L_self(N3,a3,b3,mur);
C3=1/(w0^2*L3); % Capacitor para ressonância

%% Cálculo das resistências

delta=1./sqrt(pi*f*mu0*mur*sigma); % Profundidade de penetração
R2=a2*N2./(sigma*b2*delta);
R3=a3*N3./(sigma*b3*delta);

%% Determinação das correntes em cada circuito com as drives e todas as
% bobinas compensadas

for i=1:length(h)

    M23=MLateral(a2,a3,b2*N2,b3*N3,N2,N3,d_23,h_23(i));
    
    for j=1:length(f)
        
        Z=[Rs+1j*w(j)*L2+1/(1j*w(j)*C2)+R2(j) 1j*w(j)*M23;...
        1j*w(j)*M23 R3(j)+1/(1j*w(j)*C3)+1j*w(j)*L3+RL];

        I=Z\V;

        I1(i,j)=I(1);
        I2(i,j)=I(2);
    end
end
Pout=RL.*abs(I2).^2/2;
Pin=real(v*conj(I1))/2;
eta=Pout./Pin*100;
Zin=v./I1-Rs;

%%
[F,H]=meshgrid(f,h);
figure(1);
mesh(F,H,1000*Pout);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[mW]';
title('Potência de saída');
xlabel('f[MHz]');
ylabel('h[m]');
zlabel('P_{out}[mW]');

figure(2);
mesh(F,H,eta);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Eficiência');
xlabel('f[MHz]');
ylabel('h[m]');
zlabel('\eta[%]');








