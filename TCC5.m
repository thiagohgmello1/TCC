%% Trabalho de conclusão de curso 1

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 11/11/2019

% Rotina teste para determinação da potência transferida para uma receptora
% através de uma transmissora
% Todas as bobinas estão com a driver
%
% Parâmetros variados: desalinhamento e distância
%
%-------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Parâmetros---------------------------------

f0=6.78e6; % Frequência de ressonância desejada
w0=2*pi*f0; % Frequência angular
RL=50; % Impedância da carga
Rs=50; % Impedância da fonte
v=10; % Tensão de pico
V=[v;0;0;0];
V2=[v;0];
d=0:0.0003:0.05;
h=0:0.0005:0.5;

% --------------Características do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magnética do vácuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade elétrica do material

% -----------------------Parâmetros construtivos--------------------------

% Desalinhamento lateral
d_12=0;
d_13=d;
d_14=d;
d_23=d;
d_24=d;
d_34=0;

% Distância vertical
h_12=0.0140;
h_13=0.0692+h;
h_14=0.0567+h;
h_23=0.0365+h;
h_24=0.0240+h;
h_34=0.0120;

% -----------------------------Drive loop--------------------------------

a1=(0.25+1.35e-3)/2; % Raio
N1=1; % Número de espiras
b1=1.35e-3; % Raio do condutor
L1=L_self(N1,a1,b1,mur);
C1=1/(w0^2*L1); % Capacitor para ressonância

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

% -------------------------------Load loop--------------------------------

a4=(0.104+0.5e-3)/2; % Raio
N4=1; % Número de espiras
b4=0.5e-3; % Raio do condutor
L4=L_self(N4,a4,b4,mur);
C4=1/(w0^2*L4); % Capacitor para ressonância

%% Cálculo das resistências

delta=1/sqrt(pi*f0*mu0*mur*sigma); % Profundidade de penetração
R1=a1*N1/(sigma*b1*delta);
R2=a2*N2/(sigma*b2*delta);
R3=a3*N3/(sigma*b3*delta);
R4=a4*N4/(sigma*b4*delta);

%% Determinação das correntes em cada circuito com as drives e todas as
% bobinas compensadas

M12=MLateral(a1,a2,b1,b2,N1,N2,d_12,h_12);
M34=MLateral(a3,a4,b3,b4,N3,N4,d_34,h_34);

for i=1:length(d)
    for j=1:length(h)
        M13=MLateral(a1,a3,b1,b3,N1,N3,d_13(i),h_13(j));

        M14=MLateral(a1,a4,b1,b4,N1,N4,d_14(i),h_14(j));

        M23=MLateral(a2,a3,b2,b3,N2,N3,d_23(i),h_23(j));

        M24=MLateral(a2,a4,b2,b4,N2,N4,d_24(i),h_24(j));
        
        Z1=[Rs+1j*w0*L1+1/(1j*w0*C1)+R1 1j*w0*M12 1j*w0*M13 1j*w0*M14;...
        1j*w0*M12 R2+1j*w0*L2+1/(1j*w0*C2) 1j*w0*M23 1j*w0*M24;...
        1j*w0*M13 1j*w0*M23 R3+1j*w0*L3+1/(1j*w0*C3) 1j*w0*M34;...
        1j*w0*M14 1j*w0*M24 1j*w0*M34 R4+1/(1j*w0*C4)+1j*w0*L4+RL];

        I=Z1\V;

        I11(i,j)=I(1);
        I21(i,j)=I(2);
        I31(i,j)=I(3);
        I41(i,j)=I(4);
        
        Z2=[Rs+1j*w0*L1+R1 1j*w0*M12 1j*w0*M13 1j*w0*M14;...
        1j*w0*M12 R2+1j*w0*L2+1/(1j*w0*C2) 1j*w0*M23 1j*w0*M24;...
        1j*w0*M13 1j*w0*M23 R3+1j*w0*L3+1/(1j*w0*C3) 1j*w0*M34;...
        1j*w0*M14 1j*w0*M24 1j*w0*M34 R4+1j*w0*L4+RL];

        I=Z2\V;

        I12(i,j)=I(1);
        I22(i,j)=I(2);
        I32(i,j)=I(3);
        I42(i,j)=I(4);
        
        Z3=[Rs+1j*w0*L2+1/(1j*w0*C2)+R2 1j*w0*M23;...
        1j*w0*M23 R3+1/(1j*w0*C3)+1j*w0*L3+RL];

        I=Z3\V2;

        I13(i,j)=I(1);
        I23(i,j)=I(2);
        
    end
end
Pout1=RL.*abs(I41).^2/2;
Pin1=real(v*conj(I11))/2;
eta1=Pout1./Pin1*100;
Zin1=v./I11-Rs;

Pout2=RL.*abs(I42).^2/2;
Pin2=real(v*conj(I12))/2;
eta2=Pout2./Pin2*100;
Zin2=v./I12-Rs;

Pout3=RL.*abs(I23).^2/2;
Pin3=real(v*conj(I13))/2;
eta3=Pout3./Pin3*100;
Zin3=v./I13-Rs;

%% gráficos

% [H,D]=meshgrid(h,d);
% 
% figure(1);
% mesh(H,D,Pout1);
% c=colorbar('AxisLocation','out');
% c.Label.String='P_{out}[mW]';
% title('Potência de saída');
% xlabel('h[m]');
% ylabel('d[m]');
% zlabel('P_{out}[mW]');
% 
% figure(2);
% mesh(H,D,eta1);
% c=colorbar('AxisLocation','out');
% c.Label.String='\eta[%]';
% title('Eficiência');
% xlabel('h[m]');
% ylabel('d[m]');
% zlabel('\eta[%]');
% 
% figure(3);
% mesh(H,D,eta2);
% c=colorbar('AxisLocation','out');
% c.Label.String='\eta[%]';
% title('Eficiência');
% xlabel('h[m]');
% ylabel('d[m]');
% zlabel('\eta[%]');
% 
% figure(4);
% mesh(H,D,eta3);
% c=colorbar('AxisLocation','out');
% c.Label.String='\eta[%]';
% title('Eficiência');
% xlabel('h[m]');
% ylabel('d[m]');
% zlabel('\eta[%]');

%% Gráficos de linhas

[H,D]=meshgrid(h,d);

figure(1);
contourf(H,D,eta1);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Eficiência');
xlabel('h[m]');
ylabel('d[m]');
zlabel('\eta[%]');

figure(2);
contourf(H,D,eta2);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Eficiência');
xlabel('h[m]');
ylabel('d[m]');
zlabel('\eta[%]');

figure(3);
contourf(H,D,eta3);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Eficiência');
xlabel('h[m]');
ylabel('d[m]');
zlabel('\eta[%]');

figure(4);
contourf(H,D,Pout1);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[mW]';
title('Potência de saída');
xlabel('h[m]');
ylabel('d[m]');
zlabel('P_{out}[mW]');

figure(5);
contourf(H,D,Pout2);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[mW]';
title('Potência de saída');
xlabel('h[m]');
ylabel('d[m]');
zlabel('P_{out}[mW]');

figure(6);
contourf(H,D,Pout3);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[mW]';
title('Potência de saída');
xlabel('h[m]');
ylabel('d[m]');
zlabel('P_{out}[mW]');
