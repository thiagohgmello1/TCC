%% Trabalho de conclus�o de curso 1

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 10/11/2019

% Rotina teste para determina��o da pot�ncia transferida para uma receptora
% atrav�s de uma transmissora
% Todas as bobinas est�o com a driver
% Par�metros variados: desalinhamento e frequ�ncia
%-------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Par�metros---------------------------------

f_r=6.78e6; % Frequ�ncia de resson�ncia desejada
f_min=2; % Determina��o da frequ�ncia m�nima analisada, em MHz
f_max=10; % Determina��o da frequ�ncia m�xima analisada, em MHz
delta=0.005; % Determina��o da precis�o, em MHz
f=f_min:delta:f_max;
w=2*pi*f*1e6; % Frequ�ncia angular
w0=2*pi*f_r;
RL=50; % Imped�ncia da carga
Rs=50; % Imped�ncia da fonte
v=10; % Tens�o de pico
V=[v;0;0;0];
V2=[v;0];
d=0:0.00001:0.015;

% --------------Caracter�sticas do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magn�tica do v�cuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade el�trica do material

% -----------------------Par�metros construtivos--------------------------

% Desalinhamento lateral
d_12=0;
d_13=d;
d_14=d;
d_23=d;
d_24=d;
d_34=0;

% Dist�ncia vertical
h_12=0.0140;
h_13=0.0692;
h_14=0.0567;
h_23=0.0365;
h_24=0.0240;
h_34=0.0120;

% -----------------------------Drive loop--------------------------------

a1=(0.25+1.35e-3)/2; % Raio
N1=2; % N�mero de espiras
b1=1.35e-3; % Raio do condutor
L1=L_self(N1,a1,b1,mur);
C1=1/(w0^2*L1); % Capacitor para resson�ncia

% ----------------------------Transmissora--------------------------------

a2=(0.25+1.35e-3)/2; % Raio
N2=9; % N�mero de espiras
b2=1.35e-3; % Raio do condutor
L2=L_self(N2,a2,b2,mur);
C2=1/(w0^2*L2); % Capacitor para resson�ncia

% -----------------------------Receptora----------------------------------

a3=(0.104+0.5e-3)/2; % Raio
N3=10; % N�mero de espiras
b3=0.5e-3; % Raio do condutor
L3=L_self(N3,a3,b3,mur);
C3=1/(w0^2*L3); % Capacitor para resson�ncia

% -------------------------------Load loop--------------------------------

a4=(0.104+0.5e-3)/2; % Raio
N4=2; % N�mero de espiras
b4=0.5e-3; % Raio do condutor
L4=L_self(N4,a4,b4,mur);
C4=1/(w0^2*L4); % Capacitor para resson�ncia

%% C�lculo das resist�ncias

delta=1./sqrt(pi*f*mu0*mur*sigma); % Profundidade de penetra��o
R1=a1*N1./(sigma*b1*delta);
R2=a2*N2./(sigma*b2*delta);
R3=a3*N3./(sigma*b3*delta);
R4=a4*N4./(sigma*b4*delta);

%% Determina��o das correntes em cada circuito com as drives e todas as
% bobinas compensadas

M12=MLateral(a1,a2,b1*N1,b2*N2,N1,N2,d_12,h_12);
M34=MLateral(a3,a4,b3*N3,b4*N4,N3,N4,d_34,h_34);

for i=1:length(d)
    
    kappa(i,1)=M12/sqrt(L1*L2);
    kappa(i,6)=M34/sqrt(L3*L4);

    M13=MLateral(a1,a3,b1*N1,b3*N3,N1,N3,d_13(i),h_13);
    kappa(i,2)=M13/sqrt(L1*L3);

    M14=MLateral(a1,a4,b1*N1,b4*N4,N1,N4,d_14(i),h_14);
    kappa(i,3)=M14/sqrt(L1*L4);

    M23=MLateral(a2,a3,b2*N2,b3*N3,N2,N3,d_23(i),h_23);
    kappa(i,4)=M23/sqrt(L2*L3);

    M24=MLateral(a2,a4,b2*N2,b4*N4,N2,N4,d_24(i),h_24);
    kappa(i,5)=M24/sqrt(L2*L4);
    
    for j=1:length(f)
        Z1=[Rs+1j*w(j)*L1+1/(1j*w(j)*C1)+R1(j) 1j*w(j)*M12 1j*w(j)*M13 1j*w(j)*M14;...
        1j*w(j)*M12 R2(j)+1j*w(j)*L2+1/(1j*w(j)*C2) 1j*w(j)*M23 1j*w(j)*M24;...
        1j*w(j)*M13 1j*w(j)*M23 R3(j)+1j*w(j)*L3+1/(1j*w(j)*C3) 1j*w(j)*M34;...
        1j*w(j)*M14 1j*w(j)*M24 1j*w(j)*M34 R4(j)+1/(1j*w(j)*C4)+1j*w(j)*L4+RL];

        I=Z1\V;

        I11(i,j)=I(1);
        I21(i,j)=I(2);
        I31(i,j)=I(3);
        I41(i,j)=I(4);
        
        Z2=[Rs+1j*w(j)*L1+R1(j) 1j*w(j)*M12 1j*w(j)*M13 1j*w(j)*M14;...
        1j*w(j)*M12 R2(j)+1j*w(j)*L2+1/(1j*w(j)*C2) 1j*w(j)*M23 1j*w(j)*M24;...
        1j*w(j)*M13 1j*w(j)*M23 R3(j)+1j*w(j)*L3+1/(1j*w(j)*C3) 1j*w(j)*M34;...
        1j*w(j)*M14 1j*w(j)*M24 1j*w(j)*M34 R4(j)+1j*w(j)*L4+RL];

        I=Z2\V;

        I12(i,j)=I(1);
        I22(i,j)=I(2);
        I32(i,j)=I(3);
        I42(i,j)=I(4);
        
        Z3=[Rs+1j*w(j)*L2+1/(1j*w(j)*C2)+R2(j) 1j*w(j)*M23;...
        1j*w(j)*M23 R3(j)+1/(1j*w(j)*C3)+1j*w(j)*L3+RL];

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

%%

[F,D]=meshgrid(f,d);

figure(1);
mesh(F,D,Pout1);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[mW]';
title('Pot�ncia de sa�da');
xlabel('f[MHz]');
ylabel('d[m]');
zlabel('P_{out}[mW]');

figure(2);
mesh(F,D,eta1);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Efici�ncia');
xlabel('f[MHz]');
ylabel('d[m]');
zlabel('\eta[%]');

figure(3);
mesh(F,D,eta2);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Efici�ncia');
xlabel('f[MHz]');
ylabel('d[m]');
zlabel('\eta[%]');

figure(4);
mesh(F,D,eta3);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Efici�ncia');
xlabel('f[MHz]');
ylabel('d[m]');
zlabel('\eta[%]');






