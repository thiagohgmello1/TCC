%% Trabalho de conclusão de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 23/04/2020

% Rotina teste para determinação da eficiência e potência de saída em
% função da frequência para desalinhamentos verticais e horizontais.
% A compensação da drive transmissora é em paralelo.
% Todas as bobinas estão com a driver.
% Programa com implementação da expressão geral para as correntes.
%
% Parâmetros variados: desalinhamento lateral e frequência.
%
% ------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Parâmetros---------------------------------

f0=1.5e6; % Frequência de ressonância desejada
w0=2*pi*f0; % Frequência angular
f_min=0.6; % Determinação da frequência mínima analisada, em MHz
f_max=4; % Determinação da frequência máxima analisada, em MHz
delta_f=0.0005; % Determinação da precisão, em MHz
f=f_min:delta_f:f_max;
w=2*pi*f*1e6; % Frequência angular

Zinp=50;
RL=10; % Impedância da carga

v=34; % Tensão de pico
V=[v;0;0;0;0;0;0];

d=0.06;

h_min=0;
h_max=0.2;
delta_h=0.001;
h=h_min:delta_h:h_max;

% --------------Características do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magnética do vácuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade elétrica do material

% -----------------------Parâmetros construtivos--------------------------

% Desalinhamento lateral
d_12=0;
d_13=d;
d_14=d;
d_15=d_13;
d_16=d_14;
d_23=d;
d_24=d;
d_25=d_23;
d_26=d_24;
d_34=0;
d_35=2*d;
d_36=2*d;
d_45=2*d;
d_46=2*d;
d_56=0;

% Distância vertical
h_12=0.0140;
h_13=0.0692+h;
h_14=0.0567+h;
h_15=h_13;
h_16=h_14;
h_23=0.0365+h;
h_24=0.0240+h;
h_25=h_23;
h_26=h_24;
h_34=0.0120;
h_35=0;
h_36=h_34;
h_45=h_34;
h_46=0;
h_56=h_34;

% -----------------------------Drive loop--------------------------------

N1=2; % Número de espiras
b1=1.35e-3; % Raio do condutor
a1=0.25/2+2*b1; % Raio
L1=L_self(N1,a1,b1,mur);
C1=1/(w0^2*L1); % Capacitor para ressonância

% ----------------------------Transmissora--------------------------------

N2=9; % Número de espiras
b2=1.35e-3; % Raio do condutor
a2=0.25/2+2*b2; % Raio
L2=L_self(N2,a2,b2,mur);
C2=1/(w0^2*L2); % Capacitor para ressonância

% -----------------------------Receptora 1---------------------------------

N3=10; % Número de espiras
b3=0.5e-3; % Raio do condutor
a3=0.104/2+2*b3; % Raio
L3=L_self(N3,a3,b3,mur);
C3=1/(w0^2*L3); % Capacitor para ressonância

% -----------------------------Load loop 1---------------------------------

N4=2; % Número de espiras
b4=0.5e-3; % Raio do condutor
a4=0.104/2+2*b4; % Raio
L4=L_self(N4,a4,b4,mur);
C4=1/(w0^2*L4); % Capacitor para ressonância

% -----------------------------Receptora 2---------------------------------

N5=10; % Número de espiras
b5=0.5e-3; % Raio do condutor
a5=0.104/2+2*b3; % Raio
L5=L_self(N5,a5,b5,mur);
C5=1/(w0^2*L5); % Capacitor para ressonância

% -----------------------------Load loop 2---------------------------------

N6=2; % Número de espiras
b6=0.5e-3; % Raio do condutor
a6=0.104/2+2*b6; % Raio
L6=L_self(N6,a6,b6,mur);
C6=1/(w0^2*L6); % Capacitor para ressonância

%% Cálculo das resistências

delta_f=1./sqrt(pi*f*1e6*mu0*mur*sigma); % Profundidade de penetração
R1=a1*N1./(sigma*b1*delta_f);
R2=a2*N2./(sigma*b2*delta_f);
R3=a3*N3./(sigma*b3*delta_f);
R4=a4*N4./(sigma*b4*delta_f);
R5=a5*N5./(sigma*b5*delta_f);
R6=a6*N6./(sigma*b6*delta_f);

%% Determinação das correntes em cada circuito com as drives e todas as
% bobinas compensadas

Isource = zeros(length(w),length(d));
I11=zeros(length(w),length(d));
I12=zeros(length(w),length(d));
I13=zeros(length(w),length(d));
I14=zeros(length(w),length(d));
I15=zeros(length(w),length(d));
I16=zeros(length(w),length(d));

M12=MLateral(a1,a2,b1,b2,N1,N2,d_12,h_12);
M34=MLateral(a3,a4,b3,b4,N3,N4,d_34,h_34);
M35=MLateral(a3,a5,b3,b5,N3,N5,d_35,h_35);
M36=MLateral(a3,a6,b3,b6,N3,N6,d_36,h_36);
M45=MLateral(a4,a5,b4,b5,N4,N5,d_45,h_45);
M46=MLateral(a4,a6,b4,b6,N4,N6,d_46,h_46);
M56=MLateral(a5,a6,b5,b6,N5,N6,d_56,h_56);

for j=1:length(h)
    
        M13=MLateral(a1,a3,b1,b3,N1,N3,d_13,h_13(j));
        M14=MLateral(a1,a4,b1,b4,N1,N4,d_14,h_14(j));
        M15=MLateral(a1,a5,b1,b5,N1,N5,d_15,h_15(j));
        M16=MLateral(a1,a6,b1,b6,N1,N6,d_16,h_16(j));
        M23=MLateral(a2,a3,b2,b3,N2,N3,d_23,h_23(j));
        M24=MLateral(a2,a4,b2,b4,N2,N4,d_24,h_24(j));
        M25=MLateral(a2,a5,b2,b5,N2,N5,d_25,h_25(j));
        M26=MLateral(a2,a6,b2,b6,N2,N6,d_26,h_26(j));
        
    for i=1:length(w)
        
        Z=[Zinp+1/(1j*w(i)*C1) -1/(1j*w(i)*C1) 0 0 0 0 0;
            -1/(1j*w(i)*C1) 1j*w(i)*L1+1/(1j*w(i)*C1)+R1(i) 1j*w(i)*M12 1j*w(i)*M13 1j*w(i)*M14 1j*w(i)*M15 1j*w(i)*M16;...
            0 1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2) 1j*w(i)*M23 1j*w(i)*M24 1j*w(i)*M25 1j*w(i)*M26;...
            0 1j*w(i)*M13 1j*w(i)*M23 R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3) 1j*w(i)*M34 1j*w(i)*M35 1j*w(i)*M36;...
            0 1j*w(i)*M14 1j*w(i)*M24 1j*w(i)*M34 R4(i)+1/(1j*w(i)*C4)+1j*w(i)*L4+RL 1j*w(i)*M45 1j*w(i)*M46;...
            0 1j*w(i)*M15 1j*w(i)*M25 1j*w(i)*M35 1j*w(i)*M45 R5(i)+1/(1j*w(i)*C5)+1j*w(i)*L5 1j*w(i)*M56;...
            0 1j*w(i)*M16 1j*w(i)*M26 1j*w(i)*M36 1j*w(i)*M46 1j*w(i)*M56 R6(i)+1/(1j*w(i)*C6)+1j*w(i)*L6+RL];

        I=Z\V;
        
        Isource(i,j)=I(1);
        I11(i,j)=I(2);
        I12(i,j)=I(3);
        I13(i,j)=I(4);
        I14(i,j)=I(5);
        I15(i,j)=I(6);
        I16(i,j)=I(7);
        
    end
end

%%
Pout=RL*abs(I14).^2/2;
Pin=real(v*conj(Isource))/2-real(Zinp*abs(Isource).^2)/2;
eta=Pout./Pin*100;
Zin=v./I11-Zinp;

%% Gráficos de linhas

[H,F]=meshgrid(h,f);

figure(1);
contourf(H,F,eta);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Eficiência');
xlabel('d[m]');
ylabel('f[MHz]');
zlabel('\eta[%]');

figure(2);
contourf(H,F,Pout);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[W]';
title('Potência de saída');
xlabel('d[m]');
ylabel('f[MHz]');
zlabel('P_{out}[W]');

figure(3)
s=surf(H,F,eta);
hold on;
[~,i_eta] = max(eta(:));
st = scatter3(H(i_eta),F(i_eta),eta(i_eta),'r','filled');
st.SizeData = 100;
hold off;
s.EdgeColor='none';
title('Eficiência');
xlabel('d[m]');
ylabel('f[MHz]');
zlabel('\eta[%]');

figure(4)
s=surf(H,F,Pout);
hold on;
[~,i_Pout] = max(Pout(:));
st = scatter3(H(i_Pout),F(i_Pout),Pout((i_Pout)),'r','filled');
st.SizeData = 100;
hold off;
s.EdgeColor='none';
title('Potência de saída');
xlabel('d[m]');
ylabel('f[MHz]');
zlabel('P_{out}[W]');





