%% Trabalho de conclusão de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 21/04/2020

% Rotina teste para determinação da eficiência e potência de saída em
% funçãoda frequência para desalinhamentos verticais e horizontais fixos.
% A compensação da drive transmissora é em paralelo.
% Todas as bobinas estão com a driver.
% Programa com implementação de transformação de fonte na entrada.
%
% Parâmetro variado: frequência.
%
% ------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Parâmetros---------------------------------

f0=1e6; % Frequência de ressonância desejada
w0=2*pi*f0; % Frequência angular
f_min=0.6; % Determinação da frequência mínima analisada, em MHz
f_max=1.3; % Determinação da frequência máxima analisada, em MHz
delta_f=0.0005; % Determinação da precisão, em MHz
f=f_min:delta_f:f_max;
w=2*pi*f*1e6; % Frequência angular

RL=50; % Impedância da carga
Rs=50; % Impedância da fonte

v=10; % Tensão de pico
Vs=[v;0;0;0];

d=0.06;

h=0.1;

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
N1=2; % Número de espiras
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
N4=2; % Número de espiras
b4=0.5e-3; % Raio do condutor
L4=L_self(N4,a4,b4,mur);
C4=1/(w0^2*L4); % Capacitor para ressonância

%% Cálculo das resistências

delta_f=1./sqrt(pi*f*1e6*mu0*mur*sigma); % Profundidade de penetração
R1=a1*N1./(sigma*b1*delta_f);
R2=a2*N2./(sigma*b2*delta_f);
R3=a3*N3./(sigma*b3*delta_f);
R4=a4*N4./(sigma*b4*delta_f);

%% Determinação das correntes em cada circuito com as drives e todas as
% bobinas compensadas

M12=MLateral(a1,a2,b1,b2,N1,N2,d_12,h_12);
M34=MLateral(a3,a4,b3,b4,N3,N4,d_34,h_34);
M13=MLateral(a1,a3,b1,b3,N1,N3,d_13,h_13);
M14=MLateral(a1,a4,b1,b4,N1,N4,d_14,h_14);
M23=MLateral(a2,a3,b2,b3,N2,N3,d_23,h_23);
M24=MLateral(a2,a4,b2,b4,N2,N4,d_24,h_24);

I11p=zeros(1,length(w));
I12p=zeros(1,length(w));
I13p=zeros(1,length(w));
I14p=zeros(1,length(w));

I11s=zeros(1,length(w));
I12s=zeros(1,length(w));
I13s=zeros(1,length(w));
I14s=zeros(1,length(w));

vp = v./(1+1j*w*C1*Rs);
Zeqp = Rs./(1j*w*C1*Rs+1);

for i=1:length(w)
    
    Vp = [vp(i);0;0;0];
    Zp=[Zeqp(i)+1j*w(i)*L1+R1(i) 1j*w(i)*M12 1j*w(i)*M13 ...
        1j*w(i)*M14;1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2)...
        1j*w(i)*M23 1j*w(i)*M24;1j*w(i)*M13 1j*w(i)*M23...
        R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3) 1j*w(i)*M34;1j*w(i)*M14...
        1j*w(i)*M24 1j*w(i)*M34 R4(i)+1/(1j*w(i)*C4)+1j*w(i)*L4+RL];
    
    Zs=[Rs+1j*w(i)*L1+1/(1j*w(i)*C1)+R1(i) 1j*w(i)*M12 1j*w(i)*M13 ...
        1j*w(i)*M14;1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2)...
        1j*w(i)*M23 1j*w(i)*M24;1j*w(i)*M13 1j*w(i)*M23...
        R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3) 1j*w(i)*M34;1j*w(i)*M14...
        1j*w(i)*M24 1j*w(i)*M34 R4(i)+1/(1j*w(i)*C4)+1j*w(i)*L4+RL];

    Ip=Zp\Vp;
    
    Is = Zs\Vs;

    I11p(i)=Ip(1);
    I12p(i)=Ip(2);
    I13p(i)=Ip(3);
    I14p(i)=Ip(4);

    I11s(i)=Is(1);
    I12s(i)=Is(2);
    I13s(i)=Is(3);
    I14s(i)=Is(4);
    
end
vin = vp-I11p.*Zeqp;
Poutp=RL*abs(I14p).^2/2;
Pinp=real(vin.*conj(I11p))/2;
etap=Poutp./Pinp*100;

Pouts=RL*abs(I14s).^2/2;
Pins=real(v*conj(I11s))/2;
etas=Pouts./Pins*100;

%%

figure(1);
plot(f,etap);
hold on;
plot(f,etas);
hold off;
legend('Paralelo','Série');
title('Eficiência(%)');
xlabel('f[MHz]');
ylabel('\eta[%]');

figure(2);
plot(f,Poutp);
hold on;
plot(f,Pouts);
hold off;
legend('Paralelo','Série');
title('Potência de saída');
xlabel('f[MHz]');
ylabel('P_{out}[W]');


