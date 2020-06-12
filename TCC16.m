%% Trabalho de conclus�o de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 18/05/2020

% Rotina teste para determina��o da efici�ncia e pot�ncia de sa�da em
% fun��o da frequ�ncia para desalinhamentos verticais e horizontais fixos.
% A compensa��o da transmissora � em s�rie e em paralelo.
% Todas as bobinas est�o sem a driver.
% Programa com 3 bobinas
%
% Par�metro variado: frequ�ncia.
%
% ------------------------------------------------------------------------

clear all;
clc;
close all;

% -----------------------------Par�metros---------------------------------

f0=1e6; % Frequ�ncia de resson�ncia desejada
w0=2*pi*f0; % Frequ�ncia angular
f_min=0.6; % Determina��o da frequ�ncia m�nima analisada, em MHz
f_max=5; % Determina��o da frequ�ncia m�xima analisada, em MHz
delta_f=0.0005; % Determina��o da precis�o, em MHz
f=f_min:delta_f:f_max;
w=2*pi*f*1e6; % Frequ�ncia angular

RL=10; % Imped�ncia da carga
Rs=50; % Imped�ncia da fonte
Zinp=50;

v=37.9; % Tens�o de pico
Vs=[v;0;0];

d=0.06;

h=0;

% --------------Caracter�sticas do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magn�tica do v�cuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade el�trica do material

% -----------------------Par�metros construtivos--------------------------

% Desalinhamento lateral
d_12=0;
d_13=d;
d_23=d;

% Dist�ncia vertical
h_12=h;
h_13=h;
h_23=0;

% ----------------------------- Transmissora -----------------------------

N1=2; % N�mero de espiras
b1=1.35e-3; % Raio do condutor
a1=0.25/2+2*b1; % Raio
L1=L_self(N1,a1,b1,mur);
C1=1/(w0^2*L1); % Capacitor para resson�ncia

% ---------------------------- Receptora 1--------------------------------

N2=10; % N�mero de espiras
b2=0.5e-3; % Raio do condutor
a2=0.104/2+2*b2; % Raio
L2=L_self(N2,a2,b2,mur);
C2=1/(w0^2*L2); % Capacitor para resson�ncia

% -----------------------------Receptora 2---------------------------------

N3=10; % N�mero de espiras
b3=0.5e-3; % Raio do condutor
a3=0.104/2+2*b3; % Raio
L3=L_self(N3,a3,b3,mur);
C3=1/(w0^2*L3); % Capacitor para resson�ncia

%% C�lculo das resist�ncias

delta_f=1./sqrt(pi*f*1e6*mu0*mur*sigma); % Profundidade de penetra��o
R1=a1*N1./(sigma*b1*delta_f);
R2=a2*N2./(sigma*b2*delta_f);
R3=a3*N3./(sigma*b3*delta_f);

%% Determina��o das correntes em cada circuito com as drives e todas as
% bobinas compensadas

M12=MLateral(a1,a2,b1,b2,N1,N2,d_12,h_12);
M13=MLateral(a1,a3,b1,b3,N1,N3,d_13,h_13);
M23=MLateral(a2,a3,b2,b3,N2,N3,d_23,h_23);

kappa(1) = M12/sqrt(L1*L2);
kappa(2) = M13/sqrt(L1*L3);
kappa(3) = M23/sqrt(L2*L3);

Isource = zeros(1,length(w));
I1p=zeros(1,length(w));
I2p=zeros(1,length(w));
I3p=zeros(1,length(w));

I1s=zeros(1,length(w));
I2s=zeros(1,length(w));
I3s=zeros(1,length(w));

Vp = [v;0;0;0];

for i=1:length(w)
    
    Zp=[Zinp+1/(1j*w(i)*C1) -1/(1j*w(i)*C1) 0 0;
        -1/(1j*w(i)*C1) 1j*w(i)*L1+1/(1j*w(i)*C1)+R1(i) 1j*w(i)*M12 1j*w(i)*M13;...
        0 1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2)+RL 1j*w(i)*M23;...
        0 1j*w(i)*M13 1j*w(i)*M23 R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3)+RL];
    
    Zs=[Rs+1j*w(i)*L1+1/(1j*w(i)*C1)+R1(i) 1j*w(i)*M12 1j*w(i)*M13;...
        1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2)+RL 1j*w(i)*M23;...
        1j*w(i)*M13 1j*w(i)*M23 R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3)+RL];

    Ip=Zp\Vp;
    
    Is = Zs\Vs;

    Isource(i) = Ip(1);
    I1p(i)=Ip(2);
    I2p(i)=Ip(3);
    I3p(i)=Ip(4);

    I1s(i)=Is(1);
    I2s(i)=Is(2);
    I3s(i)=Is(3);
    
end

%%
vp = v-Isource*Zinp;
Poutp=RL*abs(I2p).^2/2+RL*abs(I3p).^2/2;
Pinp=real(vp.*conj(I1p))/2;
etap=Poutp./Pinp*100;

Pouts=RL*abs(I2s).^2/2+RL*abs(I3s).^2/2;
Pins=real(v*conj(I1s))/2;
etas=Pouts./Pins*100;

%%

figure(1);
plot(f,etap);
hold on;
plot(f,etas);
hold off;
legend('Paralelo','S�rie');
title('Efici�ncia(%)');
xlabel('f[MHz]');
ylabel('\eta[%]');

figure(2);
plot(f,Poutp);
hold on;
plot(f,Pouts);
hold off;
legend('Paralelo','S�rie');
title('Pot�ncia de sa�da');
xlabel('f[MHz]');
ylabel('P_{out}[W]');

figure(3);
plot(f,abs(I1p));
hold on;
plot(f,abs(I1s));
hold off;
legend('Paralelo','S�rie');
title('Corrente na drive transmissora');
xlabel('f[MHz]');
ylabel('I[A]');

figure(4);
plot(f,abs(I3p));
hold on;
plot(f,abs(I3s));
hold off;
legend('Paralelo','S�rie');
title('Corrente na carga');
xlabel('f[MHz]');
ylabel('I[A]');
