%% Trabalho de conclus�o de curso 2

% Desenvolvedor: Thiago Henrique G. Mello
% Data: 23/04/2020

% Rotina teste para determina��o da efici�ncia e pot�ncia de sa�da em
% fun��o da frequ�ncia para desalinhamentos verticais e horizontais.
% A compensa��o da drive transmissora � em paralelo.
% Todas as bobinas est�o com a driver.
% Programa com implementa��o da express�o geral para as correntes.
%
% Par�metros variados: dist�ncia e frequ�ncia.
%
% ------------------------------------------------------------------------

clear all;
clc;

% -----------------------------Par�metros---------------------------------

f0=1e6; % Frequ�ncia de resson�ncia desejada
w0=2*pi*f0; % Frequ�ncia angular
f_min=0.6; % Determina��o da frequ�ncia m�nima analisada, em MHz
f_max=4; % Determina��o da frequ�ncia m�xima analisada, em MHz
delta_f=0.0005; % Determina��o da precis�o, em MHz
f=f_min:delta_f:f_max;
w=2*pi*f*1e6; % Frequ�ncia angular

Zinp=50;
RL=10; % Imped�ncia da carga

v=34; % Tens�o de pico
V=[v;0;0;0;0];

d0=0;
d=0;
h_min=0;
h_max=0.2;
delta_h=0.001;
h=h_min:delta_h:h_max;

% --------------Caracter�sticas do material e ambiente--------------------

mu0=4*pi*1e-7; % Permeabilidade magn�tica do v�cuo
mur=1; % Permeabilidade relativa
sigma=5.8e7; % Condutividade el�trica do material

% -----------------------Par�metros construtivos--------------------------

% Desalinhamento lateral
d_12=0;
d_13=d0+d;
d_14=d0+d;
d_23=d0+d;
d_24=d0+d;
d_34=0;

% Dist�ncia vertical
h_12=0.0140;
h_13=0.0692+h;
h_14=0.0567+h;
h_23=0.0365+h;
h_24=0.0240+h;
h_34=0.0120;

% -----------------------------Drive loop--------------------------------

N1=2; % N�mero de espiras
b1=1.35e-3; % Raio do condutor
a1=0.25/2+2*b1; % Raio
L1=L_self(N1,a1,b1,mur);
C1=1/(w0^2*L1); % Capacitor para resson�ncia

% ----------------------------Transmissora--------------------------------

N2=9; % N�mero de espiras
b2=1.35e-3; % Raio do condutor
a2=0.25/2+2*b2; % Raio
L2=L_self(N2,a2,b2,mur);
C2=1/(w0^2*L2); % Capacitor para resson�ncia

% -----------------------------Receptora----------------------------------


N3=10; % N�mero de espiras
b3=0.5e-3; % Raio do condutor
a3=0.104/2+2*b3; % Raio
L3=L_self(N3,a3,b3,mur);
C3=1/(w0^2*L3); % Capacitor para resson�ncia

% -------------------------------Load loop--------------------------------

N4=2; % N�mero de espiras
b4=0.5e-3; % Raio do condutor
a4=0.104/2+2*b4; % Raio
L4=L_self(N4,a4,b4,mur);
C4=1/(w0^2*L4); % Capacitor para resson�ncia

%% C�lculo das resist�ncias

delta_f=1./sqrt(pi*f*1e6*mu0*mur*sigma); % Profundidade de penetra��o
R1=a1*N1./(sigma*b1*delta_f);
R2=a2*N2./(sigma*b2*delta_f);
R3=a3*N3./(sigma*b3*delta_f);
R4=a4*N4./(sigma*b4*delta_f);

%% Determina��o das correntes em cada circuito com as drives e todas as
% bobinas compensadas

Isource = zeros(length(w),length(h));
I11=zeros(length(w),length(h));
I12=zeros(length(w),length(h));
I13=zeros(length(w),length(h));
I14=zeros(length(w),length(h));

M12=MLateral(a1,a2,b1,b2,N1,N2,d_12,h_12);
M34=MLateral(a3,a4,b3,b4,N3,N4,d_34,h_34);

for j=1:length(h)
    
        M13=MLateral(a1,a3,b1,b3,N1,N3,d_13,h_13(j));
        M14=MLateral(a1,a4,b1,b4,N1,N4,d_14,h_14(j));
        M23=MLateral(a2,a3,b2,b3,N2,N3,d_23,h_23(j));
        M24=MLateral(a2,a4,b2,b4,N2,N4,d_24,h_24(j));
        
    for i=1:length(w)
        
        Z=[Zinp+1/(1j*w(i)*C1) -1/(1j*w(i)*C1) 0 0 0;
            -1/(1j*w(i)*C1) 1j*w(i)*L1+1/(1j*w(i)*C1)+R1(i) 1j*w(i)*M12 1j*w(i)*M13 1j*w(i)*M14;...
            0 1j*w(i)*M12 R2(i)+1j*w(i)*L2+1/(1j*w(i)*C2) 1j*w(i)*M23 1j*w(i)*M24;...
            0 1j*w(i)*M13 1j*w(i)*M23 R3(i)+1j*w(i)*L3+1/(1j*w(i)*C3) 1j*w(i)*M34;...
            0 1j*w(i)*M14 1j*w(i)*M24 1j*w(i)*M34 R4(i)+1/(1j*w(i)*C4)+1j*w(i)*L4+RL];

        I=Z\V;
        
        Isource(i,j)=I(1);
        I11(i,j)=I(2);
        I12(i,j)=I(3);
        I13(i,j)=I(4);
        I14(i,j)=I(5);
        
    end
end

%%
Pout=RL*abs(I14).^2/2;
Pin=real(v*conj(Isource))/2-real(Zinp*abs(Isource).^2)/2;
eta=Pout./Pin*100;
Zin=v./I11-Zinp;

%% Gr�ficos de linhas

[H,F]=meshgrid(h,f);

figure(1);
contourf(H,F,eta);
c=colorbar('AxisLocation','out');
c.Label.String='\eta[%]';
title('Efici�ncia');
xlabel('h[m]');
ylabel('f[MHz]');
zlabel('\eta[%]');

figure(2);
contourf(H,F,Pout);
c=colorbar('AxisLocation','out');
c.Label.String='P_{out}[W]';
title('Pot�ncia de sa�da');
xlabel('h[m]');
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
title('Efici�ncia');
xlabel('h[m]');
ylabel('f[MHz]');
zlabel('\eta[%]');

%%
figure(4)
s=surf(H,F,Pout);
hold on;
[~,i_Pout] = max(Pout(:));
st = scatter3(H(i_Pout),F(i_Pout),Pout((i_Pout)),'r','filled');
st.SizeData = 100;
hold off;
s.EdgeColor='none';
title('Pot�ncia de sa�da');
xlabel('h[m]');
ylabel('f[MHz]');
zlabel('P_{out}[W]');





