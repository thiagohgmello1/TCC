%% 
% Algorítmo para plotar a característica da capacitância em função do
% coeficiente de acoplamento para as quatro topologias de compensação

clear all;
clc;

Q=5;
R1=0.1;
R2=0.1;
Rs=50;
L1=1e-6;
L2=1e-6;
k=0:0.00001:0.8;
M=k*sqrt(L1*L2);
fo=1e7;
wo=2*pi*fo;
C2=1/(wo^2*L2);
RL=wo*L1/Q;
Vs=10;

C1SS1=1/(wo^2*L1);
C1SP1=((RL*R2+wo^2*L2^2)^2+(wo*L2*R2)^2)./(wo^2*(L1*((RL*R2+wo^2*L2^2)^2+(wo*L2*R2)^2)-wo^4*M.^2*L2^3));
C1PS1=L1*(R2+RL)^2./(wo^2*L1^2*(R2+RL)^2+((R2+RL)*R1+wo^2*M.^2).^2);
C1PP1=(wo*L1-wo^2*M.^2*(wo*L2)^3/((RL*R2+wo^2*L2^2)^2+(wo*L2*R2)^2))./(wo*((R1+(wo^2*M.^2*(RL^2*R2+(wo*L2)^2*(RL+R2)))/((RL*R2+wo^2*L2^2)^2+(wo*L2*R2)^2)).^2+(wo*L1-wo^2*M.^2*(wo*L2)^3/((RL*R2+wo^2*L2^2)^2+(wo*L2*R2)^2))).^2);

C1SP2=((wo^2*L2^2)^2)./(wo^2*((L1*(wo^2*L2^2)^2)-wo^4*M.^2*L2^3));
C1PS2=L1*(RL)^2./(wo^2*L1^2*(RL)^2+(wo^2*M.^2).^2);
C1PP2=(wo*L1-wo^2*M.^2*(wo*L2)^3/((wo^2*L2^2)^2))./(wo*(((wo^2*M.^2*((wo*L2)^2*(RL)))/((wo^2*L2^2)^2)).^2+(wo*L1-wo^2*M.^2*(wo*L2)^3/((wo^2*L2^2)^2))).^2);

eSS=(k./k-k./k);
eSP=(C1SP1-C1SP2)./C1SP1*100;
ePS=(C1PS1-C1PS2)./C1PS1*100;
ePP=(C1PP1-C1PP2)./C1PP1*100;

figure(1);
plot(k,abs(eSS));
hold on;
plot(k,abs(eSP));
hold on;
plot(k,abs(ePS));
hold on;
plot(k,abs(ePP));
title('Erro(%)');
xlabel('\kappa');
ylabel('e(%)');
legend('SS','SP','PS','PP','Location','Southwest');

figure(2);
plot(k,k./k);
hold on;
plot(k,C1SP1/C1SS1);
hold on;
plot(k,C1PS1/C1SS1);
hold on;
plot(k,C1PP1/C1SS1);
title('C_1 para várias topologias de compensação');
legend('C_{SS}','C_{SP}','C_{PS}','C_{PP}','Location','northwest');
xlabel('Coeficiente de acoplamento magnético \kappa');
ylabel('C_1 normalizado');

figure(3);
plot(k,k./k);
hold on;
plot(k,C1SP2/C1SS1);
hold on;
plot(k,C1PS2/C1SS1);
hold on;
plot(k,C1PP2/C1SS1);
title('C_1 para várias topologias de compensação');
legend('C_{SS}','C_{SP}','C_{PS}','C_{PP}','Location','northwest');
xlabel('Coeficiente de acoplamento magnético \kappa');
ylabel('C_1 normalizado');

%%

RL=0:0.001:10000;
M=0.1e-6;
Pout=RL*wo^2*M^2*Vs^2./(2*R1^2*(R2+RL+wo^2*M^2/R1).^2);
plot(RL,Pout);

%%
C = (L1*L2-M.^2)*C2*L2^2./(M.^4*C2*RL/L2+(L1*L2-M.^2).^2);

%%
plot(k,C/C1SS1);
hold on;
plot(k,C1PP1/C1SS1);
hold off;
legend('Artigo','Calculado')

