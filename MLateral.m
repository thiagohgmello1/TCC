%% Corrigida

function [M] = MLateral(rp,rs,zp,zs,Np,Ns,dc,H)
% Calculate the mutual inductance of two coils with lateral
% misalignment
% [M] = MLateral(rp,rs,hp,hs,Np,Ns,dc,H) returns de mutual
% inductance between solenoids with ratios rp and rs, hp and hs coils
% height, N1 and N2 turns, height H and misaligneds by dc meters. M is the
% mean mutual inductance.
% ----------------------------- Verified ---------------------------------

% Parâmetro
mu0=4*pi*1e-7;

M=0;
h=zeros(Np,Ns);
for i=1:Np
    for j=1:Ns
        % h(i,j)=((2*i-1)*bp/2-(2*j-1)*bs/2);
        h(i,j)=H-zp*(2*i-1)+zs*(2*j-1);
        
        % Seção necessária para cálculo através da indutância mútua média
        % rmin=sqrt(4*rp*(rs-dc)/((rs+rp-dc)^2+h(i,j)^2));
        % [K,E]=ellipke(rmin^2);
        % Gmin=(2/rmin-rmin)*K-(2/rmin)*E;
        % ML_min=mu0*rp*rs/sqrt(rp*(rs+dc))*Gmin;
        
        % Seção necessária para cálculo através da indutância mútua média e
        % da fórmula alternativa
        rmax=sqrt(4*rp*(rs+dc)/((rs+rp+dc)^2+h(i,j)^2));
        [K,E]=ellipke(rmax.^2);
        Gmax=(2/rmax-rmax)*K-(2/rmax)*E;
        % ML_max=mu0*rp*rs/sqrt(rp*(rs-dc))*Gmax;
        
        % Indutância mútua média
        % M=M+(ML_min+ML_max)/2;
        M=M+mu0*rp*rs/sqrt(rp*(rs+dc))*Gmax; % Fórmula alternativa
    end
end
end


%% Fórmula anterior

% function [M] = MLateral(rp,rs,hp,hs,Np,Ns,dc,H)
% % Calculate the mutual inductance of two coils with lateral
% % misalignment
% % [M] = MLateral(rp,rs,hp,hs,Np,Ns,dc,H) returns de mutual
% % inductance between solenoids with ratios rp and rs, hp and hs coils
% % height, N1 and N2 turns, height H and misaligneds by dc meters. M is the
% % mean mutual inductance.
% % ----------------------------- Verified ---------------------------------
% 
% % Parâmetro
% mu0=4*pi*1e-7;
% 
% zp=hp/Np; % z position 1
% zs=hs/Ns; % z position 2
% 
% M=0;
% h=zeros(Np,Ns);
% for i=1:Np
%     for j=1:Ns
%         % h(i,j)=((2*i-1)*bp/2-(2*j-1)*bs/2);
%         h(i,j)=H+zp*i-zs*j;
%         
%         % Seção necessária para cálculo através da indutância mútua média
%         % rmin=sqrt(4*rp*(rs-dc)/((rs+rp-dc)^2+h(i,j)^2));
%         % [K,E]=ellipke(rmin^2);
%         % Gmin=(2/rmin-rmin)*K-(2/rmin)*E;
%         % ML_min=mu0*rp*rs/sqrt(rp*(rs+dc))*Gmin;
%         
%         % Seção necessária para cálculo através da indutância mútua média e
%         % da fórmula alternativa
%         rmax=sqrt(4*rp*(rs+dc)/((rs+rp+dc)^2+h(i,j)^2));
%         [K,E]=ellipke(rmax.^2);
%         Gmax=(2/rmax-rmax)*K-(2/rmax)*E;
%         % ML_max=mu0*rp*rs/sqrt(rp*(rs-dc))*Gmax;
%         
%         % Indutância mútua média
%         % M=M+(ML_min+ML_max)/2;
%         M=M+mu0*rp*rs/sqrt(rp*(rs+dc))*Gmax; % Fórmula alternativa
%     end
% end
% end