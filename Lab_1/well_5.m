%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buca TRIANGOLARE - WELL_5.m                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
%close all
PhysConstants();
SetGraphics();
%% DEFINIZIONE DEL PROBLEMA
L  = 20e-9;            %[m] larghezza totale dominio
F  = 13e8;              %[Vm-1] campo elettrico
V0 = -3*q;             %[J] profondita'buca
n  = 30;               %[1] numero autovalori
dx = 5e-11;            %[m] passo discretizzazione
x  = -L/2:dx:L/2';     %[m] asse x
%% COSTRUZIONE DEL PROFILO di POTENZIALE 1D
V = V0 + q*F*(x+L/2);
%% RISOLUZIONE
[E,psi] = es(x,V,n);
%% ...GRAFICA
% Schrodinger
figure(1);subplot(1,2,1)
plot(x*1e9,V/q,'k')
xlabel('x [nm]'); ylabel('V [eV], \Psi_{n} [a.u.]');
hold on
plot(x*1e9,(repmat(E',length(x),1)/q-psi))
axis([1e9*[-L/2 L/2] [min(V) max(V)]/q])
% Autovalori
subplot(2,2,2)
plot((E-V0)/q,'-ok')
xlabel('n'); ylabel('E_n-V_0 [eV]');
hold on
subplot(2,2,4)
loglog((E-V0)/q,'-ok')
xlabel('n'); ylabel('E_n-V_0 [eV]')
hold on
%% Analisi
% lsq fit su scala loglog per verificare la dipendenza di En-E0 da n-1
p = polyfit(log((1:n)),log((E(1:end)-V0)/q)',1);
xx = 1:40;
figure(2);plot(xx,exp(p(1)*log(xx)+p(2)),'--m')
text(20,4e-1,sprintf('slope = %.2f',p(1)))

figure(4);plot(1:1:numel(E),(E-V0)/q,'-or'); hold on;xlabel('n'); ylabel('E_n-V_0 [eV]');