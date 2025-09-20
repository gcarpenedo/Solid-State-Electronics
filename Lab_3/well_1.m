%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buca QUADRATA SINGOLA - WELL_1.m                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
PhysConstants();
SetGraphics();
%% DEFINIZIONE DEL PROBLEMA
a  = 3e-9;             %[m] larghezza della buca
L  = 10e-9;            %[m] larghezza totale dominio
V0 = 0.1*q;            %[J] profondita' buca
n  = 3;                %[1] numero autovalori
dx = 1e-11;            %[m] passo discretizzazione
x  = -L/2:dx:L/2';     %[m] asse x
%% COSTRUZIONE DEL PROFILO di POTENZIALE 1D
V = zeros(size(x));
V = V - V0*(x>=-a/2 & x<=a/2);
%% RISOLUZIONE
[E,psi] = es(x,V,n);
%%  ...GRAFICA
% Schrodinger
subplot(1,2,1)
plot(x*1e9,V/q,'k')
xlabel('x [nm]'); ylabel('V [eV], \Psi_{n} [a.u.]');
hold on
plot(x*1e9,(repmat(E',length(x),1)/q-psi))
axis([1e9*[-L/2 L/2] [min(V)*1.2 min(V)+(max(V)-min(V))*1.2]/q])
% Autovalori
subplot(2,2,2)
plot(E/q+V0/q,'-ok')
xlabel('n'); ylabel('E_n+V_0 [eV]');
hold on
subplot(2,2,4)
plot(((E+V0)/q).^0.5,'-ok')
xlabel('n'); ylabel('(E_n+V_0)^{1/2} [eV^{1/2}]')
hold on