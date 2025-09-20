%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buca QUADRATA TRIPLA - WELL_3.m                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
%close all
PhysConstants();
SetGraphics();
%% DEFINIZIONE DEL PROBLEMA
a  = 1e-9;             %[m] larghezza della buca
L  = 20e-9;            %[m] larghezza totale dominio
d  = 1e-9;             %[m] distanza interbuca
V0 = 5*q;              %[J] profondita' buca
n  = 18;               %[1] numero autovalori
dx = 5e-11;            %[m] passo discretizzazione
x  = -L/2:dx:L/2';     %[m] asse x
%% COSTRUZIONE DEL PROFILO di POTENZIALE 1D
V = zeros(size(x));
V = V - V0*((x>=-d/2-a-d & x<=-d-a/2) | ...
            (x>=-a/2 & x<=a/2) | (x>=a/2+d & x<=a/2+d+a));
%% RISOLUZIONE
[E,psi] = es(x,V,n);
%% ...GRAFICA
% Schrodinger
subplot(1,2,1)
plot(x*1e9,V/q,'k')
xlabel('x [nm]'); ylabel('V [eV], \Psi_{n} [a.u.]');
hold on
plot(x*1e9,(repmat(E',length(x),1)/q-psi))
axis([1e9*[-L/2 L/2] [min(V)*1.2 min(V)+(max(V)-min(V))*1.2]/q])
% Autovalori
subplot(2,2,2)
plot(E/q,'-ok')
xlabel('n'); ylabel('E_n [eV]');
hold on
subplot(2,2,4)
plot((E+V0)/q,'-ok')
xlabel('n'); ylabel('E_n+V_0 [eV]')
hold on