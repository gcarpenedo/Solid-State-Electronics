%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buca ARMONICA - WELL_4.m                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
% close all
PhysConstants();
SetGraphics();
%% DEFINIZIONE DEL PROBLEMA
L  = 20e-9;            %[m] larghezza totale dominio
%a  = 5e-9;            %[m] distanza interbuca
k  = 100000e-5;        %[Nm-1] profondità buca
k2 = 100000e0;         %
V0 = 6*q;              %[J] profondita'buca
n  = 20;               %[1] numero autovalori

dx = 5e-11;            %[m] passo discretizzazione
x  = -L/2:dx:L/2';     %[m] asse x

%% COSTRUZIONE DEL PROFILO di POTENZIALE 1D
% Scelta profilo di potenziale
profilo = 'a'; % default value
press = 0;
while press == 0
    clc
    disp('Scegliere un profilo di potenziale:');
    disp(sprintf('\n [a] BUCA QUADRATICA\n [b] BUCA QUARTA \n [c] DOPPIA BUCA PARABOLICA \n'));
    profilo = input('','s');
    if (profilo == 'a' || profilo == 'b' || profilo == 'c')
        press = 1;
    end
end

if profilo == 'a'
    % BUCA QUADRATICA
    V=.5*k*x.^2;
elseif profilo == 'b'
    % BUCA QUARTA
    V=.25*k2^3*x.^4;
else
    % DOPPIA BUCA PARABOLICA
    x0=1e-9;
    V=1e-18*((x/x0).^4-2*(x/x0).^2+2);
end
%% RISOLUZIONE
[E,psi] = es(x,V,n);
%% ...GRAFICA
% Schrodinger
figure(1);subplot(1,2,1)
plot(x*1e9,V/q,'k')
xlabel('x [nm]'); ylabel('V [eV], \Psi_{n} [a.u.]');
hold on
plot(x*1e9,(repmat(E',length(x),1)/q-psi))
axis([1e9*[-L/2 L/2] 0 max(E)/q])
% Autovalori
subplot(2,2,2)
plot(E/q,'-ok')
xlabel('n'); ylabel('E_n [eV]');
hold on
subplot(2,2,4)
loglog((0:n-1),(E-min(E))/q,'-ok')
xlabel('n-1'); ylabel('E_n - E_0 [eV]')
hold on
%% Analisi
% lsq fit su scala loglog per verificare la dipendenza di En-E0 da n-1
if (profilo == 'a' || profilo == 'b')
p = polyfit(log((1:n-1)),log((E(2:end)-min(E))/q)',1);
xx = 1:40;
figure(2);plot(xx,exp(p(1)*log(xx)+p(2)),'--m')
text(20,3e-1,sprintf('slope = %.2f',p(1)))
end


figure(3);plot(1:1:numel(E),E/q,'*r');xlabel('n'); ylabel('E_n [eV]');hold on;