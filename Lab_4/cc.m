%% CC.M [Carrier Concentration]
% Lo script calcola la concentrazione di portatori al variare della
% temperatura nelle regioni I (freeze out), II (estrinseca) e III (intrinseca)
% nota la concentrazione di drogante

clear all
close all


% Concentrazione di droganti
Nd=1e18;%[cm-3] donori
Na=1e14;%[cm-3] accettori  
dEd=-0.045;%[eV] rispetto a Ec
dEa=0.045;%[eV] rispetto a Ev

% Costanti
q=1.6e-19;%[C]
hb=6.626e-34/2/pi;%[Js]
m0=9.1e-31;%[kg]
k=1.38e-23;%[JK-1]
T0=300;%[K]
Ec_Si=1.12;%[eV]
Nc_Si=3.22e19;%[cm-3]
Nv_Si=1.83e19;%[cm-3]

% Parametri liberi
Ev=0;%[eV]
Ec=Ec_Si;%[eV]
Ed=Ec+dEd;%[eV]
Ea=Ev+dEa;%[eV]


T=(1:2000)';%[K]

ni = zeros(numel(T),1);
n = zeros(numel(T),1);
p = zeros(numel(T),1);
Ef = zeros(numel(T),1);
Nd0 = zeros(numel(T),1); 
Na0 = zeros(numel(T),1);

for ik=1:numel(T)

Nc = Nc_Si*(T(ik)/T0)^1.5;%[cm-3]
Nv = Nv_Si*(T(ik)/T0)^1.5;%[cm-3]
ni(ik)=(Nc*Nv)^0.5*exp(-q*Ec/2/k/T(ik));%[cm-3]
Nc_=Nc/2*exp(-q*(Ec-Ed)/k/T(ik));
Nv_=Nv/4*exp(q*(Ev-Ea)/k/T(ik));
if Nd>=Na
     %n-I-II
     n(ik)=Nc_/2*((1+4*(Nd-Na)/Nc_)^0.5-1);
     if n(ik)/ni(ik)< 10
          %n-II-III
          n(ik)=(Nd-Na)/2+(((Nd-Na)/2)^2+ni(ik)^2)^0.5;%[cm-3]
     end
     p(ik)=ni(ik)^2/n(ik);%[cm-3]
     Ef(ik)=Ec+k*T(ik)/q*log(n(ik)/Nc);%[eV]
     Nd0(ik)=Nd/(1+exp(q*(Ed-Ef(ik))/k/T(ik))/2);%[cm-3]
else
     %p-I-II
     p(ik)=Nv_/2*((1+4*(Na-Nd)/Nv_)^0.5-1);
     if p(ik)/ni(ik)<1
          %p-II-III
          p(ik)=(Na-Nd)/2+(((Na-Nd)/2)^2+ni(ik)^2)^0.5;%[cm-3]
     end
     n(ik)=ni(ik)^2/p(ik);%[cm-3]
     Ef(ik)=Ev-k*T(ik)/q*log(p(ik)/Nv);%[eV]
     Na0(ik)=Na/(1+exp(-q*(Ea-Ef(ik))/k/T(ik))/4);%[cm-3]
end
end

%% Grafica
figure(1)

subplot(2,1,1)
if Nd>=Na
     semilogy(T,n,'b',T,p,'r',T,ni,'g',T,Nd0,'--b')
     legend('n', 'p','n_i', 'N_{D0}')
else
     semilogy(T,n,'b',T,p,'r',T,ni,'g',T,Na0,'--r')
     legend('n', 'p','n_i', 'N_{A0}')
end
axis([0 2000 1e14 1e20])
ylabel('Concentration [cm^{-3}]'), xlabel('T [K]')
hold on

subplot(2,1,2)
plot(T,Ef)
ylabel('E_F [eV]'), xlabel('T [K]')
hold on