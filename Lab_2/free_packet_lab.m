%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pacchetto di funzioni d'onda libera:                                    %
% velocita' di fase vs. velocita' gruppo                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
PhysConstants();
m = m0; %[kg]
%% Definizione problema: costruzione del pacchetto
% pacchetto centrato in k0 e relazione di dispersione:
k0=1e9;          %[m-1] (numero d'onda su cui e' centrato il pacchetto)
dk=k0/200;     %[m-1] (passo di discretizzazione spazio k)
kk=(0:dk:2*k0);   %[m-1] (asse k su cui definire il pacchetto)

% Scelta del tipo di relazione di dispersione:
% 'free', 'linear', 'sublinear'
press = 0;
while press == 0
    clc
    disp('Scegliere la relazione di dispersione:');
    disp(sprintf('\n [a] free \n [b] linear \n [c] sublinear \n'));
    type_ = input('','s');
    if (type_ == 'a' || type_ == 'b' || type_ == 'c')
        press = 1;
        switch type_
            case 'a'
                type = 'free';
            case 'b'
                type = 'linear';
            case 'c'
                type = 'sublinear';
        end
    end
end

[ww,~]=rel_disp(kk,type); % (calcolo relazione dispersione E(k))
%% Si consideri il caso di una FUNZIONE PESO g(k) GAUSSIANA [vd. lez.] e
% si focalizzi l'attenzione su tre k specifici, al fine di un confronto
alpha=1./(0.01*k0^2);         % coeff. alfa della g(k) gaussiana
gkk = exp(-alpha*(kk-k0).^2);  % distribuzione g(k)

k_l=k0-alpha^(-0.5);
k_r=k0+alpha^(-0.5);
[w0,vg,beta] =rel_disp(k0,type);
[w_l,~]=rel_disp(k_l,type);
[w_r,~]=rel_disp(k_r,type);
%% GRAFICI: relazione di dispersione e g(k)
figure(1), subplot(2,1,1)
plot(kk,ww,'b'), hold on
plot(k0,w0,'or',k_l,w_l,'ok',k_r,w_r,'og')
title('Relazione di dispersione'), xlabel('k [m^{-1}]'), ylabel('\omega [rad/s]')
subplot(2,1,2)
plot(kk,gkk,'k')
title('Funzione peso g(k)'),xlabel('k [m^{-1}]'), ylabel('g(k) [a.u.]')

%% GRAFICI nel tempo:
% Re{Psi_l(x)},Re{Psi_0(x)},Re{Psi_r(x)},Re{Psi(x)},Im{Psi(x)},|Psi(x)|^2
figure(2)
dx=0.1/k0;              %[m] (passo di discretizzazione spaziale)
x=(-0.5e3*dx:dx:2e3*dx); %[m] (definizione asse x)
dt=1/w0;              %[s] (def. passo di discretizzazione temporale)

Nt=floor(0.8*length(x)*dx/vg/dt);

for kt=1:Nt
    t=kt*dt; %[s]
    
    y=zeros(1,length(x));
    for ke=-500:500
        k=k0+ke*dk;
        [w,~]=rel_disp(k,type);
        g=exp(-alpha*(k-k0)^2);
        y=y+g*exp(i*(k*x-w*t));
    end
    

    % definizione asse spaziale per le tre componenti sotto osservazione:
    xx_l=((-0.5*pi+w_l*t)/k_l:dx:(1.5*pi+w_l*t)/k_l);
    xx_0=((-0.5*pi+w0*t)/k0:dx:(1.5*pi+w0*t)/k0);
    xx_r=((-0.5*pi+w_r*t)/k_r:dx:(1.5*pi+w_r*t)/k_r);
    % andamento spaziale delle tre componenti, fissato un certo t:
    y_l=exp(i*(k_l*xx_l-w_l*t));
    y_0=exp(i*(k0*xx_0-w0*t));
    y_r=exp(i*(k_r*xx_r-w_r*t));
    % Andamento spaziale delle tre componenti, fissato un certo t, supponendo
    % si muovano alla velocita' di gruppo vg:
    xxx=(vg*t-10*dx:dx:vg*t+10*dx);
    yy_l=exp(i*(k_l*xxx-w_l*t));
    yy_0=exp(i*(k0*xxx-w0*t));
    yy_r=exp(i*(k_r*xxx-w_r*t));
    % Grafica nell'istante t
    figure(2);
    subplot(3,1,1);plot(xx_l,real(y_l),'b',xxx,real(yy_l),'m'), axis([x(1) x(end) -2 2]), ylabel('Re\{\Psi_l(x)\}'); 
    subplot(3,1,2);plot(xx_0,real(y_0),'k',xxx,real(yy_0),'m'), axis([x(1) x(end) -2 2]), ylabel('Re\{\Psi_0(x)\}');
    subplot(3,1,3);plot(xx_r,real(y_r),'g',xxx,real(yy_r),'m'), axis([x(1) x(end) -2 2]), ylabel('Re\{\Psi_r(x)\}');
    
    figure(3);
    subplot(3,1,1);plot(x,real(y),'m'), axis([x(1) x(end) -40 40]), ylabel('Re\{\Psi(x)\}');
    subplot(3,1,2);plot(x,imag(y),'m'), axis([x(1) x(end) -40 40]), ylabel('Im\{\Psi(x)\}');
    subplot(3,1,3);plot(x,abs(y).^2,'m'), axis([x(1) x(end) 0 1400]),xlabel('x [m]'),ylabel('|\Psi(x)|^2');

    pause(0.01)
end
