%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pacchetto di funzioni d'onda libera + Barriera di Potenziale (E0 <V0)   %
%                                                                         %
%                                 ######################            -     %
%        *                        #                    #            |     %
%  E0  *   *   *   -->            #                    #            | V0  %
%           *                     #                    #            |     %
%                     #############         a          ###########  -     %
%                                  <------------------->                  %
%      wave packet                   potential barrier                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
PhysConstants();
m=m0;           %[Kg]
%% Definizione problema: profilo di potenziale e pacchetto
E0=55e-3;       %[eV] (energia del pacchetto)
V0=60e-3;       %[eV] (altezza della barriera)
a=4e-9;          %[m] (larghezza della barriera)

k0=(2*m*q*E0)^0.5/hb;                 %[m-1]
k0_=k0*(2*q*V0*m/hb^2/k0^2-1)^.5;     %[m-1]
w0=hb*k0^2/2/m;                       %[rad/s]
dk=k0/200;
alp=1./(0.01*k0^2);

dx=0.05/k0;                           %[m]
dx_nm = dx*1e9;                       %[nm]
Ndx_lr = 2.5e3;                       %[1]
x_l=(-Ndx_lr*dx:dx:0);                %[m]
x_0=(0:dx:a);                         %[m]
x_r=(a:dx:Ndx_lr*dx);                 %[m]
dt=0.5/w0;                            %[s]
x = [x_l x_0 x_r];
V = zeros(size(x)) + V0*(x>=0 & x<=a)*1e4; % [a.u.]


for kt=-90:90
 t=kt*dt;                             %[s]
 y_l=zeros(1,length(x_l));
 y_0=zeros(1,length(x_0));
 y_r=zeros(1,length(x_r));
    
    for ke=-40:40
     k=k0+ke*dk;
     k_=k*(2*q*V0*m/hb^2/k^2-1)^0.5; %[m-1]
     w=hb*k^2/2/m;                     %[rad/s]

     B=-i*(k^2+k_^2)/k_/k*sinh(k_*a)/(2*cosh(k_*a)+i*sinh(k_*a)*(k_^2-k^2)/k_/k);
     E=2*exp(-i*k*a)/(2*cosh(k_*a)+i*sinh(k_*a)*(k_^2-k^2)/k_/k);
     D=E*exp(-k_*a)*exp(i*k*a)/2*(1+i*k/k_);
     C=E*exp(k_*a)*exp(i*k*a)/2*(1-i*k/k_);

     g=exp(-alp*(k-k0)^2);
     y_l=y_l+g*(exp(i*(k*x_l))+B*exp(i*(-k*x_l)))*exp(-i*w*t);
     y_0=y_0+g*(C*exp(-k_*x_0)+D*exp(k_*x_0))*exp(-i*w*t);
     y_r=y_r+g*E*exp(i*(k*x_r))*exp(-i*w*t);
    end

 subplot(3,1,1)
 plot(x_l*1e9,real(y_l),x_0*1e9,real(y_0),x_r*1e9,real(y_r),x*1e9,V)
 axis([[-Ndx_lr Ndx_lr]*dx_nm -100 100])
 ylabel('Re\{\Psi(x)\}'), xlabel('x [nm]')
 subplot(3,1,2)
 plot(x_l*1e9,imag(y_l),x_0*1e9,imag(y_0),x_r*1e9,imag(y_r),x*1e9,V)
 axis([[-Ndx_lr Ndx_lr]*dx_nm -100 100])
 ylabel('Im\{\Psi(x)\}'), xlabel('x [nm]')
 subplot(3,1,3)
 plot(x_l*1e9,abs(y_l).^2,x_0*1e9,abs(y_0).^2,x_r*1e9,abs(y_r).^2,x*1e9,V*10)
 axis([[-Ndx_lr Ndx_lr]*dx_nm 0 2000])
 ylabel('|\Psi(x)|^2'), xlabel('x [nm]')

pause(0.1)

end
