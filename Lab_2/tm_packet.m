%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO della MATRICE di TRASFERIMENTO applicato ad una doppia barriera  %
% PACCHETTO D'ONDE                                                        %
%                                                                         %
%                incident         ######      ######                 -    %
%        *       wavepacket       #    #      #    #    T5*...       |    %
%  E0  *   *   *   -->            #    #      #    #      -->        |    %
%           *                     #    #      #    #      *          |    %
%                                 #    #      #    #    *   *   *    |    %
%       R1*...           *        #    #      #    #          *      | V0 %
%             <--      *   *   *  #    #      #    #                 |    %
%                            *    #    #      #    #                 |    %
%                                 #    #      #    #                 |    %
%                       ###########    ########    ##########        -    %
%                     --|---------|----|------|----|--------|---> x       %
%                       0        a1   a2     a3   a4       a5             %
%                                                                         %
%                                                                         %
%                                  potential barrier                      %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
PhysConstants();
m = m0; %[kg]
%% DEFINIZIONE PROBLEMA: costruzione del profilo di potenziale a tratti
V0= 80e-3; %[eV]  (altezza di barriera)

V1=0;         %[eV]
a1=100e-9;    %[m]
V2=V0;         %[eV]
a2=a1+5e-9;   %[m]
V3=0;         %[eV]
a3=a2+100e-9; %[m]
V4=V0;         %[eV]
a4=a3+5e-9;   %[m]
V5=0;         %[eV]
a5=a4+100e-9; %[m]
x = linspace(0,a5,500);
V = zeros(size(x))+V1*(x<a1)+V2*(x>=a1 & x<=a2)+V3*(x>=a2 & x<=a3)+...
     V4*(x>=a3 & x<=a4)+V5*(x>=a4 & x<=a5);
%% CALCOLO COEFFICIENTE DI TRASMISSIONE (T5) E RIFLESSIONE (R1), AL VARIARE
%  DELL' ENERGIA DEL PACCHETTO D'ONDE INCIDENTE con funzione peso assunta
%  gaussiana (vd. lez. ed es. prec.)
dE=.0005;                      %[eV] (passo discretizzazione energie)
E0=90e-3;                      %[eV]
w0=q*E0/hb;                    %[rad/s]
k10=(2*m*q*(E0-V1))^0.5/hb;    %[m-1] (k di picco del pacchetto incidente)
dk=k10/200;                    %[m-1] (passo discretizzazione k)
alp=1./(0.001*k10^2);           %[m2]  (nota: sigmak^2=0.5*(alp)^-1)

%0.01 in alp => sigmak = k10/sqrt(200)

% ciclo per esplorare tutte le componenti nello spazio k
for kk=1:80
    k1(kk)=k10+dk*(kk-40);           %[m-1]
    E(kk)=(hb*k1(kk))^2/2/m/q;       %[eV]
    k2(kk)=(2*m*q*(E(kk)-V2))^0.5/hb; %[m-1]
    k3(kk)=(2*m*q*(E(kk)-V3))^0.5/hb; %[m-1]
    k4(kk)=(2*m*q*(E(kk)-V4))^0.5/hb; %[m-1]
    k5(kk)=(2*m*q*(E(kk)-V5))^0.5/hb; %[m-1]
    M11=matr_trasf(k1(kk),a1);
    M21=matr_trasf(k2(kk),a1);
    M22=matr_trasf(k2(kk),a2);
    M32=matr_trasf(k3(kk),a2);
    M33=matr_trasf(k3(kk),a3);
    M43=matr_trasf(k4(kk),a3);
    M44=matr_trasf(k4(kk),a4);
    M54=matr_trasf(k5(kk),a4);

    M=(M54\M44)*(M43\M33)*(M32\M22)*(M21\M11);
%   M=inv(M54)*M44*inv(M43)*M33*inv(M32)*M22*inv(M21)*M11;
    T5(kk)=M(1,1)-M(1,2)*M(2,1)/M(2,2);
    R1(kk)=-M(2,1)/M(2,2);

    M_45=inv(M44)*M54;
    T4(kk)=M_45(1,1)*T5(kk); R4(kk)=M_45(2,1)*T5(kk);
    M_34=inv(M33)*M43;
    T3(kk)=M_34(1,1)*T4(kk)+M_34(1,2)*R4(kk);
    R3(kk)=M_34(2,1)*T4(kk)+M_34(2,2)*R4(kk);
    M_23=inv(M22)*M32;
    T2(kk)=M_23(1,1)*T3(kk)+M_23(1,2)*R3(kk);
    R2(kk)=M_23(2,1)*T3(kk)+M_23(2,2)*R3(kk);
end

%% Construzione del pacchetto nel dominio dello spazio e Grafica
%mov = avifile ('free_lab.avi')
dx=a1/100;
x1=(0:dx:a1)';
x2=(a1:dx:a2)';
x3=(a2:dx:a3)';
x4=(a3:dx:a4)';
x5=(a4:dx:a5)';
figure(2)
dt=2.5/w0; %[s]

for kt=1:250
    t=kt*dt; %[s]
    y1=zeros(length(x1),1);
    y2=zeros(length(x2),1);
    y3=zeros(length(x3),1);
    y4=zeros(length(x4),1);
    y5=zeros(length(x5),1);

    for kk=1:80
        g=exp(-alp*(k1(kk)-k10)^2);
        w=q*E(kk)/hb;                  %[rad/s]
        y1=y1+g*(exp(i*(k1(kk)*x1-w*t))+R1(kk)*exp(-i*(k1(kk)*x1+w*t)));
        y2=y2+g*(T2(kk)*exp(i*(k2(kk)*x2-w*t))+R2(kk)*exp(-i*(k2(kk)*x2+w*t)));
        y3=y3+g*(T3(kk)*exp(i*(k3(kk)*x3-w*t))+R3(kk)*exp(-i*(k3(kk)*x3+w*t)));
        y4=y4+g*(T4(kk)*exp(i*(k4(kk)*x4-w*t))+R4(kk)*exp(-i*(k4(kk)*x4+w*t)));
        y5=y5+g*(T5(kk)*exp(i*(k5(kk)*x5-w*t)));
    end
    
%    subplot(3,1,1)
%    plot(x1,real(y1),x2,real(y2),x3,real(y3),x4,real(y4),x5,real(y5))
%    axis([0 a5 -5 5])
%    subplot(3,1,2)
%    plot(x1,imag(y1),x2,imag(y2),x3,imag(y3),x4,imag(y4),x5,imag(y5))
%    axis([0 a5 -5 5])
%    subplot(3,1,3)
    plot(x1*1e9,10*abs(y1).^2,x2*1e9,10*abs(y2).^2,...
        x3*1e9,10*abs(y3).^2,x4*1e9,10*abs(y4).^2,...
        x5*1e9,10*abs(y5).^2,x*1e9,V*100000)
    axis([0 a5*1e9 1e-4 1e4])
    xlabel('x [nm]'),ylabel('|\Psi(x)|^2')
    pause(.1)

%    F=getframe(2);
%    mov=addframe(mov,F);
%    pause(.01)
end
%mov = close(mov);
%%
%movie(M)
