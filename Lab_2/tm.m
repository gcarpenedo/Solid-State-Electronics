%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO della MATRICE di TRASFERIMENTO applicato ad una doppia barriera  %
%                                                                         %
%                                                                         %
%                                 ######      ######                 -    %
%        *         e^(i*k1*x)     #    #      #    #  T5*e^(i*k5*x)  |    %
%  E0  *   *   *   -->            #    #      #    #      -->        |    %
%           *                     #    #      #    #      *          |    %
%                                 #    #      #    #    *   *   *    |    %
%       R1*e^(-i*k1*x)   *        #    #      #    #          *      | V0 %
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
%% DEFINIZIONE PROBLEMA: profilo di potenziale a tratti
k0=1e9;   %[m-1]
V0=0.705; %[eV]  (altezza di barriera)

V1=0;       % [eV]
a1=5e-9;    % [m]
V2=V0;      % [eV]
a2=a1+2e-9; % [m]
V3=0;       % [eV]
a3=a2+1e-9; % [m]
V4=V0;      % [eV]
a4=a3+2e-9; % [m]
V5=0;       % [eV]
a5=a4+5e-9; % [m]

x = linspace(0,a5,500);
V = zeros(size(x))+V1*(x<a1)+V2*(x>=a1 & x<=a2)+V3*(x>=a2 & x<=a3)+...
     V4*(x>=a3 & x<=a4)+V5*(x>=a4 & x<=a5);
%% CALCOLO COEFFICIENTE DI TRASMISSIONE (T5) E RIFLESSIONE (R1), AL VARIARE
%  DELL' ENERGIA DELL'ONDA INCIDENTE


dE=0.001; % Definisco degli step per calcolare il cambiamento delle proprietà dell'onda per diverse energie incidenti

for kk=1:3*floor(V0/dE)

    E0(kk)=kk*dE;
    k1=(2*m*q*(E0(kk)-V1))^0.5/hb; % [m-1] %N.B Vado a calcolare i 5 k interessanti per le 5 zone create dalla successione di due buche e lo faccio per ognuna delle energie dE
    k2=(2*m*q*(E0(kk)-V2))^0.5/hb; % [m-1]
    k3=(2*m*q*(E0(kk)-V3))^0.5/hb; % [m-1]
    k4=(2*m*q*(E0(kk)-V4))^0.5/hb; % [m-1]
    k5=(2*m*q*(E0(kk)-V5))^0.5/hb; % [m-1]
    M11=matr_trasf(k1,a1);
    M21=matr_trasf(k2,a1);
    M22=matr_trasf(k2,a2);
    M32=matr_trasf(k3,a2);
    M33=matr_trasf(k3,a3);
    M43=matr_trasf(k4,a3);
    M44=matr_trasf(k4,a4);
    M54=matr_trasf(k5,a4);

   %M=inv(M54)*M44*inv(M43)*M33*inv(M32)*M22*inv(M21)*M11;
    M = (M54\M44)*(M43\M33)*(M32\M22)*(M21\M11);      %%%%%%%Ora posso trovare T5!!!!!!!
    T5(kk)=M(1,1)-M(1,2)*M(2,1)/M(2,2);
    R1(kk)=-M(2,1)/M(2,2);
end

% Grafica
figure(1);semilogy(E0,abs(T5).^2,'-b',E0,abs(R1).^2,'-r');xlabel('Energy [eV]'), ylabel('|T_{5}|^{2}, |R_{1}|^{2}');


%% CALCOLO FUNZIONI D'ONDA PER UNA DATA E0 --> Vado a calcolare praticamente tutti coefficienti di riflessione e trasmissione e poi le funzioni d'onda viaggianti
E0=0.596; %[eV]
w=E0/hb; %[rad/s]
k1=(2*m*q*(E0-V1))^0.5/hb;%[m-1]
k2=(2*m*q*(E0-V2))^0.5/hb;%[m-1]
k3=(2*m*q*(E0-V3))^0.5/hb;%[m-1]
k4=(2*m*q*(E0-V4))^0.5/hb;%[m-1]
k5=(2*m*q*(E0-V5))^0.5/hb;%[m-1]
M11=matr_trasf(k1,a1);
M21=matr_trasf(k2,a1);
M22=matr_trasf(k2,a2);
M32=matr_trasf(k3,a2);
M33=matr_trasf(k3,a3);
M43=matr_trasf(k4,a3);
M44=matr_trasf(k4,a4);
M54=matr_trasf(k5,a4);
M=(M54\M44)*(M43\M33)*(M32\M22)*(M21\M11);
%M=inv(M54)*M44*inv(M43)*M33*inv(M32)*M22*inv(M21)*M11;

T5=M(1,1)-M(1,2)*M(2,1)/M(2,2);
R1=-M(2,1)/M(2,2);
%M_45=inv(M44)*M54;
M_45=M44\M54;
T4=M_45(1,1)*T5; R4=M_45(2,1)*T5;
%M_34=inv(M33)*M43;
M_34=M33\M43;
T3=M_34(1,1)*T4+M_34(1,2)*R4; R3=M_34(2,1)*T4+M_34(2,2)*R4;
%M_23=inv(M22)*M32;
M_23= M22\M32;
T2=M_23(1,1)*T3+M_23(1,2)*R3; R2=M_23(2,1)*T3+M_23(2,2)*R3;

dx=a1/100;       % (passo di discretizzazione spaziale)
x1=(0:dx:a1)';   % (definizione asse x a tratti: 0->1)
x2=(a1:dx:a2)';  % (definizione asse x a tratti: 1->2)
x3=(a2:dx:a3)';  % (definizione asse x a tratti: 2->3)
x4=(a3:dx:a4)';  % (definizione asse x a tratti: 3->4)
x5=(a4:dx:a5)';  % (definizione asse x a tratti: 4->5)
dt=0.1/w;  %[s]   % (passo di discretizzazione temporale)

% Grafica
figure(2)
%mov = avifile ('restun.avi')
for kt=1:200
   t=kt*dt; %[s] (istante di tempo attuale)
   % calcolo funzioni d'onda
   y1=exp(i*(k1*x1-w*t))+R1*exp(-i*(k1*x1+w*t));
   y2=T2*exp(i*(k2*x2-w*t))+R2*exp(-i*(k2*x2+w*t));
   y3=T3*exp(i*(k3*x3-w*t))+R3*exp(-i*(k3*x3+w*t));
   y4=T4*exp(i*(k4*x4-w*t))+R4*exp(-i*(k4*x4+w*t));
   y5=T5*exp(i*(k5*x5-w*t));
   figure(2),subplot(3,1,1)
   plot(x,30*V,x1,real(y1),x2,real(y2),x3,real(y3),x4,real(y4),x5,real(y5))
   axis([0 a5 -30 30]), xlabel('x [m]'), ylabel('Re\{\Psi(x)\} [a.u.]')
   figure(2),subplot(3,1,2)
   plot(x,30*V,x1,imag(y1),x2,imag(y2),x3,imag(y3),x4,imag(y4),x5,imag(y5))
   axis([0 a5 -30 30]), xlabel('x [m]'), ylabel('Im\{\Psi(x)\} [a.u.]')
   figure(2),subplot(3,1,3)
   semilogy(x,1e3*V+1e-1,x1,abs(y1).^2,x2,abs(y2).^2,x3,abs(y3).^2,...
            x4,abs(y4).^2,x5,abs(y5).^2)
   axis([0 a5 1e-2 1000]), xlabel('x [m]'), ylabel('|\Psi(x)|^2 [a.u.]')
   pause(.0001)

%     F=getframe(2);
%     mov=addframe(mov,F);

end
%mov = close(mov);
%%
%movie(M)
