%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO della MATRICE di TRASFERIMENTO applicato ad un reticolo di       %
% sei buche di potenziale e sette barriere                                %
%                                 (1)    ...   (7)                        %
%                incident        ######   \\  ######                 -    %
%        *       wavepacket      #    #   //  #    #    T15*...      |    %
%  E0  *   *   *   -->           #    #   \\  #    #      -->        |    %
%           *                    #    #   //  #    #      *          |    %
%                                #    #   \\  #    #    *   *   *    |    %
%       R1*...          *        #    #   //  #    #          *      | Vn %
%             <--     *   *   *  #    #   \\  #    #                 |    %
%                           *    #    #   //  #    #                 |    %
%                                #    #   \\  #    #                 |    %
%                      ###########    ##  //  #    #############     -    %
%                    --|---------|----|-- \\ -|----|--------|---> x       %
%                       0        a1   a2  // a13  a14      a15            %
%                                                                         %
%                                                                         %
%                                  potential barrier                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

PhysConstants();
m = m0; %[kg]
%% DEFINIZIONE PROBLEMA: profilo di potenziale a tratti (15 sezioni)
V0= 0.1;%[eV]

V1=0;   a1=10e-9;      % [eV]; [m]
V2=V0;  a2=a1+2e-9;    % [eV]; [m]
V3=0;   a3=a2+3e-9;    % [eV]; [m]
V4=V0;  a4=a3+2e-9;    % [eV]; [m]
V5=0;   a5=a4+3e-9;    % [eV]; [m]
V6=V0;  a6=a5+2e-9;    % [eV]; [m]
V7=0;   a7=a6+3e-9;    % [eV]; [m]
V8=V0;  a8=a7+2e-9;    % [eV]; [m]
V9=0;   a9=a8+3e-9;    % [eV]; [m]
V10=V0; a10=a9+2e-9;   % [eV]; [m]
V11=0;  a11=a10+3e-9;  % [eV]; [m]
V12=V0; a12=a11+2e-9;  % [eV]; [m]
V13=0;  a13=a12+3e-9;  % [eV]; [m]
V14=V0; a14=a13+2e-9;  % [eV]; [m]
V15=0;  a15=a14+10e-9; % [eV]; [m]

x = linspace(0,a15,5000);
V = zeros(size(x))+V1*(x<a1)+V2*(x>=a1 & x<a2)+V3*(x>=a2 & x<a3)+...
     V4*(x>=a3 & x<a4)+V5*(x>=a4 & x<a5)+V6*(x>=a5 & x<a6)+ ...
     V7*(x>=a6 & x<a7)+V8*(x>=a7 & x<a8)+V9*(x>=a8 & x<a9)+ ...
     V10*(x>=a9 & x<a10)+V11*(x>=a10 & x<a11)+V12*(x>=a11 & x<a12)+ ...
     V13*(x>=a12 & x<a13)+V14*(x>=a13 & x<a14)+V15*(x>=a14 & x<=a15);
 
%% CALCOLO COEFFICIENTE DI TRASMISSIONE (T15) E RIFLESSIONE (R1), AL VARIARE
%  DELL' ENERGIA DELL'ONDA INCIDENTE

dE = 4e-5; %[eV]
E0_vect = dE:dE:V0; %[eV]
T15 = zeros(numel(E0_vect),1);
R1 = zeros(numel(E0_vect),1);

for kk=1:numel(E0_vect)
    E0 = E0_vect(kk);

    % % definizione k sui vari tratti per ogni E0
    
    k1=(2*m*q*(E0-V1))^0.5/hb;   %[m-1]
    k2=(2*m*q*(E0-V2))^0.5/hb;   %[m-1]
    k3=(2*m*q*(E0-V3))^0.5/hb;   %[m-1]
    k4=(2*m*q*(E0-V4))^0.5/hb;   %[m-1]
    k5=(2*m*q*(E0-V5))^0.5/hb;   %[m-1]
    k6=(2*m*q*(E0-V6))^0.5/hb;   %[m-1]
    k7=(2*m*q*(E0-V7))^0.5/hb;   %[m-1]
    k8=(2*m*q*(E0-V8))^0.5/hb;   %[m-1]
    k9=(2*m*q*(E0-V9))^0.5/hb;   %[m-1]
    k10=(2*m*q*(E0-V10))^0.5/hb; %[m-1]
    k11=(2*m*q*(E0-V11))^0.5/hb; %[m-1]
    k12=(2*m*q*(E0-V12))^0.5/hb; %[m-1]
    k13=(2*m*q*(E0-V13))^0.5/hb; %[m-1]
    k14=(2*m*q*(E0-V14))^0.5/hb; %[m-1]
    k15=(2*m*q*(E0-V15))^0.5/hb; %[m-1]

    M1_1=matr_trasf(k1,a1);
    M2_1=matr_trasf(k2,a1);
    M2_2=matr_trasf(k2,a2);
    M3_2=matr_trasf(k3,a2);
    M3_3=matr_trasf(k3,a3);
    M4_3=matr_trasf(k4,a3);
    M4_4=matr_trasf(k4,a4);
    M5_4=matr_trasf(k5,a4);
    M5_5=matr_trasf(k5,a5);
    M6_5=matr_trasf(k6,a5);
    M6_6=matr_trasf(k6,a6);
    M7_6=matr_trasf(k7,a6);
    M7_7=matr_trasf(k7,a7);
    M8_7=matr_trasf(k8,a7);
    M8_8=matr_trasf(k8,a8);
    M9_8=matr_trasf(k9,a8);
    M9_9=matr_trasf(k9,a9);
    M10_9=matr_trasf(k10,a9);
    M10_10=matr_trasf(k10,a10);
    M11_10=matr_trasf(k11,a10);
    M11_11=matr_trasf(k11,a11);
    M12_11=matr_trasf(k12,a11);
    M12_12=matr_trasf(k12,a12);
    M13_12=matr_trasf(k13,a12);
    M13_13=matr_trasf(k13,a13);
    M14_13=matr_trasf(k14,a13);
    M14_14=matr_trasf(k14,a14);
    M15_14=matr_trasf(k15,a14);

    % calcolo matrice di trasferimento globale 
    
    M=(M5_4\M4_4)*(M4_3\M3_3)*(M3_2\M2_2)*(M2_1\M1_1);
    M=(M9_8\M8_8)*(M8_7\M7_7)*(M7_6\M6_6)*(M6_5\M5_5)*M;
    M=(M12_11\M11_11)*(M11_10\M10_10)*(M10_9\M9_9)*M;
    M=(M15_14\M14_14)*(M14_13\M13_13)*(M13_12\M12_12)*M;

    % calcolo coefficiente T15 e R1
    
    T15(kk)=M(1,1)-M(1,2)*M(2,1)/M(2,2);
    R1(kk)=-M(2,1)/M(2,2);
end

figure(1);semilogy(E0_vect,abs(T15).^2,'-r',E0_vect,abs(R1).^2,'-b');xlabel('Energy [eV]');ylabel('|T|^2, |R_1|^2');
figure(3);plot(E0_vect,abs(T15).^2,'-r',E0_vect,abs(R1).^2,'-b');xlabel('Energy [eV]');ylabel('|T|^2, |R_1|^2');
figure(4);plot(E0_vect,abs(T15).^2,'-r');xlabel('Energy [eV]');ylabel('|T|^2, |R_1|^2');
%% CALCOLO FUNZIONI D'ONDA PER UNA DATA E0
E0 = 0.077;  %[eV]
w  = E0/hb;  %[rad/s]

% definizione k sui vari tratti:
k1=(2*m*q*(E0-V1))^0.5/hb;   %[m-1]
k2=(2*m*q*(E0-V2))^0.5/hb;   %[m-1]
k3=(2*m*q*(E0-V3))^0.5/hb;   %[m-1]
k4=(2*m*q*(E0-V4))^0.5/hb;   %[m-1]
k5=(2*m*q*(E0-V5))^0.5/hb;   %[m-1]
k6=(2*m*q*(E0-V6))^0.5/hb;   %[m-1]
k7=(2*m*q*(E0-V7))^0.5/hb;   %[m-1]
k8=(2*m*q*(E0-V8))^0.5/hb;   %[m-1]
k9=(2*m*q*(E0-V9))^0.5/hb;   %[m-1]
k10=(2*m*q*(E0-V10))^0.5/hb; %[m-1]
k11=(2*m*q*(E0-V11))^0.5/hb; %[m-1]
k12=(2*m*q*(E0-V12))^0.5/hb; %[m-1]
k13=(2*m*q*(E0-V13))^0.5/hb; %[m-1]
k14=(2*m*q*(E0-V14))^0.5/hb; %[m-1]
k15=(2*m*q*(E0-V15))^0.5/hb; %[m-1]

% definizione delle matrici di trasferimento sui vari tratti:
M1_1=matr_trasf(k1,a1);
M2_1=matr_trasf(k2,a1);
M2_2=matr_trasf(k2,a2);
M3_2=matr_trasf(k3,a2);
M3_3=matr_trasf(k3,a3);
M4_3=matr_trasf(k4,a3);
M4_4=matr_trasf(k4,a4);
M5_4=matr_trasf(k5,a4);
M5_5=matr_trasf(k5,a5);
M6_5=matr_trasf(k6,a5);
M6_6=matr_trasf(k6,a6);
M7_6=matr_trasf(k7,a6);
M7_7=matr_trasf(k7,a7);
M8_7=matr_trasf(k8,a7);
M8_8=matr_trasf(k8,a8);
M9_8=matr_trasf(k9,a8);
M9_9=matr_trasf(k9,a9);
M10_9=matr_trasf(k10,a9);
M10_10=matr_trasf(k10,a10);
M11_10=matr_trasf(k11,a10);
M11_11=matr_trasf(k11,a11);
M12_11=matr_trasf(k12,a11);
M12_12=matr_trasf(k12,a12);
M13_12=matr_trasf(k13,a12);
M13_13=matr_trasf(k13,a13);
M14_13=matr_trasf(k14,a13);
M14_14=matr_trasf(k14,a14);
M15_14=matr_trasf(k15,a14);

% calcolo matrice di trasferimento globale
M=(M5_4\M4_4)*(M4_3\M3_3)*(M3_2\M2_2)*(M2_1\M1_1);  % M=inv(M5_4)*M4_4*inv(M4_3)*M3_3*inv(M3_2)*M2_2*inv(M2_1)*M1_1;
M=(M9_8\M8_8)*(M8_7\M7_7)*(M7_6\M6_6)*(M6_5\M5_5)*M; % M=inv(M9_8)*M8_8*inv(M8_7)*M7_7*inv(M7_6)*M6_6*inv(M6_5)*M5_5*M;
M=(M12_11\M11_11)*(M11_10\M10_10)*(M10_9\M9_9)*M;   % M=inv(M12_11)*M11_11*inv(M11_10)*M10_10*inv(M10_9)*M9_9*M;
M=(M15_14\M14_14)*(M14_13\M13_13)*(M13_12\M12_12)*M; % M=inv(M15_14)*M14_14*inv(M14_13)*M13_13*inv(M13_12)*M12_12*M;

T15=M(1,1)-M(1,2)*M(2,1)/M(2,2);
R1=-M(2,1)/M(2,2);

% calcolo coefficienti T e R in ogni regione di potenziale

%%%%%
M_14_15 = M14_14\M15_14;
T14 = M_14_15(1,1)*T15; 
R14 = M_14_15(2,1)*T15;
%%%%%%
M_13_14 = M13_13\M14_13;
T13 = M_13_14(1,1)*T14 + M_13_14(1,2)*R14;
R13 = M_13_14(2,1)*T14 + M_13_14(2,2)*R14;
%%%%%%
M_12_13 = M12_12\M13_12;
T12 = M_12_13(1,1)*T13 + M_12_13(1,2)*R13;
R12 = M_12_13(2,1)*T13 + M_12_13(2,2)*R13;
%%%%%%
M_11_12 = M11_11\M12_11;
T11 = M_11_12(1,1)*T12 + M_11_12(1,2)*R12;
R11 = M_11_12(2,1)*T12 + M_11_12(2,2)*R12;
%%%%%%
M_10_11 = M10_10\M11_10;
T10 = M_10_11(1,1)*T11 + M_10_11(1,2)*R11;
R10 = M_10_11(2,1)*T11 + M_10_11(2,2)*R11;
%%%%%%
M_9_10 = M9_9\M10_9;
T9 = M_9_10(1,1)*T10 + M_9_10(1,2)*R10;
R9 = M_9_10(2,1)*T10 + M_9_10(2,2)*R10;
%%%%%%
M_8_9 = M8_8\M9_8;
T8 = M_8_9(1,1)*T9 + M_8_9(1,2)*R9;
R8 = M_8_9(2,1)*T9 + M_8_9(2,2)*R9;
%%%%%%
M_7_8 = M7_7\M8_7;
T7 = M_7_8(1,1)*T8 + M_7_8(1,2)*R8;
R7 = M_7_8(2,1)*T8 + M_7_8(2,2)*R8;
%%%%%%
M_6_7 = M6_6\M7_6;
T6 = M_6_7(1,1)*T7 + M_6_7(1,2)*R7;
R6 = M_6_7(2,1)*T7 + M_6_7(2,2)*R7;
%%%%%%
M_5_6 = M5_5\M6_5;
T5 = M_5_6(1,1)*T6 + M_5_6(1,2)*R6;
R5 = M_5_6(2,1)*T6 + M_5_6(2,2)*R6;
%%%%%%
M_4_5 = M4_4\M5_4;
T4 = M_4_5(1,1)*T5 + M_4_5(1,2)*R5;
R4 = M_4_5(2,1)*T5 + M_4_5(2,2)*R5;
%%%%%%
M_3_4 = M3_3\M4_3;
T3 = M_3_4(1,1)*T4 + M_3_4(1,2)*R4;
R3 = M_3_4(2,1)*T4 + M_3_4(2,2)*R4;
%%%%%%
M_2_3 = M2_2\M3_2;
T2 = M_2_3(1,1)*T3 + M_2_3(1,2)*R3;
R2 = M_2_3(2,1)*T3 + M_2_3(2,2)*R3;

% Grafica
dx=a1/100;
x1=(0:dx:a1)';
x2=(a1:dx:a2)';
x3=(a2:dx:a3)';
x4=(a3:dx:a4)';
x5=(a4:dx:a5)';
x6=(a5:dx:a6)';
x7=(a6:dx:a7)';
x8=(a7:dx:a8)';
x9=(a8:dx:a9)';
x10=(a9:dx:a10)';
x11=(a10:dx:a11)';
x12=(a11:dx:a12)';
x13=(a12:dx:a13)';
x14=(a13:dx:a14)';
x15=(a14:dx:a15)';

dt=0.1/w; %[s]

for kt=1:100
    t=kt*dt; %[s]
    
    y1=exp(i*(k1*x1-w*t))+R1*exp(-i*(k1*x1+w*t));
    y2=T2*exp(i*(k2*x2-w*t))+R2*exp(-i*(k2*x2+w*t));
    y3=T3*exp(i*(k3*x3-w*t))+R3*exp(-i*(k3*x3+w*t));
    y4=T4*exp(i*(k4*x4-w*t))+R4*exp(-i*(k4*x4+w*t));
    y5=T5*exp(i*(k5*x5-w*t))+R5*exp(-i*(k5*x5+w*t));
    y6=T6*exp(i*(k6*x6-w*t))+R6*exp(-i*(k6*x6+w*t));
    y7=T7*exp(i*(k7*x7-w*t))+R7*exp(-i*(k7*x7+w*t));
    y8=T8*exp(i*(k8*x8-w*t))+R8*exp(-i*(k8*x8+w*t));
    y9=T9*exp(i*(k9*x9-w*t))+R9*exp(-i*(k9*x9+w*t));
    y10=T10*exp(i*(k10*x10-w*t))+R10*exp(-i*(k10*x10+w*t));
    y11=T11*exp(i*(k11*x11-w*t))+R11*exp(-i*(k11*x11+w*t));
    y12=T12*exp(i*(k12*x12-w*t))+R12*exp(-i*(k12*x12+w*t));
    y13=T13*exp(i*(k13*x13-w*t))+R13*exp(-i*(k13*x13+w*t));
    y14=T14*exp(i*(k14*x14-w*t))+R14*exp(-i*(k14*x14+w*t));
    y15=T15*exp(i*(k15*x15-w*t));
    
    figure(2), subplot(3,1,1)
    plot(x1,real(y1),x2,real(y2),x3,real(y3),x4,real(y4),x5,real(y5),x6,real(y6),x7,real(y7),x8,real(y8),x9,real(y9),x10,real(y10),x11,real(y11),x12,real(y12),x13,real(y13),x14,real(y14),x15,real(y15))
    axis([0 a15 -6 6]), xlabel('x [m]'), ylabel('Re\{\Psi(x)\} [a.u.]')
    figure(2), subplot(3,1,2)
    plot(x1,imag(y1),x2,imag(y2),x3,imag(y3),x4,imag(y4),x5,imag(y5),x6,imag(y6),x7,imag(y7),x8,imag(y8),x9,imag(y9),x10,imag(y10),x11,imag(y11),x12,imag(y12),x13,imag(y13),x14,imag(y14),x15,imag(y15))
    axis([0 a15 -6 6]), xlabel('x [m]'), ylabel('Im\{\Psi(x)\} [a.u.]')
    figure(2), subplot(3,1,3)
    plot(x1,abs(y1).^2,x2,abs(y2).^2,x3,abs(y3).^2,x4,abs(y4).^2,x5,abs(y5).^2,x6,abs(y6).^2,x7,abs(y7).^2,x8,abs(y8).^2,x9,abs(y9).^2,x10,abs(y10).^2,x11,abs(y11).^2,x12,abs(y12).^2,x13,abs(y13).^2,x14,abs(y14).^2,x15,abs(y15).^2)
    axis([0 a15 0 32]), xlabel('x [m]'), ylabel('|\Psi(x)|^2 [a.u.]')
    
%     
%     pause(0.01)
end
