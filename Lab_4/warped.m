% WARPED.M
% Riproduce la superficie di Fermi in banda di valenza per Si e Ge

clear all
close all

%costanti
q=1.6e-19;%[C]
hb=6.626e-34/2/pi;%[Js]
m0=9.1e-31;%[kg]
k=1.38e-23;%[JK-1]
Ec_Si=1.12;%[eV]
ml_Si=0.9163*m0;%[kg]
mt_Si=0.1905*m0;%[kg]
mhh_Si=0.537*m0;%[kg]
mlh_Si=0.153*m0;%[kg]
Nc_Si=3.22d19;%[cm-3]
Nv_Si=1.83d19;%[cm-3]
Ec_Ge=0.66;%[eV]
ml_Ge=1.588*m0;%[kg]
mt_Ge=0.08152*m0;%[kg]
mhh_Ge=0.347*m0;%[kg]
mlh_Ge=0.0429*m0;%[kg]
Ec_GaAs=1.42;%[eV]
ml_GaAs=0.067*m0;%[kg]
mt_GaAs=0.067*m0;%[kg]
mhh_GaAs=0.51*m0;%[kg]
mlh_GaAs=0.082*m0;%[kg]
Nc_GaAs=4.21e17;%[cm-3]
Nv_GaAs=9.52e18;%[cm-3]
a_Si=5.43e-10;%[m] lattice constant
a_Ge=5.646e-10;%[m] lattice constant
a_GaAs=5.6533e-10;%[m] lattice constant

%parametri Ge (4K)
A_Ge=-13*hb^2.d0/2/m0;%[Jm2]
B_Ge=8.9*hb^2.d0/2/m0;%[Jm2]
C_Ge=10.3*hb^2/2/m0;%[Jm2]

%parametri Si (4K)
A_Si=-4.1*hb^2/2/m0;%[Jm2]
B_Si=1.6*hb^2/2/m0;%[Jm2]
C_Si=3.3*hb^2/2/m0;%[Jm2]

%parametri liberi
Ev=0;%[eV]
Ec=Ec_Si;%[eV]
T=300;%[K]
Ef=Ec/2;%[eV]
a=a_Si;%[m]
b=2*pi/a;%[m-1]
A=A_Si;
B=B_Si;
C=C_Si;

kx=(-b:b/100:b)';
ky=(-b:b/100:b)';

aa=size(kx);Na=aa(1);

for ikx=1:Na
for iky=1:Na
E_lh(ikx,iky)=A*(kx(ikx)^2+ky(iky)^2)+(B^2*(kx(ikx)^2+ky(iky)^2)^2+C^2*kx(ikx)^2*ky(iky)^2)^0.5;
E_hh(ikx,iky)=A*(kx(ikx)^2+ky(iky)^2)-(B^2*(kx(ikx)^2+ky(iky)^2)^2+C^2*kx(ikx)^2*ky(iky)^2)^0.5;
end
end

figure
subplot(2,1,1)
surf(kx,ky,E_lh)
title('Light Hole'), xlabel('k_x'), ylabel('k_y'), shading flat
subplot(2,1,2)
surf(kx,ky,E_hh)
title('Heavy Hole'), xlabel('k_x'), ylabel('k_y'), shading flat