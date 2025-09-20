%% BS_SI.M [SILICIO - BANDE E POPOLAMENTO]
% riproduce la popolazione degli stati nella prima zona di Brilluoin
clear all
close all

% Definizione delle costanti
q=1.6e-19;%[C]
hb=6.626e-34/2/pi;%[Js]
m0=9.1e-31;%[kg]
k=1.38e-23;%[JK-1]
ml_Si=0.9163*m0;%[kg]
mt_Si=0.1905*m0;%[kg]
mhh_Si=0.537*m0;%[kg]
mlh_Si=0.153*m0;%[kg]
Nc_Si=3.22e19;%[cm-3]
Nv_Si=1.83e19;%[cm-3]
ml_Ge=1.588*m0;%[kg]
mt_Ge=0.08152*m0;%[kg]
mhh_Ge=0.347*m0;%[kg]
mlh_Ge=0.0429*m0;%[kg]
m_GaAs=0.067*m0;%[kg]
mhh_GaAs=0.51*m0;%[kg]
mlh_GaAs=0.082*m0;%[kg]
a_Si=5.43e-10;%[m] lattice constant
a_Ge=5.646e-10;%[m] lattice constant
a_GaAs=5.6533e-10;%[m] lattice constant
Ev=0;%[eV]
Ec=1.12;%[eV] Si gap

% Parametri liberi 

% scelgo parametri silicio
T=300;%[K] 
Ef= 0.56;%[eV]
Nc=Nc_Si;
Nv=Nv_Si;
Ntot=1e22/(Nc*Nv)^0.5/exp(-q*Ec/2/k/T)/exp(abs(Ef-Ec/2)*q/k/T);%[1] normalizzazione
ml=ml_Si;%[kg] 
mt=mt_Si;%[kg]
mhh=mhh_Si;%[kg]
mlh=mlh_Si;%[kg]
a=a_Si;%[m]
b=2*pi/a;%[m-1]
kx0=0.8*2*pi/a;%[m-1] coordinata kx del minimo CB-Si
ky0=0;%[m-1] coordinata ky del minimo CB-Si
kz0=0;%[m-1] coordinata kz del minimo CB-Si
CBmin=[kx0 0 0
-kx0 0 0
0 kx0 0
0 -kx0 0
0 0 kx0
0 0 -kx0];%[m-1] coordinate minimi CB

%% Grafica BZ
figure(1)
plot3([b b],[b/2 0],[0 b/2],'k')
hold on
plot3([b b],[-b/2 0],[0 -b/2],'k')
plot3([b b],[b/2 0],[0 -b/2],'k')
plot3([b b],[-b/2 0],[0 b/2],'k')

plot3([-b -b],[b/2 0],[0 b/2],'k')
plot3([-b -b],[-b/2 0],[0 -b/2],'k')
plot3([-b -b],[b/2 0],[0 -b/2],'k')
plot3([-b -b],[-b/2 0],[0 b/2],'k')

plot3([b/2 0],[b b],[0 b/2],'k')
plot3([-b/2 0],[b b],[0 -b/2],'k')
plot3([b/2 0],[b b],[0 -b/2],'k')
plot3([-b/2 0],[b b],[0 b/2],'k')

plot3([b/2 0],[-b -b],[0 b/2],'k')
plot3([-b/2 0],[-b -b],[0 -b/2],'k')
plot3([b/2 0],[-b -b],[0 -b/2],'k')
plot3([-b/2 0],[-b -b],[0 b/2],'k')

plot3([b/2 0],[0 b/2],[b b],'k')
plot3([-b/2 0],[0 -b/2],[b b],'k')
plot3([b/2 0],[0 -b/2],[b b],'k')
plot3([-b/2 0],[0 b/2],[b b],'k')

plot3([b/2 0],[0 b/2],[-b -b],'k')
plot3([-b/2 0],[0 -b/2],[-b -b],'k')
plot3([b/2 0],[0 -b/2],[-b -b],'k')
plot3([-b/2 0],[0 b/2],[-b -b],'k')

plot3([b b/2],[b/2 b],[0 0],'k')
plot3([-b -b/2],[-b/2 -b],[0 0],'k')
plot3([b b/2],[-b/2 -b],[0 0],'k')
plot3([-b -b/2],[b/2 b],[0 0],'k')

plot3([0 0],[b b/2],[b/2 b],'k')
plot3([0 0],[-b -b/2],[-b/2 -b],'k')
plot3([0 0],[b b/2],[-b/2 -b],'k')
plot3([0 0],[-b -b/2],[b/2 b],'k')

plot3([b/2 b],[0 0],[b b/2],'k')
plot3([-b/2 -b],[0 0],[-b -b/2],'k')
plot3([-b/2 -b],[0 0],[b b/2],'k')
plot3([b/2 b],[0 0],[-b -b/2],'k')

%plot3([CBmin(1,1) CBmin(2,1)],[CBmin(1,2) CBmin(2,2)],[CBmin(1,3) CBmin(2,3)],'r')
%plot3([CBmin(3,1) CBmin(4,1)],[CBmin(3,2) CBmin(4,2)],[CBmin(3,3) CBmin(4,3)],'r')
%plot3([CBmin(5,1) CBmin(6,1)],[CBmin(5,2) CBmin(6,2)],[CBmin(5,3) CBmin(6,3)],'r')

%% Calcolo - Banda di Conduzione
f_c=exp(-q*(Ec-Ef)/k/T);%[1]
f_v=exp(q*(Ev-Ef)/k/T);%[1]
Ne=100;%[1] numero di step in energia
dE=k*T/q/10;%[eV]

for ik=1:Ne
E(ik)=(ik-0.5)*dE;%[eV] rispetto alla CB

% Banda di Conduzione
f_cb=exp(-q*(E(ik)+Ec-Ef)/k/T);%[1]
n=f_cb*Ntot;
ak=(2*ml*q*E(ik))^0.5/hb;%[m-1] asse kx dell'ellissoide
bk=(2*mt*q*E(ik))^0.5/hb;%[m-1] asse ky-kz dell'ellissoide
[ ex, ey, ez ] = ellipsoidrand (ak,bk,bk,n);
plot3(ex+kx0,ey,ez,'.b')
[ ex, ey, ez ] = ellipsoidrand (ak,bk,bk,n);
plot3(ex-kx0,ey,ez,'.b')
[ ex, ey, ez ] = ellipsoidrand (bk,ak,bk,n);
plot3(ex,ey+kx0,ez,'.b')
[ ex, ey, ez ] = ellipsoidrand (bk,ak,bk,n);
plot3(ex,ey-kx0,ez,'.b')
[ ex, ey, ez ] = ellipsoidrand (bk,bk,ak,n);
plot3(ex,ey,ez+kx0,'.b')
[ ex, ey, ez ] = ellipsoidrand (bk,bk,ak,n);
plot3(ex,ey,ez-kx0,'.b')

% Banda di Valenza
f_vb=exp(q*(-E(ik)+Ev-Ef)/k/T);%[1]
n=f_vb*Ntot;
ak=(2*mhh*q*E(ik))^0.5/hb;%[m-1] asse kx-ky-kz dell'ellissoide (sfera)
[ ex, ey, ez ] = ellipsoidrand (ak,ak,ak,n);
plot3(ex,ey,ez,'.r')  %lacune pesanti
ak=(2*mlh*q*E(ik))^0.5/hb;%[m-1] asse kx-ky-kz dell'ellissoide (sfera)
[ ex, ey, ez ] = ellipsoidrand (ak,ak,ak,n);
plot3(ex,ey,ez,'.m') %lacune leggere

end

xlabel('k_x')
ylabel('k_y')
zlabel('k_z')