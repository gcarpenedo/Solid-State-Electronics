%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELLO di KRONIG-PENNEY
%
% Reticolo di potenziale periodico: Studio della struttura a bande
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
PhysConstants();
%% DEFINIZIONE DEL PROBLEMA: reticolo periodico
V0 = 0.1;                 %[eV]  altezza barriere
m  = m0;                  %[kg]
Nde = 10;                 %[-]
dE = V0/Nde;              %[eV]  step in energia

a=1e-10;                  %[m]   larghezza buca (w)
b=4e-10;                  %[m]   larghezza barriera --> periodo del reticolo (a+b)
a0=(2*m*q*V0)^0.5/hb;     %[m-1] alpha_0

N = floor(V0/dE);         % numero di energie studiate nel caso E<V
Nfin = 250;               % N_fin*N*dE = 25 eV è la max energia presa in considerazione

%% SOLUZIONE NUMERICA per K&P

% > Soluzioni per E<V0 (0<eta<1)
for k1=1:N
    
    E(k1)  = dE*k1;
    eta(k1)= E(k1)/V0;
    y(k1)  = cos(a0*a*(eta(k1)).^0.5).*cosh(a0*b*(1-eta(k1)).^0.5)+(1-2*eta(k1))/2./(eta(k1).*(1-eta(k1))).^0.5.*sin(a0*a*(eta(k1)).^0.5).*sinh(a0*b*(1-eta(k1)).^0.5);
    k(k1)  = 1e11;
    
    if abs(y(k1))<1  % -1 <y(k1)< 1
        k(k1)= acos(y(k1))/(a+b); % cos(k*(a+b)) = y(k1) --> k = arc_cos(y(k1))/(a+b)
        alpha(k1)=(2*m*q*E(k1))^0.5/hb;
        beta(k1)=(2*m*q*(E(k1)-V0))^0.5/hb;
        
        % sostituendo le condizioni al contorno alle condizioni di
        % periodicità, si ottiene il seguente sistema
        
        % A_a *a11 + B_a *a12 = 0;
        % A_a *a21 + B_a *a22 = 0;
        
        %dove
        
        a11(k1)= beta(k1)*sin(alpha(k1)*a)+alpha(k1)*exp(i*k(k1)*(a+b))*sin(beta(k1)*b);
        a12(k1)= beta(k1)*cos(alpha(k1)*a)- beta(k1)*exp(i*k(k1)*(a+b))*cos(beta(k1)*b);
        a21(k1)= alpha(k1)*cos(alpha(k1)*a)-alpha(k1)*exp(i*k(k1)*(a+b))*cos(beta(k1)*b);
        a22(k1)=-alpha(k1)*sin(alpha(k1)*a)-beta(k1)*exp(i*k(k1)*(a+b))*sin(beta(k1)*b);
        
        A_a(k1)=1; % tutte le altre costanti di integrazione sono espresse in funzione di A_a --> imponiamo per semplicità A_a = 1;
        B_a(k1)=-A_a(k1)*a11(k1)/a12(k1); % soluzione prima equazione del sistema
        A_b(k1)=alpha(k1)/beta(k1)*A_a(k1); % dalla seconda CC
        B_b(k1)=B_a(k1); % dalla prima CC
    end
end

% > Soluzioni per E>V0 (eta>1)
for k1=N:Nfin*N
    E(k1)=dE*k1;    %[eV]
    eta(k1)=E(k1)/V0;
    y(k1)=cos(a0*a*(eta(k1)).^0.5).*cos(a0*b*(eta(k1)-1).^0.5);
    y(k1)=y(k1)-(2*eta(k1)-1)/2./(eta(k1).*(eta(k1)-1)).^0.5.*sin(a0*a*(eta(k1)).^0.5).*sin(a0*b*(eta(k1)-1).^0.5);
    k(k1)=1e11;
    
    if abs(y(k1))<1 % -1 <y(k1)< 1
        k(k1)=acos(y(k1))/(a+b); % cos(k*(a+b)) = y(k1) --> k = arc_cos(y(k1))/(a+b)
        alpha(k1)=(2*m*q*E(k1))^0.5/hb;
        beta(k1)=(2*m*q*(E(k1)-V0))^0.5/hb;
        
        % sostituendo le condizioni al contorno alle condizioni di
        % periodicità, si ottiene il seguente sistema
        
        % A_a *a11 + B_a *a12 = 0;
        % A_a*a21 + B_a*a22 = 0;
        
        %dove
        
        a11(k1)= beta(k1)*sin(alpha(k1)*a)+alpha(k1)*exp(i*k(k1)*(a+b))*sin(beta(k1)*b);
        a12(k1)= beta(k1)*cos(alpha(k1)*a)- beta(k1)*exp(i*k(k1)*(a+b))*cos(beta(k1)*b);
        a21(k1)= alpha(k1)*cos(alpha(k1)*a)-alpha(k1)*exp(i*k(k1)*(a+b))*cos(beta(k1)*b);
        a22(k1)=-alpha(k1)*sin(alpha(k1)*a)-beta(k1)*exp(i*k(k1)*(a+b))*sin(beta(k1)*b);
        
        A_a(k1)=1; % tutte le altre costanti di integrazione sono espresse in funzione di A_a --> imponiamo per semplicità A_a = 1;
        B_a(k1)=-A_a(k1)*a11(k1)/a12(k1); % soluzione prima equazione del sistema
        A_b(k1)=alpha(k1)/beta(k1)*A_a(k1); % dalla seconda CC
        B_b(k1)=B_a(k1); % dalla prima CC
    end
end

%% RAPPRESENTAZIONE GRAFICA delle soluzioni
% Soluzione grafica (1)
figure(1)
subplot(2,1,1)
plot(eta,y)
axis([0 Nfin -1 1])
xlabel('\eta','FontSize',24), ylabel('cos(ka)','FontSize',24)
subplot(2,1,2)
plot(k,'ok'),ylabel('k','FontSize',24)
ylim([0 pi/(a+b)])

%% RAPPRESENTAZIONE GRAFICA: funzioni d'onda
figure(4)
N_a=100;
dx_a=a/N_a; %passo discretizzazione nella buca
x_a=(0:dx_a:a);
V_a=zeros(1,N_a+1);
dx_b=b/N_a; %passo discretizzazione nella barriera
x_b=(-b:dx_b:0);

% definizione potenziale periodico
V_b=ones(1,N_a+1)*V0;
Vdo=[x_b x_a x_b+a+b x_a+a+b x_b+2*(a+b) x_a+2*(a+b) x_b+3*(a+b) x_a+3*(a+b) x_b+4*(a+b) x_a+4*(a+b) x_b+5*(a+b) x_a+5*(a+b) x_b+6*(a+b) x_a+6*(a+b) x_b+7*(a+b) x_a+7*(a+b) x_b+8*(a+b) x_a+8*(a+b) x_b+9*(a+b) x_a+9*(a+b) x_b+10*(a+b) x_a+10*(a+b)
    V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a V_b V_a];
Vdo=Vdo';
%%

for k1 = 65
    
    psi_a = A_a(k1)*sin(alpha(k1)*x_a)+B_a(k1)*cos(alpha(k1)*x_a);
    psi_b = A_b(k1)*sin(beta(k1)*x_b)+B_b(k1)*cos(beta(k1)*x_b);
    
    % definizione funzione d'onda psi regione per regione mediante applicazione teorema di Bloch
    psi=[x_b x_a x_b+a+b x_a+a+b x_b+2*(a+b) x_a+2*(a+b) x_b+3*(a+b) x_a+3*(a+b) x_b+4*(a+b) x_a+4*(a+b) x_b+5*(a+b) x_a+5*(a+b) x_b+6*(a+b) x_a+6*(a+b) x_b+7*(a+b) x_a+7*(a+b) x_b+8*(a+b) x_a+8*(a+b) x_b+9*(a+b) x_a+9*(a+b) x_b+10*(a+b) x_a+10*(a+b)
        psi_b psi_a psi_b*exp(i*k(k1)*(a+b)) psi_a*exp(i*k(k1)*(a+b)) psi_b*exp(2*i*k(k1)*(a+b)) psi_a*exp(2*i*k(k1)*(a+b)) psi_b*exp(3*i*k(k1)*(a+b)) psi_a*exp(3*i*k(k1)*(a+b)) psi_b*exp(4*i*k(k1)*(a+b)) psi_a*exp(4*i*k(k1)*(a+b)) psi_b*exp(5*i*k(k1)*(a+b)) psi_a*exp(5*i*k(k1)*(a+b)) psi_b*exp(6*i*k(k1)*(a+b)) psi_a*exp(6*i*k(k1)*(a+b)) psi_b*exp(7*i*k(k1)*(a+b)) psi_a*exp(7*i*k(k1)*(a+b)) psi_b*exp(8*i*k(k1)*(a+b)) psi_a*exp(8*i*k(k1)*(a+b)) psi_b*exp(9*i*k(k1)*(a+b)) psi_a*exp(9*i*k(k1)*(a+b)) psi_b*exp(10*i*k(k1)*(a+b)) psi_a*exp(10*i*k(k1)*(a+b))];
    
    figure(4)
    subplot(1,3,1)
    plot(k,E,'.k')
    hold on
    plot(-k,E,'.k')
    plot(k(k1),E(k1),'*r')
    axis([-pi/(a+b) pi/(a+b) 0 Nfin*V0])
    xlabel('k'), ylabel('E')
    hold off
    subplot(4,3,2)
    plot(Vdo(:,1),Vdo(:,2))
    axis([-b 6*a+5*b -V0/5 2*V0]), title('\psi_k(x)')
    subplot(4,3,5)
    plot(psi(1,:),real(psi(2,:))), title('Re[\psi_k(x)]')
    axis([-b 6*a+5*b -4 4])
    subplot(4,3,8)
    plot(psi(1,:),imag(psi(2,:))), title('Im[\psi_k(x)]')
    axis([-b 6*a+5*b -4 4])
    subplot(4,3,11)
    plot(psi(1,:),abs(psi(2,:)).^2), title('|\psi_k(x)|^2')
    axis([-b 6*a+5*b 0 10]) 
    subplot(4,3,3)
    plot(Vdo(:,1),Vdo(:,2))
    axis([-b 6*a+5*b -V0/5 2*V0]), title('u_k(x), e^{ikx}')
    subplot(4,3,6)
    plot(psi(1,:),real(psi(2,:).*exp(-i*k(k1)*psi(1,:))),'b',psi(1,:),real(exp(i*k(k1)*psi(1,:))),'r')
    axis([-b 6*a+5*b -4 4])
    subplot(4,3,9)
    plot(psi(1,:),imag(psi(2,:).*exp(-i*k(k1)*psi(1,:))),'b',psi(1,:),imag(exp(i*k(k1)*psi(1,:))),'r')
    axis([-b 6*a+5*b -4 4]) 
    subplot(4,3,12)
    plot(psi(1,:),abs(psi(2,:).*exp(-i*k(k1)*psi(1,:))).^2,'b',psi(1,:),abs(exp(i*k(k1)*psi(1,:))).^2,'r')
    axis([-b 6*a+5*b 0 10])
    
    pause(0.01)
end
