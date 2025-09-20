function [E,psi] = es(x,V,Neig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funzione per la risoluzione dell'equazione di Schrodinger 1D per profili
% continui di potenziale, con condizione di annullamento al contorno
% > INPUT
% x    = vettore spaziale equispaziato              %[m]
% V    = vettore potenziale                         %[J]
% Neig = numero autovalori/autovettori richiesti    %[1]
% > OUTPUT
% E   = autovalori calcolati                        %[J]
% psi = autofunzioni normalizzate                   %[m-1/2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definizione delle costanti fisiche e variabili locali
PhysConstants();
a = x(2)-x(1);        %[m]  (passo di discretizzazione spaziale)
m = m0;               %[Kg]
N = length(V);
%% Costruzione della Matrice per la risoluzione dell'Eq. (costruisco tre matrici diagonali e le sommo)
H = diag([(-(h/2/pi/a)^2/2/m)*ones(1,N-2) 0],-1) + ...
     diag([0 (-(h/2/pi/a)^2/2/m)*ones(1,N-2)],1) + ...
      diag(((h/2/pi/a)^2/m+V));  %[J]

% Condizioni al contorno (Boundary conditions --> possible different choices)
H(1,1)=(h/2/pi/a)^2/m;
H(N,N)=(h/2/pi/a)^2/m;
%% Calcolo autovalori ed autovettori
[F,D] = eig(H);
W = diag(D);                                    % [J]
%% Ricerca/ordinamento degli stati confinati
[E, kk] = sort(W,'ascend');
E = E(1:Neig);
psi = F(:,kk);
psi = psi(:,1:Neig);