function Mij=matr_trasf(ki,aj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcola la matrice Mij nel TRANSFER MATRIX METHOD                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mij=[exp(i*ki*aj)    exp(-i*ki*aj) 
     ki*exp(i*ki*aj) -ki*exp(-i*ki*aj)];