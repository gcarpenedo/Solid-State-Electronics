function [w,vg,beta] = rel_disp(k,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcola la pulsazione angolare a patire da k                            %
% <>> INPUT                                                               %
%  > k   : valori di k su cui calcolare la dispersione [m-1]              %
%  > type: specifica la relazione di dispersione da considerarsi:         %
%    - 'free'      : relaz. di dispers. di elettrone libero;              %
%    - 'linear'    : relaz. di dispers. lineare;                          %
%    - 'sublinear' : relaz. di dispers. sub-lineare.                      %
% <>> OUTPUT                                                              %
%  > w   : pulsazione angolare (k)                                        %
%  > vg  : velocita' di gruppo (k)                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysConstants();
m = m0;
if(strcmp(type,'free'))          % relaz. di dispers. di elettrone libero

    w = hb*k.^2/m/2.d0;
    vg = hb*k/m; %[ms-1]
    beta = hb/2/m;

elseif(strcmp(type,'linear'))    % relaz. di dispers. lineare
    
    c=1.d5;
    w=c*k;
    vg=c*ones(size(k));
    beta = 0;
    
elseif(strcmp(type,'sublinear')) % relaz. di dispers. sub-lineare
    
    c=1.d5;
    d=5e-5;
    w=c*k-d*k.^2/2;
    vg=c-d*k;
    beta = -d/2;
    
end


