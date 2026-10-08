

species = {'Pyr','OAA','AcCoA','Cit','CisAcon','IsoCit','aKG', ...
           'SucCoA','Succ','Fum','Mal','CoA','AcetylP','NAD','NADH'};

rxns = {'PK','PDH','CS','ACO1','ACO2','IDH','AKGDH', ...
        'SCS','SDH','F','MDH','PTA','ACKA','REDOX'};

S = [...
%   PK  PDH  CS  ACO1 ACO2 IDH AKGDH SCS SDH  F   MDH PTA ACKA Redox
     1  -1   0   0    0    0    0    0   0   0    0   0   0     0;  % Pyr
     0   0  -1   0    0    0    0    0   0   0    1   0   0     0;  % OAA
     0   1  -1   0    0    0    0    0   0   0    0  -1   0     0;  % AcCoA
     0   0   1  -1    0    0    0    0   0   0    0   0   0     0;  % Cit
     0   0   0   1   -1    0    0    0   0   0    0   0   0     0;  % CisAcon
     0   0   0   0    1   -1    0    0   0   0    0   0   0     0;  % IsoCit
     0   0   0   0    0    1   -1    0   0   0    0   0   0     0;  % aKG
     0   0   0   0    0    0    1   -1   0   0    0   0   0     0;  % SucCoA
     0   0   0   0    0    0    0    1  -1   0    0   0   0     0;  % Succ
     0   0   0   0    0    0    0    0   1  -1    0   0   0     0;  % Fum
     0   0   0   0    0    0    0    0   0   1   -1   0   0     0;  % Mal
     0  -1   1   0    0    0   -1    1   0   0    0   1   0     0;  % CoA
     0   0   0   0    0    0    0    0   0   0    0   1  -1     0;  % AcetylP
     0  -1   0   0    0   -1   -1    0   0   0   -1   0   0     1;  % NAD+
     0   1   0   0    0    1    1    0   0   0    1   0   0    -1]; % NADH



% left null space basis (conserved quantities)
L = null(S','r');


disp('Left null space basis L (columns are conserved quantities):')
disp(L)
