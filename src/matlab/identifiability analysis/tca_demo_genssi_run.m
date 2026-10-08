clear;       
clc;          
close all; 


syms enzyme_pdh enzyme_pk enzyme_cs
syms enzyme_aco enzyme_akd
syms enzyme_id enzyme_scs enzyme_sd acetate
syms enzyme_fumarase enzyme_md enzyme_pta k alp enzyme_ack

%enzyme_ack  acetate

%qh2 q

% qh2 q enzyme_pta enzyme_ack

% Confirm execution
genssiAskForConfirmation(3);


result = genssiMain('tca_demo_genssi', 3); 
        
