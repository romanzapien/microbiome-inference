% Symbolic parameters for identifiability analysis
syms gR1 mR1 dR1 gR2 mR2 dR2 gR3 mR3 dR3 N

% Options
options.reportCompTime = true;

% Structural identifiability analysis
diary logistic_sto_3_rel_abund.txt
genssiMain('logistic_sto_3_rel_abund',7,[gR1;mR1;dR1;gR2;mR2;dR2;gR3;mR3;dR3;N],options);
diary off