% Symbolic parameters for identifiability analysis
syms gR1 B11 A12 B13 gR2 B21 B22 A23 gR3 A31 B32 B33

% Options
options.reportCompTime = true;

% Structural identifiability analysis
diary LV_sto_3_rel_abund.txt
genssiMain('LV_sto_3_rel_abund',7,[gR1;B11;A12;B13;gR2;B21;B22;A23;gR3;A31;B32;B33],options);
diary off