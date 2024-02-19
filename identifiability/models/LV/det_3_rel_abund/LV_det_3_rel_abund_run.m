% Symbolic parameters for identifiability analysis
syms gR1 I11 I12 I13 gR2 I21 I22 I23 gR3 I31 I32 I33

% Options
options.reportCompTime = true;

% Structural identifiability analysis
diary LV_det_3_rel_abund.txt
genssiMain('LV_det_3_rel_abund',7,[gR1;I11;I12;I13;gR2;I21;I22;I23;gR3;I31;I32;I33],options);
diary off