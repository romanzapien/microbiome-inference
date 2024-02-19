function model = logistic_sto_3types()

	% Symbolic variables
	syms n1 n11 n12 n13 n2 n22 n23 n3 n33
	syms gR1 mR1 dR1 gR2 mR2 dR2 gR3 mR3 dR3 N

	% Parameters
	model.sym.p = [gR1;mR1;dR1;gR2;mR2;dR2;gR3;mR3;dR3;N];

	% State variables
	model.sym.x = [n1;n11;n12;n13;n2;n22;n23;n3;n33];

	% Control vectors (g)
	model.sym.g = [0;0;0;0;0;0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	gR1*(n1-n11/N-n12/N-n13/N)+mR1*(1-n1/N-n2/N-n3/N)-dR1*n1/N;

		gR1*(n1-n11/N-n12/N-n13/N+2*(n11-n11*n1/N-n11*n2/N-n11*n3/N))+mR1*(1-n1/N-n2/N-n3/N+2*(n1-n11/N-n12/N-n13/N))+dR1/N*(n1-2*n11);

		(gR1+gR2)*(n12-n12*n1/N-n12*n2/N-n12*n3/N)+mR1*(n2-n12/N-n22/N-n23/N)+mR2*(n1-n11/N-n12/N-n13/N)-(dR1+dR2)/N*n12;

		(gR1+gR3)*(n13-n13*n1/N-n13*n2/N-n13*n3/N)+mR1*(n3-n13/N-n23/N-n33/N)+mR3*(n1-n11/N-n12/N-n13/N)-(dR1+dR3)/N*n13;

		gR2*(n2-n12/N-n22/N-n23/N)+mR2*(1-n1/N-n2/N-n3/N)-dR2*n2/N;

		gR2*(n2-n12/N-n22/N-n23/N+2*(n22-n22*n1/N-n22*n2/N-n22*n3/N))+mR2*(1-n1/N-n2/N-n3/N+2*(n2-n12/N-n22/N-n23/N))+dR2/N*(n2-2*n22);

		(gR2+gR3)*(n23-n23*n1/N-n23*n2/N-n23*n3/N)+mR2*(n3-n13/N-n23/N-n33/N)+mR3*(n2-n12/N-n22/N-n23/N)-(dR2+dR3)/N*n23;

		gR3*(n3-n13/N-n23/N-n33/N)+mR3*(1-n1/N-n2/N-n3/N)-dR3*n3/N;

		gR3*(n3-n13/N-n23/N-n33/N+2*(n33-n33*n1/N-n33*n2/N-n33*n3/N))+mR3*(1-n1/N-n2/N-n3/N+2*(n3-n13/N-n23/N-n33/N))+dR3/N*(n3-2*n33);
];

	% Initial conditions
	model.sym.x0 = [2800.0000;7840000.0000;3920000.0000;1960000.0000;1400.0000;1960000.0000;980000.0000;700.0000;490000.0000];

	% Observables
	model.sym.y = [n1;n11;n12;n13;n2;n22;n23;n3;n33];

end