function model = logistic_det_3types()

	% Symbolic variables
	syms n1 n2 n3
	syms gR1 mR1 dR1 gR2 mR2 dR2 gR3 mR3 dR3 N

	% Parameters
	model.sym.p = [gR1;mR1;dR1;gR2;mR2;dR2;gR3;mR3;dR3;N];

	% State variables
	model.sym.x = [n1;n2;n3];

	% Control vectors (g)
	model.sym.g = [0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(((N-n1-n2-n3)*(gR1*n1+mR1)-dR1*n1)/N);

		(((N-n1-n2-n3)*(gR2*n2+mR2)-dR2*n2)/N);

		(((N-n1-n2-n3)*(gR3*n3+mR3)-dR3*n3)/N);
];

	% Initial conditions
	model.sym.x0 = [2800.0000;1400.0000;700.0000];

	% Observables
	model.sym.y = [n1;n2;n3];

end