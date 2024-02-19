function model = LV_det_3types()

	% Symbolic variables
	syms n1 n2 n3
	syms gR1 I11 I12 I13 gR2 I21 I22 I23 gR3 I31 I32 I33

	% Parameters
	model.sym.p = [gR1;I11;I12;I13;gR2;I21;I22;I23;gR3;I31;I32;I33];

	% State variables
	model.sym.x = [n1;n2;n3];

	% Control vectors (g)
	model.sym.g = [0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(n1*(gR1+I11*n1+I12*n2+I13*n3));

		(n2*(gR2+I21*n1+I22*n2+I23*n3));

		(n3*(gR3+I31*n1+I32*n2+I33*n3));
];

	% Initial conditions
	model.sym.x0 = [700.0000;7000.0000;12000.0000];

	% Observables
	model.sym.y = [n1;n2;n3];

end