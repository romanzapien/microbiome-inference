function model = LV_sto_3types()

	% Symbolic variables
	syms n1 n11 n12 n13 n2 n22 n23 n3 n33
	syms gR1 B11 A12 B13 gR2 B21 B22 A23 gR3 A31 B32 B33

	% Parameters
	model.sym.p = [gR1;B11;A12;B13;gR2;B21;B22;A23;gR3;A31;B32;B33];

	% State variables
	model.sym.x = [n1;n11;n12;n13;n2;n22;n23;n3;n33];

	% Control vectors (g)
	model.sym.g = [0;0;0;0;0;0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	gR1*n1+(0-B11)*n11+(A12-0)*n12+(0-B13)*n13;

		gR1*(n1+2*n11)+(0+B11)*n11+(A12+0)*n12+(0+B13)*n13+2*(0-B11)*n11*n1+2*(A12-0)*n12*n1+2*(0-B13)*n13*n1;

		(gR1+gR2)*n12+(0-B11+0-B21)*n12*n1+(A12-0+0-B22)*n12*n2+(0-B13+A23-0)*n12*n3;

		(gR1+gR3)*n13+(0-B11+A31-0)*n13*n1+(A12-0+0-B32)*n13*n2+(0-B13+0-B33)*n13*n3;

		gR2*n2+(0-B21)*n12+(0-B22)*n22+(A23-0)*n23;

		gR2*(n2+2*n22)+(0+B21)*n12+(0+B22)*n22+(A23+0)*n23+2*(0-B21)*n12*n2+2*(0-B22)*n22*n2+2*(A23-0)*n23*n2;

		(gR2+gR3)*n23+(0-B21+A31-0)*n23*n1+(0-B22+0-B32)*n23*n2+(A23-0+0-B33)*n23*n3;

		gR3*n3+(A31-0)*n13+(0-B32)*n23+(0-B33)*n33;

		gR3*(n3+2*n33)+(A31+0)*n13+(0+B32)*n23+(0+B33)*n33+2*(A31-0)*n13*n3+2*(0-B32)*n23*n3+2*(0-B33)*n33*n3;
];

	% Initial conditions
	model.sym.x0 = [700.0000;490000.0000;4900000.0000;8400000.0000;7000.0000;49000000.0000;84000000.0000;12000.0000;144000000.0000];

	% Observables
	model.sym.y = [n1;n11;n12;n13;n2;n22;n23;n3;n33];

end