function model = LV_sto_3types()

	% Symbolic variables
	syms x1 x11 x12 x13 x2 x22 x23 x3 x33 nS
	syms gR1 B11 A12 B13 gR2 B21 B22 A23 gR3 A31 B32 B33

	% Parameters
	model.sym.p = [gR1;B11;A12;B13;gR2;B21;B22;A23;gR3;A31;B32;B33];

	% State variables
	model.sym.x = [x1;x11;x12;x13;x2;x22;x23;x3;x33;nS];

	% Control vectors (g)
	model.sym.g = [0;0;0;0;0;0;0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(1/nS)*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2-x1*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*(gR1*(x1*nS+2*x11*nS.^2)+(0+B11)*x11*nS.^2+(A12+0)*x12*nS.^2+(0+B13)*x13*nS.^2+2*(0-B11)*x11*x1*nS.^3+2*(A12-0)*x12*x1*nS.^3+2*(0-B13)*x13*x1*nS.^3-2*x11*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*((gR1+gR2)*x12*nS.^2+(0-B11+0-B21)*x12*x1*nS.^3+(A12-0+0-B22)*x12*x2*nS.^3+(0-B13+A23-0)*x12*x3*nS.^3-2*x12*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*((gR1+gR3)*x13*nS.^2+(0-B11+A31-0)*x13*x1*nS.^3+(A12-0+0-B32)*x13*x2*nS.^3+(0-B13+0-B33)*x13*x3*nS.^3-2*x13*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS)*(gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2-x2*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*(gR2*(x2*nS+2*x22*nS.^2)+(0+B21)*x12*nS.^2+(0+B22)*x22*nS.^2+(A23+0)*x23*nS.^2+2*(0-B21)*x12*x2*nS.^3+2*(0-B22)*x22*x2*nS.^3+2*(A23-0)*x23*x2*nS.^3-2*x22*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*((gR2+gR3)*x23*nS.^2+(0-B21+A31-0)*x23*x1*nS.^3+(0-B22+0-B32)*x23*x2*nS.^3+(A23-0+0-B33)*x23*x3*nS.^3-2*x23*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS)*(gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2-x3*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

		(1/nS.^2)*(gR3*(x3*nS+2*x33*nS.^2)+(A31+0)*x13*nS.^2+(0+B32)*x23*nS.^2+(0+B33)*x33*nS.^2+2*(A31-0)*x13*x3*nS.^3+2*(0-B32)*x23*x3*nS.^3+2*(0-B33)*x33*x3*nS.^3-2*x33*nS*(gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2));

	gR1*x1*nS+(0-B11)*x11*nS.^2+(A12-0)*x12*nS.^2+(0-B13)*x13*nS.^2+gR2*x2*nS+(0-B21)*x12*nS.^2+(0-B22)*x22*nS.^2+(A23-0)*x23*nS.^2+gR3*x3*nS+(A31-0)*x13*nS.^2+(0-B32)*x23*nS.^2+(0-B33)*x33*nS.^2];

	% Initial conditions
	model.sym.x0 = [0.0355;0.0013;0.0126;0.0216;0.3553;0.1263;0.2164;0.6091;0.3710;19700.0000];

	% Observables
	model.sym.y = [x1;x11;x12;x13;x2;x22;x23;x3;x33];

end