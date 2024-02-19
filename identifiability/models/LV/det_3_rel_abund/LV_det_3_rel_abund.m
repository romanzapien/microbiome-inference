function model = LV_det_3types()

	% Symbolic variables
	syms x1 x2 x3 nS
	syms gR1 I11 I12 I13 gR2 I21 I22 I23 gR3 I31 I32 I33

	% Parameters
	model.sym.p = [gR1;I11;I12;I13;gR2;I21;I22;I23;gR3;I31;I32;I33];

	% State variables
	model.sym.x = [x1;x2;x3;nS];

	% Control vectors (g)
	model.sym.g = [0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(1/nS)*(x1*nS*(gR1+I11*x1*nS+I12*x2*nS+I13*x3*nS)-x1*(x1*nS*(gR1+I11*x1*nS+I12*x2*nS+I13*x3*nS)+x2*nS*(gR2+I21*x1*nS+I22*x2*nS+I23*x3*nS)+x3*nS*(gR3+I31*x1*nS+I32*x2*nS+I33*x3*nS)));

		(1/nS)*(x2*nS*(gR2+I21*x1*nS+I22*x2*nS+I23*x3*nS)-x2*(x1*nS*(gR1+I11*x1*nS+I12*x2*nS+I13*x3*nS)+x2*nS*(gR2+I21*x1*nS+I22*x2*nS+I23*x3*nS)+x3*nS*(gR3+I31*x1*nS+I32*x2*nS+I33*x3*nS)));

		(1/nS)*(x3*nS*(gR3+I31*x1*nS+I32*x2*nS+I33*x3*nS)-x3*(x1*nS*(gR1+I11*x1*nS+I12*x2*nS+I13*x3*nS)+x2*nS*(gR2+I21*x1*nS+I22*x2*nS+I23*x3*nS)+x3*nS*(gR3+I31*x1*nS+I32*x2*nS+I33*x3*nS)));

	x1*nS*(gR1+I11*x1*nS+I12*x2*nS+I13*x3*nS)+x2*nS*(gR2+I21*x1*nS+I22*x2*nS+I23*x3*nS)+x3*nS*(gR3+I31*x1*nS+I32*x2*nS+I33*x3*nS)];

	% Initial conditions
	model.sym.x0 = [0.0355;0.3553;0.6091;19700.0000];

	% Observables
	model.sym.y = [x1;x2;x3];

end