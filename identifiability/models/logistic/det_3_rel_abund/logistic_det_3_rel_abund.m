function model = logistic_det_3types()

	% Symbolic variables
	syms x1 x2 x3 nS
	syms gR1 mR1 dR1 gR2 mR2 dR2 gR3 mR3 dR3 N

	% Parameters
	model.sym.p = [gR1;mR1;dR1;gR2;mR2;dR2;gR3;mR3;dR3;N];

	% State variables
	model.sym.x = [x1;x2;x3;nS];

	% Control vectors (g)
	model.sym.g = [0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(1/nS)*(((N-x1*nS-x2*nS-x3*nS)*(gR1*x1*nS+mR1)-dR1*x1*nS)/N-x1*(((N-x1*nS-x2*nS-x3*nS)*(gR1*x1*nS+mR1)-dR1*x1*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR2*x2*nS+mR2)-dR2*x2*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR3*x3*nS+mR3)-dR3*x3*nS)/N));

		(1/nS)*(((N-x1*nS-x2*nS-x3*nS)*(gR2*x2*nS+mR2)-dR2*x2*nS)/N-x2*(((N-x1*nS-x2*nS-x3*nS)*(gR1*x1*nS+mR1)-dR1*x1*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR2*x2*nS+mR2)-dR2*x2*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR3*x3*nS+mR3)-dR3*x3*nS)/N));

		(1/nS)*(((N-x1*nS-x2*nS-x3*nS)*(gR3*x3*nS+mR3)-dR3*x3*nS)/N-x3*(((N-x1*nS-x2*nS-x3*nS)*(gR1*x1*nS+mR1)-dR1*x1*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR2*x2*nS+mR2)-dR2*x2*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR3*x3*nS+mR3)-dR3*x3*nS)/N));

	((N-x1*nS-x2*nS-x3*nS)*(gR1*x1*nS+mR1)-dR1*x1*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR2*x2*nS+mR2)-dR2*x2*nS)/N+((N-x1*nS-x2*nS-x3*nS)*(gR3*x3*nS+mR3)-dR3*x3*nS)/N];

	% Initial conditions
	model.sym.x0 = [0.5714;0.2857;0.1429;4900.0000];

	% Observables
	model.sym.y = [x1;x2;x3];

end