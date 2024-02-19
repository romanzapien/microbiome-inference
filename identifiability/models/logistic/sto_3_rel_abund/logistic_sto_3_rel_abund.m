function model = logistic_sto_3types()

	% Symbolic variables
	syms x1 x11 x12 x13 x2 x22 x23 x3 x33 nS
	syms gR1 mR1 dR1 gR2 mR2 dR2 gR3 mR3 dR3 N

	% Parameters
	model.sym.p = [gR1;mR1;dR1;gR2;mR2;dR2;gR3;mR3;dR3;N];

	% State variables
	model.sym.x = [x1;x11;x12;x13;x2;x22;x23;x3;x33;nS];

	% Control vectors (g)
	model.sym.g = [0;0;0;0;0;0;0;0;0;0];

	% Autonomous dynamics
	 model.sym.xdot = [	(1/nS)*(gR1*(x1*nS-x11*nS.^2/N-x12*nS.^2/N-x13*nS.^2/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N-x1*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*(gR1*(x1*nS-x11*nS.^2/N-x12*nS.^2/N-x13*nS.^2/N+2*(x11*nS.^2-x11*x1*nS.^3/N-x11*x2*nS.^3/N-x11*x3*nS.^3/N))+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N+2*(x1*nS-x11*nS.^2/N-x12*nS.^2/N-x13*nS.^2/N))+dR1/N*(x1*nS-2*x11*nS.^2)-2*x11*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*((gR1+gR2)*(x12*nS.^2-x12*x1*nS.^3/N-x12*x2*nS.^3/N-x12*x3*nS.^3/N)+mR1*(x2*nS-x12*nS.^2/N-x22*nS.^2/N-x23*nS.^2/N)+mR2*(x1*nS-x11*nS.^2/N-x12*nS.^2/N-x13*nS.^2/N)-(dR1+dR2)/N*x12*nS.^2-2*x12*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*((gR1+gR3)*(x13*nS.^2-x13*x1*nS.^3/N-x13*x2*nS.^3/N-x13*x3*nS.^3/N)+mR1*(x3*nS-x13*nS.^2/N-x23*nS.^2/N-x33*nS.^2/N)+mR3*(x1*nS-x11*nS.^2/N-x12*nS.^2/N-x13*nS.^2/N)-(dR1+dR3)/N*x13*nS.^2-2*x13*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS)*(gR2*(x2*nS-x12*nS.^2/N-x22*nS.^2/N-x23*nS.^2/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N-x2*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*(gR2*(x2*nS-x12*nS.^2/N-x22*nS.^2/N-x23*nS.^2/N+2*(x22*nS.^2-x22*x1*nS.^3/N-x22*x2*nS.^3/N-x22*x3*nS.^3/N))+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N+2*(x2*nS-x12*nS.^2/N-x22*nS.^2/N-x23*nS.^2/N))+dR2/N*(x2*nS-2*x22*nS.^2)-2*x22*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*((gR2+gR3)*(x23*nS.^2-x23*x1*nS.^3/N-x23*x2*nS.^3/N-x23*x3*nS.^3/N)+mR2*(x3*nS-x13*nS.^2/N-x23*nS.^2/N-x33*nS.^2/N)+mR3*(x2*nS-x12*nS.^2/N-x22*nS.^2/N-x23*nS.^2/N)-(dR2+dR3)/N*x23*nS.^2-2*x23*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS)*(gR3*(x3*nS-x13*nS.^2/N-x23*nS.^2/N-x33*nS.^2/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N-x3*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

		(1/nS.^2)*(gR3*(x3*nS-x13*nS.^2/N-x23*nS.^2/N-x33*nS.^2/N+2*(x33*nS.^2-x33*x1*nS.^3/N-x33*x2*nS.^3/N-x33*x3*nS.^3/N))+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N+2*(x3*nS-x13*nS.^2/N-x23*nS.^2/N-x33*nS.^2/N))+dR3/N*(x3*nS-2*x33*nS.^2)-2*x33*nS*(gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N));

	gR1*(x1*nS-(x11*nS.^2)/N-(x12*nS.^2)/N-(x13*nS.^2)/N)+mR1*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR1*x1*nS/N+gR2*(x2*nS-(x12*nS.^2)/N-(x22*nS.^2)/N-(x23*nS.^2)/N)+mR2*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR2*x2*nS/N+gR3*(x3*nS-(x13*nS.^2)/N-(x23*nS.^2)/N-(x33*nS.^2)/N)+mR3*(1-x1*nS/N-x2*nS/N-x3*nS/N)-dR3*x3*nS/N];

	% Initial conditions
	model.sym.x0 = [0.5714;0.3265;0.1633;0.0816;0.2857;0.0816;0.0408;0.1429;0.0204;4900.0000];

	% Observables
	model.sym.y = [x1;x11;x12;x13;x2;x22;x23;x3;x33];

end