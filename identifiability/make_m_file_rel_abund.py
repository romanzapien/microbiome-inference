## Make input file for identifiability analysis using GENSSI (in Matlab)

# Choose model
model = 'logistic'#logistic,LV
dynamics = 'stochastic'#deterministic,stochastic
n_types = 3

# Write file containing the chosen model
with open('%s_%s_%i_rel_abund.m'%(model,dynamics[:3],n_types), 'w') as f:
	f.write('function model = %s_%s_%itypes()'%(model,dynamics[:3],n_types))
	## vector of variables
	var_ = []
	for i in range(1,n_types+1): 
		var_.append('x%i'%i)
		if dynamics == 'stochastic':
			for j in range(1,n_types+1):
				if j>=i: var_.append('x%i%i'%(i,j))
	var_.append('nS')
	f.write('\n\n\t%% Symbolic variables\n\tsyms %s'%' '.join(var_))

	## vector of parameters
	par_ = []
	if model == 'logistic':
		for i in range(1,n_types+1): 
			par_.append('gR%i'%i)
			par_.append('mR%i'%i)
			par_.append('dR%i'%i)
		par_.append('N')
	if model == 'LV':
		for i in range(1,n_types+1): 
			par_.append('gR%i'%i)
			if dynamics == 'deterministic':
				for j in range(1,n_types+1): par_.append('I%i%i'%(i,j))
			else:
				for j in range(1,n_types+1): 
					par_.append('A%i%i'%(i,j))
					par_.append('B%i%i'%(i,j))
	f.write('\n\tsyms %s'%' '.join(par_))

	f.write('\n\n\t% Parameters')
	f.write('\n\tmodel.sym.p = [%s];'%';'.join(par_))

	f.write('\n\n\t% State variables')
	f.write('\n\tmodel.sym.x = [%s];'%';'.join(var_))

	f.write('\n\n\t% Control vectors (g)')
	f.write('\n\tmodel.sym.g = [%s];'%';'.join(len(var_)*['0']))#[];')#

	## model equations
	df_ = []
	nS_ = ''

	# scaling factor
	for i in range(1,n_types+1):
		
		if dynamics == 'deterministic':
			
			if model == 'logistic':
				nS_ += '((N'
				for j in range(1,n_types+1): nS_ += '-x%i*nS'%j
				nS_ += ')*(gR%i*x%i*nS+mR%i)-dR%i*x%i*nS)/N+'%(i,i,i,i,i)
			if model == 'LV':
				nS_ += 'x%i*nS*(gR%i'%(i,i)
				for j in range(1,n_types+1): nS_ += '+I%i%i*x%i*nS'%(i,j,j)
				nS_ += ')+'

		if dynamics == 'stochastic':

			if model == 'logistic':
				nS_ += 'gR%i*(x%i*nS'%(i,i)
				for j in range(1,n_types+1): 
					ij_s = sorted([i,j])
					nS_ += '-(x%i%i*nS.^2)/N'%(ij_s[0],ij_s[1])
				nS_ += ')+mR%i*(1'%i
				for j in range(1,n_types+1): nS_ += '-x%i*nS/N'%j
				nS_ += ')-dR%i*x%i*nS/N+'%(i,i)
			if model == 'LV':
				nS_ += 'gR%i*x%i*nS'%(i,i)
				for j in range(1,n_types+1):
					ij_s = sorted([i,j])
					nS_ += '+(A%i%i-B%i%i)*x%i%i*nS.^2'%(i,j,i,j,ij_s[0],ij_s[1])
				nS_ += '+'

	nS_ = nS_[:-1]

	for i in range(1,n_types+1):
		
		if dynamics == 'deterministic':
			dfi = '\t(1/nS)*('
			
			if model == 'logistic':
				dfi += '((N'
				for j in range(1,n_types+1): dfi += '-x%i*nS'%j
				dfi += ')*(gR%i*x%i*nS+mR%i)-dR%i*x%i*nS)/N-x%i*(%s));\n'%(i,i,i,i,i,i,nS_)
			
			if model == 'LV':
				dfi += 'x%i*nS*(gR%i'%(i,i)
				for j in range(1,n_types+1): dfi += '+I%i%i*x%i*nS'%(i,j,j)
				dfi += ')-x%i*(%s));\n'%(i,nS_)
			
			df_.append(dfi)
		
		if dynamics == 'stochastic':
			
			if model == 'logistic':
				# first moments
				dfi = '\t(1/nS)*(gR%i*(x%i*nS'%(i,i)
				for j in range(1,n_types+1): 
					ij_s = sorted([i,j])
					dfi += '-x%i%i*nS.^2/N'%(ij_s[0],ij_s[1])
				dfi += ')+mR%i*(1'%i
				for j in range(1,n_types+1): dfi += '-x%i*nS/N'%j
				dfi += ')-dR%i*x%i*nS/N-x%i*(%s));\n'%(i,i,i,nS_)
				df_.append(dfi)

				for j in range(1,n_types+1):
					# second moments
					if i == j:
						dfi = '\t(1/nS.^2)*(gR%i*(x%i*nS'%(i,i)
						for k in range(1,n_types+1):
							ik_s = sorted([i,k])
							dfi += '-x%i%i*nS.^2/N'%(ik_s[0],ik_s[1])
						dfi += '+2*(x%i%i*nS.^2'%(i,i)
						for k in range(1,n_types+1): dfi += '-x%i%i*x%i*nS.^3/N'%(i,i,k)
						dfi += '))+mR%i*(1'%i
						for k in range(1,n_types+1): dfi += '-x%i*nS/N'%(k)
						dfi += '+2*(x%i*nS'%(i)
						for k in range(1,n_types+1):
							ik_s = sorted([i,k])
							dfi += '-x%i%i*nS.^2/N'%(ik_s[0],ik_s[1])
						dfi += '))+dR%i/N*(x%i*nS-2*x%i%i*nS.^2)-2*x%i%i*nS*(%s));\n'%(i,i,i,i,i,i,nS_)
						df_.append(dfi)
					# co-moments
					else:
						if j>=i:
							ij_s = sorted([i,j])
							dfi = '\t(1/nS.^2)*((gR%i+gR%i)*(x%i%i*nS.^2'%(i,j,ij_s[0],ij_s[1])
							for k in range(1,n_types+1): dfi += '-x%i%i*x%i*nS.^3/N'%(ij_s[0],ij_s[1],k)
							dfi += ')+mR%i*(x%i*nS'%(i,j)
							for k in range(1,n_types+1):
								jk_s = sorted([j,k])
								dfi += '-x%i%i*nS.^2/N'%(jk_s[0],jk_s[1])
							dfi += ')+mR%i*(x%i*nS'%(j,i)
							for k in range(1,n_types+1): 
								ik_s = sorted([i,k])
								dfi += '-x%i%i*nS.^2/N'%(ik_s[0],ik_s[1])
							dfi += ')-(dR%i+dR%i)/N*x%i%i*nS.^2-2*x%i%i*nS*(%s));\n'%(i,j,ij_s[0],ij_s[1],ij_s[0],ij_s[1],nS_)
							df_.append(dfi)
			
			if model == 'LV':
				# first moments
				dfi = '\t(1/nS)*(gR%i*x%i*nS'%(i,i)
				for j in range(1,n_types+1):
					ij_s = sorted([i,j])
					dfi += '+(A%i%i-B%i%i)*x%i%i*nS.^2'%(i,j,i,j,ij_s[0],ij_s[1])
				dfi += '-x%i*(%s));\n'%(i,nS_)
				df_.append(dfi)
				for j in range(1,n_types+1):
					# second moments
					if i == j:
						dfi = '\t(1/nS.^2)*(gR%i*(x%i*nS+2*x%i%i*nS.^2)'%(i,i,i,i)
						for j in range(1,n_types+1):
							ij_s = sorted([i,j])
							dfi += '+(A%i%i+B%i%i)*x%i%i*nS.^2'%(i,j,i,j,ij_s[0],ij_s[1])
						for j in range(1,n_types+1):
							ij_s = sorted([i,j])
							dfi += '+2*(A%i%i-B%i%i)*x%i%i*x%i*nS.^3'%(i,j,i,j,ij_s[0],ij_s[1],i)
						dfi += '-2*x%i%i*nS*(%s));\n'%(i,i,nS_)
						df_.append(dfi)
					# co-moments
					else:
						if j>=i:
							ij_s = sorted([i,j])
							dfi = '\t(1/nS.^2)*((gR%i+gR%i)*x%i%i*nS.^2'%(i,j,ij_s[0],ij_s[1])
							for k in range(1,n_types+1): dfi += '+(A%i%i-B%i%i+A%i%i-B%i%i)*x%i%i*x%i*nS.^3'%(i,k,i,k,j,k,j,k,ij_s[0],ij_s[1],k)
							dfi += '-2*x%i%i*nS*(%s));\n'%(ij_s[0],ij_s[1],nS_)
							df_.append(dfi)
	df_.append(nS_)

	f.write('\n\n\t% Autonomous dynamics')
	f.write('\n\t model.sym.xdot = [%s];'%'\n\t'.join(df_))

	f.write('\n\n\t% Initial conditions')
	if n_types == 1: init = [4900]
	if n_types == 2: init = [2000, 2900]
	if n_types == 3:
		if model == 'logistic': init = [2800,1400,700]
		if model == 'LV': init = [700,7000,12000]

	nS_init = sum(init)
	init = [i/nS_init for i in init]
		
	if dynamics == 'deterministic': init_ = init
	else:
		init_ = []
		for i in range(n_types): 
			init_.append(init[i])
			for j in range(n_types):
				if j>=i: 
					init_.append(init[i]*init[j])

	init_.append(nS_init)
	init_ = ['%.4f'%i for i in init_]
	f.write('\n\tmodel.sym.x0 = [%s];'%';'.join(init_))

	f.write('\n\n\t% Observables')
	f.write('\n\tmodel.sym.y = [%s];'%';'.join(var_[:-1]))

	## end file
	f.write('\n\nend')


# Write file to run model 
with open('%s_%s_%i_rel_abund_run.m'%(model,dynamics[:3],n_types), 'w') as f:
	f.write('% Symbolic parameters for identifiability analysis')
	f.write('\nsyms %s'%' '.join(par_))

	f.write('\n\n% Options')
	f.write('\noptions.reportCompTime = true;')

	f.write('\n\n% Structural identifiability analysis')
	f.write('\ndiary %s_%s_%i_rel_abund.txt'%(model,dynamics[:3],n_types))
	f.write("\ngenssiMain('%s_%s_%i_rel_abund',7,[%s],options);"%(model,dynamics[:3],n_types,';'.join(par_)))
	f.write('\ndiary off')