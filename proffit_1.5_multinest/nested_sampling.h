// Tools to run Multinest from Proffit

bool passstat;
double lnlike(double *pars){
	int npar,iflag;
	double *gin;
	double f;
	if (passstat){ //chi2
		fcncounts(npar,gin,f,pars,iflag);
	}
	else { //c-stat
		fcnlikehcounts(npar,gin,f,pars,iflag);
	}
    return -0.5*f; //logL instead of -2logL
}

double *passboundlow;
double *passboundhigh;
int passnpar;
int passnpfunc;
bool *passfix;
double *passpars;
bool *passislog;

void lnprior(double *pars){
  	/* convert unit Cube to actual parameter values */
	for (int i=0;i<passnpar;i++){
		pars[i]=passboundlow[i]+pars[i]*(passboundhigh[i]-passboundlow[i]);
	}
}

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	// Load boundaries for active parameters
	lnprior(Cube);
	// Pass total parameter variable (including fixed params) to the likelihood function
	double totCube[passnpfunc];
	int tpar=0;
	for (int i=0;i<passnpfunc;i++){
		if (passfix[i]){
			totCube[i]=passpars[i];
			//printf("Fix i, totCube: %d %g\n",i+1,totCube[i]);
		}
		else {
			if (passislog[i]){
				totCube[i]=pow(10.,Cube[tpar]);
			}
			else {
				totCube[i]=Cube[tpar];
			}
			//printf("Free i, totCube: %d %g\n",i+1,totCube[i]);
			tpar++;
		}
		//printf("totCube[%d]: %g\n",i,totCube[i]);
	}
	lnew=lnlike(totCube);
}

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

void run_multinest(TF1 *model,double *boundlow,double *boundhigh,bool *fix,bool *islogpar,char *statmet,char *dirname){
			// set the MultiNest sampling parameters
	
	
	// set the MultiNest sampling parameters
	int npfunc=model->GetNpar();
	int ndims=0; 					// dimensionality (no. of free parameters)
	passpars=new double[npfunc];
	double pars[npfunc];
	for (int i=0;i<npfunc;i++){
		passpars[i]=model->GetParameter(i);
		if (!fix[i]){
			pars[ndims]=0.5; //middle of the accepted range
			ndims++;
			/*if (islogpar[i]){
				printf("parameter %d  is log\n",i+1);
			}*/
		}
	}

	passnpar=ndims;
	passnpfunc=npfunc;
	passboundlow=boundlow;
	passboundhigh=boundhigh;
	passfix=fix;
	//passpars=pars;
	if (!strcmp(statmet,"chi2")) {
		passstat=true;				 
	}
	else {
		passstat=false;
	}
	passislog=islogpar;

	double llike;
	LogLike(pars,ndims,ndims,llike,0);
	printf("Original log likelihood: %g\n",llike);
	
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 0.8;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int nPar = ndims;					// total no. of parameters including free & derived parameters
	int nClsPar = ndims;				// no. of parameters to do mode separation on
	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
	char file[100];			// root for output files
	sprintf(file,"%s_",dirname);
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 0;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	
	// calling MultiNest

	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, file, seed, pWrap, fb, resume, outfile, initMPI,logZero, maxiter, LogLike, dumper, context);
}

void load_chains(char *name,int npar,bool *fix,int &npt,double **chains){
	do{
	char fname[200];
	sprintf(fname,"%s_post_equal_weights.dat",name);
	npt=line_num(fname);
	if (npt<1){
		printf("   Error: chains could not be read\n");
		break;
	}
	for (int np=0;np<npar;np++){
		if (!fix[np]){
			chains[np]=new double[npt];
		}
	}
	FILE *fchains=fopen(fname,"r");
	for (int i=0;i<npt;i++){
		for (int np=0;np<npar;np++){
			if (!fix[np]){
				fscanf(fchains,"%lf",&chains[np][i]);
			}
		}
		double likeh;
		fscanf(fchains,"%lf\n",&likeh);
	}
	fclose(fchains);
	}while(0);
}

void calc_envelope(double **chains,int npt,TF1 *model,double *bins,int nbin,bool *fix,bool *islogpar,double *vals,double *evals){
	double **allvals=new double*[nbin];
	for (int i=0;i<nbin;i++){
		allvals[i]=new double[npt];
	}
	int npar=model->GetNpar();
	for (int nc=0;nc<npt;nc++){
		for (int np=0;np<npar;np++){
			if (!fix[np]){
				if (!islogpar[np]){
					model->SetParameter(np,chains[np][nc]);
				}
				else {
					model->SetParameter(np,TMath::Power(10.,chains[np][nc]));
				}
			}
		}
		for (int i=0;i<nbin;i++){
			allvals[i][nc]=model->Eval(bins[i]);
		}
	}
	for (int i=0;i<nbin;i++){
		vals[i]=TMath::Mean(npt,allvals[i]);
		evals[i]=TMath::RMS(npt,allvals[i]);
	}
	for (int i=0;i<nbin;i++){
		delete [] allvals[i];
	}
	delete [] allvals;
}
