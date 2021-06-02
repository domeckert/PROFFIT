/*
 *  models.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

double betaprofile(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double amp=par[2];
	double cst=par[3];
	double expo=-3.*beta+1./2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
	return amp*xx+cst;
}

double cst(double *x, double *par){
    double b=x[0]*0.0;
	double cc=par[0];
	return cc+b;
}

double cstint(double *x, double *par){
	double cc=par[0];
	return cc*2.*TMath::Pi()*x[0];
}

double powerlaw(double *x, double *par){
	double slope=-par[0];
    double rs=par[1];
	double amp=par[2];
	double cst=par[3];
	double xx=amp*pow(x[0]/rs,slope)+cst;
	return xx;
}

double powerlawint(double *x, double *par){
	double slope=-par[0];
    double rs=par[1];
    double amp=par[2];
    double cst=par[3];
	double xx=amp*pow(x[0]/rs,slope)+cst;
	return xx*2.*TMath::Pi()*x[0];
}

double betaint(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double amp=par[2];
	double cst=par[3];
	double expo=-3.*beta+1./2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
	return (amp*xx+cst)*2.*TMath::Pi()*x[0];
}

double doublebeta(double *x, double *par){
	double beta=par[0];
	double rc1=par[1];
	double rc2=par[2];
	double rat=par[3];
	double amp=par[4];
	double cst=par[5];
	double expo=-3*beta+0.5;
	double base=1.+(x[0]/rc1)*(x[0]/rc1);
	double xx1=pow(base,expo);
	base=1.+(x[0]/rc2)*(x[0]/rc2);
	double xx2=pow(base,expo);
	return amp*(xx1+rat*xx2)+cst;
}

double dbetaint(double *x, double *par){
	double beta=par[0];
	double rc1=par[1];
	double rc2=par[2];
	double rat=par[3];
	double amp=par[4];
	double cst=par[5];
	double expo=-3*beta+0.5;
	double base=1.+(x[0]/rc1)*(x[0]/rc1);
	double xx1=pow(base,expo);
	base=1.+(x[0]/rc2)*(x[0]/rc2);
	double xx2=pow(base,expo);
	return (amp*(xx1+rat*xx2)+cst)*2.*TMath::Pi()*x[0];
}

double cuspbeta(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double slope=par[2];
    double rs=par[3];
	double amp=par[4];
	double cst=par[5];
	double expo=-3.*beta+1./2.+slope/2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
	double cusp=pow(x[0]/rs,-slope);
	return amp*cusp*xx+cst;
}

double cuspbetaint(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double slope=par[2];
    double rs=par[3];
	double amp=par[4];
	double cst=par[5];
	double expo=-3.*beta+1./2.+slope/2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
	double cusp=pow(x[0]/rs,-slope);
	return (amp*cusp*xx+cst)*2.*TMath::Pi()*x[0];
}

double gausbeta(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double amp=par[2];
    double mug=par[3];
	double sigma=par[4];
    double ampg=par[5];
	double cst=par[6];
	double expo=-3.*beta+1./2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
	return amp*xx+ampg/sqrt(2*TMath::Pi())/sigma*exp(-(x[0]-mug)*(x[0]-mug)/2/sigma/sigma)+cst;
}

double gausbetaint(double *x, double *par){
	double beta=par[0];
	double rc=par[1];
	double amp=par[2];
    double mug=par[3];
    double sigma=par[4];
    double ampg=par[5];
    double cst=par[6];
	double expo=-3.*beta+1./2.;
	double base=1.+x[0]*x[0]/rc/rc;
	double xx=pow(base,expo);
    return (amp*xx+ampg/sqrt(2*TMath::Pi())/sigma*exp(-(x[0]-mug)*(x[0]-mug)/2/sigma/sigma)+cst)*2.*TMath::Pi()*x[0];
}

double gausdbeta(double *x, double *par){
	double beta=par[0];
	double rc1=par[1];
	double rc2=par[2];
	double rat=par[3];
	double amp=par[4];
    double mug=par[5];
    double sigma=par[6];
    double ampg=par[7];
	double cst=par[8];
	double expo=-3*beta+0.5;
	double base=1.+(x[0]/rc1)*(x[0]/rc1);
	double xx1=pow(base,expo);
	base=1.+(x[0]/rc2)*(x[0]/rc2);
	double xx2=pow(base,expo);
	return amp*(xx1+rat*xx2)+ampg/sqrt(2*TMath::Pi())/sigma*exp(-(x[0]-mug)*(x[0]-mug)/2/sigma/sigma)+cst;
}

double gausdbetaint(double *x, double *par){
	double beta=par[0];
	double rc1=par[1];
	double rc2=par[2];
    double rat=par[3];
    double amp=par[4];
    double mug=par[5];
    double sigma=par[6];
    double ampg=par[7];
    double cst=par[8];
	double expo=-3*beta+0.5;
	double base=1.+(x[0]/rc1)*(x[0]/rc1);
	double xx1=pow(base,expo);
	base=1.+(x[0]/rc2)*(x[0]/rc2);
	double xx2=pow(base,expo);
	return (amp*(xx1+rat*xx2)+ampg/sqrt(2*TMath::Pi())/sigma*exp(-(x[0]-mug)*(x[0]-mug)/2/sigma/sigma)+cst)*2.*TMath::Pi()*x[0];
}

int getbin2(double x,double *bins,double *ebins){
	int in=0;
	while (bins[in]+ebins[in]<x) {
		in++;
	}
	return in;
}

double *passback=NULL;
double *passbins;
double *passebins;

double funcback(double *x, double *par){
	double cst=par[0];
	double cback=par[1];
	int bin=getbin2(x[0],passbins,passebins);
	double fc=cst+cback*passback[bin];
	return fc;
}

double funcbackint(double *x, double *par){
	double cst=par[0];
	double cback=par[1];
	int bin=getbin2(x[0],passbins,passebins);
	double fc=cst+cback*passback[bin];
	return fc*2.*TMath::Pi()*x[0];
}

double bknpl(double *x, double *par){
	double sl1=-par[0];
	double sl2=-par[1];
	double rc1=par[2];
	double jump=par[4];
	double amp=par[3];
    double amp2=amp*pow(rc1,sl1-sl2)/jump;
    double xx;
    if (x[0]<rc1) {
        xx=amp*pow(x[0],sl1);
    }
    else {
        xx=amp2*pow(x[0],sl2);
    }
	return xx;
}

TF1* mkftemp(char *modname,char *funcname){
	TF1 *ftemp=NULL;
	if (!strcmp(modname,"beta")){
		ftemp=new TF1(funcname,betaprofile,1e-4,1e4,4);
	}
	if (!strcmp(modname,"doublebeta")){
		ftemp=new TF1(funcname,doublebeta,1e-4,1e4,6);
	}
	if (!strcmp(modname,"cuspbeta")){
		ftemp=new TF1(funcname,cuspbeta,1e-4,1e4,6);
	}
	if (!strcmp(modname,"const")){
		ftemp=new TF1(funcname,cst,1e-4,1e4,1);
	}
	if (!strcmp(modname,"power")){
		ftemp=new TF1(funcname,powerlaw,1e-4,1e4,4);
	}
	if (!strcmp(modname,"gausbeta")){
		ftemp=new TF1(funcname,gausbeta,1e-4,1e4,7);
	}
	if (!strcmp(modname,"gausdbeta")){
		ftemp=new TF1(funcname,gausdbeta,1e-4,1e4,9);
	}
	if (!strcmp(modname,"bknpow")){
		ftemp=new TF1(funcname,bknpow,1e-4,1e4,6);
	}
	if (!strcmp(modname,"bknbeta")){
		ftemp=new TF1(funcname,bknbeta,1e-4,1e4,7);
	}
	if (!strcmp(modname,"backfit")){
		ftemp=new TF1(funcname,funcback,1e-4,1e4,2);
	}	
	if (!strcmp(modname,"triplepl")){
		ftemp=new TF1(funcname,triplepl,1e-4,1e4,7);
	}	
	if (!strcmp(modname,"triplebkn")){
		ftemp=new TF1(funcname,triplebkn,1e-4,1e4,7);
	}	
    if (!strcmp(modname,"bknpowgauss")){
        ftemp=new TF1(funcname,bknpowgauss,1e-4,1e4,9);
    }
	return ftemp;
}	

TF1* mkfint(char *modname,char *funcname){
	TF1 *ftemp=NULL;
	if (!strcmp(modname,"beta")){
		ftemp=new TF1(funcname,betaint,1e-4,1e4,4);
	}
	if (!strcmp(modname,"doublebeta")){
		ftemp=new TF1(funcname,dbetaint,1e-4,1e4,6);
	}
	if (!strcmp(modname,"cuspbeta")){
		ftemp=new TF1(funcname,cuspbetaint,1e-4,1e4,6);
	}
	if (!strcmp(modname,"const")){
		ftemp=new TF1(funcname,cstint,1e-4,1e4,1);
	}
	if (!strcmp(modname,"power")){
		ftemp=new TF1(funcname,powerlawint,1e-4,1e4,4);
	}
	if (!strcmp(modname,"gausbeta")){
		ftemp=new TF1(funcname,gausbetaint,1e-4,1e4,7);
	}
	if (!strcmp(modname,"gausdbeta")){
		ftemp=new TF1(funcname,gausdbetaint,1e-4,1e4,9);
	}
	if (!strcmp(modname,"bknpow")){
		ftemp=new TF1(funcname,bknpowint,1e-4,1e4,6);
	}
	if (!strcmp(modname,"bknbeta")){
		ftemp=new TF1(funcname,bknbetaint,1e-4,1e4,7);
	}
	if (!strcmp(modname,"backfit")){
		ftemp=new TF1(funcname,funcbackint,1e-4,1e4,2);
	}	
	if (!strcmp(modname,"triplepl")){
		ftemp=new TF1(funcname,tripleplint,1e-4,1e4,7);
	}	
	if (!strcmp(modname,"triplebkn")){
		ftemp=new TF1(funcname,triplebknint,1e-4,1e4,7);
	}	
    if (!strcmp(modname,"bknpowgauss")){
        ftemp=new TF1(funcname,bknpowgaussint,1e-4,1e4,9);
    }
	return ftemp;
}	


