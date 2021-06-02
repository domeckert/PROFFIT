/*
 *  bknpow.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

double intfunc(double omega,double rf,double alpha,double xmin,double xmax){
	int nb=100;
	double *x=new double[nb+1];
	double logmin=log10(xmin);
	double logmax=log10(xmax);
	for (int i=0;i<nb+1;i++){
		double expon=logmin+i*(logmax-logmin)/nb;
		x[i]=TMath::Power(10.,expon);
	}
	double intot=0.0;
	for (int i=0;i<nb;i++){
		double z=(x[i]+x[i+1])/2.0;
		double width=x[i+1]-x[i];
		double term1=(omega*omega+z*z)/rf/rf;
		double term2=TMath::Power(term1,-alpha);
		intot+=term2*width;
	}
	delete [] x;
	return intot;
}

double intfuncbeta(double omega,double rc, double beta,double xmin,double xmax){
	int nb=100;
	double *x=new double[nb+1];
	double logmin=log10(xmin);
	double logmax=log10(xmax);
	for (int i=0;i<nb+1;i++){
		double expon=logmin+i*(logmax-logmin)/nb;
		x[i]=TMath::Power(10.,expon);
	}
	double intot=0.0;
	for (int i=0;i<nb;i++){
		double z=(x[i]+x[i+1])/2.0;
		double width=x[i+1]-x[i];
		double term1=1.+(omega*omega+z*z)/(rc*rc);
		double term2=TMath::Power(term1,-3*beta);
		intot+=term2*width;
	}
	delete [] x;
	return intot;
}

// Function for the projected surface-brightness profile expected from a beta and power-law density profile
double bknbeta(double *x,double *par){
	double beta=par[0];
	double rc=par[1];
	double alpha2=par[2];
	double rf=par[3];
	double A1=par[4];
	double jump=par[5];
	double A2=A1/jump/jump*TMath::Power(1.+rf*rf/rc/rc,-3*beta);
	double r=x[0]/rf;
	double sx;
	if (r<1.0){
		double term1=intfuncbeta(x[0],rc,beta,0.01,sqrt(rf*rf-x[0]*x[0]));
		double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
		sx=A1*term1+A2*term2;
	}
	else {
		double term=intfunc(x[0],rf,alpha2,0.01,1e3);
		sx=A2*term;
	}
	
	double cst=par[6];
	return sx+cst;	
}

// Function for the projected surface-brightness profile expected from a broken power-law density profile
double bknpow(double *x,double *par){
	double alpha1=par[0];
	double alpha2=par[1];
	double rf=par[2];
	double A1=par[3];
	double jump=par[4];
	double A2=A1/jump/jump;
	double r=x[0]/rf;
	double sx;
	if (r<1.0){
		double term1=intfunc(x[0],rf,alpha1,0.01,sqrt(rf*rf-x[0]*x[0]));
		double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
		sx=A1*term1+A2*term2;
	}
	else {
		double term=intfunc(x[0],rf,alpha2,0.01,1e3);
		sx=A2*term;
	}
	
	double cst=par[5];
	return sx+cst;	
}

double triplepl(double *x, double *par){
	double sl1=-par[0];
	double sl2=-par[1];
	double sl3=-par[2];
	double rc1=par[3];
	double rc2=par[4];
	double amp=par[5];
	double cst=par[6];
    double amp2=amp*pow(rc1,sl1-sl2);
    double amp3=amp2*pow(rc2,sl2-sl3);
    double xx;
    if (x[0]<rc1) {
        xx=amp*pow(x[0],sl1)+cst;
    }
    if (x[0]>=rc1 && x[0]<rc2) {
        xx=amp2*pow(x[0],sl2)+cst;
    }
    if (x[0]>=rc2) {
        xx=amp3*pow(x[0],sl3)+cst;
    }
	return xx;
}

double tripleplint(double *x, double *par){
	double sl1=-par[0];
	double sl2=-par[1];
	double sl3=-par[2];
	double rc1=par[3];
	double rc2=par[4];
	double amp=par[5];
	double cst=par[6];
    double amp2=amp*pow(rc1,sl1-sl2);
    double amp3=amp2*pow(rc2,sl2-sl3);
    double xx;
    if (x[0]<rc1) {
        xx=amp*pow(x[0],sl1)+cst;
    }
    if (x[0]>=rc1 && x[0]<rc2) {
        xx=amp2*pow(x[0],sl2)+cst;
    }
    if (x[0]>=rc2) {
        xx=amp3*pow(x[0],sl3)+cst;
    }
	return xx;
	return xx*2.*TMath::Pi()*x[0];
}


double intfunctriple(double omega,double sl1,double sl2,double sl3,double rc1,double rc2,double amp,double xmin,double xmax){
	int nb=200;
    double pars[7]={sl1,sl2,sl3,rc1,rc2,amp,0.0};
	double *x=new double[nb+1];
	double logmin=log10(xmin);
	double logmax=log10(xmax);
	for (int i=0;i<nb+1;i++){
		double expon=logmin+i*(logmax-logmin)/nb;
		x[i]=TMath::Power(10.,expon);
	}
	double intot=0.0;
	for (int i=0;i<nb;i++){
		double z=(x[i]+x[i+1])/2.0;
		double width=x[i+1]-x[i];
		double term1=sqrt(omega*omega+z*z);
		double term2=triplepl(&term1,pars);
		intot+=term2*term2*width;
	}
	delete [] x;
	return intot;
}


// Function for the projected surface-brightness profile expected from a broken power-law density profile
double triplebkn(double *x,double *par){
	double alpha1=par[0];
	double alpha2=par[1];
    double alpha3=par[2];
	double rf1=par[3];
    double rf2=par[4];
	double A1=par[5];
	double sx=intfunctriple(x[0],alpha1,alpha2,alpha3,rf1,rf2,A1,0.01,1e3);	
	double cst=par[6];
	return sx+cst;	
}

double triplebknint(double *x,double *par){
	double alpha1=par[0];
	double alpha2=par[1];
    double alpha3=par[2];
	double rf1=par[3];
    double rf2=par[4];
	double A1=par[5];
	double sx=intfunctriple(x[0],alpha1,alpha2,alpha3,rf1,rf2,A1,0.01,1e3);	
	double cst=par[6];
	return (sx+cst)*2.*TMath::Pi()*x[0];	
}


double bknbetaint(double *x,double *par){
	double beta=par[0];
	double rc=par[1];
	double alpha2=par[2];
	double rf=par[3];
	double A1=par[4];
	double jump=par[5];
	double A2=A1/jump/jump*TMath::Power(1.+rf*rf/rc/rc,-3*beta);
	double r=x[0]/rf;
	double sx;
	if (r<1.0){
		double term1=intfuncbeta(x[0],rc,beta,0.01,sqrt(rf*rf-x[0]*x[0]));
		double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
		sx=A1*term1+A2*term2;
	}
	else {
		double term=intfunc(x[0],rf,alpha2,0.01,1e3);
		sx=A2*term;
	}
	
	double cst=par[6];
	return (sx+cst)*2.*TMath::Pi()*x[0];	
}

double bknpowint(double *x,double *par){
	double alpha1=par[0];
	double alpha2=par[1];
	double rf=par[2];
	double A1=par[3];
	double jump=par[4];
	double A2=A1/jump/jump;
	double r=x[0]/rf;
	double sx;
	if (r<1.0){
		double term1=intfunc(x[0],rf,alpha1,0.01,sqrt(rf*rf-x[0]*x[0]));
		double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
		sx=A1*term1+A2*term2;
	}
	else {
		double term=intfunc(x[0],rf,alpha2,0.01,1e3);
		sx=A2*term;
	}
	
	double cst=par[5];
	return (sx+cst)*2.*TMath::Pi()*x[0];	
}

// Function for the projected surface-brightness profile expected from a broken power-law density profile
double bknpowgauss(double *x,double *par){
    double alpha1=par[0];
    double alpha2=par[1];
    double rf=par[2];
    double A1=par[3];
    double jump=par[4];
    double Ng=par[5];
    double mug=par[6];
    double sigmag=par[7];
    double A2=A1/jump/jump;
    double r=x[0]/rf;
    double sx;
    if (r<1.0){
        double term1=intfunc(x[0],rf,alpha1,0.01,sqrt(rf*rf-x[0]*x[0]));
        double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
        sx=A1*term1+A2*term2;
    }
    else {
        double term=intfunc(x[0],rf,alpha2,0.01,1e3);
        sx=A2*term;
    }
    double gaus=Ng/sqrt(2.*TMath::Pi()*sigmag*sigmag)*exp(-(x[0]-mug)*(x[0]-mug)/2./sigmag/sigmag);
    double cst=par[8];
    return sx+gaus+cst;
}

double bknpowgaussint(double *x,double *par){
    double alpha1=par[0];
    double alpha2=par[1];
    double rf=par[2];
    double A1=par[3];
    double jump=par[4];
    double Ng=par[5];
    double mug=par[6];
    double sigmag=par[7];
    double A2=A1/jump/jump;
    double r=x[0]/rf;
    double sx;
    if (r<1.0){
        double term1=intfunc(x[0],rf,alpha1,0.01,sqrt(rf*rf-x[0]*x[0]));
        double term2=intfunc(x[0],rf,alpha2,sqrt(rf*rf-x[0]*x[0]),1e3);
        sx=A1*term1+A2*term2;
    }
    else {
        double term=intfunc(x[0],rf,alpha2,0.01,1e3);
        sx=A2*term;
    }
    double gaus=Ng/sqrt(2.*TMath::Pi()*sigmag*sigmag)*exp(-(x[0]-mug)*(x[0]-mug)/2./sigmag/sigmag);
    double cst=par[8];
    return (sx+gaus+cst)*2.*TMath::Pi()*x[0];
}

