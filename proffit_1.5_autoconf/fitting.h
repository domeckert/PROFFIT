
double passfitlow,passfithigh;
TVirtualFitter *gFitter=0;
TH2F *histpsfmat;
bool ispsf=false;
TF1 *passfmod;
TH1F *passhist;
double *passbck;
bool isback=false;
double syserr=0.0;

double fastint(int np,double xm,double xp,TF1 *func,double *par){
	double width=(xp-xm)/np;
	double integ=0.0;
	for (int i=0; i<np; i++) {
		double center=xm+width/2.+width*i;
		double mnp=func->EvalPar(&center,par);
		integ+=mnp*width;
	}
	return integ;
}

void fcnstd(int &npar,double *gin,double &f,double *par,int iflag){
	Double_t x;
	
	passfmod->InitArgs(&x,par);
	npar = passfmod->GetNpar();
	f = 0;
	
	Int_t npfit = 0;
	Int_t nPoints=passhist->GetNbinsX();
	Double_t *df=new Double_t[nPoints];
	for (Int_t i=0;i<nPoints;i++) {
        double xmin=passhist->GetBinLowEdge(i+1);
        double xmax=passhist->GetBinLowEdge(i+2);
        double xminp=pow(xmin,3./2.);
        double xmaxp=pow(xmax,3./2.);
		x     = pow(0.5*(xminp+xmaxp),2./3.);//EM-weighted radius, McLaughlin 1999
		if (x>passfitlow && x<passfithigh) {
			TF1::RejectPoint(kFALSE);
			if (!ispsf){
				/*double xm=passhist->GetBinLowEdge(i+1);
				double xp=passhist->GetBinLowEdge(i+2);
				double area=TMath::Pi()*(xp*xp-xm*xm);//arcmin2
				//double mm=passfmod->Integral(xm,xp);
				double mm=fastint(10,xm,xp,passfmod,par);
				df[i] = mm-passhist->GetBinContent(i+1)*area;*/
				df[i] = passfmod->EvalPar(&x,par)-passhist->GetBinContent(i+1);
			}
			else {
				Double_t ll=0.0;
				for (int j=0; j<nPoints; j++) {
					if (!TF1::RejectedPoint()){
						x = passhist->GetBinCenter(j+1);
						ll+=passfmod->EvalPar(&x,par)*histpsfmat->GetBinContent(i+1,j+1);
					}
				}
				df[i]=ll-passhist->GetBinContent(i+1);
			}
			if (TF1::RejectedPoint()) df[i]=0.0;
			else npfit++;
		}
	}
	for (Int_t i=0;i<nPoints;i++) {
		x     = passhist->GetBinCenter(i+1);
		if (x>passfitlow && x<passfithigh) {
			/*double xm=passhist->GetBinLowEdge(i+1);
			double xp=passhist->GetBinLowEdge(i+2);
			double area=TMath::Pi()*(xp*xp-xm*xm);//arcmin2*/
			Double_t err=passhist->GetBinError(i+1);
			f += df[i]*df[i]/err/err;
		}
	}
	delete[] df;
	passfmod->SetNumberFitPoints(npfit);
    iflag=0;
    gin=NULL;
}


double passmaxexp;

void fcnstdlikeh(int &npar,double *gin,double &f,double *par,int iflag){
	Double_t x;
	
	passfmod->InitArgs(&x,par);
	npar = passfmod->GetNpar();
	f = 0;
	
	Int_t npfit = 0;
	Int_t nPoints=passhist->GetNbinsX();
	Double_t *df=new Double_t[nPoints];
	for (Int_t i=0;i<nPoints;i++) {
        double xmin=passhist->GetBinLowEdge(i+1);
        double xmax=passhist->GetBinLowEdge(i+2);
        double xminp=pow(xmin,3./2.);
        double xmaxp=pow(xmax,3./2.);
		x     = pow(0.5*(xminp+xmaxp),2./3.);//EM-weighted radius, McLaughlin 1999
		if (x>passfitlow && x<passfithigh) {
			TF1::RejectPoint(kFALSE);
			if (!ispsf){
				double mm=passfmod->EvalPar(&x,par)*passmaxexp;
                double sb=passhist->GetBinContent(i+1);
                if (sb>0.0) {
                    f += 2*(mm-sb*passmaxexp*log(mm)-sb*passmaxexp+sb*passmaxexp*log(sb*passmaxexp));
                }
                else {
                    f += 2*mm;
                }				
			}
			else {
				Double_t ll=0.0;
				for (int j=0; j<nPoints; j++) {
					if (!TF1::RejectedPoint()){
						x = passhist->GetBinCenter(j+1);
						ll+=passfmod->EvalPar(&x,par)*histpsfmat->GetBinContent(i+1,j+1);
					}
				}
				ll*=passmaxexp;
				f += 2*(ll-passhist->GetBinContent(i+1)*passmaxexp*log(ll));
			}
			if (TF1::RejectedPoint()) df[i]=0.0;
			else npfit++;
		}
	}
	passfmod->SetNumberFitPoints(npfit);
    iflag=0;
    gin=NULL;
}


void fcncounts(int &npar,double *gin,double &f,double *par,int iflag){
	Double_t x;
	
	passfmod->InitArgs(&x,par);
	npar = passfmod->GetNpar();
	f = 0;
	
	Int_t npfit = 0;
	Int_t nPoints=passhist->GetNbinsX();
	Double_t *df=new Double_t[nPoints];
	for (Int_t i=0;i<nPoints;i++) {
        double xmin=passhist->GetBinLowEdge(i+1);
        double xmax=passhist->GetBinLowEdge(i+2);
        double xminp=pow(xmin,3./2.);
        double xmaxp=pow(xmax,3./2.);
		x     = pow(0.5*(xminp+xmaxp),2./3.);//EM-weighted radius, McLaughlin 1999
		if (x>passfitlow && x<passfithigh) {
			TF1::RejectPoint(kFALSE);
			double cf=passhist->GetBinError(i+1);
			if (!ispsf){
				if (isback) {
					double bc=passbck[i];
					df[i] = passfmod->EvalPar(&x,par)*cf-passhist->GetBinContent(i+1)+bc;					
				}
				else {
					df[i] = passfmod->EvalPar(&x,par)*cf-passhist->GetBinContent(i+1);
				}
			}
			else {
				Double_t ll=0.0;
				for (int j=0; j<nPoints; j++) {
					if (!TF1::RejectedPoint()){
						x = passhist->GetBinCenter(j+1);
						ll+=passfmod->EvalPar(&x,par)*cf*histpsfmat->GetBinContent(i+1,j+1);
					}
				}
				if (isback) {
					double bc=passbck[i];
					df[i] = ll-passhist->GetBinContent(i+1)+bc;					
				}
				else {
					df[i] = ll-passhist->GetBinContent(i+1);
				}
			}
			if (TF1::RejectedPoint()) df[i]=0.0;
			else npfit++;
		}
	}
	for (Int_t i=0;i<nPoints;i++) {
		x     = passhist->GetBinCenter(i+1);
		if (x>passfitlow && x<passfithigh) {
			Double_t err=sqrt(passhist->GetBinContent(i+1)+syserr/100.*passhist->GetBinContent(i+1)*syserr/100.*passhist->GetBinContent(i+1));
			f += df[i]*df[i]/err/err;
		}
	}
	delete[] df;
	passfmod->SetNumberFitPoints(npfit);
    iflag=0;
    gin=NULL;
}


void fcnlikehcounts(int &npar,double *gin,double &f,double *par,int iflag){
	Double_t x;
	
	passfmod->InitArgs(&x,par);
	npar = passfmod->GetNpar();
	f = 0;
	
	Int_t npfit = 0;
	Int_t nPoints=passhist->GetNbinsX();
	Double_t *df=new Double_t[nPoints];
	for (Int_t i=0;i<nPoints;i++) {
        double xmin=passhist->GetBinLowEdge(i+1);
        double xmax=passhist->GetBinLowEdge(i+2);
        double xminp=pow(xmin,3./2.);
        double xmaxp=pow(xmax,3./2.);
		x     = pow(0.5*(xminp+xmaxp),2./3.);//EM-weighted radius, McLaughlin 1999
		if (x>passfitlow && x<passfithigh) {
			TF1::RejectPoint(kFALSE);
			double cf=passhist->GetBinError(i+1);
			if (!ispsf){
				if (isback) {
					double bc=passbck[i];
					double mm=passfmod->EvalPar(&x,par)*cf+bc;
                    int nc=passhist->GetBinContent(i+1);
                    if (nc>0.0) {
                        f += 2*(mm-nc*log(mm)-nc+nc*log(nc));
                    }
                    else {
                        f += 2*mm;
                    }
				}
				else {
					double mm=passfmod->EvalPar(&x,par)*cf;
                    int nc=passhist->GetBinContent(i+1);
                    if (nc>0.0) {
                        f += 2*(mm-nc*log(mm)-nc+nc*log(nc));
                    }
                    else {
                        f += 2*mm;
                    }
				}
			}
			else {
				Double_t ll=0.0;
				for (int j=0; j<nPoints; j++) {
					if (!TF1::RejectedPoint()){
						x = passhist->GetBinCenter(j+1);
						ll+=passfmod->EvalPar(&x,par)*cf*histpsfmat->GetBinContent(i+1,j+1);
					}
				}
				if (isback) {
					double bc=passbck[i];
					ll+=bc;
                    int nc=passhist->GetBinContent(i+1);
                    if (nc>0.0) {
                        f += 2*(ll-nc*log(ll)-nc+nc*log(nc));
                    }
                    else {
                        f += 2*ll;
                    }
				}
				else {
                    int nc=passhist->GetBinContent(i+1);
                    if (nc>0.0) {
                        f += 2*(ll-nc*log(ll)-nc+nc*log(nc));
                    }
                    else {
                        f += 2*ll;
                    }
				}
			}
			if (TF1::RejectedPoint()) df[i]=0.0;
			else npfit++;
		}
	}
	passfmod->SetNumberFitPoints(npfit);
    iflag=0;
    gin=NULL;
}

