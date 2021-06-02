/*
 *  miscellaneous.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */
const double omegam=0.3; //Omega_m
const double omegal=0.7; //Omega_lambda
const double omegar=2.48e-5; //Omega_radiation (CMB)
const double c=2.998E10; // speed of light in cm/s
const double Mpc=3.0856776e+24; //1 Mpc in cm
const double h0=70.*1e5; //cm/s/Mpc
const double xmmpsf_r0=5.16; //r0 of the XMM/EPIC PSF at the centre, in arcsec
const double xmmpsf_alpha=1.52; //alpha of the XMM/EPIC PSF at the centre


double gausint(double *x, double *par){
	double ampg=par[0];
	double sigma=par[1];
	return ampg*exp(-x[0]*x[0]/2/sigma/sigma)*2.*TMath::Pi()*x[0];
}

double cosmofunc(double *x, double *par){
    par=NULL;
	double unpz=1+x[0];
	double pp=omegar*pow(unpz,4.)+omegam*pow(unpz,3.)+omegal;
	double xx=1./sqrt(pp);
	return xx;
}

double kingfunction(double *x, double *par){
	double amp=par[0];
	double r0=par[1];
	double alpha=par[2];
	double xx=1.+(x[0]/r0)*(x[0]/r0);
	double yy=pow(xx,-alpha);
	return amp*yy;
}

double angdist(double z){
	TF1* ff=new TF1("ff",cosmofunc,0.,z+0.2,1);
	double resint=ff->Integral(0.,z);
	double d=c/h0/(1.+z)*resint;
	delete ff;
	return d;
}

int getbin(double dist,double *bins,double *ebins,int nbin){
    int tbin=0;
    while (dist>bins[tbin]+ebins[tbin] && tbin<nbin){
        tbin++;
    }
    return tbin;
}

int line_num(char *filename){
	FILE *f=fopen(filename, "r");;
	char c;
	int lines = 0;
	if(f == NULL) return 0;
	while((c = fgetc(f)) != EOF){
		if(c == '\n') lines++;
	}
	fclose(f);
	return lines;
}

double log12(double x){
	double tt=log(x)/log(1.2);
	return tt;
}

void backsub(double back,double eback,int nbin,double *profile,double *eprof){
	for (int i=0; i<nbin; i++) {
		double pr=profile[i];
		double epr=eprof[i];
		profile[i]=pr-back;
		eprof[i]=sqrt(epr*epr+eback*eback);
	}
}


// routine to calculate the nearest zero of f1-f2 using Newton technique 
double findrad(TF1* f1,TF1* f2,double root){
	double tol=1e-4;
	double grad=1.0;
	double diff=100.0;
	int niter=0;
	double pt1=1e-6;
	double rold;
	int zero=0;
	while (fabs(diff)>tol){
		double x1=f1->Eval(root);
		double x2=f1->Eval(root+pt1);
		double y1=f2->Eval(root);
		double y2=f2->Eval(root+pt1);
		grad=(x2-y2-x1+y1)/pt1;
		niter++;
		//      if (niter>100) break;
		if (grad!=0.0){
			rold=root;
			root=root-(x1-y1)/grad;
		}
		else break;
		if (root<0.0) root=-root;
		if (root<tol) {
			if (zero==0) {
				root=20.;
				zero++;
			}
			else {
				return root;
			}
		}
		diff=fabs(root-rold);
		//      printf("niter, root, grad, diff: %d %f %f %f\n",niter,root,grad,diff);
	}
	return root;
}

void getflux(double *profile,double  *eprof,double *bins, double *ebins,double bkg, double ebkg, double rmin,double rmax,double &flux,double &eflux){
	int i=0;
	double f=0.0;
	double ef=0.0;
	double r=bins[i]+ebins[i];
	while (r<rmin) {
		i++;
		r=bins[i]+ebins[i];
	}
	while (r<rmax){
		double area=(profile[i]-bkg)*2*TMath::Pi()*bins[i]*2*ebins[i];
		double earea=sqrt(eprof[i]*eprof[i]+ebkg*ebkg)*2*TMath::Pi()*bins[i]*2*ebins[i];
		ef+=earea*earea;
		f+=area;
		i++;
		r=bins[i]+ebins[i];
	}
	eflux=sqrt(ef);
	flux=f;
	
}


double backestimate(double *img,double *exposure,double rin,double rout,long *axes){
	double mean=0.0;
	int npix=0;
	double centerx=(double)axes[0]/2.0;
	double centery=(double)axes[1]/2.0;
	double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
	for (int i=0; i<axes[0]; i++) {
		for (int j=0; j<axes[1]; j++) {
			double dist=sqrt((i-centerx)*(i-centerx)+(j-centery)*(j-centery));//pixel
			if ((dist<rout) && (dist>rin) && (exposure[j*axes[0]+i]>0.25*maxexp)) {
				mean+=img[j*axes[0]+i];
				npix++;
			}
		}
	}
	mean/=npix;
	return mean;
}

double maxinrange(TH1F* hh,double xmin,double xmax){
	double mm=hh->GetBinContent(1);
	int npoints=hh->GetNbinsX();
	for (int i=0; i<npoints; i++) {
		double xx=hh->GetBinCenter(i+1);
		double val=hh->GetBinContent(i+1);
		if (val>mm && xx>xmin && xx<xmax){
			mm=val;
		}
	}
	return mm;
}

double mininrange(TH1F* hh,double xmin,double xmax){
	double mm=hh->GetBinContent(1);
	int npoints=hh->GetNbinsX();
	for (int i=0; i<npoints; i++) {
		double xx=hh->GetBinCenter(i+1);
		double val=hh->GetBinContent(i+1);
		if (val<mm && xx>xmin && xx<xmax){
			mm=val;
		}
	}
	return mm;
}

double meanx(double* img,double *exposure,double cra,double cdec,double rmax,long *axes,double angle){
	double cx=0.0;
	double ftot=0.0;
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double dist=sqrt((i-cra)*(i-cra)+(j-cdec)*(j-cdec));
			if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
				double xtil=cos(angle)*(i-cra)+sin(angle)*(j-cdec);
				cx+=img[j*axes[0]+i]*xtil;
				ftot+=img[j*axes[0]+i];
			}
		}
	}
	cx/=ftot;
	return cx;
}

double meany(double* img,double *exposure,double cra,double cdec,double rmax,long *axes,double angle){
	double cx=0.0;
	double ftot=0.0;
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double dist=sqrt((i-cra)*(i-cra)+(j-cdec)*(j-cdec));
			if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
				double ytil=-sin(angle)*(i-cra)+cos(angle)*(j-cdec);
				cx+=img[j*axes[0]+i]*ytil;
				ftot+=img[j*axes[0]+i];
			}
		}
	}
	cx/=ftot;
	return cx;
}

double meanx2(double* img,double *exposure,double meanx,double cra,double cdec,double rmax,long *axes,double angle){
	double cx=0.0;
	double ftot=0.0;
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double dist=sqrt((i-cra)*(i-cra)+(j-cdec)*(j-cdec));
			if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
				double xtil=cos(angle)*(i-cra)+sin(angle)*(j-cdec);
				cx+=img[j*axes[0]+i]*xtil*xtil;
				ftot+=img[j*axes[0]+i];
			}
		}
	}
	cx/=ftot;
	double mx2=cx-meanx*meanx;
	return mx2;
}

double meany2(double* img,double *exposure,double meany,double cra,double cdec,double rmax,long *axes,double angle){
	double cx=0.0;
	double ftot=0.0;
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double dist=sqrt((i-cra)*(i-cra)+(j-cdec)*(j-cdec));
			if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
				double ytil=-sin(angle)*(i-cra)+cos(angle)*(j-cdec);
				cx+=img[j*axes[0]+i]*ytil*ytil;
				ftot+=img[j*axes[0]+i];
			}
		}
	}
	cx/=ftot;
	double my2=cx-meany*meany;
	return my2;
}

double meanxy(double* img,double *exposure,double meanx,double meany,double cra,double cdec,double rmax,long *axes,double angle){
	double cx=0.0;
	double ftot=0.0;
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double dist=sqrt((i-cra)*(i-cra)+(j-cdec)*(j-cdec));
			if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
				double xtil=cos(angle)*(i-cra)+sin(angle)*(j-cdec);
				double ytil=-sin(angle)*(i-cra)+cos(angle)*(j-cdec);
				cx+=img[j*axes[0]+i]*xtil*ytil;
				ftot+=img[j*axes[0]+i];
			}
		}
	}
	cx/=ftot;
	double mxy=cx-meanx*meany;
	return mxy;
}

void invertp(double *inp,double *outp,int np){
	for (int i=0; i<np; i++) {
		int k=np-i-1;
		outp[k]=inp[i];
	}
}

double vij(double *binsh,int i,int j){
	double fact=0.0;
	if (j>i){
		double cd=binsh[i]*binsh[i]-binsh[j]*binsh[j];
		fact=4./3.*TMath::Pi()*TMath::Power(cd,3./2.);
	}
	else {
		fact=0.0;
	}
	return fact;
}

double* medsmooth(int nbin,double *profile,int width){
    if (nbin+1<width) {
        printf("Error: Size of smoothing window must be smaller than array size\n");
        return profile;
    }
    else if (width%2==0){
        printf("Error: Window width must be an odd integer\n");
        return profile;
    }
    else {
        double xx[width];
        double *smoothed=new double[nbin];
        int nm=width/2-1;
        for (int i=0; i<nbin; i++) {
            int tn=0;
            for (int j=i-nm; j<i+nm+1; j++) {
                if (j>-1 && j<nbin) {
                    xx[tn]=profile[j];
                    tn++;
                }
            }
            smoothed[i]=TMath::Median(tn,xx);
        }
        //Fix the ends
        double Y0=3.*profile[0]-2.*profile[1];
        xx[0]=Y0;
        xx[1]=profile[0];
        xx[2]=profile[1];
        smoothed[0]=TMath::Median(3,xx);
        Y0=3.*profile[nbin-1]-2.*profile[nbin-2];
        xx[0]=Y0;
        xx[1]=profile[nbin-1];
        xx[2]=profile[nbin-2];
        smoothed[nbin-1]=TMath::Median(3,xx);
        return smoothed;
    }
}

void deproject(double z, double* profile, double *eprof,int nbin, double *bins, double *ebins, double *deprof, double *edeprof){
	double dist=angdist(z)*Mpc; //cm
	double kpc3=pow(Mpc/1e3,3.0);
	double b=(TMath::Pi()*dist/10800.)*(TMath::Pi()*dist/10800.);
	double *binsh=new double[nbin+1];
	double *invpr=new double[nbin];
	double *einvpr=new double[nbin];
	invertp(profile,invpr,nbin);
	invertp(eprof,einvpr,nbin);
	binsh[0]=0.1;
	for (int i=0;i<nbin;i++){
		double bh=bins[i]+ebins[i];
		binsh[i+1]=bh/60.*TMath::Pi()/180.*dist; //radius in cm
	}
	double *radius=new double[nbin+1];
	invertp(binsh,radius,nbin+1);
	
	double *indpr=new double[nbin];
	double *eindpr=new double[nbin];
    //edge correction
    double Rm=binsh[nbin];
    double Rm1=binsh[nbin-1];
    double rin=Rm1;
    double rout=Rm;
    double term1=1.0;
    double term2=rout/rin*TMath::ACos(rin/Rm);
    double term3=rout/Rm*sqrt(1-rin*rin/Rm/Rm);
    double term4=rout/rin-1.;
    double f=term1*(1.-2./TMath::Pi()*(term2-term3)/term4);    
	double anb=TMath::Pi()*(radius[0]*radius[0]-radius[1]*radius[1]);
    double Nout=invpr[0]*anb/b/(vij(radius,0,1)-vij(radius,1,1)-vij(radius,0,0)+vij(radius,1,0));
	indpr[0]=Nout*(1.-f);
	eindpr[0]=einvpr[0]*(anb/b/(vij(radius,0,1)-vij(radius,1,1)-vij(radius,0,0)+vij(radius,1,0)))*(1.-f);
	for (int m=1; m<nbin; m++) {
		anb=TMath::Pi()*(radius[m]*radius[m]-radius[m+1]*radius[m+1]);
		double sum=0.0;
		double es2=0.0;
		for (int i=0; i<m; i++) {
            double vol=vij(radius,i,m+1)-vij(radius,i+1,m+1)-vij(radius,i,m)+vij(radius,i+1,m);
			sum+=indpr[i]*vol;
			es2+=eindpr[i]*vol*eindpr[i]*vol;
		}
        rin=radius[m+1];
        rout=radius[m];
        term1=(Rm1+Rm)*Rm*Rm1/((rin+rout)*rin*rout);
        term2=(rout/rin*TMath::ACos(rin/Rm)-TMath::ACos(rout/Rm));
        term3=rout/Rm*(sqrt(1-rin*rin/Rm/Rm)-sqrt(1-rout*rout/Rm/Rm));
        term4=rout/rin-1.;
        f=term1*(1.-2./TMath::Pi()*(term2-term3)/term4);        
        if (m==nbin-1) {
            f=0.0;
        }
		indpr[m]=(invpr[m]*anb/b-sum)/(vij(radius,m,m+1)-vij(radius,m+1,m+1)-vij(radius,m,m)+vij(radius,m+1,m))-f*Nout;
		eindpr[m]=sqrt((einvpr[m]*anb/b)*(einvpr[m]*anb/b)+es2)/(vij(radius,m,m+1)-vij(radius,m+1,m+1)-vij(radius,m,m)+vij(radius,m+1,m));
	}
	for (int i=0; i<nbin; i++) {
		indpr[i]*=kpc3;
		eindpr[i]*=kpc3;
		//printf("rad, invpr, einvpr, indpr, eindpr: %g %g %g %g %g\n",radius[i]/Mpc*1e3,invpr[i],einvpr[i],indpr[i],eindpr[i]);
	}
	invertp(indpr,deprof,nbin);
	invertp(eindpr,edeprof,nbin);
    
    
	delete [] invpr;
	delete [] einvpr;
	delete [] indpr;
	delete [] radius;
	delete [] eindpr;
	delete [] binsh;
}

void mcdeproj(double z,double *profile,double *eprof,int nbin,double *bins,double *ebins,double *deprof,double *edeprof){
	int nsim=1e4;
	TRandom3 *g=new TRandom3(0);
	double **dv=new double*[nbin];
	for (int i=0; i<nbin; i++) {
		dv[i]=new double[nsim];
	}
    double *psm=medsmooth(nbin,profile,5);
	for (int i=0;i<nsim;i++){
		//create randomized profile
		double *npr=new double[nbin];
		double *enpr=new double[nbin];
		for (int nb=0;nb<nbin;nb++){
			npr[nb]=g->Gaus(psm[nb],eprof[nb]);
			enpr[nb]=eprof[nb];
		}
		
		//deproject
		double *deptmp=new double[nbin];
		double *edeptmp=new double[nbin];
		deproject(z,npr,enpr,nbin,bins,ebins,deptmp,edeptmp);
		double *depsm=medsmooth(nbin,deptmp,5);
		for (int nb=0; nb<nbin; nb++) {
			dv[nb][i]=depsm[nb];
		}
		delete [] npr;
		delete [] enpr;
		delete [] deptmp;
		delete [] edeptmp;
        delete [] depsm;
	}
	for (int i=0; i<nbin; i++) {
		double val=TMath::Mean(nsim,dv[i]);
		double eval=TMath::RMS(nsim,dv[i]);
		deprof[i]=val;
		edeprof[i]=eval;
	}
    //deprof=medsmooth(nbin,psm,5);
	for (int i=0;i<nbin;i++){
		delete [] dv[i];
	}
	delete [] dv;
    delete [] psm;
	delete g;
}

	
void logbinning(double binsize,double maxrad,int nbin,double *bins,double *ebins,int &newnbin){
	do {
		newnbin=0;
		double binedge=binsize/60./2.;
		double db=0.0;
		int i=0;
		while (db<binsize/60.) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;
			double bb=nbin/log10(maxrad/binsize*60.)*(log10(binedge)-log10(binsize/60.));
			//double base=log10(binsize/60.)+log10(maxrad/binsize*60.)*(bb-0.5)/nbin;
			double base1=log10(binsize/60.)+log10(maxrad/binsize*60.)*(bb-1)/nbin;
			double base2=log10(binsize/60.)+log10(maxrad/binsize*60.)*bb/nbin;
			double ebe=TMath::Power(10.,base2)-TMath::Power(10.,base1);
			db=ebe;
			binedge=bins[i];
			i++;
			if (i>nbin) {
				break;
			}
		}
		binedge+=binsize/60./2.;
		double thisbin=binedge;
		int b2=1;
		while (thisbin<maxrad) {
			double base1=log10(binedge)+log10(maxrad/binsize*60.)*(b2-1)/nbin;
			double base2=log10(binedge)+log10(maxrad/binsize*60.)*b2/nbin;
			thisbin=1./2.*(TMath::Power(10.,base1)+TMath::Power(10.,base2));
			double ebe=1./2.*(TMath::Power(10.,base2)-TMath::Power(10.,base1));
			bins[i]=thisbin;
			ebins[i]=ebe;
			b2++;
			i++;
			if (i>nbin) {
				break;
			}
		}
		newnbin=i-1;
	}while (0);
}

void sbpeak(double *img,double *exposure,long *axes,double &centroid_ra,double &centroid_dec){
    double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
    double maximg=-1e5;
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            if (exposure[j*axes[0]+i]>0.25*maxexp){
                double expocor=img[j*axes[0]+i]*maxexp/exposure[j*axes[0]+i];
                if (expocor>maximg){
                    centroid_ra=i;
                    centroid_dec=j;
                    maximg=expocor;
                }
            }
        }
    }
}

void centroid(double *img,double *exposure,long *axes,double rmax,double pixsize,double &centroid_ra,double &centroid_dec){
    double pra,pdec;
    sbpeak(img,exposure,axes,pra,pdec);
    double ftot=0.0;
    centroid_ra=0.0;
    centroid_dec=0.0;
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double dist=sqrt((i-pra)*(i-pra)+(j-pdec)*(j-pdec))*pixsize*60.;//arcmin
            if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
                centroid_ra+=img[j*axes[0]+i]/exposure[j*axes[0]+i]*i;
                centroid_dec+=img[j*axes[0]+i]/exposure[j*axes[0]+i]*j;
                ftot+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
            }
        }
    }
    centroid_ra/=ftot;
    centroid_dec/=ftot;
}

double medianerr(int nsect,int nmc,double *allprofs,double *allerrs){
    TRandom3 *g=new TRandom3(0);
    double *vals=new double[nmc];
    double *tempmed=new double[nsect];
    for (int i=0; i<nmc; i++) {
        for (int ns=0; ns<nsect; ns++) {
            tempmed[ns]=g->Gaus(allprofs[ns],allerrs[ns]);
        }
        vals[i]=TMath::Median(nsect,tempmed);
    }
    double err=TMath::RMS(nmc,vals);
    delete [] vals;
    delete [] tempmed;
    return err;
}

void centroid_guess(double *img,double *exposure,long *axes,double gra,double gdec,double rmax,double pixsize,double &centroid_ra,double &centroid_dec){
    double ftot=0.0;
    centroid_ra=0.0;
    centroid_dec=0.0;
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double dist=sqrt((i-gra)*(i-gra)+(j-gdec)*(j-gdec))*pixsize*60.;//arcmin
            if (exposure[j*axes[0]+i]>0.0 && dist<rmax){
                centroid_ra+=img[j*axes[0]+i]/exposure[j*axes[0]+i]*i;
                centroid_dec+=img[j*axes[0]+i]/exposure[j*axes[0]+i]*j;
                ftot+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
            }
        }
    }
    centroid_ra/=ftot;
    centroid_dec/=ftot;
}

double centroid_shift(double *img,double *exposure,double centroid_ra,double centroid_dec,long *axes,double pixsize,double rmax,int niter){
    double ramax,decmax;
    centroid_guess(img,exposure,axes,centroid_ra,centroid_dec,rmax,pixsize,ramax,decmax);
    double w=0.0;
    for (int it=1; it<niter; it++) {
        double trad=rmax/(niter-1)*it;
        double cra,cdec;
        centroid_guess(img,exposure,axes,centroid_ra,centroid_dec,trad,pixsize,cra,cdec);
        w+=((cra-ramax)*(cra-ramax)+(cdec-decmax)*(cdec-decmax))*pixsize*pixsize*60.*60.;
    }
    double wn=sqrt(w/(niter-1))/rmax;
    return wn;
}

double *mkbinsh(int nbin,double *bins,double *ebins,bool isr200,double r200obs,bool isr500,double r500obs,bool iskpc,double kpcobs){
    double *binsh=new double[nbin+1];
    binsh[0]=1e-2;
    for (int i=0;i<nbin;i++){
        if (isr500) {
            binsh[0]=1e-2/r500obs;
            binsh[i+1]=(bins[i]+ebins[i])/r500obs;
        }
        else if (isr200){
            binsh[0]=1e-2/r200obs;
            binsh[i+1]=(bins[i]+ebins[i])/r200obs;
        }
        else if (iskpc){
            binsh[0]=1e-2/kpcobs;
            binsh[i+1]=(bins[i]+ebins[i])/kpcobs;
        }
        else {
            binsh[i+1]=bins[i]+ebins[i];
        }
    }
    return binsh;
}

void region(char *fnam,double *exposure,long *axes){
    int nlin=line_num(fnam);
    if (nlin==0) {
        printf("    Unable to read file %s\n",fnam);
    }
    else {
        FILE *ff=fopen(fnam,"r");
        double *xpos=new double[nlin-3];
        double *ypos=new double[nlin-3];
        double *rad=new double[nlin-3];
        double *r2=new double[nlin-3];
        double *angles=new double[nlin-3];
        char c;
        // Go to the 3rd line
        int nn=0;
        while (nn<2){
            c = fgetc(ff);
            if(c == '\n') nn++;
        }
        char temp[200];
        fscanf(ff,"%s\n",temp);
        if (!strcmp(temp,"physical")||!strcmp(temp,"wcs")){
            printf("    Region file must be in image coordinates\n");
        }
        else if (strcmp(temp,"image")){
            printf("    Invalid region file %s\n",temp);
        }
        else {
            printf("    %d regions will be ignored\n",nlin-3);
            char type[100];
            for (int i=0;i<nlin-3;i++){
                fscanf(ff,"%[^\n]\n",type);
                char *pch=strtok(type,"-(,)");
                if (!strcmp(pch,"circle")) {
                    pch = strtok (NULL, "(,)");
                    xpos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    ypos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    rad[i]=atof(pch);
                    printf("Excluding region: circle(%g,%g,%g)\n",xpos[i]+1,ypos[i]+1,rad[i]);
                    r2[i]=rad[i];
                    angles[i]=0.0;
                }
                else if (!strcmp(pch,"ellipse")) {
                    pch = strtok (NULL, "(,)");
                    xpos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    ypos[i]=atof(pch)-1;
                    pch = strtok (NULL, "(,)");
                    rad[i]=atof(pch);
                    pch = strtok (NULL, "(,)");
                    r2[i]=atof(pch);
                    pch = strtok (NULL, "(,)");
                    angles[i]=atof(pch)*TMath::Pi()/180.+TMath::Pi()/2.;
                    printf("Excluding region: ellipse(%g,%g,%g,%g,%g)\n",xpos[i]+1,ypos[i]+1,rad[i],r2[i],atof(pch));
                }
            }
            fclose(ff);
            // Modified 01-24-2017 ; improved performance
            for (int ns=0;ns<nlin-3;ns++){
                int boxsize=(int)round(rad[ns]+0.5);
                if (r2[ns]>rad[ns]) {
                    boxsize=(int)round(r2[ns]+0.5);
                }
                int cx=(int)round(xpos[ns]);
                int cy=(int)round(ypos[ns]);
                for (int i=cx-boxsize; i<cx+boxsize+1; i++) {
                    for (int j=cy-boxsize; j<cy+boxsize+1; j++) {
                        if (i>=0 && j>=0 && i<axes[0] && j<axes[1]) {
                            double posx=(i-xpos[ns]);
                            double posy=(j-ypos[ns]);
                            double xtil=cos(angles[ns])*posx+sin(angles[ns])*posy;
                            double ytil=-sin(angles[ns])*posx+cos(angles[ns])*posy;
                            double aoverb=rad[ns]/r2[ns];
                            double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
                            if (dist<rad[ns]){
                                exposure[j*axes[0]+i]=0.0;
                            }
                        }
                    }
                }
            }
            /*
            for (int i=0;i<axes[0];i++){
                for (int j=0;j<axes[1];j++){
                    for (int ns=0;ns<nlin-3;ns++){
                        double posx=(i-xpos[ns]);
                        double posy=(j-ypos[ns]);
                        double xtil=cos(angles[ns])*posx+sin(angles[ns])*posy;
                        double ytil=-sin(angles[ns])*posx+cos(angles[ns])*posy;
                        double aoverb=rad[ns]/r2[ns];
                        double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
                        if (dist<rad[ns]){
                            exposure[j*axes[0]+i]=0.0;
                        }
                    }
                }
            }*/
            delete [] xpos;
            delete [] ypos;
            delete [] rad;
            delete [] r2;
            delete [] angles;
        }
    }
}

double csb(double z,char *modname,TF1 *model){
    int npp=model->GetNpar();
    TF1 *ftt=mkfint(modname,(char *)"ftemp");
    char *pars=new char[200];
    for (int k=0;k<npp;k++){
        pars=(char *)model->GetParName(k);
        if (!strcmp(pars,"const")){
            ftt->SetParameter(k,0.0);
            ftt->SetParName(k,"const");
            model->SetParName(k,"const");
        }
        else {
            double tpar=model->GetParameter(k);
            ftt->SetParameter(k,tpar);
        }
    }
    double dist=angdist(z)*1e3; //angular diameter distance in kpc
    double kpcobs=(1./dist)*180./TMath::Pi()*60.;//1 kpc in arcmin
    double kpc40=40.*kpcobs;
    double kpc400=400.*kpcobs;
    printf("    Integrating the flux between radii of 40 and 400 kpc (%g and %g arcmin)\n",kpc40,kpc400);
    double tcsb=ftt->Integral(0.,kpc40)/ftt->Integral(0.,kpc400);
    delete ftt;
    return tcsb;
}

double medianval(int bin,int nval,int nmc,double **allvals,double **allerrs){
    TRandom3 *g=new TRandom3(0);
    double *vals=new double[nmc];
    double *tempmed=new double[nval];
    for (int i=0; i<nmc; i++) {
        for (int ns=0; ns<nval; ns++) {
            tempmed[ns]=g->Gaus(allvals[bin][ns],allerrs[bin][ns]);
        }
        vals[i]=TMath::Median(nval,tempmed);
    }
    double err=TMath::RMS(nmc,vals);
    delete [] vals;
    delete [] tempmed;
    return err;
}

