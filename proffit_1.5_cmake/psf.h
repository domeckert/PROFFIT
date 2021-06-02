/*
 *  psf.h
 *  
 *
 *  Created by Dominique Eckert on 05.10.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

//PSF convolution code

//Gaussian PSF
double intphi(double r0,double r,double sigpsf){
	int nb=100;
	double intot=0.0;
	double width=2*TMath::Pi()/nb;
	double normpsf=1./sqrt(2*TMath::Pi())/sigpsf;
	for (int i=0;i<nb;i++){
		double phi=i*2*TMath::Pi()/(nb-1);
		double d=r0*sqrt((r/r0-cos(phi))*(r/r0-cos(phi))+sin(phi)*sin(phi));
		double psfd=normpsf*exp(-d*d/2./sigpsf/sigpsf);
		intot+=width*psfd;
	}
	if (intot<1e-30)intot=0.0;
	return 1e3*intot;
}

double funcphi(double *x,double *par){
	double sigpsf=par[1];
	double r=par[0];
	double fc=intphi(x[0],r,sigpsf);
	return fc*x[0];
}

void psfconv(double *bins,double *ebins,int nbin,double sigpsf,double *psfmat){
	int np=100;
	double *xgl=new double[np];
	double *w=new double[np];
	for (int i=0;i<nbin;i++){
		double xc=bins[i];
		double sumprob=0.0;
		for (int j=0; j<nbin; j++) {
			double bl=bins[j]-ebins[j];
			double bh=bins[j]+ebins[j];
			TF1 *fphi=new TF1("fphi",funcphi,bl,bh,2);
			fphi->SetParameter(1,sigpsf/60.);
			fphi->SetParameter(0,xc);
			fphi->CalcGaussLegendreSamplingPoints(np,xgl,w,1e-16);
			double prob=fphi->IntegralFast(np,xgl,w,bl,bh);
			psfmat[i*nbin+j]=prob;
			sumprob+=prob;
			delete fphi;
		}
		for (int j=0; j<nbin; j++) {
			psfmat[i*nbin+j]/=sumprob;
		}
	}
	delete [] xgl;
	delete [] w;
}

//King profile
double intphiking(double r0,double r,double rc,double alpha){
	int nb=100;
	double intot=0.0;
	double width=2*TMath::Pi()/nb;
	for (int i=0;i<nb;i++){
		double phi=i*2*TMath::Pi()/(nb-1);
		double d=r0*sqrt((r/r0-cos(phi))*(r/r0-cos(phi))+sin(phi)*sin(phi));
		double base=1.+(d/rc)*(d/rc);
		double psfd=TMath::Power(base,-alpha);
		intot+=width*psfd;
	}
	if (intot<1e-30)intot=0.0;
	return 1e3*intot;
}

double funcphiking(double *x,double *par){
	double rc=par[1];
	double alpha=par[2];
	double r=par[0];
	double fc=intphiking(x[0],r,rc,alpha);
	return fc*x[0];
}

void psfconvking(double *bins,double *ebins,int nbin,double rc,double alpha,double *psfmat){
	int np=100;
	double *xgl=new double[np];
	double *w=new double[np];
	for (int i=0;i<nbin;i++){
		double xc=bins[i];
		double sumprob=0.0;
		for (int j=0; j<nbin; j++) {
			double bl=bins[j]-ebins[j];
			double bh=bins[j]+ebins[j];
			TF1 *fphi=new TF1("fphi",funcphiking,bl,bh,3);
			fphi->SetParameter(1,rc/60.);
			fphi->SetParameter(0,xc);
			fphi->SetParameter(2,alpha);
			fphi->CalcGaussLegendreSamplingPoints(np,xgl,w,1e-16);
			double prob=fphi->IntegralFast(np,xgl,w,bl,bh);
			psfmat[i*nbin+j]=prob;
			sumprob+=prob;
			delete fphi;
		}
		for (int j=0; j<nbin; j++) {
			psfmat[i*nbin+j]/=sumprob;
		}
	}
	delete [] xgl;
	delete [] w;
}

double calcking(double dist,double rc,double alpha){
    //double dist=sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]));
    double base=1.+(dist/rc)*(dist/rc);
    double king=TMath::Power(base,-alpha);
    return king;
}

// 2D Gaussian PSF model
double psfgauss2d(double *x,double *pars){
    double tx=x[0];
    double ty=x[1];
    double sigma=pars[0];
    double tr2=tx*tx+ty*ty;
    double bofr=1./sqrt(2*TMath::Pi())/sigma*exp(-tr2/2./sigma/sigma);
    return bofr;
}

TH2F* hpsf_lin(int nbh,double xhigh,TF2 *fpsf){
    if (nbh%2==1) {
        nbh=nbh+1;
    }
    TH2F *hist=new TH2F("hpsf","hpsf",nbh,-xhigh,xhigh,nbh,-xhigh,xhigh);
    //TH2F *hist=new TH2F("hpsf","hpsf",nbh,-xhigh,xhigh,nbh,-xhigh,xhigh);
    for (int i=0;i<nbh;i++){
        for (int j=0;j<nbh;j++){
            double tx=hist->GetXaxis()->GetBinCenter(i+1);
            double ty=hist->GetYaxis()->GetBinCenter(j+1);
            double tpsf=fpsf->Eval(tx,ty);
            hist->SetBinContent(i+1,j+1,tpsf);
        }
    }
    return hist;
}

//Ray-tracing-like code
void psfgaussnew(int nphot,double *bins,double *ebins,double *exposure,int nbin,long *axes,double centroid_ra,double centroid_dec,double sigpsf,double pixsize,double *psfmat){
    //if(gRandom) delete gRandom;
    //TRandom3 *gr = new TRandom3(0);
    gRandom->SetSeed(0);
    TH2F *photdist=new TH2F("pdist","pdist",axes[0],0,axes[0],axes[1],0,axes[1]);
    /*double *fbins=new double[1001];
    fbins[0]=0.0;
    for (int i=1; i<1e3+1; i++) {
        double logmax=log10(axes[0]);
        double expon=-1.+(logmax+1.)*i/1e3;
        fbins[i]=pow(10.,expon);
    }
    TH1F *fking=new TH1F("fk1","fk1",1e3,fbins);
    double sigpix=sigpsf/60./60./pixsize;
    for (int i=1; i<1e3; i++) {
        double vc=fking->GetBinCenter(i);
        double dr=2.*(vc-fking->GetBinLowEdge(i));
        double vk=-vc*vc/2./sigpix/sigpix;
        double tvk=1./sqrt(2.*TMath::Pi())/sigpix*exp(vk)*2.*TMath::Pi()*vc*dr;
        fking->SetBinContent(i,tvk);
    }
    delete [] fbins;*/
    TF2 *fpsf=new TF2("fpsf",psfgauss2d,-axes[0]/2.,axes[0]/2.,-axes[1]/2.,axes[1]/2.,1);
    double sigpix=sigpsf/60./60./pixsize;
    fpsf->SetParameter(0,sigpix);
    int nbh=200;
    TH2F *hpsf=hpsf_lin(nbh,axes[0]/2.,fpsf);
    double maxrad=bins[nbin-1]+ebins[nbin-1];
    for (int nb=0; nb<nbin; nb++) {
        for (int i=0; i<axes[0]*axes[1]; i++) {
            int ix=i%axes[0];
            int iy=i/axes[0];
            double dist=sqrt((ix-centroid_ra)*(ix-centroid_ra)+(iy-centroid_dec)*(iy-centroid_dec))*pixsize*60.;
            if (dist>=bins[nb]-ebins[nb] && dist<=bins[nb]+ebins[nb]) {
                photdist->SetBinContent(ix,iy,exposure[iy*axes[0]+ix]);
            }
            else {
                photdist->SetBinContent(ix,iy,0.0);
            }
        }
        double *rtrc=new double[nbin];
        for (int nbt=0; nbt<nbin; nbt++) {
            rtrc[nbt]=0.0;
        }
        int photx,photy;
        int npt=0;
        while (npt<nphot) {
            double px,py;
            photdist->GetRandom2(px,py);
            double xr,yr;
            /* double tr=fking->GetRandom();
            double theta=gr->Rndm()*2.*TMath::Pi();
            xr=tr*cos(theta);
            yr=tr*sin(theta);*/
            hpsf->GetRandom2(xr,yr);
            photx=(int)floor(px+xr+0.5);
            photy=(int)floor(py+yr+0.5);
            if (photx>=0 && photy>=0 && photx<axes[0] && photy<axes[1]) {
                if (exposure[photy*axes[0]+photx]>0.0) {
                    double dist=sqrt((photx-centroid_ra)*(photx-centroid_ra)+(photy-centroid_dec)*(photy-centroid_dec))*pixsize*60.;
                    if (dist<=maxrad) {
                        int tb=getbin(dist,bins,ebins,nbin);
                        rtrc[tb]+=1.0;
                        npt++;                   
                    }
                }
            }
        }
        for (int j=0;j<nbin;j++){
            double pij=rtrc[j]/nphot;
            psfmat[nb*nbin+j]=pij;
        }
        delete [] rtrc;
    }
    fpsf->Delete();
    //gr->Delete();
    hpsf->Delete();
    photdist->Delete();
}

void psfkingnew(int nphot,double *bins,double *ebins,double *exposure,int nbin,long *axes,double centroid_ra,double centroid_dec,double rc,double alpha,double pixsize,double *psfmat){
    //if(gRandom) delete gRandom;
    TRandom3 *gr = new TRandom3(0);
    gRandom->SetSeed(0);
    TH2F *photdist=new TH2F("pdist","pdist",axes[0],0,axes[0],axes[1],0,axes[1]);
    double *fbins=new double[1001];
    fbins[0]=0.0;
    for (int i=1; i<1e3+1; i++) {
        double logmax=log10(axes[0]);
        double expon=-1.+(logmax+1.)*i/1e3;
        fbins[i]=pow(10.,expon);
    }
    TH1F *fking=new TH1F("fk1","fk1",1e3,fbins);
    double trc=rc/60./60./pixsize;
    for (int i=1; i<1e3; i++) {
        double vc=fking->GetBinCenter(i);
        double dr=2.*(vc-fking->GetBinLowEdge(i));
        double vk=(1.+(vc/trc)*(vc/trc));
        double tvk=pow(vk,-alpha)*2.*TMath::Pi()*vc*dr;
        fking->SetBinContent(i,tvk);
    }
    delete [] fbins;
    double maxrad=bins[nbin-1]+ebins[nbin-1];
    for (int nb=0; nb<nbin; nb++) {
        for (int i=0; i<axes[0]*axes[1]; i++) {
            int ix=i%axes[0];
            int iy=i/axes[0];
            double dist=sqrt((ix-centroid_ra)*(ix-centroid_ra)+(iy-centroid_dec)*(iy-centroid_dec))*pixsize*60.;
            if (dist>=bins[nb]-ebins[nb] && dist<=bins[nb]+ebins[nb]) {
                photdist->SetBinContent(ix,iy,exposure[iy*axes[0]+ix]);
            }
            else {
                photdist->SetBinContent(ix,iy,0.0);
            }
        }
        double *rtrc=new double[nbin];
        for (int nbt=0; nbt<nbin; nbt++) {
            rtrc[nbt]=0.0;
        }
        int photx,photy;
        int npt=0;
        while (npt<nphot) {
            double px,py;
            photdist->GetRandom2(px,py);
            double tr,xr,yr;
            tr=fking->GetRandom();
            double theta=gr->Rndm()*2.*TMath::Pi();
            xr=tr*cos(theta);
            yr=tr*sin(theta);
            photx=(int)floor(px+xr+0.5);
            photy=(int)floor(py+yr+0.5);
            //printf("px,py,xr,yr: %g %g %g %g\n",px,py,xr,yr);
            if (photx>=0 && photy>=0 && photx<axes[0] && photy<axes[1]) {
                if (exposure[photy*axes[0]+photx]>0.0){
                    double dist=sqrt((photx-centroid_ra)*(photx-centroid_ra)+(photy-centroid_dec)*(photy-centroid_dec))*pixsize*60.;
                    if (dist<=maxrad) {
                        int tb=getbin(dist,bins,ebins,nbin);
                        rtrc[tb]+=1.0;
                        npt++;                   
                    }
                }
            }
        }
        //printf("hello %d 2\n",nb);
        for (int j=0;j<nbin;j++){
            double pij=rtrc[j]/nphot;
            psfmat[nb*nbin+j]=pij;
        }
        delete [] rtrc;
    }
    delete fking;
    delete gr;
    delete photdist;
}

void psfkingan(double *bins,double *area,int nbin,double rc,double alpha,double *psfmat){
    double trc=rc/60.;
    for (int i=0; i<nbin; i++) {
        double *vals=new double[nbin];
        double vtot=0.0;
        for (int j=0;j<nbin;j++){
            double dist=fabs(bins[i]-bins[j]);
            double vk=(1.+(dist/trc)*(dist/trc));
            double psfval=pow(vk,-alpha);
            //printf("i, j, psfval, vals: %d %d %g %g\n",i,j,psfval,psfval*area[j]);
            vals[j]=psfval*area[j];
            vtot+=vals[j];
        }
        for (int j=0;j<nbin;j++){
            psfmat[i*nbin+j]=vals[j]/vtot;
        }
        delete [] vals;
    }
}
