/*
 *  sbprofs.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */

int getbinsb(double dist,double *bins,double *ebins,int nbin){
    int tbin=0;
    while (dist>bins[tbin]+ebins[tbin] && tbin<nbin){
        tbin++;
    }
    return tbin;
}


void mk_mod_img(TF1* fmod,double *exposure,double *outimg,long *axes,double centroid_ra,double centroid_dec,
				double pixsize, double ellang,double aoverb,bool ellipse){
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist;
			if (ellipse){
				double xtil=cos(ellang)*posx+sin(ellang)*posy;
				double ytil=-sin(ellang)*posx+cos(ellang)*posy;
				dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);		  
			}
			else {
				dist=sqrt(posx*posx+posy*posy);
			}
			outimg[j*axes[0]+i]=fmod->Eval(dist)*pixsize*60*pixsize*60*exposure[j*axes[0]+i];
		}
	}
}

void mk_scat_img(double *bins,int nbin,int nsect,double *profile,double *eprof,double **allprofs,double **allerrs,double *outimg,long *axes,double centroid_ra,double centroid_dec,double pixsize, double ellang,double aoverb,bool ellipse,double maxrad){
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
            double theta;
            if (posx==0.0 && posy==0.0)theta=0.0;
            else theta=atan(posy/posx);
            //Set all angles between 0 and 2pi
			if (posx<0.0){
				theta+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>=0.0)){
				theta+=2.*TMath::Pi();
			}
			double dist;
			if (ellipse){
				double xtil=cos(ellang)*posx+sin(ellang)*posy;
				double ytil=-sin(ellang)*posx+cos(ellang)*posy;
				dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);		  
			}
			else {
				dist=sqrt(posx*posx+posy*posy);
			}
            if (dist>maxrad) {
                outimg[j*axes[0]+i]=0.0;
            }
            else {
                int tbin=0;
                double mind=fabs(dist-bins[0]);
                for (int nb=0; nb<nbin; nb++) {
                    double td=fabs(dist-bins[nb]);
                    if (td<mind) {
                        mind=td;
                        tbin=nb;
                    }
                }
                int tsect=(int)floor(theta/2./TMath::Pi()*nsect);
                //if (tsect>=nsect || tsect<0)printf("tsect, theta, posx, posy: %d %g %g %g\n",tsect, theta,posx,posy);
                //if (tbin>=nbin || tbin<0)printf("tbin: %d\n",tbin); 
                if (allerrs[tsect][tbin]>0.0 && allprofs[tsect][tbin]!=0.0) {
                    double toterr=sqrt(allerrs[tsect][tbin]*allerrs[tsect][tbin]+eprof[tbin]*eprof[tbin]);
                    outimg[j*axes[0]+i]=(allprofs[tsect][tbin]-profile[tbin])/toterr;                  
                }
                else outimg[j*axes[0]+i]=0.0;
            }
		}
	}
}


void mk_counts_profile(double *img,double *exposure,double *profile,int *numbers,double *effexp,int &nbin,long *axes,
					   double centroid_ra,double centroid_dec,double pixsize,double binsize,double maxrad,bool islogbin){
	for (int i=0;i<nbin;i++){
		profile[i]=0.0;
		effexp[i]=0.0;
		numbers[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		double *bbt=new double[nbin];
		double *ebbt=new double[nbin];
		logbinning(binsize,maxrad,nbin,bbt,ebbt,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bbt[i-1]+ebbt[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		//nbin=newbin;
		delete [] tb;
		delete [] bbt;
		delete [] ebbt;
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					numbers[bin]++;
					effexp[bin]+=exposure[j*axes[0]+i];
				}
			}
		}
	}
	for (int i=0; i<nbin; i++) {
		effexp[i]/=numbers[i];
	}
	delete hbins;
}

void mk_modprof(double *img,double *exposure,double *bins,double *ebins,double *profile,int nbin,bool sector,double angh,double angl,
				long *axes,double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,double mincounts){
	for (int i=0;i<nbin;i++){
		profile[i]=0.0;
	}
	double anglow=0.0;
	double anghigh=2*TMath::Pi()+1e-10;
	if (sector){
		anghigh=angh;
		anglow=angl;
		if (angh<angl){//We cross the zero
			anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
		}
		else {
			anghigh-=anglow;
		}
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (mincounts>0.0){
					bin=getbinsb(dist,bins,ebins,nbin);
				}
				else {
					bin=(int)floor(dist/(binsize/60.)); //goes to the nearest integer
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
				}
			}
		}
	}  
}

void mk_sector_counts(double *img,double *exposure,double anglow,double anghigh,double *profile,int *numbers,double *effexp,int &nbin,long *axes,
					  double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin){
	for (int i=0;i<nbin;i++){
        profile[i]=0.0;
        numbers[i]=0;
        effexp[i]=0.0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		double *bbt=new double[nbin];
		double *ebbt=new double[nbin];
		logbinning(binsize,maxrad,nbin,bbt,ebbt,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bbt[i-1]+ebbt[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		//nbin=newbin;
		delete [] tb;
		delete [] bbt;
		delete [] ebbt;
	}
	double angh=anghigh;
	double angl=anglow;
	if (angh<angl){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					numbers[bin]++;
					effexp[bin]+=exposure[j*axes[0]+i];
				}
			}
		}
	} 
	for (int i=0; i<nbin; i++) {
		effexp[i]/=numbers[i];
	}
	delete hbins;
}

void mk_counts_box(double *img,double *exposure,double *profile,int *numbers,double *effexp,int &nbin,long *axes,
                       double centroid_ra,double centroid_dec,double pixsize,double binsize,double maxrad,double ellang,double width,bool islogbin){
    for (int i=0;i<nbin;i++){
        profile[i]=0.0;
        effexp[i]=0.0;
        numbers[i]=0;
    }
    TH1F *hbins=NULL;
    if (islogbin) {
        int newbin=0;
        double *bbt=new double[nbin];
        double *ebbt=new double[nbin];
        logbinning(binsize,maxrad,nbin,bbt,ebbt,newbin);
        double *tb=new double[newbin+1];
        tb[0]=0.0;
        for (int i=1; i<newbin+1; i++) {
            tb[i]=bbt[i-1]+ebbt[i-1];
        }
        hbins=new TH1F("ttt","ttt",newbin,tb);
        //nbin=newbin;
        delete [] tb;
        delete [] bbt;
        delete [] ebbt;
    }
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double posx=(i-centroid_ra)*pixsize*60;//arcmin
            double posy=(j-centroid_dec)*pixsize*60;
            double emod=ellang+TMath::Pi()/2.;
            double xtil=cos(emod)*posx+sin(emod)*posy;
            double ytil=-sin(emod)*posx+cos(emod)*posy;
            if ((xtil<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(fabs(ytil)<width/2.)){
                int bin;
                if (!islogbin) {
                    bin=(int)floor(xtil/(binsize/60.)); //left-inclusive
                }
                else {
                    bin=hbins->FindBin(xtil)-1;
                }
                if ((bin<nbin)&&(bin>=0)){
                    profile[bin]+=img[j*axes[0]+i];
                    numbers[bin]++;
                    effexp[bin]+=exposure[j*axes[0]+i];
                }
            }
        }
    }
    for (int i=0; i<nbin; i++) {
        effexp[i]/=numbers[i];
    }
    delete hbins;
}

void mk_sectell_counts(double *img,double *exposure,double anglow,double anghigh,double *profile,int *numbers,double *effexp,int &nbin,long *axes,
					   double centroid_ra,double centroid_dec,double pixsize,double ellang,double aoverb,double maxrad,double binsize,bool islogbin){
	for (int i=0;i<nbin;i++){
		profile[i]=0.0;
		numbers[i]=0;
		effexp[i]=0.0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		double *bbt=new double[nbin];
		double *ebbt=new double[nbin];
		logbinning(binsize,maxrad,nbin,bbt,ebbt,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bbt[i-1]+ebbt[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		//nbin=newbin;
		delete [] tb;
		delete [] bbt;
		delete [] ebbt;
	}
	double angh=anghigh;
	double angl=anglow;
	if (angh<angl){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
    double maxrad_pixel=maxrad/pixsize/60./aoverb;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double xtil=cos(ellang)*posx+sin(ellang)*posy;
			double ytil=-sin(ellang)*posx+cos(ellang)*posy;
			double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					numbers[bin]++;
					effexp[bin]+=exposure[j*axes[0]+i];
				}
			}
		}
	}  
	for (int i=0; i<nbin; i++) {
		effexp[i]/=numbers[i];
	}
	delete hbins;
}

void mk_ellipse_counts(double *img,double *exposure,double *profile,int *numbers,double *effexp,int &nbin,long *axes,double centroid_ra,
					   double centroid_dec,double pixsize,double ellang,double aoverb,double maxrad,double binsize,bool islogbin){
	for (int i=0;i<nbin;i++){
		profile[i]=0.0;
		numbers[i]=0;
		effexp[i]=0.0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		double *bbt=new double[nbin];
		double *ebbt=new double[nbin];
		logbinning(binsize,maxrad,nbin,bbt,ebbt,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bbt[i-1]+ebbt[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		//nbin=newbin;
		delete [] tb;
		delete [] bbt;
		delete [] ebbt;
	}
    double maxrad_pixel=maxrad/pixsize/60./aoverb;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double xtil=cos(ellang)*posx+sin(ellang)*posy;
			double ytil=-sin(ellang)*posx+cos(ellang)*posy;
			double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					numbers[bin]++;
					effexp[bin]+=exposure[j*axes[0]+i];
				}
			}
		}
	} 
	for (int i=0; i<nbin; i++) {
		effexp[i]/=numbers[i];
	}
	delete hbins;
}

void mk_sb_profile(double *img,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
				   double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
					eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			profile[i]/=pixsize*pixsize*60*60;
			eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}

void mk_sb_sig(double *img,double *sig,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
				   double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					eprof[bin]+=sig[j*axes[0]+i]*sig[j*axes[0]+i];
					nnn[bin]++;						
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]/=nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}


void mk_sector(double *img,double *exposure,double anglow,double anghigh,double *profile,double *eprof,double *bins,double *ebins,int &nbin,
			   long *axes,double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
	double angh=anghigh;
	double angl=anglow;
	if (angh<angl){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
    //printf("anglow, anghigh: %g , %g\n",anglow*180./TMath::Pi(),anghigh*180./TMath::Pi());
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
					eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			profile[i]/=pixsize*pixsize*60*60;
			eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}

void mk_sector_sig(double *img,double *sig,double *exposure,double anglow,double anghigh,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
			   double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
	double angh=anghigh;
	double angl=anglow;
	if (angh<angl){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i];
					eprof[bin]+=sig[j*axes[0]+i]*sig[j*axes[0]+i];
					nnn[bin]++;						
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			//profile[i]/=pixsize*pixsize*60*60;
			//eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}


void mk_growth_curve(double *img,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
					 double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		profile[i]=0.0;
		eprof[i]=0.0;
		nnn[i]=0;
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin=(int)floor(dist/(binsize/60.)); //left-inclusive
				if (bin<nbin){
					for (int bb=bin;bb<nbin;bb++){
						profile[bb]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
						eprof[bb]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					}
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (nnn[i]>0){
			eprof[i]=sqrt(eprof[i]);
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	for (int i=0;i<nbin;i++){
		bins[i]=(i+0.5)*binsize/60.;
		ebins[i]=binsize/60./2.;
	}
}

void mk_profile_ellipse(double *img,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
						double centroid_ra,double centroid_dec,double pixsize,double ellang,double aoverb,double maxrad,double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
    double maxrad_pixel=maxrad/pixsize/60./aoverb;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double xtil=cos(ellang)*posx+sin(ellang)*posy;
			double ytil=-sin(ellang)*posx+cos(ellang)*posy;
			double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
					eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			profile[i]/=pixsize*pixsize*60*60;
			eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}

void mk_ellipse_sector(double *img,double *exposure,double anglow,double anghigh,double *profile,double *eprof,double *bins,double *ebins,
                       int &nbin,long *axes,double centroid_ra,double centroid_dec,double pixsize,double ellang,double aoverb,double maxrad,
                       double binsize,bool islogbin){
    int nnn[nbin];
    for (int i=0;i<nbin;i++){
        eprof[i]=0.0;
        profile[i]=0.0;
        nnn[i]=0;
    }
    TH1F *hbins=NULL;
    if (islogbin) {
        int newbin=0;
        logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
        double *tb=new double[newbin+1];
        tb[0]=0.0;
        for (int i=1; i<newbin+1; i++) {
            tb[i]=bins[i-1]+ebins[i-1];
        }
        hbins=new TH1F("ttt","ttt",newbin,tb);
        nbin=newbin;
        delete [] tb;
    }
    
    // Calculate the opening angle of the sector.
    // At the end of this, ANGHIGH becomes the
    // opening angle of the sector; it is no longer
    // the ANGHIGH passed by user when defining the
    // elliptical sector.
    double angh=anghigh;
    double angl=anglow;
    if (angh < angl)
    {
        anghigh += 2*TMath::Pi() - anglow;
    }
    else {
        anghigh -= anglow;
    }
    
    double maxrad_pixel=maxrad/pixsize/60./aoverb;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
        for (int j=jmin;j<jmax+1;j++){
            double posx=(i-centroid_ra)*pixsize*60;//arcmin
            double posy=(j-centroid_dec)*pixsize*60;
            double xtil=cos(ellang)*posx+sin(ellang)*posy;
            double ytil=-sin(ellang)*posx+cos(ellang)*posy;
            double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
            
            double angle = atan2(posy, posx);
            if (angle<0)
                angle += 2*TMath::Pi();
            
            double low = anglow + ellang+0.5*TMath::Pi();
            double high = low + anghigh;
            
            if (low>2*TMath::Pi())
                low = fmod(low+8.*TMath::Pi(),2*TMath::Pi());
            if (high>2*TMath::Pi())
                high = fmod(high+8.*TMath::Pi(),2*TMath::Pi());
            
            if (high<low)
            {
                high += 2*TMath::Pi();
                if (angle < low)
                    angle += 2*TMath::Pi();
            }
            if((dist < maxrad) && (exposure[j*axes[0]+i]>0.0) && (angle>low) && (angle<high)) {
                int bin;
                //   printf("pixel %d, %d, val = %g\n",i+1,j+1,img[j*axes[0]+i]);
                if (!islogbin) {
                    bin=(int)floor(dist/(binsize/60.)); //left-inclusive
                }
                else {
                    bin=hbins->FindBin(dist)-1;
                }
                if (bin<nbin){
                    img[j*axes[0]+i] = TMath::Max(img[j*axes[0]+i], 0.);
                    profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
                    eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
                    nnn[bin]++;
                }
            }
        }
    }
    for (int i=0;i<nbin;i++){
        if (!islogbin) {
            bins[i]=(i+0.5)*binsize/60.;
            ebins[i]=binsize/60./2.;			
        }
        if (nnn[i]>0){
            profile[i]=profile[i]/nnn[i];
            eprof[i]=sqrt(eprof[i])/nnn[i];
            profile[i]/=pixsize*pixsize*60*60; // pix to arcmin**2
            eprof[i]/=pixsize*pixsize*60*60;
        }
        else {
            profile[i]=0.0;
            eprof[i]=0.0;
        }
        //   printf("sx profile, bin %d = %g\n",i,profile[i]);
    }
    delete hbins;
}

/*void mk_ellipse_sector(double *img,double *exposure,double anglow,double anghigh,double *profile,double *eprof,double *bins,double *ebins,
					   int &nbin,long *axes,double centroid_ra,double centroid_dec,double pixsize,double ellang,double aoverb,double maxrad,
					   double binsize,bool islogbin){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	TH1F *hbins=NULL;
	if (islogbin) {
		int newbin=0;
		logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
		double *tb=new double[newbin+1];
		tb[0]=0.0;
		for (int i=1; i<newbin+1; i++) {
			tb[i]=bins[i-1]+ebins[i-1];
		}
		hbins=new TH1F("ttt","ttt",newbin,tb);
		nbin=newbin;
		delete [] tb;
	}
	double angh=anghigh;
	double angl=anglow;
	if (anghigh<anglow){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
	for (int i=0;i<axes[0];i++){
		for (int j=0;j<axes[1];j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			double xtil=cos(ellang)*posx+sin(ellang)*posy;
			double ytil=-sin(ellang)*posx+cos(ellang)*posy;
			double dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin;
				if (!islogbin) {
					bin=(int)floor(dist/(binsize/60.)); //left-inclusive					
				}
				else {
					bin=hbins->FindBin(dist)-1;
				}
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
					eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (!islogbin) {
			bins[i]=(i+0.5)*binsize/60.;
			ebins[i]=binsize/60./2.;			
		}
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			profile[i]/=pixsize*pixsize*60*60;
			eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete hbins;
}*/

void mk_sb_box(double *img,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
                   double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,double ellang,double width,bool islogbin){
    int nnn[nbin];
    for (int i=0;i<nbin;i++){
        eprof[i]=0.0;
        profile[i]=0.0;
        nnn[i]=0;
    }
    TH1F *hbins=NULL;
    if (islogbin) {
        int newbin=0;
        logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
        double *tb=new double[newbin+1];
        tb[0]=0.0;
        for (int i=1; i<newbin+1; i++) {
            tb[i]=bins[i-1]+ebins[i-1];
        }
        hbins=new TH1F("ttt","ttt",newbin,tb);
        nbin=newbin;
        delete [] tb;
    }
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double posx=(i-centroid_ra)*pixsize*60;//arcmin
            double posy=(j-centroid_dec)*pixsize*60;
            double emod=ellang+TMath::Pi()/2.;
            double xtil=cos(emod)*posx+sin(emod)*posy;
            double ytil=-sin(emod)*posx+cos(emod)*posy;
            if ((xtil<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(fabs(ytil)<width/2.)){
                int bin;
                if (!islogbin) {
                    bin=(int)floor(xtil/(binsize/60.)); //left-inclusive
                }
                else {
                    bin=hbins->FindBin(xtil)-1;
                }
                if ((bin<nbin)&&(bin>=0)){
                    profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
                    eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
                    nnn[bin]++;
                }
            }
        }
    }
    for (int i=0;i<nbin;i++){
        if (!islogbin) {
            bins[i]=(i+0.5)*binsize/60.;
            ebins[i]=binsize/60./2.;			
        }
        if (nnn[i]>0){
            profile[i]=profile[i]/nnn[i];
            eprof[i]=sqrt(eprof[i])/nnn[i];
            profile[i]/=pixsize*pixsize*60*60;
            eprof[i]/=pixsize*pixsize*60*60;
        }
        else {
            profile[i]=0.0;
            eprof[i]=0.0;
        }
    }
    delete hbins;
}

void mk_median_prof(double *img,double *error,double *profile,double *eprof,double *bins,double *ebins,int nbin,
                    long *axes,double centroid_ra,double centroid_dec,double pixsize,double maxrad){
    int nnn[nbin];
    for (int i=0;i<nbin;i++){
        eprof[i]=0.0;
        profile[i]=0.0;
        nnn[i]=0;
    }
    double *tb=mkbinsh(nbin,bins,ebins,false,0.,false,0.,false,0.);
    TH1F *hbins=new TH1F("ttt","ttt",nbin,tb);
    double **allvals=new double*[nbin];
    double **allerrs=new double*[nbin];
    for (int i=0; i<nbin; i++) {
        allvals[i]=new double[axes[0]*axes[1]];
        allerrs[i]=new double[axes[0]*axes[1]];
    }
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double posx=(i-centroid_ra)*pixsize*60.;//arcmin
            double posy=(j-centroid_dec)*pixsize*60.;
            double dist=sqrt(posx*posx+posy*posy);
            if ((dist<maxrad)&&(error[j*axes[0]+i]>0.0)){
                int bin=hbins->FindBin(dist)-1;
                if (bin>=0 && bin<nbin){
                    int tn=nnn[bin];
                    if (bin<nbin){
                        allvals[bin][tn]=img[j*axes[0]+i];
                        allerrs[bin][tn]=error[j*axes[0]+i];
                        nnn[bin]++;
                    }
                }
            }
        }
    }
    for (int i=0;i<nbin;i++){
        if (nnn[i]>0){
            profile[i]=TMath::Median(nnn[i],allvals[i]);
            eprof[i]=medianval(i,nnn[i],1e3,allvals,allerrs);
        }
        else {
            profile[i]=0.0;
            eprof[i]=0.0;
        }
    }
    for (int i=0; i<nbin; i++) {
        delete [] allvals[i];
        delete [] allerrs[i];
    }
    delete [] allvals;
    delete [] tb;
    hbins->Delete();
    delete [] allerrs;
}

void mk_grouping(double *cprof,double *bins,double *ebins,double *profile,double *eprof,double *area,double *effexp,int &nbin,double mincounts){
	double *newprof=new double[nbin];
	double *newbins=new double[nbin];
	double *neweb=new double[nbin];
	double *newpr=new double[nbin];
	double *newepr=new double[nbin];
	double *newarea=new double[nbin];
	double *neweffexp=new double[nbin];
	int k=0;
	for (int i=0;i<nbin;i++){
		newprof[k]=cprof[i];
		newbins[k]=bins[i];
		neweb[k]=ebins[i];
		newpr[k]=profile[i];
		newepr[k]=eprof[i];
		newarea[k]=area[i];
		neweffexp[k]=effexp[i];
		if (newprof[k]>mincounts){
			k++;
		}
		else {
			if (i<nbin-1){
				int l=0;
				double tpr,tv;
				tpr=profile[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];//counts s-1
				tv=eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i]*eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				while ((newprof[k]<mincounts)&&(i+l<nbin-1)){
					l++;
					newprof[k]+=cprof[i+l];
					newbins[k]+=bins[i+l];
					neweb[k]+=ebins[i+l];
					newarea[k]+=area[i+l];
					neweffexp[k]+=effexp[i+l];
					tpr+=profile[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tv+=eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l]*eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
				}
				if (i+l<nbin){
					newbins[k]/=(l+1);
					neweffexp[k]/=(l+1);
				}
				newpr[k]=tpr/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);//counts s-1 arcmin-2
				newepr[k]=sqrt(tv)/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				i+=l;
				k++;
			}
		}
	}
	nbin=k;
	for (int i=0;i<nbin;i++){
		cprof[i]=newprof[i];
		bins[i]=newbins[i];
		ebins[i]=neweb[i];
		profile[i]=newpr[i];
		eprof[i]=newepr[i];
		area[i]=newarea[i];
		effexp[i]=neweffexp[i];
	}
	delete [] newprof;
	delete [] newbins;
	delete [] neweb;
	delete [] newepr;
	delete [] newpr;
	delete [] neweffexp;
	delete [] newarea;
}

void mk_group_isback(double *cprof,double *bins,double *ebins,double*profile,double*eprof,double *area,double *effexp,int &nbin,double mincounts,double *backprof,double *backcounts){
	double *newprof=new double[nbin];
	double *newbins=new double[nbin];
	double *neweb=new double[nbin];
	double *newpr=new double[nbin];
	double *newepr=new double[nbin];
	double *newback=new double[nbin];
	double *newarea=new double[nbin];
	double *neweffexp=new double[nbin];
	double *newbck=new double[nbin];
	int k=0;
	for (int i=0;i<nbin;i++){
		newprof[k]=cprof[i];
		newbins[k]=bins[i];
		neweb[k]=ebins[i];
		newpr[k]=profile[i];
		newepr[k]=eprof[i];
		newback[k]=backprof[i];
		newarea[k]=area[i];
		neweffexp[k]=effexp[i];
		newbck[k]=backcounts[i];
		if (newprof[k]>mincounts){
			k++;
		}
		else {
			if (i<nbin-1){
				int l=0;
				double tpr,tv,tb;
				tpr=profile[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];//counts s-1
				tv=eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i]*eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				tb=backprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				while ((newprof[k]<mincounts)&&(i+l<nbin-1)){
					l++;
					newprof[k]+=cprof[i+l];
					newbins[k]+=bins[i+l];
					neweb[k]+=ebins[i+l];
					newarea[k]+=area[i+l];
					neweffexp[k]+=effexp[i+l];
					newbck[k]+=backcounts[i+l];
					tpr+=profile[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tv+=eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l]*eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tb+=backprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
				}
				if (i+l<nbin){
					newbins[k]/=(l+1);
					neweffexp[k]/=(l+1);
				}
				/*else {
					newbins[k]/=l;
				}*/
				newpr[k]=tpr/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);//counts s-1 arcmin-2
				newepr[k]=sqrt(tv)/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				newback[k]=tb/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				i+=l;
				k++;
			}
			/*else {
				double tpr=newpr[k-1]*2.*TMath::Pi()*bins[k-1]*2.*ebins[k-1]+profile[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				double tv=newepr[k-1]*2.*TMath::Pi()*bins[k-1]*2.*ebins[k-1]*newepr[k-1]*2.*TMath::Pi()*bins[k-1]*2.*ebins[k-1]+eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i]*eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				double tb=newback[k-1]*2.*TMath::Pi()*bins[k-1]*2.*ebins[k-1]+backprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				newprof[k-1]+=cprof[i];
				newbins[k-1]+=bins[i];
				newbins[k-1]/=2.0;
				neweb[k-1]+=ebins[i];
				newpr[k-1]=tpr/(2.*TMath::Pi()*newbins[k-1]*2.*neweb[k-1]);
				newepr[k-1]=sqrt(tv)/(2.*TMath::Pi()*newbins[k-1]*2.*neweb[k-1]);
				newback[k-1]=tb/(2.*TMath::Pi()*newbins[k-1]*2.*neweb[k-1]);
			}*/
		}
	}
	nbin=k;
	for (int i=0;i<nbin;i++){
		cprof[i]=newprof[i];
		bins[i]=newbins[i];
		ebins[i]=neweb[i];
		profile[i]=newpr[i];
		eprof[i]=newepr[i];
		backprof[i]=newback[i];
		area[i]=newarea[i];
		effexp[i]=neweffexp[i];
		backcounts[i]=newbck[i];
	}
	delete [] newprof;
	delete [] newbins;
	delete [] neweb;
	delete [] newepr;
	delete [] newpr;
	delete [] newback;
	delete [] neweffexp;
	delete [] newarea;
	delete [] newbck;
}


void mk_group_counts(double *cprof,double *bins,double *ebins,double *area,double *effexp,int nbtot,int &nbin,double mincounts){
	double *newprof=new double[nbtot];
	double *newbins=new double[nbtot];
	double *neweb=new double[nbtot];
	double *newarea=new double[nbin];
	double *neweffexp=new double[nbin];
	int k=0;
	for (int i=0;i<nbtot;i++){
		newprof[k]=cprof[i];
		newbins[k]=bins[i];
		neweb[k]=ebins[i];
		newarea[k]=area[i];
		neweffexp[k]=effexp[i];
		if (newprof[k]>mincounts){
			k++;
		}
		else {
			if (i<nbtot-1){
				int l=0;
				while ((newprof[k]<mincounts)&&(i+l<nbtot-1)){
					l++;
					newprof[k]+=cprof[i+l];
					newbins[k]+=bins[i+l];
					neweb[k]+=ebins[i+l];
					newarea[k]+=area[i+l];
					neweffexp[k]+=effexp[i+l];
				}
				if (i+l<nbtot){
					newbins[k]/=(l+1);
					neweffexp[k]/=(l+1);
				}
				/*else {
					newbins[k]/=l;
				}*/
				i+=l;
				k++;
			}
			/*else {
				newprof[k-1]+=cprof[i];
				newbins[k-1]+=bins[i];
				newbins[k-1]/=2.0;
				neweb[k-1]+=ebins[i];
			}*/
		}
	}
	for (int i=0;i<nbin;i++){
		cprof[i]=newprof[i];
		bins[i]=newbins[i];
		ebins[i]=neweb[i];
		area[i]=newarea[i];
		effexp[i]=neweffexp[i];
	}
	delete [] newprof;
	delete [] newbins;
	delete [] neweb;
	delete [] neweffexp;
	delete [] newarea;
}

void mk_group_sn(double *cprof,double *bins,double *ebins,double*profile,double*eprof,double *area,double *effexp,int &nbin,double minsn){
	double *newprof=new double[nbin];
	double *newbins=new double[nbin];
	double *neweb=new double[nbin];
	double *newpr=new double[nbin];
	double *newepr=new double[nbin];
	double *newarea=new double[nbin];
	double *neweffexp=new double[nbin];
	int k=0;
	for (int i=0;i<nbin;i++){
		newprof[k]=cprof[i];
		newbins[k]=bins[i];
		neweb[k]=ebins[i];
		newpr[k]=profile[i];
		newepr[k]=eprof[i];
		newarea[k]=area[i];
		neweffexp[k]=effexp[i];
		double sn=newpr[k]/newepr[k];
		if (sn>=minsn){
			k++;
		}
		else {
			if (i<nbin-1){
				int l=0;
				double tpr,tv;
				tpr=profile[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];//counts s-1
				tv=eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i]*eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				while ((sn<minsn)&&(i+l<nbin-1)){
					l++;
					newprof[k]+=cprof[i+l];
					newbins[k]+=bins[i+l];
					neweb[k]+=ebins[i+l];
					newarea[k]+=area[i+l];
					neweffexp[k]+=effexp[i+l];
					tpr+=profile[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tv+=eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l]*eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					sn=tpr/sqrt(tv);
				}
				if (i+l<nbin){
					newbins[k]/=(l+1);
					neweffexp[k]/=(l+1);
				}
				/*else {
					newbins[k]/=l;
				}*/
				newpr[k]=tpr/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);//counts s-1 arcmin-2
				newepr[k]=sqrt(tv)/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				i+=l;
				k++;
			}
		}
	}
	nbin=k;
	for (int i=0;i<nbin;i++){
		cprof[i]=newprof[i];
		bins[i]=newbins[i];
		ebins[i]=neweb[i];
		profile[i]=newpr[i];
		eprof[i]=newepr[i];
		area[i]=newarea[i];
		effexp[i]=neweffexp[i];
	}
	delete [] newprof;
	delete [] newbins;
	delete [] neweb;
	delete [] newepr;
	delete [] newpr;
	delete [] neweffexp;
	delete [] newarea;
}

void mk_group_isback_sn(double *cprof,double *bins,double *ebins,double*profile,double*eprof,double *area,double *effexp,int &nbin,double minsn,double *backprof,double *backcounts){
	double *newprof=new double[nbin];
	double *newbins=new double[nbin];
	double *neweb=new double[nbin];
	double *newpr=new double[nbin];
	double *newepr=new double[nbin];
	double *newback=new double[nbin];
	double *newarea=new double[nbin];
	double *neweffexp=new double[nbin];
	double *newbck=new double[nbin];
	int k=0;
	for (int i=0;i<nbin;i++){
		newprof[k]=cprof[i];
		newbins[k]=bins[i];
		neweb[k]=ebins[i];
		newpr[k]=profile[i];
		newepr[k]=eprof[i];
		newback[k]=backprof[i];
		newarea[k]=area[i];
		neweffexp[k]=effexp[i];
		newbck[k]=backcounts[i];
		double sn=newpr[k]/newepr[k];
		if (sn>=minsn){
			k++;
		}
		else {
			if (i<nbin-1){
				int l=0;
				double tpr,tv,tb;
				tpr=profile[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];//counts s-1
				tv=eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i]*eprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				tb=backprof[i]*2.*TMath::Pi()*bins[i]*2.*ebins[i];
				while ((sn<minsn)&&(i+l<nbin-1)){
					l++;
					newprof[k]+=cprof[i+l];
					newbins[k]+=bins[i+l];
					neweb[k]+=ebins[i+l];
					newarea[k]+=area[i+l];
					neweffexp[k]+=effexp[i+l];
					newbck[k]+=backcounts[i+l];
					tpr+=profile[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tv+=eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l]*eprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					tb+=backprof[i+l]*2.*TMath::Pi()*bins[i+l]*2.*ebins[i+l];
					sn=tpr/sqrt(tv);
				}
				if (i+l<nbin){
					newbins[k]/=(l+1);
					neweffexp[k]/=(l+1);
				}
				/*else {
					newbins[k]/=l;
				}*/
				newpr[k]=tpr/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);//counts s-1 arcmin-2
				newepr[k]=sqrt(tv)/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				newback[k]=tb/(2.*TMath::Pi()*newbins[k]*2.*neweb[k]);
				i+=l;
				k++;
			}
		}
	}
	nbin=k;
	for (int i=0;i<nbin;i++){
		cprof[i]=newprof[i];
		bins[i]=newbins[i];
		ebins[i]=neweb[i];
		profile[i]=newpr[i];
		eprof[i]=newepr[i];
		backprof[i]=newback[i];
		area[i]=newarea[i];
		effexp[i]=neweffexp[i];
		backcounts[i]=newbck[i];
	}
	delete [] newprof;
	delete [] newbins;
	delete [] neweb;
	delete [] newepr;
	delete [] newpr;
	delete [] newback;
	delete [] neweffexp;
	delete [] newarea;
	delete [] newbck;
}


void mk_sector_scat(double *img,double *exposure,double anglow,double anghigh,double *bins,double *ebins,int nbin,
					double *profile,double *eprof,long *axes,double centroid_ra,double centroid_dec,double pixsize,double maxrad){
	int nnn[nbin];
	for (int i=0;i<nbin;i++){
		eprof[i]=0.0;
		profile[i]=0.0;
		nnn[i]=0;
	}
	double *binsh=new double[nbin+1];
	binsh[0]=0.0;
	for (int i=0;i<nbin;i++){
		binsh[i+1]=bins[i]+ebins[i];
	}
	TH1F* hh=new TH1F("hh","hh",nbin,binsh);
	double angh=anghigh;
	double angl=anglow;
	if (angh<angl){//We cross the zero
		anghigh+=2*TMath::Pi()-anglow;//angle with respect to anglow
	}
	else {
		anghigh-=anglow;
	}
    double maxrad_pixel=maxrad/pixsize/60.;
    int imin=(int)round(centroid_ra-maxrad_pixel-1.);
    int imax=(int)round(centroid_ra+maxrad_pixel+1.);
    int jmin=(int)round(centroid_dec-maxrad_pixel-1.);
    int jmax=(int)round(centroid_dec+maxrad_pixel+1.);
    if (imin<0) imin=0;
    if (jmin<0) jmin=0;
    if (imax>axes[0]-1) imax=axes[0]-1;
    if (jmax>axes[1]-1) jmax=axes[1]-1;
    for (int i=imin;i<imax+1;i++){
		for (int j=jmin;j<jmax+1;j++){
			double posx=(i-centroid_ra)*pixsize*60;//arcmin
			double posy=(j-centroid_dec)*pixsize*60;
			double angle=atan(posy/posx);
			//Set all angles between 0 and 2pi
			if (posx<0.0){
				angle+=TMath::Pi();
			}
			if ((posy<0.0)&&(posx>0.0)){
				angle+=2.*TMath::Pi();
			}
			if (angh<angl && angle<angl){// we cross the zero
				angle+=2*TMath::Pi()-anglow;
			}
			else angle-=anglow; //Set the origin at anglow
			double dist=sqrt(posx*posx+posy*posy);
			if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
				int bin=hh->FindBin(dist)-1; //left-inclusive
				if (bin<nbin){
					profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
					eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
					nnn[bin]++;
				}
			}
		}
	}  
	for (int i=0;i<nbin;i++){
		if (nnn[i]>0){
			profile[i]=profile[i]/nnn[i];
			eprof[i]=sqrt(eprof[i])/nnn[i];
			profile[i]/=pixsize*pixsize*60*60;
			eprof[i]/=pixsize*pixsize*60*60;
		}
		else {
			profile[i]=0.0;
			eprof[i]=0.0;
		}
	}
	delete [] binsh;
	delete hh;
}
