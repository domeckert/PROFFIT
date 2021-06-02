//Input-output functions

int getaxes(char *temp,long *axes){
	int status=0;
	char *sbmap=new char[200];
	sprintf(sbmap,"%s[0]",temp);
	fitsfile *x;
	fits_open_file(&x,sbmap,READONLY,&status);
	int bitpix,naxis;
	fits_get_img_param(x,2,&bitpix,&naxis,axes,&status);
	fits_close_file(x,&status);
	return status;
}

void initwcs(struct wcsprm *m_wcs, double* crpix,
             double* crval, double* cdelt,char **ctype)
{
    // Call wcsini() in WCSLIB.
    int naxis = 2;
    m_wcs->flag = -1;
    wcsini(1, naxis, m_wcs);
    
    for( int i=0; i<naxis; i++){
        strcpy(m_wcs->ctype[i], ctype[i] ); //axis type
        m_wcs->crval[i] = crval[i]; // reference value
        m_wcs->crpix[i] = crpix[i]; // pixel coordinate
        m_wcs->cdelt[i] = cdelt[i]; // scale factor
    }
    int status=0;
    if ((status = wcsset(m_wcs))) {
        fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
    }
}


int readimg(char *temp,double *img,long *axes,double &pixsize,double &cdelt1,double &crval1,double &crval2,
			 double &crpix1,double &crpix2,struct wcsprm *wcs_inp){
	int status=0;
	char *sbmap=new char[200];
	sprintf(sbmap,"%s",temp);
	fitsfile *x;
	fits_open_file(&x,sbmap,READONLY,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	int bitpix,naxis;
	fits_get_img_param(x,2,&bitpix,&naxis,axes,&status);
	long start[2]={1,1};
	int anynul;
	fits_read_pix(x,TDOUBLE,start,axes[0]*axes[1],NULL,img,&anynul,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
    char **ctype=new char*[2];
    ctype[0]=new char[20];
    ctype[1]=new char[20];
	fits_read_key(x,TDOUBLE,(char *)"CDELT2",&pixsize,NULL,&status);
    if (status!=0) {
        status=0;
        fits_read_key(x,TDOUBLE,(char *)"CD2_2",&pixsize,NULL,&status);
        if (status!=0) {
            printf("Warning: Unable to read CDELT parameter\n");
            pixsize=0.000694444;
            cdelt1=-0.000694444;
            status=0;
        }
        else {
            fits_read_key(x,TDOUBLE,(char *)"CD1_1",&cdelt1,NULL,&status);
        }
    }
    else {
        fits_read_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,NULL,&status);
    }
	fits_read_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRVAL1 could not be read. WCS conversions will not work properly.\n");
        crval1=0.0;
        status=0;
    }
	fits_read_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRVAL2 could not be read. WCS conversions will not work properly.\n");
        crval2=0.0;
        status=0;
    }
	fits_read_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRPIX1 could not be read. WCS conversions will not work properly.\n");
        crpix1=axes[0]/2.;
        status=0;
    }
	fits_read_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRPIX2 could not be read. WCS conversions will not work properly.\n");
        crpix2=axes[1]/2.;
        status=0;
    }
	fits_read_key(x,TSTRING,(char *)"CTYPE1",ctype[0],NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CTYPE1 could not be read. WCS conversions will not work properly.\n");
        sprintf(ctype[0],"RA---TAN");
        status=0;
    }
	fits_read_key(x,TSTRING,(char *)"CTYPE2",ctype[1],NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CTYPE2 could not be read. WCS conversions will not work properly.\n");
        sprintf(ctype[1],"DEC--TAN");
        status=0;
    }
	char telescope[200];
    int ts=0;
	fits_read_key(x,TSTRING,(char *)"TELESCOP",telescope,NULL,&ts);
	char instrument[200];
	fits_read_key(x,TSTRING,(char *)"INSTRUME",instrument,NULL,&ts);
    if (ts==0)	printf("    The instrument is %s/%s \n",telescope,instrument);
    double crpix[2]={crpix1,crpix2};
    double cdelt[2]={cdelt1,pixsize};
    double crval[2]={crval1,crval2};
    initwcs(wcs_inp,crpix,crval,cdelt,ctype);
    fits_close_file(x,&status);
	delete [] sbmap;
	return status;
}

int readexp(char *temp,long *axes,double *exposure){
	int status=0;
	char *expmap=new char[200];
	sprintf(expmap,"%s",temp);
	fitsfile *x;
	fits_open_file(&x,expmap,READONLY,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	long axexp[2];
	int bitpix,naxis;
	fits_get_img_param(x,2,&bitpix,&naxis,axexp,&status);
	if (axexp[0]!=axes[0] || axexp[1]!=axes[1]) {
		printf("    Error: image file and exposure map have different size\n");
		return status;
	}
	long start[2]={1,1};
	int anynul;
	fits_read_pix(x,TDOUBLE,start,axes[0]*axes[1],NULL,exposure,&anynul,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	fits_close_file(x,&status);
	delete [] expmap; 					
	return status;
}

int getaxes_expo(char *temp,long *axes,double &pixsize,double &cdelt1,double &crval1,double &crval2,
                 double &crpix1,double &crpix2,struct wcsprm *wcs_inp){
    int status=0;
    fitsfile *x;
    fits_open_file(&x,temp,READONLY,&status);
    int bitpix,naxis;
    fits_get_img_param(x,2,&bitpix,&naxis,axes,&status);
    char **ctype=new char*[2];
    ctype[0]=new char[20];
    ctype[1]=new char[20];
    fits_read_key(x,TDOUBLE,(char *)"CDELT2",&pixsize,NULL,&status);
    if (status!=0) {
        status=0;
        fits_read_key(x,TDOUBLE,(char *)"CD2_2",&pixsize,NULL,&status);
        if (status!=0) {
            printf("Warning: Unable to read CDELT parameter\n");
            pixsize=0.000694444;
            cdelt1=-0.000694444;
            status=0;
        }
        else {
            fits_read_key(x,TDOUBLE,(char *)"CD1_1",&cdelt1,NULL,&status);
        }
    }
    else {
        fits_read_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,NULL,&status);
    }
    fits_read_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRVAL1 could not be read. WCS conversions will not work properly.\n");
        crval1=0.0;
        status=0;
    }
    fits_read_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRVAL2 could not be read. WCS conversions will not work properly.\n");
        crval2=0.0;
        status=0;
    }
    fits_read_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRPIX1 could not be read. WCS conversions will not work properly.\n");
        crpix1=axes[0]/2.;
        status=0;
    }
    fits_read_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CRPIX2 could not be read. WCS conversions will not work properly.\n");
        crpix2=axes[1]/2.;
        status=0;
    }
    fits_read_key(x,TSTRING,(char *)"CTYPE1",ctype[0],NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CTYPE1 could not be read. WCS conversions will not work properly.\n");
        sprintf(ctype[0],"RA---TAN");
        status=0;
    }
    fits_read_key(x,TSTRING,(char *)"CTYPE2",ctype[1],NULL,&status);
    if (status!=0) {
        printf("Warning: Keyword CTYPE2 could not be read. WCS conversions will not work properly.\n");
        sprintf(ctype[1],"DEC--TAN");
        status=0;
    }
    char telescope[200];
    int ts=0;
    fits_read_key(x,TSTRING,(char *)"TELESCOP",telescope,NULL,&ts);
    char instrument[200];
    fits_read_key(x,TSTRING,(char *)"INSTRUME",instrument,NULL,&ts);
    if (ts==0)	printf("    The instrument is %s/%s \n",telescope,instrument);
    double crpix[2]={crpix1,crpix2};
    double cdelt[2]={cdelt1,pixsize};
    double crval[2]={crval1,crval2};
    initwcs(wcs_inp,crpix,crval,cdelt,ctype);
    fits_close_file(x,&status);
    return status;
}

int readback(char *temp,long *axes,double *backmap){
	int status=0;
	char *bmap=new char[200];
	sprintf(bmap,"%s",temp);
	fitsfile *x;
	fits_open_file(&x,bmap,READONLY,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	long axexp[2];
	int bitpix,naxis;
	fits_get_img_param(x,2,&bitpix,&naxis,axexp,&status);
	if (axexp[0]!=axes[0] || axexp[1]!=axes[1]) {
		printf("    Error: image file and background map have different size\n");
		return status;
	}
	long start[2]={1,1};
	int anynul;
	fits_read_pix(x,TDOUBLE,start,axes[0]*axes[1],NULL,backmap,&anynul,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	fits_close_file(x,&status);
	delete [] bmap;
	return status;
}

int save_fits_basic(fitsfile *x,int nbin,double *bins,double *ebins,double *profile,double *eprof,double *cprof,double *area,double *effexp){
	int status=0;
	char *ttype[7]={(char *)"RADIUS",(char *)"WIDTH",(char *)"SB",(char *)"ERR_SB",(char *)"COUNTS",(char *)"AREA",(char *)"EXPOSURE"};
	char *tform[7]={(char *)"1E",(char *)"1E",(char *)"1E",(char *)"1E",(char *)"1E",(char *)"1E",(char *)"1E"};
	char *tunit[7]={(char *)"arcmin",(char *)"arcmin",(char *)"counts s-1 arcmin-2",(char *)"counts s-1 arcmin-2",(char *)"counts",(char *)"arcmin2",(char *)"s"};
	int tfields=7;
	char extname[]="DATA";
	fits_create_tbl(x,BINARY_TBL,nbin,tfields,(char **)ttype,(char **)tform,(char **)tunit,extname,&status);
	fits_write_col(x,TDOUBLE,1,1,1,nbin,bins,&status);
	fits_write_col(x,TDOUBLE,2,1,1,nbin,ebins,&status);
	fits_write_col(x,TDOUBLE,3,1,1,nbin,profile,&status);
	fits_write_col(x,TDOUBLE,4,1,1,nbin,eprof,&status);
	fits_write_col(x,TDOUBLE,5,1,1,nbin,cprof,&status);
	fits_write_col(x,TDOUBLE,6,1,1,nbin,area,&status);				
	fits_write_col(x,TDOUBLE,7,1,1,nbin,effexp,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
	}
	return status;
}

int save_fits_isback(fitsfile *x,int nbin,double *backprof,double *backcounts,int &tfields){
	int status=0;
	char keyword[100];
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"BKG",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,backprof,&status);
	fits_make_keyn("TUNIT",tfields,keyword,&status);
	fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 arcmin-2",(char *)"physical unit of field",&status);
    
    tfields+=1;
    fits_insert_col(x,tfields,(char *)"BKGCOUNT",(char *)"1E",&status);
    fits_write_col(x,TDOUBLE,tfields,1,1,nbin,backcounts,&status);
    fits_make_keyn("TUNIT",tfields,keyword,&status);
    fits_write_key(x,TSTRING,keyword,(char *)"counts",(char *)"physical unit of field",&status);
	return status;
}

int save_fits_isdepr(fitsfile *x,int nbin,double *deprof,double *edeprof,int &tfields){
	int status=0;
	char keyword[100];
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"DEPR",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,deprof,&status);
	fits_make_keyn("TUNIT",tfields,keyword,&status);
	fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 kpc-3",(char *)"physical unit of field",&status);					
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"ERR_DEPR",(char *)"1E",&status);
	fits_make_keyn("TUNIT",tfields,keyword,&status);
	fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 kpc-3",(char *)"physical unit of field",&status);					
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,edeprof,&status);
	return status;
}

int save_fits_isdens(fitsfile *x,int nbin,double *dens,double *edens,int &tfields){
    int status=0;
    char keyword[100];
    tfields+=1;
    fits_insert_col(x,tfields,(char *)"DENSITY",(char *)"1E",&status);
    fits_write_col(x,TDOUBLE,tfields,1,1,nbin,dens,&status);
    fits_make_keyn("TUNIT",tfields,keyword,&status);
    fits_write_key(x,TSTRING,keyword,(char *)"cm-3",(char *)"physical unit of field",&status);
    tfields+=1;
    fits_insert_col(x,tfields,(char *)"ERR_DENS",(char *)"1E",&status);
    fits_make_keyn("TUNIT",tfields,keyword,&status);
    fits_write_key(x,TSTRING,keyword,(char *)"cm-3",(char *)"physical unit of field",&status);
    fits_write_col(x,TDOUBLE,tfields,1,1,nbin,edens,&status);
    return status;
}

int save_fits_ismed(fitsfile *x,int nbin,double *medprof,double *emedprof,int &tfields){
    int status=0;
    char keyword[100];
    tfields+=1;
    fits_insert_col(x,tfields,(char *)"MEDIANSB",(char *)"1E",&status);
    fits_write_col(x,TDOUBLE,tfields,1,1,nbin,medprof,&status);
    fits_make_keyn("TUNIT",tfields,keyword,&status);
    fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 arcmin-2",(char *)"physical unit of field",&status);
    tfields+=1;
    fits_insert_col(x,tfields,(char *)"ERR_MED",(char *)"1E",&status);
    fits_make_keyn("TUNIT",tfields,keyword,&status);
    fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 arcmin-2",(char *)"physical unit of field",&status);
    fits_write_col(x,TDOUBLE,tfields,1,1,nbin,emedprof,&status);
    return status;
}

int save_fits_ismod1(fitsfile *x,int nbin,double *bins,TF1 *model,double *profile,double *eprof,int &tfields){
	int status=0;
	char keyword[100];
	double *cmod=new double[nbin];
	double *deltachi=new double[nbin];
	for (int i=0; i<nbin; i++) {
		double rad=bins[i];
		cmod[i]=model->Eval(rad);
		deltachi[i]=(profile[i]-cmod[i])/eprof[i];
	}
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"MODEL",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,cmod,&status);
	fits_make_keyn("TUNIT",tfields,keyword,&status);
	fits_write_key(x,TSTRING,keyword,(char *)"counts s-1 arcmin-2",(char *)"physical unit of field",&status);					
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"DELTACHI",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,deltachi,&status);
	fits_make_keyn("TUNIT",tfields,keyword,&status);
	fits_write_key(x,TSTRING,keyword,(char *)"sigma",(char *)"physical unit of field",&status);					
	delete [] cmod;
	delete [] deltachi;
	return status;
}

int save_fits_isscat(fitsfile *x,int nbin,double *scat,double *escat,int &tfields){
	int status=0;
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"SCATTER",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,scat,&status);
	tfields+=1;
	fits_insert_col(x,tfields,(char *)"ERR_SCATTER",(char *)"1E",&status);
	fits_write_col(x,TDOUBLE,tfields,1,1,nbin,escat,&status);
	return status;
}

int save_fits_keys(fitsfile *x,long *axes,double *exposure,double centroid_ra,double centroid_dec,struct wcsprm *wcs_inp,bool ellipse,double ellang,double aoverb,bool sector,double angl,double angh,char *imgfile,char *expfile,char *backfile,bool isback){
	int status=0;
	double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
	fits_write_key(x,TDOUBLE,(char *)"EXPOSURE",&maxexp,(char *)"Total effective on-axis exposure",&status);
    double pixcrd[2],imgcrd[2],world[2];
    double phi,theta;
    int stat;
    pixcrd[0]=centroid_ra+1.;
    pixcrd[1]=centroid_dec+1.;
    status = wcsp2s(wcs_inp, 1, 2, pixcrd, imgcrd, &phi,&theta, world,&stat);
	fits_write_key(x,TDOUBLE,(char *)"RA_C",&world[0],(char *)"Right ascension of the centre in degrees",&status);
	fits_write_key(x,TDOUBLE,(char *)"DEC_C",&world[1],(char *)"Declination of the centre in degrees",&status);
	fits_write_key(x,TDOUBLE,(char *)"X_C",&pixcrd[0],(char *)"X axis fine position of the centre in pixels",&status);
	fits_write_key(x,TDOUBLE,(char *)"Y_C",&pixcrd[1],(char *)"Y axis fine position of the centre in pixels",&status);
    if (ellipse) {
        double angdeg=ellang*180./TMath::Pi();
        fits_write_key(x,TDOUBLE,(char *)"ELLANG",&angdeg,(char *)"Angle between RA axis and major axis (degrees)",&status);
        fits_write_key(x,TDOUBLE,(char *)"AOVERB",&aoverb,(char *)"Ratio between major and minor axis",&status);
    }
    if (sector) {
        double angldeg=angl*180./TMath::Pi();
        double anghdeg=angh*180./TMath::Pi();
        fits_write_key(x,TDOUBLE,(char *)"ANGLOW",&angldeg,(char *)"Lower angle of the sector (deg)",&status);
        fits_write_key(x,TDOUBLE,(char *)"ANGHIGH",&anghdeg,(char *)"Higher angle of the sector (deg)",&status);
    }
    fits_write_key(x,TSTRING,(char *)"IMGFILE",imgfile,(char *)"Link to image file",&status);
    fits_write_key(x,TSTRING,(char *)"EXPFILE",expfile,(char *)"Link to exposure map",&status);
    if (isback){
        fits_write_key(x,TSTRING,(char *)"BACKFILE",backfile,(char *)"Link to background map",&status);
    }
    return status;
}	

int save_fits_ismod2(fitsfile *x,TF1 *model,char *modname,char **names){
	int status=0;
	int npp=model->GetNpar();
	char *ttym[4]={(char *)"PAR",(char *)"NAME",(char *)"VALUE",(char *)"ERROR"};
	char *tfm[4]={(char *)"1I",(char *)"16A",(char *)"1E",(char *)"1E"};
	char extn[]="MODEL";
	fits_create_tbl(x,BINARY_TBL,npp,4,(char **)ttym,(char **)tfm,NULL,extn,&status);
	if (status!=0) {
		printf("    Error %d\n",status);
		return status;
	}
	int *parn=new int[npp];
	double *parval=new double[npp];
	double *parerr=new double[npp];
	for (int i=0; i<npp; i++) {
		parn[i]=i+1;
		parval[i]=model->GetParameter(i);
		parerr[i]=model->GetParError(i);
	}
	fits_write_col(x,TINT,1,1,1,npp,parn,&status);
	fits_write_col(x,TSTRING,2,1,1,npp,names,&status);
	fits_write_col(x,TDOUBLE,3,1,1,npp,parval,&status);
	fits_write_col(x,TDOUBLE,4,1,1,npp,parerr,&status);
    fits_write_key(x,TSTRING,(char *)"MODNAME",modname,(char *)"Model name",&status);
	delete [] parn;
	delete [] parval;
	delete [] parerr;
	return status;
}	

int save_fits_ispsf(fitsfile *x,int nbin,double *psfmat){
	int status=0;
	long nb2[2]={nbin,nbin};					
	fits_create_img(x,-32,2,nb2,&status);
	long start[2]={1,1};
	fits_write_pix(x,TDOUBLE,start,nbin*nbin,psfmat,&status);
	fits_write_key(x,TSTRING,(char *)"EXTNAME",(char *)"PSF",(char *)"PSF convolution matrix",&status);
	return status;
}

int save_img(char *temp,double *modimg,long *axes,double crpix1,double crval1,double crpix2,double crval2,double cdelt1,double pixsize){
	int status=0;
	fitsfile *x;
	char nfn[200];
	sprintf(nfn,"!%s",temp);
	fits_create_file(&x,nfn,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	fits_create_img(x,-32,2,axes,&status);
	long start[2]={1,1};
	fits_write_pix(x,TDOUBLE,start,axes[0]*axes[1],modimg,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	fits_write_key(x,TSTRING,(char *)"CREATOR",(void *)"proffit 1.5",NULL,&status);
	fits_write_key(x,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
	fits_write_key(x,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
	fits_write_key(x,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
	fits_write_key(x,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
	fits_write_key(x,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
	fits_write_key(x,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
	fits_write_key(x,TDOUBLE,(char *)"CDELT2",&pixsize,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
	fits_write_key(x,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
	fits_write_key(x,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
	fits_close_file(x,&status);
	if (status!=0){
		printf("    Error %d\n",status);
		return status;
	}
	else {
		printf("    Image succesfully written\n");
	}
	return status;
}	

int load_fits_structure(char *infile,int &nbin,long *axes,bool &isback,bool &ispsf,bool &ismod,char *modname){
    fitsfile *x;
    int status=0;
    fits_open_file(&x,infile,READONLY,&status);
    if (status!=0){
        printf("    Error %d\n",status);
        return status;
    }
    int nhdu,hdutype;
    char extname[200];
    status=fits_get_num_hdus(x,&nhdu,&status);
    bool isdata=false;
    for (int i=1;i<nhdu;i++){
        status=fits_movabs_hdu(x,i+1,&hdutype,&status);
        fits_read_key(x,TSTRING,(char *)"EXTNAME",extname,NULL,&status);
        if (!strcmp(extname,"DATA")){
            isdata=true;
            printf("    Found DATA structure\n");
            // Read image and exposure map keywords and load images
            char imgfile[200];
            fits_read_key(x,TSTRING,(char *)"IMGFILE",imgfile,NULL,&status);
            if (status!=0){
                printf("    Error: IMGFILE keyword could not be read\n");
                return status;
            }
            status=getaxes(imgfile,axes);
            char backfile[200];
            fits_read_key(x,TSTRING,(char *)"BACKFILE",backfile,NULL,&status);
            if (status==0){
                isback=true;
            }
            else {
                status=0;
            }
            long nsrctemp;
            fits_get_num_rows(x,&nsrctemp,&status);
            if (status!=0) {
                printf("    Error: DATA structure could not be read\n");
                printf("    Error %d\n",status);
                return status;
            }
            nbin=nsrctemp;
        }
        else if (!strcmp(extname,"MODEL")){
            fits_read_key(x,TSTRING,(char *)"MODNAME",modname,NULL,&status);
            if (status!=0){
                printf("    Error: Model name could not be read\n");
                return status;
            }
            ismod=true;
        }
        else if (!strcmp(extname,"PSF")){
            ispsf=true;
        }
    }
    fits_close_file(x,&status);
    if (!isdata){
        printf("    Error: No DATA structure could be found in input FITS file\n");
        return 1;
    }
    return 0;
}

int load_fits(char *infile,int nbin,long *axes,double *bins,double *ebins,double *profile,double *eprof,double *cprof,double *area,double *effexp,bool &isprofile,char *imgfile,char *expfile,double *img,bool &isimg,double *exposure,bool &isexp,double &centroid_ra,double &centroid_dec,double &pixsize,double &cdelt1,double &crval1,double &crval2,double &crpix1,double &crpix2,struct  wcsprm *wcs_inp,char *backfile,double *backmap,bool &isback,double *backprof,double *backcounts,char **names,TF1 *model,double *psfmat){
    fitsfile *x;
    int status=0;
    fits_open_file(&x,infile,READONLY,&status);
    if (status!=0){
        printf("    Error %d\n",status);
        return status;
    }
    int nhdu,hdutype;
    int colnum;
    int anynul;
    char extname[200];
    status=fits_get_num_hdus(x,&nhdu,&status);
    for (int i=1;i<nhdu;i++){
        status=fits_movabs_hdu(x,i+1,&hdutype,&status);
        fits_read_key(x,TSTRING,(char *)"EXTNAME",extname,NULL,&status);
        if (!strcmp(extname,"DATA")){
            //printf("    Found DATA structure\n");
            // Read image and exposure map keywords and load images
            fits_read_key(x,TSTRING,(char *)"IMGFILE",imgfile,NULL,&status);
            if (status!=0){
                printf("    Error: IMGFILE keyword could not be read\n");
                return status;
            }
            printf("Hello\n");
            status=readimg(imgfile,img,axes,pixsize,cdelt1,crval1,crval2,crpix1,crpix2,wcs_inp);
            if (status!=0){
                printf("    Error: image could not be read\n");
                printf("    Error %d\n",status);
                return status;
            }
            else {
                printf("    Image successfully loaded\n");
                isimg=true;
            }
            fits_read_key(x,TSTRING,(char *)"EXPFILE",expfile,NULL,&status);
            if (status!=0){
                printf("Error: EXPFILE keyword could not be read\n");
                return status;
            }
            status=readexp(expfile,axes,exposure);
            if (status!=0){
                printf("    Error: exposure map could not be read\n");
                printf("    Error %d\n",status);
                return status;
            }
            else {
                printf("    Exposure map successfully loaded\n");
                isexp=true;
            }
            fits_read_key(x,TSTRING,(char *)"BACKFILE",backfile,NULL,&status);
            if (status==0){
                status=readback(backfile,axes,backmap);
                if (status!=0){
                    printf("    Error: background map could not be read\n");
                    printf("    Error %d\n",status);
                    return status;
                }
                else {
                    printf("    Background map successfully loaded\n");
                    isback=true;
                }
            }
            else {
                isback=false;
                status=0;
            }
            double xc,yc;
            fits_read_key(x,TDOUBLE,(char *)"X_C",&xc,NULL,&status);
            if (status!=0){
                printf("    Error: Coordinates of center could not be read\n");
                return status;
            }
            centroid_ra=xc-1.;
            fits_read_key(x,TDOUBLE,(char *)"Y_C",&yc,NULL,&status);
            if (status!=0){
                printf("    Error: Coordinates of center could not be read\n");
                return status;
            }
            centroid_dec=yc-1.;
            fits_get_colnum(x,CASEINSEN,(char *)"RADIUS",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,bins,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"WIDTH",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,ebins,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"SB",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,profile,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"ERR_SB",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,eprof,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"COUNTS",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,cprof,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"AREA",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,area,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"EXPOSURE",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,effexp,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            isprofile=true;
            if (isback){
                fits_get_colnum(x,CASEINSEN,(char *)"BKG",&colnum,&status);
                if (status!=0) {
                    printf("Error %d\n",status);
                    return status;
                }
                fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,backprof,&anynul,&status);
                if (status!=0) {
                    printf("Error %d\n",status);
                    return status;
                }
                fits_get_colnum(x,CASEINSEN,(char *)"BKGCOUNT",&colnum,&status);
                if (status!=0) {
                    printf("Error %d\n",status);
                    return status;
                }
                fits_read_col(x,TDOUBLE,colnum,1,1,nbin,NULL,backcounts,&anynul,&status);
                if (status!=0) {
                    printf("Error %d\n",status);
                    return status;
                }
            }
        }
        if (!strcmp(extname,"MODEL")){
            printf("    Found MODEL structure\n");
            int npar=model->GetNpar();
            long npartemp;
            fits_get_num_rows(x,&npartemp,&status);
            if (status!=0) {
                printf("    Error: MODEL structure could not be read\n");
                printf("    Error %d\n",status);
                return status;
            }
            if (npar!=npartemp){
                printf("    Error: The number of model parameters in MODEL structure does not match the required number\n");
                return 1;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"NAME",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TSTRING,colnum,1,1,npar,NULL,names,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            double *vals=new double[npar];
            double *evals=new double[npar];
            fits_get_colnum(x,CASEINSEN,(char *)"VALUE",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,npar,NULL,vals,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_get_colnum(x,CASEINSEN,(char *)"ERROR",&colnum,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            fits_read_col(x,TDOUBLE,colnum,1,1,npar,NULL,evals,&anynul,&status);
            if (status!=0) {
                printf("Error %d\n",status);
                return status;
            }
            for (int i=0;i<npar;i++){
                model->SetParameter(i,vals[i]);
                model->SetParError(i,vals[i]);
                model->SetParName(i,names[i]);
            }
            delete [] vals;
            delete [] evals;
        }
        if (!strcmp(extname,"PSF")){
            printf("    Found PSF structure\n");
            int bitpix,naxis;
            fits_get_img_param(x,2,&bitpix,&naxis,axes,&status);
            if ((axes[0]!=nbin)||(axes[1]!=nbin)){
                printf("    Error: image size not equal to number of bins in profile\n");
                return 1;
            }
            long start[2]={1,1};
            int anynul;
            fits_read_pix(x,TDOUBLE,start,nbin*nbin,NULL,psfmat,&anynul,&status);
            if (status!=0) {
                printf("    Error %d\n",status);
                return status;
            }
        }
    }
    fits_close_file(x,&status);
    return status;
}
