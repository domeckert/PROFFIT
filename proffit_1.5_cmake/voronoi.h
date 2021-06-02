
const double Pi=3.14159265358979312;

double maxel(int npix,double *array){
    double tmax=array[0];
    for (int i=0; i<npix; i++) {
        if (array[i]>tmax) {
            tmax=array[i];
        }
    }
    return tmax;
}

void parametrization(int sbox,int n,int &px,int &py){
    int next=2*sbox+1;
    int nint=2*(sbox-1)+1;
    int npix=next*next-nint*nint;
    if (sbox<1 || n>=npix || n<0) {
        px=0;
        py=0;
    }
    if (n>=0 && n<next) {
        py=-sbox;
        px=n-sbox;
    }
    if (n>=next && n<next+nint) {
        px=sbox;
        py=n-next-(sbox-1);
    }
    if (n>=next+nint && n<2*next+nint) {
        py=sbox;
        px=n-next-nint-sbox;
    }
    if (n>=2*next+nint && n<npix) {
        px=-sbox;
        py=n-2*next-nint-(sbox-1);
    }
}

int getclosest(long *axes,int *binning, double cx, double cy){//find the next closest pixel
    int tpix=0;
    int ix=(int)floor(cx+0.5);
    int iy=(int)floor(cy+0.5);
    //printf("ix, iy: %d  %d\n",ix,iy);
    if (ix<0 || iy<0 || ix>axes[0]-1 || iy>axes[1]-1) {
        tpix=-1;
        return tpix;
    }
    if (binning[iy*axes[0]+ix]<1. && binning[iy*axes[0]+ix]>-1) {
        tpix=iy*axes[0]+ix;
    }
    else {
        int sbox=1;
        double dist=1e10;
        bool cont=true;
        int nit=0;
        while (cont) {
            int next=2*sbox+1;
            int nint=2*(sbox-1)+1;
            int npix=next*next-nint*nint;
            int nimg=0;
            for (int n=0;n<npix;n++){
                int px,py,tx,ty;
                parametrization(sbox,n,tx,ty);
                px=ix+tx;
                py=iy+ty;
                int pix=py*axes[0]+px;
                if (px>=0 && py>=0 && px<axes[0] && py<axes[1]) {
                    if (binning[pix]<1 && binning[pix]>-1) {
                        double td=(px-cx)*(px-cx)+(py-cy)*(py-cy);
                        if (td<dist) {
                            dist=td;
                            tpix=pix;
                        }
                        else if (td==dist){
                            double gr=rand()%1;
                            if (gr<0.5) {//switch half the time
                                dist=td;
                                tpix=pix;                                
                            }
                        }
                    }
                    nit++;
                    nimg++;
                }
            }
            if (nimg==0) {
                nit++;
            }
            if (dist<1e10) {
                cont=false;
            }
            if (nit>=axes[0]*axes[1]) {
                cont=false;
                tpix=-1;
            }
            else {
                sbox++;
            }
        }
        //printf("dist: %d\n",dist);
    }
    return tpix;
}

bool totok(int npix,int *binning){
    bool tok=true;
    int tpix=0;
    while (tok && tpix<npix) {
        if (binning[tpix]<1. && binning[tpix]>-1) {
            tok=false;
        }
        else {
            tpix++;
        }
    }
    return tok;
}


void getmean(int npixbin,long *axes,int *tpix,double &cra,double &cdec){//geometrical centroid
    cra=0.0;
    cdec=0.0;
    for (int i=0;i<npixbin;i++){
        int pixx=tpix[i]%axes[0];
        int pixy=tpix[i]/axes[0];
        cra+=pixx;
        cdec+=pixy;
    }
    if (npixbin>0) {
        cra/=npixbin;
        cdec/=npixbin;        
    }
}


double totcounts(int npixbin,int *tpix,double *img){
    int tc=0;
    for (int i=0; i<npixbin; i++) {
        tc+=img[tpix[i]];
    }
    return tc;
}

double totsnr(int npixbin,int *tpix,double *img,double *backimg){
    int tc=0;
    double tbkg=0;
    for (int i=0; i<npixbin; i++) {
        tc+=img[tpix[i]];
        tbkg+=backimg[tpix[i]];
    }
    double snr=tc/sqrt(tc+tbkg);
    return snr;
}


bool topologicalcrit(int pixel,int npixbin,int *tpix,long *axes){//topological criterion
    int pixx=pixel%axes[0];
    int pixy=pixel/axes[0];
    bool isadj=false;
    int i=0;
    while (!isadj && i<npixbin) {
        int ix=tpix[i]%axes[0];
        int iy=tpix[i]/axes[0];
        double dist=(ix-pixx)*(ix-pixx)+(iy-pixy)*(iy-pixy);
        if (dist<2.0) {//adjacent bin
            isadj=true;
        }
        i++;
    }
    return isadj;
}

bool morphologicalcrit(int pixel,int npixbin,int *tpix,long *axes){//morphological criterion
    int pixx=pixel%axes[0];
    int pixy=pixel/axes[0];
    //geometrical centroid
    double cra=pixx;
    double cdec=pixy;
    for (int i=0;i<npixbin;i++){
        int px=tpix[i]%axes[0];
        int py=tpix[i]/axes[0];
        cra+=px;
        cdec+=py;
    }
    cra/=npixbin+1;
    cdec/=npixbin+1;

    double rmax=(cra-pixx)*(cra-pixx)+(cdec-pixy)*(cdec-pixy);
    for (int i=0; i<npixbin; i++) {
        int ix=tpix[i]%axes[0];
        int iy=tpix[i]/axes[0];
        double dist=(ix-cra)*(ix-cra)+(iy-cdec)*(iy-cdec);
        if (dist>rmax) {
            rmax=dist;
        }
    }
    double reff=sqrt((npixbin+1.)/Pi);
    double R=sqrt(rmax)/reff-1.;
    //printf("rmax, reff, R: %g  %g  %g\n",sqrt(rmax),reff,R);
    bool isround=false;
    if (R<0.3) {//criterion to be "round"
        isround=true;
    }
    return isround;
}

bool uniformitycrit(int pixel,int npixbin,int *tpix,double *img,int mincounts){//uniformity criterion
    bool isuniform;
    double ncold=totcounts(npixbin,tpix,img);
    double ncnew=ncold+img[pixel];
    double devold=(ncold-mincounts)*(ncold-mincounts);
    double devnew=(ncnew-mincounts)*(ncnew-mincounts);
    if (devnew<=devold) {
        isuniform=true;
    }
    else isuniform=false;
    return isuniform;
}

bool uniformitycrit_snr(int pixel,int npixbin,int *tpix,double *img,double *backimg,double minsnr){//uniformity criterion
    bool isuniform;
    double snrold=totsnr(npixbin,tpix,img,backimg);
    int *tpp=new int[npixbin+1];
    for (int i=0; i<npixbin; i++) {
        tpp[i]=tpix[i];
    }
    tpp[npixbin]=pixel;
    double snrnew=totsnr(npixbin+1,tpp,img,backimg);
    double devold=(snrold-minsnr)*(snrold-minsnr);
    double devnew=(snrnew-minsnr)*(snrnew-minsnr);
    if (devnew<=devold) {
        isuniform=true;
    }
    else isuniform=false;
    delete [] tpp;
    return isuniform;
}


int findbin(int pixel,long *axes,int nbin,bool *isgood,double *cra,double *cdec){//determine which is the closest bins for unbinned pixels
    int pixx=pixel%axes[0];
    int pixy=pixel/axes[0];
    int bin=0;
    double dist=(cra[0]-pixx)*(cra[0]-pixx)+(cdec[0]-pixy)*(cdec[0]-pixy);
    for (int i=0; i<nbin; i++) {
        if (isgood[i]) {
            double td=(cra[i]-pixx)*(cra[i]-pixx)+(cdec[i]-pixy)*(cdec[i]-pixy);
            if (td<dist) {
                dist=td;
                bin=i;
            }            
        }
    }
    return bin;
}

int findbin_good(int pixel,long *axes,int nbin,double *scale,double *cra,double *cdec){//determine which is the closest bins for unbinned pixels
    int pixx=pixel%axes[0];
    int pixy=pixel/axes[0];
    int bin=0;
    double dist=((cra[0]-pixx)*(cra[0]-pixx)+(cdec[0]-pixy)*(cdec[0]-pixy))/scale[0]/scale[0];
    for (int i=0; i<nbin; i++) {
        double td=((cra[i]-pixx)*(cra[i]-pixx)+(cdec[i]-pixy)*(cdec[i]-pixy))/scale[i]/scale[i];
        if (td<dist) {
            dist=td;
            bin=i;
        }            
    }
    return bin;
}

void getcentroid(int npixbin,long *axes,int *tpix,double *img,double &cra,double &cdec){//img^2-weighted centroid
    cra=0.0;
    cdec=0.0;
    double ftot=0.0;
    for (int i=0;i<npixbin;i++){
        int pixx=tpix[i]%axes[0];
        int pixy=tpix[i]/axes[0];
        cra+=img[tpix[i]]*img[tpix[i]]*pixx;
        cdec+=img[tpix[i]]*img[tpix[i]]*pixy;
        ftot+=img[tpix[i]]*img[tpix[i]];
    }
    cra/=ftot;
    cdec/=ftot;
}


double compute_voronoi(long *axes,double *img,int *binning,int nbin,double *scale,bool wvt,double *cra,double *cdec,double *cpb){
    double *xnode=new double[nbin];
    double *ynode=new double[nbin];
    int *npixbin=new int[nbin];
    int **pix=new int*[nbin];
    for (int i=0; i<nbin; i++) {
        npixbin[i]=0;
        pix[i]=new int[1];
    }
    for (int i=0; i<axes[0]*axes[1]; i++) {
        if (binning[i]>-1) {
            int tbin=findbin_good(i,axes,nbin,scale,cra,cdec);
            binning[i]=tbin+1;
            int tnb=npixbin[tbin];
            int *transarr=new int[tnb+1];
            for (int j=0; j<tnb; j++) {
                transarr[j]=pix[tbin][j];
            }
            transarr[tnb]=i;
            delete [] pix[tbin];
            pix[tbin]=transarr;
            npixbin[tbin]++;
        }
    }
    //compute centroids weighted by ncounts^2
    for (int i=0; i<nbin; i++) {
        if (npixbin[i]>0) {
            getcentroid(npixbin[i],axes,pix[i],img,xnode[i],ynode[i]);
            cpb[i]=totcounts(npixbin[i],pix[i],img);
        }
        else {
            xnode[i]=cra[i];
            ynode[i]=cdec[i];
        }
    }
    if (wvt) {//Diehl method
        //adapt scale factors
        for (int i=0; i<nbin; i++) {
            if (npixbin[i]>0.0) {
                scale[i]=sqrt(npixbin[i]/cpb[i]);            
            }
            else {
                scale[i]=1.0;
            }
        }        
    }
    
    double totdiff=0.0;
    for (int i=0; i<nbin; i++) {
        totdiff+=(cra[i]-xnode[i])*(cra[i]-xnode[i])+(cdec[i]-ynode[i])*(cdec[i]-ynode[i]);
        cra[i]=xnode[i];
        cdec[i]=ynode[i];
    }
    delete [] npixbin;
    for (int i=0; i<nbin; i++) {
        delete [] pix[i];
    }
    delete [] pix;
    delete [] xnode;
    delete [] ynode;
    return totdiff;
}

double compute_voronoi_snr(long *axes,double *img,double *backimg,int *binning,int nbin,double *scale,bool wvt,double *cra,double *cdec,double *cpb){
    double *xnode=new double[nbin];
    double *ynode=new double[nbin];
    int *npixbin=new int[nbin];
    int **pix=new int*[nbin];
    for (int i=0; i<nbin; i++) {
        npixbin[i]=0;
        pix[i]=new int[1];
    }
    for (int i=0; i<axes[0]*axes[1]; i++) {
        if (binning[i]>-1) {
            int tbin=findbin_good(i,axes,nbin,scale,cra,cdec);
            binning[i]=tbin+1;
            int tnb=npixbin[tbin];
            int *transarr=new int[tnb+1];
            for (int j=0; j<tnb; j++) {
                transarr[j]=pix[tbin][j];
            }
            transarr[tnb]=i;
            delete [] pix[tbin];
            pix[tbin]=transarr;
            npixbin[tbin]++;
        }
    }
    //compute centroids weighted by ncounts^2
    for (int i=0; i<nbin; i++) {
        if (npixbin[i]>0) {
            getcentroid(npixbin[i],axes,pix[i],img,xnode[i],ynode[i]);
            cpb[i]=totsnr(npixbin[i],pix[i],img,backimg);
        }
        else {
            xnode[i]=cra[i];
            ynode[i]=cdec[i];
        }
    }
    if (wvt) {//Diehl method
        //adapt scale factors
        for (int i=0; i<nbin; i++) {
            if (npixbin[i]>0.0) {
                scale[i]=sqrt(npixbin[i]/cpb[i]);
            }
            else {
                scale[i]=1.0;
            }
        }
    }
    
    double totdiff=0.0;
    for (int i=0; i<nbin; i++) {
        totdiff+=(cra[i]-xnode[i])*(cra[i]-xnode[i])+(cdec[i]-ynode[i])*(cdec[i]-ynode[i]);
        cra[i]=xnode[i];
        cdec[i]=ynode[i];
    }
    delete [] npixbin;
    for (int i=0; i<nbin; i++) {
        delete [] pix[i];
    }
    delete [] pix;
    delete [] xnode;
    delete [] ynode;
    return totdiff;
}

void compute_sb(int npixbin,int *pix,double pixsize,double *img,double *expo,
                bool isexp,double *backimg,bool isback,double back,double eback,
                bool isbacksub,double &sb,double &esb){
    double profile=0.0;
    double vprof=0.0;
    double bp=0.0;
    for (int i=0; i<npixbin; i++) {
        int tpix=pix[i];
        double texp=1.0;
        if (isexp) {
            texp=expo[tpix];
        }
        profile+=img[tpix]/texp;
        vprof+=img[tpix]/texp/texp;
        if (isback) {
            bp+=backimg[tpix]/texp;
        }
    }
    if (npixbin>0) {
        sb=(profile-bp)/npixbin/pixsize/pixsize*60.*60.;
        //esb=sqrt(vprof)/npixbin/pixsize/pixsize*60.*60.;        
        esb=sqrt(vprof/npixbin)/pixsize/pixsize*60.*60.;        
        if (isbacksub) {
            sb-=back;
            double tesb=sqrt(esb*esb+eback*eback);//Error on the SB of a single pixel assuming the same error for all
            //the total error on the SB in the bin is this value divided by sqrt(npixbin)
            esb=tesb;
        }
    }
    else {
        sb=0.0;
        esb=0.0;
    }
}

void compute_outimg(int *binning,int nbin,long *axes,double *cra,double *cdec,
                    double pixsize,double *img,double *expo,bool isexp,double *backimg,
                    bool isback,double back,double eback,bool isbacksub,double *outimg,
                    double *error,long npix,double *scale,double *sb,double *esb,int *npixbin){
    
    int **pix=new int*[nbin];
    for (int i=0; i<nbin; i++) {
        npixbin[i]=0;
        pix[i]=new int[1];
    }
    for (int i=0; i<axes[0]*axes[1]; i++) {
        if (binning[i]>-1) {
            int tbin=findbin_good(i,axes,nbin,scale,cra,cdec);
            binning[i]=tbin;
            int tnb=npixbin[tbin];
            int *transarr=new int[tnb+1];
            for (int j=0; j<tnb; j++) {
                transarr[j]=pix[tbin][j];
            }
            transarr[tnb]=i;
            delete [] pix[tbin];
            pix[tbin]=transarr;
            npixbin[tbin]++;
        }
    }
    for (int i=0; i<nbin; i++) {
        compute_sb(npixbin[i],pix[i],pixsize,img,expo,isexp,backimg,isback,back,eback,isbacksub,sb[i],esb[i]);
    }
    for (int i=0; i<npix; i++) {
        if (binning[i]<0) {
            outimg[i]=0.0;
            error[i]=0.0;
        }
        else {
            int tbin=binning[i];
            outimg[i]=sb[tbin];
            error[i]=esb[tbin];
        }
    }
    for (int i=0; i<nbin; i++) {
        delete [] pix[i];
    }
    
}

void voronoi(double *img,long *axes,int mincounts,double pixsize,double crpix1,double crpix2,double crval1,double crval2,double cdelt1,double *expo,bool isback,double *backimg,char *outfile,double *outimg,double *error){
    do {
    pixsize*=3600.; //arcsec/
    bool isexp=true;
    long npix=axes[0]*axes[1];
    double cx,cy;//image pixels of pixel with largest number of counts
    double maximg=maxel(npix,img);
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            //if (img[j*axes[0]+i]==maximg){
            if (img[j*axes[0]+i]==maximg){
                if (isexp) {
                    if (expo[j*axes[0]+i]>0.0) {
                        cx=i;
                        cy=j;
                    }
                }
                else {
                    cx=i;
                    cy=j;
                }
            }
        }
    }
    int *binning=new int[npix];//mapping of pixels to bins
    bool *isbinned=new bool[npix];
    long npt=0;
    if (!isexp) {
        npt=npix;
    }
    for (int i=0; i<npix; i++) {
        isbinned[i]=false;
        if (isexp) {
            if (expo[i]>0.0) {
                binning[i]=0;
                npt++;
            }
            else {
                binning[i]=-1;
            }
        }
        else {
            binning[i]=0;
        }
    }
    printf("Number of pixels: %ld\n",npt);
    int cpix=cy*axes[0]+cx;
    //printf("Step 0: %g  %g\n",cx,cy);
    int **tpix=new int*[npix];//array to store pixels present in each bin
    int *npixbin=new int[npix];//number of pixels per bin
    int *binnedpix=new int[npix];//array to store all pixels already binned
    double *cra=new double[npix];//centroids of the bins
    double *cdec=new double[npix];
    bool *isgood=new bool[npix];
    long nok=0;
    long nbinned=0;
    int bin=0;
    int *tarray=new int[npix];
    printf("Creating initial condition...\n");
    while (nok<npt-2) {
        //printf("Working in bin %d...\n",bin);
        tarray[0]=cpix;
        npixbin[bin]=1;
        int nbin=bin+1;
        binning[cpix]=nbin;
        double ncb=totcounts(npixbin[bin],tarray,img);//number of counts in the current pixel
        getmean(npixbin[bin],axes,tarray,cra[bin],cdec[bin]);//Centroid of one pixel
        bool cont=true;
        while (cont) {
            binnedpix[nok]=cpix;
            nok++;
            if (nok%10000==0) {
                printf("Working with pixel %ld...\n",nok);
            }
            cpix=getclosest(axes,binning,cra[bin],cdec[bin]);
            if (cpix==-1) break;
            bool crit1=topologicalcrit(cpix,npixbin[bin],tarray,axes);
            bool crit2=morphologicalcrit(cpix,npixbin[bin],tarray,axes);//Estimate the roundness of the possible new bin
            bool crit3=uniformitycrit(cpix,npixbin[bin],tarray,img,mincounts);
            if (!crit1 || !crit2 || !crit3){
                /*if (!crit1) {
                 printf("Crit 1 failed\n");
                 }
                 if (!crit2) {
                 printf("Crit 2 failed\n");
                 }
                 if (!crit3) {
                 printf("Crit 3 failed\n");
                 }*/
                ncb=totcounts(npixbin[bin],tarray,img);
                if (ncb>0.8*mincounts) {//The bin has reached the target number of counts
                    //printf("Bin %d has reached the target number of counts\n",nbin);
                    for (int i=0; i<npixbin[bin]; i++) {
                        int pix=tarray[i];
                        isbinned[pix]=true;
                        nbinned++;
                    }
                    isgood[bin]=true;
                }
                else {
                    isgood[bin]=false;
                }
                tpix[bin]=new int[npixbin[bin]];
                for (int i=0; i<npixbin[bin]; i++) {
                    tpix[bin][i]=tarray[i];
                }
                //printf("npixbin, ncounts: %d  %g\n",npixbin[bin],ncb);
                bin++;
                cont=false;
            }
            else {//we accrete
                int tnb=npixbin[bin];
                tarray[tnb]=cpix;
                npixbin[bin]++;
                binning[cpix]=nbin;
                //Update centroid of the current bin
                getmean(npixbin[bin],axes,tarray,cra[bin],cdec[bin]);
            }
        }
        //Calculate new centroid of binned pixels and determine first pixel of next bin
        double tcx,tcy;
        getmean(nok,axes,binnedpix,tcx,tcy);
        /*if (tcx==cx && tcy==cy) {
         nok++;
         }*/
        cx=tcx;
        cy=tcy;
        cpix=getclosest(axes,binning,cx,cy);
        if (cpix==-1) break;
    }
    
    printf("Adding unsuccessfully binned pixels...\n");
    for (int i=0; i<npix; i++) {
        if (!isbinned[i] && binning[i]>-1) {
            int tbin=findbin(i,axes,bin,isgood,cra,cdec);
            binning[i]=tbin+1;
            int tnb=npixbin[tbin];
            int *transarr=new int[tnb+1];
            for (int j=0; j<tnb; j++) {
                transarr[j]=tpix[tbin][j];
            }
            transarr[tnb]=i;
            delete [] tpix[tbin];
            tpix[tbin]=transarr;
            npixbin[tbin]++;
        }
    }
    printf("Recomputing centroids...\n");
    int nbgood=0;
    for (int i=0; i<bin; i++) {
        if (isgood[i]) {
            getmean(npixbin[i],axes,tpix[i],cra[i],cdec[i]);
            nbgood++;
        }
    }
    printf("Number of bins for Voronoi tessellation: %d\n",nbgood);
    
    //Sort good bins
    double *cxn=new double[nbgood];
    double *cyn=new double[nbgood];
    int ng=0;
    for (int i=0; i<bin; i++) {
        if (isgood[i]) {
            cxn[ng]=cra[i];
            cyn[ng]=cdec[i];
            ng++;
        }
    }
    delete [] cra;
    delete [] cdec;
    delete [] isbinned;
    delete [] tarray;
    for (int i=0; i<npix; i++) {
        delete [] tpix[i];
    }
    delete [] tpix;
    delete [] npixbin;
    delete [] isgood;
    printf("Compute Voronoi tessellation...\n");
    //Scale factor for the use of the Diehl & Statler 2006 method
    //If wvt=false these are fixed to 1.0
    double *scale=new double[nbgood];
    for (int i=0; i<nbgood; i++) {
        scale[i]=1.0;
    }
    double *cpb=new double[nbgood];
    double *sb=new double[nbgood];
    double *esb=new double[nbgood];
    double totdiff=1e10;
    int maxiter=20;
    int iter=0;
    while (totdiff>0.2 && iter<maxiter) {//threshold for fixed point
        totdiff=compute_voronoi(axes,img,binning,nbgood,scale,true,cxn,cyn,cpb)/nbgood;
        printf("Iteration %d, difference %g\n",iter+1,totdiff);
        iter++;
    }
    
    if (ng==nbgood) {
        printf("No empty bins\n");
    }
    else {
        printf("Some bins are empty\n");
    }
    
    printf("Binning completed, now we compute the surface-brightness map...\n");
    int *npb=new int[nbgood];
    double back=0.0;
    double eback=0.0;
    bool isbacksub=false;
    compute_outimg(binning,nbgood,axes,cxn,cyn,pixsize,img,expo,isexp,backimg,isback,
                   back,eback,isbacksub,outimg,error,npix,scale,sb,esb,npb);
    
    
    /*FILE *fc=fopen("centroids.txt","w");
    ng=0;
    for (int i=0; i<nbgood; i++) {
        fprintf(fc,"%d  %g  %g  %g  %g  %g  %d  %g\n",i,cxn[i],cyn[i],cpb[i],sb[i],esb[i],npb[i],scale[i]);
        if (cpb[i]>0.0) {
            ng++;
        }
    }
    fclose(fc);*/
    
    
    int status=0;
    fitsfile *x2;
    char fnam[200];
    sprintf(fnam,"!%s",outfile);
    fits_create_file(&x2,fnam,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        break;
    }
    long start[2]={1,1};
    fits_create_img(x2,-32,2,axes,&status);
    fits_write_pix(x2,TDOUBLE,start,axes[0]*axes[1],outimg,&status);
    if (status!=0){
        printf("Error %d\n",status);
        break;
    }
    double tps=pixsize/3600.;
    fits_write_key(x2,TSTRING,(char *)"CREATOR",(void *)"proffit 1.5",NULL,&status);
    fits_write_key(x2,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
    fits_write_key(x2,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
    fits_write_key(x2,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
    fits_write_key(x2,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CDELT2",&tps,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
    fits_write_key(x2,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
    fits_write_key(x2,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
    fits_create_img(x2,-32,2,axes,&status);
    if (status!=0) {
        printf("Error %d\n",status);
        break;
    }
    fits_write_pix(x2,TDOUBLE,start,npix,error,&status);
    if (status!=0) {
        printf("    Error %d\n",status);
    }
    fits_write_key(x2,TSTRING,(char *)"CREATOR",(void *)"proffit 1.5",NULL,&status);
    fits_write_key(x2,TSTRING,(char *)"CTYPE1",(void *)"RA---TAN",(char *)"LONGPROJ where LONG can be RA, GLON, ELON and PROJ can be CAR, TAN or AIT",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRPIX1",&crpix1,(char *)"Pixel at reference point",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRVAL1",&crval1,(char *)"LONG at the reference value",&status);
    fits_write_key(x2,TSTRING,(char *)"CUNIT1",(void *)"deg",(char *)"Physical units of axis 1",&status);
    fits_write_key(x2,TSTRING,(char *)"CTYPE2",(void *)"DEC--TAN",(char *)"LAT-PROJ where LAT can be DEC, GLAT, ELAT and PROJ can be CAR, TAN or AIT",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRPIX2",&crpix2,(char *)"Pixel at reference point",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CRVAL2",&crval2,(char *)"LAT at the reference value",&status);
    fits_write_key(x2,TSTRING,(char *)"CUNIT2",(void *)"deg",(char *)"Physical units of axis 2",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CDELT1",&cdelt1,(char *)"Element (1,1) of coordinate transf. matrix (default 1)",&status);
    fits_write_key(x2,TDOUBLE,(char *)"CDELT2",&pixsize,(char *)"Element (2,2) of coordinate transf. matrix (default 1)",&status);
    fits_write_key(x2,TSTRING,(char *)"RADECSYS",(void *)"FK5",(char *)"Stellar reference frame",&status);
    fits_write_key(x2,TSTRING,(char *)"EQUINOX",(void *)"2000.0",(char *)"Coordinate system equinox",&status);
    fits_close_file(x2,&status);
    if (status!=0) {
        printf("Program exiting with status %d\n",status);
        break;
    }
    
    printf("Program exited successfully\n");
    
    delete [] binning;
    delete [] cpb;
    delete [] cxn;
    delete [] cyn;
    delete [] scale;
    delete [] sb;
    delete [] esb;

    }
    while (0);
}
