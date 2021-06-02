void fakeimg(double *exposure,long *axes,TF1 *model,double pixsize,double centroid_ra,double centroid_dec,double newexp,double *outimg,bool isback,double *backmap,bool ellipse,double ellang,double aoverb){
    TRandom3 *gr=new TRandom3(0);
    double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
    for (int i=0; i<axes[0]; i++) {
        for (int j=0; j<axes[1]; j++) {
            exposure[j*axes[0]+i]*=newexp/maxexp;
            double posx=(i-centroid_ra)*pixsize*60;//arcmin
            double posy=(j-centroid_dec)*pixsize*60;
            double dist;
            if (!ellipse) {
                dist=sqrt(posx*posx+posy*posy);
            }
            else {
                double xtil=cos(ellang)*posx+sin(ellang)*posy;
                double ytil=-sin(ellang)*posx+cos(ellang)*posy;
                dist=aoverb*sqrt(xtil*xtil+ytil*ytil/aoverb/aoverb);
            }
            double ncp=model->Eval(dist)*pixsize*pixsize*60.*60.*exposure[j*axes[0]+i];
            if (isback && exposure[j*axes[0]+i]>0.0) {
                ncp+=backmap[j*axes[0]+i];
            }
            if (ncp<=0.0) {
                outimg[j*axes[0]+i]=0.0;
            }
            else {
                outimg[j*axes[0]+i]=gr->Poisson(ncp);
                //printf("Exposure, model, ncp, outimg: %g %g %g %g\n",exposure[j*axes[0]+i],model->Eval(dist),ncp,outimg[j*axes[0]+i]);
            }
        }
    }
    gr->Delete();
}
