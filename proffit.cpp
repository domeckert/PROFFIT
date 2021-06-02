#include <stdlib.h>
#include <stdio.h>
#include "fitsio.h"
#include "TMath.h"
#include <iostream>
#include "TH2F.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TRint.h"
#include <TROOT.h>
#include <TApplication.h>
#include "TFile.h"
#include "TTree.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFitter.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include <wcslib/wcs.h>
#include "math.h"
#include "io.h"
#include "bknpow.h"
#include "models.h"
#include "miscellaneous.h"
#include "help.h"
#include "sbprofs.h"
#include "psf.h"
#include "fitting.h"
#include "TLegend.h"
#include "fakeit.h"
#include "multinest.h"
#include "nested_sampling.h"
#include "voronoi.h"

extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };

TROOT root("Rint","The ROOT Interactive Interface", initfuncs);

using namespace std;


void proffit(int argc, char **argv){
do {
  int status=0;
  printf(" \n");
  printf("*******************************\n");
  printf("*** Welcome to proffit v1.5 ***\n");
  printf("*******************************\n");  
  printf(" \n");
  printf("Give me a command (for a list of commands, type help)\n");
  printf(" \n");
  
  // Initialize variables
  
  double *img=NULL;
  bool isimg=false;
  double *exposure=NULL;
  bool isexp=false;
  double *sig=NULL;
  bool issig=false;
  double *backmap=NULL;
  double *cprof=NULL;
  bool isprofile=false;
  double *profile=NULL;
  double *bins=NULL;
  double *eprof=NULL;
  double *ebins=NULL;
	double *effexp=NULL;
	double *area=NULL;
	double *backcounts=NULL;
  double *ecp=NULL;
  double *deltachi=NULL;
  double *backprof=NULL;
  double *deprof=NULL;
  double *edeprof=NULL;
    double *dens=NULL;
    double *edens=NULL;
  double *scat=NULL;
  double *escat=NULL;
  bool isscat=false;
  bool isdepr=false;
  bool isdens=false;
    double *medprof=NULL;
    double *emedprof=NULL;
    bool ismed=false;
    bool scripting=false;
    double *voronoimap=NULL;
    double *voronoierr=NULL;
    bool isvoronoi=false;
  TF1 *model=NULL;
  bool ismod=false;
  char **names=new char*[20];
  for (int i=0;i<20;i++){
    names[i]=new char[200];
  }
  bool *fix=new bool[20];
  for (int i=0;i<20;i++){
    fix[i]=false;
  }
  char *modname=new char[30];
  double *growth=NULL;
  double *egr=NULL;
  bool isgrowth=false;
  bool islimits=false;
  bool logx=true;
  bool logy=true;
  double crpix1,crpix2,crval1,crval2,cdelt1;
  struct wcsprm *wcs_inp=new struct wcsprm;
  double binsize,maxrad,pixsize,centroid_ra,centroid_dec,rad2pix,bin2pix;
  long *axes=new long[2];
  int nbin,nbtot;
  double fitlow,fithigh;
  double ellang,aoverb;
  double mincounts=0.0;
  double angl,angh;
  bool sector=false;
  bool ellipse=false;
  double *psfmat=NULL;
    char imgfile[200];
    char expfile[200];
    char backfile[200];
  bool isr200=false;
  bool isr500=false;
  bool iskpc=false;
	bool fitc=false;
	bool isnested=false;
	char *nestname=new char[200];
	bool *islogpar=NULL;
	double **chains=NULL;
	int npchain;
	double r200obs,r500obs,kpcobs;
	
  // End initialize variables
  
  TApplication theApp("name", &argc, argv);
  TCanvas *c1=new TCanvas("c1","c1");
  c1->SetLogx();
  c1->SetLogy();
  c1->SetWindowSize(689,600);
  gROOT->SetStyle("Plain");
  c1->SetFillColor(kWhite);
  c1->SetFrameBorderMode(0);
  TMultiGraph *m1=new TMultiGraph("m1","m1");
  char *temp=new char[200];
  char *statmet=new char[200];
  sprintf(statmet,"chi2");
  char comm[200];
  int ncom=1;
  bool cont=true;
  while (cont){
	  sprintf(comm,"%d > ",ncom);
	  cout << comm;
	  cin >> temp;
	  ncom++;
	  if (!strcmp(temp,"help")){
		commlist();
		sprintf(comm,"    Command > ");
		printf("%s",comm);
		cin >> temp;		
		help(temp);
    }
      else if (!strcmp(temp,"scripting")){
          sprintf(comm,"    Run in scripting mode? (y/n) > ");
          cout << comm;
          cin >> temp;
          if (!strcmp(temp,"y")||!strcmp(temp,"Y")||!strcmp(temp,"yes")||!strcmp(temp,"true")) {
              scripting=true;
          }
      }
      else if (!strcmp(temp,"readimg")){
          do {
              img=NULL;
              sprintf(comm,"    File name > ");
              cout << comm;
              cin >> imgfile;
              status=getaxes(imgfile,axes);
              img=new double[axes[0]*axes[1]];
              status=readimg(imgfile,img,axes,pixsize,cdelt1,crval1,crval2,crpix1,crpix2,wcs_inp);
              if (status!=0) {
                  printf("    Error %d\n",status);
                  break;
              }
              else {
                  printf("    Image succesfully loaded\n");
                  isimg=true;
              }
              double *cd_inp=new double[4];
              cd_inp[0]=cdelt1;
              cd_inp[1]=0.0;
              cd_inp[2]=0.0;
              cd_inp[3]=pixsize;
          }
          while (0);
      }
      else if (!strcmp(temp,"fakeit")){
          do {
              if (!ismod){
                  printf("    No model loaded\n");
                  break;
              }
              if (!isexp){
                  printf("    Exposure map not loaded\n");
                  break;
              }
              sprintf(comm,"    Total exposure for simulation > ");
              double newexp;
              cout << comm;
              cin >> temp;
              newexp=atof(temp);
              img=new double[axes[0]*axes[1]];
              cout << "    Center (1: image coord, 2: FK5) > ";
              cin >> temp;
              int center=atoi(temp);
              if ((center!=1)&&(center!=2)){
                  printf("    Invalid option %s\n",temp);
                  break;
              }
              centroid_ra=0.0;
              centroid_dec=0.0;
              if (center==1){
                  cout << "       X center (image coord) > ";
                  cin >> temp;
                  centroid_ra=atof(temp);
                  cout << "       Y center (image coord) > ";
                  cin >> temp;
                  centroid_dec=atof(temp);
                  centroid_ra-=1;
                  centroid_dec-=1;
              }
              if (center==2){
                  double pixcrd[2],imgcrd[2],world[2];
                  double inra,indec;
                  int stat;
                  cout << "       RA center (J2000.0) > ";
                  cin >> temp;
                  world[0]=atof(temp);
                  cout << "       Dec center (J2000.0) > ";
                  cin >> temp;
                  world[1]=atof(temp);
                  status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                  centroid_ra=pixcrd[0];
                  centroid_dec=pixcrd[1];
                  printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
                  centroid_ra-=1;
                  centroid_dec-=1;
              }
              img=NULL;
              img=new double[axes[0]*axes[1]];
              fakeimg(exposure,axes,model,pixsize,centroid_ra,centroid_dec,newexp,img,isback,backmap,ellipse,ellang,aoverb);
              printf("    Image succesfully simulated\n");
              isimg=true;
              cout << "    Save simulated image? (y/n) > ";
              cin >> temp;
              if (!strcmp(temp,"y")||!strcmp(temp,"Y")||!strcmp(temp,"yes")||!strcmp(temp,"true")) {
                  cout << "      Output FITS file name > ";
                  cin >> temp;
                  save_img(temp,img,axes,crpix1,crval1,crpix2,crval2,cdelt1,pixsize);
              }
          }
          while (0);
      }
      else if (!strcmp(temp,"readexp")){
          do {
              exposure=NULL;
              sprintf(comm,"    File name > ");
              cout << comm;
              cin >> expfile;
              if (!strcmp(expfile,"none")){
                  exposure=new double[axes[0]*axes[1]];
                  for (int i=0; i<axes[0]*axes[1]; i++) {
                      exposure[i]=1.0;
                  }
                  isexp=true;
              }
              else {
                  if (!isimg) {
                      status=getaxes_expo(expfile,axes,pixsize,cdelt1,crval1,crval2,crpix1,crpix2,wcs_inp);
                      if (status!=0) {
                          printf("Error: unable to read parameters of exposure map\n");
                          break;
                      }
                  }
                  exposure=new double[axes[0]*axes[1]];
                  status=readexp(expfile,axes,exposure);
                  if (status!=0) {
                      break;
                  }
                  else {
                      printf("    Exposure map succesfully loaded\n");
                      isexp=true;
                  }
              }
          }
          while (0);
      }
      else if (!strcmp(temp,"readback")){
          do {
              backmap=NULL;
              sprintf(comm,"    File name > ");
              cout << comm;
              cin >> backfile;
              backmap=new double[axes[0]*axes[1]];
              status=readback(backfile,axes,backmap);
              if (status!=0) {
                  printf("    Error %d\n",status);
                  break;
              }
              else {
                  printf("    Background map succesfully loaded\n");
                  isback=true;
              }
          }
          while (0);
      }
		else if (!strcmp(temp,"model")){
			cout << "    model (type help for a list) > ";
			cin >> modname;
			if (!strcmp(modname,"help")){
				allmodels();
				cout << "    Model > ";
				cin >> temp;
				helpmodels(temp);
			}
			else if (!strcmp(modname,"beta")){
				model=NULL;
				model=new TF1("model",betaprofile,0.01,1e4,4);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc,amp,cc;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc > ";
				cin >> temp;
				rc=atof(temp);
				model->SetParameter(1,rc);
				model->SetParName(1,"rc");
				cout << "      norm > ";
				cin >> temp;
				amp=atof(temp);
				model->SetParameter(2,amp);
				model->SetParName(2,"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(3,cc);
				model->SetParName(3,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc");
				sprintf(names[2],"norm");
				sprintf(names[3],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"doublebeta")){
				model=NULL;
				model=new TF1("model",doublebeta,0.01,1e4,6);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc1,rc2,amp1,amp2,cc;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc1 > ";
				cin >> temp;
				rc1=atof(temp);
				model->SetParameter(1,rc1);
				model->SetParName(1,"rc1");
				cout << "      rc2 > ";
				cin >> temp;
				rc2=atof(temp);
				model->SetParameter(2,rc2);
				model->SetParName(2,"rc2");
				cout << "      ratio > ";
				cin >> temp;
				amp1=atof(temp);
				model->SetParameter(3,amp1);
				model->SetParName(3,"ratio");
				cout << "      norm > ";
				cin >> temp;
				amp2=atof(temp);
				model->SetParameter(4,amp2);
				model->SetParName(4,"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(5,cc);       
				model->SetParName(5,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc1");
				sprintf(names[2],"rc2");
				sprintf(names[3],"ratio");
				sprintf(names[4],"norm");
				sprintf(names[5],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"cuspbeta")){
				model=NULL;
				model=new TF1("model",cuspbeta,0.01,1e4,6);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc,alpha,amp,cc,rs;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc > ";
				cin >> temp;
				rc=atof(temp);
				model->SetParameter(1,rc);
				model->SetParName(1,"rc");
				cout << "      alpha > ";
				cin >> temp;
				alpha=atof(temp);
                model->SetParameter(2,alpha);
				model->SetParName(2,"alpha");
                cout << "      rs > ";
                cin >> temp;
                rs=atof(temp);
                model->SetParameter(3,rs);
                model->SetParName(3,"rs");
				cout << "      norm > ";
				cin >> temp;
				amp=atof(temp);
				model->SetParameter(4,amp);
				model->SetParName(4,"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(5,cc);
				model->SetParName(5,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc");
				sprintf(names[2],"alpha");
                sprintf(names[3],"rs");
				sprintf(names[4],"norm");
				sprintf(names[5],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"const")){
				model=NULL;
				model=new TF1("model",cst,0.01,1e4,1);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				double cc;
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(0,cc);
				model->SetParName(0,"const");
				sprintf(names[0],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"power")){
				model=NULL;
				model=new TF1("model",powerlaw,0.01,1e4,4);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				double norm,alpha,cc,rs;
				cout << "      alpha > ";
				cin >> temp;
				alpha=atof(temp);
				model->SetParameter(0,alpha);
				model->SetParName(0,"alpha");
				sprintf(names[0],"alpha");
                cout << "      rs > ";
                cin >> temp;
                rs=atof(temp);
                model->SetParameter(1,rs);
                model->SetParName(1,"rs");
                sprintf(names[1],"rs");
				cout << "      norm > ";
				cin >> temp;
				norm=atof(temp);
				model->SetParameter(2,norm);
				model->SetParName(2,"norm");
				sprintf(names[2],"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(3,cc);
				model->SetParName(3,"const");
				sprintf(names[3],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"gausbeta")){
				model=NULL;
				model=new TF1("model",gausbeta,0.01,1e4,7);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc,norm,normg,sigma,cc,mu;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc > ";
				cin >> temp;
				rc=atof(temp);
				model->SetParameter(1,rc);
				model->SetParName(1,"rc");
				cout << "      norm > ";
				cin >> temp;
				norm=atof(temp);
				model->SetParameter(2,norm);
				model->SetParName(2,"norm");
                cout << "      mu > ";
                cin >> temp;
                mu=atof(temp);
                model->SetParameter(3,mu);
                model->SetParName(3,"mu");
				cout << "      sigma > ";
				cin >> temp;
				sigma=atof(temp);
				model->SetParameter(4,sigma);
				model->SetParName(4,"sigma");
                cout << "      normg > ";
                cin >> temp;
                normg=atof(temp);
                model->SetParameter(5,normg);
                model->SetParName(5,"normg");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(6,cc);
				model->SetParName(6,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc");
				sprintf(names[2],"norm");
                sprintf(names[3],"mu");
				sprintf(names[4],"sigma");
                sprintf(names[5],"normg");
				sprintf(names[6],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"gausdbeta")){
				model=NULL;
				model=new TF1("model",gausdbeta,0.01,1e4,9);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc1,rc2,rat,norm,normg,sigma,cc,mu;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc1 > ";
				cin >> temp;
				rc1=atof(temp);
				model->SetParameter(1,rc1);
				model->SetParName(1,"rc1");
				cout << "      rc2 > ";
				cin >> temp;
				rc2=atof(temp);
				model->SetParameter(2,rc2);
				model->SetParName(2,"rc2");
				cout << "      ratio > ";
				cin >> temp;
				rat=atof(temp);
				model->SetParameter(3,rat);
				model->SetParName(3,"ratio");
				cout << "      norm > ";
				cin >> temp;
				norm=atof(temp);
				model->SetParameter(4,norm);
				model->SetParName(4,"norm");
                cout << "      mu > ";
                cin >> temp;
                mu=atof(temp);
                model->SetParameter(5,mu);
                model->SetParName(5,"mu");
				cout << "      sigma > ";
				cin >> temp;
				sigma=atof(temp);
				model->SetParameter(6,sigma);
				model->SetParName(6,"sigma");
                cout << "      normg > ";
                cin >> temp;
                normg=atof(temp);
                model->SetParameter(7,normg);
                model->SetParName(7,"normg");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(8,cc);
				model->SetParName(8,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc1");
				sprintf(names[2],"rc2");
				sprintf(names[3],"ratio");
				sprintf(names[4],"norm");
                sprintf(names[5],"mu");
				sprintf(names[6],"sigma");
                sprintf(names[7],"normg");
				sprintf(names[8],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"bknpow")){
				model=NULL;
				model=new TF1("model",bknpow,0.01,1e4,6);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      alpha1 > ";
				double alpha1,alpha2,cutrad,norm1,jump,cc;
				cin >> temp;
				alpha1=atof(temp);
				model->SetParameter(0,alpha1);
				model->SetParName(0,"alpha1");
				cout << "      alpha2 > ";
				cin >> temp;
				alpha2=atof(temp);
				model->SetParameter(1,alpha2);
				model->SetParName(1,"alpha2");
				cout << "      cutrad > ";
				cin >> temp;
				cutrad=atof(temp);
				model->SetParameter(2,cutrad);
				model->SetParName(2,"cutrad");
				cout << "      norm > ";
				cin >> temp;
				norm1=atof(temp);
				model->SetParameter(3,norm1);
				model->SetParName(3,"norm");
				cout << "      jump > ";
				cin >> temp;
				jump=atof(temp);
				model->SetParameter(4,jump);
				model->SetParName(4,"jump");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(5,cc);       
				model->SetParName(5,"const");
				sprintf(names[0],"alpha1");
				sprintf(names[1],"alpha2");
				sprintf(names[2],"cutrad");
				sprintf(names[3],"norm");
				sprintf(names[4],"jump");
				sprintf(names[5],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"bknbeta")){
				model=NULL;
				model=new TF1("model",bknbeta,0.01,1e4,7);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				cout << "      beta > ";
				double beta,rc,alpha2,cutrad,norm1,jump,cc;
				cin >> temp;
				beta=atof(temp);
				model->SetParameter(0,beta);
				model->SetParName(0,"beta");
				cout << "      rc > ";
				cin >> temp;
				rc=atof(temp);
				model->SetParameter(1,rc);
				model->SetParName(1,"rc");
				cout << "      alpha2 > ";
				cin >> temp;
				alpha2=atof(temp);
				model->SetParameter(2,alpha2);
				model->SetParName(2,"alpha2");
				cout << "      cutrad > ";
				cin >> temp;
				cutrad=atof(temp);
				model->SetParameter(3,cutrad);
				model->SetParName(3,"cutrad");
				cout << "      norm > ";
				cin >> temp;
				norm1=atof(temp);
				model->SetParameter(4,norm1);
				model->SetParName(4,"norm");
				cout << "      jump > ";
				cin >> temp;
				jump=atof(temp);
				model->SetParameter(5,jump);
				model->SetParName(5,"jump");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(6,cc);       
				model->SetParName(6,"const");
				sprintf(names[0],"beta");
				sprintf(names[1],"rc");
				sprintf(names[2],"alpha2");
				sprintf(names[3],"cutrad");
				sprintf(names[4],"norm");
				sprintf(names[5],"jump");
				sprintf(names[6],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"backfit")){
				do {
					if (!isback) {
						printf("    No background model loaded/n");
						break;
					}
					if (!isprofile) {
						printf("    No profile loaded\n");
						break;
					}
					model=NULL;
					passback=backprof;
					passbins=bins;
					passebins=ebins;
					model=new TF1("model",funcback,0.01,1e4,2);
					model->SetLineWidth(2);
					model->SetLineColor(kBlue);
					double cc,cb;
					cout << "      const > ";
					cin >> temp;
					cc=atof(temp);
					model->SetParameter(0,cc);
					model->SetParName(0,"const");					
					cout << "      back > ";
					cin >> temp;
					cb=atof(temp);
					model->SetParameter(1,cb);
					model->SetParName(1,"back");
					sprintf(names[0],"const");
					sprintf(names[1],"back");
					ismod=true;
				} while (0);
			}
			else if (!strcmp(modname,"triplepl")){
				model=NULL;
				model=new TF1("model",triplepl,0.01,1e4,7);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				double norm,sl1,sl2,sl3,rc1,rc2,cc;
				cout << "      alpha 1 > ";
				cin >> temp;
				sl1=atof(temp);
				model->SetParameter(0,sl1);
				model->SetParName(0,"alpha1");
				sprintf(names[0],"alpha1");
				cout << "      alpha 2 > ";
				cin >> temp;
				sl2=atof(temp);
				model->SetParameter(1,sl2);
				model->SetParName(1,"alpha2");
				sprintf(names[1],"alpha2");
				cout << "      alpha 3 > ";
				cin >> temp;
				sl3=atof(temp);
				model->SetParameter(2,sl3);
				model->SetParName(2,"alpha3");
				sprintf(names[0],"alpha3");
				cout << "      rc1 > ";
				cin >> temp;
				rc1=atof(temp);
				model->SetParameter(3,rc1);
				model->SetParName(3,"rc1");
				sprintf(names[3],"rc1");
				cout << "      rc2 > ";
				cin >> temp;
				rc2=atof(temp);
				model->SetParameter(4,rc2);
				model->SetParName(4,"rc2");
				sprintf(names[4],"rc2");
				cout << "      norm > ";
				cin >> temp;
				norm=atof(temp);
				model->SetParameter(5,norm);
				model->SetParName(5,"norm");
				sprintf(names[5],"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(6,cc);
				model->SetParName(6,"const");
				sprintf(names[6],"const");
				ismod=true;
			}
			else if (!strcmp(modname,"triplebkn")){
				model=NULL;
				model=new TF1("model",triplebkn,0.01,1e4,7);
				model->SetLineWidth(2);
				model->SetLineColor(kBlue);
				double norm,sl1,sl2,sl3,rc1,rc2,cc;
				cout << "      alpha 1 > ";
				cin >> temp;
				sl1=atof(temp);
				model->SetParameter(0,sl1);
				model->SetParName(0,"alpha1");
				sprintf(names[0],"alpha1");
				cout << "      alpha 2 > ";
				cin >> temp;
				sl2=atof(temp);
				model->SetParameter(1,sl2);
				model->SetParName(1,"alpha2");
				sprintf(names[1],"alpha2");
				cout << "      alpha 3 > ";
				cin >> temp;
				sl3=atof(temp);
				model->SetParameter(2,sl3);
				model->SetParName(2,"alpha3");
				sprintf(names[2],"alpha3");
				cout << "      rc1 > ";
				cin >> temp;
				rc1=atof(temp);
				model->SetParameter(3,rc1);
				model->SetParName(3,"rc1");
				sprintf(names[3],"rc1");
				cout << "      rc2 > ";
				cin >> temp;
				rc2=atof(temp);
				model->SetParameter(4,rc2);
				model->SetParName(4,"rc2");
				sprintf(names[4],"rc2");
				cout << "      norm > ";
				cin >> temp;
				norm=atof(temp);
				model->SetParameter(5,norm);
				model->SetParName(5,"norm");
				sprintf(names[5],"norm");
				cout << "      const > ";
				cin >> temp;
				cc=atof(temp);
				model->SetParameter(6,cc);
				model->SetParName(6,"const");
				sprintf(names[6],"const");
				ismod=true;
			}
            else if (!strcmp(modname,"bknpowgauss")){
                model=NULL;
                model=new TF1("model",bknpowgauss,0.01,1e4,9);
                model->SetLineWidth(2);
                model->SetLineColor(kBlue);
                cout << "      alpha1 > ";
                double alpha1,alpha2,cutrad,norm1,jump,Ng,mug,sigmag,cc;
                cin >> temp;
                alpha1=atof(temp);
                model->SetParameter(0,alpha1);
                model->SetParName(0,"alpha1");
                cout << "      alpha2 > ";
                cin >> temp;
                alpha2=atof(temp);
                model->SetParameter(1,alpha2);
                model->SetParName(1,"alpha2");
                cout << "      cutrad > ";
                cin >> temp;
                cutrad=atof(temp);
                model->SetParameter(2,cutrad);
                model->SetParName(2,"cutrad");
                cout << "      norm > ";
                cin >> temp;
                norm1=atof(temp);
                model->SetParameter(3,norm1);
                model->SetParName(3,"norm");
                cout << "      jump > ";
                cin >> temp;
                jump=atof(temp);
                model->SetParameter(4,jump);
                model->SetParName(4,"jump");
                cout << "      normg > ";
                cin >> temp;
                Ng=atof(temp);
                model->SetParameter(5,Ng);
                model->SetParName(5,"normg");
                cout << "      mu > ";
                cin >> temp;
                mug=atof(temp);
                model->SetParameter(6,mug);
                model->SetParName(6,"mu");
                cout << "      sigma > ";
                cin >> temp;
                sigmag=atof(temp);
                model->SetParameter(7,sigmag);
                model->SetParName(7,"sigma");
                cout << "      const > ";
                cin >> temp;
                cc=atof(temp);
                model->SetParameter(8,cc);
                model->SetParName(8,"const");
                sprintf(names[0],"alpha1");
                sprintf(names[1],"alpha2");
                sprintf(names[2],"cutrad");
                sprintf(names[3],"norm");
                sprintf(names[4],"jump");
                sprintf(names[5],"normg");
                sprintf(names[6],"mu");
                sprintf(names[7],"sigma");
                sprintf(names[8],"const");
                ismod=true;
            }
			else {
				printf("    Unknown model %s\n",modname);
				printf("    For the list of available models type help\n");
			}
			for (int i=0; i<20; i++) {
				fix[i]=false;
			}
		}
		else if (!strcmp(temp,"profile")){
		 do {
			 ellipse=false;
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 cout << "    Center (1: centroid, 2: sb peak, 3: user input (image coord), 4: user input (J2000.0 coord)) > ";
			 cin >> temp;
			 int center=atoi(temp);
			 if ((center!=1)&&(center!=2)&&(center!=3)&&(center!=4)){
				 printf("    Invalid option %s\n",temp);
				 break;
			 }
			 double rmax;
			 if (center==1) {
				 cout << "       Radius of the region to evaluate (arcmin) > ";
				 cin >> temp;
				 rmax=atof(temp)/pixsize/60.;//pixel
				 if (rmax<=0.0) {
					 printf("    Invalid value\n");
					 break;
				 }
			 }
			 centroid_ra=0.0;
			 centroid_dec=0.0;
			 if (center==3){
				 cout << "       X center (image coord) > ";
				 cin >> temp;
				 centroid_ra=atof(temp);
				 cout << "       Y center (image coord) > ";
				 cin >> temp;
				 centroid_dec=atof(temp);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 if (center==4){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
 				 cout << "       RA center (J2000.0) > ";
				 cin >> temp;
				 world[0]=atof(temp);
				 cout << "       Dec center (J2000.0) > ";
				 cin >> temp;
				 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
				 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 cout << "    Bin size (arcsec) > ";
			 cin >> temp;
			 binsize=atof(temp);
			 if (binsize<pixsize*60.*60.){
				 printf("    Error: bin size is smaller than pixel size\n");
				 break;
			 }
			 cout << "    Maximal radius (arcmin) > ";
			 cin >> temp;
			 maxrad=atof(temp);
			 if (maxrad<1e-10){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 cout << "    Logarithmic binning? (y/n) > ";
			 cin >> temp;
			 bool islogbin=false;
			 if (!strcmp(temp,"y")) {
				 islogbin=true;
			 }
			 double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
			 passmaxexp=maxexp;
			 if (center==1){
                 centroid(img,exposure,axes,rmax,pixsize,centroid_ra,centroid_dec);
                 double pixcrd[2],imgcrd[2],world[2];
                 double phi,theta;
                 int stat;
                 pixcrd[0]=centroid_ra+1.;
                 pixcrd[1]=centroid_dec+1.;
                 status = wcsp2s(wcs_inp, 1, 2, pixcrd, imgcrd, &phi,&theta, world,&stat);
				 printf("    Centroid of the image (J2000.0): %g %g\n",world[0],world[1]);
			 }
			 if (center==2){
                 sbpeak(img,exposure,axes,centroid_ra,centroid_dec);
                 double pixcrd[2],imgcrd[2],world[2];
                 double phi,theta;
                 int stat;
                 pixcrd[0]=centroid_ra+1.;
                 pixcrd[1]=centroid_dec+1.;
                 status = wcsp2s(wcs_inp, 1, 2, pixcrd, imgcrd, &phi,&theta, world,&stat);
				 printf("    Surface-brightness peak of the cluster (J2000.0): %g %g\n",world[0],world[1]);
			 }
			 int cra=(int)floor(centroid_ra);
			 int cdec=(int)floor(centroid_dec);
			 if (exposure[cdec*axes[0]+cra]==0.0){
				 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
			 }
			 printf("    Exposure at the centre: %g sec\n",maxexp);
			 bin2pix=binsize/pixsize/3600.;
			 rad2pix=maxrad/pixsize/60.;	
			 nbin=(int)floor(rad2pix/bin2pix);				 				 
			 /*if (islogbin) {
				 nbin/=2;
			 }*/
			 nbtot=nbin;
			 mincounts=0.0;
			 cprof=NULL;
			 profile=NULL;
			 bins=NULL;
			 eprof=NULL;
			 ebins=NULL;
			 area=NULL;
			 effexp=NULL;
			 cprof=new double[nbin];
			 profile=new double[nbin];
			 bins=new double[nbin];
			 eprof=new double[nbin];
			 ebins=new double[nbin];
			 int *numbers=new int[nbin];
			 effexp=new double[nbin];
			 area=NULL;
			 area=new double[nbin];
			 if (!issig) {
				 mk_counts_profile(img,exposure,cprof,numbers,effexp,nbin,axes,centroid_ra,centroid_dec,pixsize,binsize,maxrad,islogbin);
				 mk_sb_profile(img,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
			 }
			 else {
				 mk_sb_sig(img,sig,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
			 }
			 
			 if ((isback) && (!issig)) {
				 nbin=(int)floor(rad2pix/bin2pix);
				 backprof=NULL;
				 backcounts=NULL;
				 bins=NULL;
				 ebins=NULL;
				 backprof=new double[nbin];
				 bins=new double[nbin];
				 ebins=new double[nbin];
				 backcounts=new double[nbin];
				 int *nnn=new int[nbin];
				 double *ebp=new double[nbin];
				 double *eee=new double[nbin];
				 mk_counts_profile(backmap,exposure,backcounts,nnn,eee,nbin,axes,centroid_ra,centroid_dec,pixsize,binsize,maxrad,islogbin);
				 mk_sb_profile(backmap,exposure,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
				 cout << "    Do you want to subtract the background profile? (y/n) ";
				 cin >> temp;
				 if (!strcmp(temp,"y")) {
					 for (int i=0; i<nbin; i++) {
						 profile[i]-=backprof[i];
					 }
				 }
				 delete [] ebp;
				 delete [] eee;
				 delete [] nnn;
			 }
			 for (int i=0; i<nbin; i++) {
				 area[i]=numbers[i]*pixsize*60.*pixsize*60.;
			 }
			 delete [] numbers;
			 isprofile=true;
			 isdepr=false;
			 sector=false;
			 ellipse=false;
			 ispsf=false;
		 }
		 while (0);
     }
     else if (!strcmp(temp,"ellipse")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 cout << "    Center (1: centroid, 2: sb peak, 3: user input (image coord), 4: user input (J2000.0 coord)) > ";
			 cin >> temp;
			 int center=atoi(temp);
			 if ((center!=1)&&(center!=2)&&(center!=3)&&(center!=4)){
				 printf("    Invalid option %s\n",temp);
				 break;
			 }
			 double rmax;
			 if (center==1) {
				 cout << "       Radius of the region to evaluate (arcmin) > ";
				 cin >> temp;
				 rmax=atof(temp)/pixsize/60.;//pixel
				 if (rmax<=0.0) {
					 printf("    Invalid value\n");
					 break;
				 }
			 }
			 centroid_ra=0.0;
			 centroid_dec=0.0;
			 if (center==3){
				 cout << "       X center (image coord) > ";
				 cin >> temp;
				 centroid_ra=atof(temp);
				 cout << "       Y center (image coord) > ";
				 cin >> temp;
				 centroid_dec=atof(temp);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 if (center==4){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
 				 cout << "       RA center (J2000.0) > ";
				 cin >> temp;
				 world[0]=atof(temp);
				 cout << "       Dec center (J2000.0) > ";
				 cin >> temp;
				 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
				 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 cout << "    Bin size (arcsec) > ";
			 cin >> temp;
			 binsize=atof(temp);
			 if (binsize<pixsize*60.*60.){
				 printf("    Error: bin size is smaller than pixel size\n");
				 break;
			 }
			 cout << "    Maximal radius (arcmin) > ";
			 cin >> temp;
			 maxrad=atof(temp);
			 if (maxrad<1e-10){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 cout << "    Angle between RA axis and major axis (degrees) > ";
			 cin >> temp;
			 double tta=atof(temp)-90.;
			 if (tta<-90 || tta>270){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 ellang=tta*TMath::Pi()/180.;
			 cout << "    Ratio between major and minor axis > ";
			 cin >> temp;
			 aoverb=atof(temp);
			 if (aoverb<1.0){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 else if (aoverb==1.0) {
				 printf("    You could have just as well used the \"profile\" command, you idiot\n");
			 }
			 cout << "    Logarithmic binning? (y/n) > ";
			 cin >> temp;
			 bool islogbin=false;
			 if (!strcmp(temp,"y")) {
				 islogbin=true;
			 }
			 double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
			 passmaxexp=maxexp;
			 if (center==1){
                 centroid(img,exposure,axes,rmax,pixsize,centroid_ra,centroid_dec);
                 double pixcrd[2],imgcrd[2],world[2];
                 double phi,theta;
                 int stat;
                 pixcrd[0]=centroid_ra+1.;
                 pixcrd[1]=centroid_dec+1.;
                 status = wcsp2s(wcs_inp, 1, 2, pixcrd, imgcrd, &phi,&theta, world,&stat);
				 printf("    Centroid of the image (J2000.0): %g %g\n",world[0],world[1]);
			 }
			 if (center==2){
                 sbpeak(img,exposure,axes,centroid_ra,centroid_dec);
                 double pixcrd[2],imgcrd[2],world[2];
                 double phi,theta;
                 int stat;
                 pixcrd[0]=centroid_ra+1.;
                 pixcrd[1]=centroid_dec+1.;
                 status = wcsp2s(wcs_inp, 1, 2, pixcrd, imgcrd, &phi,&theta, world,&stat);
				 printf("    Surface-brightness peak of the cluster (J2000.0): %g %g\n",world[0],world[1]);
			 }
			 int cra=(int)floor(centroid_ra);
			 int cdec=(int)floor(centroid_dec);
			 if (exposure[cdec*axes[0]+cra]==0.0){
				 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
			 }
			 printf("    Exposure at the centre: %g sec\n",maxexp);
			 bin2pix=binsize/pixsize/3600;
			 rad2pix=maxrad/pixsize/60;
			 nbin=(int)floor(rad2pix/bin2pix);				 				 
			 nbtot=nbin;
			 mincounts=0.0;
			 cprof=NULL;
			 profile=NULL;
			 bins=NULL;
			 eprof=NULL;
			 ebins=NULL;
			 cprof=new double[nbin];
			 profile=new double[nbin];
			 bins=new double[nbin];
			 eprof=new double[nbin];
			 ebins=new double[nbin];
			 effexp=NULL;
			 int *numbers=new int[nbin];
			 effexp=new double[nbin];
			 area=NULL;
			 area=new double[nbin];
			 if (!issig) {
				 mk_ellipse_counts(img,exposure,cprof,numbers,effexp,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 mk_profile_ellipse(img,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
			 }
			 else {
				 mk_sb_sig(img,sig,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
			 }
			 
			 if ((isback) && (!issig)) {
				 nbin=(int)floor(rad2pix/bin2pix);
				 backprof=NULL;
				 bins=NULL;
				 ebins=NULL;
				 backprof=new double[nbin];
				 bins=new double[nbin];
				 ebins=new double[nbin];
				 backcounts=NULL;
				 backcounts=new double[nbin];
				 int *nnn=new int[nbin];
				 double *eee=new double[nbin];
				 double *ebp=new double[nbin];
				 mk_ellipse_counts(backmap,exposure,backcounts,nnn,eee,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 mk_profile_ellipse(backmap,exposure,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 cout << "    Do you want to subtract the background profile? (y/n) ";
				 cin >> temp;
				 if (!strcmp(temp,"y")) {
					 for (int i=0; i<nbin; i++) {
						 profile[i]-=backprof[i];
					 }
				 }
				 delete [] nnn;
				 delete [] eee;
				 delete [] ebp;
			 }
			 for (int i=0; i<nbin; i++) {
				 area[i]=numbers[i]*pixsize*60.*pixsize*60.;
			 }
			 delete [] numbers;
			 isprofile=true;
			 ellipse=true;
			 isdepr=false;
			 sector=false;
			 ispsf=false;
		 }
		 while (0);
     }
     else if (!strcmp(temp,"fitcounts")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model loaded\n");
				 break;
			 }
			 char *pars=new char[200];
			 int npfunc=model->GetNpar();
			 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			 int modw=model->GetLineWidth();
			 int modc=model->GetLineColor();
			 ftemp->SetLineWidth(modw);
			 ftemp->SetLineColor(modc);
			 for (int k=0;k<npfunc;k++){
				 pars=(char *)model->GetParName(k);
				 ftemp->SetParName(k,pars);
				 double tpar=model->GetParameter(k);
				 if (fix[k]){
					 ftemp->FixParameter(k,tpar);
				 }
				 else {
					 ftemp->SetParameter(k,tpar);
				 }
			 }
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F *hc=new TH1F("hc","hc",nbin,binsh);
			 TH1F *hh=new TH1F("hh","hh",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 hc->SetBinContent(i+1,cprof[i]);
				 hc->SetBinError(i+1,effexp[i]*area[i]);
				 hh->SetBinContent(i+1,profile[i]);
				 double toterr=sqrt(eprof[i]*eprof[i]+syserr/100.*profile[i]*syserr/100.*profile[i]);
				 hh->SetBinError(i+1,toterr);
			 }
			 histpsfmat=NULL;
			 if (ispsf) {
				 histpsfmat=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
				 for (int i=0; i<nbin; i++) {
					 for (int j=0; j<nbin; j++) {
						 histpsfmat->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
					 }
				 }
			 }
			 if (!islimits){
				 fitlow=0.0;
				 fithigh=maxrad;
			 }
			 passfitlow=fitlow;
			 passfithigh=fithigh;
			 passfmod=ftemp;
			 passhist=hc;
			 passbck=backcounts;
			 TFitter* minimizer=new TFitter(npfunc);
			 double p1 = -1;
			 minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
			 if (!strcmp(statmet,"chi2")) {
				 minimizer->SetFCN(fcncounts);				 
			 }
			 else {
				 minimizer->SetFCN(fcnlikehcounts);
			 }
			 for (int i=0; i<npfunc; i++) {
				 pars=(char *)ftemp->GetParName(i);
				 double tpar=ftemp->GetParameter(i);
				 minimizer->SetParameter(i,pars,tpar,1,0,0);
				 if (fix[i]) {
					 minimizer->FixParameter(i);
				 }
			 }
			 minimizer->ExecuteCommand("MIGRAD",0,0);
			 int npoints=ftemp->GetNumberFitPoints();
			 int npused=minimizer->GetNumberFreeParameters();
			 Double_t amin,edm,errdef;
			 Int_t nvpar,nparx,icstat;
			 TMinuit *gm=minimizer->GetMinuit();
			 gm->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
			 gm->mnprin(3,amin);
			 for (int i=0; i<npfunc; i++) {
				 double pp=minimizer->GetParameter(i);
				 double epp=minimizer->GetParError(i);
				 ftemp->SetParameter(i,pp);
				 ftemp->SetParError(i,epp);
				 model->SetParameter(i,pp);
				 model->SetParError(i,epp);
			 }
			 if (!strcmp(statmet,"chi2")){
				 int nlib=npoints-npused;
				 double prob=TMath::Prob(amin,nlib);
				 printf("    Chi-squared = %g for %d d.o.f\n",amin,nlib);
				 printf("    Reduced chi-squared = %g\n",amin/nlib);
				 printf("    Null-hypothesis probability = %g\n",prob);
			 }
             else {
                 int nlib=npoints-npused;
                 printf("    Minimum C-statistic = %g for %d d.o.f\n",amin,nlib);
                 printf("    Reduced C-statistic = %g\n",amin/nlib);
             }
			 c1->Clear();
			 TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
			 if (logx){
				 pad1->SetLogx();
			 }
			 if (logy){
				 pad1->SetLogy();
			 }
			 pad1->SetBottomMargin(0);
			 pad1->SetLeftMargin(0.15);
			 pad1->SetTopMargin(0.05);
			 pad1->SetFillColor(kWhite);
			 pad1->SetTickx();
			 pad1->SetTicky();
			 pad1->Draw();
			 pad1->cd();
			 double maxy,miny;
			 maxy=maxinrange(hh,fitlow,fithigh);
			 miny=mininrange(hh,fitlow,fithigh);
			 hh->SetTitle("");
			 hh->SetStats(false);
			 hh->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
			 hh->GetYaxis()->CenterTitle();
			 hh->SetAxisRange(fitlow,fithigh);
			 hh->SetLabelSize(0.05,"Y");
			 hh->SetTitleSize(0.05,"Y");
			 hh->SetTitleOffset(1.1,"Y");
			 hh->Draw();
			 double *modpsf=new double[nbin];
			 if (!ispsf){
				 model->Draw("same");
				 if (!strcmp(modname,"doublebeta")){
					 TF1 *comp1=new TF1("comp1",betaprofile,0.,1e4,4);
					 TF1 *comp2=new TF1("comp2",betaprofile,0.,1e4,4);
					 double betafit=model->GetParameter(0);
					 comp1->SetParameter(0,betafit);
					 comp2->SetParameter(0,betafit);
					 double rc1=model->GetParameter(1);
					 comp1->SetParameter(1,rc1);
					 double rc2=model->GetParameter(2);
					 comp2->SetParameter(1,rc2);
					 double ratio=model->GetParameter(3);
					 double norm=model->GetParameter(4);
					 double amp2=ratio*norm;
					 comp1->SetParameter(2,norm);
					 comp2->SetParameter(2,amp2);
					 comp1->SetParameter(3,0.0);
					 comp2->SetParameter(3,0.0);
					 comp1->SetLineStyle(2);
					 comp1->SetLineWidth(1);
					 comp2->SetLineStyle(2);
					 comp2->SetLineWidth(2);
					 comp1->SetLineColor(kBlack);
					 comp2->SetLineColor(kBlack);
					 comp1->Draw("same");
					 comp2->Draw("same");
				 }
			 }
			 else {
				 for (int i=0; i<nbin; i++) {
					 double ll=0.0;
					 for (int j=0; j<nbin; j++) {
						 double x = bins[j];
						 ll+=model->Eval(x)*histpsfmat->GetBinContent(i+1,j+1);
					 }
					 modpsf[i]=ll;					 
				 }
				 TGraph *ggm=new TGraph(nbin,bins,modpsf);
				 ggm->SetLineColor(kBlue);
				 ggm->SetLineWidth(2);
				 ggm->Draw("CP");
			 }
			 c1->cd();
			 TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
			 if (logx){
				 pad2->SetLogx();
			 }
			 pad2->SetLeftMargin(0.15);
			 pad2->SetBottomMargin(0.28); 
			 pad2->SetTopMargin(0);
			 pad2->SetFillColor(kWhite);
			 pad2->SetTickx();
			 pad2->SetTicky();
			 pad2->Draw();
			 pad2->cd();
			 TH1F *hhd=new TH1F("hhd","hhd",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 double bin=hh->GetBinCenter(i+1);
				 double dat=hh->GetBinContent(i+1);
				 double edat=hh->GetBinError(i+1);
				 if (edat>0.0){
					 double dc;
					 if (!ispsf) {
						 dc=(dat-model->Eval(bin))/edat;						 
					 }
					 else {
						 dc=(dat-modpsf[i])/edat;
					 }
					 hhd->SetBinContent(i+1,dc);
					 hhd->SetBinError(i+1,1.0);
				 }
				 else {
					 double dc=dat-model->Eval(bin);
					 hhd->SetBinContent(i+1,dc);
					 hhd->SetBinError(i+1,1.0);
				 }
			 }
			 hhd->SetTitle("");
			 hhd->SetStats(false);
			 hhd->SetYTitle("#chi");
			 hhd->GetYaxis()->CenterTitle();
			 hhd->SetXTitle("Distance [arcmin]");
			 if (isr500){
				 hhd->SetXTitle("r/r_{500}");
			 }
			 if (isr200){
				 hhd->SetXTitle("r/r_{200}");
			 }
             if (iskpc) {
                 hhd->SetXTitle("Radius [kpc]");
             }
			 hhd->GetXaxis()->CenterTitle();
			 hhd->SetAxisRange(fitlow,fithigh);
			 hhd->SetLabelSize(0.09,"X");
			 hhd->SetTitleSize(0.09,"X");
			 hhd->SetTitleOffset(1.25,"X");
			 hhd->SetLabelSize(0.09,"Y");
			 hhd->SetTitleSize(0.12,"Y");
			 hhd->SetTitleOffset(0.35,"Y");
			 hhd->Draw();
			 double xmin=hhd->GetXaxis()->GetXmin();
			 double xmax=hhd->GetXaxis()->GetXmax();
			 if (islimits){
				 xmin=fitlow;
				 xmax=fithigh;
			 }			 
			 TLine *ll=new TLine(xmin,0.,xmax,0.);
			 ll->Draw();
			 c1->cd();
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
             else {
		     c1->SaveAs("fitcounts.pdf");
	     }
			 if (ispsf) {
				 histpsfmat->Delete();
			 }
			 delete ftemp;
			 delete hh;
			 delete hc;
			 delete hhd;
			 fitc=true;
		 }while(0);
     }
     else if (!strcmp(temp,"fit")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model yet loaded\n");
				 break;
			 }
			 char *pars=new char[200];
			 int npfunc=model->GetNpar();
			 //TF1 *ftemp=mkfint(modname,(char *)"ftemp");;
			 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");;
			 int modw=model->GetLineWidth();
			 int modc=model->GetLineColor();
			 ftemp->SetLineWidth(modw);
			 ftemp->SetLineColor(modc);
			 for (int k=0;k<npfunc;k++){
				 pars=(char *)model->GetParName(k);
				 ftemp->SetParName(k,pars);
				 double tpar=model->GetParameter(k);
				 if (fix[k]){
					 ftemp->FixParameter(k,tpar);
				 }
				 else {
					 ftemp->SetParameter(k,tpar);
				 }
			 }
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F *hh=new TH1F("hh","hh",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				hh->SetBinContent(i+1,profile[i]);
				double toterr=sqrt(eprof[i]*eprof[i]+syserr/100.*profile[i]*syserr/100.*profile[i]);
				hh->SetBinError(i+1,toterr);
			 }
			 histpsfmat=NULL;
			 if (ispsf) {
				 histpsfmat=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
				 for (int i=0; i<nbin; i++) {
					 for (int j=0; j<nbin; j++) {
						 histpsfmat->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
					 }
				 }
			 }
			 if (!islimits){
				 fitlow=0.0;
				 fithigh=maxrad;
			 }
			 passfitlow=fitlow;
			 passfithigh=fithigh;
			 passfmod=ftemp;
			 passhist=hh;
			 TFitter* minimizer=new TFitter(npfunc);
			 double p1 = -1;
			 minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
			 if (!strcmp(statmet,"chi2")) {
				 minimizer->SetFCN(fcnstd);				 
			 }
			 else {
				 minimizer->SetFCN(fcnstdlikeh);
			 }
			 for (int i=0; i<npfunc; i++) {
				pars=(char *)ftemp->GetParName(i);
				double tpar=ftemp->GetParameter(i);
				minimizer->SetParameter(i,pars,tpar,1,0,0);
				if (fix[i]) {
					minimizer->FixParameter(i);
				}
			 }
			 minimizer->ExecuteCommand("MIGRAD",0,0);
			 int npoints=ftemp->GetNumberFitPoints();
			 int npused=minimizer->GetNumberFreeParameters();
			 Double_t amin,edm,errdef;
			 Int_t nvpar,nparx,icstat;
			 TMinuit *gm=minimizer->GetMinuit();
			 gm->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
			 gm->mnprin(3,amin);
			 for (int i=0; i<npfunc; i++) {
				double pp=minimizer->GetParameter(i);
				double epp=minimizer->GetParError(i);
				ftemp->SetParameter(i,pp);
				ftemp->SetParError(i,epp);
				model->SetParameter(i,pp);
				model->SetParError(i,epp);
			 }
			 if (!strcmp(statmet,"chi2")){
				 int nlib=npoints-npused;
				 double prob=TMath::Prob(amin,nlib);
				 printf("    Chi-squared = %g for %d d.o.f\n",amin,nlib);
				 printf("    Reduced chi-squared = %g\n",amin/nlib);
				 printf("    Null-hypothesis probability = %g\n",prob);
			 }
             else {
                 int nlib=npoints-npused;
                 printf("    Minimum C-statistic = %g for %d d.o.f\n",amin,nlib);
                 printf("    Reduced C-statistic = %g\n",amin/nlib);
             }
			 c1->Clear();
			 TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
			 if (logx){
				 pad1->SetLogx();
			 }
			 if (logy){
				 pad1->SetLogy();
			 }
			 pad1->SetBottomMargin(0);
			 pad1->SetLeftMargin(0.15);
			 pad1->SetTopMargin(0.05);
			 pad1->SetFillColor(kWhite);
			 pad1->SetTickx();
			 pad1->SetTicky();
			 pad1->Draw();
			 pad1->cd();
			 double maxy,miny;
			 maxy=maxinrange(hh,fitlow,fithigh);
			 miny=mininrange(hh,fitlow,fithigh);
			 hh->SetTitle("");
			 hh->SetStats(false);
			 hh->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
			 hh->GetYaxis()->CenterTitle();
			 hh->SetAxisRange(fitlow,fithigh);
			 hh->SetLabelSize(0.05,"Y");
			 hh->SetTitleSize(0.05,"Y");
			 hh->SetTitleOffset(1.1,"Y");
			 hh->Draw();
			 double *modpsf=new double[nbin];
			 if (!ispsf){
				 model->Draw("same");
				 if (!strcmp(modname,"doublebeta")){
					 TF1 *comp1=new TF1("comp1",betaprofile,0.,1e4,4);
					 TF1 *comp2=new TF1("comp2",betaprofile,0.,1e4,4);
					 double betafit=model->GetParameter(0);
					 comp1->SetParameter(0,betafit);
					 comp2->SetParameter(0,betafit);
					 double rc1=model->GetParameter(1);
					 comp1->SetParameter(1,rc1);
					 double rc2=model->GetParameter(2);
					 comp2->SetParameter(1,rc2);
					 double ratio=model->GetParameter(3);
					 double norm=model->GetParameter(4);
					 double amp2=ratio*norm;
					 comp1->SetParameter(2,norm);
					 comp2->SetParameter(2,amp2);
					 comp1->SetParameter(3,0.0);
					 comp2->SetParameter(3,0.0);
					 comp1->SetLineStyle(2);
					 comp1->SetLineWidth(1);
					 comp2->SetLineStyle(2);
					 comp2->SetLineWidth(2);
					 comp1->SetLineColor(kBlack);
					 comp2->SetLineColor(kBlack);
					 comp1->Draw("same");
					 comp2->Draw("same");
				 }
			 }
			 else {
				 for (int i=0; i<nbin; i++) {
					 double ll=0.0;
					 for (int j=0; j<nbin; j++) {
						 double x = bins[j];
						 ll+=model->Eval(x)*histpsfmat->GetBinContent(i+1,j+1);
					 }
					 modpsf[i]=ll;					 
				 }
				 TGraph *ggm=new TGraph(nbin,bins,modpsf);
				 ggm->SetLineColor(kBlue);
				 ggm->SetLineWidth(2);
				 ggm->Draw("CP");
			 }
			 c1->cd();
             if (!strcmp(modname,"bknpow")) {
                 TPad *p3=new TPad("p3","p3",0.6,0.7,0.87,0.95);
                 p3->SetLogx();
                 p3->SetLogy();
                 p3->SetTickx();
                 p3->SetTicky();
                 p3->Draw();
                 p3->cd();
                 TF1 *fne=new TF1("",bknpl,fitlow,fithigh,5);
                 for (int i=0; i<5; i++) {
                     double tp=model->GetParameter(i);
                     fne->SetParameter(i,tp);
                 }
                 fne->SetLineWidth(2);
                 fne->SetLineColor(kBlue);
                 fne->GetXaxis()->SetLabelSize(0.0);
                 fne->GetYaxis()->SetLabelSize(0.0);
                 fne->Draw();
                 c1->cd();
             }
             if (!strcmp(modname,"triplebkn")) {
                 TPad *p3=new TPad("p3","p3",0.6,0.7,0.87,0.95);
                 p3->SetLogx();
                 p3->SetLogy();
                 p3->SetTickx();
                 p3->SetTicky();
                 p3->Draw();
                 p3->cd();
                 TF1 *fne=new TF1("",triplepl,fitlow,fithigh,7);
                 for (int i=0; i<6; i++) {
                     double tp=model->GetParameter(i);
                     fne->SetParameter(i,tp);
                 }
                 fne->SetLineWidth(2);
                 fne->SetParameter(6,0.0);
                 fne->SetLineColor(kBlue);
                 fne->GetXaxis()->SetLabelSize(0.0);
                 fne->GetYaxis()->SetLabelSize(0.0);
                 fne->Draw();
                 c1->cd();                 
             }
			 TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
			 if (logx){
				 pad2->SetLogx();
			 }
			 pad2->SetLeftMargin(0.15);
			 pad2->SetBottomMargin(0.28); 
			 pad2->SetTopMargin(0);
			 pad2->SetFillColor(kWhite);
			 pad2->SetTickx();
			 pad2->SetTicky();
			 pad2->Draw();
			 pad2->cd();
			 TH1F *hhd=new TH1F("hhd","hhd",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 double bin=hh->GetBinCenter(i+1);
				 double dat=hh->GetBinContent(i+1);
				 double edat=hh->GetBinError(i+1);
				 if (edat>0.0){
					 double dc;
					 if (!ispsf) {
						 dc=(dat-model->Eval(bin))/edat;						 
					 }
					 else {
						 dc=(dat-modpsf[i])/edat;
					 }
					 hhd->SetBinContent(i+1,dc);
					 hhd->SetBinError(i+1,1.0);
				 }
				 else {
					 double dc=dat-model->Eval(bin);
					 hhd->SetBinContent(i+1,dc);
					 hhd->SetBinError(i+1,1.0);
				 }
			 }
			 hhd->SetTitle("");
			 hhd->SetStats(false);
			 hhd->SetYTitle("#chi");
			 hhd->GetYaxis()->CenterTitle();
			 hhd->SetXTitle("Distance [arcmin]");
			 if (isr500){
				 hhd->SetXTitle("r/r_{500}");
			 }
			 if (isr200){
				 hhd->SetXTitle("r/r_{200}");
			 }
             if (iskpc) {
                 hhd->SetXTitle("Radius [kpc]");
             }
			 hhd->GetXaxis()->CenterTitle();
			 hhd->SetAxisRange(fitlow,fithigh);
			 hhd->SetLabelSize(0.09,"X");
			 hhd->SetTitleSize(0.09,"X");
			 hhd->SetTitleOffset(1.25,"X");
			 hhd->SetLabelSize(0.09,"Y");
			 hhd->SetTitleSize(0.12,"Y");
			 hhd->SetTitleOffset(0.35,"Y");
			 hhd->Draw();
			 double xmin=hhd->GetXaxis()->GetXmin();
			 double xmax=hhd->GetXaxis()->GetXmax();
			 if (islimits){
				 xmin=fitlow;
				 xmax=fithigh;
			 }			 
			 TLine *ll=new TLine(xmin,0.,xmax,0.);
			 ll->Draw();
			 c1->cd();
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
             else {
		     c1->SaveAs("fit.pdf");
	     }
			 if (ispsf) {
				 histpsfmat->Delete();
			 }
			 delete [] binsh;
			 delete ftemp;
			 delete hh;
			 delete hhd;
			 fitc=false;
		 }while(0);
     }
	 else if (!strcmp(temp,"contour")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model yet loaded\n");
				 break;
			 }
			 int par1,par2;
			 char *pars=new char[200];
			 int npfunc=model->GetNpar();
			 cout << "       Parameter 1  > ";
			 cin >> temp;
			 par1=atoi(temp);
			 par1--;
			 if (par1<0 || par1>=npfunc){
				 printf("    Invalid parameter %d\n",par1);
				 break;
			 }
			 cout << "       Parameter 2  > ";
			 cin >> temp;
			 par2=atoi(temp);
			 par2--;
			 if (par2<0 || par1>=npfunc){
				 printf("    Invalid parameter %d\n",par2);
				 break;
			 }
			 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");;
			 int modw=model->GetLineWidth();
			 int modc=model->GetLineColor();
			 ftemp->SetLineWidth(modw);
			 ftemp->SetLineColor(modc);
			 for (int k=0;k<npfunc;k++){
				 pars=(char *)model->GetParName(k);
				 ftemp->SetParName(k,pars);
				 double tpar=model->GetParameter(k);
				 if (fix[k]){
					 ftemp->FixParameter(k,tpar);
				 }
				 else {
					 ftemp->SetParameter(k,tpar);
				 }
			 }
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F *hc=new TH1F("hc","hc",nbin,binsh);
			 TH1F *hh=new TH1F("hh","hh",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 hc->SetBinContent(i+1,cprof[i]);
				 hc->SetBinError(i+1,effexp[i]*area[i]);
				 hh->SetBinContent(i+1,profile[i]);
				 double toterr=sqrt(eprof[i]*eprof[i]+syserr/100.*profile[i]*syserr/100.*profile[i]);
				 hh->SetBinError(i+1,toterr);
			 }
			 histpsfmat=NULL;
			 if (ispsf) {
				 histpsfmat=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
				 for (int i=0; i<nbin; i++) {
					 for (int j=0; j<nbin; j++) {
						 histpsfmat->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
					 }
				 }
			 }
			 if (!islimits){
				 fitlow=0.0;
				 fithigh=maxrad;
			 }
			 passfitlow=fitlow;
			 passfithigh=fithigh;
			 passfmod=ftemp;
			 if (fitc) {
				 passhist=hc;
			 }
			 else {
				 passhist=hh;
			 }
			 TFitter* minimizer=new TFitter(npfunc);
			 double p1 = -1;
			 minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
			 if (!strcmp(statmet,"chi2")&&(!fitc)) {
				 minimizer->SetFCN(fcnstd);				 
			 }
			 else if (!strcmp(statmet,"cash")&&(!fitc)) {
				 minimizer->SetFCN(fcnstdlikeh);
			 }
			 else if (!strcmp(statmet,"chi2")&&(fitc)) {
				 minimizer->SetFCN(fcncounts);
			 }
			 else {
				 minimizer->SetFCN(fcnlikehcounts);
			 }
			 for (int i=0; i<npfunc; i++) {
				 pars=(char *)ftemp->GetParName(i);
				 double tpar=ftemp->GetParameter(i);
				 minimizer->SetParameter(i,pars,tpar,1,0,0);
				 if (fix[i]) {
					 minimizer->FixParameter(i);
				 }
			 }
			 minimizer->ExecuteCommand("MIGRAD",0,0);
			 TMinuit *gm=minimizer->GetMinuit();
			 m1=NULL;
			 TMultiGraph *m1=new TMultiGraph("m1","m1");
			 m1->SetTitle();
			 gm->SetErrorDef(9.0);
			 TGraph *gg3=(TGraph*)gm->Contour(50,par1,par2);
			 int stat=gMinuit->GetStatus();
			 if (stat==-1 || stat==1){
				 printf("   Impossible to make contours\n");
				 delete ftemp;
				 delete hh;
				 break;
			 }
			 gg3->SetMarkerColor(kGreen);
			 gg3->SetLineColor(kGreen);
			 m1->Add(gg3,"CP");
			 m1->Draw("A");
			 gm->SetErrorDef(4.0);
			 TGraph *gg2=(TGraph*)gm->Contour(50,par1,par2);
			 stat=gMinuit->GetStatus();
			 if (stat==-1 || stat==1){
				 printf("   Impossible to make contours\n");
				 delete ftemp;
				 delete hh;
				 break;
			 }
			 gg2->SetMarkerColor(kRed);
			 gg2->SetLineColor(kRed);
			 m1->Add(gg2,"CP");
			 m1->Draw("A");
			 gm->SetErrorDef(1.0);
			 TGraph *gg=(TGraph*)gm->Contour(50,par1,par2);
			 stat=gMinuit->GetStatus();
			 if (stat==-1 || stat==1){
				 printf("   Impossible to make contours\n");
				 delete ftemp;
				 delete hh;
				 break;
			 }
			 gg->SetMarkerColor(kBlue);
			 gg->SetLineColor(kBlue);
			 c1->Clear();
			 c1->SetLogx(0);
			 c1->SetLogy(0);
			 m1->Add(gg,"CP");
			 m1->Draw("A");
			 char nn[100];
			 sprintf(nn,"%s",names[par1]);
			 m1->GetXaxis()->SetTitle(nn);
			 sprintf(nn,"%s",names[par2]);
			 m1->GetYaxis()->SetTitle(nn);
			 m1->Draw("A");
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
			 if (logx) c1->SetLogx();
			 if (logy) c1->SetLogy();
			 if (ispsf) {
				 histpsfmat->Delete();
			 }
			 delete [] binsh;
			 delete ftemp;
			 delete hh;
			 delete hc;
			 delete gg;
			 delete gg2;
		 }while(0);
		 
	 }
	 else if (!strcmp(temp,"error")){
		 do {
			 if (!isimg && !isnested){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp && !isnested){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 if (!isprofile && !isnested){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model yet loaded\n");
				 break;
			 }
			 int par1;
			 char *pars=new char[200];
			 int npfunc=model->GetNpar();
			 cout << "       Parameter  > ";
			 cin >> temp;
			 par1=atoi(temp);
			 if (par1<1 || par1>npfunc){
				 printf("    Invalid parameter %d\n",par1);
				 break;
			 }
			 cout << "       Confidence level (%)  > ";
			 cin >> temp;
			 double cper=atof(temp);
			 if (cper==0.0){
				 printf("    Invalid parameter %d\n",par1);
				 break;
			 }
			if (!isnested){
				 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");;
				 int modw=model->GetLineWidth();
				 int modc=model->GetLineColor();
				 ftemp->SetLineWidth(modw);
				 ftemp->SetLineColor(modc);
				 for (int k=0;k<npfunc;k++){
					 pars=(char *)model->GetParName(k);
					 ftemp->SetParName(k,pars);
					 double tpar=model->GetParameter(k);
					 if (fix[k]){
						 ftemp->FixParameter(k,tpar);
					 }
					 else {
						 ftemp->SetParameter(k,tpar);
					 }
				 }
				 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
				 TH1F *hc=new TH1F("hc","hc",nbin,binsh);
				 TH1F *hh=new TH1F("hh","hh",nbin,binsh);
				 for (int i=0;i<nbin;i++){
					 hc->SetBinContent(i+1,cprof[i]);
					 hc->SetBinError(i+1,effexp[i]*area[i]);
					 hh->SetBinContent(i+1,profile[i]);
					 double toterr=sqrt(eprof[i]*eprof[i]+syserr/100.*profile[i]*syserr/100.*profile[i]);
					 hh->SetBinError(i+1,toterr);
				 }
				 histpsfmat=NULL;
				 if (ispsf) {
					 histpsfmat=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
					 for (int i=0; i<nbin; i++) {
						 for (int j=0; j<nbin; j++) {
							 histpsfmat->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
						 }
					 }
				 }
				 if (!islimits){
					 fitlow=0.0;
					 fithigh=maxrad;
				 }
				 passfitlow=fitlow;
				 passfithigh=fithigh;
				 passfmod=ftemp;
				 if (fitc) {
					 passhist=hc;
				 }
				 else {
					 passhist=hh;
				 }
				 TFitter* minimizer=new TFitter(npfunc);
				 double p1 = -1;
				 minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
				 if (!strcmp(statmet,"chi2")&&(!fitc)) {
					 minimizer->SetFCN(fcnstd);				 
				 }
				 else if (!strcmp(statmet,"cash")&&(!fitc)) {
					 minimizer->SetFCN(fcnstdlikeh);
				 }
				 else if (!strcmp(statmet,"chi2")&&(fitc)) {
					 minimizer->SetFCN(fcncounts);
				 }
				 else {
					 minimizer->SetFCN(fcnlikehcounts);
				 }
				 for (int i=0; i<npfunc; i++) {
					 pars=(char *)ftemp->GetParName(i);
					 double tpar=ftemp->GetParameter(i);
					 minimizer->SetParameter(i,pars,tpar,1,0,0);
					 if (fix[i]) {
						 minimizer->FixParameter(i);
					 }
				 }
				 minimizer->ExecuteCommand("MIGRAD",0,0);
				 TMinuit *gm=minimizer->GetMinuit();
				 double errdef=sqrt(2.)*TMath::ErfInverse(cper/100.);
				 gm->SetErrorDef(errdef*errdef);
				 double pp=minimizer->GetParameter(par1-1);
				 double eplus,eminus,eparab,gcc;
				 char nn[20];
				 sprintf(nn,"MINOS 500 %d",par1);
				 gm->Command(nn);
				 gm->mnerrs(par1-1,eplus,eminus,eparab,gcc);
				 printf("    Confidence interval (%g%%) : %g ( %g , %g)\n",cper,pp,eminus,eplus);
				 if (ispsf) {
					 histpsfmat->Delete();
				 }
				 delete [] binsh;
				 delete ftemp;
				 delete hh;
				 delete hc;
			}
			else {
				if (fix[par1-1]){
					printf("    Parameter %d is frozen\n",par1);
					break;
				}
				int id[npchain];
				TMath::Sort(npchain,chains[par1-1],id,0);
				int ilow=(int)round(npchain*(0.5-cper/100./2.)-1);
				int ihigh=(int)round(npchain*(0.5+cper/100./2.)-1);
				double pp,eminus,eplus;
				if (islogpar[par1-1]){
					pp=pow(10.,TMath::Median(npchain,chains[par1-1]));
					eminus=pp-pow(10.,chains[par1-1][id[ilow]]);
					eplus=pow(10.,chains[par1-1][id[ihigh]])-pp;
				}
				else {
					pp=TMath::Median(npchain,chains[par1-1]);
					eminus=pp-chains[par1-1][id[ilow]];
					eplus=chains[par1-1][id[ihigh]]-pp;
				}
				printf("    Errors calculated from chains\n");
				printf("    Confidence interval (%g%%) : %g ( %g , %g)\n",cper,pp,eminus,eplus);
			}
		 }while(0);
		 
	 }
     else if (!strcmp(temp,"growth")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 cout << "    Center (1: centroid, 2: sb peak, 3: user input (image coord), 4: user input (J2000.0 coord)) > ";
			 cin >> temp;
			 int center=atoi(temp);
			 if ((center!=1)&&(center!=2)&&(center!=3)&&(center!=4)){
				 printf("    Invalid option %s\n",temp);
				 break;
			 }
			 centroid_ra=0.0;
			 centroid_dec=0.0;
			 if (center==3){
				 cout << "       X center (image coord) > ";
				 cin >> temp;
				 centroid_ra=atof(temp);
				 if ((centroid_ra<1+1e-10)||(centroid_ra>axes[0])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 cout << "       Y center (image coord) > ";
				 cin >> temp;
				 centroid_dec=atof(temp);
				 if ((centroid_dec<1+1e-10)||(centroid_ra>axes[1])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 if (center==4){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
 				 cout << "       RA center (J2000.0) > ";
				 cin >> temp;
				 world[0]=atof(temp);
				 cout << "       Dec center (J2000.0) > ";
				 cin >> temp;
				 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
				 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 
			 cout << "    Bin size (arcsec) > ";
			 cin >> temp;
			 binsize=atof(temp);
			 if (binsize<pixsize*60.*60.){
				 printf("    Error: bin size is smaller than pixel size\n");
				 break;
			 }
			 cout << "    Maximal radius (arcmin) > ";
			 cin >> temp;
			 maxrad=atof(temp);
			 if (maxrad<1e-10){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 double ftot=0.0;
			 double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
			 passmaxexp=maxexp;
			 double *expocor=new double[axes[0]*axes[1]];
			 for (int i=0;i<axes[0]*axes[1];i++){
				 if (exposure[i]>0.25*maxexp){
					 expocor[i]=img[i]*maxexp/exposure[i];
				 }
				 else {
					 expocor[i]=0.0;
				 }
			 }
			 double maximg=TMath::MaxElement(axes[0]*axes[1],expocor);
			 if (center==1){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]>0.0){
							 centroid_ra+=img[j*axes[0]+i]*i;
							 centroid_dec+=img[j*axes[0]+i]*j;
							 ftot+=img[j*axes[0]+i];
						 }
					 }
				 }
				 centroid_ra/=ftot;
				 centroid_dec/=ftot;
				 printf("    Centroid of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 if (center==2){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]==maximg){
							 centroid_ra=i;
							 centroid_dec=j;
						 }
					 }
				 }
				 printf("    Surface-brightness peak of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 int cra=(int)floor(centroid_ra);
			 int cdec=(int)floor(centroid_dec);
			 if (exposure[cdec*axes[0]+cra]==0.0){
				 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
			 }
			 printf("    Exposure at the centre: %g sec\n",maxexp);
			 bin2pix=binsize/pixsize/3600;
			 rad2pix=maxrad/pixsize/60;
			 nbin=(int)floor(rad2pix/bin2pix);
			 nbtot=nbin;
			 growth=NULL;
			 egr=NULL;
			 bins=NULL;
			 ebins=NULL;
			 growth=new double[nbin];
			 egr=new double[nbin];
			 bins=new double[nbin];
			 ebins=new double[nbin];
			 mk_growth_curve(img,exposure,growth,egr,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize);
			 printf("    Total enclosed count rate: %g [counts/sec]\n",growth[nbin-1]);
			 isgrowth=true;
			 delete [] expocor;
		 }
		 while (0);
     }
     else if (!strcmp(temp,"box")){
         do {
             if (!isimg){
                 printf("    Image file not yet loaded\n");
                 break;
             }
             if (!isexp){
                 printf("    Exposure file not yet loaded\n");
                 break;
             }
             cout << "    Starting point (1: cuser input (image coord), 2: user input (FK5)) > ";
             cin >> temp;
             int center=atoi(temp);
             if ((center!=1)&&(center!=2)){
                 printf("    Invalid option %s\n",temp);
                 break;
             }
             centroid_ra=0.0;
             centroid_dec=0.0;
             if (center==1){
                 cout << "       X center (image coord) > ";
                 cin >> temp;
                 centroid_ra=atof(temp);
                 cout << "       Y center (image coord) > ";
                 cin >> temp;
                 centroid_dec=atof(temp);
                 centroid_ra-=1;
                 centroid_dec-=1;
             }
             if (center==2){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
                 cout << "       RA center (J2000.0) > ";
                 cin >> temp;
                 world[0]=atof(temp);
                 cout << "       Dec center (J2000.0) > ";
                 cin >> temp;
                 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
                 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
                 centroid_ra-=1;
                 centroid_dec-=1;
             }
             cout << "    Bin size (arcsec) > ";
             cin >> temp;
             binsize=atof(temp);
             if (binsize<pixsize*60.*60.){
                 printf("    Error: bin size is smaller than pixel size\n");
                 break;
             }
             cout << "    Maximal radius (arcmin) > ";
             cin >> temp;
             maxrad=atof(temp);
             if (maxrad<1e-10){
                 printf("    Invalid value %s\n",temp);
                 break;
             }
             cout << "    Box rotation angle (degrees) > ";
             cin >> temp;
             double tta=atof(temp)-90.;
             if (tta<-90 || tta>270){
                 printf("    Invalid value %s\n",temp);
                 break;
             }
             ellang=tta*TMath::Pi()/180.;
             cout << "    Box width (arcmin) > ";
             cin >> temp;
             double width=atof(temp);
             if (width<=0.0 || width>axes[0] || width>axes[1]){
                 printf("    Invalid value %s\n",temp);
                 break;
             }
             cout << "    Logarithmic binning? (y/n) > ";
             cin >> temp;
             bool islogbin=false;
             if (!strcmp(temp,"y")) {
                 islogbin=true;
             }
             double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
             passmaxexp=maxexp;
             int cra=(int)floor(centroid_ra);
             int cdec=(int)floor(centroid_dec);
             if (exposure[cdec*axes[0]+cra]==0.0){
                 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
             }
             printf("    Exposure at the centre: %g sec\n",maxexp);
             bin2pix=binsize/pixsize/3600;
             rad2pix=maxrad/pixsize/60;
             nbin=(int)floor(rad2pix/bin2pix);
             nbtot=nbin;
             mincounts=0.0;
             cprof=NULL;
             profile=NULL;
             bins=NULL;
             eprof=NULL;
             ebins=NULL;
             cprof=new double[nbin];
             profile=new double[nbin];
             bins=new double[nbin];
             eprof=new double[nbin];
             ebins=new double[nbin];
             effexp=NULL;
             int *numbers=new int[nbin];
             effexp=new double[nbin];
             area=NULL;
             area=new double[nbin];
             if (!issig) {
                 mk_counts_box(img,exposure,cprof,numbers,effexp,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,ellang,width,islogbin);
                 mk_sb_box(img,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,ellang,width,islogbin);
             }
             else {
                 mk_sb_sig(img,sig,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
             }
             
             if ((isback) && (!issig)) {
                 nbin=(int)floor(rad2pix/bin2pix);
                 backprof=NULL;
                 bins=NULL;
                 ebins=NULL;
                 backprof=new double[nbin];
                 bins=new double[nbin];
                 ebins=new double[nbin];
                 backcounts=NULL;
                 backcounts=new double[nbin];
                 int *nnn=new int[nbin];
                 double *eee=new double[nbin];
                 double *ebp=new double[nbin];
                 mk_counts_box(backmap,exposure,backcounts,nnn,eee,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,ellang,width,islogbin);
                 mk_sb_box(backmap,exposure,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,ellang,width,islogbin);
                 cout << "    Do you want to subtract the background profile? (y/n) ";
                 cin >> temp;
                 if (!strcmp(temp,"y")) {
                     for (int i=0; i<nbin; i++) {
                         profile[i]-=backprof[i];
                     }
                 }
                 delete [] nnn;
                 delete [] eee;
                 delete [] ebp;
             }
             for (int i=0; i<nbin; i++) {
                 area[i]=numbers[i]*pixsize*60.*pixsize*60.;
             }
             delete [] numbers;
             isprofile=true;
             ellipse=false;
             isdepr=false;
             sector=false;
             ispsf=false;
         }
         while (0);
     }
     else if (!strcmp(temp,"sector")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 cout << "    Center (1: centroid, 2: sb peak, 3: user input (image coord), 4: user input (J2000.0 coord)) > ";
			 cin >> temp;
			 int center=atoi(temp);
			 if ((center!=1)&&(center!=2)&&(center!=3)&&(center!=4)){
				 printf("    Invalid option %s\n",temp);
				 break;
			 }
			 centroid_ra=0.0;
			 centroid_dec=0.0;
			 if (center==3){
				 cout << "       X center (image coord) > ";
				 cin >> temp;
				 centroid_ra=atof(temp);
				 if ((centroid_ra<1+1e-10)||(centroid_ra>axes[0])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 cout << "       Y center (image coord) > ";
				 cin >> temp;
				 centroid_dec=atof(temp);
				 if ((centroid_dec<1+1e-10)||(centroid_dec>axes[1])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 if (center==4){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
 				 cout << "       RA center (J2000.0) > ";
				 cin >> temp;
				 world[0]=atof(temp);
				 cout << "       Dec center (J2000.0) > ";
				 cin >> temp;
				 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
				 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 
			 cout << "    Bin size (arcsec) > ";
			 cin >> temp;
			 binsize=atof(temp);
			 if (binsize<pixsize*60.*60.){
				 printf("    Error: bin size is smaller than pixel size\n");
				 break;
			 }
			 cout << "    Maximal radius (arcmin) > ";
			 cin >> temp;
			 maxrad=atof(temp);
			 if (maxrad<1e-10){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 cout << "    Lower angle of the sector (deg) > ";
			 cin >> temp;
			 angl=atof(temp);
			 if ((angl<0.0)&&(angl>360.0)){
				 printf("    Error: angle must be between 0 and 360 deg\n");
				 break;
			 }
			 angl*=TMath::Pi()/180.;
			 cout << "    Higher angle of the sector (deg) > ";
			 cin >> temp;
			 angh=atof(temp);
			 if ((angh<0.0)&&(angh>360.0)){
				 printf("    Error: angle must be between 0 and 360 deg\n");
				 break;
			 }
			 angh*=TMath::Pi()/180.;
			 cout << "    Logarithmic binning? (y/n) > ";
			 cin >> temp;
			 bool islogbin=false;
			 if (!strcmp(temp,"y")) {
				 islogbin=true;
			 }
			 double ftot=0.0;
			 double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
			 passmaxexp=maxexp;
			 double *expocor=new double[axes[0]*axes[1]];
			 for (int i=0;i<axes[0]*axes[1];i++){
				 if (exposure[i]>0.25*maxexp){
					 expocor[i]=img[i]*maxexp/exposure[i];
				 }
				 else {
					 expocor[i]=0.0;
				 }
			 }
			 double maximg=TMath::MaxElement(axes[0]*axes[1],expocor);
			 if (center==1){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]>0.0){
							 centroid_ra+=img[j*axes[0]+i]*i;
							 centroid_dec+=img[j*axes[0]+i]*j;
							 ftot+=img[j*axes[0]+i];
						 }
					 }
				 }
				 centroid_ra/=ftot;
				 centroid_dec/=ftot;
				 printf("    Centroid of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 if (center==2){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]==maximg){
							 centroid_ra=i;
							 centroid_dec=j;
						 }
					 }
				 }
				 printf("    Surface-brightness peak of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 int cra=(int)floor(centroid_ra);
			 int cdec=(int)floor(centroid_dec);
			 if (exposure[cdec*axes[0]+cra]==0.0){
				 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
			 }
			 printf("    Exposure at the centre: %g sec\n",maxexp);
			 bin2pix=binsize/pixsize/3600;
			 rad2pix=maxrad/pixsize/60;
			 nbin=(int)floor(rad2pix/bin2pix);				 				 
			 nbtot=nbin;
			 mincounts=0.0;
			 cprof=NULL;
			 profile=NULL;
			 bins=NULL;
			 eprof=NULL;
			 ebins=NULL;
			 cprof=new double[nbin];
			 profile=new double[nbin];
			 bins=new double[nbin];
			 eprof=new double[nbin];
			 ebins=new double[nbin];
			 effexp=NULL;
			 int *numbers=new int[nbin];
			 effexp=new double[nbin];
			 area=NULL;
			 area=new double[nbin];
			 if (!issig) {
				 mk_sector_counts(img,exposure,angl,angh,cprof,numbers,effexp,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
				 mk_sector(img,exposure,angl,angh,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);				 
			 }
			 else {
				 mk_sector_sig(img,sig,exposure,angl,angh,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
			 }
			 
			 if ((isback) && (!issig)) {
				 nbin=(int)floor(rad2pix/bin2pix);
				 backprof=NULL;
				 bins=NULL;
				 ebins=NULL;
				 backprof=new double[nbin];
				 bins=new double[nbin];
				 ebins=new double[nbin];
				 backcounts=NULL;
				 backcounts=new double[nbin];
				 int *nnn=new int[nbin];
				 double *eee=new double[nbin];
				 double *ebp=new double[nbin];
				 mk_sector_counts(backmap,exposure,angl,angh,backcounts,nnn,eee,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
				 mk_sector(backmap,exposure,angl,angh,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
				 cout << "    Do you want to subtract the background profile? (y/n) ";
				 cin >> temp;
				 if (!strcmp(temp,"y")) {
					 for (int i=0; i<nbin; i++) {
						 profile[i]-=backprof[i];
					 }
				 }
				 delete [] ebp;
				 delete [] eee;
				 delete [] nnn;
			 }
			 for (int i=0; i<nbin; i++) {
				 area[i]=numbers[i]*pixsize*60.*pixsize*60.;
			 }
			 delete [] numbers;
			 delete [] expocor;
			 isprofile=true;
			 sector=true;
			 isdepr=false;
			 ellipse=false;
			 ispsf=false;
		 }
		 while (0);
     }
     else if (!strcmp(temp,"sectorellipse")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 cout << "    Center (1: centroid, 2: sb peak, 3: user input (image coord), 4: user input (J2000.0 coord)) > ";
			 cin >> temp;
			 int center=atoi(temp);
			 if ((center!=1)&&(center!=2)&&(center!=3)&&(center!=4)){
				 printf("    Invalid option %s\n",temp);
				 break;
			 }
			 centroid_ra=0.0;
			 centroid_dec=0.0;
			 if (center==3){
				 cout << "       X center (image coord) > ";
				 cin >> temp;
				 centroid_ra=atof(temp);
				 if ((centroid_ra<1+1e-10)||(centroid_ra>axes[0])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 cout << "       Y center (image coord) > ";
				 cin >> temp;
				 centroid_dec=atof(temp);
				 if ((centroid_dec<1+1e-10)||(centroid_dec>axes[1])){
					 printf("    Invalid value %s\n",temp);
					 break;
				 }
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 if (center==4){
                 double pixcrd[2],imgcrd[2],world[2];
                 double inra,indec;
                 int stat;
 				 cout << "       RA center (J2000.0) > ";
				 cin >> temp;
				 world[0]=atof(temp);
				 cout << "       Dec center (J2000.0) > ";
				 cin >> temp;
				 world[1]=atof(temp);
                 status = wcss2p(wcs_inp, 1, 2, world, &inra,&indec,imgcrd, pixcrd,&stat);
                 centroid_ra=pixcrd[0];
                 centroid_dec=pixcrd[1];
				 printf("       Corresponding pixel coordinates : %g %g\n",centroid_ra,centroid_dec);
				 centroid_ra-=1;
				 centroid_dec-=1;
			 }
			 
			 cout << "    Bin size (arcsec) > ";
			 cin >> temp;
			 binsize=atof(temp);
			 if (binsize<pixsize*60.*60.){
				 printf("    Error: bin size is smaller than pixel size\n");
				 break;
			 }
			 cout << "    Maximal radius (arcmin) > ";
			 cin >> temp;
			 maxrad=atof(temp);
			 if (maxrad<1e-10){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 cout << "    Lower angle of the sector (deg) > ";
			 cin >> temp;
			 angl=atof(temp);
			 if ((angl<0.0)||(angl>360.0)){
				 printf("    Error: angle must be between 0 and 360 deg\n");
				 break;
			 }
			 angl*=TMath::Pi()/180.;
			 cout << "    Higher angle of the sector (deg) > ";
			 cin >> temp;
			 angh=atof(temp);
			 if ((angh<0.0)||(angh>360.0)){
				 printf("    Error: angle must be between 0 and 360 deg\n");
				 break;
			 }
			 angh*=TMath::Pi()/180.;
			 cout << "    Angle between RA axis and major axis (degrees) > ";
			 cin >> temp;
			 double tta=atof(temp)-90.;
			 if (tta<-90 || tta>270){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 ellang=tta*TMath::Pi()/180.;
			 cout << "    Ratio between major and minor axis > ";
			 cin >> temp;
			 aoverb=atof(temp);
			 if (aoverb<1.0){
				 printf("    Invalid value %s\n",temp);
				 break;
			 }
			 else if (aoverb==1.0) {
				 printf("    You could have just as well used the \"sector\" command...\n");
			 }
			 cout << "    Logarithmic binning? (y/n) > ";
			 cin >> temp;
			 bool islogbin=false;
			 if (!strcmp(temp,"y")) {
				 islogbin=true;
			 }
			 double ftot=0.0;
			 double maxexp=TMath::MaxElement(axes[0]*axes[1],exposure);
			 passmaxexp=maxexp;
			 double *expocor=new double[axes[0]*axes[1]];
			 for (int i=0;i<axes[0]*axes[1];i++){
				 if (exposure[i]>0.25*maxexp){
					 expocor[i]=img[i]*maxexp/exposure[i];
				 }
				 else {
					 expocor[i]=0.0;
				 }
			 }
			 double maximg=TMath::MaxElement(axes[0]*axes[1],expocor);
			 if (center==1){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]>0.0){
							 centroid_ra+=img[j*axes[0]+i]*i;
							 centroid_dec+=img[j*axes[0]+i]*j;
							 ftot+=img[j*axes[0]+i];
						 }
					 }
				 }
				 centroid_ra/=ftot;
				 centroid_dec/=ftot;
				 printf("    Centroid of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 if (center==2){
				 for (int i=0;i<axes[0];i++){
					 for (int j=0;j<axes[1];j++){
						 if (expocor[j*axes[0]+i]==maximg){
							 centroid_ra=i;
							 centroid_dec=j;
						 }
					 }
				 }
				 printf("    Surface-brightness peak of the cluster (image coordinates): %g %g\n",centroid_ra+1,centroid_dec+1);
			 }
			 int cra=(int)floor(centroid_ra);
			 int cdec=(int)floor(centroid_dec);
			 if (exposure[cdec*axes[0]+cra]==0.0){
				 printf("    WARNING: Exposure is 0 at the central position. The fitting procedure may fail.\n");
			 }
			 printf("    Exposure at the centre: %g sec\n",maxexp);
			 bin2pix=binsize/pixsize/3600;
			 rad2pix=maxrad/pixsize/60;
			 nbin=(int)floor(rad2pix/bin2pix);				 				 
			 nbtot=nbin;
			 mincounts=0.0;
			 cprof=NULL;
			 profile=NULL;
			 bins=NULL;
			 eprof=NULL;
			 ebins=NULL;
			 cprof=new double[nbin];
			 profile=new double[nbin];
			 bins=new double[nbin];
			 eprof=new double[nbin];
			 ebins=new double[nbin];
			 effexp=NULL;
			 int *numbers=new int[nbin];
			 effexp=new double[nbin];
			 area=NULL;
			 area=new double[nbin];
			 if (!issig) {
				 mk_sectell_counts(img,exposure,angl,angh,cprof,numbers,effexp,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 mk_ellipse_sector(img,exposure,angl,angh,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
			 }
			 else {
				 mk_sector_sig(img,sig,exposure,angl,angh,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin);
			 }
			 
			 if ((isback) && (!issig)) {
				 nbin=(int)floor(rad2pix/bin2pix);
				 backprof=NULL;
				 bins=NULL;
				 ebins=NULL;
				 backprof=new double[nbin];
				 bins=new double[nbin];
				 ebins=new double[nbin];
				 double *ebp=new double[nbin];
				 backcounts=NULL;
				 backcounts=new double[nbin];
				 int *nnn=new int[nbin];
				 double *eee=new double[nbin];
				 mk_sectell_counts(backmap,exposure,angl,angh,backcounts,nnn,eee,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 mk_ellipse_sector(backmap,exposure,angl,angh,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,maxrad,binsize,islogbin);
				 cout << "    Do you want to subtract the background profile? (y/n) ";
				 cin >> temp;
				 if (!strcmp(temp,"y")) {
					 for (int i=0; i<nbin; i++) {
						 profile[i]-=backprof[i];
					 }
				 }
				 delete [] ebp;
				 delete [] eee;
				 delete [] nnn;
			 }
			 for (int i=0; i<nbin; i++) {
				 area[i]=numbers[i]*pixsize*60.*pixsize*60.;
			 }
			 delete [] numbers;
			 delete [] expocor;
			 isprofile=true;
			 sector=true;
			 ellipse=true;
			 isdepr=false;
			 ispsf=false;
		 }
		 while (0);
     }
	 else if (!strcmp(temp,"scatter")){
		 do {
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod) {
				 printf("    No model loaded\n");
				 break;
			 }
			 double back,eback;
			 int npar=model->GetNpar();
			 char *pname=new char[20];
			 for (int i=0; i<npar; i++) {
				 pname=(char *)model->GetParName(i);
				 if (!strcmp(pname,"const")) {
					 back=model->GetParameter(i);
					 eback=model->GetParError(i);
				 }
			 }
			 cout << "    Number of sectors > ";
			 cin >> temp;
			 int nsect=atoi(temp);
			 if (nsect<2) {
				 printf("    Not enough sectors\n");
				 break;
			 }
			 double angles[nsect+1];
			 angles[0]=0.0;
			 for (int i=1; i<nsect+1; i++) {
				 double ai=2.*TMath::Pi()/nsect*i;
				 angles[i]=ai;
			 }			 
			 cout << "    Subtract the background? (y/n) > ";
			 cin >> temp;
			 double **allprofs=new double*[nsect];
			 double **allerrs=new double*[nsect];
			 for (int i=0; i<nsect; i++) {
				 allprofs[i]=new double[nbin];
				 allerrs[i]=new double [nbin];
			 }
			 for (int i=0; i<nsect; i++) {
				 double angl=angles[i];
				 double angh=angles[i+1];
				 mk_sector_scat(img,exposure,angl,angh,bins,ebins,nbin,allprofs[i],allerrs[i],axes,centroid_ra,centroid_dec,pixsize,maxrad);
				 if (!strcmp(temp,"y")) {
					 if (isback) {
						 double *backsect=new double[nbin];
						 double *dummyerr=new double[nbin];
						 mk_sector_scat(backmap,exposure,angl,angh,bins,ebins,nbin,backsect,dummyerr,axes,centroid_ra,centroid_dec,pixsize,maxrad);
						 for (int nb=0; nb<nbin; nb++) {
							 allprofs[i][nb]-=backsect[nb];
						 }
						 delete [] backsect;
						 delete [] dummyerr;
					 }
					 backsub(back,eback,nbin,allprofs[i],allerrs[i]);
				 }
			 }
			 scat=NULL;
			 escat=NULL;
			 scat=new double[nbin];
			 escat=new double[nbin];
			 TRandom3 *g=new TRandom3(0);
			 int nsim=1e3;
			 for (int i=0; i<nbin; i++) {
				 double statscat=0.0;
				 for (int ns=0; ns<nsect; ns++) {
					 statscat+=allerrs[ns][i]*allerrs[ns][i]/profile[i]/profile[i];					 
				 }
				 statscat/=nsect;
				 double *vals=new double[nsim];
				 for (int sim=0; sim<nsim; sim++) {
					 double mv=g->Gaus(profile[i],eprof[i]);
					 double ts=0.0;
					 double statscat=0.0;
					 for (int ns=0; ns<nsect; ns++) {
						 double vns=g->Gaus(allprofs[ns][i],allerrs[ns][i]);
						 ts+=(vns-mv)*(vns-mv)/mv/mv;
					 }
					 if (ts/nsect>statscat) {
						 vals[sim]=sqrt(ts/nsect-statscat);						 
					 }
					 else {
						 vals[sim]=0.0;
					 }
				 }
				 double sc=0.0;
				 for (int ns=0; ns<nsect; ns++) {
					 sc+=(allprofs[ns][i]-profile[i])*(allprofs[ns][i]-profile[i])/profile[i]/profile[i];					 
				 }
				 if (sc/nsect>statscat) {
					 scat[i]=sqrt(sc/nsect-statscat);					 
				 }
				 else {
					 scat[i]=0.0;
				 }
				 escat[i]=TMath::RMS(nsim,vals);
				 delete [] vals;
			 }
			 
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F* hh=new TH1F("hh","hh",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 hh->SetBinContent(i+1,scat[i]);
				 hh->SetBinError(i+1,escat[i]);
			 }
			 hh->SetTitle("");
			 hh->SetStats(false);
			 hh->SetXTitle("Radius [arcmin]");
			 if (isr500){
				 hh->SetXTitle("r/r_{500}");
			 }
			 if (isr200){
				 hh->SetXTitle("r/r_{200}");
			 }
             if (iskpc) {
                 hh->SetXTitle("Radius [kpc]");
             }
			 if (!logy){
				 hh->GetYaxis()->SetTitleOffset(1.8);
			 }
			 else {
				 hh->GetYaxis()->SetTitleOffset(1.5);
			 }
			 hh->GetXaxis()->SetTitleOffset(1.2);
			 hh->SetYTitle("Scatter");
			 hh->Draw();			 
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
			 delete [] binsh;
			 for (int i=0; i<nsect; i++) {
				 delete [] allprofs[i];
				 delete [] allerrs[i];
			 }
			 delete [] allprofs;
			 delete [] allerrs;
			 delete hh;
			 isscat=true;
			 
		 }while (0);
	 }
	 else if (!strcmp(temp,"allsectors")){
		 do {
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod) {
				 printf("    No model loaded\n");
				 break;
			 }
			 double back,eback;
			 int npar=model->GetNpar();
			 char *pname=new char[20];
			 for (int i=0; i<npar; i++) {
				 pname=(char *)model->GetParName(i);
				 if (!strcmp(pname,"const")) {
					 back=model->GetParameter(i);
					 eback=model->GetParError(i);
				 }
			 }
			 cout << "    Number of sectors > ";
			 cin >> temp;
			 int nsect=atoi(temp);
			 if (nsect<2) {
				 printf("    Not enough sectors\n");
				 break;
			 }
			 double angles[nsect+1];
			 angles[0]=0.0;
			 for (int i=1; i<nsect+1; i++) {
				 double ai=2.*TMath::Pi()/nsect*i;
				 angles[i]=ai;
			 }			 
			 cout << "    Subtract the background? (y/n) > ";
			 cin >> temp;
			 double **allprofs=new double*[nsect];
			 double **allerrs=new double*[nsect];
			 for (int i=0; i<nsect; i++) {
				 allprofs[i]=new double[nbin];
				 allerrs[i]=new double [nbin];
			 }
			 for (int i=0; i<nsect; i++) {
				 double angl=angles[i];
				 double angh=angles[i+1];
				 mk_sector_scat(img,exposure,angl,angh,bins,ebins,nbin,allprofs[i],allerrs[i],axes,centroid_ra,centroid_dec,pixsize,maxrad);
				 if (!strcmp(temp,"y")) {
					 if (isback) {
						 double *backsect=new double[nbin];
						 double *dummyerr=new double[nbin];
						 mk_sector_scat(backmap,exposure,angl,angh,bins,ebins,nbin,backsect,dummyerr,axes,centroid_ra,centroid_dec,pixsize,maxrad);
						 for (int nb=0; nb<nbin; nb++) {
							 allprofs[i][nb]-=backsect[nb];
						 }
						 delete [] backsect;
						 delete [] dummyerr;
					 }
					 backsub(back,eback,nbin,allprofs[i],allerrs[i]);
				 }
			 }
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F **hh=new TH1F*[nsect];
			 char hname[100];
			 for (int i=0; i<nsect; i++) {
				 sprintf(hname,"h%d",i+1);
				 hh[i]=new TH1F(hname,hname,nbin,binsh);				 
				 for (int nb=0;nb<nbin;nb++){
					 hh[i]->SetBinContent(nb+1,allprofs[i][nb]);
					 hh[i]->SetBinError(nb+1,allerrs[i][nb]);
				 }
				 hh[i]->SetLineColor(i+1);
				 hh[i]->SetMarkerColor(i+1);
				 hh[i]->SetStats(0);
			 }
			 hh[0]->SetTitle("");
			 hh[0]->SetXTitle("Radius [arcmin]");
			 if (isr500){
				 hh[0]->SetXTitle("r/r_{500}");
			 }
			 if (isr200){
				 hh[0]->SetXTitle("r/r_{200}");
			 }
             if (iskpc) {
                 hh[0]->SetXTitle("Radius [kpc]");
             }
			 if (!logy){
				 hh[0]->GetYaxis()->SetTitleOffset(1.8);
			 }
			 else {
				 hh[0]->GetYaxis()->SetTitleOffset(1.5);
			 }
			 hh[0]->GetXaxis()->SetTitleOffset(1.2);
			 hh[0]->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
			 hh[0]->Draw();
			 for (int i=1; i<nsect; i++) {
				 hh[i]->Draw("same");
			 }
             TLegend* l=new TLegend(0.65,0.7,0.85,0.9);
             char legtext[200];
             for (int i=0; i<nsect; i++) {
 				 double angl=angles[i]*180./TMath::Pi();
				 double angh=angles[i+1]*180./TMath::Pi();
                 sprintf(legtext,"Sector %g-%g",angl,angh);
				 sprintf(hname,"h%d",i+1);
                 l->AddEntry(hname,legtext,"L");
             }
             l->SetHeader("All sectors");
             l->SetFillColor(kWhite);
             l->Draw();
			 FILE *fall=fopen("allsectors.txt","w");
			 for (int i=0; i<nbin; i++) {
				 fprintf(fall,"%g  %g ",binsh[i],binsh[i+1]);
				 for (int ns=0; ns<nsect-1; ns++) {
					 fprintf(fall," %g  %g ",allprofs[ns][i],allerrs[ns][i]);
				 }
				 fprintf(fall," %g  %g\n",allprofs[nsect-1][i],allerrs[nsect-1][i]);
			 }
			 fclose(fall);
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
			 delete [] binsh;
			 for (int i=0; i<nsect; i++) {
				 delete [] allprofs[i];
				 delete [] allerrs[i];
				 delete hh[i];
			 }
			 delete [] allprofs;
			 delete [] allerrs;
			 delete hh;
			 isscat=true;
			 
		 }while (0);
	 }
      else if (!strcmp(temp,"mediansb")){
          //Median in sectors
          do {
              if (!ismod) {
                  printf("    No model loaded\n");
                  break;
              }
              double back,eback;
              int npar=model->GetNpar();
              char *pname=new char[20];
              for (int i=0; i<npar; i++) {
                  pname=(char *)model->GetParName(i);
                  if (!strcmp(pname,"const")) {
                      back=model->GetParameter(i);
                      eback=model->GetParError(i);
                  }
              }
              medprof=NULL;
              emedprof=NULL;
              medprof=new double[nbin];
              emedprof=new double[nbin];
              if (!isvoronoi){
                  cout << "    Number of sectors > ";
                  cin >> temp;
                  int nsect=atoi(temp);
                  if (nsect<2) {
                      printf("    Not enough sectors\n");
                      break;
                  }
                  double angles[nsect+1];
                  angles[0]=0.0;
                  for (int i=1; i<nsect+1; i++) {
                      double ai=2.*TMath::Pi()/nsect*i;
                      angles[i]=ai;
                  }
                  cout << "    Subtract the background? (y/n) > ";
                  cin >> temp;
                  double **allprofs=new double*[nsect];
                  double **allerrs=new double*[nsect];
                  for (int i=0; i<nsect; i++) {
                      allprofs[i]=new double[nbin];
                      allerrs[i]=new double [nbin];
                  }
                  for (int i=0; i<nsect; i++) {
                      double angl=angles[i];
                      double angh=angles[i+1];
                      mk_sector_scat(img,exposure,angl,angh,bins,ebins,nbin,allprofs[i],allerrs[i],axes,centroid_ra,centroid_dec,pixsize,maxrad);
                      if (!strcmp(temp,"y")) {
                          if (isback) {
                              double *backsect=new double[nbin];
                              double *dummyerr=new double[nbin];
                              mk_sector_scat(backmap,exposure,angl,angh,bins,ebins,nbin,backsect,dummyerr,axes,centroid_ra,centroid_dec,pixsize,maxrad);
                              for (int nb=0; nb<nbin; nb++) {
                                  allprofs[i][nb]-=backsect[nb];
                              }
                              delete [] backsect;
                              delete [] dummyerr;
                          }
                          backsub(back,eback,nbin,allprofs[i],allerrs[i]);
                      }
                  }
                  for (int i=0; i<nbin; i++) {
                      int nsectobs=0;
                      for (int ns=0; ns<nsect; ns++) {
                          if (allerrs[ns][i]>0.0) {
                              nsectobs++;
                          }
                      }
                      double *vals=new double[nsectobs];
                      double *evals=new double[nsectobs];
                      int is=0;
                      for (int ns=0; ns<nsect; ns++) {
                          if (allerrs[ns][i]>0.0) {
                              vals[is]=allprofs[ns][i];
                              evals[is]=allerrs[ns][i];
                              is++;
                          }
                      }
                      medprof[i]=TMath::Median(nsectobs,vals);
                      emedprof[i]=medianerr(nsectobs,1e3,vals,evals);
                      delete [] vals;
                      delete [] evals;
                      for (int i=0; i<nsect; i++) {
                          delete [] allprofs[i];
                          delete [] allerrs[i];
                      }
                      delete [] allprofs;
                      delete [] allerrs;
                  }
              }
              else {
                  mk_median_prof(voronoimap,voronoierr,medprof,emedprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad);
              }
              double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
              TH1F *hhmed=new TH1F("median","median",nbin,binsh);
              for (int i=0;i<nbin;i++){
                  hhmed->SetBinContent(i+1,medprof[i]);
                  hhmed->SetBinError(i+1,emedprof[i]);
              }
              c1->Clear();
              if (!logy){
                  c1->SetLeftMargin(0.15);
              }
              else {
                  c1->SetLeftMargin(0.12);
              }
              c1->SetBottomMargin(0.12);
              hhmed->SetTitle("");
              hhmed->SetStats(false);
              hhmed->SetXTitle("Radius [arcmin]");
              if (isr500) {
                  hhmed->SetXTitle("r/r_{500}");
              }
              if (isr200) {
                  hhmed->SetXTitle("r/r_{200}");
              }
              if (iskpc) {
                  hhmed->SetXTitle("Radius [kpc]");
              }
              if (!logy){
                  hhmed->GetYaxis()->SetTitleOffset(1.8);
              }
              else {
                  hhmed->GetYaxis()->SetTitleOffset(1.5);
              }
              hhmed->GetXaxis()->SetTitleOffset(1.2);
              hhmed->GetXaxis()->CenterTitle();
              hhmed->GetYaxis()->CenterTitle();
              hhmed->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
              hhmed->Draw();
              c1->Update();
              ismed=true;
              if (!scripting) {
                  theApp.Run(kTRUE);
              }
              delete [] binsh;
              delete hhmed;
          }
          while (0);
      }
	 /*else if (!strcmp(temp,"deviations")){
		 do {
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod) {
				 printf("    No model loaded\n");
				 break;
			 }
			 double back,eback;
			 int npar=model->GetNpar();
			 char *pname=new char[20];
			 for (int i=0; i<npar; i++) {
				 pname=(char *)model->GetParName(i);
				 if (!strcmp(pname,"const")) {
					 back=model->GetParameter(i);
					 eback=model->GetParError(i);
				 }
			 }
			 cout << "    Number of sectors > ";
			 cin >> temp;
			 int nsect=atoi(temp);
			 if (nsect<2) {
				 printf("    Not enough sectors\n");
				 break;
			 }
			 double angles[nsect+1];
			 angles[0]=0.0;
			 for (int i=1; i<nsect+1; i++) {
				 double ai=2.*TMath::Pi()/nsect*i;
				 angles[i]=ai;
			 }			 
			 cout << "    Subtract the background? (y/n) > ";
			 cin >> temp;
			 double **allprofs=new double*[nsect];
			 double **allerrs=new double*[nsect];
			 for (int i=0; i<nsect; i++) {
				 allprofs[i]=new double[nbin];
				 allerrs[i]=new double [nbin];
			 }
			 for (int i=0; i<nsect; i++) {
				 double angl=angles[i];
				 double angh=angles[i+1];
				 mk_sector_scat(img,exposure,angl,angh,bins,ebins,nbin,allprofs[i],allerrs[i],axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize);				 
				 if (!strcmp(temp,"y")) {
					 if (isback) {
						 double *backsect=new double[nbin];
						 double *dummyerr=new double[nbin];
						 mk_sector_scat(backmap,exposure,angl,angh,bins,ebins,nbin,backsect,dummyerr,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize);
						 for (int nb=0; nb<nbin; nb++) {
							 allprofs[i][nb]-=backsect[nb];
						 }
						 delete [] backsect;
						 delete [] dummyerr;
					 }
					 backsub(back,eback,nbin,allprofs[i],allerrs[i]);
				 }
			 }
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH2F *hdev=new TH2F("hdev","",nsect,0.,360.,nbin,binsh);
             for (int i=0; i<nsect; i++) {
                 for (int nb=0; nb<nbin; nb++) {
                     double toterr=sqrt(eprof[nb]*eprof[nb]+allerrs[i][nb]*allerrs[i][nb]);
                     double dev=(allprofs[i][nb]-profile[nb])/toterr;
                     hdev->SetBinContent(i+1,nb+1,dev);
                 }
             }
             hdev->GetXaxis()->SetTitle("Angle [deg]");
             hdev->Draw("cont4 pol");
             c1->SetLogx(1);
             c1->SetLogy(1);
			 c1->Update();
			 theApp.Run(kTRUE);
			 delete [] binsh;
			 delete [] allprofs;
			 delete [] allerrs;
			 delete hdev;			 
		 }while (0);
	 }*/
     else if (!strcmp(temp,"plot")){
		 do{
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 m1=NULL;
			 c1->Clear();
			 c1->SetTickx();
			 c1->SetTicky();
			 if (!logy){
				 c1->SetLeftMargin(0.15);
			 }
			 else {
				 c1->SetLeftMargin(0.12);
			 }
			 c1->SetBottomMargin(0.12);
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F* hh=new TH1F("hh","hh",nbin,binsh);
			 for (int i=0;i<nbin;i++){
			   hh->SetBinContent(i+1,profile[i]);
			   hh->SetBinError(i+1,eprof[i]);
			 }
			 hh->SetTitle("");
			 hh->SetStats(false);
			 hh->SetXTitle("Radius [arcmin]");
			 if (isr500){
				 hh->SetXTitle("r/r_{500}");
			 }
			 if (isr200){
				 hh->SetXTitle("r/r_{200}");
			 }
             if (iskpc) {
                 hh->SetXTitle("Radius [kpc]");
             }
			 if (!logy){
				 hh->GetYaxis()->SetTitleOffset(1.8);
			 }
			 else {
				 hh->GetYaxis()->SetTitleOffset(1.5);
			 }
			 hh->GetXaxis()->SetTitleOffset(1.2);
             hh->GetXaxis()->CenterTitle();
             hh->GetYaxis()->CenterTitle();
			 hh->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
			 hh->Draw();
			 TH1F *hback=NULL;
			 if (isback) {
				 hback=(TH1F*)hh->Clone("hback");
				 for (int i=0; i<nbin; i++) {
					 hback->SetBinContent(i+1,backprof[i]);
					 hback->SetBinError(i+1,0.0);
				 }
				 hback->SetLineColor(kBlue);
				 double minback=hback->GetMinimum();
				 double minprof=hh->GetMinimum();
				 double maxprof=hh->GetMaximum();
				 if (minback>0.0 && minback<minprof) {
					 hh->GetYaxis()->SetRangeUser(minback/1.2,maxprof*1.2);					 
				 }
				 hh->Draw();
				 hback->Draw("same");
			 }
             if (ismod) {
                 model->Draw("same");
             }
			 if (isnested){
				TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
				for (int k=0;k<model->GetNpar();k++){
					double tpar=model->GetParameter(k);
					ftemp->SetParameter(k,tpar);
				}
				double *vnest=new double[nbin];
				double *evnest=new double[nbin];
				calc_envelope(chains,npchain,ftemp,bins,nbin,fix,islogpar,vnest,evnest);
				double *zer=new double[nbin];
				for (int i=0;i<nbin;i++){
					zer[i]=0.0;
					//printf("bins, vnest, evnest: %g  %g  %g\n",bins[i],vnest[i],evnest[i]);
				}
				TGraphErrors *gnest=new TGraphErrors(nbin,bins,vnest,zer,evnest);
				gnest->SetFillColor(2);
				gnest->Draw("3");
				ftemp->Delete();
			 }
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
             else {
		     c1->SaveAs("plot.pdf");
	     }
             delete [] binsh;
			 delete hh;
			 delete hback;
       }while(0);
     }
     else if (!strcmp(temp,"plotcounts")){
		 do{
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 m1=NULL;
			 c1->Clear();
			 c1->SetTickx();
			 c1->SetTicky();
			 c1->SetLeftMargin(0.12);
			 c1->SetBottomMargin(0.12);
			 ecp=NULL;
			 ecp=new double[nbin];
			 for (int i=0;i<nbin;i++){
				 ecp[i]=sqrt(cprof[i]);
			 }
			 TGraphErrors *gg=new TGraphErrors(nbin,bins,cprof,ebins,ecp);
			 TMultiGraph *m1=new TMultiGraph("m1","m1");
			 m1->Add(gg);
			 m1->Draw("AP");
			 m1->SetTitle("");
			 m1->GetXaxis()->SetTitleOffset(1.2);
			 m1->GetXaxis()->SetTitle("Radius [arcmin]");
			 m1->GetYaxis()->SetTitleOffset(1.5);
			 m1->GetYaxis()->SetTitle("Counts profile");
			 m1->SetTitle("");
			 m1->Draw("AP");
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
             else {
		     c1->SaveAs("plotcounts.pdf");
	     }
		 }while(0);
	 }
	 else if (!strcmp(temp,"plotgrowth")){
		 do{
			 if (!isgrowth){
				 printf("    Growth curve not yet extracted\n");
				 break;
			 }
			 m1=NULL;
			 c1->Clear();
			 c1->SetTickx();
			 c1->SetTicky();
			 c1->SetLeftMargin(0.14);
			 c1->SetBottomMargin(0.12);
			 TGraphErrors *gg=new TGraphErrors(nbin,bins,growth,ebins,egr);
			 TMultiGraph *m1=new TMultiGraph("m1","m1");
			 m1->Add(gg);
			 m1->Draw("AP");
			 m1->SetTitle("");
			 m1->GetXaxis()->SetTitleOffset(1.2);
			 m1->GetXaxis()->SetTitle("Radius [arcmin]");
			 m1->GetYaxis()->SetTitleOffset(1.5);
			 m1->GetYaxis()->SetTitle("Enclosed brightness [counts/sec]");
			 m1->SetTitle("");
			 m1->Draw("AP");
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
             else {
		     c1->SaveAs("plotgrowth.pdf");
	     }
		 }while(0);
	 }
	 else if (!strcmp(temp,"plotgrmod")){
		 do{
			 if (!isgrowth){
				 printf("    Growth curve not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model yet loaded\n");
				 break;
			 }
			 char *pars=new char[200];
			 int npp=model->GetNpar();
			 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			 for (int k=0;k<npp;k++){
				 pars=(char *)model->GetParName(k);
				 if ((!strcmp(pars,"const"))&&(strcmp(modname,"const"))){
					 ftemp->SetParameter(k,0.0);
					 ftemp->SetParName(k,"const");
				 }
				 else {
					 double tpar=model->GetParameter(k);
					 ftemp->SetParameter(k,tpar);
				 }
			 }
			 m1=NULL;
			 c1->Clear();
			 c1->SetTickx();
			 c1->SetTicky();
			 c1->SetLeftMargin(0.14);
			 c1->SetBottomMargin(0.12);
			 TGraphErrors *gg=new TGraphErrors(nbin,bins,growth,ebins,egr);
			 TMultiGraph *m1=new TMultiGraph("m1","m1");
			 m1->Add(gg,"P");
			 m1->Draw("AP");
			 m1->SetTitle("");
			 m1->GetXaxis()->SetTitleOffset(1.2);
			 m1->GetXaxis()->SetTitle("Radius [arcmin]");
			 m1->GetYaxis()->SetTitleOffset(1.5);
			 m1->GetYaxis()->SetTitle("Enclosed brightness [counts/sec]");
			 m1->SetTitle("");
			 m1->Draw("AP");
			 double *modimg=new double[axes[0]*axes[1]];
			 double *modgr=new double[nbin];
			 double *emp=new double[nbin];
			 mk_mod_img(ftemp,exposure,modimg,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse);
			 mk_growth_curve(modimg,exposure,modgr,emp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize);
			 TGraph *gm=new TGraph(nbin,bins,modgr);
			 gm->SetLineColor(kBlue);
			 gm->SetLineWidth(2);
			 m1->Add(gm,"CP");
			 m1->Draw("A");	
			 c1->Update();
             if (!scripting) {
                 theApp.Run(kTRUE);
             }
			 delete [] modimg;
			 delete [] modgr;
			 delete [] emp;
			 delete ftemp;
		 }while(0);
	 }
	 else if (!strcmp(temp,"plotmod")){
		do{
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			c1->Clear();
			c1->SetTickx();
			c1->SetTicky();
			if (!logy){
				c1->SetLeftMargin(0.15);
			}
			else {
				c1->SetLeftMargin(0.12);
			}
			c1->SetBottomMargin(0.12);
			model->Draw();
			model->GetXaxis()->SetTitleOffset(1.5);
			model->GetXaxis()->SetTitle("Radius [arcmin]");
			if (!logy){
				model->GetYaxis()->SetTitleOffset(1.8);
			}
			else {
				model->GetYaxis()->SetTitleOffset(1.8);
			}
			model->GetYaxis()->SetTitle("Model profile [counts/sec/arcmin^{2}]");
			c1->Update();
            if (!scripting) {
                theApp.Run(kTRUE);
            }
		}while(0);
    }
	else if (!strcmp(temp,"plotmodcounts")){
		do{
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			if (!isexp){
				printf("    Exposure file not yet loaded\n");
				break;
			}	
			if (!isprofile){
				printf("    Profile not yet extracted\n");
				break;
			}
			int npp=model->GetNpar();
			TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			for (int k=0;k<npp;k++){
				double tpar=model->GetParameter(k);
				ftemp->SetParameter(k,tpar);
			}
			m1=NULL;
			c1->Clear();
			c1->SetLeftMargin(0.12);
			c1->SetBottomMargin(0.12);
			double *modimg=new double[axes[0]*axes[1]];
			double *modprof=new double[nbin];
			int *nnn=new int[nbin];
			mk_mod_img(ftemp,exposure,modimg,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse);
			mk_modprof(modimg,exposure,bins,ebins,modprof,nbin,sector,angh,angl,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,mincounts);
			double *dummyerr=new double[nbin];
			for (int i=0;i<nbin;i++){
				dummyerr[i]=0.0;
			}
			TGraphErrors *gg=new TGraphErrors(nbin,bins,modprof,ebins,dummyerr);
			TMultiGraph *m1=new TMultiGraph("m1","m1");
			//gg->SetMarkerColor(kBlue);
			m1->Add(gg);
			m1->Draw("AP");
			m1->SetTitle("");
			m1->GetXaxis()->SetTitleOffset(1.2);
			m1->GetXaxis()->SetTitle("Radius [arcmin]");
			m1->GetYaxis()->SetTitleOffset(1.5);
			m1->GetYaxis()->SetTitle("Model counts profile");
			m1->SetTitle("");
			m1->Draw("AP");
			c1->Update();
            if (!scripting) {
                theApp.Run(kTRUE);
            }
			delete [] modimg;
			delete [] modprof;
			delete [] nnn;
			delete [] dummyerr;
			delete ftemp;
		}while(0);
    }
    else if (!strcmp(temp,"logx")){
		cout << "    Log scale on X axis? (y/n) > ";
		cin >> temp;
		if (!strcmp(temp,"y")) logx=true;
		if (!strcmp(temp,"n")) logx=false;
		if (logx){
			c1->SetLogx();
		}
		else {
			c1->SetLogx(0);
		}
	}
	else if (!strcmp(temp,"logy")){
		cout << "    Log scale on Y axis? (y/n) > ";
		cin >> temp;
		if (!strcmp(temp,"y")) logy=true;
		if (!strcmp(temp,"n")) logy=false;
		if (logy){
			c1->SetLogy();
		}
		else {
			c1->SetLogy(0);
		}
	}
    else if (!strcmp(temp,"showmod")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			model->Print();
		}while(0);
    }
    else if (!strcmp(temp,"region")){
		do {
			if (!isexp){
				printf("    Exposure file not yet loaded\n");
				break;
			}
			cout << "    Region file > ";
			cin >> temp;
            region(temp,exposure,axes);
		}while(0);
    }
    else if (!strcmp(temp,"group")){
		do {
			if (!isprofile){
				printf("    Profile not yet extracted\n");
				break;
			}
			cout << "    What do you want to group (minimum counts: 1, minimum S/N: 2) > ";
			cin >> temp;
			int nn=atoi(temp);
			if (nn==1){
				cout << "    Minimum counts per bin > ";
				cin >> temp;
				mincounts=atof(temp);
				if (!isback){
					mk_grouping(cprof,bins,ebins,profile,eprof,area,effexp,nbin,mincounts);
				}
				else {
					mk_group_isback(cprof,bins,ebins,profile,eprof,area,effexp,nbin,mincounts,backprof,backcounts);
				}
			}
			else if (nn==2) {
				cout << "    Minimum S/N > ";
				cin >> temp;
				double minsn=atof(temp);
				if (!isback){
					mk_group_sn(cprof,bins,ebins,profile,eprof,area,effexp,nbin,minsn);
				}
				else {
					mk_group_isback_sn(cprof,bins,ebins,profile,eprof,area,effexp,nbin,minsn,backprof,backcounts);
				}				
			}
			else {
				printf("    Invalid entry\n");
				break;
			}
		}while(0);
    }
    else if (!strcmp(temp,"newpar")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			int npar;
			double nval;
			cin >> npar;
			cin >> nval;
			if (npar>0){
				model->SetParameter(npar-1,nval);
			}
			else {
				printf("    Invalid parameter %d\n",npar);
			}
		}while(0);
    }
    else if (!strcmp(temp,"fixpar")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			int npar;
			cin >> npar;
			if (npar>0){
				fix[npar-1]=true;
				double llll=model->GetParameter(npar-1);
				model->FixParameter(npar-1,llll);
			}
			else {
				printf("    Invalid parameter %d\n",npar);
			}
		}while(0);
    }
    else if (!strcmp(temp,"thawpar")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			int npar;
			cin >> npar;
			if (npar>0){
				fix[npar-1]=false;
				double llll=model->GetParameter(npar-1);
				model->SetParameter(npar-1,llll);
			}
			else {
				printf("    Invalid parameter %d\n",npar);
			}
		}while(0);
    }
    else if (!strcmp(temp,"limits")){
		do {
			cin >> temp;
			fitlow=atof(temp);
			cin >> temp;
			fithigh=atof(temp);
			islimits=true;
		}while(0);
    }
    else if (!strcmp(temp,"flux")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			char *pars=new char[200];	
			double fluxlow,fluxhigh;
			cin >> temp;
			fluxlow=atof(temp);
			cin >> temp;
			fluxhigh=atof(temp);
			int npp=model->GetNpar();
			//printf("fluxlow, fluxhigh: %g %g\n",fluxlow,fluxhigh);
			TF1 *ftt=mkfint(modname,(char *)"ftemp");
			for (int k=0;k<npp;k++){
				pars=(char *)model->GetParName(k);
				//printf("temp: %s\n",pars);
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
			double flux=ftt->Integral(fluxlow,fluxhigh);
			printf("    Model count rate: %g counts/sec\n",flux);
			delete ftt;
		}while(0);
    }
    else if (!strcmp(temp,"flobs")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			if (!isprofile){
				printf("    Profile not yet extracted\n");
				break;
			}	
			double fluxlow,fluxhigh;
			cin >> temp;
			fluxlow=atof(temp);
			cin >> temp;
			fluxhigh=atof(temp);
			TF1 *fbkg=new TF1("fbkg","[0]*1.0",0.,1e4);
			double bkg,ebkg;
			char *pars=new char[200];	
			int npp=model->GetNpar();
			TF1 *ftt=mkftemp(modname,(char *)"ftemp");
			for (int k=0;k<npp;k++){
				pars=(char *)model->GetParName(k);
				if (!strcmp(pars,"const")){
					bkg=model->GetParameter(k);
					ebkg=model->GetParError(k);
					ftt->SetParameter(k,0.0);
					ftt->SetParName(k,"const");
					model->SetParName(k,"const");
				}
				else {
					double tpar=model->GetParameter(k);
					ftt->SetParameter(k,tpar);
				}
			}
			fbkg->SetParameter(0,bkg);
			double raddet=findrad(ftt,fbkg,2.);
			if (raddet>bins[nbin-1]) raddet=bins[nbin-1];
			printf("    Maximum detection radius (source intensity similar to bkg): %g\n",raddet);
			double flux,eflux;
			cout << "    Subtract the background? (y/n)";
			cin >> temp;
			if (!strcmp(temp,"y")) {
				getflux(profile,eprof,bins,ebins,bkg,ebkg,fluxlow,fluxhigh,flux,eflux);				
			}
			else {
				getflux(profile,eprof,bins,ebins,0.0,0.0,fluxlow,fluxhigh,flux,eflux);
			}
			printf("    Integrated count rate from %g to %g: %g +/- %g counts/sec\n",fluxlow,fluxhigh,flux,eflux);
			delete ftt;
		}while(0);
    }
    else if (!strcmp(temp,"flcomp")){
		do {
			if (!ismod){
				printf("    No model yet loaded\n");
				break;
			}
			double fluxlow,fluxhigh;
			cin >> temp;
			fluxlow=atof(temp);
			cin >> temp;
			fluxhigh=atof(temp);
			//printf("fluxlow, fluxhigh: %g %g\n",fluxlow,fluxhigh);
			TF1 *ftt=NULL;
			if (!strcmp(modname,"doublebeta")){
				ftt=new TF1("ftt",betaint,0.01,1e4,4);
				//Component 1:
				double tpar=model->GetParameter(0);
				ftt->SetParameter(0,tpar); //beta
				tpar=model->GetParameter(1);
				ftt->SetParameter(1,tpar); //rc1
				tpar=model->GetParameter(4);
				ftt->SetParameter(2,tpar); //norm
				ftt->SetParameter(3,0.0);//const
				double flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 1: %g counts/sec\n",flux);
				//Component 2:
				tpar=model->GetParameter(2);
				ftt->SetParameter(1,tpar); //rc2
				double rat=model->GetParameter(3);
                tpar=model->GetParameter(4);
				ftt->SetParameter(2,tpar*rat); //norm*ratio
				flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 2: %g counts/sec\n",flux);
			}
			else if (!strcmp(modname,"gausbeta")){
				//Component 1:
				ftt=new TF1("ftt",betaint,0.01,1e4,4);
				double tpar=model->GetParameter(0);
				ftt->SetParameter(0,tpar); //beta
				tpar=model->GetParameter(1);
				ftt->SetParameter(1,tpar); //rc1
				tpar=model->GetParameter(2);
				ftt->SetParameter(2,tpar); //amp1
				ftt->SetParameter(3,0.0);//const
				double flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 1 (beta) : %g counts/sec\n",flux);
				//Component 2:
				ftt=NULL;
				ftt=new TF1("ftt",gausint,0.01,1e4,2);
				tpar=model->GetParameter(5);
				ftt->SetParameter(0,tpar); //ampg
				tpar=model->GetParameter(4);
				ftt->SetParameter(1,tpar); //sigma	  
				flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 2 (gauss) : %g counts/sec\n",flux);
			}
			else if (!strcmp(modname,"gausdbeta")){
				ftt=new TF1("ftt",betaint,0.01,1e4,4);
				//Component 1:
				double tpar=model->GetParameter(0);
				ftt->SetParameter(0,tpar); //beta
				tpar=model->GetParameter(1);
				ftt->SetParameter(1,tpar); //rc1
                tpar=model->GetParameter(4);
                ftt->SetParameter(2,tpar); //norm
				ftt->SetParameter(3,0.0);//const
				double flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 1 (beta 1): %g counts/sec\n",flux);
				//Component 2:
				tpar=model->GetParameter(2);
				ftt->SetParameter(1,tpar); //rc2
                double rat=model->GetParameter(3);
                tpar=model->GetParameter(4);
                ftt->SetParameter(2,tpar*rat); //norm*ratio
				flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 2 (beta 2): %g counts/sec\n",flux);
				//Component 3:
				ftt=NULL;
				ftt=new TF1("ftt",gausint,0.01,1e4,2);
				tpar=model->GetParameter(7);
				ftt->SetParameter(0,tpar); //ampg
				tpar=model->GetParameter(6);
				ftt->SetParameter(1,tpar); //sigma	  
				flux=ftt->Integral(fluxlow,fluxhigh);
				printf("    Model count rate for component 3 (gauss) : %g counts/sec\n",flux);
			}
			else {
				printf("    This feature only works for the models with multiple components\n");
				break;
			}
			delete ftt;
		}while(0);
    }
    else if (!strcmp(temp,"ls")){
		system("ls");
    }
    else if (!strcmp(temp,"psf")){
		do{
			if (!isprofile){
				printf("    Profile not yet extracted\n");
				break;
			}
			printf("    Warning: This may take a while\n");
			cout << "    PSF profile (gaus/king) > ";
			cin >> temp;
			if (!strcmp(temp,"gaus")) {
				cout << "    Number of ray-tracing photons per bin > ";
				cin >> temp;
				double ttp;
				ttp=atof(temp);
				int nphot=(int)floor(ttp);
				if (nphot<=0){
					printf("    Invalid entry\n");
					break;
				}
				cout << "    Sigma of the PSF (arcsec) > ";
				cin >> temp;
				double sigpsf=atof(temp);
				if (sigpsf<=0.0){
					printf("    Invalid PSF value\n");
					break;
				}
				psfmat=NULL;
				psfmat=new double[nbin*nbin];
				psfgaussnew(nphot,bins,ebins,exposure,nbin,axes,centroid_ra,centroid_dec,sigpsf,pixsize,psfmat);
				ispsf=true;				
			}
			else if (!strcmp(temp,"king")){
				cout << "    Number of ray-tracing photons per bin > ";
				cin >> temp;
				double ttp;
				ttp=atof(temp);
				int nphot=(int)floor(ttp);
				if (nphot<=0){
					printf("    Invalid entry\n");
					break;
				}
				cout << "    R0 (arcsec) > ";
				cin >> temp;
				double r0=atof(temp);
				if (r0<=0.0){
					printf("    Invalid PSF value\n");
					break;
				}
				cout << "    Slope > ";
				cin >> temp;
				double alpha=atof(temp);
				if (alpha<=0.0){
					printf("    Invalid PSF value\n");
					break;
				}
				psfmat=NULL;
				psfmat=new double[nbin*nbin];
				psfkingnew(nphot,bins,ebins,exposure,nbin,axes,centroid_ra,centroid_dec,r0,alpha,pixsize,psfmat);
				ispsf=true;				
			}
			else {
				printf("    Unknown profile %s\n",temp);
				break;
			}
		}while (0);
    }
    else if (!strcmp(temp,"save")){
		do{
			if (!isprofile){
				printf("    Profile not yet extracted\n");
				break;
			}	
			cout << "    Format (txt/root/fits) > ";
			cin >> temp;
			double *esbs=new double[nbin];
			for (int i=0; i<nbin; i++) {
				esbs[i]=sqrt(eprof[i]*eprof[i]+syserr/100.*profile[i]*syserr/100.*profile[i]);
			}
			if (!strcmp(temp,"txt")){
				if (!ismod){
					printf("    No model loaded\n");
					break;
				}
				cout << "    Output file name > ";
				cin >> temp;
				FILE *ff=fopen(temp,"w");
				fprintf(ff,"  Bin  Ebin  Counts  Profile  Error\n");
				for (int i=0;i<nbin;i++){
					fprintf(ff,"%f  %g  %g  %g  %g\n",bins[i],ebins[i],cprof[i],profile[i],esbs[i]);
				}
				fclose(ff);
			}
			else if (!strcmp(temp,"root")){
				cout << "    Output file name > ";
				cin >> temp;
				TFile *ff=new TFile(temp,"recreate");
				TTree *t=new TTree("t","t");
				double tb,tpr,tep,teb,tcpr,tdep,etdep,tdens,etdens,tmed,etmed;
				TH2F *pps=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
				if (ispsf) {
					for (int i=0; i<nbin; i++) {
						for (int j=0; j<nbin; j++) {
							pps->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
						}
					}
				}
				t->Branch("bins",&tb,"bins/D");
				t->Branch("profile",&tpr,"profile/D");
				t->Branch("eprof",&tep,"eprof/D");
				t->Branch("ebins",&teb,"ebins/D");
				t->Branch("cprof",&tcpr,"cprof/D");
				if (isdepr) {
					t->Branch("deprof",&tdep,"deprof/D");
					t->Branch("edeprof",&etdep,"edeprof/D");
				}
                if (isdens) {
                    t->Branch("density",&tdens,"density/D");
                    t->Branch("edensity",&etdens,"edensity/D");
                }
                if (ismed) {
                    t->Branch("mediansb",&tmed,"mediansb/D");
                    t->Branch("emedian",&etmed,"emedian/D");
                }
				for (int i=0;i<nbin;i++){
					tb=bins[i];
					tpr=profile[i];
					tep=esbs[i];
					teb=ebins[i];
					tcpr=cprof[i];
					if (isdepr) {
						tdep=deprof[i];
						etdep=edeprof[i];
					}
                    if (isdens) {
                        tdens=dens[i];
                        etdep=edens[i];
                    }
                    if (ismed) {
                        tmed=medprof[i];
                        etmed=emedprof[i];
                    }
					t->Fill();
				}
				t->Write();
                double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
				TH1F *hh=new TH1F("profile","profile",nbin,binsh);
				for (int i=0;i<nbin;i++){
					hh->SetBinContent(i+1,profile[i]);
					hh->SetBinError(i+1,esbs[i]);
				}
				if (ismod){
					TH1F *hhd=new TH1F("chi","chi",nbin,binsh);
					for (int i=0;i<nbin;i++){
						double bin=hh->GetBinCenter(i+1);
						double dat=hh->GetBinContent(i+1);
						double edat=hh->GetBinError(i+1);
						if (edat>0.0){
							double dc=(dat-model->Eval(bin))/edat;
							hhd->SetBinContent(i+1,dc);
							hhd->SetBinError(i+1,1.0);
						}
						else {
							double dc=dat-model->Eval(bin);
							hhd->SetBinContent(i+1,dc);
							hhd->SetBinError(i+1,1.0);
						}
					}
					hhd->Write();
				}
				if (isdepr) {
					TH1F* hdep=new TH1F("hdep","hdep",nbin,binsh);
					for (int i=0;i<nbin;i++){
						hdep->SetBinContent(i+1,deprof[i]);
						hdep->SetBinError(i+1,edeprof[i]);
					}
					hh->SetXTitle("Radius [arcmin]");
					if (isr500){
						hdep->SetXTitle("r/r_{500}");
					}
					if (isr200){
						hdep->SetXTitle("r/r_{200}");
					}
                    if (iskpc) {
                        hdep->SetXTitle("Radius [kpc]");
                    }
					hdep->SetYTitle("Volume emission density [counts s^{-1} kpc^{-3}]");
					hdep->Write();
				}
                if (isdens) {
                    TH1F* hdens=new TH1F("hdens","hdens",nbin,binsh);
                    for (int i=0;i<nbin;i++){
                        hdens->SetBinContent(i+1,dens[i]);
                        hdens->SetBinError(i+1,edens[i]);
                    }
                    hh->SetXTitle("Radius [arcmin]");
                    if (isr500){
                        hdens->SetXTitle("r/r_{500}");
                    }
                    if (isr200){
                        hdens->SetXTitle("r/r_{200}");
                    }
                    if (iskpc) {
                        hdens->SetXTitle("Radius [kpc]");
                    }
                    hdens->SetYTitle("Electron number density [cm^{-3}]");
                    hdens->Write();
                }
				if (isback){
					TH1F *hback=new TH1F("back","back",nbin,binsh);
					for (int i=0; i<nbin; i++) {
						hback->SetBinContent(i+1,backprof[i]);
						hback->SetBinError(i+1,0.0);
					}
					hback->Write();
				}
				if (isscat) {
					TH1F* hsc=new TH1F("scat","scat",nbin,binsh);
					for (int i=0;i<nbin;i++){
						hsc->SetBinContent(i+1,scat[i]);
						hsc->SetBinError(i+1,escat[i]);
					}
					hsc->SetTitle("");
					hsc->SetStats(false);
					hsc->SetXTitle("Radius [arcmin]");
					if (isr500){
						hsc->SetXTitle("r/r_{500}");
					}
					if (isr200){
						hsc->SetXTitle("r/r_{200}");
					}
                    if (iskpc) {
                        hsc->SetXTitle("Radius [kpc]");
                    }
					hsc->SetYTitle("Scatter");
					hsc->Write();
				}
                if (ismed) {
                    TH1F* hmed=new TH1F("mediansb","mediansb",nbin,binsh);
                    for (int i=0;i<nbin;i++){
                        hmed->SetBinContent(i+1,medprof[i]);
                        hmed->SetBinError(i+1,emedprof[i]);
                    }
                    hmed->SetTitle("");
                    hmed->SetStats(false);
                    hmed->SetXTitle("Radius [arcmin]");
                    if (isr500){
                        hmed->SetXTitle("r/r_{500}");
                    }
                    if (isr200){
                        hmed->SetXTitle("r/r_{200}");
                    }
                    if (iskpc) {
                        hmed->SetXTitle("Radius [kpc]");
                    }
                    hmed->SetYTitle("SB [counts s^{-1} arcmin^{-2}]");
                    hmed->Write();
                }
				hh->Write();
				if (ismod) model->Write();
				if (ispsf) pps->Write();
				ff->Close();
				//delete pps;
				//delete hh;
				//delete hhd;
				delete [] binsh;
			}
			else if (!strcmp(temp,"fits")){
				cout << "    Output file name > ";
				cin >> temp;
				char nfn[200];
				sprintf(nfn,"!%s",temp);
				int status=0;
				fitsfile *x;
				fits_create_file(&x,nfn,&status);
				status=save_fits_basic(x,nbin,bins,ebins,profile,esbs,cprof,area,effexp);
				if (status!=0) break; 
				int tfields=7;
				if (isback) {
					status=save_fits_isback(x,nbin,backprof,backcounts,tfields);
				}
				if (isdepr) {
					status=save_fits_isdepr(x,nbin,deprof,edeprof,tfields);
				}
                if (isdens) {
                    status=save_fits_isdens(x,nbin,dens,edens,tfields);
                }
				if (ismod) {
					status=save_fits_ismod1(x,nbin,bins,model,profile,esbs,tfields);
				}
				if (isscat) {
					status=save_fits_isscat(x,nbin,scat,escat,tfields);
				}
                if (ismed){
                    status=save_fits_ismed(x,nbin,medprof,emedprof,tfields);
                }
				if (status!=0) {
					printf("    Error %d\n",status);
					break;
				}
				status=save_fits_keys(x,axes,exposure,centroid_ra,centroid_dec,wcs_inp,ellipse,ellang,aoverb,sector,angl,angh,imgfile,expfile,backfile,isback);
				if (ismod) {
					status=save_fits_ismod2(x,model,modname,names);
					if (status!=0) {
						printf("    Error %d\n",status);
						break;
					}
				}
				if (ispsf) {
					status=save_fits_ispsf(x,nbin,psfmat);
					if (status!=0) {
						printf("    Error %d\n",status);
						break;
					}
				}
				fits_close_file(x,&status);
				if (status!=0){
					printf("    Error %d\n",status);
					break;
				}
				else {
					printf("    Results saved in file %s\n",temp);
				}
			}
			else {
				printf("    Unknown format %s\n",temp);
			}
			delete [] esbs;
		}while(0);
    }
      else if (!strcmp(temp,"loadfits")){
          cout << "    Input FITS file > ";
          cin >> temp;
          int status=load_fits_structure(temp,nbin,axes,isback,ispsf,ismod,modname);
          if (status==0){
              img=NULL;
              img=new double[axes[0]*axes[1]];
              exposure=NULL;
              exposure=new double[axes[0]*axes[1]];
              bins=NULL;
              ebins=NULL;
              profile=NULL;
              eprof=NULL;
              cprof=NULL;
              area=NULL;
              effexp=NULL;
              bins=new double[nbin];
              ebins=new double[nbin];
              profile=new double[nbin];
              eprof=new double[nbin];
              cprof=new double[nbin];
              area=new double[nbin];
              effexp=new double[nbin];
              if (isback){
                  backmap=NULL;
                  backmap=new double[axes[0]*axes[1]];
                  backprof=NULL;
                  backprof=new double[nbin];
                  backprof=new double[nbin];
                  backcounts=new double[nbin];
              }
              if (ismod){
                  model=mkftemp(modname,(char *)"model");
                  model->SetLineWidth(2);
                  model->SetLineColor(kBlue);
              }
              if (ispsf){
                  psfmat=NULL;
                  psfmat=new double[nbin*nbin];
              }
              status=load_fits(temp,nbin, axes,bins,ebins,profile,eprof,cprof,area,effexp,isprofile,imgfile,expfile,img,isimg,exposure,isexp,centroid_ra,centroid_dec,pixsize,cdelt1,crval1,crval2,crpix1,crpix2,wcs_inp,backfile,backmap,isback,backprof,backcounts,names,model,psfmat);
              maxrad=bins[nbin-1]+ebins[nbin-1];
              if (status==0){
                  printf("    Session recovered!\n");
              }
          }
      }
    else if (!strcmp(temp,"statistics")){
		cout << "    Statistics (chi2/cash) (default=chi2) > ";
		cin >> statmet;
		if ((strcmp(statmet,"chi2"))&&(strcmp(statmet,"cash"))){
			printf("    Unknown statistics %s\n",statmet);
			printf("    Resetting to chi2\n");
			sprintf(statmet,"chi2");
		}
    }
    else if (!strcmp(temp,"r500")){
		// Uses scaling relations from Arnaud et al. 2005
		double z,kt;
		cout << "    Redshift > ";
		cin >> temp;
		z=atof(temp);
		cout << "    Temperature > ";
		cin >> temp;
		kt=atof(temp);
		double dist=angdist(z)*1e3; //angular diameter distance in kpc
		double hz=sqrt(omegam*pow((1.+z),3.)+omegal);
		double beta=0.0;
		double b500=0.0;
		if (kt>3.5){
			beta=0.5;
			b500=1129.;
		}
		else {
			beta=0.57;
			b500=1104.;
		}
		double r500=b500*pow(kt/5.,beta)/hz;
		r500obs=(r500/dist)*180./TMath::Pi()*60.;
		double frac=maxrad/r500obs;
		if (frac>1.0)frac=1.0;
		printf("    Angular diameter distance to the object: %g Mpc\n",dist/1e3);
		printf("    R500 = %g kpc (%g arcmin)\n",r500,r500obs);
		printf("    Fraction of R500 considered: %g\n",frac);
		cout << "    Plot the results rescaled by r500? (y/n)";
		cin >> temp;
		if (!strcmp(temp,"y")){
			isr500=true;
            isr200=false;
            iskpc=false;
		}
		else {
			isr500=false;
		}
    }
    else if (!strcmp(temp,"r200")){
		// Uses scaling relations from Arnaud et al. 2005
		double z,kt;
		cout << "    Redshift > ";
		cin >> temp;
		z=atof(temp);
		cout << "    Temperature > ";
		cin >> temp;
		kt=atof(temp);
		double dist=angdist(z)*1e3; //angular diameter distance in kpc
		double hz=sqrt(omegam*pow((1.+z),3.)+omegal);
		double beta=0.0;
		double b200=0.0;
		if (kt>3.5){
			beta=0.5;
			b200=1714.;
		}
		else {
			beta=0.57;
			b200=1674.;
		}
		double r200=b200*pow(kt/5.,beta)/hz;
		r200obs=(r200/dist)*180./TMath::Pi()*60.;
		double frac=maxrad/r200obs;
		if (frac>1.0)frac=1.0;
		printf("    Angular diameter distance to the object: %g Mpc\n",dist/1e3);
		printf("    R200 = %g kpc (%g arcmin)\n",r200,r200obs);
		printf("    Fraction of R200 considered: %g\n",frac);
		cout << "    Plot the results rescaled by r200? (y/n)";
		cin >> temp;
		if (!strcmp(temp,"y")){
			isr200=true;
            isr500=false;
            iskpc=false;
		}
		else {
			isr200=false;
		}
    }
    else if (!strcmp(temp,"kpc")){
		double z;
		cout << "    Redshift > ";
		cin >> temp;
		z=atof(temp);
		double dist=angdist(z)*1e3; //angular diameter distance in kpc
		kpcobs=(1./dist)*180./TMath::Pi()*60.;//1 kpc in arcmin
		printf("    Angular diameter distance to the object: %g Mpc\n",dist/1e3);
		printf("    1 arcmin = %g kpc\n",1./kpcobs);
		cout << "    Plot the results in kpc unit? (y/n)";
		cin >> temp;
		if (!strcmp(temp,"y")){
			iskpc=true;
            isr500=false;
            isr200=false;
		}
		else {
			iskpc=false;
		}
    }
    else if (!strcmp(temp,"ellipticity")){
		do {
			// Uses definition of Hashimoto et al. 2007
			if (!isprofile) {
				printf("    Profile not yet extracted\n");
				break;
			}
			cout << "       Maximum radius (arcmin) > ";
			cin >> temp;
			double rmax=atof(temp)/pixsize/60.;//pixel
			if (rmax<=0.0 || rmax*pixsize*60.>maxrad) {
				printf("    Invalid values\n");
				break;
			}
			int nang=36; //5 deg increment
			double A=0.0;
			double B=0.0;
			double ellang=0.0;
			for (int i=0; i<nang; i++) {
				double angle=TMath::Pi()*i/nang;
				double cx=meanx(img,exposure,centroid_ra,centroid_dec,rmax,axes,angle);
				double cy=meany(img,exposure,centroid_ra,centroid_dec,rmax,axes,angle);
				double mx2=meanx2(img,exposure,cx,centroid_ra,centroid_dec,rmax,axes,angle);
				double my2=meany2(img,exposure,cy,centroid_ra,centroid_dec,rmax,axes,angle);
				double mxy=meanxy(img,exposure,cx,cy,centroid_ra,centroid_dec,rmax,axes,angle);
				double ta=sqrt((mx2+my2)/2.0+sqrt((mx2-my2)/2.0*(mx2-my2)/2.0+mxy*mxy));
				double tb=sqrt((mx2+my2)/2.0-sqrt((mx2-my2)/2.0*(mx2-my2)/2.0+mxy*mxy));
				if (ta>A) {
					A=ta;
					B=tb;
					ellang=angle;
				}
			}
			printf("    Major axis: %g arcmin\n",A*pixsize*60.);
			printf("    Minor axis: %g arcmin\n",B*pixsize*60.);
			printf("    Angle of largest ellipticity: %g deg\n",ellang*180./TMath::Pi());
			printf("    Ellipticity: %g\n",1.-B/A);
		}
		while (0);
    }
    else if (!strcmp(temp,"centroidshift")){
		do {
			// Uses definition of Rasia et al. 2012
			if (!isprofile) {
				printf("    Profile not yet extracted\n");
				break;
			}
			cout << "       Maximum radius (arcmin) > ";
			cin >> temp;
			double rmax=atof(temp);//arcmin
			if (rmax<=0.0) {
				printf("    Invalid value\n");
				break;
			}
            int niter=(int)floor(rmax/(5.*pixsize*60.));//5-pixel separation between annuli
            double w=centroid_shift(img,exposure,centroid_ra,centroid_dec,axes,pixsize,rmax,niter);
            printf("    Centroid shift: %g\n",w);
		}
		while (0);
    }
    else if (!strcmp(temp,"csb")){
        do {
            // Uses definition of Santos et al. 2008
            if (!ismod) {
                printf("    No model loaded\n");
                break;
            }
            cout << "       Cluster redshift > ";
            cin >> temp;
            double z=atof(temp);
            if (z<=0.0) {
                printf("    Invalid value\n");
                break;
            }
            double tcsb=csb(z,modname,model);
            printf("    Surface-brightness concentration: %g\n",tcsb);
        }
        while (0);
    }
    else if (!strcmp(temp,"backsub")){
		do {
			if (!isprofile) {
				printf("    Profile not yet extracted\n");
				break;
			}
			if (!ismod) {
				printf("    No model loaded\n");
				break;
			}
			double back,eback;
			int npar=model->GetNpar();
			char *pname=new char[20];
			for (int i=0; i<npar; i++) {
				pname=(char *)model->GetParName(i);
				if (!strcmp(pname,"const")) {
					back=model->GetParameter(i);
					eback=model->GetParError(i);
				}
			}
			for (int i=0; i<nbin; i++) {
				double pr=profile[i];
				double epr=eprof[i];
				double sysp=syserr*pr/100.;
				if (!strcmp(modname,"backfit")) {
					back=model->GetParameter(0)+model->GetParameter(1)*backprof[i];
					eback=sqrt(model->GetParError(0)*model->GetParError(0)+model->GetParError(1)*model->GetParError(1)*backprof[i]*backprof[i]);
				}					
				profile[i]=pr-back;
				eprof[i]=sqrt(epr*epr+eback*eback+sysp*sysp);
			}
		}
		while (0);
    }
    else if (!strcmp(temp,"angle2dist")){
      double z,angle;
      cout << "    Redshift > ";
      cin >> temp;
      z=atof(temp);
      cout << "    Angle (arcmin) > ";
      cin >> temp;
      angle=atof(temp);
      double dist=angdist(z)*1e3; //angular diameter distance in kpc
      printf("    Angular diameter distance to the object: %g Mpc\n",dist/1e3);
      double angrad=angle/60.*TMath::Pi()/180.; //angle in radians
      double radius=dist*angrad; //physical radius in kpc
      printf("    Physical distance corresponding to the angle: %g kpc\n",radius);
    }
    else if (!strcmp(temp,"deproject")){
		if (!isprofile) {
			printf("    Profile not yet extracted\n");
			break;
		}
		double z;
		cout << "    Redshift > ";
		cin >> temp;
		z=atof(temp);
		cout << "    Error estimation method (1: error propagation; 2: Monte Carlo) > ";
		cin >> temp;
		int method=atoi(temp);
		if (method<1 && method>2) {
			printf("    Unknow method %s\n",temp);
			break;
		}
		deprof=NULL;
		deprof=new double[nbin];
		edeprof=NULL;
		edeprof=new double[nbin];
		if (method==1) {
			deproject(z,profile,eprof,nbin,bins,ebins,deprof,edeprof);			
		}
		else {
			mcdeproj(z,profile,eprof,nbin,bins,ebins,deprof,edeprof);
		}

		c1->Clear();
		c1->SetTickx();
		c1->SetTicky();
		if (!logy){
			c1->SetLeftMargin(0.15);
		}
		else {
			c1->SetLeftMargin(0.12);
		}
		c1->SetBottomMargin(0.12);
        double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
		TH1F* hh=new TH1F("hh","hh",nbin,binsh);
		for (int i=0;i<nbin;i++){
			hh->SetBinContent(i+1,deprof[i]);
			hh->SetBinError(i+1,edeprof[i]);
		}
		hh->SetTitle("");
		hh->SetStats(false);
		hh->SetXTitle("Radius [arcmin]");
		if (isr500){
			hh->SetXTitle("r/r_{500}");
		}
		if (isr200){
			hh->SetXTitle("r/r_{200}");
		}
        if (iskpc) {
            hh->SetXTitle("Radius [kpc]");
        }
		if (!logy){
			hh->GetYaxis()->SetTitleOffset(1.8);
		}
		else {
			hh->GetYaxis()->SetTitleOffset(1.5);
		}
		hh->GetXaxis()->SetTitleOffset(1.2);
		hh->SetYTitle("Volume emission density [counts s^{-1} kpc^{-3}]");
		hh->Draw();
		c1->Update();
        if (!scripting) {
            theApp.Run(kTRUE);
        }
		
		isdepr=true;
        delete [] binsh;
		delete hh;
    }
    else if (!strcmp(temp,"density")){
        if (!isdepr) {
            printf("    Deprojection not performed yet\n");
            break;
        }
        double z,conv;
        cout << "    Redshift > ";
        cin >> temp;
        z=atof(temp);
        double ad=angdist(z)*Mpc; //cm
        double kpc3=TMath::Power(Mpc/1e3,3.0);
        double numfactmekal=1e-14/(4.*TMath::Pi()*ad*ad*(1.+z)*(1.+z)); //cm-2
        double nhc=1.21; // ne = 1.21*nH
        cout << "    CR to EM conversion > ";
        cin >> temp;
        conv=atof(temp);
        double cf=1./conv;
        dens=NULL;
        edens=NULL;
        dens=new double[nbin];
        edens=new double[nbin];
        for (int i=0; i<nbin; i++) {
            double pp=deprof[i]/kpc3*cf/numfactmekal/nhc; //nh^2, [cm-6]
            double epp=edeprof[i]/kpc3*cf/numfactmekal/nhc; //nh^2, [cm-6]
            if (pp>0.0) {
                dens[i]=sqrt(pp);
                edens[i]=epp/2/sqrt(pp);
            }
            else {
                dens[i]=-sqrt(-pp);
                edens[i]=epp/2/sqrt(-pp);
            }
        }
        double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
        TH1F* hh=new TH1F("hh","hh",nbin,binsh);
        for (int i=0;i<nbin;i++){
            hh->SetBinContent(i+1,dens[i]);
            hh->SetBinError(i+1,edens[i]);
        }
        hh->SetTitle("");
        hh->SetStats(false);
        hh->SetXTitle("Radius [arcmin]");
        if (isr500){
            hh->SetXTitle("r/r_{500}");
        }
        if (isr200){
            hh->SetXTitle("r/r_{200}");
        }
        if (iskpc) {
            hh->SetXTitle("Radius [kpc]");
        }
        if (!logy){
            hh->GetYaxis()->SetTitleOffset(1.8);
        }
        else {
            hh->GetYaxis()->SetTitleOffset(1.5);
        }
        hh->GetXaxis()->SetTitleOffset(1.2);
        hh->SetYTitle("Proton number density [cm^{-3}]");
        hh->Draw();
        c1->Update();
        if (!scripting) {
            theApp.Run(kTRUE);
        }
        
        isdens=true;
        delete [] binsh;
        delete hh;
    }
    else if (!strcmp(temp,"quit")){
		printf("    Good bye\n");
		cont=false;
    }
    else if (!strcmp(temp,"savemodimg")){
		do {
			if (!ismod){
				printf("    No model loaded\n");
				break;
			} 	
			sprintf(comm,"    Apply vignetting correction? (default=y) > ");
			cout << comm;
			cin >> temp;
            bool isvig=true;
            if (!strcmp(temp,"n")||!strcmp(temp,"N")) {
                isvig=false;
            }
			int npp=model->GetNpar();
			TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			for (int k=0;k<npp;k++){
				double tpar=model->GetParameter(k);
				ftemp->SetParameter(k,tpar);
			}
			double *modimg=new double[axes[0]*axes[1]];
			if (!isvig) {
                sprintf(comm,"    Include pixles with 0 exposure? (default=y) > ");
                cout << comm;
                cin >> temp;
                bool iszero=true;
                if (!strcmp(temp,"n")||!strcmp(temp,"N")) {
                    iszero=false;
                }
                double me=TMath::MaxElement(axes[0]*axes[1],exposure);
                double *fakeexp=new double[axes[0]*axes[1]];
                for (int i=0; i<axes[0]*axes[1]; i++) {
                    if ((exposure[i]<0.05*me) && (!iszero)) {
                        fakeexp[i]=0.0;
                    }
                    else {
                        fakeexp[i]=me;
                    }
                }
                mk_mod_img(ftemp,fakeexp,modimg,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse);
                delete [] fakeexp;
            }
            else {
                mk_mod_img(ftemp,exposure,modimg,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse);
            }
			sprintf(comm,"    File name > ");
			cout << comm;
			cin >> temp;
			status=save_img(temp,modimg,axes,crpix1,crval1,crpix2,crval2,cdelt1,pixsize);
			delete [] modimg;
			delete ftemp;
		}
		while (0);
    }
	else if (!strcmp(temp,"syst")) {
		cout << "    Systematic error (%) > ";
		cin >> temp;
		syserr=atof(temp);
	}
      else if (!strcmp(temp,"voronoi")) {
          do {
              if (!isimg){
                  printf("    Image file not yet loaded\n");
                  break;
              }
              if (!isexp){
                  printf("    Exposure map not yet loaded\n");
                  break;
              }
              cout << "    Target number of counts per bin > ";
              cin >> temp;
              int mincounts=atoi(temp);
              if (mincounts==0){
                  printf("    Invalid number of counts %s\n",temp);
                  break;
              }
              char *outvoronoi=new char[200];
              cout << "    Output image > ";
              cin >> outvoronoi;
              voronoimap=NULL;
              voronoierr=NULL;
              voronoimap=new double[axes[0]*axes[1]];
              voronoierr=new double[axes[0]*axes[1]];
              voronoi(img,axes,mincounts,pixsize,crpix1,crpix2,crval1,crval2,cdelt1,exposure,isback,backmap,outvoronoi,voronoimap,voronoierr);
              isvoronoi=true;
          }
          while (0);
      }
    else if (!strcmp(temp,"savedeviations")){
		do {
			if (!isimg){
				printf("    Image file not yet loaded\n");
				break;
			}
			if (!ismod){
				printf("    No model loaded\n");
				break;
			} 	
			int npp=model->GetNpar();
			TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			for (int k=0;k<npp;k++){
				double tpar=model->GetParameter(k);
				ftemp->SetParameter(k,tpar);
			}
			double *modimg=new double[axes[0]*axes[1]];
			double *devmod=new double[axes[0]*axes[1]];
			mk_mod_img(ftemp,exposure,modimg,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse);
			for (int i=0;i<axes[0]*axes[1];i++){
                if (isback) {
                    if (img[i]>0){
                        devmod[i]=(img[i]-modimg[i]-backmap[i])/sqrt(img[i]);
                    }
                    else {
                        devmod[i]=-modimg[i]-backmap[i];
                    }
                }
                else {
                    if (img[i]>0){
                        devmod[i]=(img[i]-modimg[i])/sqrt(img[i]);
                    }
                    else {
                        devmod[i]=-modimg[i];
                    }
                }
				if (exposure[i]==0.0){
					devmod[i]=0.0;
				}
			}
			sprintf(comm,"    File name > ");
			cout << comm;
			cin >> temp;
			status=save_img(temp,devmod,axes,crpix1,crval1,crpix2,crval2,cdelt1,pixsize);
			delete [] modimg;
			delete [] devmod;
			delete ftemp;
		}
		while (0);
    }
    else if (!strcmp(temp,"savescat")){
		do {
			if (!isimg){
				printf("    Image file not yet loaded\n");
				break;
			}
            if (!isprofile){
                printf("    Profile not yet extracted\n");
                break;
            }
            if (!ismod) {
                printf("    No model loaded\n");
                break;
            }
            double *devmod=new double[axes[0]*axes[1]];
            double back,eback;
            int npar=model->GetNpar();
            char *pname=new char[20];
            for (int i=0; i<npar; i++) {
                pname=(char *)model->GetParName(i);
                if (!strcmp(pname,"const")) {
                    back=model->GetParameter(i);
                    eback=model->GetParError(i);
                }
            }
            cout << "    Number of sectors > ";
            cin >> temp;
            int nsect=atoi(temp);
            if (nsect<2) {
                printf("    Not enough sectors\n");
                break;
            }
            double angles[nsect+1];
            angles[0]=0.0;
            for (int i=1; i<nsect+1; i++) {
                double ai=2.*TMath::Pi()/nsect*i;
                angles[i]=ai;
            }			 
            cout << "    Subtract the background? (y/n) > ";
            cin >> temp;
            double **allprofs=new double*[nsect];
            double **allerrs=new double*[nsect];
            for (int i=0; i<nsect; i++) {
                allprofs[i]=new double[nbin];
                allerrs[i]=new double [nbin];
            }
            for (int i=0; i<nsect; i++) {
                double angl=angles[i];
                double angh=angles[i+1];
                mk_sector_scat(img,exposure,angl,angh,bins,ebins,nbin,allprofs[i],allerrs[i],axes,centroid_ra,centroid_dec,pixsize,maxrad);
                if (!strcmp(temp,"y")) {
                    if (isback) {
                        double *backsect=new double[nbin];
                        double *dummyerr=new double[nbin];
                        mk_sector_scat(backmap,exposure,angl,angh,bins,ebins,nbin,backsect,dummyerr,axes,centroid_ra,centroid_dec,pixsize,maxrad);
                        for (int nb=0; nb<nbin; nb++) {
                            allprofs[i][nb]-=backsect[nb];
                        }
                        delete [] backsect;
                        delete [] dummyerr;
                    }
                    backsub(back,eback,nbin,allprofs[i],allerrs[i]);
                }
            }
			sprintf(comm,"    File name > ");
			cout << comm;
			cin >> temp;
            mk_scat_img(bins,nbin,nsect,profile,eprof,allprofs,allerrs,devmod,axes,centroid_ra,centroid_dec,pixsize,ellang,aoverb,ellipse,maxrad);
            status=save_img(temp,devmod,axes,crpix1,crval1,crpix2,crval2,cdelt1,pixsize);
            for (int i=0; i<nsect; i++) {
                delete [] allprofs[i];
                delete [] allerrs[i];
            }
            delete [] allprofs;
            delete [] allerrs;
        }
        while (0);
    }
     else if (!strcmp(temp,"multinest")){
		 do {
			 if (!isimg){
				 printf("    Image file not yet loaded\n");
				 break;
			 }
			 if (!isexp){
				 printf("    Exposure file not yet loaded\n");
				 break;
			 }
			 if (!isprofile){
				 printf("    Profile not yet extracted\n");
				 break;
			 }
			 if (!ismod){
				 printf("    No model loaded\n");
				 break;
			 }
			 char *parname=new char[200];
			 int npfunc=model->GetNpar();
			 int npfit=0;
			 for (int i=0;i<npfunc;i++){
				if (!fix[i]){
					npfit++;
				}
			 }
             printf("npfit: %d\n",npfit);
			 cout << "    Name of MultiNest output files > ";
			 cin >> nestname;
			 double *boundlow=new double[npfit];
			 double *boundhigh=new double[npfit];
			 islogpar=new bool[npfit];
			 TF1 *ftemp=mkftemp(modname,(char *)"ftemp");
			 int tp=0;
			 sprintf(temp,"%s_params.txt",nestname);
			 FILE *saveparams=fopen(temp,"w");
			 for (int k=0;k<npfunc;k++){
				 parname=(char *)model->GetParName(k);
				 ftemp->SetParName(k,parname);
				 double tpar=model->GetParameter(k);
				 if (fix[k]){
					 ftemp->FixParameter(k,tpar);
					 printf("    Parameter %d is frozen\n",k+1);
				 }
				 else {
					 ftemp->SetParameter(k,tpar);
					 islogpar[tp]=false;
					 cout << "    Set log prior on parameter " << k+1 << "? (y/n; default=n) > ";
					 cin >> temp;
					 if (!strcmp(temp,"y")||!strcmp(temp,"Y")||!strcmp(temp,"yes")){
						islogpar[tp]=true;
					 }
					 cout << "    Set lower and upper boundaries for parameter " << k+1 << " > ";
					 cin >> temp;
					 boundlow[tp]=atof(temp);
					 cin >> temp;
					 boundhigh[tp]=atof(temp);
					 if (islogpar[tp]){
					 	fprintf(saveparams,"%d %g log %g %g\n",k+1,tpar,boundlow[tp],boundhigh[tp]);
						boundlow[tp]=log10(boundlow[tp]);
						boundhigh[tp]=log10(boundhigh[tp]);
					 }
					 else {
					 	fprintf(saveparams,"%d %g flat %g %g\n",k+1,tpar,boundlow[tp],boundhigh[tp]);
					 }
					 tp++;
				 }
			 }
			 fclose(saveparams);
			 double *binsh=mkbinsh(nbin,bins,ebins,isr200,r200obs,isr500,r500obs,iskpc,kpcobs);
			 TH1F *hc=new TH1F("hc","hc",nbin,binsh);
			 for (int i=0;i<nbin;i++){
				 hc->SetBinContent(i+1,cprof[i]);
				 hc->SetBinError(i+1,effexp[i]*area[i]);
			 }
			 histpsfmat=NULL;
			 if (ispsf) {
				 histpsfmat=new TH2F("psf","psf",nbin,0.0,nbin*1.0,nbin,0.0,nbin*1.0);
				 for (int i=0; i<nbin; i++) {
					 for (int j=0; j<nbin; j++) {
						 histpsfmat->SetBinContent(i+1,j+1,psfmat[i*nbin+j]);
					 }
				 }
			 }
			 if (!islimits){
				 fitlow=0.0;
				 fithigh=maxrad;
			 }
			 passfitlow=fitlow;
			 passfithigh=fithigh;
			 passfmod=ftemp;
			 passhist=hc;
			 passbck=backcounts;
			sprintf(temp,"%s_params.txt",nestname);
			FILE *fparams=fopen(temp,"r");
			for (int k=0;k<npfit;k++){
				int tp;
				double vp,blp,bhp;
				fscanf(fparams,"%d %lf %s %lf %lf\n",&tp,&vp,temp,&blp,&bhp);
				if (!strcmp(temp,"log")){
					islogpar[tp-1]=true;
				}
				else {
					islogpar[tp-1]=false;
				}
			}
			 run_multinest(ftemp,boundlow,boundhigh,fix,islogpar,statmet,nestname);
			 chains=NULL;
			 chains=new double*[npfunc];
			 load_chains(nestname,npfunc,fix,npchain,chains);
			 if (ispsf) {
				 histpsfmat->Delete();
			 }
			 ftemp->Delete();
			 hc->Delete();
			 delete [] boundlow;
			 delete [] boundhigh;
			isnested=true;
		 }while(0);
     }
	else if (!strcmp(temp,"loadchains")){
		do{
		if (!ismod){
			printf("    No model loaded\n");
			break;
		}
		cout << "    Name of MultiNest output files > ";
		cin >> nestname;
		sprintf(temp,"%s_params.txt",nestname);
		int npfit=line_num(temp);
		FILE *fparams=fopen(temp,"r");
		int npfunc=model->GetNpar();
		islogpar=new bool[npfit];
		for (int k=0;k<npfit;k++){
			int tp;
			double vp,blp,bhp;
			fscanf(fparams,"%d %lf %s %lf %lf\n",&tp,&vp,temp,&blp,&bhp);
			if (!strcmp(temp,"log")){
				islogpar[tp-1]=true;
			}
			else {
				islogpar[tp-1]=false;
			}
		}
		chains=NULL;
		chains=new double*[npfunc];
		load_chains(nestname,npfunc,fix,npchain,chains);
		if (npchain>0){
			isnested=true;			
		}
		}while(0);
	}
	else if (!strcmp(temp,"margin")){
		do {
			if (!isnested){
				printf("    No active chains\n");
				break;
			}
			 if (!ismod){
				 printf("    No model loaded\n");
				 break;
			 }
			int nbpost=100;
			int nb2d=20;
			int npfit=0;
			int npfunc=model->GetNpar();
			for (int i=0;i<npfunc;i++){
				if (!fix[i]) npfit++;
			}
			c1->Clear();
			c1->SetTickx();
			c1->SetTicky();
			c1->Divide(npfit,npfit);
			TF1 *fdiv=new TF1("fdiv","[0]",0.,1e3);
			fdiv->SetParameter(0,npchain);
			int npad=1;
			char *parnam=new char[200];
			for (int i=0;i<npfunc;i++){
				for (int j=0;j<npfunc;j++){
					if (!fix[i]&&!fix[j]){
						c1->cd(npad);
						gPad->SetTickx();
						gPad->SetTicky();
						gPad->SetLeftMargin(0.15);
						gPad->SetBottomMargin(0.15);
						gPad->SetTopMargin(0.05);
						gPad->SetRightMargin(0.05);
						if (i==j){
							double minval=TMath::MinElement(npchain,chains[i]);
							double maxval=TMath::MaxElement(npchain,chains[i]);
							TH1F *hpost;
							sprintf(temp,"post%d",i);
							if (islogpar[i]){
								hpost=new TH1F(temp,"",nbpost,pow(10.,minval),pow(10.,maxval));
							}
							else {
								hpost=new TH1F(temp,"",nbpost,minval,maxval);
							}
							for (int np=0;np<npchain;np++){
								if (islogpar[i]){
									hpost->Fill(pow(10.,chains[i][np]));
									gPad->SetLogx();
								}
								else {
									hpost->Fill(chains[i][np]);
									gPad->SetLogx(0);
								}
							}
							hpost->Divide(fdiv);
							parnam=(char *)model->GetParName(i);
							hpost->SetXTitle(parnam);
							//hpost->SetYTitle("Probability");
							hpost->GetXaxis()->SetTitleOffset(1.2);
							hpost->GetXaxis()->SetLabelSize(0.06);
		      			  	hpost->GetXaxis()->CenterTitle();
							hpost->GetXaxis()->SetTitleSize(0.06);
							hpost->GetYaxis()->SetTitleOffset(2.0);
							hpost->GetYaxis()->SetLabelSize(0.06);
		     			   	hpost->GetYaxis()->CenterTitle();
							hpost->GetYaxis()->SetTitleSize(0.06);
							hpost->SetStats(0);
							hpost->SetFillColor(kCyan);
							hpost->Draw();
						}
						else if (i<j){
							
							double minvx=TMath::MinElement(npchain,chains[i]);
							double maxvx=TMath::MaxElement(npchain,chains[i]);
							double minvy=TMath::MinElement(npchain,chains[j]);
							double maxvy=TMath::MaxElement(npchain,chains[j]);
							TH2F *hxy;
							sprintf(temp,"hxy%d%d",i,j);
							if (islogpar[i]&&islogpar[j]){
								hxy=new TH2F(temp,"",nb2d,pow(10.,minvx),pow(10.,maxvx),nb2d,pow(10.,minvy),pow(10.,maxvy));
							}
							else if (islogpar[i]&&!islogpar[j]){
								hxy=new TH2F(temp,"",nb2d,pow(10.,minvx),pow(10.,maxvx),nb2d,minvy,maxvy);
							}
							else if (!islogpar[i]&&islogpar[j]){
								hxy=new TH2F(temp,"",nb2d,minvx,maxvx,nb2d,pow(10.,minvy),pow(10.,maxvy));
							}
							else {
								hxy=new TH2F(temp,"",nb2d,minvx,maxvx,nb2d,minvy,maxvy);
							}
							for (int np=0;np<npchain;np++){
								double vx,vy;
								if (islogpar[i]){
									vx=pow(10.,chains[i][np]);
									gPad->SetLogx();
								}
								else {
									vx=chains[i][np];
									gPad->SetLogx(0);
								}
								if (islogpar[j]){
									vy=pow(10.,chains[j][np]);
									gPad->SetLogy();
								}
								else {
									vy=chains[j][np];
									gPad->SetLogy(0);
								}
								hxy->Fill(vx,vy);
							}
							
							hxy->Divide(fdiv);
							hxy->SetStats(0);
							parnam=(char *)model->GetParName(i);
							hxy->SetXTitle(parnam);
							parnam=(char *)model->GetParName(j);
							hxy->SetYTitle(parnam);
							hxy->GetXaxis()->SetTitleOffset(1.2);
							hxy->GetXaxis()->SetTitleSize(0.06);
							hxy->GetXaxis()->SetLabelSize(0.06);
		      			  	hxy->GetXaxis()->CenterTitle();
							hxy->GetYaxis()->SetTitleOffset(1.2);
							hxy->GetYaxis()->SetLabelSize(0.06);
							hxy->GetYaxis()->SetTitleSize(0.06);
		     			   	hxy->GetYaxis()->CenterTitle();
							hxy->Draw("cont0");
						}
						npad++;
					}
				}
			}
			c1->cd();		
			c1->Update();
        	if (!scripting) {
			theApp.Run(kTRUE);
        	}
        	else {
			sprintf(temp,"%s_margin.pdf",nestname);
			c1->SaveAs(temp);
		}
			fdiv->Delete();
		}while(0);
	}
	else if (!strcmp(temp,"posterior")){
		do {
			if (!isnested){
				printf("    No active chains\n");
				break;
			}
			if (!ismod){
				printf("    No model loaded\n");
				break;
			}
       		cout << "    Parameter > ";
        	cin >> temp;
			int par=atoi(temp)-1;
			int npfunc=model->GetNpar();
			if (par<0||par>npfunc){
				printf("   Invalid parameter %s\n",temp);
				break;
			}
			if (fix[par]){
				printf("    Parameter %d is frozen\n",par+1);
				break;
			}
      		cout << "    Number of bins in output histogram > ";
        	cin >> temp;
			int nbpost=atoi(temp);
			if (nbpost<1){
				printf("    Invalid number of bins %s\n",temp);
				break;
			}
			int npfit=0;
			for (int i=0;i<npfunc;i++){
				if (!fix[i]) npfit++;
			}
			c1->Clear();
			c1->SetTickx();
			c1->SetTicky();
			TF1 *fdiv=new TF1("fdiv","[0]",0.,1e3);
			fdiv->SetParameter(0,npchain);
			char *parnam=new char[200];
			double minval=TMath::MinElement(npchain,chains[par]);
			double maxval=TMath::MaxElement(npchain,chains[par]);
			TH1F *hpost;
			if (islogpar[par]){
				hpost=new TH1F("post","",nbpost,pow(10.,minval),pow(10.,maxval));
			}
			else {
				hpost=new TH1F("post","",nbpost,minval,maxval);
			}
			for (int np=0;np<npchain;np++){
				if (islogpar[par]){
					hpost->Fill(pow(10.,chains[par][np]));
					c1->SetLogx();
				}
				else {
					hpost->Fill(chains[par][np]);
					c1->SetLogx(0);
				}
			}
			hpost->Divide(fdiv);
			parnam=(char *)model->GetParName(par);
			hpost->SetXTitle(parnam);
			hpost->SetYTitle("Probability");
			hpost->GetXaxis()->SetTitleOffset(1.2);
 			hpost->GetXaxis()->CenterTitle();
			hpost->GetYaxis()->SetTitleOffset(1.2);
		    hpost->GetYaxis()->CenterTitle();
			hpost->SetStats(0);
			hpost->Draw();
			c1->Update();
        	if (!scripting) {
            	theApp.Run(kTRUE);
        	}
		else {
			sprintf(temp,"%s_posterior.pdf",nestname);
			c1->SaveAs(temp);
		}
			fdiv->Delete();
			hpost->Delete();
		}while(0);
	}
    else {
       printf("   Unknown command %s\n",temp);
       printf("   For a list of commands type help\n");
    }
  }
  delete [] img;
  delete [] exposure;
  delete [] sig;
  delete [] backmap;
  delete [] profile;
  delete [] bins;
  delete [] ebins;
  delete [] eprof;
  delete [] ecp;
	if (isnested){
		int npmod=model->GetNpar();
		for (int i=0;i<npmod;i++){
			delete [] chains[i];
		}
		delete [] chains;
	}
  delete model;
 // c1->Delete();
  for (int i=0;i<20;i++){
    delete [] names[i];
  }
  delete [] names;
  delete [] deltachi;
  delete [] egr;
  delete [] axes;
	delete [] backprof;
	delete [] wcs_inp;
	delete [] deprof;
	delete [] edeprof;
	delete [] scat;
	delete [] escat;
	delete [] effexp;
    delete [] medprof;
    delete [] emedprof;
	delete [] fix;
	delete [] islogpar;
}
while (0);
}

int main(int argc, char **argv) {
  proffit(argc,argv);
}
