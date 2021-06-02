/*
 *  help.h
 *  
 *
 *  Created by Dominique Eckert on 23.09.10.
 *  Copyright 2010 INAF/IASF-Milano. All rights reserved.
 *
 */
void commlist(){
    int nhelp=77;
	char **allcoms=new char*[nhelp];
	for (int i=0;i<nhelp;i++){
		allcoms[i]=new char[nhelp];
	}
    int n=0;
    allcoms[n]=(char *)"allsectors    ";n++;
	allcoms[n]=(char *)"angle2dist    ";n++;
	allcoms[n]=(char *)"backsub       ";n++;
    allcoms[n]=(char *)"box           ";n++;
    allcoms[n]=(char *)"centroidshift ";n++;
	allcoms[n]=(char *)"contour       ";n++;
    allcoms[n]=(char *)"csb           ";n++;
    allcoms[n]=(char *)"density       ";n++;
	allcoms[n]=(char *)"deproject     ";n++;
	allcoms[n]=(char *)"ellipse       ";n++;
	allcoms[n]=(char *)"ellipticity   ";n++;
	allcoms[n]=(char *)"error         ";n++;
	allcoms[n]=(char *)"fakeit        ";n++;
    allcoms[n]=(char *)"fit           ";n++;
	allcoms[n]=(char *)"fitcounts     ";n++;
	allcoms[n]=(char *)"fixpar        ";n++;
	allcoms[n]=(char *)"flcomp        ";n++;
	allcoms[n]=(char *)"flobs         ";n++;
	allcoms[n]=(char *)"flux          ";n++;
	allcoms[n]=(char *)"group         ";n++;
	allcoms[n]=(char *)"growth        ";n++;
	allcoms[n]=(char *)"help          ";n++;
    allcoms[n]=(char *)"kpc           ";n++;
	allcoms[n]=(char *)"limits        ";n++;
    allcoms[n]=(char *)"loadchains    ";n++;
    allcoms[n]=(char *)"loadfits      ";n++;
	allcoms[n]=(char *)"logx          ";n++;
	allcoms[n]=(char *)"logy          ";n++;
	allcoms[n]=(char *)"ls            ";n++;
    allcoms[n]=(char *)"margin        ";n++;
    allcoms[n]=(char *)"mediansb      ";n++;
	allcoms[n]=(char *)"model         ";n++;
    allcoms[n]=(char *)"multinest     ";n++;
	allcoms[n]=(char *)"newpar        ";n++;
	allcoms[n]=(char *)"plot          ";n++;
	allcoms[n]=(char *)"plotcounts    ";n++;
	allcoms[n]=(char *)"plotgrowth    ";n++;
	allcoms[n]=(char *)"plotgrmod     ";n++;
	allcoms[n]=(char *)"plotmod       ";n++;
	allcoms[n]=(char *)"plotmodcounts ";n++;
    allcoms[n]=(char *)"posterior     ";n++;
	allcoms[n]=(char *)"psf           ";n++;
	allcoms[n]=(char *)"profile       ";n++;
	allcoms[n]=(char *)"quit          ";n++;
	allcoms[n]=(char *)"r200          ";n++;
	allcoms[n]=(char *)"r500          ";n++;
	allcoms[n]=(char *)"readback      ";n++;
	allcoms[n]=(char *)"readexp       ";n++;
	allcoms[n]=(char *)"readimg       ";n++;
	allcoms[n]=(char *)"region        ";n++;
	allcoms[n]=(char *)"save          ";n++;
	allcoms[n]=(char *)"savedeviations";n++;
	allcoms[n]=(char *)"savemodimg    ";n++;
    allcoms[n]=(char *)"savescat      ";n++;
	allcoms[n]=(char *)"scatter       ";n++;
    allcoms[n]=(char *)"scripting     ";n++;
	allcoms[n]=(char *)"sector        ";n++;
	allcoms[n]=(char *)"sectorellipse ";n++;
	allcoms[n]=(char *)"showmod       ";n++;
	allcoms[n]=(char *)"statistics    ";n++;
	allcoms[n]=(char *)"syst          ";n++;
	allcoms[n]=(char *)"thawpar       ";n++;
    allcoms[n]=(char *)"voronoi       ";n++;
	int nlins=n/6+1;
	for (int i=n;i<nlins*6;i++){
		allcoms[i]=(char *)"              ";
	}
	printf("    List of commands: \n");
	for (int i=0;i<nlins;i++){
		printf("    %s %s %s %s %s %s\n",allcoms[0+i],allcoms[nlins+i],allcoms[2*nlins+i],allcoms[3*nlins+i],allcoms[4*nlins+i],allcoms[5*nlins+i]);
	}
	/*for (int i=0;i<39;i++){
		delete [] allcoms[i];
	}
	delete [] allcoms;*/
}

void help(char *temp){
    if (!strcmp(temp,"allsectors")){
        printf("    Plot surface-brightness profiles side by side in sectors covering the entire azimuth. The number of sectors can be adjusted.\n");
    }
	else if (!strcmp(temp,"angle2dist")){
		printf("    For a given redshift, convert an angular distance on the sky into a physical distance. Uses H0=72, Omega_m=0.27, Omega_lambda=0.73.\n");
	}
	else if (!strcmp(temp,"backsub")){
		printf("    Subtract the background from the profile. Uses the best-fit results for the background.\n");
	}
    else if (!strcmp(temp,"box")){
        printf("    Extract a brightness profile from the image in a box-shaped region. The coordinate of the center should correspond to the inner edge of the box. The width of the box and the maximum radius should be given.\n");
    }
    else if (!strcmp(temp,"centroidshift")){
        printf("    Calculate the centroid shift within a given aperture using the definition of Rasia et al. 2012.\n");
    }
	else if (!strcmp(temp,"contour")){
		printf("    Produce a 2-dimensional contour-plot (1, 2 and 3 sigma) of the chi2 function around the minimum for 2 parameters.\n");
	}	
    else if (!strcmp(temp,"csb")){
        printf("    Calculate the surface-brightness concentration parameter, which is defined as the ratio between the flux within 40 kpc radius to the flux within 400 kpc radius (see Santos et al. 2008, A&A 483, 35).\n");
    }	
    else if (!strcmp(temp,"density")){
        printf("    Compute the density profile from the deprojected profile. The conversion between count rate and emission measure must be provided by the user.\n");
    }	
	else if (!strcmp(temp,"deproject")){
		printf("    Deproject the surface-brightness profile using the procedure of Kriss et al. 1983, ApJ 272, 439, applying edge correction following McLaughlin (1999) and median smoothing to avoid jumping profiles. To estimate the error bars, simple error propagation (fast) or Monte Carlo (lengthy, but much more accurate) can be used.\n");
	}	
	else if (!strcmp(temp,"ellipse")){
		printf("    Extract the surface-brightness profile in ellipse-shaped annuli. The center can be either the image centroid (1), the surface-brightness peak (2), or any user-given input (3 and 4).\n");
	}
	else if (!strcmp(temp,"ellipticity")){
		printf("    Calculate the ellipticity of the emission using the definition of Hashimoto et al. 2007, A&A 467, 485.\n");
	}
	else if (!strcmp(temp,"error")){
		printf("    Compute the errors (upper and lower) at a given confidence level for a parameter using the MINOS algorithm.\n");
	}
    else if (!strcmp(temp,"fakeit")){
        printf("    Simulate an image using the current model. The exposure map is used to take the vignetting curve into account. In case the profile was extracted in elliptical annuli, an elliptical image is simulated. The simulated can the be analyzed in the same way as a true image.\n");
    }
	else if (!strcmp(temp,"fit")){
		printf("    Fit the model to the surface-brightness profile (see also fitcounts). To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"fitcounts")){
		printf("    Fit the model folded through the instrument response to the counts image (see also fit). To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"fixpar")){
		printf("    Fix the model parameter number n while fitting.\n");
		printf("    Usage: fixpar n\n");
	}
	else if (!strcmp(temp,"flcomp")){
		printf("    For the models with multiple components only, compute the background-subtracted model flux (in counts/sec) from radius r1 to r2 for both components.\n");
		printf("    Usage: flcomp r1 r2\n");
	}
	else if (!strcmp(temp,"flobs")){
		printf("    Compute the integrated count rate of the profile from r1 to r2.\n");
		printf("    Usage: flobs r1 r2\n");
	}
	else if (!strcmp(temp,"flux")){
		printf("    Compute the background-subtracted model flux (in counts/sec) from radius r1 to r2.\n");
		printf("    Usage: flux r1 r2\n");
	}
	else if (!strcmp(temp,"group")){
		printf("    Group the bins in the surface-brightness profile to achieve a minimum number of counts per bin or a minimum S/N per bin.\n");
	} 
	else if (!strcmp(temp,"growth")){
		printf("    Extract the growth curve of the image. The center can be either the image centroid (1), the surface-brightness peak (2), or any user-given input (3 and 4).\n");
	}
	else if (!strcmp(temp,"help")){
		printf("    Display this help.\n");
	}
    else if (!strcmp(temp,"kpc")){
        printf("    Plot the surface-brightness profile in physical distance (kiloparsec). Uses H0=72, Omega_m=0.27, Omega_lambda=0.73.\n");
    }
	else if (!strcmp(temp,"limits")){
		printf("    Set the limits for the fitting from r1 to r2 to fit only some part of the profile. By default, the whole surface-brightness profile is used.\n");
		printf("    Usage: limits r1 r2\n");
	}
    else if (!strcmp(temp,"loadfits")){
        printf("    Load a FITS file saved from a previous PROFFIT session and recover the saved data in the current session. Requires a FITS file saved with PROFFIT version 1.5 or greater.\n");
    }
	else if (!strcmp(temp,"logx")){
		printf("    Set/unset log scale on X axis (default is logarithmic).\n");
	}
	else if (!strcmp(temp,"logy")){
		printf("    Set/unset log scale on Y axis (default is logarithmic).\n");
	}
	else if (!strcmp(temp,"ls")){
		printf("    Shell command \"ls\".\n");
	}
	else if (!strcmp(temp,"model")){
		printf("    Set the surface-brightness model and give inital parameters.\n");
	}
	else if (!strcmp(temp,"newpar")){
		printf("    Set the initial value for parameter n.\n");
		printf("    Usage: newpar n value\n");
	}
	else if (!strcmp(temp,"plot")){
		printf("    Plot the surface-brightness profile. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"plotcounts")){
		printf("    Plot the number of counts per bin. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"plotgrowth")){
		printf("    Plot the growth curve. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"plotgrmod")){
		printf("    Plot the growth curve and the current background-subtracted model. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"plotmod")){
		printf("    Plot the current model. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"plotmodcounts")){
		printf("    Plot the counts profile predicted by the model. To return to the command prompt, click File->Quit ROOT.\n");
	}
	else if (!strcmp(temp,"psf")){
		printf("    Model the PSF of the instrument (Gaussian or King profiles are aloud). A probability matrix is extracted using a ray-tracing approach and convolved with the model while fitting. The number of ray-tracing photons determines the accuracy in the PSF matrix. See Appendix C of Eckert et al. 2016, arXiv:1512.03814.\n");
	}
	else if (!strcmp(temp,"profile")){
		printf("    Extract the surface-brightness profile from the image. The center can be either the image centroid (1), the surface-brightness peak (2), or any user-given input (3 and 4).\n");
	}
	else if (!strcmp(temp,"quit")){
		printf("    Quit proffit.\n");
	}
    else if (!strcmp(temp,"mediansb")){
        printf("    Calculate median the surface brightness profiles (Eckert et al. 2015, MNRAS 447,2198). If Voronoi binning has been extracted, the algorithm is the same as described in Eckert et al. 2015. Otherwise, surface brightness profiles are extracted in sectors covering the whole azimuth and the median value of the various sectors is taken at each radius.\n");
    }
	else if (!strcmp(temp,"r200")){
		printf("    Compute R_200 for a given redshift and temperature using the scaling relations of Arnaud et al. 2005, A&A 441, 893. Uses H0=70, Omega_m=0.3, Omega_lambda=0.7.\n");
	}
	else if (!strcmp(temp,"r500")){
		printf("    Compute R_500 for a given redshift and temperature using the scaling relations of Arnaud et al. 2005, A&A 441, 893. Uses H0=70, Omega_m=0.3, Omega_lambda=0.7.\n");
	}
	else if (!strcmp(temp,"readimg")){
		printf("    Read counts image file (FITS format).\n");
	}	
	else if (!strcmp(temp,"readexp")){
		printf("    Read exposure map file (FITS format).\n");
	}
	else if (!strcmp(temp,"readback")){
		printf("    Read background map file (FITS format).\n");
	}
	else if (!strcmp(temp,"region")){
		printf("    Filter out some parts of the image. The input file should be a SAODS9-compatible region file in \"image\" format. Only circular regions can be used.\n");
	}
	else if (!strcmp(temp,"save")){
		printf("    Save results of the analysis in an output file, in txt, FITS or ROOT format.\n");
	}
	else if (!strcmp(temp,"savedeviations")){
		printf("    Save the deviations of the image from the model in an output FITS file.\n");
	}
	else if (!strcmp(temp,"savemodimg")){
		printf("    Save the model image in an output FITS file.\n");
	}
	else if (!strcmp(temp,"scatter")){
		printf("    Compute the azimuthal scatter of the current profile in any number of sectors, following the definition of Vazza et al. 2010.\n");
	}
    else if (!strcmp(temp,"scripting")){
        printf("    Run PROFFIT in scripting mode. In this case, the interactive plotting window is switched off.\n");
    }	
    else if (!strcmp(temp,"sector")){
		printf("    Extract the surface-brightness profile in a sector. Zero angle corresponds to the RA axis.\n");
	}
	else if (!strcmp(temp,"sectorellipse")){
		printf("    Extract the surface-brightness profile in a sector in ellipse-shaped annuli. Zero angle corresponds to the RA axis.\n");
	}
	else if (!strcmp(temp,"showmod")){
		printf("    Show the current model.\n");
	}
	else if (!strcmp(temp,"statistics")){
		printf("    Set the statistical method for fitting (chi2 or cash, default=chi2).\n");
	}
	else if (!strcmp(temp,"syst")){
		printf("    Set the percentage of systematic error (default=0).\n");
		printf("    Usage: syst value\n");
	}
	else if (!strcmp(temp,"thawpar")){
		printf("    Set the model parameter number n free while fitting.\n");
		printf("    Usage: thawpar n\n");
	}
    else if (!strcmp(temp,"savescat")){
        printf("    Create an image of the deviations in sectors from the azimuthal average and save it to a FITS image. The number of sectors can be adjusted.\n");
    }
    else if (!strcmp(temp,"multinest")){
        printf("    Run a Bayesian analysis using the nested sampling algorithm MultiNest. For each free parameter, the task requires the user to decide between flat and logarithmic prios and provide upper and lower boundaries on the parameter values. The standard MultiNest output files will be stored with a custom prefix.\n");
    }
    else if (!strcmp(temp,"margin")){
        printf("    Draw an upper triangular plot of marginalized posterior probabilities and covariances from the current MultiNest chain.\n");
    }
    else if (!strcmp(temp,"posterior")){
        printf("    Draw the posterior probability distribution for a given parameter and a customizable number of bins.\n");
    }
    else if (!strcmp(temp,"loadchains")){
        printf("    Reload the results of a previous MultiNest run, using the prefix defined when running the multinest task.\n");
    }
    else if (!strcmp(temp,"voronoi")){
        printf("    Extract a Voronoi-binned surface brightness image using the Cappellari & Copin (2003) algorithm, with a target number of counts to be defined by the user. The Voronoi-binned image can be saved in an output FITS image and is then loaded in the current session for analysis with the mediansb task.\n");
    }
	else {
		printf("    Command %s does not exist.\n",temp);
	}
}

void allmodels(){
	printf("    Available models:\n");
	char **allmod=new char*[50];
	allmod[0]=(char *)"backfit    ";
	allmod[1]=(char *)"beta       ";
	allmod[2]=(char *)"bknbeta    ";
	allmod[3]=(char *)"bknpow     ";
	allmod[4]=(char *)"const      ";
	allmod[5]=(char *)"cuspbeta   ";
	allmod[6]=(char *)"doublebeta ";
	allmod[7]=(char *)"gausbeta   ";
	allmod[8]=(char *)"gausdbeta  ";
	allmod[9]=(char *)"power      ";
	allmod[10]=(char *)"triplebkn   ";
    allmod[11]=(char *)"triplepl   ";
	int nlins=12/5+1;
	for (int i=12;i<nlins*5;i++){
		allmod[i]=(char *)"              ";
	}
	printf("    List of commands: \n");
	for (int i=0;i<nlins;i++){
		printf("    %s %s %s %s %s\n",allmod[0+i],allmod[nlins+i],allmod[2*nlins+i],allmod[3*nlins+i],allmod[4*nlins+i]);
	}
	printf("    Other simple models can be used through the ROOT interactive user interface.\n");
}

void helpmodels(char *temp){
	if (!strcmp(temp,"backfit")){
		printf("    When a background map is provided, model to adjust the background profile (backprof) to the data:\n");
		printf("    S(r)=const+back*backprof\n");
	}
	else if (!strcmp(temp,"beta")){
		printf("    Standard beta model, plus constant for background fitting:\n");
		printf("    S(r)=norm*(1+(r/rc)^2)^(-3*beta+0.5)+const\n");
	}
	else if (!strcmp(temp,"doublebeta")){
		printf("    Double beta model, using a single beta parameter for the 2 components (see beta), plus constant for background fitting:\n");
		printf("    S(r)=norm*[(1+(r/rc1)^2)^(-3*beta+0.5)+ratio*(1+(r/rc2)^2)^(-3*beta+0.5)]+const\n");
	}
	else if (!strcmp(temp,"cuspbeta")){
		printf("    Cusp beta model, plus constant for background fitting:\n");
		printf("    S(r)=norm*(r/rs)^(-alpha)*(1+(r/rc)^2)^(-3*beta+0.5)+const\n");
	}
	else if (!strcmp(temp,"const")){
		printf("    Simple constant for background fitting\n");
	}
	else if (!strcmp(temp,"power")){
		printf("    Simple power law, plus constant for background fitting:\n");
		printf("    S(r)=norm*(r/rs)^(-alpha)+const\n");
	}
	else if (!strcmp(temp,"triplepl")){
		printf("    Triple power law joined continuously at the edges, plus constant for background fitting:\n");
        printf("    S(r)=norm*r^(-alpha1)+const, r<rc1; S(r)~r^(-alpha2)+const, rc1<r<rc2; S(r)~r^(-alpha3)+const, r>rc2.\n");
	}
	else if (!strcmp(temp,"triplebkn")){
		printf("    Same as triplepl, but projected along the line of sight. See Appendix A of Rossetti et al. 2013, A&A 556, 44.\n");
	}
	else if (!strcmp(temp,"gausbeta")){
		printf("    Beta model with a central Gaussian for fitting a central point source, plus constant for background fitting:\n");
		printf("    S(r)=norm*(1+(r/rc)^2)^(-3*beta+0.5)+normg*exp(-r^2/2*sigma^2)+const\n");
	}
	else if (!strcmp(temp,"gausdbeta")){
		printf("    Double beta model with a central Gaussian for fitting a central point source, plus constant for background fitting:\n");
		printf("    S(r)=norm*[(1+(r/rc1)^2)^(-3*beta+0.5)+ratio*(1+(r/rc2)^2)^(-3*beta+0.5)]+normg*exp(-r^2/2*sigma^2)+const\n");
	}
	else if (!strcmp(temp,"bknbeta")){
		printf("    Model for fitting density discontinuities, following Rossetti et al. 2007, A&A 463, 839\n");
		printf("    This model assumes a beta model inside the discontinuity and a power law outside. Parameters:\n");
		printf("      beta: beta parameter of the inner beta model\n");
		printf("      rc: core radius of the inner beta model\n");
		printf("      alpha2: power-law slope of the outer density profile\n");
		printf("      cutrad: discontinuity radius (in arcmin)\n");
		printf("      norm: normalization of the beta model\n");
		printf("      jump: density jump at the discontinuity\n");
		printf("      const: constant for background fitting\n");
	}
	else if (!strcmp(temp,"bknpow")){
		printf("    Model for fitting density discontinuities, following Owers et al. 2009, ApJ 704, 1349\n");
		printf("    This model assumes power-law density profiles inside and outside the discontinuity. Parameters:\n");
		printf("      alpha1: power-law slope of the inner density profile\n");
		printf("      alpha2: power-law slope of the outer density profile\n");
		printf("      cutrad: discontinuity radius (in arcmin)\n");
		printf("      norm: normalization of the inner component\n");
		printf("      jump: density jump at the discontinuity\n");
		printf("      const: constant for background fitting\n");
	}
	else {
		printf("    Model %s does not exist\n",temp);
	}
}
