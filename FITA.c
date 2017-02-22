/*
 * FIA.c
 * Fluorescence Intermittency Trajectory Analyzer
 * Copyright 2017 Anthony Ruth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */
#include <stdio.h>
#include <math.h>
#include "aruthmath.h"
#include "aruthio.h"
#include <gtk/gtk.h>
#include <fftw3.h>

void bin(double *f, int numinputpoints,  double *binlocations,int *output,int numoutputpoints);
void logarithmicaverage(double *x, double *y, int numinputpoints,  double *binlocations, double *output, int numoutputpoints);
void logarithmicbin(double *f, int numinputpoints,  double *binlocations,double *output,int numoutputpoints,int dividebybinsize);


double *datax;
double *datay;
int numdatapoints = 0;

double *trajectoryx;
double *trajectoryy;
int numtrajectorypoints = 0;

//int *onTimeBins;
//int *offTimeBins;
double *onTimeBins;
double *offTimeBins;
double *onTimes;
double *offTimes;
int numonoffbins;


double *intensities;
int *intensityBins;
int numintensitybins;

double *powerspectrum;
double *psdfrequencies;
int numpsdbins;

char *savedirectory = ".";

GtkWidget *TopLevel;
GtkWidget *directorylabel;

//Global
GtkWidget *trimLowerSpinButton;
GtkWidget *trimUpperSpinButton;
GtkWidget *rebinningSpinButton;
GtkWidget *scalingSpinButton;

//Threshold
GtkWidget *thresholdValueSpinButton;
GtkWidget *onoffBinsSpinButton;

//Intensity
GtkWidget *intensityBinsSpinButton;
GtkWidget *intensityFitFunctionComboBox;
GtkWidget *intensityPeaksSpinButton;
GtkWidget *displayIntensityFitsCheckButton;
GtkWidget *intensityFitQualityOutput;

//Power Spectrum
GtkWidget *psdBinsSpinButton;
GtkWidget *fourierFilterLowerSpinButton;
GtkWidget *fourierFilterUpperSpinButton;


int numspinbuttons = 11;
gulong *spin_handlers;
GtkWidget **spinbuttons;


FILE *trajectoryAndDistributionPlotPipe;
FILE *onoffDistributionPlotPipe;
FILE *PSDPlotPipe;


void globalCalcs()
{
	if(trajectoryx != NULL)
	{
		free(trajectoryx);
		free(trajectoryy);
	}
	int trimlower = gtk_spin_button_get_value(trimLowerSpinButton);
	int trimupper = gtk_spin_button_get_value(trimUpperSpinButton);
	int rebinning =  gtk_spin_button_get_value(rebinningSpinButton);
	double scaling = gtk_spin_button_get_value(scalingSpinButton);
	
	numtrajectorypoints = ceil((trimupper-trimlower)/(double)rebinning);
	trajectoryx = malloc(numtrajectorypoints*sizeof(double));
	trajectoryy = malloc(numtrajectorypoints*sizeof(double));
	
	int i;
	int j;
	double starttime = *(datax+trimlower);
	for(i=0;i<numtrajectorypoints;i++)
	{
		*(trajectoryx+i) = scaling*(*(datax+i*rebinning)-starttime);
		*(trajectoryy+i) = 0;
		for(j = 0;j<rebinning;j++)
		*(trajectoryy+i) += *(datay+trimlower+i*rebinning+j);
	}
	
}

void thresholdCalcs()
{
	if(onTimeBins != NULL)
	{
		free(onTimeBins);
		free(offTimeBins);
		free(onTimes);
		free(offTimes);

	}
	double thresholdValue = gtk_spin_button_get_value(thresholdValueSpinButton);
	numonoffbins = gtk_spin_button_get_value(onoffBinsSpinButton);
//	onTimeBins = malloc(numonoffbins*sizeof(int));
//	offTimeBins = malloc(numonoffbins*sizeof(int));
	onTimeBins = malloc(numonoffbins*sizeof(double));
	offTimeBins = malloc(numonoffbins*sizeof(double));

	onTimes = malloc(numonoffbins*sizeof(double));
	offTimes = malloc(numonoffbins*sizeof(double));
	
	int duration=1;
	int state; //1 on, 0 off
	int oldstate;
	
	double ons[numtrajectorypoints];
	double offs[numtrajectorypoints];
	int numontimes=0;
	int numofftimes=0;
	
	int i;
	for(i = 0;i<numtrajectorypoints;i++)
	{
		state = *(trajectoryy+i) >= thresholdValue;
		if(i != 0)
		{
		//printf("Old state %d new state %d\n",oldstate,state);
		if(state == oldstate)
		duration++;
		else
		{
		if(oldstate == 0)
		{
		offs[numofftimes] = *(trajectoryx+i) - *(trajectoryx+i-duration);
		//printf("Off time %g\n",offs[numofftimes]);
		duration=1;
		numofftimes++;
		}
		if(oldstate == 1)
		{
		ons[numontimes] = *(trajectoryx+i) - *(trajectoryx+i-duration);
		//printf("On time %g\n",ons[numofftimes]);
		duration=1;
		numontimes++;
		}	
		}
		}
		oldstate = state;
	}
	//bin(ons,numontimes,onTimes,onTimeBins,numonoffbins);
	//bin(offs,numofftimes,offTimes,offTimeBins,numonoffbins);
	logarithmicbin(ons,numontimes,onTimes,onTimeBins,numonoffbins,1);
	logarithmicbin(offs,numofftimes,offTimes,offTimeBins,numonoffbins,1);
}
FILE *openpipe()
    {
    
    FILE *pipe = popen("gnuplot","w");    
    return pipe;
    }

void trajectoryAndDistributionPlot()
{
	char *basefile = 
"set multiplot layout 1, 2 title 'Trajectory and Intensity Distribution'\n"
"set xlabel 'Time (s)'\n"
//"unset ytics\n"
"set key noautotitles\n"
"set rmargin 0\n"
"ntics = 4\n"
"stats '%s/trajectory.dat' using 1 name 'x' nooutput\n"
"set xtics x_max/ntics\n"
"plot '%s/trajectory.dat' w lines lc rgb 'black', %g w lines\n"
"set lmargin 1\n"
"unset rmargin\n"
"set xlabel 'Occurance'\n"
"set y2label 'Counts/bin'\n"
"set y2tics\n"
"unset ytics\n"
"stats '%s/intensitydistribution.dat' using 2 name 'v' nooutput\n"
"set logscale x\n"
"set xtics 1,v_max**(1.0/ntics),v_max\n"
"plot '%s/intensitydistribution.dat' u 2:1 axes x1y2 w lines lc rgb 'black',%g axes x1y2 w lines\n"
"reset\n";

		double thresholdValue = gtk_spin_button_get_value(thresholdValueSpinButton);
	char filename[1000];
	sprintf(filename,"%s/trajectoryanddistributionplotscript",savedirectory);
	FILE *outfile = fopen(filename,"w");
	fprintf(outfile,basefile,savedirectory,savedirectory,thresholdValue,savedirectory,savedirectory,thresholdValue);
	fclose(outfile);
	char command[1000];
	sprintf(command,"load '%s'\n",filename);
	fprintf(trajectoryAndDistributionPlotPipe,command);
	fflush(trajectoryAndDistributionPlotPipe);
	
	
	}
	
	
	void onoffDistributionPlot()
{
	char *basefile = 
"set logscale x\n"
"set logscale y\n"
"set xlabel 'On/Off Time'\n"
"set ylabel 'Occurance'\n"
"b=-1\n"
"d=-1\n"
"on(x) = a*x + b\n"
"off(x) = c*x + d\n"
"fit on(x) '%s/ondistribution.dat' using (log($1)):(log($2)) via a, b\n"
"fit off(x) '%s/offdistribution.dat' using (log($1)):(log($2)) via c, d\n"
"on(x) = exp(b)*x**a\n"
"off(x) = exp(d)*x**c\n"
"t1 = 'On Fit m='.sprintf('%%.3g',a)\n"
"t2 = 'Off Fit m='.sprintf('%%.3g',c)\n"
"plot '%s/ondistribution.dat' title 'On Distribution','%s/offdistribution.dat' title 'Off Distribution', on(x) title t1, off(x) title t2\n"
"reset\n";

	char filename[1000];
	sprintf(filename,"%s/onoffdistributionplotscript",savedirectory);
	FILE *outfile = fopen(filename,"w");
	fprintf(outfile,basefile,savedirectory,savedirectory,savedirectory,savedirectory);
	fclose(outfile);
	char command[1000];
	sprintf(command,"load '%s'\n",filename);
	fprintf(onoffDistributionPlotPipe,command);
	fflush(onoffDistributionPlotPipe);
	
	
	}

void intensityDistributionCalcs()
{
	numintensitybins = gtk_spin_button_get_value(intensityBinsSpinButton);
	intensities = malloc(numintensitybins*sizeof(double));
	intensityBins = malloc(numintensitybins*sizeof(int));
	bin(trajectoryy,numtrajectorypoints,intensities,intensityBins,numintensitybins);
	
}

void psdCalcs()
{
	fftw_complex *out;
    fftw_plan p;
    
    
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numtrajectorypoints);
    double in[numtrajectorypoints];
    int i;
	for(i=0;i<numtrajectorypoints;i++)
	{
	in[i]=trajectoryy[i];	
	//printf("Input point %d is %g\n",i,in[i]);
	}	
    
    p = fftw_plan_dft_r2c_1d(numtrajectorypoints, in, out, FFTW_ESTIMATE);
    
    fftw_execute(p); /* repeat as needed */
    
    //Fourier Filtering should go here
    
    double powerspectruminitial[numtrajectorypoints/2];
    double psdfrequenciesinitial[numtrajectorypoints/2];
    
    double basefrequency = 1/ (*(trajectoryx+numtrajectorypoints-1)-*(trajectoryx))/2.0;
	
    for(i=0;i<numtrajectorypoints/2;i++)
    {
		//printf("Point %d components %g +  i*%g\n",i,out[0][i],out[1][i]);
		//printf("Input point %d is %g\n",i,in[i]);
		*(powerspectruminitial+i) = pow(out[0][i],2) + pow(out[1][i],2);//The square of the real and imaginary parts
		*(psdfrequenciesinitial+i) = i*basefrequency;
	}
    
    numpsdbins = gtk_spin_button_get_value(psdBinsSpinButton);

    powerspectrum = malloc(numpsdbins*sizeof(double));
    psdfrequencies = malloc(numpsdbins*sizeof(double));
    logarithmicaverage(psdfrequenciesinitial,powerspectruminitial,numtrajectorypoints/2,psdfrequencies,powerspectrum,numpsdbins);
    
    
    fftw_destroy_plan(p);
    fftw_free(out);
    
}

	void psdPlot()
{
	char *basefile = 
"set logscale x\n"
"set logscale y\n"
"set xlabel 'Frequency'\n"
"set ylabel 'PSD (counts^2/hz)'\n"
"b=-1\n"
"psd(x) = a*x + b\n"
"fit psd(x) '%s/psd.dat' using (log($1)):(log($2)) via a, b\n"
"psd(x) = exp(b)*x**a\n"
"t1 = 'PSD Fit m='.sprintf('%%.3g',a)\n"
"plot '%s/psd.dat' title 'Calculated PSD', psd(x) title t1\n"
"reset\n";

	char filename[1000];
	sprintf(filename,"%s/psdplotscript",savedirectory);
	FILE *outfile = fopen(filename,"w");
	fprintf(outfile,basefile,savedirectory,savedirectory);
	fclose(outfile);
	char command[1000];
	sprintf(command,"load '%s'\n",filename);
	fprintf(PSDPlotPipe,command);
	fflush(PSDPlotPipe);
	
	}


void bin(double *f, int numinputpoints,  double *binlocations,int *output,int numoutputpoints)
{
	//This function is passed an array of values and it creates a distribution function of that data
	//The array does not need to be sorted.
	//The first and last bin will always have at least one entry
	double epsilon = 1e-6;
	double lowerbound = listmin(f,numinputpoints);
	double stepsize = (listmax(f,numinputpoints)-lowerbound)/(numoutputpoints-1);
	int i,j;
	
	printf("Binning %d data points into %d bins Smallest value %g, largest value %g\n",numinputpoints,numoutputpoints,listmin(f,numinputpoints),listmax(f,numinputpoints));
	for(j = 0;j<numoutputpoints;j++)
	{
		*(output+j)=0;
		*(binlocations+j)=lowerbound + stepsize*j;
		//printf("Looking for points between %g and %g\n",lowerbound + stepsize*j,lowerbound + stepsize*(j+1));
		for(i=0;i<numinputpoints;i++)
		if(lowerbound-epsilon + stepsize*j<*(f+i) && *(f+i) < lowerbound-epsilon + stepsize*(j+1))
		(*(output+j))++;
	}
}

void logarithmicbin(double *f, int numinputpoints,  double *binlocations,double *output,int numoutputpoints,int dividebybinsize)
{
	//This function is passed an array of values and it creates a distribution function of that data
	//The array does not need to be sorted.
	//The first and last bin will always have at least one entry
	//This function expects to take in integer-valued function f. If it is not integer-valued, then the binsize factor should be modified.
	double epsilon = 1e-6;
	double lowerbound = listmin(f,numinputpoints);
	double upperbound = listmax(f,numinputpoints);
	double stepsize = exp(log(upperbound/lowerbound)/(numoutputpoints-1)); //lowerbound * stepsize ^(n-1) = upperbound
	int i,j;
	
	printf("Binning %d data points into %d bins Smallest value %g, largest value %g\n",numinputpoints,numoutputpoints,lowerbound,upperbound);
	for(j = 0;j<numoutputpoints;j++)
	{
		*(output+j)=0;
		*(binlocations+j)=lowerbound * pow(stepsize,j);
		//printf("Looking for points between %g and %g\n",lowerbound + stepsize*j,lowerbound + stepsize*(j+1));
		for(i=0;i<numinputpoints;i++)
		{
		if(*(binlocations+j)-epsilon<*(f+i) && *(f+i) <= *(binlocations+j)*stepsize-epsilon)
		{
			//printf("binsize factor for bin %d which is from %g to %g is %g\n",i,*(binlocations+j)-epsilon,*(binlocations+j)*stepsize-epsilon,1.0/ (floor(*(binlocations+j)*stepsize-epsilon) - floor(*(binlocations+j)-epsilon)));
		if(dividebybinsize)
		*(output+j) += 1.0/ (floor(*(binlocations+j)*stepsize-epsilon) - floor(*(binlocations+j)-epsilon));
		else
		(*(output+j))++;
		}
	}
	}
}

void logarithmicaverage(double *x, double *y, int numinputpoints,  double *binlocations, double *output, int numoutputpoints)
{
	//This function is passed an array of values and it logarithmically averages the values y along the axis x
	//The array does need to be sorted.
	//The first and last bin will always have at least one entry
	//printf("In logarithmic average function\n");
	double epsilon = 1e-6;
	//double lowerbound = listmin(x,numinputpoints);
	double lowerbound = *(x+1);//This should be nonzero
	double upperbound = listmax(x,numinputpoints);
	double stepsize = exp(log(upperbound/lowerbound)/(numoutputpoints-1)); //lowerbound * stepsize ^(n-1) = upperbound
	int i,j;
	
	printf("Averaging %d data points into %d bins Smallest x %g, largest x %g\n",numinputpoints,numoutputpoints,lowerbound,upperbound);
	int numpointsdetected;
	for(j = 0;j<numoutputpoints;j++)
	{
		*(output+j)=0;
		*(binlocations+j)=lowerbound * pow(stepsize,j);
		//printf("Looking for points between %g and %g\n",lowerbound + stepsize*j,lowerbound + stepsize*(j+1));
		numpointsdetected = 0;
		for(i=0;i<numinputpoints;i++)
		if(*(binlocations+j)-epsilon<*(x+i) && *(x+i) <= *(binlocations+j)*stepsize-epsilon)
		{
		(*(output+j)) += *(y+i);
		numpointsdetected++;
		}
		if(numpointsdetected != 0)
		(*(output+j)) /= numpointsdetected;
	}
}

/*void plottrajectory(char *filename,FILE *pipe)
    {
    
    char *basefile = 
    
    
    
    fprintf(pipe,command);
    fflush(pipe);
    }
*/
void calc()
{
	printf("Running all calculations\n");
	globalCalcs();
	psdCalcs();
	thresholdCalcs();
	intensityDistributionCalcs();
	
	int i;
	
	char filename[1000];
	sprintf(filename,"%s/trajectory.dat",savedirectory);
	printf("Saveing %s\n",filename);
	FILE *trajectoryfile = fopen(filename,"w");
	for(i=0;i<numtrajectorypoints;i++)
	fprintf(trajectoryfile,"%g %g\n",*(trajectoryx+i),*(trajectoryy+i));
	fclose(trajectoryfile);
	
	sprintf(filename,"%s/intensitydistribution.dat",savedirectory);
	printf("Saveing %s\n",filename);
	FILE *intensityfile = fopen(filename,"w");
	for(i=0;i<numintensitybins;i++)
	fprintf(intensityfile,"%g %d\n",*(intensities+i),*(intensityBins+i));
	fclose(intensityfile);
	
	sprintf(filename,"%s/ondistribution.dat",savedirectory);
	printf("Saveing %s\n",filename);
	FILE *onfile = fopen(filename,"w");
	for(i=0;i<numonoffbins;i++)
	fprintf(onfile,"%g %g\n",*(onTimes+i),*(onTimeBins+i));
	fclose(onfile);
	
	sprintf(filename,"%s/offdistribution.dat",savedirectory);
	printf("Saveing %s\n",filename);
	FILE *offfile = fopen(filename,"w");
	for(i=0;i<numonoffbins;i++)
	fprintf(offfile,"%g %g\n",*(offTimes+i),*(offTimeBins+i));
	fclose(offfile);
	
	sprintf(filename,"%s/psd.dat",savedirectory);
	printf("Saveing %s\n",filename);
	FILE *psdfile = fopen(filename,"w");
	for(i=0;i<numpsdbins;i++)
	fprintf(psdfile,"%g %g\n",*(psdfrequencies+i),*(powerspectrum+i));
	fclose(psdfile);
	
	
}

gboolean runAll()
{
	printf("Running all\n");
	calc();
	printf("Finished all calcs now plotting\n");
	trajectoryAndDistributionPlot();
	onoffDistributionPlot();
	psdPlot();
	return FALSE;
}	

GtkWidget * global()
{
	//Returns the frame for global settings
	double T_max=0;
	
	
	
	GtkWidget *trimLowerBox = gtk_hbox_new(FALSE,0);
	trimLowerSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0, 0, numdatapoints, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(trimLowerSpinButton,"Cuts Off Data Points Prior to This");
	GtkWidget *trimLowerText = gtk_label_new("Lower Trim: "); 
	gtk_box_pack_start (GTK_BOX(trimLowerBox),trimLowerText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(trimLowerBox),trimLowerSpinButton, FALSE, TRUE, 10);

	GtkWidget *trimUpperBox = gtk_hbox_new(FALSE,0);
	trimUpperSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (numdatapoints, 0, numdatapoints, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(trimUpperSpinButton,"Cuts Off Data Points After This");
	GtkWidget *trimUpperText = gtk_label_new("Upper Trim: "); 
	gtk_box_pack_start (GTK_BOX(trimUpperBox),trimUpperText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(trimUpperBox),trimUpperSpinButton, FALSE, TRUE, 10);

	GtkWidget *rebinningBox = gtk_hbox_new(FALSE,0);
	rebinningSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (1, 0, numdatapoints, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(rebinningSpinButton,"Rebins the data by summing over this many bins");
	GtkWidget *rebinningText = gtk_label_new("Rebinning Size: "); 
	gtk_box_pack_start (GTK_BOX(rebinningBox),rebinningText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(rebinningBox),rebinningSpinButton, FALSE, TRUE, 10);

	GtkWidget *scalingBox = gtk_hbox_new(FALSE,0);
	scalingSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (1, 0, 1000000000, 0.00001, 0.000001, 0.0)),5,0);
	gtk_widget_set_tooltip_text(scalingSpinButton,"Scales the x axis by this number. This should be the time/bin in the original data.");
	GtkWidget *scalingText = gtk_label_new("Time Scaling (s/bin): "); 
	gtk_box_pack_start (GTK_BOX(scalingBox),scalingText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(scalingBox),scalingSpinButton, FALSE, TRUE, 10);

	GtkWidget *globalBox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(globalBox),trimLowerBox, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(globalBox),trimUpperBox, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(globalBox),rebinningBox, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(globalBox),scalingBox, FALSE, TRUE, 10);

	return globalBox;
}

GtkWidget * threshold()
{
	double I_max=0;
	
	GtkWidget *thresholdValueBox = gtk_hbox_new(FALSE,0);
	thresholdValueSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (I_max/2, 0, I_max, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(thresholdValueSpinButton,"Threshold to use for on/off analysis");
	GtkWidget *thresholdValueText = gtk_label_new("Threshold: "); 
	gtk_box_pack_start (GTK_BOX(thresholdValueBox),thresholdValueText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(thresholdValueBox),thresholdValueSpinButton, FALSE, TRUE, 10);

	GtkWidget *onoffBinsBox = gtk_hbox_new(FALSE,0);
	onoffBinsSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (pow(numdatapoints,0.5), 0, numdatapoints, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(onoffBinsSpinButton,"Number of bins to use in on/off times distribution function");
	GtkWidget *onoffBinsText = gtk_label_new("On/Off Time Bins"); 
	gtk_box_pack_start (GTK_BOX(onoffBinsBox),onoffBinsText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(onoffBinsBox),onoffBinsSpinButton, FALSE, TRUE, 10);

	GtkWidget *thresholdBox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(thresholdBox),thresholdValueBox, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(thresholdBox),onoffBinsBox, FALSE, TRUE, 10);
	
	return thresholdBox;
}

GtkWidget * intensityDistribution()
{
	double I_max = 0;
	
	GtkWidget *intensityBinsBox = gtk_hbox_new(FALSE,0);
	intensityBinsSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (sqrt(I_max), 0, I_max, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(intensityBinsSpinButton,"Number of bins to use for intensity distribution");
	GtkWidget *intensityBinsText = gtk_label_new("Intensity Bins: "); 
	gtk_box_pack_start (GTK_BOX(intensityBinsBox),intensityBinsText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(intensityBinsBox),intensityBinsSpinButton, FALSE, TRUE, 10);

	GtkWidget *intensityFitFunctionBox = gtk_hbox_new(FALSE,0);
	intensityFitFunctionComboBox = gtk_combo_box_text_new ();
	gtk_combo_box_text_append (intensityFitFunctionComboBox,NULL,"Gaussian");
	gtk_combo_box_text_append (intensityFitFunctionComboBox,NULL,"Poisson");
	gtk_combo_box_set_active(intensityFitFunctionComboBox,0);
	gtk_widget_set_tooltip_text(intensityFitFunctionComboBox,"Function to use for peak fitting");
	GtkWidget *intensityFitFunctionText = gtk_label_new("Fit Function: "); 
	gtk_box_pack_start (GTK_BOX(intensityFitFunctionBox),intensityFitFunctionText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(intensityFitFunctionBox),intensityFitFunctionComboBox, FALSE, TRUE, 10);

	GtkWidget *intensityPeaksBox = gtk_hbox_new(FALSE,0);
	intensityPeaksSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (2, 0, 10, 1, 1, 0.0)),1,0);
	gtk_widget_set_tooltip_text(intensityPeaksSpinButton,"Number of peaks to use in fit");
	GtkWidget *intensityPeaksText = gtk_label_new("Peaks: "); 
	gtk_box_pack_start (GTK_BOX(intensityPeaksBox),intensityPeaksText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(intensityPeaksBox),intensityPeaksSpinButton, FALSE, TRUE, 10);

	displayIntensityFitsCheckButton = gtk_check_button_new_with_label("Display Fit");
	gtk_toggle_button_set_active(displayIntensityFitsCheckButton,TRUE);
	
	GtkWidget *intensityDistributionBox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(intensityDistributionBox),intensityBinsBox, FALSE, TRUE, 10);
	//gtk_box_pack_start (GTK_BOX(intensityDistributionBox),intensityFitFunctionBox, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(intensityDistributionBox),intensityPeaksBox, FALSE, TRUE, 10);
	//gtk_box_pack_start (GTK_BOX(intensityDistributionBox),displayIntensityFitsCheckButton, FALSE, TRUE, 10);
	//These two have not been connected yet


	return intensityDistributionBox;
}

GtkWidget * psd()
{
	double T_max = 0;

	if(datax != NULL)
	T_max = *(datax+numdatapoints-1);

	GtkWidget *psdBinsBox = gtk_hbox_new(FALSE,0);
	psdBinsSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (pow(numdatapoints,0.5), 0, numdatapoints, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(psdBinsSpinButton,"Number of bins to use in psd calculation");
	GtkWidget *psdBinsText = gtk_label_new("PSD Bins"); 
	gtk_box_pack_start (GTK_BOX(psdBinsBox),psdBinsText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(psdBinsBox),psdBinsSpinButton, FALSE, TRUE, 10);


	GtkWidget *fourierFilterLowerBox = gtk_hbox_new(FALSE,0);
	fourierFilterLowerSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0, 0, T_max, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(fourierFilterLowerSpinButton,"Low Frequency Fourier Filtering Cutoff");
	GtkWidget *fourierFilterLowerText = gtk_label_new("Low Cut:"); 
	gtk_box_pack_start (GTK_BOX(fourierFilterLowerBox),fourierFilterLowerText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(fourierFilterLowerBox),fourierFilterLowerSpinButton, FALSE, TRUE, 10);

	GtkWidget *fourierFilterUpperBox = gtk_hbox_new(FALSE,0);
	fourierFilterUpperSpinButton =  gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (T_max, 0, T_max, 1, 10, 0.0)),1,0);
	gtk_widget_set_tooltip_text(fourierFilterLowerSpinButton,"High Frequency Fourier Filtering Cutoff");
	GtkWidget *fourierFilterUpperText = gtk_label_new("High Cut:"); 
	gtk_box_pack_start (GTK_BOX(fourierFilterUpperBox),fourierFilterUpperText, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(fourierFilterUpperBox),fourierFilterUpperSpinButton, FALSE, TRUE, 10);

	GtkWidget *powerspectrum = gtk_vbox_new(FALSE,0);
	//gtk_box_pack_start(GTK_BOX(powerspectrum),fourierFilterLowerBox,FALSE,TRUE,10);
	//gtk_box_pack_start(GTK_BOX(powerspectrum),fourierFilterUpperBox,FALSE,TRUE,10);
	gtk_box_pack_start(GTK_BOX(powerspectrum),psdBinsBox,FALSE,TRUE,10);
//Fourier Filtering has not been enabled yet
	
	return powerspectrum;
}

void readjust()
{
	double I_max = 0;
	double T_max = 0;
	if(datay != NULL)
	I_max = listmax(datay,numdatapoints);
	
	if(datax != NULL)
	T_max = *(datax+numdatapoints-1);
	
	printf("Starting Adjustments\n");
	int i;
	GtkWidget *spin[] =
	{trimLowerSpinButton,trimUpperSpinButton,rebinningSpinButton,scalingSpinButton,thresholdValueSpinButton,onoffBinsSpinButton,
		intensityBinsSpinButton,intensityPeaksSpinButton,fourierFilterLowerSpinButton,fourierFilterUpperSpinButton,psdBinsSpinButton};
	
	for(i=0;i<numspinbuttons;i++)
	 g_signal_handler_block (spin[i], spin_handlers[i]);
	 printf("Here\n");
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(trimLowerSpinButton),gtk_adjustment_new (0, 0, numdatapoints, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(trimUpperSpinButton),gtk_adjustment_new (numdatapoints, 0, numdatapoints, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(rebinningSpinButton),gtk_adjustment_new (1, 0, numdatapoints, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(scalingSpinButton),gtk_adjustment_new (1, 0, 1000000000, 0.00001, 0.000001, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(thresholdValueSpinButton),gtk_adjustment_new (I_max/2, 0, I_max, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(onoffBinsSpinButton),gtk_adjustment_new (pow(numdatapoints,0.5), 0, numdatapoints, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(intensityBinsSpinButton),gtk_adjustment_new (sqrt(I_max), 0, I_max, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(intensityPeaksSpinButton),gtk_adjustment_new (2, 0, 10, 1, 1, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(psdBinsSpinButton),gtk_adjustment_new (pow(numdatapoints,0.5), 0, numdatapoints, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(fourierFilterLowerSpinButton),gtk_adjustment_new (0, 0, T_max, 1, 10, 0.0));
	 gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(fourierFilterUpperSpinButton),gtk_adjustment_new (T_max, 0, T_max, 1, 10, 0.0));
	for(i=0;i<numspinbuttons;i++)
	 g_signal_handler_unblock (spin[i], spin_handlers[i]);

	printf("Done With Adjustments\n");
}

void readInputFile(char *filename)
{
	if(datax != NULL)
	free(datax);
	if(datay != NULL)
	free(datay);
	
	numdatapoints = countlines(filename);
	int twocolumn = 0;
	if(countcolumns(filename) == 2)
	twocolumn = 1;
	
	datax = malloc(numdatapoints*sizeof(double));
	datay = malloc(numdatapoints*sizeof(double));
	
	FILE *infile = fopen(filename,"r");
	int index;
	for(index = 0;index<numdatapoints;index++)
	{
		if(twocolumn)
		{
		fscanf(infile,"%lf %lf",datax+index,datay+index);
		}
		else
		{
		fscanf(infile,"%lf",datay+index);
		*(datax+index) = index;
		}
	}
}

	gboolean load()
{
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	gint res;

dialog = gtk_file_chooser_dialog_new ("Load Data From File",
                                      TopLevel,
                                     GTK_FILE_CHOOSER_ACTION_OPEN,
                                      "Cancel",
                                      GTK_RESPONSE_CANCEL,
                                      "Load",
                                      GTK_RESPONSE_ACCEPT,
                                      NULL);
chooser = GTK_FILE_CHOOSER (dialog);

gtk_file_chooser_set_do_overwrite_confirmation (chooser, TRUE);


res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    filename = gtk_file_chooser_get_filename (chooser);
	readInputFile(filename);
	readjust();
	//calc();
	runAll();
  }

	gtk_widget_destroy (dialog);
	return TRUE;
}

	gboolean save()
{
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	gint res;

dialog = gtk_file_chooser_dialog_new ("Load Data From File",
                                      TopLevel,
                                     GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER,
                                      "Cancel",
                                      GTK_RESPONSE_CANCEL,
                                      "Save",
                                      GTK_RESPONSE_ACCEPT,
                                      NULL);
chooser = GTK_FILE_CHOOSER (dialog);

gtk_file_chooser_set_do_overwrite_confirmation (chooser, TRUE);


res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
  {
    savedirectory = gtk_file_chooser_get_filename (chooser);
	gtk_label_set_text(directorylabel,savedirectory);
  }

	gtk_widget_destroy (dialog);
	return FALSE;
}


void connectAll()
{
	//printf("Location of trimLowerSpinButton %d\n",trimLowerSpinButton);
	//spinbuttons = 
	GtkWidget *spin[] =
	{trimLowerSpinButton,trimUpperSpinButton,rebinningSpinButton,scalingSpinButton,thresholdValueSpinButton,onoffBinsSpinButton,
		intensityBinsSpinButton,intensityPeaksSpinButton,fourierFilterLowerSpinButton,fourierFilterUpperSpinButton,psdBinsSpinButton};
	spinbuttons = spin;
	
		//printf("Location of spinbutton[0] %d\n",(spinbuttons[0]));	
		spin_handlers = malloc(numspinbuttons*sizeof(gulong));
		
	int i;
		for(i=0;i<numspinbuttons;i++)
		spin_handlers[i] = g_signal_connect_swapped(spinbuttons[i], "value-changed", G_CALLBACK (runAll),NULL);
	
}

int main(int argc, char **argv)
{
	gtk_init(NULL,NULL);

	trajectoryAndDistributionPlotPipe = openpipe();
	onoffDistributionPlotPipe = openpipe();
	PSDPlotPipe = openpipe();


	GtkWidget *notebook = gtk_notebook_new();
	gtk_notebook_append_page(notebook,global(),gtk_label_new("Global"));
	gtk_notebook_append_page(notebook,threshold(),gtk_label_new("Thresholding"));
	gtk_notebook_append_page(notebook,intensityDistribution(),gtk_label_new("Intensity Distribution"));
	gtk_notebook_append_page(notebook,psd(),gtk_label_new("Power Spectrum"));

	
	GtkWidget *loadbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(loadbutton),gtk_image_new_from_stock(GTK_STOCK_OPEN,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(loadbutton, "clicked", G_CALLBACK (load),NULL);

	GtkWidget *savebutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(savebutton),gtk_image_new_from_stock(GTK_STOCK_SAVE,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(savebutton, "clicked", G_CALLBACK (save),NULL);

	GtkWidget *savetext = gtk_label_new("Save Directory:");
	directorylabel = gtk_label_new(savedirectory);

	GtkWidget *topBox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start(topBox,loadbutton,FALSE,TRUE,10);
	gtk_box_pack_start(topBox,savebutton,FALSE,TRUE,10);
	gtk_box_pack_start(topBox,savetext,FALSE,TRUE,10);
	gtk_box_pack_start(topBox,directorylabel,TRUE,TRUE,10);

	
	GtkWidget *mainBox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(mainBox),topBox,TRUE,TRUE,10);
	gtk_box_pack_start(GTK_BOX(mainBox),notebook,TRUE,TRUE,10);
	
	TopLevel = gtk_window_new (GTK_WINDOW_TOPLEVEL); 
	gtk_container_add (GTK_CONTAINER (TopLevel), mainBox);
	gtk_widget_show_all(TopLevel);

//GtkWidget *intensityFitFunctionComboBox;
//GtkWidget *displayIntensityFitsCheckButton;
	
	connectAll();
	
	gtk_main();
}

