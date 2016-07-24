/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUBFILT: $Revision: 1.22 $ ; $Date: 2012/11/28 22:13:13 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "recstalta.c"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                               ",
" SUSTALTA - Apply STA/LTA Algorithm                         	",
" Guided by :SUBFILT - apply Butterworth bandpass filter        ",
"                                   							",
" sustalta <stdin >stdout [optional parameters]			        ",
"                            							        ",
" Required parameters:		                    				",
" 	- if dt is not set in header, then dt is mandatory	        ",
" 	- if ns is not correct in header, then ns is mandatory      ",
"                                                               ",
" Optional parameters: (nyquist calculated internally)		    ",
" 	verbose=0		=1 for filter design info 	                ",
"                                                               ",
" 	dt = (from header)	time sampling interval (sec)	        ",
"       ns = (from header)      number of samples per trace     ",
"       nlta =  ns / 10         length of long time avg.        ",
"                               in samples                      ",
"       nsta =  nlta / 10       length of short time avg.       ",
"                               in samples                      ",
NULL};

/* Credits:
 *	CWP: Dave Hale c. 1993 for bf.c subs and test drivers
 *	CWP: Jack K. Cohen for su wrapper c. 1993
 *      SEAM Project: Bruce Verwest 2009 added explicit pole option
 *                    in a program called "subfiltpole"
 *      CWP: John Stockwell (2012) combined Bruce Verwests changes
 *           into the original subfilt.
 *
 * Caveat: zerophase will not do good if trace has a spike near
 *	   the end.  One could make a try at getting the "effective"
 *	   length of the causal filter, but padding the traces seems
 *	   painful in an already expensive algorithm.
 *
 *
 * Theory:
 * The 
 *
 * Trace header fields accessed: ns, dt, trid
 */
/**************** end self doc ***********************************/

/* prototypes of functions used internally */

void recstalta(double *a, double *charfct, int ndat, int nsta, int nlta);
void fputdata(FILE *fileptr, FILE *headerptr, double *outdata, int nt);
void fputdata3c(FILE *fileptr, FILE *headerptr, double **outdata3c, int nt);

segy tr;

int
main(int argc, char **argv)
{

        /* OUTPUT FILE POINTERS */
    FILE *stalta=NULL;
    FILE *headerfp=NULL; /* temporary file for trace headers */
        
    char *file=NULL;         /* base of output file name(s) */
    char *fname=NULL;        /* complete output file name */
        
    int nt;			/* number of samples in one trace	*/
    int nsta;               /* number of samples for short term window */
    int nlta;               /* number of samples for long term window  */
    int verbose;		/* design info flag 			*/
    int it;                 /* index for time sample in main loop */


    double *data_trace;      /* individual trace data                */
   	double dt;		/* sample spacing, sec			*/
//  	double nyq;		/* nyquist frequency			*/
//    double sta=0;              /* short term avg. value              */
//    double lta=1.0e-99;              /* long term avg. value           */
    double trigger;          /* threshold value for detection        */
    double *charfct=NULL; /* output sta/lta trace                */

	cwp_Bool seismic;	/* is this seismic data?		*/



	
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get info from first trace */ 
	if (!gettr(&tr))  err("can't get first trace");
	seismic = ISSEISMIC(tr.trid); 
		 
	if (!seismic)
		warn("input is not seismic data, trid=%d", tr.trid);
        if (!getparint("ns", &nt))	nt = tr.ns;
	if (!getpardouble("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt) err("dt field is zero and not getparred");


	/* Get Parameters */
	if (!getparint("verbose", &verbose))	verbose = 0;
	if (!getparint("nlta", &nlta))	nlta = nt/ 10;
	if (!getparint("nsta", &nsta))	nsta = nlta / 10;
    if (!getparstring("file", &file)) file="sta_lta";	
	
    /* allocate space for input data and analysis vectors */
    data_trace = ealloc1double(nt); data_trace-=1;


    /* open temporary file for trace headers */
    headerfp = etmpfile();

    /* set filenames and open files */
    fname = malloc( strlen(file)+7 );    
    sprintf(fname, "%s.su", file);    stalta = efopen(fname, "w");	    
    free(fname);    
    
    /* allocate and zero out space for output data */
    charfct = ealloc1double(nt); 
    memset((void *) charfct, 0, nt*FSIZE);    

        /* Main loop over traces */
        

	do {
	    /* store trace header in temporary file and read data */
        efwrite(&tr, HDRBYTES, 1, headerfp);
      	memcpy((void *) data_trace, (const void *) &tr.data, nt*FSIZE);        


        /* STA/LTA trace generation */
        recstalta(data_trace, charfct, nt, nsta, nlta);           
        fputdata(stalta, headerfp, charfct, nt);


	} while (gettr(&tr));

   /* close files */
   efclose(headerfp);
   efclose(stalta);
	return(CWP_Exit());
}


/**********************************************************************/
/* Functions used internally                                          */
/**********************************************************************/


/**********************************************************************/
/* Functions for data output                                          */
/**********************************************************************/


/* write one-component data into file */

void fputdata(FILE *fileptr, FILE *headerptr, double *outdata, int nt)
{    
    efread(&tr, 1, HDRBYTES, headerptr);
    erewind(headerptr);
    
    memcpy((void *)tr.data, (const void *) outdata, nt*FSIZE);

    fputtr(fileptr, &tr);
}

/* write three-component data into file */

void fputdata3c(FILE *fileptr, FILE *headerptr, double **outdata3c, int nt)
{
    int i;
   
    for(i=1;i<=3;i++) {
        efread(&tr, 1, HDRBYTES, headerptr);
        
        memcpy((void *)tr.data, (const void *) outdata3c[i], nt*FSIZE);

        fputtr(fileptr, &tr);
    }
    erewind(headerptr);
}

/* END OF FILE */    
