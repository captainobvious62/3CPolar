/* Modified from SUPOLAR, Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                                                */

/* SUPOLAR_PS: $Revision: 1.9 $ ; $Date: Sat, 23 Jul 2016 00:10:26 -0700 */

#include <omp.h>
#include <stdio.h>

#include "su.h"
#include "segy.h"
#include "header.h"


/*********************** self documentation *****************************/
char *sdoc[] = {
"                                                                       ",
" SUMATRIX   - Sandbox to produce and play with SU data in              ",
"            - matrix format                                            ",
"            - to provide P and S arrivals                              ",
"            - Modification of SUPOLAR                                  ",
"                                                                       ",
" sumatrix <stdin ntr= [optional parameters]                               ",
"                                                                       ",
" Required parameters:                                                  ",
"    -ntr (number of traces) needs to be specified in the SU file       ",
"                                                                       ",
" Optional parameters:                                                  ",
"    dt=(from header)  time sampling interval in seconds                ",
"    wl=0.1            correlation window length in seconds             ",
"    win=boxcar        correlation window shape, choose \"boxcar\",     ",
"                      \"hanning\", \"bartlett\", or \"welsh\"          ",
"    file=polar        base of output file name(s)                      ",
"    rl=1              1 = rectilinearity evaluating 2 eigenvalues,     ",
"                      2, 3 = rectilinearity evaluating 3 eigenvalues   ",
"    rlq=1.0           contrast parameter for rectilinearity            ",
"    dir=1             1 = 3 components of direction of polarization    ",
"                          (the only three-component output file)       ",
NULL};

/* 
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 */
/**************** end self doc *******************************************/

/* prototypes of functions used internally */
void calc_window(float *w, int iwl, int iwin);
void fputdata(FILE *fileptr, FILE *headerptr, float outdata[], int nt);
void fputdata3c(FILE *fileptr, FILE *headerptr, float **outdata3c, int nt);


/* window shape identifiers */
#define WBOXCAR 0
#define WBARTLETT 1
#define WHANNING 2
#define WWELSH 3


segy tr;                /* SEG-Y record (trace and header) */

int main(int argc, char **argv)
{
    /* output file pointers */    
    FILE *fileout=NULL;
    FILE *headerfp=NULL;

    char *file=NULL;         /* base of output file name(s) */
    char *fname=NULL;        /* complete output file name */
    char *win=NULL;          /* shape of used time window */        
    /* flags (see selfdoc) */
    int verbose;    
    int i,j;      /* indices for components (in loops) */
  //  int it;             /* index for time sample in main loop */
    int iwl;            /* correlation window length in samples */
    int nt;             /* number of time samples in one trace */
    int iwin=0;           /* time window shape identifier */
    int ntr;               /* number of traces */
//    int kwl;            /* kurtosis window length in seconds */    


    float **data_matrix;     /* three-component data ([1..3][0..nt-1]) */
    float dt;           /* sampling interval in seconds */    
    float wl;           /* correlation window length in seconds */
    float *w;

    /* initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* get number of traces */
    if (!getparint("ntr", &ntr)) err("Please specify number of traces");
    
    /* get info from first trace */
    if(!gettr(&tr)) err("can't get first trace");
    nt = tr.ns;
    warn("%d",tr.ns);
    /* get parameters ... */
    if (!getparstring("file", &file)) file="polar";
    if (!getparfloat("wl", &wl)) wl = 0.1;
    if (!getparstring("win", &win)) win="boxcar";
    if (!getparfloat("dt", &dt)) dt = ((double) tr.dt)/1000000.0;
    if (!getparint("verbose", &verbose)) verbose = 0;
//    if (!getparint("kwl", &kwl)) kwl = 5 * ((int) 1/dt);    

    checkpars();


    /* get time window shape */
    if      (STREQ(win, "boxcar"))   iwin=WBOXCAR;
    else if (STREQ(win, "bartlett")) iwin=WBARTLETT;
    else if (STREQ(win, "hanning"))  iwin=WHANNING;
    else if (STREQ(win, "welsh"))    iwin=WWELSH;
    else err("unknown win=%s", win);

     /* convert seconds to samples */
    if (!dt) {
        dt = 0.004;
        warn("dt not set, assuming dt=0.004");
    }

    iwl = NINT(wl/dt);
        
    /* data validation */
    if (iwl<1) err("wl=%g must be positive", wl);
    if (iwl>nt) err("wl=%g too long for trace", wl);


    /* open temporary file for trace headers */
    headerfp = etmpfile();

    /* set filenames and open files */
    fname = malloc( strlen(file)+7 );
    sprintf(fname, "%s.txt", file); 
    fileout = efopen(fname, "w");
  //  sprintf(fname, "%s.pkur", file); if (rl && theta) pkur = efopen(fname, "w");
  //  sprintf(fname, "%s.skur", file); if (rl && theta) skur = efopen(fname, "w");
    free(fname);

    /* allocate space for input data and analysis matrices */
    /* index ranges used here: data3c[1..3][0..nt-1],      */
    /* a[1..3][1..3], v[1..3][1..3], d[1..3]               */

    data_matrix = ealloc2float(nt,ntr); data_matrix-=1;
    warn("Trace Start Time: %d %d %d %d %d", tr.year, tr.day, tr.hour, tr.minute, tr.sec);    



/* ************************ BEGIN CALCULATION ******************************* */    

    /* loop over traces */
    // Need to convert this do while loop into a for loop so as to be easier
    // to parallelize

    for (i=0; i<ntr; i++)
    {
        efwrite(&tr, HDRBYTES, 1, headerfp);
        memcpy((void *)data_matrix[i], (const void *) tr.data, nt*FSIZE);
        /* store trace header in temporary file and read data */
       // efwrite(&tr, HDRBYTES, 1, headerfp);
        for (j=0; j<nt; j++)
        {   
            fprintf(fileout, "%f", data_matrix[i][j]);
        }

    }
/* *************************** END CALCULATION ****************************** */
   /* end loop over traces */
    



/* ***************************** BEGIN WRITE ******************************** */


/* ****************************** END WRITE ********************************* */    
 
    return(CWP_Exit());
}





/**********************************************************************/
/* Functions used internally                                          */
/**********************************************************************/

/* calculate time window weights for a smooth time window,        */
/* after Sheriff, 1991, p. 335; tapered windows are suggested by  */
/* Jurkevics, 1988                                                */

void calc_window(float *w, int iwl, int iwin)
{
    int i,j;      /* loop indices within window */
    float m;      /* half of time window (iwl/2) */
    
    m = (float) iwl/2;
    
    switch (iwin) {
        case WBOXCAR :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 1.0;
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;
        case WBARTLETT :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 1.0 - (fabs((float) i) / m);
                
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;
        case WHANNING :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 0.5 + 0.5*cos(PI*fabs((float) i) / m);
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;            
        case WWELSH :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = fabs( fabs((float) i) / m - 1.0);
                w[j] *= w[j];            
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;            
    }
}

/* Note: The SIGN function is introduced to resolve  */
/* the 180 deg ambiguity by taking the positive      */
/* vertical component of v[1][1] (Jurkevics, 1988).  */
#define VSIGN ( (v[1][1]<0.0) ? -1.0 : 1.0 )

float calc_phi(float **v, int opt)
{
    float phi;
    
    switch (opt) {
        default:
            /* case 1: */
            /* definition after Kanasewich, 1981 */
            /* interval -pi/2 <= phi <= pi/2 */ 
            if (v[1][2]) {
                phi = atan( v[3][1] / v[2][1] );
            }
            else {
                phi = (v[3][1]>0.0) ? 0.5*PI : -0.5*PI;
            }
            break;
        case 2:
        case 3:
            /* definitions after Jurkevics, 1988 */
            /* interval -pi <= phi <= pi */
            if (v[2][1]) {
                phi = atan2( v[3][1]*VSIGN, v[2][1]*VSIGN);
            }
            else {
                phi = (v[3][1]>0.0) ? 0.5*PI*VSIGN : -0.5*PI*VSIGN;
            }
            
            /* interval 0.0 <= phi <= 2*pi */
            if (phi<0.0 && opt==3) phi += 2.0*PI;
            break;
    }   
    return phi;
}
#undef VSIGN

/**********************************************************************/
/* Functions for data output                                          */
/**********************************************************************/


/* write one-component data into file */

void fputdata(FILE *fileptr, FILE *headerptr, float *outdata, int nt)
{    
    efread(&tr, 1, HDRBYTES, headerptr);
    erewind(headerptr);
    
    memcpy((void *)tr.data, (const void *) outdata, nt*FSIZE);

    fputtr(fileptr, &tr);
}

/* write three-component data into file */

void fputdata3c(FILE *fileptr, FILE *headerptr, float **outdata3c, int nt)
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
