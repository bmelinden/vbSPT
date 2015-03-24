/*
 HMM_multiViterbi_log.c
 =========================================================================
 
 Copyright (C) 2014 Martin Lindén, E-mail: bmelinden@gmail.com

 This program is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or any later
 version.   
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 Public License for more details.

 Additional permission under GNU GPL version 3 section 7
 
 If you modify this Program, or any covered work, by linking or combining it
 with Matlab or any Matlab toolbox, the licensors of this Program grant you 
 additional permission to convey the resulting work.
 
 You should have received a copy of the GNU General Public License along
 with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*
 Matlab function definition:
 function s=HMM_multiViterbi_log(lnQ,lnH,iEnd) 

 Finds Viterbi path (most likely hidden state sequence) in an HMM, where Q 
 is the transition matrix, not necessarily normalized, lnQ=log(Q), and 
 H(t,s) is the emission likelihood at time t, lnH=log(H). 
 
 M.L. 2011-12-20

 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

/* input Arguments */
#define	LNQ_IN  prhs[0]
#define	LNH_IN  prhs[1]
#define	IEND_IN prhs[2]

#define	S_OUT plhs[0]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    double *lnQ,*lnH,*iEnd; /* inputparameters */
    int *s;          /* output parameter */
    int T,t,Mends,Nends,n,tStart,tEnd;                       /* temporary variables */
    int ss,N,j,k;
    double Z;
    double *lnp0, *lnp1,*lnpp;
    int *MaxPrev;
    mxArray *mxMaxPrev,*mxpp, *mxp1, *mxp0; /* temporary arrays */
    
    /* check number of input/output arguments */
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required.");
    if (nrhs != 3)
        mexErrMsgTxt("Three input arguments required.");
    
    /* size of input variables */
    T = mxGetM(LNH_IN); /* number of rows = number of time points*/
    N = mxGetN(LNH_IN); /* number of columns = number of states*/
    
    /* check that the other variables have consistent sizes */
    if ( mxGetM(LNQ_IN) != N || mxGetN(LNQ_IN) != N)
        mexErrMsgTxt("Q is not an N by N matrix (N = # columns in lnH).");
    
    /* require iEnd to be a 1 x N matrix */
    Mends = mxGetM(IEND_IN); /* number of rows */
    Nends = mxGetN(IEND_IN); /* number of rows */
    /* printf("input iEnds is [ %d %d ] \n",Mends,Nends); */
    if( (Mends==1) && (Nends >=1) ){ /* then do nothing */
    }else if( (Nends == 1) && (Mends >=1) ){ /* a N x 1 is also OK */
        Nends=Mends;
        Mends=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("iEnds must be a 1-by-N or N-by-1 matrix.");
    }

    /* retrieve input data */
    lnQ=mxGetPr(LNQ_IN);
    lnH=mxGetPr(LNH_IN);
    iEnd=mxGetPr(IEND_IN);
    
    /* Create an mxArray for the output data */
    S_OUT = mxCreateNumericMatrix(T,1,mxINT32_CLASS, mxREAL);
    s= (int*) mxGetPr(S_OUT); /* pointer to output array*/
    
    /* allocate temporary arrays */
    mxp0=mxCreateDoubleMatrix(1,N, mxREAL);
    mxp1=mxCreateDoubleMatrix(1,N, mxREAL);
    lnp0=mxGetPr(mxp0);
    lnp1=mxGetPr(mxp1);
    
    mxMaxPrev=mxCreateNumericMatrix(T,N,mxINT32_CLASS,mxREAL);
    MaxPrev=(int*)mxGetPr(mxMaxPrev);
    mxpp=mxCreateDoubleMatrix(N,1, mxREAL);
    lnpp=mxGetPr(mxpp);
    
    /* start of actual algorithm */
    tStart=0;
    for(n=0;n<Nends;n++){
        tEnd=(int)(iEnd[n]);   

        /*lnP1=lnH(tStart,:); % initial probability, not normalized*/
        for(k=0;k<N;k++){
            lnp1[k]=lnH[tStart+k*T];}
        
        /*for tV=(tStart+1):tEnd  */
        for(t=tStart+1;t<tEnd;t++){
            Z=0.0;
            for(k=0;k<N;k++){/*lnP0=lnP1-mean(lnP1); lnP1=zeros(1,N);*/
                lnp0[k]=lnp1[k];
                Z=Z+lnp0[k]/N;
                lnp1[k]=0.0;
            }
            for(k=0;k<N;k++){
                lnp0[k]=lnp0[k]-Z;
            }
                    
            for(j=0;j<N;j++){
                for(k=0;k<N;k++){
                    /* lnPP(kV)=lnP0(kV)+lnQ(kV,jV)+lnH(tV,jV); */
                    lnpp[k]=lnp0[k]+lnQ[k+j*N]+lnH[t+j*T];
                }
                /* probability of previous state before ending up at jV.
                 * [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP); */
                lnp1[j]=lnpp[0];
                MaxPrev[t+j*T]=0;
                for(k=1;k<N;k++){
                    if(lnpp[k]>lnp1[j]){
                        lnp1[j]=lnpp[k];
                        MaxPrev[t+j*T]=k;
                    }
                }
            }
        }
        /* [~,S(tEnd)]=max(lnP1);*/
        s[tEnd-1]=0;
        Z=lnp1[0];
        for(k=1;k<N;k++){
            if(lnp1[k]>Z){
                Z=lnp1[k];
                s[tEnd-1]=k;
            }
        }
        /*  for tV=tEnd-1:-1:tStart*/
        /*      S(tV)=MaxPrev(tV+1,S(tV+1)); end*/
        for(t=tEnd-2;t>=tStart;t--){
            ss=s[t+1];
            s[t]=MaxPrev[t+1+T*ss];
        }
        tStart=tEnd;
    }
    
    /* finished actual algorithm */
    /* transform back to matlab's index convention */
    for(t=0;t<T;t++){
        s[t]=s[t]+1;
    }
    /* destroy temporary arrays */
    mxDestroyArray(mxp0);
    mxDestroyArray(mxp1);
    mxDestroyArray(mxpp);
    mxDestroyArray(mxMaxPrev);
}
