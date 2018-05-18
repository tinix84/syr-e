#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fsolver.h"
#include "femmcomplex.h"
#include "spars.h"
#include "mesh.h"

#include "mex.h"

/*
 * mexfsolver.cpp
 *
 * gateway to fsolver program
 *
 * This is a MEX-file for MATLAB.
 *
 */

// extern void _main();
void voidmexPrintF(const char *message);

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    FSolver SolveObj;
    std::string FilePath;
    char *inputbuf;
    mwSize buflen;
    int status;

    //(void) plhs;    /* unused parameters */

    /* Check for proper number of input and output arguments */
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    }

    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    /* Check for proper input type */
    if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )  {
        mexErrMsgTxt("Input argument must be a string.");
    }

    /* Find out how long the input string is.  Allocate enough memory
       to hold the converted string.  NOTE: MATLAB stores characters
       as 2 byte unicode ( 16 bit ASCII) on machines with multi-byte
       character sets.  You should use mxChar to ensure enough space
       is allocated to hold the string */

    buflen = mxGetN(prhs[0])*sizeof(mxChar)+1;
    inputbuf = (char *) mxMalloc(buflen);

    /* Copy the string data into buf. */
    status = mxGetString(prhs[0], inputbuf, buflen);
    //mexPrintf("Solved file: %s\n",inputbuf);

    // Tell FSolver the location of the mesh and fem files this should be
    // an extension free base name
    SolveObj.PathName = inputbuf;

    // free the input buffer
    mxFree(inputbuf);

    // make the warning message function point to mexPrintf
    SolveObj.WarnMessage = &voidmexPrintF;

    status = SolveObj.LoadProblemFile ();
	if (status != TRUE){
		mexPrintf("problem loading .fem file\n");
		plhs[0] = mxCreateDoubleScalar (1.0);
        return;
	}

	// load mesh
	status = SolveObj.LoadMesh();
	if (status != 0){

	    mexPrintf("problem loading mesh\n");
	    plhs[0] = mxCreateDoubleScalar(status);

	    switch (status)
	    {
	        case ( BADEDGEFILE ):

                mexErrMsgTxt("Could not open .edge file.\n");
                break;

	        case ( BADELEMENTFILE ):
                mexErrMsgTxt("Could not open .ele file.\n");
                break;

	        case( BADFEMFILE ):
                mexErrMsgTxt("Could not open .fem file.\n");
                break;

	        case( BADNODEFILE ):
                mexErrMsgTxt("Could not open .node file.\n");
                break;

	        case( BADPBCFILE ):
                mexErrMsgTxt("Could not open .pbc file.\n");
                break;

	        case( MISSINGMATPROPS ):
                mexErrMsgTxt("Material properties have not been defined for all regions.\n");
                break;

	        default:
                mexErrMsgTxt("AN unknown error occured.\n");
                break;
	    }

        return;
		//return -1;
	}

	// renumber using Cuthill-McKee
    //mexPrintf("renumbering nodes");
	if (SolveObj.Cuthill() != TRUE){
		mexPrintf("problem renumbering node points\n");
		plhs[0] = mxCreateDoubleScalar(3.0);
        return;
	}

    //mexPrintf("solving...");

	//mexPrintf("Problem Statistics:\n%i nodes\n%i elements\nPrecision: %3.2e\n",
	//		SolveObj.NumNodes,SolveObj.NumEls,SolveObj.Precision);
//
//     mexPrintf(outstr);

	//double mr = (8.*((double) SolveObj.NumNodes)*((double) SolveObj.BandWidth)) / 1.e06;

	if(SolveObj.Frequency == 0){

		CBigLinProb L;

		L.Precision = SolveObj.Precision;

		// initialize the problem, allocating the space required to solve it.
		if (L.Create(SolveObj.NumNodes, SolveObj.BandWidth) == FALSE)
		{
            mexPrintf("couldn't allocate enough space for matrices\n");
			plhs[0] = mxCreateDoubleScalar(4.0);
            return;
		}

		// Create element matrices and solve the problem;
		if (SolveObj.ProblemType == FALSE)
		{
			if (SolveObj.Static2D(L) == FALSE)
			{
				mexPrintf("Couldn't solve the problem\n");
				plhs[0] = mxCreateDoubleScalar(5.0);
                return;
			}
			//mexPrintf("Static 2-D problem solved\n");
		}
		else{

			if (SolveObj.StaticAxisymmetric(L) == FALSE)
			{
				mexPrintf("Couldn't solve the problem\n");
				plhs[0] = mxCreateDoubleScalar(5.0);
                return;
			}
            mexPrintf("Static axisymmetric problem solved\n");
		}

		if (SolveObj.WriteStatic2D(L) == FALSE)
		{
			mexPrintf("couldn't write results to disk\n");
			plhs[0] = mxCreateDoubleScalar(6.0);
            return;
		}
        //mexPrintf("results written to disk\n");
	}
	else
	{
		CBigComplexLinProb L;

		L.Precision = SolveObj.Precision;

		// initialize the problem, allocating the space required to solve it.

		if (L.Create(SolveObj.NumNodes+SolveObj.NumCircProps, SolveObj.BandWidth, SolveObj.NumNodes) == FALSE)
		{
            mexPrintf("couldn't allocate enough space for matrices\n");
			plhs[0] = mxCreateDoubleScalar(4.0);
            return;
		}

		// Create element matrices and solve the problem;
		if (SolveObj.ProblemType == FALSE)
		{
			if (SolveObj.Harmonic2D(L) == FALSE)
			{
				mexPrintf("Couldn't solve the problem\n");
				plhs[0] = mxCreateDoubleScalar(5.0);
                return;
			}
			mexPrintf("Harmonic 2-D problem solved\n");
		}
		else
		{
            if (SolveObj.HarmonicAxisymmetric(L) == FALSE)
            {
                mexPrintf("Couldn't solve the problem\n");
                plhs[0] = mxCreateDoubleScalar(5.0);
                return;
            }
            mexPrintf("Harmonic axisymmetric problem solved\n");
		}


		if (SolveObj.WriteHarmonic2D(L)==FALSE)
		{
			mexPrintf("couldn't write results to disk\n");
			plhs[0] = mxCreateDoubleScalar(6.0);
            return;
		}
        //mexPrintf("results written to disk.\n");
	}

	//plhs[0] = mxCreateDoubleScalar(0.0);
    plhs[0] = NULL;
}


void voidmexPrintF(const char *message)
{
    mexPrintf(message);
}
