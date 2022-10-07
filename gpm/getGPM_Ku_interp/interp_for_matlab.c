/* interpolate data onto regular grid */

/* this is the guts of the original PRInterp.c code written by
   Stacy Brodzik for interpolating the TRMM PR data.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include "interp_for_matlab.h"
#include "mex.h"

/* Computational code */

/* void interp(int argv, char** argc)  { */
void interp_for_matlab(int size2D, int gridColumns, int gridRows, int gridAlts,
	    float RADIUS, char* interp_type, double* latitude, double* longitude,
	    double* inputArray, double* gridLatitude, double* gridLongitude, 
        double* gridNoSwath, double* outputArray)  {

  /* declaration of variables */
  int ilevel, i, j, index;
  int count = 0, count2 = 0;
  double *latitude_2d,  *longitude_2d, *inputArray_2d;
  double *weighting, *differences, *products, *distance, *filteredDistance;
  double interpolatedvalue, total, neighbour, proximity;
  int *indices, *indices2;
	
  /* allocate space for arrays */
  latitude_2d = (double *)malloc(size2D*sizeof(double));
  longitude_2d = (double *)malloc(size2D*sizeof(double));
  inputArray_2d = (double *)malloc(size2D*sizeof(double));
  if (!latitude_2d || !longitude_2d || !inputArray_2d)  {
    fprintf(stdout,"Memory allocation error for one of 2d arrays.\n");
    exit(1);
  }

  /* check input params */
  fprintf(stdout,"TOP OF INTERP: gridColumns=%d gridRows=%d gridAlts=%d interp_type=%s\n",
	  gridColumns, gridRows, gridAlts, interp_type);

  /* interpolate each level separately */
  for (ilevel=0; ilevel<gridAlts; ilevel++)  { 
    if (ilevel%10 == 0){ 
      fprintf(stdout,"ilevel = %d\n", ilevel);
    }

    /* grab inputArray, latitude, and longitude for this level */
    for(j=0;j<size2D;j++)  {
      *(latitude_2d + j) = *(latitude + (ilevel*size2D) + j);
      *(longitude_2d + j) = *(longitude + (ilevel*size2D) + j);
      *(inputArray_2d + j) = *(inputArray + (ilevel*size2D) + j);
    }

    for (i=0; i<gridColumns*gridRows; i++)  {
      if (gridNoSwath[i]==1.0)  {
        outputArray[(ilevel*gridRows*gridColumns + i)] = NODATA;
      } 
      else {

        /* allocate space for working arrays */
        weighting = (double *) malloc(size2D*sizeof(double)); 
        indices = (int *) malloc(size2D*sizeof(double));
        indices2 = (int *) malloc(size2D*sizeof(double));
        differences = (double *) malloc(size2D*sizeof(double)); 
        products = (double *) malloc(size2D*sizeof(double));
        distance = (double *) malloc(size2D*sizeof(double));
        filteredDistance = (double *) malloc(size2D*sizeof(double));
        if (!weighting || !indices || !indices2  || !differences || 
                !products || !distance)  {
            fprintf(stdout,"Memory allocation error for one of working arrays.\n");
            exit(1);
        } 

        /* initialize indices arrays */
        for(j=0;j<size2D;j++) {
            indices2[j] = NODATA;
            indices[j]  = NODATA;
        }

        /* find points within latitude/longitude box */
        searcher2(latitude_2d, longitude_2d, gridLatitude, gridLongitude, 
            inputArray_2d, indices2, &count2, size2D, i);

        for(j=0; j<count2; j++)  {
            distance[indices2[j]] = 
                degrees2distance(gridLatitude[i], latitude_2d[indices2[j]], 
                    gridLongitude[i], longitude_2d[indices2[j]]);
        }

        /* find points within RADIUS */
        searcher(distance, indices, indices2, &count, &count2, RADIUS);

        if (!count)  {
            /* no points to interpolate */
            outputArray[(ilevel*gridRows*gridColumns + i)] = NODATA;
        } 

        else  { 

            /* find distance to each point within radius */
            for(j=0; j<count; j++)  {
                filteredDistance[j] = distance[indices[j]];
            }

            /* find nearest neighbor */
            proximity = RADIUS;
            for(j=0;j<count;j++)  {
                if (filteredDistance[j] < proximity)  {
                    proximity = filteredDistance[j];
                    index = j;
                }
            }
            neighbour = inputArray_2d[indices[index]];
            total = neighbour;

            if ( !strcmp(interp_type,"weighted") ) { 
            
                /* find weighting factors */
                for(j=0;j<count;j++)  {
                    weighting[j] = cressman(filteredDistance[j],RADIUS);
                }  

                /* compute weighted average */
                for(j=0;j<count;j++)  {
                    differences[indices[j]] = neighbour - inputArray_2d[indices[j]];
                    products[j] = weighting[j]*differences[indices[j]];
                    interpolatedvalue = products[j] + total;
                    total = interpolatedvalue;
                    /*
                    interpolatedvalue = 10 * log10(interpolatedvalue);
                    if (interpolatedvalue < MINDBZ) {
                        interpolatedvalue = NODATA;
                    }
                    */ 
                    outputArray[(ilevel*gridRows*gridColumns + i)] = interpolatedvalue;
                }
                
            } else if ( !strcmp(interp_type,"nn") ) {
                outputArray[(ilevel*gridRows*gridColumns + i)] = neighbour;
            } else {
                fprintf(stdout,"interp_type must be either weighted or nn.  Exiting ...\n");
                exit(-1);
            }

        } /* end else (!count) */

        /* deallocate memory for working arrays */
        free(products);
        free(differences); 
        free(weighting);
        free(indices);
        free(indices2);
        free(distance);
        free(filteredDistance);

      } /* end else (allocate space for working arrays) */

    } /* for each gridpoint */

  } /* end for level... */

  /* deallocate memory */
  free(latitude_2d);
  free(longitude_2d);
  free(inputArray_2d);

  /*****************************************************************************/
  /*****************************************************************************/

  /* For all non-NODATA values:
     1. convert linear values to log values and 
     2. filter out reflectivity values < MINDBZ
     3. Print all values to output file
  */

  /*
  for(j=0;j<gridAlts*gridRows*gridColumns;j++)  {
    if (interpolation[j] != NODATA)  {
      interpolation[j] = 10 * log10(interpolation[j]);
      if (interpolation[j] < MINDBZ )  {
	interpolation[j] = NODATA;
      }
    }
    outputArray[j] = (float)interpolation[j];
    fprintf(fpGridReflectivity, "%lf ", interpolation[j]);
  }
  fprintf(fpGridReflectivity, "\n");
  */

}

/*******************************************************************************/

double degrees2distance(double gridLatitude, double latitude, 
			double gridLongitude, double longitude)  {
  double difflong = (DEG2RAD*longitude) - (DEG2RAD*gridLongitude);
  double absdifflong = fabs(difflong);

  return (111.24 * RAD2DEG * 
	  (acos(((sin(DEG2RAD*gridLatitude) * sin(DEG2RAD*latitude)) + 
		 (cos(DEG2RAD*gridLatitude) * cos(DEG2RAD*latitude) * 
		  cos(absdifflong))))));
}

/*******************************************************************************/

double cressman(double d2gp,float RADIUS)  {
  /* calculate weighting factor */
  return ((-1)*(((RADIUS*RADIUS)-(d2gp*d2gp))/((RADIUS*RADIUS)+(d2gp*d2gp))));
}

/*******************************************************************************/

void searcher(double *distance, int indices[], int indices2[], 
	      int *count, int *count2, float RADIUS)  {
  int i;
  int subi;

  *count = 0;
  subi = 0;

  for(i=0;i<*count2;i++)  {
    
    if (((*(distance + indices2[i])) < RADIUS) && 
	((*(distance + indices2[i])) >= 0.0))  {
      indices[subi] = (*(indices2 + i));
      (*count)++;
      subi++;
    }
  }
}

/*******************************************************************************/

void searcher2(double *latitude_2d, double *longitude_2d, 
	       double *gridLatitude, double *gridLongitude, 
	       double *inputArray_2d, int indices2[], 
	       int *count2, int size2D, int i)  {

  int j;
  *count2 = 0;

  double temp1=(*(gridLatitude + i));
  double temp2=(*(gridLongitude + i));

  for(j=0;j<size2D;j++)  {
   		
    if (fabs(*(latitude_2d + j))<(fabs(temp1) + DEGREES) && 
	fabs(*(latitude_2d + j))>(fabs(temp1) - DEGREES) && 
	fabs(*(longitude_2d + j)<(fabs(temp2) + DEGREES)) && 
	fabs(*(longitude_2d + j))>(fabs(temp2) - DEGREES))  {
	 
      /* if (inputArray_2d[j]>(-666))  { WHY WAS THIS HARDCODED IN? */
      if (inputArray_2d[j] != NODATA)  {
	indices2[*count2] = j;
	(*count2)++;
      }						
    }
  }
}

/*******************************************************************************/

/* Gateway function */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize size2D;                  /* size of 2D array */
    mwSize gridColumns;             /* number of grid lons */
    mwSize gridRows;                /* number of grid lats */
    mwSize gridAlts;                /* number of grid alts */
    double RADIUS;                  /* interpolation radius in km */
    char *interp_type;              /* "weighted" or "nn" */
    double *latitude;               /* 1X(size2D*LEVELS) input latitude array */
    double *longitude;              /* 1X(size2D*LEVELS) input longitude array */
    double *inputArray;             /* 1X(size2D*LEVELS) input array */
    double *gridLatitude;           /* 1X(gridRows*gridColumns*gridAlts) array */
    double *gridLongitude;          /* 1X(gridRows*gridColumns*gridAlts) array */
    double *gridNoSwath;            /* 1X(gridRows*gridColumns*gridAlts) array */
    double *outputArray;            /* output gridded matrix */

    /* check for proper number of arguments */
    if(nrhs!=12) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:nrhs","Twelve inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notScalar","Input size2D must be a scalar.");
    }
    
    /* make sure the second input argument is scalar */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notScalar","Input gridColumns must be a scalar.");
    }
    
    /* make sure the third input argument is scalar */
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notScalar","Input gridRpws must be a scalar.");
    }
    
    /* make sure the fourth input argument is scalar */
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3]) ||
         mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notScalar","Input gridAlts must be a scalar.");
    }
    
    /* make sure the fifth input argument is scalar */
    if( !mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4]) ||
         mxGetNumberOfElements(prhs[4])!=1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notScalar","Input RADIUS must be a scalar.");
    }
    
    /* make sure the sixth input argument is string */
    if( mxIsChar(prhs[5]) != 1 ) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notString","Input interp_type must be a string.");
    }
    
    /* make sure the seventh input argument is type double */
    if( !mxIsDouble(prhs[6]) || 
         mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input latitude must be type double.");
    }
    
    /* make sure the eighth input argument is type double */
    if( !mxIsDouble(prhs[7]) || 
         mxIsComplex(prhs[7])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input longitude must be type double.");
    }
    
    /* make sure the ninth input argument is type double */
    if( !mxIsDouble(prhs[8]) || 
         mxIsComplex(prhs[8])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input inputArray must be type double.");
    }
    
    /* make sure the tenth input argument is type double */
    if( !mxIsDouble(prhs[9]) || 
         mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input gridLatitude must be type double.");
    }
    
    /* make sure the eleventh input argument is type double */
    if( !mxIsDouble(prhs[10]) || 
         mxIsComplex(prhs[10])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input gridLongitude must be type double.");
    }
    
    /* make sure the twelvth input argument is type double */
    if( !mxIsDouble(prhs[11]) || 
         mxIsComplex(prhs[11])) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notDouble","Input gridNoSwath must be type double.");
    }
    
    /* check that number of rows in seventh thru twelvth input arguments is 1 */
    if(mxGetM(prhs[6])!=1 || mxGetM(prhs[7])!=1 || mxGetM(prhs[8])!=1 || 
            mxGetM(prhs[9])!=1 || mxGetM(prhs[10])!=1 || mxGetM(prhs[11])!=1) {
        mexErrMsgIdAndTxt("getGPM:interp_for_matlab:notRowVector","Input arrays must be a row vectors.");
    }
    
    /* get the value of the first 5 scalar inputs  */
    size2D = mxGetScalar(prhs[0]);
    gridColumns = mxGetScalar(prhs[1]);
    gridRows = mxGetScalar(prhs[2]);
    gridAlts = mxGetScalar(prhs[3]);
    RADIUS = mxGetScalar(prhs[4]);

    /* get the value of the input string */
    interp_type = mxArrayToString(prhs[5]);
    if(interp_type == NULL) 
      mexErrMsgIdAndTxt( "getGPM:interp_for_matlab:conversionFailed",
              "Could not convert input to string.");
    
    /* create a pointer to the real data in the input matrices  */
    latitude = mxGetPr(prhs[6]);
    longitude = mxGetPr(prhs[7]);
    inputArray = mxGetPr(prhs[8]);
    gridLatitude = mxGetPr(prhs[9]);
    gridLongitude = mxGetPr(prhs[10]);
    gridNoSwath = mxGetPr(prhs[11]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)(gridColumns*gridRows*gridAlts),mxREAL);

    /* get a pointer to the real data in the output matrix */
    outputArray = mxGetPr(plhs[0]);

    /* call the computational routine */
    interp_for_matlab(size2D,gridColumns,gridRows,gridAlts,RADIUS,
            interp_type,latitude,longitude,inputArray,gridLatitude,
            gridLongitude,gridNoSwath,outputArray);
}
