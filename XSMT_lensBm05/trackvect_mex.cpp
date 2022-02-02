#include "mex.h"
#include "matrix.h"
#include <math.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, oi, oj, k;/* some loop variable*/

	int imw, imh, imd;/* size of the input sample image matrix*/
    int tpw, tph, tpd;/* size of the input reference image matrix*/
    int outw, outh, outd; /* real size of output (taking into account 'shape')*/
	
    int resw, resh;/* working resolution for future programming*/

	int rgw, rgh;/* how far to correlate*/

	double *im, *tp;
    double *mat;
    mxArray *matc;
    double *res;
    mwSize nDims;
    mwSize const *pSize;
    mwSize out_dims[3];
    double corrf = 0;
    
    int maxI, maxJ;
    double f0, f1, f2, snr;
    double xDisp, yDisp;
    
    int _intLengthX, _intLengthY;

    
    //check stuffs
    if( nrhs < 2 ) {
		mexErrMsgTxt("result = trackvecMEX(SAMPLE_STACK, REF_STACK,range)");
	} else if( !( mxIsDouble(prhs[0])&&mxIsDouble(prhs[1])) ) {
		mexErrMsgTxt("STACKS must be of type double");
	} else if( nrhs > 2 && !mxIsClass(prhs[2], "int8" )) {
		mexErrMsgTxt("RANGE parameter must be a scalar value");
	}
    
    
    if( nrhs>2 ) {
        rgw = mxGetScalar( prhs[2]);
        rgh = mxGetScalar( prhs[2]);
		if( rgw > 41){
			rgw = 41;    rgh = 41;}
		else if( rgw < 2){
			rgw = 2;    rgh = 2;}
     }else {
			rgw = 5;    rgh = 5;
     }
    
    
    im = mxGetPr(prhs[0]);
    nDims = mxGetNumberOfDimensions(prhs[0]);
    if( nDims != 3) {mexErrMsgTxt("Stack must be of dimension 3");}
    pSize = mxGetDimensions(prhs[0]);
    
	imw = pSize[0] ;
	imh = pSize[1] ;
    imd = pSize[2] ;

	tp = mxGetPr(prhs[1]);
	nDims = mxGetNumberOfDimensions(prhs[1]);
	if( nDims != 3) {mexErrMsgTxt("Stack must be of dimension 3");}
    pSize = mxGetDimensions(prhs[1]);
    
    tpw = pSize[0];
    tph = pSize[1];
    tpd = pSize[2];
    
    
    if( tpw != imw || tph != imh || tpd != imd ) {
		mexErrMsgTxt("Stacks must have the same size");
	} 
	// size of the ouput image will be reduced by the size of the interrogation window (no info on the edge)
    outw = imw -2*rgw-1;
    outh = imh -2*rgh-1;
    
    
        
    /* Allocate memory for the new_dims array on the fly 
    out_dims = mxMalloc(nDims * sizeof(*out_dims));*/
    out_dims[0] = outw; out_dims[1] = outh; out_dims[2] = 3;
    plhs[0] = mxCreateNumericArray( nDims, out_dims, mxDOUBLE_CLASS,mxREAL );
    
    
    res = mxGetPr(plhs[0]);// the output 3D array of size [imw -2*rgw-1,imw -2*rgw-1, 5]
    matc = mxCreateDoubleMatrix(2*rgw+1,2*rgh+1,mxREAL);
	mat = mxGetPr(matc);// a matrix which the pixel correlation map
    
    // printf("so fat so good %d",tpw);
    // oi,oj are the pixel where we calculate the peak position
	// we use i, j for the neighbouring pixels	
	// k for the third dimension of the stack
    for (oi = rgw+1; oi < (imw -rgw); oi++)
    {    
        for (oj = rgh+1; oj < (imh -rgh); oj++)
        {
            //calc correlation coef of the sample vectors with the interrogation window ref vectors
            for (i = -rgw; i < rgw+1; i++)
            {
                
                for (j = -rgh; j < rgh+1; j++)
                {
                    corrf = 0;
                    for (k = 1; k < imd+1; k++)
                    {
                       corrf = corrf +  im[(oi-1) + (oj-1)*imw + (k-1)*imw *imh] * tp[(oi +i-1) + (oj +j-1)*imw + (k-1)*imw *imh];            
                    }
                    
                    mat[(i+rgw) + (j+rgh)*(2*rgw+1)] = corrf/imd;
                }
            }         
            /* we have now a matrix with a peak of correlation for the pixel oi,oj*/
            

            /*========================pixel accuracy==============================*/    

            _intLengthX = 2*rgh+1;
            _intLengthY = 2*rgw+1;


            /* First, find the maximum value in the correlation map*/
            double max = -1000.0;
            for (i = 1; i < _intLengthX-1; i++)
            {
                for (j = 1; j < _intLengthY-1; j++)
                {
                    if (mat[_intLengthX*i + j] > max)
                    {
                        max = mat[_intLengthX*i + j];
                        maxI = i;
                        maxJ = j;
                    }
                }
            }

            /*========================subpixel accuracy==============================*/
            				
			
			
        /*     printf("wdith %d\n", _intLengthX);*/
        /*     printf("length %d\n",_intLengthY);*/
        /*     printf("maxI %d\n",maxI);*/
        /*     printf("maxJ %d\n",maxJ);*/
            /* Need to calculate the average to which to compare the signal*/
            double mean = 0.0;
            int count = 0;
            for (i = 0; i < _intLengthX; i++)
            {
                for (j = 0; j < _intLengthY; j++)
                {
                    if (!(i < (maxI+1) && i > (maxI-1) && j < (maxJ+1) && j < (maxJ-1)))
                    {
                        mean += fabs(mat[_intLengthX*i + j]);
                        count++;
                    }
                }
            }
			
			// 
            mean = mean / double(count);
            /* Assign the signal-to-noise (snr) value*/
            snr = double(max / mean);
			
			
// 			/*======================method 1 which fit 2 1D gaussian peaks===============*/
// 			
// 			
//             /* Gaussian fitting to the peak in the correlation map*/
//             f0 = log(fabs(mat[_intLengthX*maxI + maxJ]));
// 
//             f1 = log(fabs(mat[_intLengthX*(maxI-1) + maxJ]));
//             f2 = log(fabs(mat[_intLengthX*(maxI+1) + maxJ]));
//             yDisp = double(_intLengthY/2) - (double(maxI+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
//             /* Assign the vertical displacment*/
// 
// 
//             f1 = log(fabs(mat[_intLengthX*maxI + maxJ - 1]));
//             f2 = log(fabs(mat[_intLengthX*maxI + maxJ + 1]));
//             xDisp = double(_intLengthX/2) - (double(maxJ+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
//              /* Assign the horizontal displacement*/
// 
//         /*    plhs[1] = mxCreateDoubleMatrix( 5, 1, mxREAL );
//              pic = mxGetPr(plhs[1]);
//              printf("xDisp %f\n",xDisp);
//              printf("yDisp %f\n",yDisp);
//              printf("snr %f\n",snr);*/

            /*===================method 2: 1 2d peak==================*/

            double dx = (mat[_intLengthX*(maxI+1) + maxJ] - mat[_intLengthX*(maxI-1) + maxJ]) / 2.f;
            double dy = (mat[_intLengthX*maxI + maxJ + 1] - mat[_intLengthX*maxI + maxJ - 1]) / 2.f;
            double dxx = (mat[_intLengthX*(maxI+1) + maxJ] + mat[_intLengthX*(maxI-1) + maxJ] - 2 * mat[_intLengthX*maxI + maxJ]);
            double dyy = (mat[_intLengthX*maxI + maxJ + 1] + mat[_intLengthX*maxI + maxJ - 1] - 2 * mat[_intLengthX*maxI + maxJ]);
            double dxy = (mat[_intLengthX*(maxI + 1) + maxJ + 1] - mat[_intLengthX*(maxI + 1) + maxJ - 1] - mat[_intLengthX*(maxI - 1) + maxJ + 1] + mat[_intLengthX*(maxI - 1) + maxJ - 1]) / 4.f;

            double det = 1.f/(dxx*dyy - dxy*dxy);

            double ix = - (dyy*dx - dxy*dy) * det + maxI + 1 - (_intLengthY/2) ;  
            double iy = - (dxx*dy - dxy*dx) * det + maxJ + 1 - (_intLengthX/2) ;   

            /*     printf("xDisp %f\n",ix );*/
            /*     printf("yDisp %f\n",iy );*/

            res[oi-rgw-1 + (oj-rgh-1)*outw ] = iy;//-xDisp;
            res[oi-rgw-1 + (oj-rgh-1)*outw + outw * outh] = ix;//-yDisp;
            res[oi-rgw-1 + (oj-rgh-1)*outw + 2*outw *outh] = max;//iy;
            //res[oi-rgw-1 + (oj-rgh-1)*outw + 3*outw *outh] = ix;
            //res[oi-rgw-1 + (oj-rgh-1)*outw + 4*outw *outh] = max;//snr;
        
             
       }
    }

	
	
};