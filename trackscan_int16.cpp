#include "mex.h"
#include "matrix.h"
#include <math.h>

// INPUT 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, oi, oj, k;/* some loop variable*/

	int imw, imh, imd;/* size of the input sample image matrix*/
    int tpw, tph, tpd, tpl1, tpl2;/* size of the input reference image matrix*/
    int outw, outh, outd; /* real size of output (taking into account 'shape')*/
	
    int resw, resh;/* working resolution for future programming*/

	int rgw, rgh;/* how far to correlate*/

	unsigned short *tp;
    float *im, *mat;
    mxArray *matc;
    float *res;
    mwSize nDims;
    mwSize const *pSize;
    mwSize out_dims[3];
    float corrf = 0;
    
    int maxI, maxJ;
    float f0, f1, f2, snr, sum1, sum2, std1, std2;
    float xDisp, yDisp;
    
    int _intLengthX, _intLengthY;

    
    //check stuffs
    if( nrhs < 1 ) {
		mexErrMsgTxt("result = trackvecMEX(SAMPLE_STACK, REF_STACK)");
	} else if( !( mxIsNumeric(prhs[0])&&mxIsNumeric(prhs[1])) ) {
		mexErrMsgTxt("STACKS must be of type uin1t16");
	} 
    
    im = (float*)mxGetPr(prhs[0]);
    nDims = mxGetNumberOfDimensions(prhs[0]);
    if( nDims != 3) {mexErrMsgTxt("Stack must be of dimension 3");}
    pSize = mxGetDimensions(prhs[0]);
    
	imw = pSize[0] ;
	imh = pSize[1] ;
    imd = pSize[2] ;

	tp = (unsigned short *)mxGetPr(prhs[1]);
	nDims = mxGetNumberOfDimensions(prhs[1]);
	if( nDims != 5) {mexErrMsgTxt("Stack must be of dimension 5");}
    pSize = mxGetDimensions(prhs[1]);
    
    tpw = pSize[0];
    tph = pSize[1];
    tpd = pSize[4];
    tpl1 = pSize[2];
    tpl2 = pSize[3];
    
    if( tpw != imw || tph != imh || tpd != imd ) {
		mexErrMsgTxt("Stacks must have the same size in the first three dimensions");
	} 
        
    /* Allocate memory for the new_dims array on the fly 
    out_dims = mxMalloc(nDims * sizeof(*out_dims));*/
    out_dims[0] = imw; 
    out_dims[1] = imh; 
    out_dims[2] = 3;
    
    plhs[0] = mxCreateNumericArray( 3, out_dims, mxSINGLE_CLASS,mxREAL );
    
    
    res = (float*)mxGetPr(plhs[0]);// the output 3D array of size [imw -2*rgw-1,imw -2*rgw-1, 5]
    matc = mxCreateNumericMatrix(tpl1, tpl2,mxSINGLE_CLASS,mxREAL);
	mat = (float*)mxGetPr(matc);// a matrix which the pixel correlation map
    
    
    //printf("Size little square %d x %d",tpl1,tpl2);
    //printf("Size image %d x %d",tpw,tph);
    //    printf("Element 10 %d x %f",tp[10],(float)tp[10]);

    
    // oi,oj are the pixel where we calculate the peak position
	// we use i, j for the neighbouring pixels	
	// k for the third dimension of the stack
    for (oi = 0; oi < imw ; oi++)
    {    
        for (oj = 0; oj < imh; oj++)
        {
            
            //normalize sample stack vector
            sum1 = 0;   std1 = 1;
            for (k = 0; k < imd; k++){
                sum1 += (float)im[oi + oj*imw + k*imw*imh];
            }
            sum1 = sum1/(float)tpd;//mean
            for (k = 0; k < imd; k++){
                        std1 += pow((float)im[oi + oj*imw + k*imw *imh] - sum1, 2);
            }
            std1 = sqrt(std1/((float)tpd -1));//unbiased standart deviation
            
            //for (k = 0; k < imd; k++){
            //    im[oi + oj*imw + k*imw *imh] = (im[oi + oj*imw + k*imw *imh] - sum1)/std1;//normalization of the vector elements
            //}
            
            //calc correlation coef of the sample vectors with the interrogation window ref vectors
            for (i = 0; i < tpl1; i++)
            {
                
                for (j = 0; j < tpl2; j++)
                {
                   corrf = 0;  sum2 = 0;  std2 = 1;
                   for (k = 0; k < imd; k++){
                        sum2 += (float)tp[oi + oj*imw + i*imw*imh + j*imw*imh*tpl1 + k*imw *imh*tpl1*tpl2];
                    }
                    sum2 = sum2/(float)tpd;
    
                    for (k = 0; k < imd; k++){
                        std2 += pow((float)tp[oi + oj*imw + i*imw *imh + j *imw *imh *tpl1 + k* imw *imh* tpl1 *tpl2] - sum2, 2);
                    }
                    std2 = sqrt(std2/((float)tpd -1));

                    for (k = 0; k < imd; k++){
                        corrf +=  ((float)im[oi + oj*imw + k*imw *imh] - sum1) * ((float)tp[oi + oj*imw + i*imw*imh + j*imw*imh*tpl1 + k*imw *imh*tpl1*tpl2] - sum2 )/(std2 * std1);            
                    }
                    
                    mat[i + j*tpl1] = corrf /imd ;
                }
            }         
            /* we have now a matrix with a peak of correlation for the pixel oi,oj*/
            

            /*========================pixel accuracy==============================*/    

            
            _intLengthX = tpl1;
            _intLengthY = tpl2;


            /* First, find the maximum value in the correlation map*/
            float max = -1000.0;
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
            float mean = 0.0;
            int count = 0;
            for (i = 0; i < _intLengthX; i++)
            {
                for (j = 0; j < _intLengthY; j++)
                {
                    if (!(i < (maxI+1) && i > (maxI-1) && j < (maxJ+1) && j > (maxJ-1)))
                    {
                        mean += fabs(mat[_intLengthX*i + j]);
                        count++;
                    }
                }
            }
			
			// 
            mean = mean / float(count);
            /* Assign the signal-to-noise (snr) value*/
            snr = float(max / mean);
			
			
			/*======================method 1 which fit 2 1D gaussian peaks===============*/
			
			
            /* Gaussian fitting to the peak in the correlation map*/
            //f0 = log(fabs(mat[_intLengthX*maxI + maxJ]));

            //f1 = log(fabs(mat[_intLengthX*(maxI-1) + maxJ]));
            //f2 = log(fabs(mat[_intLengthX*(maxI+1) + maxJ]));
            //yDisp = double(_intLengthY/2) - (double(maxI+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
            /* Assign the vertical displacment*/


            //f1 = log(fabs(mat[_intLengthX*maxI + maxJ - 1]));
            //f2 = log(fabs(mat[_intLengthX*maxI + maxJ + 1]));
            //xDisp = double(_intLengthX/2) - (double(maxJ+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
             /* Assign the horizontal displacement*/

            /*    plhs[1] = mxCreateDoubleMatrix( 5, 1, mxREAL );
             pic = mxGetPr(plhs[1]);
             printf("xDisp %f\n",xDisp);
             printf("yDisp %f\n",yDisp);
             printf("snr %f\n",snr);*/

            /*===================method 2: 1 2d peak==================*/

            float dx = (mat[_intLengthX*(maxI+1) + maxJ] - mat[_intLengthX*(maxI-1) + maxJ]) / 2.f;
            float dy = (mat[_intLengthX*maxI + maxJ + 1] - mat[_intLengthX*maxI + maxJ - 1]) / 2.f;
            float dxx = (mat[_intLengthX*(maxI+1) + maxJ] + mat[_intLengthX*(maxI-1) + maxJ] - 2 * mat[_intLengthX*maxI + maxJ]);
            float dyy = (mat[_intLengthX*maxI + maxJ + 1] + mat[_intLengthX*maxI + maxJ - 1] - 2 * mat[_intLengthX*maxI + maxJ]);
            float dxy = (mat[_intLengthX*(maxI + 1) + maxJ + 1] - mat[_intLengthX*(maxI + 1) + maxJ - 1] - mat[_intLengthX*(maxI - 1) + maxJ + 1] + mat[_intLengthX*(maxI - 1) + maxJ - 1]) / 4.f;

            float det = 1.f/(dxx*dyy - dxy*dxy);

            float ix = - (dyy*dx - dxy*dy) * det + maxI + 1 - (_intLengthY/2) ;  
            float iy = - (dxx*dy - dxy*dx) * det + maxJ + 1 - (_intLengthX/2) ;   

            /*     printf("xDisp %f\n",ix );*/
            /*     printf("yDisp %f\n",iy );*/
            //
            res[oi + oj*imw ] = iy;//-xDisp;
            res[oi + oj*imw + 1*imw *imh] = ix;//-yDisp;
            res[oi + oj*imw + 2*imw *imh] = max;//iy;
            //res[oi + oj*imw + 3*imw *imh] = ix;
            //res[oi + oj*imw + 4*imw *imh] = max;//snr;
        
             
       }
        //printf("ok");
    }

	
	
};



