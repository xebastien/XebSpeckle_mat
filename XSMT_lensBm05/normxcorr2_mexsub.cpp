#include "mex.h"
#include "cv_src/cv.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, oi, oj, k;

	int imw, imh;
	int tpw, tph;
	int resw, resh;

	int outw, outh; // real size of output (taking into account 'shape')
	int offw, offh;

	double* im;
	double* tp;
	double* res;
    double* pic;
    
	IplImage* cvIm;
	IplImage* cvTp;
	IplImage* cvRes;

	int shapeTypeCharSize = 6;
	char* shapeTypeChar;
	int shapeType = 1; 
	/* 1, full
	   2, valid
	   3, same
	   */

    int _intLengthX;
    int _intLengthY;

    int maxI, maxJ;
    double f0, f1, f2, snr;
    double xDisp, yDisp;
    
	float* cvImP;
	float* cvTpP;
	float* cvResP;

	if( nrhs < 2 ) {
		mexErrMsgTxt("result = normxcorr2Mex(TEMPLATE, IMAGE)");
	} else if( !( mxIsDouble(prhs[0])&&mxIsDouble(prhs[1])) ) {
		mexErrMsgTxt("IMAGE and TEMPLATE must be of type double");
	} else if( nrhs > 2 && !mxIsChar(prhs[2]) ) {
		mexErrMsgTxt("SHAPE parameter must be a character string");
	}

	shapeTypeChar = (char*)malloc(shapeTypeCharSize);
	if( nrhs>2 ) {
		mxGetString( prhs[2], shapeTypeChar, shapeTypeCharSize );
		for( int ci=0; ci<strlen(shapeTypeChar); ci++ )
			shapeTypeChar[ci] = tolower(shapeTypeChar[ci]);

		if( !strcmp(shapeTypeChar, "full") )
			shapeType = 1;
		else if( !strcmp(shapeTypeChar, "valid") )
			shapeType = 2;
		else if( !strcmp(shapeTypeChar, "same") )
			shapeType = 3;
		else
			shapeType = 0;
	}

	free(shapeTypeChar);
	if( !shapeType )
		mexErrMsgTxt("unknown shape type");

	im = mxGetPr(prhs[1]);
	imw = mxGetN(prhs[1]);
	imh = mxGetM(prhs[1]);

	tp = mxGetPr(prhs[0]);
	tpw = mxGetN(prhs[0]);
	tph = mxGetM(prhs[0]);

	resw = imw-tpw+1;
	resh = imh-tph+1;

	if( resw<=0 || resh<=0 ) {
		mexErrMsgTxt("size(image) < size(template)");
	} 

	cvIm = cvCreateImage(cvSize(imw,imh), IPL_DEPTH_32F, 1);
	cvTp = cvCreateImage(cvSize(tpw,tph), IPL_DEPTH_32F, 1);
	cvRes = cvCreateImage(cvSize(resw,resh), IPL_DEPTH_32F, 1);

	cvImP = (float*)cvIm->imageData;
	cvTpP = (float*)cvTp->imageData;
	cvResP = (float*)cvRes->imageData;

	if( shapeType==1 ) { //full
		outw = imw+tpw-1;
		outh = imh+tph-1;
		offw = tpw-1;
		offh = tph-1;
	} else if( shapeType==2 ) { //valid
		outw = resw;
		outh = resh;
		offw = 0;
		offh = 0;
	} else if( shapeType==3 ) { // same
		outw = imw;
		outh = imh;
		offw = ceil((float)(tpw-1)/2.0f);
		offh = ceil((float)(tph-1)/2.0f);
	}

	plhs[0] = mxCreateDoubleMatrix( outh, outw, mxREAL );
	res = mxGetPr(plhs[0]);

	// we don't need to worry about alignment since we're using 32f

	k = 0;
	for( i=0; i<imw; i++ ) {
		for( j=0; j<imh; j++ ) {
			cvImP[ j*imw+i ] = im[ k ];
			k++;
		}
	}

	k = 0;
	for( i=0; i<tpw; i++ ) {
		for( j=0; j<tph; j++ ) {
			cvTpP[ j*tpw+i ] = tp[ k ];
			k++;
		}
	}

	cvMatchTemplate( cvIm, cvTp, cvRes, CV_TM_CCOEFF_NORMED );

	for( i=0, oi=offw; i<resw; i++, oi++ ) {
		for( j=0, oj=offh; j<resh; j++, oj++ ) {
			res[ oi*outh + oj ] = cvResP[ j*resw + i ];
		}
	}

	cvReleaseImage(&cvIm);
	cvReleaseImage(&cvTp);
	cvReleaseImage(&cvRes);
    
    
//========================subpixel accuracy==============================//    
    
    _intLengthX = outh;
    _intLengthY = outw;


    // First, find the maximum value in the correlation map
    double max = -1000.0;
    for (i = 0; i < _intLengthY; i++)
    {
        for (j = 0; j < _intLengthX; j++)
        {
            if (res[_intLengthX*i + j] > max)
            {
                max = res[_intLengthX*i + j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    
//     printf("wdith %d\n", _intLengthX);
//     printf("length %d\n",_intLengthY);
//     printf("maxI %d\n",maxI);
//     printf("maxJ %d\n",maxJ);
    // Need to calculate the average to which to compare the signal
    double mean = 0.0;
    int count = 0;
    for (i = 0; i < _intLengthY; i++)
    {
        for (j = 0; j < _intLengthX; j++)
        {
            if (!(i < (maxI+1) && i > (maxI-1) && j < (maxJ+1) && j < (maxJ-1)))
            {
                mean += fabs(res[_intLengthX*i + j]);
                count++;
            }
        }
    }
    mean = mean / double(count);
    // Assign the signal-to-noise (snr) value
    snr = double(max / mean);

    // Gaussian fitting to the peak in the correlation map
    f0 = log(fabs(res[_intLengthX*maxI + maxJ]));

    f1 = log(fabs(res[_intLengthX*(maxI-1) + maxJ]));
    f2 = log(fabs(res[_intLengthX*(maxI+1) + maxJ]));
    yDisp = double(_intLengthY/2) - (double(maxI+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
    // Assign the vertical displacment


    f1 = log(fabs(res[_intLengthX*maxI + maxJ - 1]));
    f2 = log(fabs(res[_intLengthX*maxI + maxJ + 1]));
    xDisp = double(_intLengthX/2) - (double(maxJ+1) + (f1 - f2) / (2*f1 - 4*f0 + 2*f2));
    // Assign the horizontal displacement

    plhs[1] = mxCreateDoubleMatrix( 5, 1, mxREAL );
    pic = mxGetPr(plhs[1]);
//     printf("xDisp %f\n",xDisp);
//     printf("yDisp %f\n",yDisp);
//     printf("snr %f\n",snr);
    pic[0] = -xDisp;
    pic[1] = -yDisp;
    pic[4] = snr;
    
    //=====================================
    
    float dx = (res[_intLengthX*(maxI+1) + maxJ] - res[_intLengthX*(maxI-1) + maxJ]) / 2.f;
    float dy = (res[_intLengthX*maxI + maxJ + 1] - res[_intLengthX*maxI + maxJ - 1]) / 2.f;
    float dxx = (res[_intLengthX*(maxI+1) + maxJ] + res[_intLengthX*(maxI-1) + maxJ] - 2 * res[_intLengthX*maxI + maxJ]);
    float dyy = (res[_intLengthX*maxI + maxJ + 1] + res[_intLengthX*maxI + maxJ - 1] - 2 * res[_intLengthX*maxI + maxJ]);
    float dxy = (res[_intLengthX*(maxI + 1) + maxJ + 1] - res[_intLengthX*(maxI + 1) + maxJ - 1] - res[_intLengthX*(maxI - 1) + maxJ + 1] + res[_intLengthX*(maxI - 1) + maxJ - 1]) / 4.f;

    float det = 1.f/(dxx*dyy - dxy*dxy);

    float ix = - (dyy*dx - dxy*dy) * det + maxI + 1 - double(_intLengthY/2) ;  
    float iy = - (dxx*dy - dxy*dx) * det + maxJ + 1 - double(_intLengthX/2) ;   
    
    //     printf("xDisp %f\n",ix );
    //     printf("yDisp %f\n",iy );

    pic[2] = iy;
    pic[3] = ix;
    // Note that at this point we know nothing about the location (x,y) of the data point.
    // That value is assigned elsewhere.
};