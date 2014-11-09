#include <string.h>
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <cmath>

/**
 * [max_channel get the max values between channels]
 * @param  ch1 [Channel 1 value]
 * @param  ch2 [Channel 2 value]
 * @param  ch3 [Channel 3 value]
 * @return     [Max value between channels]
 */
double max_channel(const double ch1, const double ch2, const double ch3) {
  double max = ch1;
  if (ch2 > max) max = ch2;
  if (ch3 > max) max = ch3;
  return max;
}

/**
 * [filter_response description]
 * @param  img    [Image vector]
 * @param  im_i   [Image current index]
 * @param  idims  [Image dimention]
 * @param  filter [Filter vector]
 * @param  fdims  [Filter dimention]
 * @return        [Convolution value]
 */
double filter_response(double *img, const int im_i, const int *idims,
                       double *g, const int *fdims) {
  int m_fil = fdims[0];
  int n_fil = fdims[1];
  int m_img = idims[0];
  int n_img = idims[1];

  double acc = 0, f_val, i_val, i_val_ch1, i_val_ch2, i_val_ch3;
  for (int i = 0; i < m_fil*n_fil; i++) {
    // filter access pointer
    f_val = std::abs(*(g+i));
    // images access pointer
    double *im = img + im_i;
    // image channel 1
    i_val_ch1 = std::abs(*(im + (i%m_fil) + (int)(i/m_fil)*m_img));
    // image channel 2
    i_val_ch2 = std::abs(*(im + m_img*n_img + (i%m_fil) + (int)(i/m_fil)*m_img));
    // image channel 3
    i_val_ch3 = std::abs(*(im + m_img*n_img + (i%m_fil) + (int)(i/m_fil)*m_img));
    // best response
    i_val = max_channel(i_val_ch1, i_val_ch2, i_val_ch3);
    acc += f_val * i_val;
  }
  return acc;
}

/**
 * [gabor description]
 * @param  img    [Imagen vector]
 * @param  idims  [Image dimension]
 * @param  filter [Filter vector]
 * @param  fdims  [Filter dimension]
 * @return        [Gabors vector]
 */
mxArray *gabor(double *img, const int *idims, const mxArray *filter, const int *fdims) {
  int uv = fdims[0]*fdims[1];
  int out[3];
  out[0] = idims[0];
  out[1] = idims[1];
  out[2] = uv+1;
  mxArray *mxGabor = mxCreateNumericArray(3, out, mxDOUBLE_CLASS, mxREAL);
  float *gabors = (float *)mxGetPr(mxGabor);

  // loop through each filter
  for (int fi = 0; fi < uv; fi++) {
    // Get current filter from filter pool
    mxArray *current_filter = mxGetCell(filter, fi);
    // Pointer for current filter
    double *cf = (double *)mxGetPr(current_filter);
    // Filter dimention
    const int *fcdims = mxGetDimensions(current_filter);
    // center of the filter
    int center_x = (int)(fcdims[0]/2);
    int center_y = (int)(fcdims[1]/2);
    double resp;
    int img_size = (idims[0]*idims[1])-2;
    for (int i = 0; i < img_size; i++) {
      resp = filter_response(img, i, idims, cf, fcdims);
      int offset = i + center_x + (center_y*idims[0]) + img_size*fi;
      *(gabors + offset) = resp;
    }
  }
  return mxGabor;
}

mxArray *process(const mxArray *mximage, const mxArray *filterArray) {

  double *im = (double *)mxGetPr(mximage);
  const int *imdims = mxGetDimensions(mximage);

  if (mxGetNumberOfDimensions(mximage) != 3 || mxGetClassID(mximage) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Invalid input");

  const int *fdims = mxGetDimensions(filterArray);
  const mxArray *filter = filterArray;
  int uv[2];
  uv[0] = fdims[0];
  uv[1] = fdims[1];
  return gabor(im, imdims, filter, fdims);
}

// matlab entry point
// F = features(image, bin)
// image should be color with double values
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 2)
    mexErrMsgTxt("Wrong number of inputs");
  if (nlhs != 1)
    mexErrMsgTxt("Wrong number of outputs");
  plhs[0] = process(prhs[0], prhs[1]);
}