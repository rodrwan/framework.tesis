#include <string.h>
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <cmath>
// small value, used to avoid division by zero
#define eps 0.0001

// unit vectors used to compute gradient orientation
double uu[9] = {1.0000,
    0.9397,
    0.7660,
    0.500,
    0.1736,
    -0.1736,
    -0.5000,
    -0.7660,
    -0.9397};
double vv[9] = {0.0000,
    0.3420,
    0.6428,
    0.8660,
    0.9848,
    0.9848,
    0.8660,
    0.6428,
    0.3420};

static inline double min(double x, double y) { return (x <= y ? x : y); }
// static inline double max(double x, double y) { return (x <= y ? y : x); }

static inline int min(int x, int y) { return (x <= y ? x : y); }
static inline int max(int x, int y) { return (x <= y ? y : x); }

/**
 * [max_channel get the max values between channels]
 * @param  ch1 [Ch1 value]
 * @param  ch2 [Ch2 value]
 * @param  ch3 [Ch3 value]
 * @return     [Max value]
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

/**
 * [final_response description]
 * @param  gabors   [Array with gabor convolution]
 * @param  sbin     [Bin]
 * @param  idims    [Image dimention]
 * @param  uv       [Filter count]
 * @return          [descriptor]
 */
mxArray *final_response(mxArray *gabors, const int sbin, const int *idims, const int *puv) {
  int out[3];
  out[0] = (int)round((double)idims[0]/(double)sbin);
  out[1] = (int)round((double)idims[1]/(double)sbin);
  out[2] = puv[0]*puv[1];
  mxArray *mxfeat = mxCreateNumericArray(3, out, mxDOUBLE_CLASS, mxREAL);
  float *feat = (float *)mxGetPr(mxfeat);

  float *g_response = (float *)mxGetPr(gabors);
  int Mi = idims[0];
  int Mc = out[0];
  int Nc = out[1];
  int cell_size = Mc*Nc;
  int uv, i, it, offset;

  for (int ii = 0; ii < cell_size; ii++) {
    float max = 0.0;
    uv = 0;
    // index transformation.
    i = ii + (int)(ii/Mc)*(Mi-Mc);
    it = ((int)(i/Nc))%Nc + ((int)(i/(Mi*Mc)) * ((int)(Mi/Nc)));
    // loop through each filter position
    for (int j = 0; j < out[2]; j ++) {
      float gabor_value = *(g_response + i + j*cell_size);
      if (gabor_value > max) {
        max = gabor_value;
        uv = j;
      }
    }
    offset = it + uv*cell_size;
    float *dst = feat + offset;
    *(dst) += max;
  }
  return mxfeat;
}


/**
 * [hog descriptor]
 * @param  im   [Image vector]
 * @param  sbin [Bin size]
 * @param  dims [Image dimention]
 * @return      [Descriptor for image input]
 */
mxArray *hog(double *im, const int sbin, const int *dims) {
  // memory for caching orientation histograms & their norms
  int cells[2];
  cells[0] = (int)round((double)dims[0]/(double)sbin);
  cells[1] = (int)round((double)dims[1]/(double)sbin);

  float *hist = (float *)mxCalloc(cells[0]*cells[1]*18, sizeof(float));
  float *norm = (float *)mxCalloc(cells[0]*cells[1], sizeof(float));

  // memory for HOG features
  int out[3];
  out[0] = max(cells[0]-2, 0);
  out[1] = max(cells[1]-2, 0);
  out[2] = 27+4+1;
  mxArray *mxfeat = mxCreateNumericArray(3, out, mxSINGLE_CLASS, mxREAL);
  float *feat = (float *)mxGetPr(mxfeat);

  int visible[2];
  visible[0] = cells[0]*sbin;
  visible[1] = cells[1]*sbin;

  for (int x = 1; x < visible[1]-1; x++) {
    for (int y = 1; y < visible[0]-1; y++) {
      // first color channel
      double *s = im + min(x, dims[1]-2)*dims[0] + min(y, dims[0]-2);
      double dy = *(s+1) - *(s-1);
      double dx = *(s+dims[0]) - *(s-dims[0]);
      double v = dx*dx + dy*dy;

      // second color channel
      s += dims[0]*dims[1];
      double dy2 = *(s+1) - *(s-1);
      double dx2 = *(s+dims[0]) - *(s-dims[0]);
      double v2 = dx2*dx2 + dy2*dy2;

      // third color channel
      s += dims[0]*dims[1];
      double dy3 = *(s+1) - *(s-1);
      double dx3 = *(s+dims[0]) - *(s-dims[0]);
      double v3 = dx3*dx3 + dy3*dy3;

      // pick channel with strongest gradient
      if (v2 > v) {
        v = v2;
        dx = dx2;
        dy = dy2;
      }
      if (v3 > v) {
        v = v3;
        dx = dx3;
        dy = dy3;
      }

      // snap to one of 18 orientations
      double best_dot = 0;
      int best_o = 0;
      for (int o = 0; o < 9; o++) {
        double dot = uu[o]*dx + vv[o]*dy;
        if (dot > best_dot) {
          best_dot = dot;
          best_o = o;
        } else if (-dot > best_dot) {
          best_dot = -dot;
          best_o = o+9;
        }
      }

      // add to 4 histograms around pixel using bilinear interpolation
      double xp = ((double)x+0.5)/(double)sbin - 0.5;
      double yp = ((double)y+0.5)/(double)sbin - 0.5;
      int ixp = (int)floor(xp);
      int iyp = (int)floor(yp);
      double vx0 = xp-ixp;
      double vy0 = yp-iyp;
      double vx1 = 1.0-vx0;
      double vy1 = 1.0-vy0;
      v = sqrt(v);

      if (ixp >= 0 && iyp >= 0) {
        *(hist + ixp*cells[0] + iyp + best_o*cells[0]*cells[1]) +=
          vx1*vy1*v;
      }

      if (ixp+1 < cells[1] && iyp >= 0) {
        *(hist + (ixp+1)*cells[0] + iyp + best_o*cells[0]*cells[1]) +=
          vx0*vy1*v;
      }

      if (ixp >= 0 && iyp+1 < cells[0]) {
        *(hist + ixp*cells[0] + (iyp+1) + best_o*cells[0]*cells[1]) +=
          vx1*vy0*v;
      }

      if (ixp+1 < cells[1] && iyp+1 < cells[0]) {
        *(hist + (ixp+1)*cells[0] + (iyp+1) + best_o*cells[0]*cells[1]) +=
          vx0*vy0*v;
      }
    }
  }

  // compute energy in each block by summing over orientations
  for (int o = 0; o < 9; o++) {
    float *src1 = hist + o*cells[0]*cells[1];
    float *src2 = hist + (o+9)*cells[0]*cells[1];
    float *dst = norm;
    float *end = norm + cells[1]*cells[0];
    while (dst < end) {
      *(dst++) += (*src1 + *src2) * (*src1 + *src2);
      src1++;
      src2++;
    }
  }

  // compute features
  for (int x = 0; x < out[1]; x++) {
    for (int y = 0; y < out[0]; y++) {
      float *dst = feat + x*out[0] + y;
      float *src, *p, n1, n2, n3, n4;

      p = norm + (x+1)*cells[0] + y+1;
      n1 = 1.0 / sqrt(*p + *(p+1) + *(p+cells[0]) + *(p+cells[0]+1) + eps);
      p = norm + (x+1)*cells[0] + y;
      n2 = 1.0 / sqrt(*p + *(p+1) + *(p+cells[0]) + *(p+cells[0]+1) + eps);
      p = norm + x*cells[0] + y+1;
      n3 = 1.0 / sqrt(*p + *(p+1) + *(p+cells[0]) + *(p+cells[0]+1) + eps);
      p = norm + x*cells[0] + y;
      n4 = 1.0 / sqrt(*p + *(p+1) + *(p+cells[0]) + *(p+cells[0]+1) + eps);

      float t1 = 0;
      float t2 = 0;
      float t3 = 0;
      float t4 = 0;

      // contrast-sensitive features
      src = hist + (x+1)*cells[0] + (y+1);
      for (int o = 0; o < 18; o++) {
        float h1 = min(*src * n1, 0.2);
        float h2 = min(*src * n2, 0.2);
        float h3 = min(*src * n3, 0.2);
        float h4 = min(*src * n4, 0.2);
        *dst = 0.5 * (h1 + h2 + h3 + h4);
        t1 += h1;
        t2 += h2;
        t3 += h3;
        t4 += h4;
        dst += out[0]*out[1];
        src += cells[0]*cells[1];
      }

      // contrast-insensitive features
      src = hist + (x+1)*cells[0] + (y+1);
      for (int o = 0; o < 9; o++) {
        float sum = *src + *(src + 9*cells[0]*cells[1]);
        float h1 = min(sum * n1, 0.2);
        float h2 = min(sum * n2, 0.2);
        float h3 = min(sum * n3, 0.2);
        float h4 = min(sum * n4, 0.2);
        *dst = 0.5 * (h1 + h2 + h3 + h4);
        dst += out[0]*out[1];
        src += cells[0]*cells[1];
      }

      // texture features
      *dst = 0.2357 * t1;
      dst += out[0]*out[1];
      *dst = 0.2357 * t2;
      dst += out[0]*out[1];
      *dst = 0.2357 * t3;
      dst += out[0]*out[1];
      *dst = 0.2357 * t4;

      // truncation feature
      dst += out[0]*out[1];
      *dst = 0;
    }
  }

  mxFree(hist);
  mxFree(norm);
  return mxfeat;
}
// main function:
// takes a double color image and a bin size
// returns HOG features
mxArray *process(const mxArray *mximage, const mxArray *mxsbin,
                 const mxArray *filterArray, const mxArray *type,
                 const mxArray *verbose) {

  double *im = (double *)mxGetPr(mximage);
  const int *imdims = mxGetDimensions(mximage);
  int sbin = (int)mxGetScalar(mxsbin);
  char *feat_type = mxArrayToString(type);

  if (strcmp(feat_type, "hog") == 0) {
    if (mxGetNumberOfDimensions(mximage) != 3 || imdims[2] != 3 || mxGetClassID(mximage) != mxDOUBLE_CLASS)
      mexErrMsgTxt("Invalid input");

    return hog(im, sbin, imdims);
  } else {
    if (mxGetNumberOfDimensions(mximage) != 3 || mxGetClassID(mximage) != mxDOUBLE_CLASS)
      mexErrMsgTxt("Invalid input");
    const int *fdims = mxGetDimensions(filterArray);
    const mxArray *filter = filterArray;
    int uv[2];
    uv[0] = fdims[0];
    uv[1] = fdims[1];
    mxArray *gabor_response = gabor(im, imdims, filter, fdims);
    return final_response(gabor_response, sbin, imdims, uv);
  }
}

// matlab entry point
// F = features(image, bin)
// image should be color with double values
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 5)
    mexErrMsgTxt("Wrong number of inputs");
  if (nlhs != 1)
    mexErrMsgTxt("Wrong number of outputs");
  plhs[0] = process(prhs[0], prhs[1], prhs[2], prhs[3], prhs[4]);
}



