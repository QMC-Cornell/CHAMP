#ifndef MULTI_BSPLINE_EVAL_CUDA_H
#define MULTI_BSPLINE_EVAL_CUDA_H

#include "multi_bspline_structs_cuda.h"

extern "C" void
eval_multi_multi_UBspline_3d_s_cuda 
(multi_UBspline_3d_s_cuda *spline, float *pos_d, float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_s_vgh_cuda 
(multi_UBspline_3d_s_cuda *spline,
 float *pos_d, float *vals_d[], float *grads_d[], float *hess_d[], int num);


extern "C" void
eval_multi_multi_UBspline_3d_c_cuda 
(multi_UBspline_3d_c_cuda *spline, 
 float *pos_d, complex_float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_c_vgh_cuda 
(multi_UBspline_3d_c_cuda *spline, float *pos_d, 
 complex_float *vals_d[], complex_float *grads_d[], 
 complex_float *hess_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_s_vgl_cuda
(multi_UBspline_3d_s_cuda *spline, float *pos_d, float *Linv_d, 
 float *vals_d[], float *grad_lapl_d[], int num, int row_stride);

#endif
