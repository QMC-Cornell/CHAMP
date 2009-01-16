#ifndef MULTI_BSPLINE_CREATE_CUDA_H
#define MULTI_BSPLINE_CREATE_CUDA_H

#include "multi_bspline_structs_cuda.h"

extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda (multi_UBspline_3d_s* spline);

extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda_conv (multi_UBspline_3d_d* spline);

extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda (multi_UBspline_3d_c* spline);

extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda_conv (multi_UBspline_3d_z* spline);

#endif
