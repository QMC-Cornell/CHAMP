#include "multi_bspline.h"
#include "multi_bspline_structs_cuda.h"

__constant__ float Acuda[48];

#include "multi_bspline_cuda_s_impl.h"
#include "multi_bspline_cuda_c_impl.h"

// typedef struct
// {
//   float *coefs;
//   uint3 stride;
//   float3 gridInv;
//   int num_splines;
// } multi_UBspline_3d_s_cuda;

// typedef struct
// {
//   float *coefs_real, *coefs_imag;
//   uint3 stride;
//   float3 gridInv;
//   int num_splines;
// } multi_UBspline_3d_c_cuda;


extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda (multi_UBspline_3d_c* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  multi_UBspline_3d_c_cuda *cuda_spline =
    (multi_UBspline_3d_c_cuda*) malloc (sizeof (multi_UBspline_3d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%SPLINE_BLOCK_SIZE) != 0)
    N += 64 - (N%SPLINE_BLOCK_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs_real), size);
  cudaMalloc((void**)&(cuda_spline->coefs_imag), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].real();
	}
  cudaMemcpy(cuda_spline->coefs_real, spline_buff, size, cudaMemcpyHostToDevice);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].imag();
	}
  cudaMemcpy(cuda_spline->coefs_imag, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}


extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda_conv (multi_UBspline_3d_z* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  multi_UBspline_3d_c_cuda *cuda_spline =
    (multi_UBspline_3d_c_cuda*) malloc (sizeof (multi_UBspline_3d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%SPLINE_BLOCK_SIZE) != 0)
    N += 64 - (N%SPLINE_BLOCK_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs_real), size);
  cudaMalloc((void**)&(cuda_spline->coefs_imag), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    (float)spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].real();
	}
  cudaMemcpy(cuda_spline->coefs_real, spline_buff, size, cudaMemcpyHostToDevice);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    (float)spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp].imag();
	}
  cudaMemcpy(cuda_spline->coefs_imag, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}




extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda (multi_UBspline_3d_s* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  multi_UBspline_3d_s_cuda *cuda_spline =
    (multi_UBspline_3d_s_cuda*) malloc (sizeof (multi_UBspline_3d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%SPLINE_BLOCK_SIZE) != 0)
    N += 64 - (N%SPLINE_BLOCK_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}



extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda_conv (multi_UBspline_3d_d* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  multi_UBspline_3d_s_cuda *cuda_spline =
    (multi_UBspline_3d_s_cuda*) malloc (sizeof (multi_UBspline_3d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%SPLINE_BLOCK_SIZE) != 0)
    N += 64 - (N%SPLINE_BLOCK_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}




