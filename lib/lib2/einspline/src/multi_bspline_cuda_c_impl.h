#ifndef MULTI_BSPLINE_CUDA_C_IMPL_H
#define MULTI_BSPLINE_CUDA_C_IMPL_H

#include "multi_bspline.h"
#include "multi_bspline_create_cuda.h"

__global__ static void
eval_multi_multi_UBspline_3d_c_cuda (float *pos, float3 drInv, 
				     float *coefs_real, float *coefs_imag,
				     float *vals[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;

  __shared__ float *myval;
  __shared__ float abc[64], coefs[2*SPLINE_BLOCK_SIZE];

  // __shared__ float pos_s[SPLINE_BLOCK_SIZE];
  // int ir1 = (ir >> 4)*64;
  // int ir2 = (ir & 15)*4;
  // pos_s[thr] = pos[ir1+thr];
  // __syncthreads();
  // float3 r;
  // r.x = pos_s[ir2+0];
  // r.y = pos_s[ir2+1];
  // r.z = pos_s[ir2+2];
  __shared__ float3 r;
  if (thr == 0) {
    r.x = pos[4*ir+0];
    r.y = pos[4*ir+1];
    r.z = pos[4*ir+2];
    myval = vals[ir];
  }
  __syncthreads();
  
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];

  s = r.x * drInv.x;
  sf = floor(s);
  index.x = (int)sf;
  t.x = s - sf;

  s = r.y * drInv.y;
  sf = floor(s);
  index.y = (int)sf;
  t.y = s - sf;

  s = r.z * drInv.z;
  sf = floor(s);
  index.z = (int)sf;
  t.z = s - sf;
  
  tp[0] = make_float4(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = make_float4(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = make_float4(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);

  __shared__ float a[4], b[4], c[4];
  if (thr < 4) {
    a[thr] = Acuda[4*thr+0]*tp[0].x + Acuda[4*thr+1]*tp[0].y + Acuda[4*thr+2]*tp[0].z + Acuda[4*thr+3]*tp[0].w;
    b[thr] = Acuda[4*thr+0]*tp[1].x + Acuda[4*thr+1]*tp[1].y + Acuda[4*thr+2]*tp[1].z + Acuda[4*thr+3]*tp[1].w;
    c[thr] = Acuda[4*thr+0]*tp[2].x + Acuda[4*thr+1]*tp[2].y + Acuda[4*thr+2]*tp[2].z + Acuda[4*thr+3]*tp[2].w;
  }
  __syncthreads();

  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);
  
  abc[thr] = a[i]*b[j]*c[k];
  __syncthreads();


  float val_real = 0.0;
  float val_imag = 0.0;
  val_real = val_imag = 0.0;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      for (int k=0; k<4; k++) {
  	val_real += abc[16*i+4*j+k] * base_real[off+k*strides.z];
  	val_imag += abc[16*i+4*j+k] * base_imag[off+k*strides.z];
      }
    }
  }
  __shared__ float buff[2*SPLINE_BLOCK_SIZE];
  buff[2*thr+0] = val_real;
  buff[2*thr+1] = val_imag;
  __syncthreads();
  myval[off] = buff[thr];
  myval[off+SPLINE_BLOCK_SIZE] = buff[thr+SPLINE_BLOCK_SIZE];
}



__global__ static void
eval_multi_multi_UBspline_3d_c_vgh_cuda (float *pos, float3 drInv, 
					 float *coefs_real, float *coefs_imag,
					 float *vals[], float *grads[], 
					 float *hess[], uint3 strides)
{
  int block = blockIdx.x;
  int thr   = threadIdx.x;
  int ir    = blockIdx.y;
  int off   = block*SPLINE_BLOCK_SIZE+thr;

  __shared__ float *myval, *mygrad, *myhess;
  __shared__ float3 r;
  if (thr == 0) {
    r.x = pos[4*ir+0];
    r.y = pos[4*ir+1];
    r.z = pos[4*ir+2];
    myval  = vals[ir];
    mygrad = grads[ir];
    myhess = hess[ir];
  }
  __syncthreads();
  
  int3 index;
  float3 t;
  float s, sf;
  float4 tp[3];

  s = r.x * drInv.x;
  sf = floor(s);
  index.x = (int)sf;
  t.x = s - sf;

  s = r.y * drInv.y;
  sf = floor(s);
  index.y = (int)sf;
  t.y = s - sf;

  s = r.z * drInv.z;
  sf = floor(s);
  index.z = (int)sf;
  t.z = s - sf;
  
  tp[0] = make_float4(t.x*t.x*t.x, t.x*t.x, t.x, 1.0);
  tp[1] = make_float4(t.y*t.y*t.y, t.y*t.y, t.y, 1.0);
  tp[2] = make_float4(t.z*t.z*t.z, t.z*t.z, t.z, 1.0);

  // First 4 of a are value, second 4 are derivative, last four are
  // second derivative.
  __shared__ float a[12], b[12], c[12];
  if (thr < 12) {
    a[thr] = Acuda[4*thr+0]*tp[0].x + Acuda[4*thr+1]*tp[0].y + Acuda[4*thr+2]*tp[0].z + Acuda[4*thr+3]*tp[0].w;
    b[thr] = Acuda[4*thr+0]*tp[1].x + Acuda[4*thr+1]*tp[1].y + Acuda[4*thr+2]*tp[1].z + Acuda[4*thr+3]*tp[1].w;
    c[thr] = Acuda[4*thr+0]*tp[2].x + Acuda[4*thr+1]*tp[2].y + Acuda[4*thr+2]*tp[2].z + Acuda[4*thr+3]*tp[2].w;
  }
  __syncthreads();

  __shared__ float abc[640];
  int i = (thr>>4)&3;
  int j = (thr>>2)&3;
  int k = (thr & 3);

  abc[10*(16*i+4*j+k)+0] = a[i+0]*b[j+0]*c[k+0]; // val
  abc[10*(16*i+4*j+k)+1] = a[i+4]*b[j+0]*c[k+0]; // d/dx
  abc[10*(16*i+4*j+k)+2] = a[i+0]*b[j+4]*c[k+0]; // d/dy
  abc[10*(16*i+4*j+k)+3] = a[i+0]*b[j+0]*c[k+4]; // d/dz
  abc[10*(16*i+4*j+k)+4] = a[i+8]*b[j+0]*c[k+0]; // d2/dx2
  abc[10*(16*i+4*j+k)+5] = a[i+4]*b[j+4]*c[k+0]; // d2/dxdy
  abc[10*(16*i+4*j+k)+6] = a[i+4]*b[j+0]*c[k+4]; // d2/dxdz
  abc[10*(16*i+4*j+k)+7] = a[i+0]*b[j+8]*c[k+0]; // d2/dy2
  abc[10*(16*i+4*j+k)+8] = a[i+0]*b[j+4]*c[k+4]; // d2/dydz
  abc[10*(16*i+4*j+k)+9] = a[i+0]*b[j+0]*c[k+8]; // d2/dz2

  __syncthreads();

  float v_r = 0.0;
  float v_i = 0.0;
  float g0_r=0.0, g0_i=0.0, g1_r=0.0, g1_i=0.0, g2_r=0.0, g2_i=0.0, 
    h00_r=0.0, h00_i=0.0, h01_r=0.0, h01_i=0.0, h02_r=0.0, h02_i=0.0, 
    h11_r=0.0, h11_i=0.0, h12_r=0.0, h12_i=0.0, h22_r=0.0, h22_i=0.0;
  int n = 0;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      float *base_real = coefs_real + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;
      float *base_imag = coefs_imag + (index.x+i)*strides.x + (index.y+j)*strides.y + index.z*strides.z;

      for (int k=0; k<4; k++) {
	float cr = base_real[off+k*strides.z];
	float ci = base_imag[off+k*strides.z];
	v_r   += abc[n+0] * cr;  v_i   += abc[n+0] * ci;
	g0_r  += abc[n+1] * cr;  g0_i  += abc[n+1] * ci;
	g1_r  += abc[n+2] * cr;  g1_i  += abc[n+2] * ci;
	g2_r  += abc[n+3] * cr;  g2_i  += abc[n+3] * ci;
	h00_r += abc[n+4] * cr;  h00_i += abc[n+4] * ci;
	h01_r += abc[n+5] * cr;  h01_i += abc[n+5] * ci;
	h02_r += abc[n+6] * cr;  h02_i += abc[n+6] * ci;
	h11_r += abc[n+7] * cr;  h11_i += abc[n+7] * ci;
	h12_r += abc[n+8] * cr;  h12_i += abc[n+8] * ci;
	h22_r += abc[n+9] * cr;  h22_i += abc[n+9] * ci; 
	n += 10;
      }
    }
  }
  g0_r *= drInv.x; g0_i *= drInv.x;
  g1_r *= drInv.y; g1_i *= drInv.y;
  g2_r *= drInv.z; g2_i *= drInv.z;

  h00_r *= drInv.x * drInv.x;  h00_i *= drInv.x * drInv.x;
  h01_r *= drInv.x * drInv.y;  h01_i *= drInv.x * drInv.y;
  h02_r *= drInv.x * drInv.z;  h02_i *= drInv.x * drInv.z;
  h11_r *= drInv.y * drInv.y;  h11_i *= drInv.y * drInv.y;
  h12_r *= drInv.y * drInv.z;  h12_i *= drInv.y * drInv.z;
  h22_r *= drInv.z * drInv.z;  h22_i *= drInv.z * drInv.z;

  
  __shared__ float buff[6*SPLINE_BLOCK_SIZE];
  // Note, we can reuse abc, by replacing buff with abc.
  
  buff[2*thr+0] = v_r;  buff[2*thr+1] = v_i;
  __syncthreads();
  myval[off] = buff[thr];    
  myval[off+SPLINE_BLOCK_SIZE] = buff[thr+SPLINE_BLOCK_SIZE];

  buff[6*thr+0] = g0_r;  buff[6*thr+1] = g0_i;
  buff[6*thr+2] = g1_r;  buff[6*thr+3] = g1_i;
  buff[6*thr+4] = g2_r;  buff[6*thr+5] = g2_i;
  __syncthreads();
  for (int i=0; i<6; i++) 
    mygrad[(6*block+i)*SPLINE_BLOCK_SIZE+thr] = buff[i*SPLINE_BLOCK_SIZE+thr]; 
  __syncthreads();

  // Write first half of Hessians
  if (thr < 32) {
    buff[12*thr+0]  = h00_r;    buff[12*thr+1]  = h00_i;
    buff[12*thr+2]  = h01_r;    buff[12*thr+3]  = h01_i;
    buff[12*thr+4]  = h02_r;    buff[12*thr+5]  = h02_i;
    buff[12*thr+6]  = h11_r;    buff[12*thr+7]  = h11_i;
    buff[12*thr+8]  = h12_r;    buff[12*thr+9]  = h12_i;
    buff[12*thr+10] = h22_r;    buff[12*thr+11] = h22_i;
  }
  __syncthreads();
  if (thr < 32) 
    for (int i=0; i<6; i++) 
      myhess[(12*block+i)*SPLINE_BLOCK_SIZE+thr] = buff[i*SPLINE_BLOCK_SIZE+thr];

  __syncthreads();
  int th2 = thr-32;
  if (thr >= 32) {
    buff[12*th2+0]  = h00_r;    buff[12*th2+1]  = h00_i;
    buff[12*th2+2]  = h01_r;    buff[12*th2+3]  = h01_i;
    buff[12*th2+4]  = h02_r;    buff[12*th2+5]  = h02_i;
    buff[12*th2+6]  = h11_r;    buff[12*th2+7]  = h11_i;
    buff[12*th2+8]  = h12_r;    buff[12*th2+9]  = h12_i;
    buff[12*th2+10] = h22_r;    buff[12*th2+11] = h22_i;
  }
  __syncthreads();
  if (thr >= 32) {
    for (int i=0; i<6; i++) 
      myhess[(12*block+i+6)*SPLINE_BLOCK_SIZE+th2] = buff[i*SPLINE_BLOCK_SIZE+th2];
  }
}

#endif
