#include <stdarg.h>


// return the address of a (fortran) routine
// not sure of correct type "long" or "long long"?
#ifdef UNDERSCORE
unsigned long long int address_( unsigned long long int f )
#else
unsigned long long int address( unsigned long long int f )
#endif
{
//  printf("c address: %lld\n", f);
  return f;
}

#define args void  *i1, void *i2, void *i3, void  *i4, void *i5, void  *i6, void *i7, void *i8

#ifdef UNDERSCORE
void exe_(void (**f)(), args)
#else
void exe(void  (**f)(), args)
#endif
{
//  (**f)(i1,i2,i3,i4,i5,i6,i7,i8); this result in a segmentation fault on some machines
  (**f)();

}

#define args_address int *i[]

#ifdef UNDERSCORE
void exe_by_address_0_(void(**f)())
#else
void exe_by_address_0(void (**f)())
#endif
{   (**f)();
}

#ifdef UNDERSCORE
void exe_by_address_1_(void(**f)(), args_address)
#else
void exe_by_address_1(void (**f)(), args_address)
#endif
{   (**f)(i[0]);
}

#ifdef UNDERSCORE
void exe_by_address_2_(void(**f)(), args_address)
#else
void exe_by_address_2(void (**f)(), args_address)
#endif
 
{  (**f)(i[0],i[1]);
      }

#ifdef UNDERSCORE
void exe_by_address_3_(void(**f)(), args_address)
#else
void exe_by_address3(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2]);
      }

#ifdef UNDERSCORE
void exe_by_address_4_(void(**f)(), args_address)
#else
void exe_by_address_4(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2],i[3]);
      }

#ifdef UNDERSCORE
void exe_by_address_5_(void(**f)(), args_address)
#else
void exe_by_address_5(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2],i[3],i[4]);
      }

#ifdef UNDERSCORE
void exe_by_address_6_(void(**f)(), args_address)
#else
void exe_by_address_6(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2],i[3],i[4],i[5]);
      }

#ifdef UNDERSCORE
void exe_by_address_7_(void(**f)(), args_address)
#else
void exe_by_address_7(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2],i[3],i[4],i[5],i[6]);
      }

#ifdef UNDERSCORE
void exe_by_address_8_(void(**f)(), args_address)
#else
void exe_by_address_8(void (**f)(), args_address)
#endif
{   (**f)(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);
      }

// address to target
#ifdef UNDERSCORE
void address_to_target_integer_c_( int** p, int* array, int* dim1 )
#else
void address_to_target_integer_c( int** p, int* array, int* dim1 )
#endif
{
  int i;
  for(i=0; i< *dim1; i++)
  {
   array[i]=*(*p+i);
  }
  return;
}

#ifdef UNDERSCORE
void address_to_target_integer_1_c_( int** p, int* array, int* dim1 )
#else
void address_to_target_integer_1_c( int** p, int* array, int* dim1 )
#endif
{
  int i;
  for(i=0; i< *dim1; i++)
  {
   array[i]=*(*p+i);
  }
  return;
}

#ifdef UNDERSCORE
void address_to_target_double_c_( double** p, double* d )
#else
void address_to_target_double_c( double** p, double* d )
#endif
{
   *d = **p;
  return;
}

#ifdef UNDERSCORE
void address_to_target_double_1_c_( double** p, double* array, int* dim1 )
#else
void address_to_target_double_1_c( double** p, double* array, int* dim1 )
#endif
{
  int i;
  for(i=0; i< *dim1; i++)
  {
   array[i]=*(*p+i);
  }
  return;
}

//write in address
#ifdef UNDERSCORE
void write_in_address_double_c_( double** p, double* d)
#else
void write_in_address_double_c( double** p, double* d)
#endif
{
   **p = *d;
  return;
}

#ifdef UNDERSCORE
void write_in_address_double_1_c_( double** p, double* array, int* dim1 )
#else
void write_in_address_double_1_c( double** p, double* array, int* dim1 )
#endif
{
  int i;
  for(i=0; i< *dim1; i++)
  {
   *(*p+i) = array[i];
  }
  return;
}


