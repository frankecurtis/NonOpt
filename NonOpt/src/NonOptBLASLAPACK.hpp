// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTBLASLAPACK_HPP__
#define __NONOPTBLASLAPACK_HPP__

extern "C"
{

  // Subroutines
  void daxpy_(int* n, double* a, double* x, int* incx, double* y, int* incy);
  void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
  void dscal_(int* n, double* a, double* x, int* incx);
  void dsymv_(char* u, int* n, double* a, double* A, int* m, double* x, int* incx, double* b, double* y, int* incy);
  void dsyr_(char* u, int* n, double* a, double* x, int* incx, double* A, int* m);
  void dsyr2_(char* u, int* n, double* a, double* x, int* incx, double* y, int* incy, double* A, int* m);
  void dsysv_(char* u, int* n, int* r, double* a, int* la, int* piv, double* b, int* lb, double* work, int* lw, int* info);
  void dsytrf_(char* u, int* n, double* a, int* m, int* o, double* w, int* l, int* f);
  void dsytrs_(char* u, int* n, int* r, double* a, int* m, int* o, double* b, int* d, int* f);

  void dtrsv_(char* u, char* t, char* d, int* n, double* A, int* m, double* x, int* incx);
  void dpotf2_(char* u, int* n, double* A, int* m, int* o);

  // Scalar functions
  double dasum_(int* n, double* x, int* incx);
  double ddot_(int* n, double* x, int* incx, double* y, int* incy);
  double dnrm2_(int* n, double* x, int* incx);
  int idamax_(int* n, double* x, int* incx);

} // end extern "C"

#endif /* __NONOPTBLASLAPACK_HPP__ */
