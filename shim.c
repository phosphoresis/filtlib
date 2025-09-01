/** @file shim.c
 * @brief Shim to shared object for Python FFI.
 * @author ph919@ic.ac.uk
 * @date 2023
 */

#include "filters.h"

double s_fir_filt(const double x, const unsigned int N, const double *b,
                  unsigned int *idx, double *x_b) {
  return fir_filt(x, N, b, idx, x_b);
}

double s_iir_dfi_filt(const double x, const unsigned int N, const double *b,
                      const double *a, unsigned int *idx, double *x_b,
                      double *y_b) {
  return iir_dfi_filt(x, N, b, a, idx, x_b, y_b);
}

double s_iir_dfii_filt(const double x, const unsigned int N, const double *b,
                       const double *a, unsigned int *idx, double *x_b,
                       double *y_b) {
  return iir_dfii_filt(x, N, b, a, idx, x_b, y_b);
}

double s_iir_dfiit_filt(const double x, const unsigned int N, const double *b,
                        const double *a, double *s) {
  return iir_dfiit_filt(x, N, b, a, s);
}

double s_sos_dfi_filt(const double x, const unsigned int Nstages,
                      const double *sos, double *s) {
  return sos_dfi_filt(x, Nstages, sos, s);
}

double s_sos_dfiit_filt(const double x, const unsigned int Nstages,
                        const double *sos, double *s) {
  return sos_dfiit_filt(x, Nstages, sos, s);
}