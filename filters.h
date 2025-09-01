/** @file filters.h
 *
 * @brief Header-only library of various digital filters in direct form.
 *
 * @par
 * Implementations of FIR, IIR direct form I, II, and II transposed, SOS form I,
 * and SOS form II transposed. Ring buffers are used for delay lines. For
 * simplicity, we assume the order of coefficients b and a are equal. You could
 * wrap the inputs in structs for cleaner instantiation.
 *
 * For more information on the theory behind FIR and IIR filters, see
 * Julius Smith, Introduction to Digital Filters with Audio Applications, 2007.
 *
 * These functions are intended to be easy to read and convenient for generic
 * prototyping on processors without SIMD. On supported platforms you
 * should try a vectorized library like CMSIS-DSP
 * (https://arm-software.github.io/CMSIS-DSP/latest/index.html).
 *
 * Todo: add a lattice all-pass filter
 *
 * @author ph919@ic.ac.uk
 * @date 2023
 */

#pragma once

/*!
 * @brief FIR Filter
 *
 * @param[in] x Input value
 * @param[in] N Order of coefficients
 * @param[in] b Pointer to coefficient array of size N
 * @param[out] idx Ring buffer index
 * @param[out] x_b Pointer to delay line buffer of size N - 1
 * @return Filtered value
 */
static inline double fir_filt(const double x, const unsigned int N,
                              const double *b, unsigned int *idx, double *x_b) {
  // Difference equation for filter order M = N - 1:
  // y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
  double y = b[0] * x;
  for (unsigned int k = 1; k < N; k++) {
    y += b[k] * x_b[(N - 1 + *idx - k) % (N - 1)];
  }
  // Update state
  x_b[*idx] = x;
  // Increment index
  *idx = (*idx + 1) % (N - 1);
  return y;
}

/*!
 * @brief Direct Form I IIR Filter
 *
 * @param[in] x Input value
 * @param[in] N Number of coefficients
 * @param[in] b Pointer to numerator coefficient array of size N
 * @param[in] a Pointer to denominator coefficient array of size N
 * @param[out] idx Ring buffer index
 * @param[out] x_b Pointer to input buffer of size N - 1
 * @param[out] y_b Pointer to output buffer of size N - 1
 * @return Filtered value
 */
static inline double iir_dfi_filt(const double x, const unsigned int N,
                                  const double *b, const double *a,
                                  unsigned int *idx, double *x_b, double *y_b) {
  // Difference equation for filter order M = N - 1:
  // a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
  //                       - a[1]*y[n-1] - ... - a[M]*y[n-M]
  double y = b[0] * x;
  for (unsigned int k = 1; k < N; k++) {
    y += b[k] * x_b[(N - 1 + *idx - k) % (N - 1)];
  }
  for (unsigned int k = 1; k < N; k++) {
    y -= a[k] * y_b[(N - 1 + *idx - k) % (N - 1)];
  }
  y /= a[0]; // Unnecessary if coefficients pre-normalized
  // Update state
  x_b[*idx] = x;
  y_b[*idx] = y;
  // Increment index
  *idx = (*idx + 1) % (N - 1);
  return y;
}

/*!
 * @brief Direct Form II IIR Filter
 *
 * @param[in] x Input value
 * @param[in] N Order of coefficients
 * @param[in] b Pointer to numerator coefficient array of size N
 * @param[in] a Pointer to denominator coefficient array of size N
 * @param[out] idx Ring buffer index
 * @param[out] x_b Pointer to input buffer of size N - 1
 * @param[out] y_b Pointer to output buffer of size N - 1
 * @return Filtered value
 */
static inline double iir_dfii_filt(const double x, const unsigned int N,
                                   const double *b, const double *a,
                                   unsigned int *idx, double *x_b,
                                   double *y_b) {
  // Difference equation for filter order M = N - 1:
  // a[0]*y[n] =             - a[1]*y[n-1] - ... - a[M]*y[n-M]
  //             + b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
  double y = 0;
  for (unsigned int k = 1; k < N; k++) {
    y -= a[k] * y_b[(N - 1 + *idx - k) % (N - 1)];
  }
  y += b[0] * x;
  for (unsigned int k = 1; k < N; k++) {
    y += b[k] * x_b[(N - 1 + *idx - k) % (N - 1)];
  }
  y /= a[0]; // Unnecessary if coefficients pre-normalized
  // Update state
  x_b[*idx] = x;
  y_b[*idx] = y;
  // Increment index
  *idx = (*idx + 1) % (N - 1);
  return y;
}

// /*!
//  * @brief Direct Form II Transposed IIR Filter
//  *
//  * @param[in] x Input value
//  * @param[in] N Number of coefficients
//  * @param[in] b Pointer to numerator coefficient array of size N
//  * @param[in] a Pointer to denominator coefficient array of size N
//  * @param[out] idx Ring buffer index
//  * @param[out] x_b Pointer to input buffer of size N - 1
//  * @param[out] y_b Pointer to output buffer of size N - 1
//  * @return Filtered value
//  */
// static inline double iir_dfiit_filt(const double x, const unsigned int N,
//                                     const double *b, const double *a,
//                                     unsigned int *idx, double *x_b,
//                                     double *y_b) {
//   // Difference equation for filter order M = N - 1:
//   // a[0]*y[n] = b[0]*x[n] +
//   //             b[1]*x[n-1] - a[1]*y[n-1] - ... + b[M]*x[n-M] - a[M]*y[n-M]
//   double y = b[0] * x;
//   for (unsigned int k = 1; k < N; k++) {
//     y += b[k] * x_b[(N - 1 + *idx - k) % (N - 1)] -
//          a[k] * y_b[(N - 1 + *idx - k) % (N - 1)];
//   }
//   y /= a[0]; // Unnecessary if coefficients pre-normalized
//   // Update state
//   x_b[*idx] = x;
//   y_b[*idx] = y;
//   // Increment index
//   *idx = (*idx + 1) % (N - 1);
//   return y;
// }

/*!
 * @brief Direct Form II Transposed IIR Filter
 *
 * @param[in] x Input value
 * @param[in] N Number of coefficients
 * @param[in] b Pointer to numerator coefficient array of size N
 * @param[in] a Pointer to denominator coefficient array of size N
 * @param[out] s Pointer to state array of size N - 1
 * @return Filtered value
 *
 * @par
 * Obtain filter coefficients easily in python with scipy.signal.butter:
 * print(scipy.signal.butter(N=N, Wn=Wn, output='ba'))
 */
static inline double iir_dfiit_filt(const double x, const unsigned int N,
                                    const double *b, const double *a,
                                    double *s) {
  // Difference equations for filter order M = N - 1:
  // a[0]*y[n]   = s[0][n-1]   + b[0]  *x[n]
  // s[0][n]     = s[1][n-1]   + b[1]  *x[n] - a[1]  *y[n]
  // s[1][n]     = s[2][n-1]   + b[2]  *x[n] - a[2]  *y[n]
  // ...
  // s[M-2][n]   = s[M-1][n-1] + b[M-1]*x[n] - a[M-1]*y[n]
  // s[M-1][n]   =               b[M]  *x[n] - a[M]  *y[n]

  // Division by a[0] unnecessary if coefficients pre-normalized
  double y = (b[0] * x + s[0]) / a[0];
  for (unsigned int k = 0; k < N - 2; k++) {
    s[k] = s[k + 1] + b[k + 1] * x - a[k + 1] * y;
  }
  s[N - 2] = b[N - 1] * x - a[N - 1] * y;
  return y;
}

// /*!
//  * @brief Cascaded Biquad Direct Form I IIR Filter
//  *
//  * @param[in] x Input value
//  * @param[in] Nstages Number of biquad stages
//  * @param[in] sos Pointer to coefficient array of size Nstages * 6
//  * @param[out] idx_b Pointer to ring buffer index array
//  * @param[out] s Pointer to state array of size Nstages * 4
//  * @return Filtered value
//  */
// static inline double sos_dfi_filt(const double x, const unsigned int Nstages,
//                                   const double *sos, unsigned int *idx_b,
//                                   double *s) {
//   double y = x;
//   for (unsigned int k = 0; k < Nstages; k++) {
//     y = iir_dfi_filt(y, 3, &sos[k * 6], &sos[3 + k * 6], &idx_b[k], &s[k *
//     4],
//                      &s[2 + k * 4]);
//   }
//   return y;
// }

/*!
 * @brief Cascaded Biquad Unrolled Direct Form I IIR Filter
 *
 * @param[in] x Input value
 * @param[in] Nstages Number of biquad stages
 * @param[in] sos Pointer to coefficient array of size Nstages * 6
 * @param[out] s Pointer to state array of size Nstages * 4
 * @return Filtered value
 *
 * @par
 * This version of sos_dfi_filt() unrolls
 * (https://en.wikipedia.org/wiki/Loop_unrolling) the second order filter
 * and updates the entire state manually, removing the need to keep track of
 * every ring buffer index.
 *
 * The format of sos is similar to the format used by scipy, but instead of an
 * array of shape (Nstages, 6), we use a single array of size Nstages * 6. Every
 * hextuple represents the coefficients of one second-order section, with the
 * first three values being the numerator coefficients and the last three being
 * the denominator coefficients (effectively just numpy.ndarray.flatten).
 * Example python output:
 * print(scipy.signal.butter(N=4, Wn=Wn, output='sos'))
 * > [[b[0][0] b[0][1] b[0][2] a[0][0] a[0][1] a[0][2]
 * >  [b[1][0] b[1][1] b[1][2] a[1][0] a[1][1] a[1][2]]
 * Corresponding sos definition:
 * const double sos[2 * 6] = {b[0][0], b[0][1], b[0][2], a[0][0], a[0][1],
 *                            a[0][2], b[1][0], b[1][1], b[1][2], a[1][0],
 *                            a[1][1], a[1][2]}
 */
static inline double sos_dfi_filt(const double x, const unsigned int Nstages,
                                  const double *sos, double *s) {
  double y = x;
  for (unsigned int k = 0; k < Nstages; k++) {
    const double *b = &sos[k * 6];
    const double *a = &sos[3 + k * 6];
    double *x_b = &s[k * 4];
    double *y_b = &s[2 + k * 4];
    // Difference equation for filter order 2 = 3 - 1:
    // a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + b[2]*x[n-2]
    //                       - a[1]*y[n-1] - a[2]*y[n-2]
    double acc = b[0] * y;
    acc += b[1] * x_b[0];
    acc += b[2] * x_b[1];
    acc -= a[1] * y_b[0];
    acc -= a[2] * y_b[1];
    acc /= a[0]; // Unnecessary if coefficients pre-normalized
    // Update state
    x_b[1] = x_b[0];
    x_b[0] = y;
    y_b[1] = y_b[0];
    y_b[0] = acc;
    y = y_b[0];
  }
  return y;
}

/*!
 * @brief Cascaded Biquad Direct Form II Transposed IIR Filter
 *
 * @param[in] x Input value
 * @param[in] Nstages Number of biquad stages
 * @param[in] sos Pointer to coefficient array of size Nstages * 6
 * @param[out] s Pointer to state array of size Nstages * 2
 * @return Filtered value
 *
 * @par
 * Similar to sos_dfi_filt(), except using direct form II transposed IIR biquads
 */
static inline double sos_dfiit_filt(const double x, const unsigned int Nstages,
                                    const double *sos, double *s) {
  double y = x;
  for (unsigned int k = 0; k < Nstages; k++) {
    y = iir_dfiit_filt(y, 3, &sos[k * 6], &sos[3 + k * 6], &s[k * 2]);
  }
  return y;
}