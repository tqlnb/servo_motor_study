/*
 * File: sumAll.c
 *
 * Code generated for Simulink model 'sumAll'.
 *
 * Model version                  : 1.6
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Thu Oct 10 21:22:43 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "sumAll.h"
#include "rtwtypes.h"
#include <string.h>
#include <emmintrin.h>
#include <float.h>
#include "rt_look.h"
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#define NumBitsPerChar                 8U
#ifndef CodeFormat
#define CodeFormat                     S-Function
#else
#undef CodeFormat
#define CodeFormat                     S-Function
#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#else
#undef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#endif

#ifndef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#else
#undef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#endif

#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

#if !defined(RTW_SFUNCTION_DEFINES)
#define RTW_SFUNCTION_DEFINES
#ifndef _RTW_COMMON_DEFINES_
#define _RTW_COMMON_DEFINES_
#endif
#endif

const real_T sumAll_RGND = 0.0;        /* real_T ground */

/* Block signals and states (default storage) */
DW rtDW;

/* Real-time model */
static RT_MODEL rtM_;
RT_MODEL *const rtM = &rtM_;
extern void rt_invd4x4_snf(const real_T u[16], real_T y[16]);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
static void rate_scheduler(void);
static real_T rtGetNaN(void);
static real32_T rtGetNaNF(void);

/*===========*
 * Constants *
 *===========*/
#define RT_PI                          3.14159265358979323846
#define RT_PIF                         3.1415927F
#define RT_LN_10                       2.30258509299404568402
#define RT_LN_10F                      2.3025851F
#define RT_LOG10E                      0.43429448190325182765
#define RT_LOG10EF                     0.43429449F
#define RT_E                           2.7182818284590452354
#define RT_EF                          2.7182817F

/*
 * UNUSED_PARAMETER(x)
 *   Used to specify that a function parameter (argument) is required but not
 *   accessed by the function body.
 */
#ifndef UNUSED_PARAMETER
#if defined(__LCC__)
#define UNUSED_PARAMETER(x)                                      /* do nothing */
#else

/*
 * This is the semi-ANSI standard way of indicating that an
 * unused function parameter is required.
 */
#define UNUSED_PARAMETER(x)            (void) (x)
#endif
#endif

#ifndef INTERP
#define INTERP(x,x1,x2,y1,y2)          ( (y1)+(((y2) - (y1))/((x2) - (x1)))*((x)-(x1)) )
#endif

#ifndef ZEROTECHNIQUE
#define ZEROTECHNIQUE

typedef enum {
  NORMAL_INTERP,
  AVERAGE_VALUE,
  MIDDLE_VALUE
} ZeroTechnique;

#endif

static int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u) ;
static real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y);

#define NOT_USING_NONFINITE_LITERALS   1

extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
static void rt_InitInfAndNaN(size_t realSize);
static boolean_T rtIsInf(real_T value);
static boolean_T rtIsInfF(real32_T value);
static boolean_T rtIsNaN(real_T value);
static boolean_T rtIsNaNF(real32_T value);
typedef struct {
  struct {
    uint32_T wordH;
    uint32_T wordL;
  } words;
} BigEndianIEEEDouble;

typedef struct {
  struct {
    uint32_T wordL;
    uint32_T wordH;
  } words;
} LittleEndianIEEEDouble;

typedef struct {
  union {
    real32_T wordLreal;
    uint32_T wordLuint;
  } wordL;
} IEEESingle;

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;
static real_T rtGetInf(void);
static real32_T rtGetInfF(void);
static real_T rtGetMinusInf(void);
static real32_T rtGetMinusInfF(void);

/*
 * Initialize rtNaN needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetNaN(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T nan = 0.0;
  if (bitsPerReal == 32U) {
    nan = rtGetNaNF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0xFFF80000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    nan = tmpVal.fltVal;
  }

  return nan;
}

/*
 * Initialize rtNaNF needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetNaNF(void)
{
  IEEESingle nanF = { { 0.0F } };

  nanF.wordL.wordLuint = 0xFFC00000U;
  return nanF.wordL.wordLreal;
}

/*
 * Routine to get the index of the input from a table using binary or
 * interpolation search.
 *
 * Inputs:
 * *x   : Pointer to table, x[0] ....x[xlen-1]
 * xlen : Number of values in xtable
 * u    : input value to look up
 *
 * Output:
 * idx  : the index into the table such that:
 * if u is negative
 * x[idx] <= u < x[idx+1]
 * else
 * x[idx] < u <= x[idx+1]
 *
 * Interpolation Search: If the table contains a large number of nearly
 * uniformly spaced entries, i.e., x[n] vs n is linear then the index
 * corresponding to the input can be found in one shot using the linear
 * interpolation formula. Therefore if you have a look-up table block with
 * many data points, using interpolation search might speed up the code.
 * Compile the generated code with the following flag:
 *
 * make_rtw OPTS=-DDOINTERPSEARCH
 *
 * to enable interpolation search.
 */
static int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u)
{
  int_T idx = 0;
  int_T bottom = 0;
  int_T top = xlen-1;
  int_T retValue = 0;
  boolean_T returnStatus = 0U;

#ifdef DOINTERPSEARCH

  real_T offset = 0;

#endif

  /*
   * Deal with the extreme cases first:
   *   if u <= x[bottom] then return idx = bottom
   *   if u >= x[top]    then return idx = top-1
   */
  if (u <= x[bottom]) {
    retValue = bottom;
    returnStatus = 1U;
  } else if (u >= x[top]) {
    retValue = top-1;
    returnStatus = 1U;
  } else {
    /* else required to ensure safe programming, even *
     * if it's expected that it will never be reached */
  }

  if (returnStatus == 0U) {
    if (u < 0) {
      /* For negative input find index such that: x[idx] <= u < x[idx+1] */
      for (;;) {

#ifdef DOINTERPSEARCH

        offset = (u-x[bottom])/(x[top]-x[bottom]);
        idx = bottom + (int_T)((top-bottom)*(offset-DBL_EPSILON));

#else

        idx = (bottom + top)/2;

#endif

        if (u < x[idx]) {
          top = idx - 1;
        } else if (u >= x[idx+1]) {
          bottom = idx + 1;
        } else {
          /* we have x[idx] <= u < x[idx+1], return idx */
          retValue = idx;
          break;
        }
      }
    } else {
      /* For non-negative input find index such that: x[idx] < u <= x[idx+1] */
      for (;;) {

#ifdef DOINTERPSEARCH

        offset = (u-x[bottom])/(x[top]-x[bottom]);
        idx = bottom + (int_T)((top-bottom)*(offset-DBL_EPSILON));

#else

        idx = (bottom + top)/2;

#endif

        if (u <= x[idx]) {
          top = idx - 1;
        } else if (u > x[idx+1]) {
          bottom = idx + 1;
        } else {
          /* we have x[idx] < u <= x[idx+1], return idx */
          retValue = idx;
          break;
        }
      }
    }
  }

  return retValue;
}

/* 1D lookup routine for data type of real_T. */
static real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y)
{
  int_T idx = rt_GetLookupIndex(x, xlen, u);
  real_T num = y[idx+1] - y[idx];
  real_T den = x[idx+1] - x[idx];

  /* Due to the way the binary search is implemented
     in rt_look.c (rt_GetLookupIndex), den cannot be
     0.  Equivalently, m cannot be inf or nan. */
  real_T m = num/den;
  return (y[idx] + (m * (u - x[idx])));
}

/*
 * Initialize the rtInf, rtMinusInf, and rtNaN needed by the
 * generated code. NaN is initialized as non-signaling. Assumes IEEE.
 */
static void rt_InitInfAndNaN(size_t realSize)
{
  (void) (realSize);
  rtNaN = rtGetNaN();
  rtNaNF = rtGetNaNF();
  rtInf = rtGetInf();
  rtInfF = rtGetInfF();
  rtMinusInf = rtGetMinusInf();
  rtMinusInfF = rtGetMinusInfF();
}

/* Test if value is infinite */
static boolean_T rtIsInf(real_T value)
{
  return (boolean_T)((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

/* Test if single-precision value is infinite */
static boolean_T rtIsInfF(real32_T value)
{
  return (boolean_T)(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

/* Test if value is not a number */
static boolean_T rtIsNaN(real_T value)
{
  boolean_T result = (boolean_T) 0;
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  if (bitsPerReal == 32U) {
    result = rtIsNaNF((real32_T)value);
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.fltVal = value;
    result = (boolean_T)((tmpVal.bitVal.words.wordH & 0x7FF00000) == 0x7FF00000 &&
                         ( (tmpVal.bitVal.words.wordH & 0x000FFFFF) != 0 ||
                          (tmpVal.bitVal.words.wordL != 0) ));
  }

  return result;
}

/* Test if single-precision value is not a number */
static boolean_T rtIsNaNF(real32_T value)
{
  IEEESingle tmp;
  tmp.wordL.wordLreal = value;
  return (boolean_T)( (tmp.wordL.wordLuint & 0x7F800000) == 0x7F800000 &&
                     (tmp.wordL.wordLuint & 0x007FFFFF) != 0 );
}

/*
 * Initialize rtInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetInf(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T inf = 0.0;
  if (bitsPerReal == 32U) {
    inf = rtGetInfF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0x7FF00000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    inf = tmpVal.fltVal;
  }

  return inf;
}

/*
 * Initialize rtInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetInfF(void)
{
  IEEESingle infF;
  infF.wordL.wordLuint = 0x7F800000U;
  return infF.wordL.wordLreal;
}

/*
 * Initialize rtMinusInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetMinusInf(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T minf = 0.0;
  if (bitsPerReal == 32U) {
    minf = rtGetMinusInfF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0xFFF00000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    minf = tmpVal.fltVal;
  }

  return minf;
}

/*
 * Initialize rtMinusInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetMinusInfF(void)
{
  IEEESingle minfF;
  minfF.wordL.wordLuint = 0xFF800000U;
  return minfF.wordL.wordLreal;
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (rtM->Timing.TaskCounters.TID[1])++;
  if ((rtM->Timing.TaskCounters.TID[1]) > 99999) {/* Sample time: [1.0s, 0.0s] */
    rtM->Timing.TaskCounters.TID[1] = 0;
  }
}

void rt_invd4x4_snf(const real_T u[16], real_T y[16])
{
  real_T x[16];
  real_T smax;
  int32_T ijA;
  int32_T ipiv_tmp;
  int32_T ix;
  int32_T iy;
  int32_T jA;
  int32_T jBcol;
  int32_T jj;
  int32_T jp1j;
  int32_T pipk;
  int8_T ipiv[4];
  int8_T p[4];
  for (ipiv_tmp = 0; ipiv_tmp < 16; ipiv_tmp++) {
    y[ipiv_tmp] = 0.0;
    x[ipiv_tmp] = u[ipiv_tmp];
  }

  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (pipk = 0; pipk < 3; pipk++) {
    jBcol = pipk * 5 + 2;
    jj = pipk * 5;
    ix = 4 - pipk;
    iy = 1;
    smax = fabs(x[jj]);
    for (jA = 2; jA <= ix; jA++) {
      real_T s;
      s = fabs(x[(jBcol + jA) - 3]);
      if (s > smax) {
        iy = jA;
        smax = s;
      }
    }

    if (x[(jBcol + iy) - 3] != 0.0) {
      if (iy - 1 != 0) {
        ipiv_tmp = pipk + iy;
        ipiv[pipk] = (int8_T)ipiv_tmp;
        smax = x[pipk];
        x[pipk] = x[ipiv_tmp - 1];
        x[ipiv_tmp - 1] = smax;
        smax = x[pipk + 4];
        x[pipk + 4] = x[ipiv_tmp + 3];
        x[ipiv_tmp + 3] = smax;
        smax = x[pipk + 8];
        x[pipk + 8] = x[ipiv_tmp + 7];
        x[ipiv_tmp + 7] = smax;
        smax = x[pipk + 12];
        x[pipk + 12] = x[ipiv_tmp + 11];
        x[ipiv_tmp + 11] = smax;
      }

      iy = jBcol - pipk;
      for (ix = jBcol; ix <= iy + 2; ix++) {
        x[ix - 1] /= x[jj];
      }
    }

    ix = 3 - pipk;
    jA = jj;
    jj += 4;
    for (jp1j = 0; jp1j < ix; jp1j++) {
      smax = x[(jp1j << 2) + jj];
      if (smax != 0.0) {
        iy = jA + 6;
        ipiv_tmp = jA - pipk;
        for (ijA = iy; ijA <= ipiv_tmp + 8; ijA++) {
          x[ijA - 1] += x[((jBcol + ijA) - jA) - 7] * -smax;
        }
      }

      jA += 4;
    }
  }

  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 4;
  if (ipiv[0] > 1) {
    pipk = p[ipiv[0] - 1];
    p[ipiv[0] - 1] = 1;
    p[0] = (int8_T)pipk;
  }

  if (ipiv[1] > 2) {
    pipk = p[ipiv[1] - 1];
    p[ipiv[1] - 1] = p[1];
    p[1] = (int8_T)pipk;
  }

  if (ipiv[2] > 3) {
    pipk = p[ipiv[2] - 1];
    p[ipiv[2] - 1] = p[2];
    p[2] = (int8_T)pipk;
  }

  for (jA = 0; jA < 4; jA++) {
    jj = (p[jA] - 1) << 2;
    y[jA + jj] = 1.0;
    for (pipk = jA + 1; pipk < 5; pipk++) {
      ipiv_tmp = (jj + pipk) - 1;
      if (y[ipiv_tmp] != 0.0) {
        for (ix = pipk + 1; ix < 5; ix++) {
          jBcol = (jj + ix) - 1;
          y[jBcol] -= x[(((pipk - 1) << 2) + ix) - 1] * y[ipiv_tmp];
        }
      }
    }
  }

  for (pipk = 0; pipk < 4; pipk++) {
    jBcol = pipk << 2;
    for (jA = 3; jA >= 0; jA--) {
      jp1j = jA << 2;
      ipiv_tmp = jA + jBcol;
      smax = y[ipiv_tmp];
      if (smax != 0.0) {
        y[ipiv_tmp] = smax / x[jA + jp1j];
        iy = jA - 1;
        for (ix = 0; ix <= iy; ix++) {
          jj = ix + jBcol;
          y[jj] -= x[ix + jp1j] * y[ipiv_tmp];
        }
      }
    }
  }
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int32_T tmp;
    int32_T tmp_0;
    if (u0 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u1 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = atan2(tmp, tmp_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != trunc(u1))) {
    real_T q;
    q = fabs(u0 / u1);
    if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
  }

  return y;
}

/* Model step function */
void sumAll_step(void)
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  __m128d tmp_6;
  __m128d tmp_7;
  real_T W43wr[16];
  real_T rtb_Lminrows24col24_0[16];
  real_T rtb_Sum5[16];
  real_T rtb_MultiportSwitch[4];
  real_T rtb_Sum2[4];
  real_T rtb_MathFunction[3];
  real_T TrigonometricFunction_o1;
  real_T TrigonometricFunction_o2;
  real_T isb_g;
  real_T rtb_Gain2_l;
  real_T rtb_Mult1_idx_0;
  real_T rtb_Mult1_idx_1;
  real_T rtb_MultiportSwitch1_f;
  real_T rtb_MultiportSwitch1_idx_2;
  real_T rtb_MultiportSwitch1_idx_3;
  real_T rtb_Phisat;
  real_T rtb_Saturation;
  real_T rtb_Switch_m_idx_0;
  real_T rtb_phimd;
  real_T rtb_u_Vb_idx_3;
  real_T rtb_unitconversion_idx_11;
  real_T rtb_unitconversion_idx_6;
  real_T rtb_unitconversion_idx_9;
  real_T vdr_c;
  real_T vds_m;
  int32_T i;
  int32_T i_0;
  int32_T rtb_Sum5_tmp;
  uint8_T rtb_Subtract2_m;
  boolean_T rtb_LogicalOperator2;

  /* Update for UnitDelay: '<S43>/fluxes' incorporates:
   *  UnitDelay: '<S45>/fluxes'
   */
  rtDW.fluxes_DSTATE_m[0] = rtDW.fluxes_DSTATE[0];
  rtDW.fluxes_DSTATE_m[1] = rtDW.fluxes_DSTATE[1];
  rtDW.fluxes_DSTATE_m[2] = rtDW.fluxes_DSTATE[2];
  rtDW.fluxes_DSTATE_m[3] = rtDW.fluxes_DSTATE[3];

  /* Switch: '<S38>/Switch' incorporates:
   *  Constant: '<S38>/Constant2'
   *  Product: '<S44>/inversion'
   */
  memcpy(&rtb_Lminrows24col24_0[0], &rtConstP.Constant2_Value[0], sizeof(real_T)
         << 4U);

  /* UnitDelay: '<S45>/fluxes' incorporates:
   *  Product: '<S38>/Product3'
   */
  rtb_Mult1_idx_1 = rtDW.fluxes_DSTATE[1];
  rtb_Gain2_l = rtDW.fluxes_DSTATE[0];
  rtb_Mult1_idx_0 = rtDW.fluxes_DSTATE[2];
  rtb_Saturation = rtDW.fluxes_DSTATE[3];
  for (i = 0; i <= 2; i += 2) {
    /* Product: '<S38>/Product3' incorporates:
     *  UnitDelay: '<S45>/fluxes'
     */
    tmp_4 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 4]);
    tmp_5 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i]);
    tmp_6 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 8]);
    tmp_7 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 12]);
    _mm_storeu_pd(&rtb_Sum2[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_4, _mm_set1_pd(rtb_Mult1_idx_1)), _mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_Gain2_l))), _mm_mul_pd(tmp_6, _mm_set1_pd(rtb_Mult1_idx_0))),
      _mm_mul_pd(tmp_7, _mm_set1_pd(rtb_Saturation))));
  }

  /* UnitDelay: '<S60>/wm_delay' */
  rtb_Phisat = rtDW.wm_delay_DSTATE;

  /* Sum: '<S60>/Sum1' incorporates:
   *  Gain: '<S60>/F2'
   *  UnitDelay: '<S60>/wm_predict'
   */
  rtb_phimd = 2.0 * rtb_Phisat - rtDW.wm_predict_DSTATE;

  /* Outputs for Enabled SubSystem: '<S41>/sin(thr),cos(thr)' incorporates:
   *  EnablePort: '<S58>/Enable'
   */
  if (rtDW.sinthrcosthr_MODE) {
    /* Disable for Trigonometry: '<S58>/Trigonometric Function' incorporates:
     *  Outport: '<S58>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o1_a = 0.0;

    /* Disable for Trigonometry: '<S58>/Trigonometric Function' incorporates:
     *  Outport: '<S58>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o2_a = 0.0;

    /* Disable for Outport: '<S58>/sin(thr),cos(thr)' incorporates:
     *  Constant: '<S58>/Constant'
     */
    rtDW.Constant_n[0] = 0.0;
    rtDW.Constant_n[1] = 0.0;

    /* Disable for Assignment: '<S58>/W(2,1)=-wr' incorporates:
     *  Outport: '<S58>/W'
     */
    memset(&rtDW.W21wr[0], 0, sizeof(real_T) << 4U);
    rtDW.sinthrcosthr_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S41>/sin(thr),cos(thr)' */

  /* Outputs for Enabled SubSystem: '<S41>/sin(thr),cos(thr)1' incorporates:
   *  EnablePort: '<S59>/Enable'
   */
  /* Assignment: '<S59>/W(4,3)=wr' incorporates:
   *  SignalConversion generated from: '<S59>/W(3,4)=-wr'
   */
  memset(&W43wr[0], 0, sizeof(real_T) << 4U);

  /* Gain: '<S59>/Gain3' incorporates:
   *  Assignment: '<S59>/W(4,3)=wr'
   */
  W43wr[14] = -rtb_phimd;

  /* Assignment: '<S59>/W(4,3)=wr' */
  W43wr[11] = rtb_phimd;

  /* Trigonometry: '<S59>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S37>/Rotor angle thetam'
   */
  TrigonometricFunction_o1 = sin(rtDW.Rotoranglethetam_DSTATE);

  /* Trigonometry: '<S59>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S37>/Rotor angle thetam'
   */
  TrigonometricFunction_o2 = cos(rtDW.Rotoranglethetam_DSTATE);

  /* End of Outputs for SubSystem: '<S41>/sin(thr),cos(thr)1' */

  /* Outputs for Enabled SubSystem: '<S40>/Rotor reference frame' incorporates:
   *  EnablePort: '<S54>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE) {
    /* Disable for Fcn: '<S54>/ira' incorporates:
     *  Outport: '<S54>/ira,irb'
     */
    rtDW.ira_d = 0.0;

    /* Disable for Fcn: '<S54>/irb' incorporates:
     *  Outport: '<S54>/ira,irb'
     */
    rtDW.irb_f = 0.0;

    /* Disable for Fcn: '<S54>/isa' incorporates:
     *  Outport: '<S54>/isa,isb'
     */
    rtDW.isa_p = 0.0;

    /* Disable for Fcn: '<S54>/isb' incorporates:
     *  Outport: '<S54>/isa,isb'
     */
    rtDW.isb_j = 0.0;
    rtDW.Rotorreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S40>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S40>/Stationary reference frame' incorporates:
   *  EnablePort: '<S55>/Enable'
   */
  /* Fcn: '<S55>/isb' */
  isb_g = -(1.7320508075688772 * rtb_Sum2[1] + rtb_Sum2[0]) / 2.0;

  /* End of Outputs for SubSystem: '<S40>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S40>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S56>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MODE) {
    /* Disable for Fcn: '<S56>/ira' incorporates:
     *  Outport: '<S56>/ira,irb'
     */
    rtDW.ira = 0.0;

    /* Disable for Fcn: '<S56>/irb' incorporates:
     *  Outport: '<S56>/ira,irb'
     */
    rtDW.irb = 0.0;

    /* Disable for Fcn: '<S56>/isa' incorporates:
     *  Outport: '<S56>/isa,isb'
     */
    rtDW.isa = 0.0;

    /* Disable for Fcn: '<S56>/isb' incorporates:
     *  Outport: '<S56>/isa,isb'
     */
    rtDW.isb = 0.0;
    rtDW.Synchronousreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S40>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S40>/Stationary reference frame' incorporates:
   *  EnablePort: '<S55>/Enable'
   */
  /* Gain: '<S40>/ib' incorporates:
   *  Fcn: '<S55>/ira'
   *  Fcn: '<S55>/irb'
   *  Fcn: '<S55>/isa'
   *  MultiPortSwitch: '<S40>/Multiport Switch1'
   *  MultiPortSwitch: '<S41>/Multiport Switch'
   */
  rtDW.ib[0] = (TrigonometricFunction_o2 * rtb_Sum2[2] -
                TrigonometricFunction_o1 * rtb_Sum2[3]) * 13.731987951966302;
  rtDW.ib[2] = 13.731987951966302 * rtb_Sum2[0];
  rtDW.ib[1] = ((-TrigonometricFunction_o2 - 1.7320508075688772 *
                 TrigonometricFunction_o1) * rtb_Sum2[2] +
                (TrigonometricFunction_o1 - 1.7320508075688772 *
                 TrigonometricFunction_o2) * rtb_Sum2[3]) / 2.0 *
    13.731987951966302;

  /* End of Outputs for SubSystem: '<S40>/Stationary reference frame' */
  rtDW.ib[3] = 13.731987951966302 * isb_g;

  /* S-Function (sfun_spssw_discc): '<S64>/State-Space' incorporates:
   *  Constant: '<S61>/DC'
   */

  /* S-Function block: <S64>/State-Space */
  {
    real_T accum;

    /* Circuit has switches */
    int_T *switch_status = (int_T*) rtDW.StateSpace_PWORK.SWITCH_STATUS;
    int_T *switch_status_init = (int_T*)
      rtDW.StateSpace_PWORK.SWITCH_STATUS_INIT;
    int_T *SwitchChange = (int_T*) rtDW.StateSpace_PWORK.SW_CHG;
    int_T *gState = (int_T*) rtDW.StateSpace_PWORK.G_STATE;
    real_T *yswitch = (real_T*)rtDW.StateSpace_PWORK.Y_SWITCH;
    int_T *switchTypes = (int_T*) rtDW.StateSpace_PWORK.SWITCH_TYPES;
    int_T *idxOutSw = (int_T*) rtDW.StateSpace_PWORK.IDX_OUT_SW;
    real_T *DxCol = (real_T*)rtDW.StateSpace_PWORK.DX_COL;
    real_T *tmp2 = (real_T*)rtDW.StateSpace_PWORK.TMP2;
    real_T *uswlast = (real_T*)rtDW.StateSpace_PWORK.USWLAST;
    int_T newState;
    int_T swChanged = 0;
    int loopsToDo = 20;
    real_T temp;

    /* keep an initial copy of switch_status*/
    memcpy(switch_status_init, switch_status, 6 * sizeof(int_T));
    memcpy(uswlast, &rtDW.StateSpace_o1[0], 6*sizeof(real_T));
    do {
      if (loopsToDo == 1) {            /* Need to reset some variables: */
        swChanged = 0;

        /* return to the original switch status*/
        {
          int_T i1;
          for (i1=0; i1 < 6; i1++) {
            swChanged = ((SwitchChange[i1] = switch_status_init[i1] -
                          switch_status[i1]) != 0) ? 1 : swChanged;
            switch_status[i1] = switch_status_init[i1];
          }
        }
      } else {
        /*
         * Compute outputs:
         * ---------------
         */
        real_T *Ds = (real_T*)rtDW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &rtDW.StateSpace_o1[0];
          for (i1=0; i1 < 8; i1++) {
            accum = 0.0;

            {
              int_T i2;
              const real_T *u0;
              for (i2=0; i2 < 6; i2++) {
                accum += *(Ds++) * 0.0;
              }

              accum += *(Ds++) * rtDW.ib[2];
              accum += *(Ds++) * rtDW.ib[3];
              accum += *(Ds++) * 500.0;
            }

            y0[i1] = accum;
          }
        }

        swChanged = 0;

        {
          int_T i1;
          real_T *y0 = &rtDW.StateSpace_o1[0];
          for (i1=0; i1 < 6; i1++) {
            newState = ((y0[i1] > 0.0) && (gState[i1] > 0)) || (y0[i1] < 0.0) ?
              1 : (((y0[i1] > 0.0) && gState[i1] == 0) ? 0 : switch_status[i1]);
            swChanged = ((SwitchChange[i1] = newState - switch_status[i1]) != 0)
              ? 1 : swChanged;
            switch_status[i1] = newState;/* Keep new state */
          }
        }
      }

      /*
       * Compute new As, Bs, Cs and Ds matrixes:
       * --------------------------------------
       */
      if (swChanged) {
        real_T *Ds = (real_T*)rtDW.StateSpace_PWORK.DS;
        real_T a1;

        {
          int_T i1;
          for (i1=0; i1 < 6; i1++) {
            if (SwitchChange[i1] != 0) {
              a1 = 1000.0*SwitchChange[i1];
              temp = 1/(1-Ds[i1*10]*a1);

              {
                int_T i2;
                for (i2=0; i2 < 8; i2++) {
                  DxCol[i2]= Ds[i2 * 9 + i1]*temp*a1;
                }
              }

              DxCol[i1] = temp;

              /* Copy row nSw of Ds into tmp2 and zero it out in Ds */
              memcpy(tmp2, &Ds[i1 * 9], 9 * sizeof(real_T));
              memset(&Ds[i1 * 9], '\0', 9 * sizeof(real_T));

              /* Cs = Cs + DxCol * tmp1, Ds = Ds + DxCol * tmp2 *******************/
              {
                int_T i2;
                for (i2=0; i2 < 8; i2++) {
                  a1 = DxCol[i2];

                  {
                    int_T i3;
                    for (i3=0; i3 < 9; i3++) {
                      Ds[i2 * 9 + i3] += a1 * tmp2[i3];
                    }
                  }
                }
              }
            }
          }
        }
      }                                /* if (swChanged) */
    } while (swChanged > 0 && --loopsToDo > 0);

    if (loopsToDo == 0) {
      real_T *Ds = (real_T*)rtDW.StateSpace_PWORK.DS;

      {
        int_T i1;
        real_T *y0 = &rtDW.StateSpace_o1[0];
        for (i1=0; i1 < 8; i1++) {
          accum = 0.0;

          {
            int_T i2;
            const real_T *u0;
            for (i2=0; i2 < 6; i2++) {
              accum += *(Ds++) * 0.0;
            }

            accum += *(Ds++) * rtDW.ib[2];
            accum += *(Ds++) * rtDW.ib[3];
            accum += *(Ds++) * 500.0;
          }

          y0[i1] = accum;
        }
      }
    }

    /* Output new switches states */
    {
      int_T i1;
      real_T *y1 = &rtDW.StateSpace_o2[0];
      for (i1=0; i1 < 6; i1++) {
        y1[i1] = (real_T)switch_status[i1];
      }
    }
  }

  /* Gain: '<S39>/1_Vb' */
  rtb_u_Vb_idx_3 = 0.0055670221426890416 * rtDW.StateSpace_o1[7];

  /* Outputs for Enabled SubSystem: '<S39>/Rotor reference frame' incorporates:
   *  EnablePort: '<S50>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE_b) {
    /* Disable for Fcn: '<S50>/vqr' incorporates:
     *  Outport: '<S50>/vqr,vdr'
     */
    rtDW.vqr_o = 0.0;

    /* Disable for Fcn: '<S50>/vdr' incorporates:
     *  Outport: '<S50>/vqr,vdr'
     */
    rtDW.vdr_g = 0.0;

    /* Disable for Fcn: '<S50>/vqs' incorporates:
     *  Outport: '<S50>/vqs,vds'
     */
    rtDW.vqs_m = 0.0;

    /* Disable for Fcn: '<S50>/vds' incorporates:
     *  Outport: '<S50>/vqs,vds'
     */
    rtDW.vds_k = 0.0;
    rtDW.Rotorreferenceframe_MODE_b = false;
  }

  /* End of Outputs for SubSystem: '<S39>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S39>/Stationary reference frame' incorporates:
   *  EnablePort: '<S51>/Enable'
   */
  /* Fcn: '<S51>/vdr' incorporates:
   *  MultiPortSwitch: '<S41>/Multiport Switch'
   */
  vdr_c = ((-TrigonometricFunction_o1 - 1.7320508075688772 *
            TrigonometricFunction_o2) * 0.0 + -2.0 * TrigonometricFunction_o1 *
           0.0) * 0.33333333333333331;

  /* Fcn: '<S51>/vds' */
  vds_m = -0.57735026918962573 * rtb_u_Vb_idx_3;

  /* Fcn: '<S51>/vqr' incorporates:
   *  MultiPortSwitch: '<S41>/Multiport Switch'
   */
  TrigonometricFunction_o1 = ((TrigonometricFunction_o2 - 1.7320508075688772 *
    TrigonometricFunction_o1) * 0.0 + 2.0 * TrigonometricFunction_o2 * 0.0) *
    0.33333333333333331;

  /* Fcn: '<S51>/vqs' incorporates:
   *  Gain: '<S39>/1_Vb'
   */
  rtb_u_Vb_idx_3 = (0.0055670221426890416 * rtDW.StateSpace_o1[6] * 2.0 +
                    rtb_u_Vb_idx_3) * 0.33333333333333331;

  /* End of Outputs for SubSystem: '<S39>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S39>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S52>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MOD_g) {
    /* Disable for Fcn: '<S52>/vqr' incorporates:
     *  Outport: '<S52>/vqr,vdr'
     */
    rtDW.vqr = 0.0;

    /* Disable for Fcn: '<S52>/vdr' incorporates:
     *  Outport: '<S52>/vqr,vdr'
     */
    rtDW.vdr = 0.0;

    /* Disable for Fcn: '<S52>/vqs' incorporates:
     *  Outport: '<S52>/vqs,vds'
     */
    rtDW.vqs = 0.0;

    /* Disable for Fcn: '<S52>/vds' incorporates:
     *  Outport: '<S52>/vqs,vds'
     */
    rtDW.vds = 0.0;
    rtDW.Synchronousreferenceframe_MOD_g = false;
  }

  /* End of Outputs for SubSystem: '<S39>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S40>/Stationary reference frame' incorporates:
   *  EnablePort: '<S55>/Enable'
   */
  /* Gain: '<S35>/unit conversion' incorporates:
   *  Fcn: '<S55>/isa'
   *  MultiPortSwitch: '<S40>/Multiport Switch1'
   *  Sum: '<S40>/Sum3'
   *  UnitDelay: '<S45>/fluxes'
   */
  rtb_unitconversion_idx_11 = ((0.0 - rtb_Sum2[0]) - isb_g) * 13.731987951966302;

  /* End of Outputs for SubSystem: '<S40>/Stationary reference frame' */
  TrigonometricFunction_o2 = 0.57177765423802906 * rtDW.fluxes_DSTATE[2];

  /* Outputs for Enabled SubSystem: '<S40>/Stationary reference frame' incorporates:
   *  EnablePort: '<S55>/Enable'
   */
  rtb_unitconversion_idx_9 = 13.731987951966302 * rtb_Sum2[0];

  /* End of Outputs for SubSystem: '<S40>/Stationary reference frame' */
  rtb_unitconversion_idx_6 = 0.57177765423802906 * rtDW.fluxes_DSTATE[3];
  isb_g *= 13.731987951966302;

  /* Switch: '<S45>/IC' incorporates:
   *  DigitalClock: '<S45>/Digital Clock'
   *  Gain: '<S49>/wbase*Ts//2 '
   *  Product: '<S45>/Product1'
   *  Product: '<S45>/Product2'
   *  Sum: '<S45>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
   *  UnitDelay: '<S45>/fluxes'
   */
  if ((((rtM->Timing.clockTick0+rtM->Timing.clockTickH0* 4294967296.0)) * 1.0E-5)
      >= 1.0E-5) {
    for (i = 0; i <= 14; i += 2) {
      tmp_4 = _mm_loadu_pd(&W43wr[i]);
      tmp_4 = _mm_mul_pd(_mm_sub_pd(_mm_sub_pd(_mm_set1_pd(0.0), tmp_4),
        _mm_loadu_pd(&rtConstP.Constant4_Value[i])), _mm_set1_pd
                         (0.0015707963267948967));
      tmp_5 = _mm_loadu_pd(&rtConstP.u5_Value_p[i]);
      _mm_storeu_pd(&rtb_Sum5[i], _mm_add_pd(tmp_5, tmp_4));
      _mm_storeu_pd(&rtb_Lminrows24col24_0[i], _mm_sub_pd(tmp_5, tmp_4));
    }

    /* Product: '<S49>/inversion' incorporates:
     *  Assignment: '<S59>/W(4,3)=wr'
     *  Constant: '<S38>/Constant4'
     *  Constant: '<S49>/u5'
     *  Gain: '<S49>/wbase*Ts//2'
     *  Gain: '<S49>/wbase*Ts//2 '
     *  MultiPortSwitch: '<S41>/Multiport Switch1'
     *  Sum: '<S38>/Sum1'
     *  Sum: '<S49>/Sum1'
     *  Sum: '<S49>/Sum5'
     *  Switch: '<S38>/Switch1'
     */
    rt_invd4x4_snf(rtb_Lminrows24col24_0, W43wr);

    /* Product: '<S49>/Product4' incorporates:
     *  Gain: '<S49>/wbase*Ts//2 '
     *  Sum: '<S49>/Sum5'
     */
    for (i = 0; i < 4; i++) {
      rtb_Sum5_tmp = i << 2;
      rtb_Mult1_idx_1 = rtb_Sum5[rtb_Sum5_tmp + 1];
      rtb_Gain2_l = rtb_Sum5[rtb_Sum5_tmp];
      rtb_Mult1_idx_0 = rtb_Sum5[rtb_Sum5_tmp + 2];
      rtb_Saturation = rtb_Sum5[rtb_Sum5_tmp + 3];
      for (i_0 = 0; i_0 <= 2; i_0 += 2) {
        tmp_4 = _mm_loadu_pd(&W43wr[i_0 + 4]);
        tmp_5 = _mm_loadu_pd(&W43wr[i_0]);
        tmp_6 = _mm_loadu_pd(&W43wr[i_0 + 8]);
        tmp_7 = _mm_loadu_pd(&W43wr[i_0 + 12]);
        _mm_storeu_pd(&rtb_Lminrows24col24_0[i_0 + rtb_Sum5_tmp], _mm_add_pd
                      (_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd
          (rtb_Mult1_idx_1), tmp_4), _mm_mul_pd(_mm_set1_pd(rtb_Gain2_l), tmp_5)),
          _mm_mul_pd(_mm_set1_pd(rtb_Mult1_idx_0), tmp_6)), _mm_mul_pd
                       (_mm_set1_pd(rtb_Saturation), tmp_7)));
      }
    }

    /* End of Product: '<S49>/Product4' */

    /* Sum: '<S45>/sum' incorporates:
     *  MultiPortSwitch: '<S39>/Multiport Switch'
     *  MultiPortSwitch: '<S39>/Multiport Switch1'
     *  UnitDelay: '<S45>/voltages'
     */
    rtb_Switch_m_idx_0 = rtb_u_Vb_idx_3 + rtDW.voltages_DSTATE[0];
    rtb_MultiportSwitch1_idx_2 = TrigonometricFunction_o1 +
      rtDW.voltages_DSTATE[2];
    rtb_MultiportSwitch1_f = vds_m + rtDW.voltages_DSTATE[1];
    rtb_MultiportSwitch1_idx_3 = vdr_c + rtDW.voltages_DSTATE[3];

    /* UnitDelay: '<S45>/fluxes' incorporates:
     *  Product: '<S45>/Product2'
     */
    rtb_Mult1_idx_1 = rtDW.fluxes_DSTATE[1];
    rtb_Gain2_l = rtDW.fluxes_DSTATE[0];
    rtb_Mult1_idx_0 = rtDW.fluxes_DSTATE[2];
    rtb_Saturation = rtDW.fluxes_DSTATE[3];
    for (i = 0; i <= 2; i += 2) {
      tmp_4 = _mm_loadu_pd(&W43wr[i + 4]);
      tmp_5 = _mm_set1_pd(0.0015707963267948967);
      tmp_6 = _mm_loadu_pd(&W43wr[i]);
      tmp_7 = _mm_loadu_pd(&W43wr[i + 8]);
      tmp = _mm_loadu_pd(&W43wr[i + 12]);
      tmp_0 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 4]);
      tmp_1 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i]);
      tmp_2 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 8]);
      tmp_3 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 12]);
      _mm_storeu_pd(&rtb_MultiportSwitch[i], _mm_add_pd(_mm_add_pd(_mm_add_pd
        (_mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_4, tmp_5), _mm_set1_pd
        (rtb_MultiportSwitch1_f)), _mm_mul_pd(_mm_mul_pd(tmp_5, tmp_6),
        _mm_set1_pd(rtb_Switch_m_idx_0))), _mm_mul_pd(_mm_mul_pd(tmp_7, tmp_5),
        _mm_set1_pd(rtb_MultiportSwitch1_idx_2))), _mm_mul_pd(_mm_mul_pd(tmp,
        tmp_5), _mm_set1_pd(rtb_MultiportSwitch1_idx_3))), _mm_add_pd(_mm_add_pd
        (_mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd(rtb_Mult1_idx_1)), _mm_mul_pd
                    (tmp_1, _mm_set1_pd(rtb_Gain2_l))), _mm_mul_pd(tmp_2,
        _mm_set1_pd(rtb_Mult1_idx_0))), _mm_mul_pd(tmp_3, _mm_set1_pd
        (rtb_Saturation)))));
    }
  } else {
    rtb_MultiportSwitch[0] = rtDW.fluxes_DSTATE[0];
    rtb_MultiportSwitch[1] = rtDW.fluxes_DSTATE[1];
    rtb_MultiportSwitch[2] = rtDW.fluxes_DSTATE[2];
    rtb_MultiportSwitch[3] = rtDW.fluxes_DSTATE[3];
  }

  /* End of Switch: '<S45>/IC' */
  for (i = 0; i <= 0; i += 2) {
    /* Gain: '<S22>/Gain1' incorporates:
     *  Gain: '<S22>/Gain3'
     */
    _mm_storeu_pd(&rtb_MathFunction[i], _mm_mul_pd(_mm_add_pd(_mm_add_pd
      (_mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i + 3]), _mm_set1_pd(isb_g)),
       _mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i]), _mm_set1_pd
                  (rtb_unitconversion_idx_9))), _mm_mul_pd(_mm_loadu_pd
      (&rtConstP.Gain3_Gain[i + 6]), _mm_set1_pd(rtb_unitconversion_idx_11))),
      _mm_set1_pd(0.66666666666666663)));
  }

  for (i = 2; i < 3; i++) {
    /* Gain: '<S22>/Gain1' incorporates:
     *  Gain: '<S22>/Gain3'
     */
    rtb_MathFunction[i] = ((rtConstP.Gain3_Gain[i + 3] * isb_g +
      rtConstP.Gain3_Gain[i] * rtb_unitconversion_idx_9) + rtConstP.Gain3_Gain[i
      + 6] * rtb_unitconversion_idx_11) * 0.66666666666666663;
  }

  /* Gain: '<S37>/1\p1' */
  rtb_unitconversion_idx_11 = 157.07963267948966 * rtb_phimd;

  /* ComplexToMagnitudeAngle: '<S1>/Complex to Magnitude-Angle' incorporates:
   *  RealImagToComplex: '<S1>/Real-Imag to Complex'
   */
  rtb_unitconversion_idx_9 = rt_atan2d_snf(TrigonometricFunction_o2,
    rtb_unitconversion_idx_6);

  /* Outputs for Enabled SubSystem: '<S27>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S30>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S21>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  /* Fcn: '<S25>/Fcn' incorporates:
   *  Fcn: '<S25>/Fcn1'
   *  Fcn: '<S30>/Fcn'
   */
  isb_g = cos(rtb_unitconversion_idx_9);
  rtb_Gain2_l = sin(rtb_unitconversion_idx_9);

  /* End of Outputs for SubSystem: '<S27>/Subsystem - pi//2 delay' */

  /* Switch: '<S21>/Switch' incorporates:
   *  Fcn: '<S25>/Fcn'
   *  Fcn: '<S25>/Fcn1'
   */
  rtb_Mult1_idx_0 = rtb_MathFunction[0] * rtb_Gain2_l - rtb_MathFunction[1] *
    isb_g;
  rtb_Mult1_idx_1 = rtb_MathFunction[0] * isb_g + rtb_MathFunction[1] *
    rtb_Gain2_l;

  /* End of Outputs for SubSystem: '<S21>/Subsystem - pi//2 delay' */

  /* Outputs for Enabled SubSystem: '<S27>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S30>/Enable'
   */
  /* Switch: '<S27>/Switch' incorporates:
   *  Fcn: '<S30>/Fcn'
   */
  rtb_Switch_m_idx_0 = rtb_unitconversion_idx_6 * rtb_Gain2_l -
    TrigonometricFunction_o2 * isb_g;

  /* End of Outputs for SubSystem: '<S27>/Subsystem - pi//2 delay' */
  if (rtM->Timing.TaskCounters.TID[1] == 0) {
    /* DigitalClock: '<S10>/Clock' */
    rtDW.Clock = ((rtM->Timing.clockTick1) * 1.0);
  }

  /* Switch: '<S10>/Switch' incorporates:
   *  Product: '<S10>/Divide'
   */
  if (rtDW.Clock > 0.0001) {
    TrigonometricFunction_o2 = rtb_Mult1_idx_1 / rtb_Mult1_idx_0;
  } else {
    TrigonometricFunction_o2 = 0.0;
  }

  /* End of Switch: '<S10>/Switch' */

  /* Gain: '<S10>/Gain' */
  rtb_unitconversion_idx_6 = 6.23416784551107 * TrigonometricFunction_o2;

  /* Gain: '<S10>/Gain1' */
  TrigonometricFunction_o2 = 2.0 * rtb_unitconversion_idx_11;

  /* Sum: '<S10>/Add' */
  rtb_unitconversion_idx_6 += TrigonometricFunction_o2;

  /* Sum: '<S4>/Sum1' incorporates:
   *  Constant: '<S1>/Constant'
   *  DiscreteIntegrator: '<S4>/Discrete-Time Integrator'
   *  Gain: '<S10>/Gain2'
   *  Gain: '<S10>/Gain3'
   *  Gain: '<S10>/Gain4'
   *  Gain: '<S10>/Gain5'
   *  Gain: '<S4>/Gain'
   *  Product: '<S10>/Product'
   *  Sum: '<S10>/Add2'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S4>/Sum'
   */
  isb_g = ((0.0 - rtb_unitconversion_idx_6 * rtb_Mult1_idx_1 *
            0.058148172073290927 * 0.060793999999999994) - 0.022361 *
           rtb_Switch_m_idx_0 * 270.56932718377169) + ((9.4040441629117311 -
    rtb_Mult1_idx_0) * 15.7034858759281 + rtDW.DiscreteTimeIntegrator_DSTATE);

  /* Saturate: '<S4>/Saturation' */
  if (isb_g > 359.2584956081995) {
    rtb_Gain2_l = 359.2584956081995;
  } else if (isb_g < -359.2584956081995) {
    rtb_Gain2_l = -359.2584956081995;
  } else {
    rtb_Gain2_l = isb_g;
  }

  /* End of Saturate: '<S4>/Saturation' */

  /* Sum: '<S4>/Sum3' incorporates:
   *  Constant: '<S1>/Constant'
   *  Gain: '<S4>/Gain2'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S4>/Sum2'
   */
  isb_g = (9.4040441629117311 - rtb_Mult1_idx_0) - (isb_g - rtb_Gain2_l) *
    0.063680128596982508;

  /* Product: '<S10>/Product1' */
  rtb_Mult1_idx_0 *= rtb_unitconversion_idx_6;

  /* Sum: '<S6>/Sum' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S6>/Gain'
   *  Sum: '<S1>/Sum'
   */
  rtb_unitconversion_idx_6 = (52.359877559829883 - rtb_unitconversion_idx_11) *
    13.859596918105916 + rtDW.DiscreteTimeIntegrator_DSTATE_b;

  /* Saturate: '<S6>/Saturation' */
  if (rtb_unitconversion_idx_6 > 48.400544337535294) {
    rtb_Saturation = 48.400544337535294;
  } else if (rtb_unitconversion_idx_6 < -48.400544337535294) {
    rtb_Saturation = -48.400544337535294;
  } else {
    rtb_Saturation = rtb_unitconversion_idx_6;
  }

  /* End of Saturate: '<S6>/Saturation' */

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S6>/Gain2'
   */
  rtb_Mult1_idx_1 = 0.61904294884986732 * rtb_Saturation - rtb_Mult1_idx_1;

  /* Sum: '<S5>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S5>/Discrete-Time Integrator'
   *  Gain: '<S10>/Gain6'
   *  Gain: '<S10>/Gain7'
   *  Gain: '<S10>/Gain8'
   *  Gain: '<S5>/Gain'
   *  Product: '<S10>/Product2'
   *  Sum: '<S10>/Add1'
   *  Sum: '<S5>/Sum'
   */
  TrigonometricFunction_o2 = (0.058148172073290927 * rtb_Mult1_idx_0 *
    0.060793999999999994 + rtb_Switch_m_idx_0 * TrigonometricFunction_o2 *
    0.97049050893180255) + (15.7034858759281 * rtb_Mult1_idx_1 +
    rtDW.DiscreteTimeIntegrator_DSTATE_l);

  /* Saturate: '<S5>/Saturation' */
  if (TrigonometricFunction_o2 > 359.2584956081995) {
    rtb_Mult1_idx_0 = 359.2584956081995;
  } else if (TrigonometricFunction_o2 < -359.2584956081995) {
    rtb_Mult1_idx_0 = -359.2584956081995;
  } else {
    rtb_Mult1_idx_0 = TrigonometricFunction_o2;
  }

  /* End of Saturate: '<S5>/Saturation' */

  /* Fcn: '<S7>/alpha' incorporates:
   *  Fcn: '<S7>/beta'
   */
  rtb_Switch_m_idx_0 = sin(1.5707963267948966 - rtb_unitconversion_idx_9);
  rtb_MultiportSwitch1_idx_2 = cos(1.5707963267948966 - rtb_unitconversion_idx_9);
  rtb_MultiportSwitch1_f = rtb_MultiportSwitch1_idx_2 * rtb_Gain2_l +
    rtb_Switch_m_idx_0 * rtb_Mult1_idx_0;

  /* Gain: '<S11>/Gain' */
  rtb_unitconversion_idx_9 = 1.7320508075688772 * rtb_MultiportSwitch1_f;

  /* Fcn: '<S7>/beta' */
  rtb_Gain2_l = -rtb_Switch_m_idx_0 * rtb_Gain2_l + rtb_MultiportSwitch1_idx_2 *
    rtb_Mult1_idx_0;

  /* Sum: '<S11>/Subtract2' incorporates:
   *  Constant: '<S17>/Constant'
   *  Constant: '<S18>/Constant'
   *  Constant: '<S19>/Constant'
   *  Gain: '<S11>/Gain2'
   *  Gain: '<S11>/Gain3'
   *  RelationalOperator: '<S17>/Compare'
   *  RelationalOperator: '<S18>/Compare'
   *  RelationalOperator: '<S19>/Compare'
   *  Sum: '<S11>/Subtract'
   *  Sum: '<S11>/Subtract1'
   */
  rtb_Subtract2_m = (uint8_T)(((uint32_T)((rtb_unitconversion_idx_9 -
    rtb_Gain2_l > 0.0) << 1) + (uint32_T)(rtb_Gain2_l > 0.0)) + (uint32_T)(((0.0
    - rtb_unitconversion_idx_9) - rtb_Gain2_l > 0.0) << 2));

  /* Gain: '<S16>/Gain' */
  rtb_Switch_m_idx_0 = 3.4641016151377545E-7 * rtb_Gain2_l;

  /* Gain: '<S16>/Gain2' */
  rtb_Gain2_l *= 0.5;

  /* Gain: '<S16>/Gain1' */
  rtb_MultiportSwitch1_f *= 0.8660254037844386;

  /* MultiPortSwitch: '<S12>/Multiport Switch' incorporates:
   *  Gain: '<S12>/Gain'
   *  Gain: '<S12>/Gain1'
   *  Gain: '<S12>/Gain2'
   *  Gain: '<S16>/Gain3'
   *  Gain: '<S16>/Gain4'
   *  Sum: '<S16>/Subtract'
   *  Sum: '<S16>/Subtract1'
   */
  switch (rtb_Subtract2_m) {
   case 1:
    rtb_unitconversion_idx_9 = (rtb_Gain2_l - rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7;

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S16>/Gain3'
     *  Gain: '<S16>/Gain4'
     *  Sum: '<S16>/Subtract'
     *  Sum: '<S16>/Subtract1'
     */
    rtb_MultiportSwitch1_f = (rtb_Gain2_l + rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7;
    break;

   case 2:
    rtb_unitconversion_idx_9 = (rtb_Gain2_l + rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7;

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S12>/Gain'
     *  Gain: '<S16>/Gain3'
     *  Sum: '<S16>/Subtract'
     */
    rtb_MultiportSwitch1_f = -rtb_Switch_m_idx_0;
    break;

   case 3:
    rtb_unitconversion_idx_9 = -((rtb_Gain2_l - rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7);

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S12>/Gain2'
     *  Gain: '<S16>/Gain4'
     *  Sum: '<S16>/Subtract1'
     */
    rtb_MultiportSwitch1_f = rtb_Switch_m_idx_0;
    break;

   case 4:
    rtb_unitconversion_idx_9 = -rtb_Switch_m_idx_0;

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S12>/Gain'
     *  Gain: '<S16>/Gain4'
     *  Sum: '<S16>/Subtract1'
     */
    rtb_MultiportSwitch1_f = (rtb_Gain2_l - rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7;
    break;

   case 5:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S12>/Gain1'
     *  Gain: '<S16>/Gain3'
     *  Sum: '<S16>/Subtract'
     */
    rtb_MultiportSwitch1_f = -((rtb_Gain2_l + rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7);
    break;

   default:
    rtb_unitconversion_idx_9 = -((rtb_Gain2_l + rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7);

    /* MultiPortSwitch: '<S12>/Multiport Switch1' incorporates:
     *  Gain: '<S12>/Gain1'
     *  Gain: '<S12>/Gain2'
     *  Gain: '<S16>/Gain3'
     *  Gain: '<S16>/Gain4'
     *  Sum: '<S16>/Subtract'
     *  Sum: '<S16>/Subtract1'
     */
    rtb_MultiportSwitch1_f = -((rtb_Gain2_l - rtb_MultiportSwitch1_f) *
      3.4641016151377545E-7);
    break;
  }

  /* End of MultiPortSwitch: '<S12>/Multiport Switch' */

  /* Gain: '<S14>/Gain' incorporates:
   *  Constant: '<S14>/Constant'
   *  Sum: '<S14>/Subtract'
   */
  rtb_Gain2_l = ((0.0001 - rtb_unitconversion_idx_9) - rtb_MultiportSwitch1_f) *
    0.25;

  /* Sum: '<S14>/Subtract1' incorporates:
   *  Gain: '<S14>/Gain1'
   */
  rtb_Switch_m_idx_0 = 0.5 * rtb_unitconversion_idx_9 + rtb_Gain2_l;

  /* Sum: '<S14>/Subtract2' incorporates:
   *  Gain: '<S14>/Gain2'
   */
  rtb_MultiportSwitch1_f = 0.5 * rtb_MultiportSwitch1_f + rtb_Switch_m_idx_0;

  /* MultiPortSwitch: '<S13>/Multiport Switch1' */
  switch (rtb_Subtract2_m) {
   case 1:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;

   case 2:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;

   case 3:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;

   case 4:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;

   case 5:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;

   default:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;
  }

  /* End of MultiPortSwitch: '<S13>/Multiport Switch1' */
  if (rtM->Timing.TaskCounters.TID[1] == 0) {
    /* Lookup: '<S20>/Look-Up Table' incorporates:
     *  DigitalClock: '<S20>/Digital Clock'
     *  Fcn: '<S20>/Fcn1'
     */
    rtDW.LookUpTable = rt_Lookup(&rtConstP.LookUpTable_XData[0], 3, rt_remd_snf
      (((rtM->Timing.clockTick1) * 1.0), 0.0001), &rtConstP.LookUpTable_YData[0]);
  }

  /* RelationalOperator: '<S15>/Relational Operator1' */
  rtb_LogicalOperator2 = (rtDW.LookUpTable > rtb_unitconversion_idx_9);

  /* DataTypeConversion: '<S15>/Cast To Double2' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[2] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S15>/Cast To Double3' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   *  Logic: '<S15>/Logical Operator1'
   */
  rtDW.VectorConcatenate[3] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S13>/Multiport Switch3' */
  switch (rtb_Subtract2_m) {
   case 1:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;

   case 2:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;

   case 3:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;

   case 4:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;

   case 5:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;

   default:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;
  }

  /* End of MultiPortSwitch: '<S13>/Multiport Switch3' */

  /* RelationalOperator: '<S15>/Relational Operator' */
  rtb_LogicalOperator2 = (rtDW.LookUpTable > rtb_unitconversion_idx_9);

  /* DataTypeConversion: '<S15>/Cast To Double' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[0] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S15>/Cast To Double1' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   *  Logic: '<S15>/Logical Operator'
   */
  rtDW.VectorConcatenate[1] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S13>/Multiport Switch2' */
  switch (rtb_Subtract2_m) {
   case 1:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;

   case 2:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;

   case 3:
    rtb_unitconversion_idx_9 = rtb_MultiportSwitch1_f;
    break;

   case 4:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;

   case 5:
    rtb_unitconversion_idx_9 = rtb_Switch_m_idx_0;
    break;

   default:
    rtb_unitconversion_idx_9 = rtb_Gain2_l;
    break;
  }

  /* End of MultiPortSwitch: '<S13>/Multiport Switch2' */

  /* RelationalOperator: '<S15>/Relational Operator2' */
  rtb_LogicalOperator2 = (rtDW.LookUpTable > rtb_unitconversion_idx_9);

  /* DataTypeConversion: '<S15>/Cast To Double4' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[4] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S15>/Cast To Double5' incorporates:
   *  Concatenate: '<S15>/Vector Concatenate'
   *  Logic: '<S15>/Logical Operator2'
   */
  rtDW.VectorConcatenate[5] = !rtb_LogicalOperator2;

  /* Gain: '<S37>/1_2H' incorporates:
   *  Gain: '<S37>/F'
   *  Gain: '<S42>/1-1'
   *  Product: '<S42>/Mult1'
   *  Sum: '<S37>/Sum'
   *  Sum: '<S42>/Sum2'
   *  UnitDelay: '<S45>/fluxes'
   */
  rtb_Gain2_l = (((rtb_Sum2[0] * rtDW.fluxes_DSTATE[1] + rtDW.fluxes_DSTATE[0] *
                   -rtb_Sum2[1]) - 1.6981581911296177) - 0.0 * rtb_phimd) *
    5.9506091980420592;

  /* DiscreteIntegrator: '<S37>/Rotor speed(wm)' */
  if (rtDW.Rotorspeedwm_SYSTEM_ENABLE != 0) {
    /* DiscreteIntegrator: '<S37>/Rotor speed(wm)' */
    rtb_unitconversion_idx_9 = rtDW.Rotorspeedwm_DSTATE;
  } else {
    /* DiscreteIntegrator: '<S37>/Rotor speed(wm)' */
    rtb_unitconversion_idx_9 = 5.0E-6 * rtb_Gain2_l + rtDW.Rotorspeedwm_DSTATE;
  }

  /* End of DiscreteIntegrator: '<S37>/Rotor speed(wm)' */

  /* Update for UnitDelay: '<S45>/fluxes' */
  rtDW.fluxes_DSTATE[0] = rtb_MultiportSwitch[0];
  rtDW.fluxes_DSTATE[1] = rtb_MultiportSwitch[1];
  rtDW.fluxes_DSTATE[2] = rtb_MultiportSwitch[2];
  rtDW.fluxes_DSTATE[3] = rtb_MultiportSwitch[3];

  /* Update for DiscreteIntegrator: '<S37>/Rotor angle thetam' incorporates:
   *  Gain: '<S37>/web_psb'
   */
  rtDW.Rotoranglethetam_DSTATE += 314.15926535897933 * rtb_phimd * 1.0E-5;

  /* Update for UnitDelay: '<S60>/wm_delay' */
  rtDW.wm_delay_DSTATE = rtb_unitconversion_idx_9;

  /* Update for UnitDelay: '<S60>/wm_predict' */
  rtDW.wm_predict_DSTATE = rtb_Phisat;

  /* Update for S-Function (sfun_spssw_discc): '<S64>/State-Space' incorporates:
   *  Constant: '<S61>/DC'
   */
  {
    int_T *gState = (int_T*)rtDW.StateSpace_PWORK.G_STATE;

    /* Store switch gates values for next step */
    {
      int_T i1;
      const real_T *u1 = &rtDW.VectorConcatenate[0];
      for (i1=0; i1 < 6; i1++) {
        *(gState++) = (int_T) u1[i1];
      }
    }
  }

  /* Update for UnitDelay: '<S45>/voltages' incorporates:
   *  MultiPortSwitch: '<S39>/Multiport Switch'
   *  MultiPortSwitch: '<S39>/Multiport Switch1'
   */
  rtDW.voltages_DSTATE[0] = rtb_u_Vb_idx_3;
  rtDW.voltages_DSTATE[2] = TrigonometricFunction_o1;
  rtDW.voltages_DSTATE[1] = vds_m;
  rtDW.voltages_DSTATE[3] = vdr_c;

  /* Update for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S4>/Gain1'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE += 2896.152729496454 * isb_g * 1.0E-5;

  /* Update for DiscreteIntegrator: '<S5>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S5>/Gain1'
   *  Gain: '<S5>/Gain2'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE_l += (rtb_Mult1_idx_1 -
    (TrigonometricFunction_o2 - rtb_Mult1_idx_0) * 0.063680128596982508) *
    2896.152729496454 * 1.0E-5;

  /* Update for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S6>/Gain1'
   *  Gain: '<S6>/Gain3'
   *  Sum: '<S1>/Sum'
   *  Sum: '<S6>/Sum1'
   *  Sum: '<S6>/Sum2'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE_b += ((52.359877559829883 -
    rtb_unitconversion_idx_11) - (rtb_unitconversion_idx_6 - rtb_Saturation) *
    0.072152170507471164) * 2462.6907165410848 * 1.0E-5;

  /* Update for DiscreteIntegrator: '<S37>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 0U;
  rtDW.Rotorspeedwm_DSTATE = 5.0E-6 * rtb_Gain2_l + rtb_unitconversion_idx_9;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The resolution of this integer timer is 1.0E-5, which is the step size
   * of the task. Size of "clockTick0" ensures timer will not overflow during the
   * application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  rtM->Timing.clockTick0++;
  if (!rtM->Timing.clockTick0) {
    rtM->Timing.clockTickH0++;
  }

  if (rtM->Timing.TaskCounters.TID[1] == 0) {
    /* Update absolute timer for sample time: [1.0s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 1.0, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     */
    rtM->Timing.clockTick1++;
  }

  rate_scheduler();
}

/* Model initialize function */
void sumAll_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* Start for S-Function (sfun_spssw_discc): '<S64>/State-Space' incorporates:
   *  Constant: '<S61>/DC'
   */

  /* S-Function block: <S64>/State-Space */
  {
    rtDW.StateSpace_PWORK.DS = (real_T*)calloc(8 * 9, sizeof(real_T));
    rtDW.StateSpace_PWORK.DX_COL = (real_T*)calloc(8, sizeof(real_T));
    rtDW.StateSpace_PWORK.TMP2 = (real_T*)calloc(9, sizeof(real_T));
    rtDW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.G_STATE = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(6, sizeof(real_T));
    rtDW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)calloc(6, sizeof(int_T));
    rtDW.StateSpace_PWORK.USWLAST = (real_T*)calloc(6, sizeof(real_T));
  }

  /* InitializeConditions for S-Function (sfun_spssw_discc): '<S64>/State-Space' incorporates:
   *  Constant: '<S61>/DC'
   */
  {
    int32_T i, j;
    real_T *Ds = (real_T*)rtDW.StateSpace_PWORK.DS;

    /* Copy and transpose D */
    for (i=0; i<8; i++) {
      for (j=0; j<9; j++)
        Ds[i*9 + j] = (rtConstP.StateSpace_DS_param[i + j*8]);
    }

    {
      /* Switches work vectors */
      int_T *switch_status = (int_T*) rtDW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *gState = (int_T*)rtDW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)rtDW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*)rtDW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*)rtDW.StateSpace_PWORK.IDX_OUT_SW;
      int_T *switch_status_init = (int_T*)
        rtDW.StateSpace_PWORK.SWITCH_STATUS_INIT;

      /* Initialize work vectors */
      switch_status[0] = 0;
      switch_status_init[0] = 0;
      gState[0] = (int_T) 0.0;
      yswitch[0] = 1/0.001;
      switchTypes[0] = (int_T)7.0;
      idxOutSw[0] = ((int_T)0.0) - 1;
      switch_status[1] = 0;
      switch_status_init[1] = 0;
      gState[1] = (int_T) 0.0;
      yswitch[1] = 1/0.001;
      switchTypes[1] = (int_T)7.0;
      idxOutSw[1] = ((int_T)0.0) - 1;
      switch_status[2] = 0;
      switch_status_init[2] = 0;
      gState[2] = (int_T) 0.0;
      yswitch[2] = 1/0.001;
      switchTypes[2] = (int_T)7.0;
      idxOutSw[2] = ((int_T)0.0) - 1;
      switch_status[3] = 0;
      switch_status_init[3] = 0;
      gState[3] = (int_T) 0.0;
      yswitch[3] = 1/0.001;
      switchTypes[3] = (int_T)7.0;
      idxOutSw[3] = ((int_T)0.0) - 1;
      switch_status[4] = 0;
      switch_status_init[4] = 0;
      gState[4] = (int_T) 0.0;
      yswitch[4] = 1/0.001;
      switchTypes[4] = (int_T)7.0;
      idxOutSw[4] = ((int_T)0.0) - 1;
      switch_status[5] = 0;
      switch_status_init[5] = 0;
      gState[5] = (int_T) 0.0;
      yswitch[5] = 1/0.001;
      switchTypes[5] = (int_T)7.0;
      idxOutSw[5] = ((int_T)0.0) - 1;
    }
  }

  /* Enable for DiscreteIntegrator: '<S37>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 1U;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
