/*
 * File: motor_model2.c
 *
 * Code generated for Simulink model 'motor_model2'.
 *
 * Model version                  : 13.20
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Mon Nov 11 19:15:56 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "motor_model2.h"
#include "rtwtypes.h"
#include <string.h>
#include <emmintrin.h>
#include <float.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#define NumBitsPerChar                 8U

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

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

const real_T motor_model2_RGND = 0.0;  /* real_T ground */

/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Real-time model */
static RT_MODEL rtM_;
RT_MODEL *const rtM = &rtM_;
extern void rt_invd4x4_snf(const real_T u[16], real_T y[16]);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
static real_T look1_binlx(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T maxIndex);
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

static real_T look1_binlx(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T iLeft;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    uint32_T bpIdx;
    uint32_T iRght;

    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
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
void motor_model2_step(void)
{
  /* local block i/o variables */
  real_T rtb_xk1[4];
  real_T rtb_wm_delay;
  real_T rtb_MultiportSwitch_p[2];
  real_T rtb_MultiportSwitch1[2];
  real_T rtb_IC[4];
  real_T rtb_u_2H;
  real_T rtb_web_psb;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  __m128d tmp_6;
  __m128d tmp_7;
  __m128d tmp_8;
  real_T W43wr[16];
  real_T rtb_Lminrows24col24_0[16];
  real_T rtb_Sum5[16];
  real_T rtb_Sum2[4];
  real_T rtb_MathFunction[3];
  real_T TrigonometricFunction_o1;
  real_T TrigonometricFunction_o2;
  real_T isb_d;
  real_T rtb_DiscreteTimeIntegrator;
  real_T rtb_Gain2_pc;
  real_T rtb_MultiportSwitch2;
  real_T rtb_MultiportSwitch_h_idx_0;
  real_T rtb_Phisat;
  real_T rtb_Product2_e;
  real_T rtb_Subtract2;
  real_T rtb_Switch_idx_0;
  real_T rtb_u_Vb_idx_3;
  real_T rtb_unitconversion_idx_9;
  real_T tmp;
  int32_T i;
  int32_T i_0;
  int32_T rtb_Sum5_tmp;
  uint8_T rtb_Subtract2_o;
  boolean_T rtb_LogicalOperator2;

  /* UnitDelay: '<S23>/fluxes' */
  rtb_xk1[0] = rtDW.fluxes_DSTATE[0];
  rtb_xk1[1] = rtDW.fluxes_DSTATE[1];
  rtb_xk1[2] = rtDW.fluxes_DSTATE[2];
  rtb_xk1[3] = rtDW.fluxes_DSTATE[3];

  /* Switch: '<S16>/Switch' incorporates:
   *  Constant: '<S16>/Constant2'
   *  Product: '<S22>/inversion'
   */
  memcpy(&rtb_Lminrows24col24_0[0], &rtConstP.Constant2_Value[0], sizeof(real_T)
         << 4U);

  /* Product: '<S16>/Product3' */
  rtb_Product2_e = rtb_xk1[1];
  rtb_DiscreteTimeIntegrator = rtb_xk1[0];
  rtb_Gain2_pc = rtb_xk1[2];
  rtb_MultiportSwitch_h_idx_0 = rtb_xk1[3];
  for (i = 0; i <= 2; i += 2) {
    /* Product: '<S16>/Product3' */
    tmp_5 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 4]);
    tmp_6 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i]);
    tmp_7 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 8]);
    tmp_8 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 12]);
    _mm_storeu_pd(&rtb_Sum2[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_5, _mm_set1_pd(rtb_Product2_e)), _mm_mul_pd(tmp_6, _mm_set1_pd
      (rtb_DiscreteTimeIntegrator))), _mm_mul_pd(tmp_7, _mm_set1_pd(rtb_Gain2_pc))),
      _mm_mul_pd(tmp_8, _mm_set1_pd(rtb_MultiportSwitch_h_idx_0))));
  }

  /* UnitDelay: '<S38>/wm_delay' */
  rtb_wm_delay = rtDW.wm_delay_DSTATE;

  /* Sum: '<S38>/Sum1' incorporates:
   *  Gain: '<S38>/F2'
   *  UnitDelay: '<S38>/wm_predict'
   */
  rtb_Phisat = 2.0 * rtb_wm_delay - rtDW.wm_predict_DSTATE;

  /* Outputs for Enabled SubSystem: '<S19>/sin(thr),cos(thr)' incorporates:
   *  EnablePort: '<S36>/Enable'
   */
  if (rtDW.sinthrcosthr_MODE) {
    /* Disable for Trigonometry: '<S36>/Trigonometric Function' incorporates:
     *  Outport: '<S36>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o1_e = 0.0;

    /* Disable for Trigonometry: '<S36>/Trigonometric Function' incorporates:
     *  Outport: '<S36>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o2_l = 0.0;

    /* Disable for Outport: '<S36>/sin(thr),cos(thr)' incorporates:
     *  Constant: '<S36>/Constant'
     */
    rtDW.Constant_f[0] = 0.0;
    rtDW.Constant_f[1] = 0.0;

    /* Disable for Assignment: '<S36>/W(2,1)=-wr' incorporates:
     *  Outport: '<S36>/W'
     */
    memset(&rtDW.W21wr[0], 0, sizeof(real_T) << 4U);
    rtDW.sinthrcosthr_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S19>/sin(thr),cos(thr)' */

  /* Outputs for Enabled SubSystem: '<S19>/sin(thr),cos(thr)1' incorporates:
   *  EnablePort: '<S37>/Enable'
   */
  /* Assignment: '<S37>/W(4,3)=wr' incorporates:
   *  SignalConversion generated from: '<S37>/W(3,4)=-wr'
   */
  memset(&W43wr[0], 0, sizeof(real_T) << 4U);

  /* Gain: '<S37>/Gain3' incorporates:
   *  Assignment: '<S37>/W(4,3)=wr'
   */
  W43wr[14] = -rtb_Phisat;

  /* Assignment: '<S37>/W(4,3)=wr' */
  W43wr[11] = rtb_Phisat;

  /* Trigonometry: '<S37>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S15>/Rotor angle thetam'
   */
  TrigonometricFunction_o1 = sin(rtDW.Rotoranglethetam_DSTATE);

  /* Trigonometry: '<S37>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S15>/Rotor angle thetam'
   */
  TrigonometricFunction_o2 = cos(rtDW.Rotoranglethetam_DSTATE);

  /* End of Outputs for SubSystem: '<S19>/sin(thr),cos(thr)1' */

  /* Outputs for Enabled SubSystem: '<S18>/Rotor reference frame' incorporates:
   *  EnablePort: '<S32>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE) {
    /* Disable for Fcn: '<S32>/ira' incorporates:
     *  Outport: '<S32>/ira,irb'
     */
    rtDW.ira_e = 0.0;

    /* Disable for Fcn: '<S32>/irb' incorporates:
     *  Outport: '<S32>/ira,irb'
     */
    rtDW.irb_f = 0.0;

    /* Disable for Fcn: '<S32>/isa' incorporates:
     *  Outport: '<S32>/isa,isb'
     */
    rtDW.isa_i = 0.0;

    /* Disable for Fcn: '<S32>/isb' incorporates:
     *  Outport: '<S32>/isa,isb'
     */
    rtDW.isb_e = 0.0;
    rtDW.Rotorreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S18>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S18>/Stationary reference frame' incorporates:
   *  EnablePort: '<S33>/Enable'
   */
  /* Fcn: '<S33>/isb' */
  isb_d = -(1.7320508075688772 * rtb_Sum2[1] + rtb_Sum2[0]) / 2.0;

  /* End of Outputs for SubSystem: '<S18>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S18>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S34>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MODE) {
    /* Disable for Fcn: '<S34>/ira' incorporates:
     *  Outport: '<S34>/ira,irb'
     */
    rtDW.ira = 0.0;

    /* Disable for Fcn: '<S34>/irb' incorporates:
     *  Outport: '<S34>/ira,irb'
     */
    rtDW.irb = 0.0;

    /* Disable for Fcn: '<S34>/isa' incorporates:
     *  Outport: '<S34>/isa,isb'
     */
    rtDW.isa = 0.0;

    /* Disable for Fcn: '<S34>/isb' incorporates:
     *  Outport: '<S34>/isa,isb'
     */
    rtDW.isb = 0.0;
    rtDW.Synchronousreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S18>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S18>/Stationary reference frame' incorporates:
   *  EnablePort: '<S33>/Enable'
   */
  /* Gain: '<S18>/ib' incorporates:
   *  Fcn: '<S33>/ira'
   *  Fcn: '<S33>/irb'
   *  Fcn: '<S33>/isa'
   *  MultiPortSwitch: '<S18>/Multiport Switch1'
   *  MultiPortSwitch: '<S19>/Multiport Switch'
   */
  rtDW.ib[0] = (TrigonometricFunction_o2 * rtb_Sum2[2] -
                TrigonometricFunction_o1 * rtb_Sum2[3]) * 13.731987951966302;
  rtDW.ib[2] = 13.731987951966302 * rtb_Sum2[0];
  rtDW.ib[1] = ((-TrigonometricFunction_o2 - 1.7320508075688772 *
                 TrigonometricFunction_o1) * rtb_Sum2[2] +
                (TrigonometricFunction_o1 - 1.7320508075688772 *
                 TrigonometricFunction_o2) * rtb_Sum2[3]) / 2.0 *
    13.731987951966302;

  /* End of Outputs for SubSystem: '<S18>/Stationary reference frame' */
  rtDW.ib[3] = 13.731987951966302 * isb_d;

  /* S-Function (sfun_spssw_discc): '<S63>/State-Space' incorporates:
   *  Constant: '<S39>/DC'
   */

  /* S-Function block: <S63>/State-Space */
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

  /* Gain: '<S17>/1_Vb' */
  rtb_u_Vb_idx_3 = 0.0055670221426890416 * rtDW.StateSpace_o1[7];

  /* Outputs for Enabled SubSystem: '<S17>/Rotor reference frame' incorporates:
   *  EnablePort: '<S28>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE_o) {
    /* Disable for Fcn: '<S28>/vqr' incorporates:
     *  Outport: '<S28>/vqr,vdr'
     */
    rtDW.vqr_p = 0.0;

    /* Disable for Fcn: '<S28>/vdr' incorporates:
     *  Outport: '<S28>/vqr,vdr'
     */
    rtDW.vdr_i = 0.0;

    /* Disable for Fcn: '<S28>/vqs' incorporates:
     *  Outport: '<S28>/vqs,vds'
     */
    rtDW.vqs_g = 0.0;

    /* Disable for Fcn: '<S28>/vds' incorporates:
     *  Outport: '<S28>/vqs,vds'
     */
    rtDW.vds_n = 0.0;
    rtDW.Rotorreferenceframe_MODE_o = false;
  }

  /* End of Outputs for SubSystem: '<S17>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S17>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S30>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MOD_a) {
    /* Disable for Fcn: '<S30>/vqr' incorporates:
     *  Outport: '<S30>/vqr,vdr'
     */
    rtDW.vqr = 0.0;

    /* Disable for Fcn: '<S30>/vdr' incorporates:
     *  Outport: '<S30>/vqr,vdr'
     */
    rtDW.vdr = 0.0;

    /* Disable for Fcn: '<S30>/vqs' incorporates:
     *  Outport: '<S30>/vqs,vds'
     */
    rtDW.vqs = 0.0;

    /* Disable for Fcn: '<S30>/vds' incorporates:
     *  Outport: '<S30>/vqs,vds'
     */
    rtDW.vds = 0.0;
    rtDW.Synchronousreferenceframe_MOD_a = false;
  }

  /* End of Outputs for SubSystem: '<S17>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S17>/Stationary reference frame' incorporates:
   *  EnablePort: '<S29>/Enable'
   */
  /* MultiPortSwitch: '<S17>/Multiport Switch' incorporates:
   *  Fcn: '<S29>/vdr'
   *  Fcn: '<S29>/vqr'
   *  MultiPortSwitch: '<S19>/Multiport Switch'
   */
  rtb_MultiportSwitch_p[0] = ((TrigonometricFunction_o2 - 1.7320508075688772 *
    TrigonometricFunction_o1) * 0.0 + 2.0 * TrigonometricFunction_o2 * 0.0) *
    0.33333333333333331;
  rtb_MultiportSwitch_p[1] = ((-TrigonometricFunction_o1 - 1.7320508075688772 *
    TrigonometricFunction_o2) * 0.0 + -2.0 * TrigonometricFunction_o1 * 0.0) *
    0.33333333333333331;

  /* MultiPortSwitch: '<S17>/Multiport Switch1' incorporates:
   *  Fcn: '<S29>/vds'
   *  Fcn: '<S29>/vqs'
   *  Gain: '<S17>/1_Vb'
   */
  rtb_MultiportSwitch1[0] = (0.0055670221426890416 * rtDW.StateSpace_o1[6] * 2.0
    + rtb_u_Vb_idx_3) * 0.33333333333333331;
  rtb_MultiportSwitch1[1] = -0.57735026918962573 * rtb_u_Vb_idx_3;

  /* End of Outputs for SubSystem: '<S17>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S18>/Stationary reference frame' incorporates:
   *  EnablePort: '<S33>/Enable'
   */
  /* Gain: '<S13>/unit conversion' incorporates:
   *  Fcn: '<S33>/isa'
   *  MultiPortSwitch: '<S18>/Multiport Switch1'
   *  Sum: '<S18>/Sum3'
   */
  rtb_u_Vb_idx_3 = ((0.0 - rtb_Sum2[0]) - isb_d) * 13.731987951966302;

  /* End of Outputs for SubSystem: '<S18>/Stationary reference frame' */
  TrigonometricFunction_o1 = 0.57177765423802906 * rtb_xk1[2];

  /* Outputs for Enabled SubSystem: '<S18>/Stationary reference frame' incorporates:
   *  EnablePort: '<S33>/Enable'
   */
  rtb_unitconversion_idx_9 = 13.731987951966302 * rtb_Sum2[0];

  /* End of Outputs for SubSystem: '<S18>/Stationary reference frame' */
  TrigonometricFunction_o2 = 0.57177765423802906 * rtb_xk1[3];
  isb_d *= 13.731987951966302;

  /* Outport: '<Root>/is_a' */
  rtY.is_a = rtb_unitconversion_idx_9;

  /* Outport: '<Root>/is_b' */
  rtY.is_b = isb_d;

  /* Outport: '<Root>/is_c' */
  rtY.is_c = rtb_u_Vb_idx_3;

  /* Outport: '<Root>/v_d' */
  rtY.v_d = TrigonometricFunction_o2;

  /* Outport: '<Root>/v_q' */
  rtY.v_q = TrigonometricFunction_o1;

  /* Switch: '<S23>/IC' incorporates:
   *  DigitalClock: '<S23>/Digital Clock'
   *  Gain: '<S27>/wbase*Ts//2 '
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Sum: '<S23>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
   */
  if (((rtM->Timing.clockTick1) * 0.0001) >= 0.0001) {
    for (i = 0; i <= 14; i += 2) {
      tmp_5 = _mm_loadu_pd(&W43wr[i]);
      tmp_5 = _mm_mul_pd(_mm_sub_pd(_mm_sub_pd(_mm_set1_pd(0.0), tmp_5),
        _mm_loadu_pd(&rtConstP.Constant4_Value[i])), _mm_set1_pd
                         (0.015707963267948967));
      tmp_6 = _mm_loadu_pd(&rtConstP.u5_Value_m[i]);
      _mm_storeu_pd(&rtb_Sum5[i], _mm_add_pd(tmp_6, tmp_5));
      _mm_storeu_pd(&rtb_Lminrows24col24_0[i], _mm_sub_pd(tmp_6, tmp_5));
    }

    /* Product: '<S27>/inversion' incorporates:
     *  Assignment: '<S37>/W(4,3)=wr'
     *  Constant: '<S16>/Constant4'
     *  Constant: '<S27>/u5'
     *  Gain: '<S27>/wbase*Ts//2'
     *  Gain: '<S27>/wbase*Ts//2 '
     *  MultiPortSwitch: '<S19>/Multiport Switch1'
     *  Sum: '<S16>/Sum1'
     *  Sum: '<S27>/Sum1'
     *  Sum: '<S27>/Sum5'
     *  Switch: '<S16>/Switch1'
     */
    rt_invd4x4_snf(rtb_Lminrows24col24_0, W43wr);

    /* Product: '<S27>/Product4' incorporates:
     *  Gain: '<S27>/wbase*Ts//2 '
     *  Sum: '<S27>/Sum5'
     */
    for (i = 0; i < 4; i++) {
      rtb_Sum5_tmp = i << 2;
      rtb_Product2_e = rtb_Sum5[rtb_Sum5_tmp + 1];
      rtb_DiscreteTimeIntegrator = rtb_Sum5[rtb_Sum5_tmp];
      rtb_Gain2_pc = rtb_Sum5[rtb_Sum5_tmp + 2];
      rtb_MultiportSwitch_h_idx_0 = rtb_Sum5[rtb_Sum5_tmp + 3];
      for (i_0 = 0; i_0 <= 2; i_0 += 2) {
        tmp_5 = _mm_loadu_pd(&W43wr[i_0 + 4]);
        tmp_6 = _mm_loadu_pd(&W43wr[i_0]);
        tmp_7 = _mm_loadu_pd(&W43wr[i_0 + 8]);
        tmp_8 = _mm_loadu_pd(&W43wr[i_0 + 12]);
        _mm_storeu_pd(&rtb_Lminrows24col24_0[i_0 + rtb_Sum5_tmp], _mm_add_pd
                      (_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd
          (rtb_Product2_e), tmp_5), _mm_mul_pd(_mm_set1_pd
          (rtb_DiscreteTimeIntegrator), tmp_6)), _mm_mul_pd(_mm_set1_pd
          (rtb_Gain2_pc), tmp_7)), _mm_mul_pd(_mm_set1_pd
          (rtb_MultiportSwitch_h_idx_0), tmp_8)));
      }
    }

    /* End of Product: '<S27>/Product4' */

    /* Sum: '<S23>/sum' incorporates:
     *  UnitDelay: '<S23>/voltages'
     */
    rtb_Switch_idx_0 = rtb_MultiportSwitch1[0] + rtDW.voltages_DSTATE[0];
    rtb_Subtract2 = rtb_MultiportSwitch_p[0] + rtDW.voltages_DSTATE[2];
    rtb_MultiportSwitch2 = rtb_MultiportSwitch1[1] + rtDW.voltages_DSTATE[1];
    tmp = rtb_MultiportSwitch_p[1] + rtDW.voltages_DSTATE[3];

    /* Product: '<S23>/Product2' */
    rtb_Product2_e = rtb_xk1[1];
    rtb_DiscreteTimeIntegrator = rtb_xk1[0];
    rtb_Gain2_pc = rtb_xk1[2];
    rtb_MultiportSwitch_h_idx_0 = rtb_xk1[3];
    for (i = 0; i <= 2; i += 2) {
      tmp_5 = _mm_loadu_pd(&W43wr[i + 4]);
      tmp_6 = _mm_set1_pd(0.015707963267948967);
      tmp_7 = _mm_loadu_pd(&W43wr[i]);
      tmp_8 = _mm_loadu_pd(&W43wr[i + 8]);
      tmp_0 = _mm_loadu_pd(&W43wr[i + 12]);
      tmp_1 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 4]);
      tmp_2 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i]);
      tmp_3 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 8]);
      tmp_4 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 12]);
      _mm_storeu_pd(&rtb_IC[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd
        (_mm_mul_pd(_mm_mul_pd(tmp_5, tmp_6), _mm_set1_pd(rtb_MultiportSwitch2)),
         _mm_mul_pd(_mm_mul_pd(tmp_6, tmp_7), _mm_set1_pd(rtb_Switch_idx_0))),
        _mm_mul_pd(_mm_mul_pd(tmp_8, tmp_6), _mm_set1_pd(rtb_Subtract2))),
        _mm_mul_pd(_mm_mul_pd(tmp_0, tmp_6), _mm_set1_pd(tmp))), _mm_add_pd
        (_mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd(rtb_Product2_e)),
        _mm_mul_pd(tmp_2, _mm_set1_pd(rtb_DiscreteTimeIntegrator))), _mm_mul_pd
                    (tmp_3, _mm_set1_pd(rtb_Gain2_pc))), _mm_mul_pd(tmp_4,
        _mm_set1_pd(rtb_MultiportSwitch_h_idx_0)))));
    }
  } else {
    /* Switch: '<S23>/IC' */
    rtb_IC[0] = rtb_xk1[0];
    rtb_IC[1] = rtb_xk1[1];
    rtb_IC[2] = rtb_xk1[2];
    rtb_IC[3] = rtb_xk1[3];
  }

  /* End of Switch: '<S23>/IC' */
  for (i = 0; i <= 0; i += 2) {
    /* Gain: '<S53>/Gain1' incorporates:
     *  Gain: '<S53>/Gain3'
     */
    _mm_storeu_pd(&rtb_MathFunction[i], _mm_mul_pd(_mm_add_pd(_mm_add_pd
      (_mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i + 3]), _mm_set1_pd(isb_d)),
       _mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i]), _mm_set1_pd
                  (rtb_unitconversion_idx_9))), _mm_mul_pd(_mm_loadu_pd
      (&rtConstP.Gain3_Gain[i + 6]), _mm_set1_pd(rtb_u_Vb_idx_3))), _mm_set1_pd
      (0.66666666666666663)));
  }

  for (i = 2; i < 3; i++) {
    /* Gain: '<S53>/Gain1' incorporates:
     *  Gain: '<S53>/Gain3'
     */
    rtb_MathFunction[i] = ((rtConstP.Gain3_Gain[i + 3] * isb_d +
      rtConstP.Gain3_Gain[i] * rtb_unitconversion_idx_9) + rtConstP.Gain3_Gain[i
      + 6] * rtb_u_Vb_idx_3) * 0.66666666666666663;
  }

  /* Gain: '<S15>/1\p1' */
  isb_d = 157.07963267948966 * rtb_Phisat;

  /* Outport: '<Root>/wm' */
  rtY.wm = isb_d;

  /* ComplexToMagnitudeAngle: '<S1>/Complex to Magnitude-Angle' incorporates:
   *  RealImagToComplex: '<S1>/Real-Imag to Complex'
   */
  rtb_u_Vb_idx_3 = rt_atan2d_snf(TrigonometricFunction_o1,
    TrigonometricFunction_o2);

  /* Outputs for Enabled SubSystem: '<S58>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S61>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S52>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S56>/Enable'
   */
  /* Fcn: '<S56>/Fcn' incorporates:
   *  Fcn: '<S56>/Fcn1'
   *  Fcn: '<S61>/Fcn'
   */
  rtb_Product2_e = cos(rtb_u_Vb_idx_3);
  rtb_DiscreteTimeIntegrator = sin(rtb_u_Vb_idx_3);

  /* End of Outputs for SubSystem: '<S58>/Subsystem - pi//2 delay' */

  /* Switch: '<S52>/Switch' incorporates:
   *  Fcn: '<S56>/Fcn'
   *  Fcn: '<S56>/Fcn1'
   */
  rtb_MultiportSwitch_h_idx_0 = rtb_MathFunction[0] * rtb_DiscreteTimeIntegrator
    - rtb_MathFunction[1] * rtb_Product2_e;
  rtb_unitconversion_idx_9 = rtb_MathFunction[0] * rtb_Product2_e +
    rtb_MathFunction[1] * rtb_DiscreteTimeIntegrator;

  /* End of Outputs for SubSystem: '<S52>/Subsystem - pi//2 delay' */

  /* Outputs for Enabled SubSystem: '<S58>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S61>/Enable'
   */
  /* Switch: '<S58>/Switch' incorporates:
   *  Fcn: '<S61>/Fcn'
   */
  rtb_Switch_idx_0 = TrigonometricFunction_o2 * rtb_DiscreteTimeIntegrator -
    TrigonometricFunction_o1 * rtb_Product2_e;

  /* End of Outputs for SubSystem: '<S58>/Subsystem - pi//2 delay' */

  /* Gain: '<S11>/Gain1' */
  TrigonometricFunction_o1 = 2.0 * isb_d;

  /* Clock: '<S11>/Clock' incorporates:
   *  Clock: '<S49>/Clock'
   */
  TrigonometricFunction_o2 = rtM->Timing.t[0];

  /* Switch: '<S11>/Switch' incorporates:
   *  Clock: '<S11>/Clock'
   *  Product: '<S11>/Divide'
   */
  if (TrigonometricFunction_o2 > 0.001) {
    rtb_Product2_e = rtb_unitconversion_idx_9 / rtb_MultiportSwitch_h_idx_0;
  } else {
    rtb_Product2_e = 0.0;
  }

  /* Sum: '<S11>/Add' incorporates:
   *  Gain: '<S11>/Gain'
   *  Switch: '<S11>/Switch'
   */
  rtb_Product2_e = 6.23416784551107 * rtb_Product2_e + TrigonometricFunction_o1;

  /* Sum: '<S2>/Sum1' incorporates:
   *  Constant: '<S1>/Constant'
   *  DiscreteIntegrator: '<S2>/Discrete-Time Integrator'
   *  Gain: '<S11>/Gain2'
   *  Gain: '<S11>/Gain3'
   *  Gain: '<S11>/Gain4'
   *  Gain: '<S11>/Gain5'
   *  Gain: '<S2>/Gain'
   *  Product: '<S11>/Product'
   *  Sum: '<S11>/Add2'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S2>/Sum'
   */
  rtb_DiscreteTimeIntegrator = ((0.0 - rtb_Product2_e * rtb_unitconversion_idx_9
    * 0.058148172073290927 * 0.060793999999999994) - 0.022361 * rtb_Switch_idx_0
    * 270.56932718377169) + ((9.4040441629117311 - rtb_MultiportSwitch_h_idx_0) *
    1.5703485875928103 + rtDW.DiscreteTimeIntegrator_DSTATE);

  /* Saturate: '<S2>/Saturation' */
  if (rtb_DiscreteTimeIntegrator > 359.2584956081995) {
    rtb_Gain2_pc = 359.2584956081995;
  } else if (rtb_DiscreteTimeIntegrator < -359.2584956081995) {
    rtb_Gain2_pc = -359.2584956081995;
  } else {
    rtb_Gain2_pc = rtb_DiscreteTimeIntegrator;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* Gain: '<S2>/Gain1' incorporates:
   *  Constant: '<S1>/Constant'
   *  Gain: '<S2>/Gain2'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S2>/Sum2'
   *  Sum: '<S2>/Sum3'
   */
  rtDW.Gain1_m = ((9.4040441629117311 - rtb_MultiportSwitch_h_idx_0) -
                  (rtb_DiscreteTimeIntegrator - rtb_Gain2_pc) *
                  0.636801285969825) * 289.61527294964543;

  /* Product: '<S11>/Product1' */
  rtb_MultiportSwitch_h_idx_0 *= rtb_Product2_e;

  /* Sum: '<S1>/Sum' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Inport: '<Root>/wish_speed'
   */
  rtb_DiscreteTimeIntegrator = 0.10471975511965977 * rtU.wish_speed - isb_d;

  /* Sum: '<S4>/Sum' incorporates:
   *  DiscreteIntegrator: '<S4>/Discrete-Time Integrator'
   *  Gain: '<S4>/Gain'
   */
  rtb_Product2_e = 1.3859596918105916 * rtb_DiscreteTimeIntegrator +
    rtDW.DiscreteTimeIntegrator_DSTATE_f;

  /* Saturate: '<S4>/Saturation' */
  if (rtb_Product2_e > 48.400544337535294) {
    isb_d = 48.400544337535294;
  } else if (rtb_Product2_e < -48.400544337535294) {
    isb_d = -48.400544337535294;
  } else {
    isb_d = rtb_Product2_e;
  }

  /* End of Saturate: '<S4>/Saturation' */

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S4>/Gain2'
   */
  rtb_unitconversion_idx_9 = 0.61904294884986732 * isb_d -
    rtb_unitconversion_idx_9;

  /* Sum: '<S3>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S3>/Discrete-Time Integrator'
   *  Gain: '<S11>/Gain6'
   *  Gain: '<S11>/Gain7'
   *  Gain: '<S11>/Gain8'
   *  Gain: '<S3>/Gain'
   *  Product: '<S11>/Product2'
   *  Sum: '<S11>/Add1'
   *  Sum: '<S3>/Sum'
   */
  TrigonometricFunction_o1 = (0.058148172073290927 * rtb_MultiportSwitch_h_idx_0
    * 0.060793999999999994 + rtb_Switch_idx_0 * TrigonometricFunction_o1 *
    0.97049050893180255) + (1.5703485875928103 * rtb_unitconversion_idx_9 +
    rtDW.DiscreteTimeIntegrator_DSTATE_p);

  /* Saturate: '<S3>/Saturation' */
  if (TrigonometricFunction_o1 > 359.2584956081995) {
    rtb_MultiportSwitch_h_idx_0 = 359.2584956081995;
  } else if (TrigonometricFunction_o1 < -359.2584956081995) {
    rtb_MultiportSwitch_h_idx_0 = -359.2584956081995;
  } else {
    rtb_MultiportSwitch_h_idx_0 = TrigonometricFunction_o1;
  }

  /* End of Saturate: '<S3>/Saturation' */

  /* Fcn: '<S7>/alpha' incorporates:
   *  Fcn: '<S7>/beta'
   */
  rtb_Switch_idx_0 = sin(1.5707963267948966 - rtb_u_Vb_idx_3);
  rtb_Subtract2 = cos(1.5707963267948966 - rtb_u_Vb_idx_3);
  rtb_u_Vb_idx_3 = rtb_Subtract2 * rtb_Gain2_pc + rtb_Switch_idx_0 *
    rtb_MultiportSwitch_h_idx_0;

  /* Gain: '<S40>/Gain' */
  rtb_MultiportSwitch2 = 1.7320508075688772 * rtb_u_Vb_idx_3;

  /* Fcn: '<S7>/beta' */
  rtb_Gain2_pc = -rtb_Switch_idx_0 * rtb_Gain2_pc + rtb_Subtract2 *
    rtb_MultiportSwitch_h_idx_0;

  /* Sum: '<S40>/Subtract2' incorporates:
   *  Constant: '<S46>/Constant'
   *  Constant: '<S47>/Constant'
   *  Constant: '<S48>/Constant'
   *  Gain: '<S40>/Gain2'
   *  Gain: '<S40>/Gain3'
   *  RelationalOperator: '<S46>/Compare'
   *  RelationalOperator: '<S47>/Compare'
   *  RelationalOperator: '<S48>/Compare'
   *  Sum: '<S40>/Subtract'
   *  Sum: '<S40>/Subtract1'
   */
  rtb_Subtract2_o = (uint8_T)(((uint32_T)((rtb_MultiportSwitch2 - rtb_Gain2_pc >
    0.0) << 1) + (uint32_T)(rtb_Gain2_pc > 0.0)) + (uint32_T)(((0.0 -
    rtb_MultiportSwitch2) - rtb_Gain2_pc > 0.0) << 2));

  /* Gain: '<S45>/Gain' */
  rtb_Switch_idx_0 = 3.4641016151377543E-6 * rtb_Gain2_pc;

  /* Gain: '<S45>/Gain2' */
  rtb_Gain2_pc *= 0.5;

  /* Gain: '<S45>/Gain1' */
  rtb_u_Vb_idx_3 *= 0.8660254037844386;

  /* MultiPortSwitch: '<S41>/Multiport Switch' incorporates:
   *  Gain: '<S41>/Gain'
   *  Gain: '<S41>/Gain1'
   *  Gain: '<S41>/Gain2'
   *  Gain: '<S45>/Gain3'
   *  Gain: '<S45>/Gain4'
   *  Sum: '<S45>/Subtract'
   *  Sum: '<S45>/Subtract1'
   */
  switch (rtb_Subtract2_o) {
   case 1:
    rtb_MultiportSwitch2 = (rtb_Gain2_pc - rtb_u_Vb_idx_3) *
      3.4641016151377543E-6;

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S45>/Gain3'
     *  Gain: '<S45>/Gain4'
     *  Sum: '<S45>/Subtract'
     *  Sum: '<S45>/Subtract1'
     */
    rtb_u_Vb_idx_3 = (rtb_Gain2_pc + rtb_u_Vb_idx_3) * 3.4641016151377543E-6;
    break;

   case 2:
    rtb_MultiportSwitch2 = (rtb_Gain2_pc + rtb_u_Vb_idx_3) *
      3.4641016151377543E-6;

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S41>/Gain'
     *  Gain: '<S45>/Gain3'
     *  Sum: '<S45>/Subtract'
     */
    rtb_u_Vb_idx_3 = -rtb_Switch_idx_0;
    break;

   case 3:
    rtb_MultiportSwitch2 = -((rtb_Gain2_pc - rtb_u_Vb_idx_3) *
      3.4641016151377543E-6);

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S41>/Gain2'
     *  Gain: '<S45>/Gain4'
     *  Sum: '<S45>/Subtract1'
     */
    rtb_u_Vb_idx_3 = rtb_Switch_idx_0;
    break;

   case 4:
    rtb_MultiportSwitch2 = -rtb_Switch_idx_0;

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S41>/Gain'
     *  Gain: '<S45>/Gain4'
     *  Sum: '<S45>/Subtract1'
     */
    rtb_u_Vb_idx_3 = (rtb_Gain2_pc - rtb_u_Vb_idx_3) * 3.4641016151377543E-6;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S41>/Gain1'
     *  Gain: '<S45>/Gain3'
     *  Sum: '<S45>/Subtract'
     */
    rtb_u_Vb_idx_3 = -((rtb_Gain2_pc + rtb_u_Vb_idx_3) * 3.4641016151377543E-6);
    break;

   default:
    rtb_MultiportSwitch2 = -((rtb_Gain2_pc + rtb_u_Vb_idx_3) *
      3.4641016151377543E-6);

    /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
     *  Gain: '<S41>/Gain1'
     *  Gain: '<S41>/Gain2'
     *  Gain: '<S45>/Gain3'
     *  Gain: '<S45>/Gain4'
     *  Sum: '<S45>/Subtract'
     *  Sum: '<S45>/Subtract1'
     */
    rtb_u_Vb_idx_3 = -((rtb_Gain2_pc - rtb_u_Vb_idx_3) * 3.4641016151377543E-6);
    break;
  }

  /* End of MultiPortSwitch: '<S41>/Multiport Switch' */

  /* Gain: '<S43>/Gain' incorporates:
   *  Constant: '<S43>/Constant'
   *  Sum: '<S43>/Subtract'
   */
  rtb_Gain2_pc = ((0.001 - rtb_MultiportSwitch2) - rtb_u_Vb_idx_3) * 0.25;

  /* Sum: '<S43>/Subtract1' incorporates:
   *  Gain: '<S43>/Gain1'
   */
  rtb_Switch_idx_0 = 0.5 * rtb_MultiportSwitch2 + rtb_Gain2_pc;

  /* Sum: '<S43>/Subtract2' incorporates:
   *  Gain: '<S43>/Gain2'
   */
  rtb_Subtract2 = 0.5 * rtb_u_Vb_idx_3 + rtb_Switch_idx_0;

  /* MultiPortSwitch: '<S42>/Multiport Switch1' */
  switch (rtb_Subtract2_o) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;
  }

  /* End of MultiPortSwitch: '<S42>/Multiport Switch1' */

  /* Lookup_n-D: '<S49>/Look-Up Table1' incorporates:
   *  Constant: '<S49>/Constant'
   *  Math: '<S49>/Math Function'
   */
  rtb_u_Vb_idx_3 = look1_binlx(rt_remd_snf(TrigonometricFunction_o2, 0.001),
    rtConstP.LookUpTable1_bp01Data, rtConstP.LookUpTable1_tableData, 2U);

  /* RelationalOperator: '<S44>/Relational Operator1' */
  rtb_LogicalOperator2 = (rtb_u_Vb_idx_3 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S44>/Cast To Double2' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[2] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S44>/Cast To Double3' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   *  Logic: '<S44>/Logical Operator1'
   */
  rtDW.VectorConcatenate[3] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S42>/Multiport Switch3' */
  switch (rtb_Subtract2_o) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;
  }

  /* End of MultiPortSwitch: '<S42>/Multiport Switch3' */

  /* RelationalOperator: '<S44>/Relational Operator' */
  rtb_LogicalOperator2 = (rtb_u_Vb_idx_3 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S44>/Cast To Double' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[0] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S44>/Cast To Double1' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   *  Logic: '<S44>/Logical Operator'
   */
  rtDW.VectorConcatenate[1] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S42>/Multiport Switch2' */
  switch (rtb_Subtract2_o) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Subtract2;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Switch_idx_0;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Gain2_pc;
    break;
  }

  /* End of MultiPortSwitch: '<S42>/Multiport Switch2' */

  /* RelationalOperator: '<S44>/Relational Operator2' */
  rtb_LogicalOperator2 = (rtb_u_Vb_idx_3 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S44>/Cast To Double4' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   */
  rtDW.VectorConcatenate[4] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S44>/Cast To Double5' incorporates:
   *  Concatenate: '<S44>/Vector Concatenate'
   *  Logic: '<S44>/Logical Operator2'
   */
  rtDW.VectorConcatenate[5] = !rtb_LogicalOperator2;

  /* Gain: '<S3>/Gain1' incorporates:
   *  Gain: '<S3>/Gain2'
   *  Sum: '<S3>/Sum2'
   *  Sum: '<S3>/Sum3'
   */
  rtDW.Gain1_k = (rtb_unitconversion_idx_9 - (TrigonometricFunction_o1 -
    rtb_MultiportSwitch_h_idx_0) * 0.636801285969825) * 289.61527294964543;

  /* Gain: '<S4>/Gain1' incorporates:
   *  Gain: '<S4>/Gain3'
   *  Sum: '<S4>/Sum1'
   *  Sum: '<S4>/Sum2'
   */
  rtDW.Gain1_f = (rtb_DiscreteTimeIntegrator - (rtb_Product2_e - isb_d) *
                  0.72152170507471169) * 24.626907165410849;

  /* Gain: '<S15>/1_2H' incorporates:
   *  Gain: '<S15>/F'
   *  Gain: '<S15>/Unit conversion'
   *  Gain: '<S20>/1-1'
   *  Inport: '<Root>/load_torque'
   *  Product: '<S20>/Mult1'
   *  Sum: '<S15>/Sum'
   *  Sum: '<S20>/Sum2'
   */
  rtb_u_2H = (((rtb_Sum2[0] * rtb_xk1[1] + rtb_xk1[0] * -rtb_Sum2[1]) -
               0.042453954778240446 * rtU.load_torque) - 0.0 * rtb_Phisat) *
    5.9506091980420592;

  /* DiscreteIntegrator: '<S15>/Rotor speed(wm)' */
  if (rtDW.Rotorspeedwm_SYSTEM_ENABLE != 0) {
    /* DiscreteIntegrator: '<S15>/Rotor speed(wm)' */
    rtDW.Rotorspeedwm = rtDW.Rotorspeedwm_DSTATE;
  } else {
    /* DiscreteIntegrator: '<S15>/Rotor speed(wm)' */
    rtDW.Rotorspeedwm = 5.0E-5 * rtb_u_2H + rtDW.Rotorspeedwm_DSTATE;
  }

  /* End of DiscreteIntegrator: '<S15>/Rotor speed(wm)' */

  /* Gain: '<S15>/web_psb' */
  rtb_web_psb = 314.15926535897933 * rtb_Phisat;

  /* Update for UnitDelay: '<S23>/fluxes' */
  rtDW.fluxes_DSTATE[0] = rtb_IC[0];

  /* Update for UnitDelay: '<S21>/fluxes' */
  rtDW.fluxes_DSTATE_d[0] = rtb_xk1[0];

  /* Update for UnitDelay: '<S23>/fluxes' */
  rtDW.fluxes_DSTATE[1] = rtb_IC[1];

  /* Update for UnitDelay: '<S21>/fluxes' */
  rtDW.fluxes_DSTATE_d[1] = rtb_xk1[1];

  /* Update for UnitDelay: '<S23>/fluxes' */
  rtDW.fluxes_DSTATE[2] = rtb_IC[2];

  /* Update for UnitDelay: '<S21>/fluxes' */
  rtDW.fluxes_DSTATE_d[2] = rtb_xk1[2];

  /* Update for UnitDelay: '<S23>/fluxes' */
  rtDW.fluxes_DSTATE[3] = rtb_IC[3];

  /* Update for UnitDelay: '<S21>/fluxes' */
  rtDW.fluxes_DSTATE_d[3] = rtb_xk1[3];

  /* Update for DiscreteIntegrator: '<S15>/Rotor angle thetam' */
  rtDW.Rotoranglethetam_DSTATE += 0.0001 * rtb_web_psb;

  /* Update for UnitDelay: '<S38>/wm_delay' */
  rtDW.wm_delay_DSTATE = rtDW.Rotorspeedwm;

  /* Update for UnitDelay: '<S38>/wm_predict' */
  rtDW.wm_predict_DSTATE = rtb_wm_delay;

  /* Update for S-Function (sfun_spssw_discc): '<S63>/State-Space' incorporates:
   *  Constant: '<S39>/DC'
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

  /* Update for UnitDelay: '<S23>/voltages' */
  rtDW.voltages_DSTATE[0] = rtb_MultiportSwitch1[0];
  rtDW.voltages_DSTATE[2] = rtb_MultiportSwitch_p[0];
  rtDW.voltages_DSTATE[1] = rtb_MultiportSwitch1[1];
  rtDW.voltages_DSTATE[3] = rtb_MultiportSwitch_p[1];

  /* Update for DiscreteIntegrator: '<S2>/Discrete-Time Integrator' */
  rtDW.DiscreteTimeIntegrator_DSTATE += 0.0001 * rtDW.Gain1_m;

  /* Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' */
  rtDW.DiscreteTimeIntegrator_DSTATE_p += 0.0001 * rtDW.Gain1_k;

  /* Update for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' */
  rtDW.DiscreteTimeIntegrator_DSTATE_f += 0.0001 * rtDW.Gain1_f;

  /* Update for DiscreteIntegrator: '<S15>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 0U;
  rtDW.Rotorspeedwm_DSTATE = 5.0E-5 * rtb_u_2H + rtDW.Rotorspeedwm;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   */
  rtM->Timing.t[0] =
    ((time_T)(++rtM->Timing.clockTick0)) * rtM->Timing.stepSize0;

  {
    /* Update absolute timer for sample time: [0.0001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.0001, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     */
    rtM->Timing.clockTick1++;
  }
}

/* Model initialize function */
void motor_model2_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&rtM->solverInfo, &rtM->Timing.simTimeStep);
    rtsiSetTPtr(&rtM->solverInfo, &rtmGetTPtr(rtM));
    rtsiSetStepSizePtr(&rtM->solverInfo, &rtM->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&rtM->solverInfo, (&rtmGetErrorStatus(rtM)));
    rtsiSetRTModelPtr(&rtM->solverInfo, rtM);
  }

  rtsiSetSimTimeStep(&rtM->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&rtM->solverInfo,"FixedStepDiscrete");
  rtmSetTPtr(rtM, &rtM->Timing.tArray[0]);
  rtM->Timing.stepSize0 = 0.0001;

  /* Start for S-Function (sfun_spssw_discc): '<S63>/State-Space' incorporates:
   *  Constant: '<S39>/DC'
   */

  /* S-Function block: <S63>/State-Space */
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

  /* InitializeConditions for S-Function (sfun_spssw_discc): '<S63>/State-Space' incorporates:
   *  Constant: '<S39>/DC'
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

  /* Enable for DiscreteIntegrator: '<S15>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 1U;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
