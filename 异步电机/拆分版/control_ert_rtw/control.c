/*
 * File: control.c
 *
 * Code generated for Simulink model 'control'.
 *
 * Model version                  : 1.2
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Thu Oct 10 17:03:55 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "control.h"
#include "rtwtypes.h"
#include <emmintrin.h>
#include <math.h>
#include <float.h>
#include "rt_look.h"
#include <stddef.h>
#define NumBitsPerChar                 8U

/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Real-time model */
static RT_MODEL rtM_;
RT_MODEL *const rtM = &rtM_;
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
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
void control_step(void)
{
  real_T rtb_Gain1[3];
  real_T rtb_ComplextoMagnitudeAngle_o2;
  real_T rtb_Gain;
  real_T rtb_Gain1_ng;
  real_T rtb_Gain3_j;
  real_T rtb_Gain4;
  real_T rtb_MultiportSwitch2;
  real_T rtb_Saturation_i;
  real_T rtb_Sum_i;
  real_T rtb_Sum_m;
  real_T rtb_Switch_idx_0;
  real_T rtb_Switch_idx_0_tmp;
  real_T rtb_Switch_idx_1;
  real_T rtb_Switch_m_idx_0;
  int32_T i;
  uint8_T rtb_Subtract2_a;
  boolean_T rtb_LogicalOperator2;
  for (i = 0; i <= 0; i += 2) {
    /* Gain: '<S20>/Gain1' incorporates:
     *  Gain: '<S20>/Gain3'
     *  Inport: '<Root>/ia'
     *  Inport: '<Root>/ib'
     *  Inport: '<Root>/ic'
     */
    _mm_storeu_pd(&rtb_Gain1[i], _mm_mul_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (_mm_loadu_pd(&rtConstP.Gain3_Gain[i + 3]), _mm_set1_pd(rtU.ib)),
      _mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i]), _mm_set1_pd(rtU.ia))),
      _mm_mul_pd(_mm_loadu_pd(&rtConstP.Gain3_Gain[i + 6]), _mm_set1_pd(rtU.ic))),
      _mm_set1_pd(0.66666666666666663)));
  }

  for (i = 2; i < 3; i++) {
    /* Gain: '<S20>/Gain1' incorporates:
     *  Gain: '<S20>/Gain3'
     *  Inport: '<Root>/ia'
     *  Inport: '<Root>/ib'
     *  Inport: '<Root>/ic'
     */
    rtb_Gain1[i] = ((rtConstP.Gain3_Gain[i + 3] * rtU.ib + rtConstP.Gain3_Gain[i]
                     * rtU.ia) + rtConstP.Gain3_Gain[i + 6] * rtU.ic) *
      0.66666666666666663;
  }

  /* ComplexToMagnitudeAngle: '<Root>/Complex to Magnitude-Angle' incorporates:
   *  Inport: '<Root>/vd'
   *  Inport: '<Root>/vq'
   */
  rtb_ComplextoMagnitudeAngle_o2 = rt_atan2d_snf(rtU.vq, rtU.vd);

  /* Outputs for Enabled SubSystem: '<S25>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S28>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S19>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S23>/Enable'
   */
  /* Fcn: '<S23>/Fcn' incorporates:
   *  Fcn: '<S23>/Fcn1'
   *  Fcn: '<S28>/Fcn'
   */
  rtb_Switch_idx_0_tmp = cos(rtb_ComplextoMagnitudeAngle_o2);
  rtb_Saturation_i = sin(rtb_ComplextoMagnitudeAngle_o2);

  /* End of Outputs for SubSystem: '<S25>/Subsystem - pi//2 delay' */

  /* Switch: '<S19>/Switch' incorporates:
   *  Fcn: '<S23>/Fcn'
   *  Fcn: '<S23>/Fcn1'
   */
  rtb_Switch_idx_0 = rtb_Gain1[0] * rtb_Saturation_i - rtb_Gain1[1] *
    rtb_Switch_idx_0_tmp;
  rtb_Switch_idx_1 = rtb_Gain1[0] * rtb_Switch_idx_0_tmp + rtb_Gain1[1] *
    rtb_Saturation_i;

  /* End of Outputs for SubSystem: '<S19>/Subsystem - pi//2 delay' */

  /* Outputs for Enabled SubSystem: '<S25>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S28>/Enable'
   */
  /* Switch: '<S25>/Switch' incorporates:
   *  Fcn: '<S28>/Fcn'
   *  Inport: '<Root>/vd'
   *  Inport: '<Root>/vq'
   */
  rtb_Switch_m_idx_0 = rtU.vd * rtb_Saturation_i - rtU.vq * rtb_Switch_idx_0_tmp;

  /* End of Outputs for SubSystem: '<S25>/Subsystem - pi//2 delay' */

  /* DigitalClock: '<S7>/Clock' incorporates:
   *  DigitalClock: '<S18>/Digital Clock'
   */
  rtb_Switch_idx_0_tmp = ((rtM->Timing.clockTick0) * 1.0);

  /* Switch: '<S7>/Switch' incorporates:
   *  DigitalClock: '<S7>/Clock'
   *  Product: '<S7>/Divide'
   */
  if (rtb_Switch_idx_0_tmp > 0.0001) {
    rtb_Saturation_i = rtb_Switch_idx_1 / rtb_Switch_idx_0;
  } else {
    rtb_Saturation_i = 0.0;
  }

  /* End of Switch: '<S7>/Switch' */

  /* Gain: '<S7>/Gain' */
  rtb_Gain = 6.23416784551107 * rtb_Saturation_i;

  /* Gain: '<S7>/Gain1' incorporates:
   *  Inport: '<Root>/wm'
   */
  rtb_Saturation_i = 2.0 * rtU.wm;

  /* Sum: '<S7>/Add' */
  rtb_Gain += rtb_Saturation_i;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Constant: '<Root>/Constant'
   *  DiscreteIntegrator: '<S1>/Discrete-Time Integrator'
   *  Gain: '<S1>/Gain'
   *  Gain: '<S7>/Gain2'
   *  Gain: '<S7>/Gain3'
   *  Gain: '<S7>/Gain4'
   *  Gain: '<S7>/Gain5'
   *  Product: '<S7>/Product'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<S1>/Sum'
   *  Sum: '<S7>/Add2'
   */
  rtb_Sum_i = ((0.0 - rtb_Gain * rtb_Switch_idx_1 * 0.058148172073290927 *
                0.060793999999999994) - 0.022361 * rtb_Switch_m_idx_0 *
               270.56932718377169) + ((9.4040441629117311 - rtb_Switch_idx_0) *
    15.7034858759281 + rtDW.DiscreteTimeIntegrator_DSTATE);

  /* Saturate: '<S1>/Saturation' */
  if (rtb_Sum_i > 359.2584956081995) {
    rtb_Gain1_ng = 359.2584956081995;
  } else if (rtb_Sum_i < -359.2584956081995) {
    rtb_Gain1_ng = -359.2584956081995;
  } else {
    rtb_Gain1_ng = rtb_Sum_i;
  }

  /* End of Saturate: '<S1>/Saturation' */

  /* Sum: '<S1>/Sum2' */
  rtb_Sum_i -= rtb_Gain1_ng;

  /* Product: '<S7>/Product1' */
  rtb_Sum_m = rtb_Gain * rtb_Switch_idx_0;

  /* Gain: '<S7>/Gain8' incorporates:
   *  Product: '<S7>/Product2'
   */
  rtb_Switch_m_idx_0 = rtb_Switch_m_idx_0 * rtb_Saturation_i *
    0.97049050893180255;

  /* Sum: '<S3>/Sum' incorporates:
   *  DiscreteIntegrator: '<S3>/Discrete-Time Integrator'
   *  Gain: '<S3>/Gain'
   *  Inport: '<Root>/wm'
   *  Sum: '<Root>/Sum'
   */
  rtb_Gain = (52.359877559829883 - rtU.wm) * 13.859596918105916 +
    rtDW.DiscreteTimeIntegrator_DSTATE_g;

  /* Saturate: '<S3>/Saturation' */
  if (rtb_Gain > 48.400544337535294) {
    rtb_Saturation_i = 48.400544337535294;
  } else if (rtb_Gain < -48.400544337535294) {
    rtb_Saturation_i = -48.400544337535294;
  } else {
    rtb_Saturation_i = rtb_Gain;
  }

  /* End of Saturate: '<S3>/Saturation' */

  /* Sum: '<Root>/Sum1' incorporates:
   *  Gain: '<S3>/Gain2'
   */
  rtb_Switch_idx_1 = 0.61904294884986732 * rtb_Saturation_i - rtb_Switch_idx_1;

  /* Sum: '<S2>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S2>/Discrete-Time Integrator'
   *  Gain: '<S2>/Gain'
   *  Gain: '<S7>/Gain6'
   *  Gain: '<S7>/Gain7'
   *  Sum: '<S2>/Sum'
   *  Sum: '<S7>/Add1'
   */
  rtb_Sum_m = (0.058148172073290927 * rtb_Sum_m * 0.060793999999999994 +
               rtb_Switch_m_idx_0) + (15.7034858759281 * rtb_Switch_idx_1 +
    rtDW.DiscreteTimeIntegrator_DSTATE_p);

  /* Saturate: '<S2>/Saturation' */
  if (rtb_Sum_m > 359.2584956081995) {
    rtb_Switch_m_idx_0 = 359.2584956081995;
  } else if (rtb_Sum_m < -359.2584956081995) {
    rtb_Switch_m_idx_0 = -359.2584956081995;
  } else {
    rtb_Switch_m_idx_0 = rtb_Sum_m;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* Fcn: '<S4>/beta' incorporates:
   *  Fcn: '<S4>/alpha'
   */
  rtb_MultiportSwitch2 = cos(1.5707963267948966 - rtb_ComplextoMagnitudeAngle_o2);
  rtb_Gain4 = sin(1.5707963267948966 - rtb_ComplextoMagnitudeAngle_o2);
  rtb_ComplextoMagnitudeAngle_o2 = -rtb_Gain4 * rtb_Gain1_ng +
    rtb_MultiportSwitch2 * rtb_Switch_m_idx_0;

  /* Fcn: '<S4>/alpha' */
  rtb_Gain1_ng = rtb_MultiportSwitch2 * rtb_Gain1_ng + rtb_Gain4 *
    rtb_Switch_m_idx_0;

  /* Gain: '<S9>/Gain' */
  rtb_MultiportSwitch2 = 1.7320508075688772 * rtb_Gain1_ng;

  /* Sum: '<S9>/Subtract2' incorporates:
   *  Constant: '<S15>/Constant'
   *  Constant: '<S16>/Constant'
   *  Constant: '<S17>/Constant'
   *  Gain: '<S9>/Gain2'
   *  Gain: '<S9>/Gain3'
   *  RelationalOperator: '<S15>/Compare'
   *  RelationalOperator: '<S16>/Compare'
   *  RelationalOperator: '<S17>/Compare'
   *  Sum: '<S9>/Subtract'
   *  Sum: '<S9>/Subtract1'
   */
  rtb_Subtract2_a = (uint8_T)(((uint32_T)((rtb_MultiportSwitch2 -
    rtb_ComplextoMagnitudeAngle_o2 > 0.0) << 1) + (uint32_T)
    (rtb_ComplextoMagnitudeAngle_o2 > 0.0)) + (uint32_T)(((0.0 -
    rtb_MultiportSwitch2) - rtb_ComplextoMagnitudeAngle_o2 > 0.0) << 2));

  /* Gain: '<S14>/Gain2' */
  rtb_MultiportSwitch2 = 0.5 * rtb_ComplextoMagnitudeAngle_o2;

  /* Gain: '<S14>/Gain1' */
  rtb_Gain1_ng *= 0.8660254037844386;

  /* Gain: '<S14>/Gain4' incorporates:
   *  Sum: '<S14>/Subtract1'
   */
  rtb_Gain4 = (rtb_MultiportSwitch2 - rtb_Gain1_ng) * 3.4641016151377545E-7;

  /* Gain: '<S14>/Gain3' incorporates:
   *  Sum: '<S14>/Subtract'
   */
  rtb_Gain3_j = (rtb_MultiportSwitch2 + rtb_Gain1_ng) * 3.4641016151377545E-7;

  /* MultiPortSwitch: '<S10>/Multiport Switch' incorporates:
   *  Gain: '<S10>/Gain'
   *  Gain: '<S10>/Gain1'
   *  Gain: '<S10>/Gain2'
   *  Gain: '<S14>/Gain'
   *  Gain: '<S14>/Gain3'
   *  Gain: '<S14>/Gain4'
   *  Sum: '<S14>/Subtract'
   *  Sum: '<S14>/Subtract1'
   */
  switch (rtb_Subtract2_a) {
   case 1:
    rtb_MultiportSwitch2 = (rtb_MultiportSwitch2 - rtb_Gain1_ng) *
      3.4641016151377545E-7;

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S14>/Gain3'
     *  Gain: '<S14>/Gain4'
     *  Sum: '<S14>/Subtract1'
     */
    rtb_ComplextoMagnitudeAngle_o2 = rtb_Gain3_j;
    break;

   case 2:
    rtb_MultiportSwitch2 = (rtb_MultiportSwitch2 + rtb_Gain1_ng) *
      3.4641016151377545E-7;

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S10>/Gain'
     *  Gain: '<S14>/Gain'
     *  Gain: '<S14>/Gain3'
     *  Sum: '<S14>/Subtract'
     */
    rtb_ComplextoMagnitudeAngle_o2 = -(3.4641016151377545E-7 *
      rtb_ComplextoMagnitudeAngle_o2);
    break;

   case 3:
    rtb_MultiportSwitch2 = -((rtb_MultiportSwitch2 - rtb_Gain1_ng) *
      3.4641016151377545E-7);

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S10>/Gain2'
     *  Gain: '<S14>/Gain'
     *  Gain: '<S14>/Gain4'
     *  Sum: '<S14>/Subtract1'
     */
    rtb_ComplextoMagnitudeAngle_o2 *= 3.4641016151377545E-7;
    break;

   case 4:
    rtb_MultiportSwitch2 = -(3.4641016151377545E-7 *
      rtb_ComplextoMagnitudeAngle_o2);

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S10>/Gain'
     *  Gain: '<S14>/Gain'
     *  Gain: '<S14>/Gain4'
     */
    rtb_ComplextoMagnitudeAngle_o2 = rtb_Gain4;
    break;

   case 5:
    rtb_MultiportSwitch2 = 3.4641016151377545E-7 *
      rtb_ComplextoMagnitudeAngle_o2;

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S10>/Gain1'
     *  Gain: '<S14>/Gain'
     *  Gain: '<S14>/Gain3'
     */
    rtb_ComplextoMagnitudeAngle_o2 = -rtb_Gain3_j;
    break;

   default:
    rtb_MultiportSwitch2 = -((rtb_MultiportSwitch2 + rtb_Gain1_ng) *
      3.4641016151377545E-7);

    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Gain: '<S10>/Gain1'
     *  Gain: '<S10>/Gain2'
     *  Gain: '<S14>/Gain3'
     *  Gain: '<S14>/Gain4'
     *  Sum: '<S14>/Subtract'
     */
    rtb_ComplextoMagnitudeAngle_o2 = -rtb_Gain4;
    break;
  }

  /* End of MultiPortSwitch: '<S10>/Multiport Switch' */

  /* Gain: '<S12>/Gain' incorporates:
   *  Constant: '<S12>/Constant'
   *  Sum: '<S12>/Subtract'
   */
  rtb_Gain1_ng = ((0.0001 - rtb_MultiportSwitch2) -
                  rtb_ComplextoMagnitudeAngle_o2) * 0.25;

  /* Sum: '<S12>/Subtract1' incorporates:
   *  Gain: '<S12>/Gain1'
   */
  rtb_Gain4 = 0.5 * rtb_MultiportSwitch2 + rtb_Gain1_ng;

  /* Sum: '<S12>/Subtract2' incorporates:
   *  Gain: '<S12>/Gain2'
   */
  rtb_Gain3_j = 0.5 * rtb_ComplextoMagnitudeAngle_o2 + rtb_Gain4;

  /* MultiPortSwitch: '<S11>/Multiport Switch3' */
  switch (rtb_Subtract2_a) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;
  }

  /* End of MultiPortSwitch: '<S11>/Multiport Switch3' */

  /* Lookup: '<S18>/Look-Up Table' incorporates:
   *  Fcn: '<S18>/Fcn1'
   */
  rtb_ComplextoMagnitudeAngle_o2 = rt_Lookup(&rtConstP.LookUpTable_XData[0], 3,
    rt_remd_snf(rtb_Switch_idx_0_tmp, 0.0001), &rtConstP.LookUpTable_YData[0]);

  /* RelationalOperator: '<S13>/Relational Operator' */
  rtb_LogicalOperator2 = (rtb_ComplextoMagnitudeAngle_o2 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S13>/Cast To Double' incorporates:
   *  Outport: '<Root>/gate'
   */
  rtY.gate[0] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S13>/Cast To Double1' incorporates:
   *  Logic: '<S13>/Logical Operator'
   *  Outport: '<Root>/gate'
   */
  rtY.gate[1] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S11>/Multiport Switch1' */
  switch (rtb_Subtract2_a) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;
  }

  /* End of MultiPortSwitch: '<S11>/Multiport Switch1' */

  /* RelationalOperator: '<S13>/Relational Operator1' */
  rtb_LogicalOperator2 = (rtb_ComplextoMagnitudeAngle_o2 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S13>/Cast To Double2' incorporates:
   *  Outport: '<Root>/gate'
   */
  rtY.gate[2] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S13>/Cast To Double3' incorporates:
   *  Logic: '<S13>/Logical Operator1'
   *  Outport: '<Root>/gate'
   */
  rtY.gate[3] = !rtb_LogicalOperator2;

  /* MultiPortSwitch: '<S11>/Multiport Switch2' */
  switch (rtb_Subtract2_a) {
   case 1:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;

   case 2:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;

   case 3:
    rtb_MultiportSwitch2 = rtb_Gain3_j;
    break;

   case 4:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;

   case 5:
    rtb_MultiportSwitch2 = rtb_Gain4;
    break;

   default:
    rtb_MultiportSwitch2 = rtb_Gain1_ng;
    break;
  }

  /* End of MultiPortSwitch: '<S11>/Multiport Switch2' */

  /* RelationalOperator: '<S13>/Relational Operator2' */
  rtb_LogicalOperator2 = (rtb_ComplextoMagnitudeAngle_o2 > rtb_MultiportSwitch2);

  /* DataTypeConversion: '<S13>/Cast To Double4' incorporates:
   *  Outport: '<Root>/gate'
   */
  rtY.gate[4] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S13>/Cast To Double5' incorporates:
   *  Logic: '<S13>/Logical Operator2'
   *  Outport: '<Root>/gate'
   */
  rtY.gate[5] = !rtb_LogicalOperator2;

  /* Update for DiscreteIntegrator: '<S1>/Discrete-Time Integrator' incorporates:
   *  Constant: '<Root>/Constant'
   *  Gain: '<S1>/Gain1'
   *  Gain: '<S1>/Gain2'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<S1>/Sum3'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE += ((9.4040441629117311 - rtb_Switch_idx_0)
    - 0.063680128596982508 * rtb_Sum_i) * 2896.152729496454;

  /* Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S3>/Gain3'
   *  Inport: '<Root>/wm'
   *  Sum: '<Root>/Sum'
   *  Sum: '<S3>/Sum1'
   *  Sum: '<S3>/Sum2'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE_g += ((52.359877559829883 - rtU.wm) -
    (rtb_Gain - rtb_Saturation_i) * 0.072152170507471164) * 2462.6907165410848;

  /* Update for DiscreteIntegrator: '<S2>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S2>/Gain2'
   *  Sum: '<S2>/Sum2'
   *  Sum: '<S2>/Sum3'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE_p += (rtb_Switch_idx_1 - (rtb_Sum_m -
    rtb_Switch_m_idx_0) * 0.063680128596982508) * 2896.152729496454;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The resolution of this integer timer is 1.0, which is the step size
   * of the task. Size of "clockTick0" ensures timer will not overflow during the
   * application lifespan selected.
   */
  rtM->Timing.clockTick0++;
}

/* Model initialize function */
void control_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
