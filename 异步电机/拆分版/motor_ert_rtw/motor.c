/*
 * File: motor.c
 *
 * Code generated for Simulink model 'motor'.
 *
 * Model version                  : 1.0
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Thu Oct 10 15:56:04 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "motor.h"
#include "rtwtypes.h"
#include <string.h>
#include <emmintrin.h>
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

const real_T motor_RGND = 0.0;         /* real_T ground */

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
static real_T rtGetNaN(void);
static real32_T rtGetNaNF(void);

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

/* Model step function */
void motor_step(void)
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
  real_T TrigonometricFunction_o1;
  real_T TrigonometricFunction_o2;
  real_T fluxes_DSTATE;
  real_T fluxes_DSTATE_0;
  real_T isb_f;
  real_T rtb_MultiportSwitch1_idx_0;
  real_T rtb_MultiportSwitch1_idx_1;
  real_T rtb_MultiportSwitch1_idx_2;
  real_T rtb_MultiportSwitch1_idx_3;
  real_T rtb_Phisat;
  real_T rtb_phimd;
  real_T rtb_u_Vb_idx_3;
  real_T vdr_p;
  real_T vds_i;
  int32_T i;
  int32_T i_0;
  int32_T rtb_Sum5_tmp;

  /* Update for UnitDelay: '<S13>/fluxes' incorporates:
   *  UnitDelay: '<S15>/fluxes'
   */
  rtDW.fluxes_DSTATE_l[0] = rtDW.fluxes_DSTATE[0];
  rtDW.fluxes_DSTATE_l[1] = rtDW.fluxes_DSTATE[1];
  rtDW.fluxes_DSTATE_l[2] = rtDW.fluxes_DSTATE[2];
  rtDW.fluxes_DSTATE_l[3] = rtDW.fluxes_DSTATE[3];

  /* Switch: '<S8>/Switch' incorporates:
   *  Constant: '<S8>/Constant2'
   *  Product: '<S14>/inversion'
   */
  memcpy(&rtb_Lminrows24col24_0[0], &rtConstP.Constant2_Value[0], sizeof(real_T)
         << 4U);

  /* UnitDelay: '<S15>/fluxes' incorporates:
   *  Product: '<S8>/Product3'
   */
  isb_f = rtDW.fluxes_DSTATE[1];
  TrigonometricFunction_o2 = rtDW.fluxes_DSTATE[0];
  fluxes_DSTATE = rtDW.fluxes_DSTATE[2];
  fluxes_DSTATE_0 = rtDW.fluxes_DSTATE[3];
  for (i = 0; i <= 2; i += 2) {
    /* Product: '<S8>/Product3' incorporates:
     *  UnitDelay: '<S15>/fluxes'
     */
    tmp_4 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 4]);
    tmp_5 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i]);
    tmp_6 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 8]);
    tmp_7 = _mm_loadu_pd(&rtb_Lminrows24col24_0[i + 12]);
    _mm_storeu_pd(&rtb_Sum2[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_4, _mm_set1_pd(isb_f)), _mm_mul_pd(tmp_5, _mm_set1_pd
      (TrigonometricFunction_o2))), _mm_mul_pd(tmp_6, _mm_set1_pd(fluxes_DSTATE))),
      _mm_mul_pd(tmp_7, _mm_set1_pd(fluxes_DSTATE_0))));
  }

  /* UnitDelay: '<S30>/wm_delay' */
  rtb_Phisat = rtDW.wm_delay_DSTATE;

  /* Sum: '<S30>/Sum1' incorporates:
   *  Gain: '<S30>/F2'
   *  UnitDelay: '<S30>/wm_predict'
   */
  rtb_phimd = 2.0 * rtb_Phisat - rtDW.wm_predict_DSTATE;

  /* Outputs for Enabled SubSystem: '<S11>/sin(thr),cos(thr)' incorporates:
   *  EnablePort: '<S28>/Enable'
   */
  if (rtDW.sinthrcosthr_MODE) {
    /* Disable for Trigonometry: '<S28>/Trigonometric Function' incorporates:
     *  Outport: '<S28>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o1_d = 0.0;

    /* Disable for Trigonometry: '<S28>/Trigonometric Function' incorporates:
     *  Outport: '<S28>/sin(thr),cos(thr)'
     */
    rtDW.TrigonometricFunction_o2_j = 0.0;

    /* Disable for Outport: '<S28>/sin(thr),cos(thr)' incorporates:
     *  Constant: '<S28>/Constant'
     */
    rtDW.Constant_e[0] = 0.0;
    rtDW.Constant_e[1] = 0.0;

    /* Disable for Assignment: '<S28>/W(2,1)=-wr' incorporates:
     *  Outport: '<S28>/W'
     */
    memset(&rtDW.W21wr[0], 0, sizeof(real_T) << 4U);
    rtDW.sinthrcosthr_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S11>/sin(thr),cos(thr)' */

  /* Outputs for Enabled SubSystem: '<S11>/sin(thr),cos(thr)1' incorporates:
   *  EnablePort: '<S29>/Enable'
   */
  /* Assignment: '<S29>/W(4,3)=wr' incorporates:
   *  SignalConversion generated from: '<S29>/W(3,4)=-wr'
   */
  memset(&W43wr[0], 0, sizeof(real_T) << 4U);

  /* Gain: '<S29>/Gain3' incorporates:
   *  Assignment: '<S29>/W(4,3)=wr'
   */
  W43wr[14] = -rtb_phimd;

  /* Assignment: '<S29>/W(4,3)=wr' */
  W43wr[11] = rtb_phimd;

  /* Trigonometry: '<S29>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S7>/Rotor angle thetam'
   */
  TrigonometricFunction_o1 = sin(rtDW.Rotoranglethetam_DSTATE);

  /* Trigonometry: '<S29>/Trigonometric Function' incorporates:
   *  DiscreteIntegrator: '<S7>/Rotor angle thetam'
   */
  TrigonometricFunction_o2 = cos(rtDW.Rotoranglethetam_DSTATE);

  /* End of Outputs for SubSystem: '<S11>/sin(thr),cos(thr)1' */

  /* Outputs for Enabled SubSystem: '<S10>/Rotor reference frame' incorporates:
   *  EnablePort: '<S24>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE) {
    /* Disable for Fcn: '<S24>/ira' incorporates:
     *  Outport: '<S24>/ira,irb'
     */
    rtDW.ira_p = 0.0;

    /* Disable for Fcn: '<S24>/irb' incorporates:
     *  Outport: '<S24>/ira,irb'
     */
    rtDW.irb_n = 0.0;

    /* Disable for Fcn: '<S24>/isa' incorporates:
     *  Outport: '<S24>/isa,isb'
     */
    rtDW.isa_m = 0.0;

    /* Disable for Fcn: '<S24>/isb' incorporates:
     *  Outport: '<S24>/isa,isb'
     */
    rtDW.isb_p = 0.0;
    rtDW.Rotorreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S10>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S10>/Stationary reference frame' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  /* Fcn: '<S25>/isb' */
  isb_f = -(1.7320508075688772 * rtb_Sum2[1] + rtb_Sum2[0]) / 2.0;

  /* End of Outputs for SubSystem: '<S10>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S10>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S26>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MODE) {
    /* Disable for Fcn: '<S26>/ira' incorporates:
     *  Outport: '<S26>/ira,irb'
     */
    rtDW.ira = 0.0;

    /* Disable for Fcn: '<S26>/irb' incorporates:
     *  Outport: '<S26>/ira,irb'
     */
    rtDW.irb = 0.0;

    /* Disable for Fcn: '<S26>/isa' incorporates:
     *  Outport: '<S26>/isa,isb'
     */
    rtDW.isa = 0.0;

    /* Disable for Fcn: '<S26>/isb' incorporates:
     *  Outport: '<S26>/isa,isb'
     */
    rtDW.isb = 0.0;
    rtDW.Synchronousreferenceframe_MODE = false;
  }

  /* End of Outputs for SubSystem: '<S10>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S10>/Stationary reference frame' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  /* Gain: '<S10>/ib' incorporates:
   *  Fcn: '<S25>/ira'
   *  Fcn: '<S25>/irb'
   *  Fcn: '<S25>/isa'
   *  MultiPortSwitch: '<S10>/Multiport Switch1'
   *  MultiPortSwitch: '<S11>/Multiport Switch'
   */
  rtDW.ib[0] = (TrigonometricFunction_o2 * rtb_Sum2[2] -
                TrigonometricFunction_o1 * rtb_Sum2[3]) * 13.731987951966302;
  rtDW.ib[2] = 13.731987951966302 * rtb_Sum2[0];
  rtDW.ib[1] = ((-TrigonometricFunction_o2 - 1.7320508075688772 *
                 TrigonometricFunction_o1) * rtb_Sum2[2] +
                (TrigonometricFunction_o1 - 1.7320508075688772 *
                 TrigonometricFunction_o2) * rtb_Sum2[3]) / 2.0 *
    13.731987951966302;

  /* End of Outputs for SubSystem: '<S10>/Stationary reference frame' */
  rtDW.ib[3] = 13.731987951966302 * isb_f;

  /* S-Function (sfun_spssw_discc): '<S34>/State-Space' incorporates:
   *  Constant: '<S31>/DC'
   *  Inport: '<Root>/gate'
   */

  /* S-Function block: <S34>/State-Space */
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

  /* Gain: '<S9>/1_Vb' */
  rtb_u_Vb_idx_3 = 0.0055670221426890416 * rtDW.StateSpace_o1[7];

  /* Outputs for Enabled SubSystem: '<S9>/Rotor reference frame' incorporates:
   *  EnablePort: '<S20>/Enable'
   */
  if (rtDW.Rotorreferenceframe_MODE_p) {
    /* Disable for Fcn: '<S20>/vqr' incorporates:
     *  Outport: '<S20>/vqr,vdr'
     */
    rtDW.vqr_j = 0.0;

    /* Disable for Fcn: '<S20>/vdr' incorporates:
     *  Outport: '<S20>/vqr,vdr'
     */
    rtDW.vdr_l = 0.0;

    /* Disable for Fcn: '<S20>/vqs' incorporates:
     *  Outport: '<S20>/vqs,vds'
     */
    rtDW.vqs_l = 0.0;

    /* Disable for Fcn: '<S20>/vds' incorporates:
     *  Outport: '<S20>/vqs,vds'
     */
    rtDW.vds_f = 0.0;
    rtDW.Rotorreferenceframe_MODE_p = false;
  }

  /* End of Outputs for SubSystem: '<S9>/Rotor reference frame' */

  /* Outputs for Enabled SubSystem: '<S9>/Stationary reference frame' incorporates:
   *  EnablePort: '<S21>/Enable'
   */
  /* Fcn: '<S21>/vdr' incorporates:
   *  MultiPortSwitch: '<S11>/Multiport Switch'
   */
  vdr_p = ((-TrigonometricFunction_o1 - 1.7320508075688772 *
            TrigonometricFunction_o2) * 0.0 + -2.0 * TrigonometricFunction_o1 *
           0.0) * 0.33333333333333331;

  /* Fcn: '<S21>/vds' */
  vds_i = -0.57735026918962573 * rtb_u_Vb_idx_3;

  /* Fcn: '<S21>/vqr' incorporates:
   *  MultiPortSwitch: '<S11>/Multiport Switch'
   */
  TrigonometricFunction_o1 = ((TrigonometricFunction_o2 - 1.7320508075688772 *
    TrigonometricFunction_o1) * 0.0 + 2.0 * TrigonometricFunction_o2 * 0.0) *
    0.33333333333333331;

  /* Fcn: '<S21>/vqs' incorporates:
   *  Gain: '<S9>/1_Vb'
   */
  rtb_u_Vb_idx_3 = (0.0055670221426890416 * rtDW.StateSpace_o1[6] * 2.0 +
                    rtb_u_Vb_idx_3) * 0.33333333333333331;

  /* End of Outputs for SubSystem: '<S9>/Stationary reference frame' */

  /* Outputs for Enabled SubSystem: '<S9>/Synchronous reference frame' incorporates:
   *  EnablePort: '<S22>/Enable'
   */
  if (rtDW.Synchronousreferenceframe_MOD_e) {
    /* Disable for Fcn: '<S22>/vqr' incorporates:
     *  Outport: '<S22>/vqr,vdr'
     */
    rtDW.vqr = 0.0;

    /* Disable for Fcn: '<S22>/vdr' incorporates:
     *  Outport: '<S22>/vqr,vdr'
     */
    rtDW.vdr = 0.0;

    /* Disable for Fcn: '<S22>/vqs' incorporates:
     *  Outport: '<S22>/vqs,vds'
     */
    rtDW.vqs = 0.0;

    /* Disable for Fcn: '<S22>/vds' incorporates:
     *  Outport: '<S22>/vqs,vds'
     */
    rtDW.vds = 0.0;
    rtDW.Synchronousreferenceframe_MOD_e = false;
  }

  /* End of Outputs for SubSystem: '<S9>/Synchronous reference frame' */

  /* Outputs for Enabled SubSystem: '<S10>/Stationary reference frame' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  /* Outport: '<Root>/ia' incorporates:
   *  Fcn: '<S25>/isa'
   *  Gain: '<S5>/unit conversion'
   */
  rtY.ia = 13.731987951966302 * rtb_Sum2[0];

  /* End of Outputs for SubSystem: '<S10>/Stationary reference frame' */

  /* Outport: '<Root>/ib' incorporates:
   *  Gain: '<S5>/unit conversion'
   *  MultiPortSwitch: '<S10>/Multiport Switch1'
   */
  rtY.ib = 13.731987951966302 * isb_f;

  /* Outputs for Enabled SubSystem: '<S10>/Stationary reference frame' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  /* Outport: '<Root>/ic' incorporates:
   *  Fcn: '<S25>/isa'
   *  Gain: '<S5>/unit conversion'
   *  MultiPortSwitch: '<S10>/Multiport Switch1'
   *  Sum: '<S10>/Sum3'
   */
  rtY.ic = ((0.0 - rtb_Sum2[0]) - isb_f) * 13.731987951966302;

  /* End of Outputs for SubSystem: '<S10>/Stationary reference frame' */

  /* Outport: '<Root>/vd' incorporates:
   *  Gain: '<S5>/unit conversion'
   *  UnitDelay: '<S15>/fluxes'
   */
  rtY.vd = 0.57177765423802906 * rtDW.fluxes_DSTATE[3];

  /* Outport: '<Root>/vq' incorporates:
   *  Gain: '<S5>/unit conversion'
   *  UnitDelay: '<S15>/fluxes'
   */
  rtY.vq = 0.57177765423802906 * rtDW.fluxes_DSTATE[2];

  /* Switch: '<S15>/IC' incorporates:
   *  DigitalClock: '<S15>/Digital Clock'
   *  Gain: '<S19>/wbase*Ts//2 '
   *  Product: '<S15>/Product1'
   *  Product: '<S15>/Product2'
   *  Sum: '<S15>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
   *  UnitDelay: '<S15>/fluxes'
   */
  if ((((rtM->Timing.clockTick0+rtM->Timing.clockTickH0* 4294967296.0)) * 1.0E-5)
      >= 1.0E-5) {
    for (i = 0; i <= 14; i += 2) {
      tmp_4 = _mm_loadu_pd(&W43wr[i]);
      tmp_4 = _mm_mul_pd(_mm_sub_pd(_mm_sub_pd(_mm_set1_pd(0.0), tmp_4),
        _mm_loadu_pd(&rtConstP.Constant4_Value[i])), _mm_set1_pd
                         (0.0015707963267948967));
      tmp_5 = _mm_loadu_pd(&rtConstP.u5_Value_l[i]);
      _mm_storeu_pd(&rtb_Sum5[i], _mm_add_pd(tmp_5, tmp_4));
      _mm_storeu_pd(&rtb_Lminrows24col24_0[i], _mm_sub_pd(tmp_5, tmp_4));
    }

    /* Product: '<S19>/inversion' incorporates:
     *  Assignment: '<S29>/W(4,3)=wr'
     *  Constant: '<S19>/u5'
     *  Constant: '<S8>/Constant4'
     *  Gain: '<S19>/wbase*Ts//2'
     *  Gain: '<S19>/wbase*Ts//2 '
     *  MultiPortSwitch: '<S11>/Multiport Switch1'
     *  Sum: '<S19>/Sum1'
     *  Sum: '<S19>/Sum5'
     *  Sum: '<S8>/Sum1'
     *  Switch: '<S8>/Switch1'
     */
    rt_invd4x4_snf(rtb_Lminrows24col24_0, W43wr);

    /* Product: '<S19>/Product4' incorporates:
     *  Gain: '<S19>/wbase*Ts//2 '
     *  Sum: '<S19>/Sum5'
     */
    for (i = 0; i < 4; i++) {
      rtb_Sum5_tmp = i << 2;
      isb_f = rtb_Sum5[rtb_Sum5_tmp + 1];
      TrigonometricFunction_o2 = rtb_Sum5[rtb_Sum5_tmp];
      fluxes_DSTATE = rtb_Sum5[rtb_Sum5_tmp + 2];
      fluxes_DSTATE_0 = rtb_Sum5[rtb_Sum5_tmp + 3];
      for (i_0 = 0; i_0 <= 2; i_0 += 2) {
        tmp_4 = _mm_loadu_pd(&W43wr[i_0 + 4]);
        tmp_5 = _mm_loadu_pd(&W43wr[i_0]);
        tmp_6 = _mm_loadu_pd(&W43wr[i_0 + 8]);
        tmp_7 = _mm_loadu_pd(&W43wr[i_0 + 12]);
        _mm_storeu_pd(&rtb_Lminrows24col24_0[i_0 + rtb_Sum5_tmp], _mm_add_pd
                      (_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd(isb_f),
          tmp_4), _mm_mul_pd(_mm_set1_pd(TrigonometricFunction_o2), tmp_5)),
          _mm_mul_pd(_mm_set1_pd(fluxes_DSTATE), tmp_6)), _mm_mul_pd(_mm_set1_pd
          (fluxes_DSTATE_0), tmp_7)));
      }
    }

    /* End of Product: '<S19>/Product4' */

    /* Sum: '<S15>/sum' incorporates:
     *  MultiPortSwitch: '<S9>/Multiport Switch'
     *  MultiPortSwitch: '<S9>/Multiport Switch1'
     *  UnitDelay: '<S15>/voltages'
     */
    rtb_MultiportSwitch1_idx_0 = rtb_u_Vb_idx_3 + rtDW.voltages_DSTATE[0];
    rtb_MultiportSwitch1_idx_2 = TrigonometricFunction_o1 +
      rtDW.voltages_DSTATE[2];
    rtb_MultiportSwitch1_idx_1 = vds_i + rtDW.voltages_DSTATE[1];
    rtb_MultiportSwitch1_idx_3 = vdr_p + rtDW.voltages_DSTATE[3];

    /* UnitDelay: '<S15>/fluxes' incorporates:
     *  Product: '<S15>/Product2'
     */
    isb_f = rtDW.fluxes_DSTATE[1];
    TrigonometricFunction_o2 = rtDW.fluxes_DSTATE[0];
    fluxes_DSTATE = rtDW.fluxes_DSTATE[2];
    fluxes_DSTATE_0 = rtDW.fluxes_DSTATE[3];
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
        (rtb_MultiportSwitch1_idx_1)), _mm_mul_pd(_mm_mul_pd(tmp_5, tmp_6),
        _mm_set1_pd(rtb_MultiportSwitch1_idx_0))), _mm_mul_pd(_mm_mul_pd(tmp_7,
        tmp_5), _mm_set1_pd(rtb_MultiportSwitch1_idx_2))), _mm_mul_pd(_mm_mul_pd
        (tmp, tmp_5), _mm_set1_pd(rtb_MultiportSwitch1_idx_3))), _mm_add_pd
        (_mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd(isb_f)), _mm_mul_pd
        (tmp_1, _mm_set1_pd(TrigonometricFunction_o2))), _mm_mul_pd(tmp_2,
        _mm_set1_pd(fluxes_DSTATE))), _mm_mul_pd(tmp_3, _mm_set1_pd
        (fluxes_DSTATE_0)))));
    }
  } else {
    rtb_MultiportSwitch[0] = rtDW.fluxes_DSTATE[0];
    rtb_MultiportSwitch[1] = rtDW.fluxes_DSTATE[1];
    rtb_MultiportSwitch[2] = rtDW.fluxes_DSTATE[2];
    rtb_MultiportSwitch[3] = rtDW.fluxes_DSTATE[3];
  }

  /* End of Switch: '<S15>/IC' */

  /* Outport: '<Root>/wm' incorporates:
   *  Gain: '<S7>/1\p1'
   */
  rtY.wm = 157.07963267948966 * rtb_phimd;

  /* Gain: '<S7>/1_2H' incorporates:
   *  Gain: '<S12>/1-1'
   *  Gain: '<S7>/F'
   *  Gain: '<S7>/Unit conversion'
   *  Inport: '<Root>/Tm'
   *  Product: '<S12>/Mult1'
   *  Sum: '<S12>/Sum2'
   *  Sum: '<S7>/Sum'
   *  UnitDelay: '<S15>/fluxes'
   */
  isb_f = (((rtb_Sum2[0] * rtDW.fluxes_DSTATE[1] + rtDW.fluxes_DSTATE[0] *
             -rtb_Sum2[1]) - 0.042453954778240446 * rtU.Tm) - 0.0 * rtb_phimd) *
    5.9506091980420592;

  /* DiscreteIntegrator: '<S7>/Rotor speed(wm)' */
  if (rtDW.Rotorspeedwm_SYSTEM_ENABLE != 0) {
    /* DiscreteIntegrator: '<S7>/Rotor speed(wm)' */
    TrigonometricFunction_o2 = rtDW.Rotorspeedwm_DSTATE;
  } else {
    /* DiscreteIntegrator: '<S7>/Rotor speed(wm)' */
    TrigonometricFunction_o2 = 5.0E-6 * isb_f + rtDW.Rotorspeedwm_DSTATE;
  }

  /* End of DiscreteIntegrator: '<S7>/Rotor speed(wm)' */

  /* Update for UnitDelay: '<S15>/fluxes' */
  rtDW.fluxes_DSTATE[0] = rtb_MultiportSwitch[0];
  rtDW.fluxes_DSTATE[1] = rtb_MultiportSwitch[1];
  rtDW.fluxes_DSTATE[2] = rtb_MultiportSwitch[2];
  rtDW.fluxes_DSTATE[3] = rtb_MultiportSwitch[3];

  /* Update for DiscreteIntegrator: '<S7>/Rotor angle thetam' incorporates:
   *  Gain: '<S7>/web_psb'
   */
  rtDW.Rotoranglethetam_DSTATE += 314.15926535897933 * rtb_phimd * 1.0E-5;

  /* Update for UnitDelay: '<S30>/wm_delay' */
  rtDW.wm_delay_DSTATE = TrigonometricFunction_o2;

  /* Update for UnitDelay: '<S30>/wm_predict' */
  rtDW.wm_predict_DSTATE = rtb_Phisat;

  /* Update for S-Function (sfun_spssw_discc): '<S34>/State-Space' incorporates:
   *  Constant: '<S31>/DC'
   *  Inport: '<Root>/gate'
   */
  {
    int_T *gState = (int_T*)rtDW.StateSpace_PWORK.G_STATE;

    /* Store switch gates values for next step */
    {
      int_T i1;
      const real_T *u1 = &rtU.gate[0];
      for (i1=0; i1 < 6; i1++) {
        *(gState++) = (int_T) u1[i1];
      }
    }
  }

  /* Update for UnitDelay: '<S15>/voltages' incorporates:
   *  MultiPortSwitch: '<S9>/Multiport Switch'
   *  MultiPortSwitch: '<S9>/Multiport Switch1'
   */
  rtDW.voltages_DSTATE[0] = rtb_u_Vb_idx_3;
  rtDW.voltages_DSTATE[2] = TrigonometricFunction_o1;
  rtDW.voltages_DSTATE[1] = vds_i;
  rtDW.voltages_DSTATE[3] = vdr_p;

  /* Update for DiscreteIntegrator: '<S7>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 0U;
  rtDW.Rotorspeedwm_DSTATE = 5.0E-6 * isb_f + TrigonometricFunction_o2;

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
}

/* Model initialize function */
void motor_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* Start for S-Function (sfun_spssw_discc): '<S34>/State-Space' incorporates:
   *  Constant: '<S31>/DC'
   *  Inport: '<Root>/gate'
   */

  /* S-Function block: <S34>/State-Space */
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

  /* InitializeConditions for S-Function (sfun_spssw_discc): '<S34>/State-Space' incorporates:
   *  Constant: '<S31>/DC'
   *  Inport: '<Root>/gate'
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

  /* Enable for DiscreteIntegrator: '<S7>/Rotor speed(wm)' */
  rtDW.Rotorspeedwm_SYSTEM_ENABLE = 1U;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
