/*
 * Asynchronous_motor_2017a.c
 *
 * Code generation for model "Asynchronous_motor_2017a".
 *
 * Model version              : 13.12
 * Simulink Coder version : 9.9 (R2023a) 19-Nov-2022
 * C source code generated on : Mon Nov 11 18:52:56 2024
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Asynchronous_motor_2017a.h"
#include "rtwtypes.h"
#include <math.h>
#include <string.h>
#include "Asynchronous_motor_2017a_private.h"
#include "rt_nonfinite.h"
#include <float.h>
#include "rt_defines.h"

const real_T Asynchronous_motor_2017a_RGND = 0.0;/* real_T ground */

/* Block signals (default storage) */
B_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_B;

/* Continuous states */
X_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_X;

/* Block states (default storage) */
DW_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_DW;

/* Real-time model */
static RT_MODEL_Asynchronous_motor_2_T Asynchronous_motor_2017a_M_;
RT_MODEL_Asynchronous_motor_2_T *const Asynchronous_motor_2017a_M =
  &Asynchronous_motor_2017a_M_;
real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
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
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

real_T look1_pbinlxpw(real_T u0, const real_T bp0[], const real_T table[],
                      uint32_T prevIndex[], uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T bpIdx;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'on'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    bpIdx = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    uint32_T found;
    uint32_T iLeft;
    uint32_T iRght;

    /* Binary Search using Previous Index */
    bpIdx = prevIndex[0U];
    iLeft = 0U;
    iRght = maxIndex;
    found = 0U;
    while (found == 0U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx - 1U;
        bpIdx = ((bpIdx + iLeft) - 1U) >> 1U;
      } else if (u0 < bp0[bpIdx + 1U]) {
        found = 1U;
      } else {
        iLeft = bpIdx + 1U;
        bpIdx = ((bpIdx + iRght) + 1U) >> 1U;
      }
    }

    frac = (u0 - bp0[bpIdx]) / (bp0[bpIdx + 1U] - bp0[bpIdx]);
  } else {
    bpIdx = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  prevIndex[0U] = bpIdx;

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[bpIdx];
  return (table[bpIdx + 1U] - yL_0d0) * frac + yL_0d0;
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 3;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Asynchronous_motor_2017a_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  Asynchronous_motor_2017a_step();
  Asynchronous_motor_2017a_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  Asynchronous_motor_2017a_step();
  Asynchronous_motor_2017a_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else if (rtIsNaN(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
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
void Asynchronous_motor_2017a_step(void)
{
  /* local block i/o variables */
  real_T rtb_xk1[4];
  real_T rtb_wm_delay;
  real_T rtb_MultiportSwitch_p[2];
  real_T rtb_MultiportSwitch1[2];
  real_T rtb_IC[4];
  real_T rtb_web_psb;
  real_T rtb_unitconversion[19];
  real_T rtb_Lminrows24col24[16];
  real_T rtb_Lminrows24col24_0[16];
  real_T rtb_inversion_0[16];
  real_T rtb_Sum2[4];
  real_T rtb_MathFunction[3];
  real_T rtb_Delay;
  real_T rtb_MultiportSwitch1_l_idx_0;
  real_T rtb_MultiportSwitch1_l_idx_1;
  real_T rtb_MultiportSwitch_h_idx_0;
  real_T rtb_MultiportSwitch_idx_0;
  real_T rtb_MultiportSwitch_idx_1;
  real_T rtb_MultiportSwitch_idx_2;
  real_T rtb_Sum3_e;
  real_T rtb_Switch_idx_0;
  real_T rtb_Switch_idx_1;
  real_T rtb_phimd;
  real_T rtb_u_Vb_idx_0;
  real_T rtb_u_Vb_idx_1;
  real_T rtb_u_Vb_idx_2;
  real_T rtb_u_Vb_idx_3;
  int32_T Linv_tmp;
  int32_T i;
  int32_T i_0;
  uint8_T rtb_Compare;
  boolean_T rtb_LogicalOperator2;
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* set solver stop time */
    if (!(Asynchronous_motor_2017a_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Asynchronous_motor_2017a_M->solverInfo,
                            ((Asynchronous_motor_2017a_M->Timing.clockTickH0 + 1)
        * Asynchronous_motor_2017a_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Asynchronous_motor_2017a_M->solverInfo,
                            ((Asynchronous_motor_2017a_M->Timing.clockTick0 + 1)
        * Asynchronous_motor_2017a_M->Timing.stepSize0 +
        Asynchronous_motor_2017a_M->Timing.clockTickH0 *
        Asynchronous_motor_2017a_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Asynchronous_motor_2017a_M)) {
    Asynchronous_motor_2017a_M->Timing.t[0] = rtsiGetT
      (&Asynchronous_motor_2017a_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* UnitDelay: '<S22>/fluxes' */
    rtb_xk1[0] = Asynchronous_motor_2017a_DW.fluxes_DSTATE[0];
    rtb_xk1[1] = Asynchronous_motor_2017a_DW.fluxes_DSTATE[1];
    rtb_xk1[2] = Asynchronous_motor_2017a_DW.fluxes_DSTATE[2];
    rtb_xk1[3] = Asynchronous_motor_2017a_DW.fluxes_DSTATE[3];

    /* Outputs for Enabled SubSystem: '<S15>/Saturation' incorporates:
     *  EnablePort: '<S21>/Enable'
     */
    /* Constant: '<S15>/Constant' */
    if (Asynchronous_motor_2017a_P.Constant_Value_a2 > 0.0) {
      /* Math: '<S23>/Math Function1' incorporates:
       *  Constant: '<S23>/u2'
       *  Math: '<S23>/Math Function'
       *  SignalConversion generated from: '<S23>/Math Function'
       *  Sum: '<S23>/Sum2'
       *  UnitDelay: '<S21>/Delay'
       *
       * About '<S23>/Math Function1':
       *  Operator: reciprocal
       *
       * About '<S23>/Math Function':
       *  Operator: reciprocal
       *
       * About SignalConversion generated from '<S23>/Math Function':
       *  Operator: reciprocal
       */
      rtb_Delay = 1.0 / ((1.0 / Asynchronous_motor_2017a_P.u2_Value[0] + 1.0 /
                          Asynchronous_motor_2017a_P.u2_Value[1]) + 1.0 /
                         Asynchronous_motor_2017a_B.Switch_m);

      /* Product: '<S25>/Product2' incorporates:
       *  Constant: '<S25>/u1'
       */
      rtb_Switch_idx_0 = Asynchronous_motor_2017a_P.u1_Value[0] * rtb_Delay;
      rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.u1_Value[1] * rtb_Delay;

      /* Sum: '<S25>/Sum1' incorporates:
       *  Gain: '<S20>/Gain'
       *  Product: '<S25>/Product1'
       *  Sum: '<S20>/Sum2'
       *  UnitDelay: '<S20>/fluxes'
       *  UnitDelay: '<S22>/fluxes'
       */
      rtb_phimd = (Asynchronous_motor_2017a_P.Gain_Gain *
                   Asynchronous_motor_2017a_DW.fluxes_DSTATE[1] -
                   Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[1]) *
        rtb_Switch_idx_0;

      /* Math: '<S21>/Math Function' incorporates:
       *  Gain: '<S20>/Gain'
       *  Product: '<S25>/Product'
       *  Product: '<S25>/Product1'
       *  Sum: '<S20>/Sum2'
       *  Sum: '<S25>/Sum1'
       *  Sum: '<S25>/Sum2'
       *  UnitDelay: '<S20>/fluxes'
       *  UnitDelay: '<S22>/fluxes'
       */
      rtb_Switch_idx_1 = rt_hypotd_snf((Asynchronous_motor_2017a_P.Gain_Gain *
        Asynchronous_motor_2017a_DW.fluxes_DSTATE[0] -
        Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[0]) * rtb_Switch_idx_0 +
        (Asynchronous_motor_2017a_P.Gain_Gain *
         Asynchronous_motor_2017a_DW.fluxes_DSTATE[2] -
         Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[2]) * rtb_Switch_idx_1,
        (Asynchronous_motor_2017a_P.Gain_Gain *
         Asynchronous_motor_2017a_DW.fluxes_DSTATE[3] -
         Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[3]) * rtb_Switch_idx_1 +
        rtb_phimd);

      /* Lookup_n-D: '<S21>/1-D Lookup Table' incorporates:
       *  Math: '<S21>/Math Function'
       */
      rtb_phimd = look1_pbinlxpw(rtb_Switch_idx_1,
        Asynchronous_motor_2017a_P.uDLookupTable_bp01Data,
        Asynchronous_motor_2017a_P.uDLookupTable_tableData,
        &Asynchronous_motor_2017a_DW.m_bpIndex, 1U);

      /* Switch: '<S21>/Switch' */
      if (rtb_phimd != 0.0) {
        /* Switch: '<S21>/Switch' incorporates:
         *  Product: '<S21>/Product'
         */
        Asynchronous_motor_2017a_B.Switch_m = rtb_Switch_idx_1 / rtb_phimd;
      } else {
        /* Switch: '<S21>/Switch' incorporates:
         *  Constant: '<S21>/Constant1'
         */
        Asynchronous_motor_2017a_B.Switch_m =
          Asynchronous_motor_2017a_P.Constant1_Value;
      }

      /* End of Switch: '<S21>/Switch' */

      /* Assignment: '<S24>/Lm in rows[2,4] & col[2,4]' incorporates:
       *  Assignment: '<S24>/Lm in rows[1,3] & col[1,3]'
       *  Constant: '<S24>/u1'
       */
      memcpy(&rtb_Lminrows24col24[0], &Asynchronous_motor_2017a_P.u1_Value_j[0],
             sizeof(real_T) << 4U);

      /* Assignment: '<S24>/Lm in rows[1,3] & col[1,3]' incorporates:
       *  Assignment: '<S24>/Lm in rows[2,4] & col[2,4]'
       */
      rtb_Lminrows24col24[0] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[2,4] & col[2,4]' */
      rtb_Lminrows24col24[5] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[1,3] & col[1,3]' incorporates:
       *  Assignment: '<S24>/Lm in rows[2,4] & col[2,4]'
       */
      rtb_Lminrows24col24[2] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[2,4] & col[2,4]' */
      rtb_Lminrows24col24[7] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[1,3] & col[1,3]' incorporates:
       *  Assignment: '<S24>/Lm in rows[2,4] & col[2,4]'
       */
      rtb_Lminrows24col24[8] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[2,4] & col[2,4]' */
      rtb_Lminrows24col24[13] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[1,3] & col[1,3]' incorporates:
       *  Assignment: '<S24>/Lm in rows[2,4] & col[2,4]'
       */
      rtb_Lminrows24col24[10] = Asynchronous_motor_2017a_B.Switch_m;

      /* Assignment: '<S24>/Lm in rows[2,4] & col[2,4]' */
      rtb_Lminrows24col24[15] = Asynchronous_motor_2017a_B.Switch_m;

      /* Sum: '<S24>/Sum2' incorporates:
       *  Assignment: '<S24>/Lm in rows[2,4] & col[2,4]'
       *  Constant: '<S24>/u5'
       */
      for (i = 0; i < 16; i++) {
        rtb_Lminrows24col24_0[i] = rtb_Lminrows24col24[i] +
          Asynchronous_motor_2017a_P.u5_Value[i];
      }

      /* Product: '<S21>/inversion' incorporates:
       *  Sum: '<S24>/Sum2'
       */
      rt_invd4x4_snf(rtb_Lminrows24col24_0, Asynchronous_motor_2017a_B.Linv);
      for (i = 0; i < 4; i++) {
        /* Product: '<S21>/Product1' incorporates:
         *  Constant: '<S21>/u1'
         *  Product: '<S21>/inversion'
         */
        Linv_tmp = i << 2;
        rtb_Switch_idx_1 = Asynchronous_motor_2017a_B.Linv[Linv_tmp + 1];
        rtb_phimd = Asynchronous_motor_2017a_B.Linv[Linv_tmp];
        rtb_Switch_idx_0 = Asynchronous_motor_2017a_B.Linv[Linv_tmp + 2];
        rtb_MultiportSwitch_idx_0 = Asynchronous_motor_2017a_B.Linv[Linv_tmp + 3];
        for (i_0 = 0; i_0 < 4; i_0++) {
          /* Product: '<S21>/Product1' incorporates:
           *  Constant: '<S21>/u1'
           */
          Asynchronous_motor_2017a_B.RLinv[i_0 + Linv_tmp] =
            ((Asynchronous_motor_2017a_P.u1_Value_e[i_0 + 4] * rtb_Switch_idx_1
              + rtb_phimd * Asynchronous_motor_2017a_P.u1_Value_e[i_0]) +
             Asynchronous_motor_2017a_P.u1_Value_e[i_0 + 8] * rtb_Switch_idx_0)
            + Asynchronous_motor_2017a_P.u1_Value_e[i_0 + 12] *
            rtb_MultiportSwitch_idx_0;
        }
      }
    }

    /* End of Constant: '<S15>/Constant' */
    /* End of Outputs for SubSystem: '<S15>/Saturation' */

    /* Switch: '<S15>/Switch' incorporates:
     *  Constant: '<S15>/Constant1'
     *  Constant: '<S15>/Constant2'
     *  Product: '<S21>/inversion'
     */
    rtb_LogicalOperator2 = (Asynchronous_motor_2017a_P.Constant1_Value_m >=
      Asynchronous_motor_2017a_P.Switch_Threshold);
    if (rtb_LogicalOperator2) {
      memcpy(&rtb_inversion_0[0], &Asynchronous_motor_2017a_B.Linv[0], sizeof
             (real_T) << 4U);
    } else {
      memcpy(&rtb_inversion_0[0], &Asynchronous_motor_2017a_P.Constant2_Value[0],
             sizeof(real_T) << 4U);
    }

    /* End of Switch: '<S15>/Switch' */

    /* Product: '<S15>/Product3' */
    rtb_Switch_idx_0 = rtb_xk1[1];
    rtb_Sum3_e = rtb_xk1[0];
    rtb_MultiportSwitch_h_idx_0 = rtb_xk1[2];
    rtb_Delay = rtb_xk1[3];
    for (i = 0; i < 4; i++) {
      rtb_Sum2[i] = ((rtb_inversion_0[i + 4] * rtb_Switch_idx_0 +
                      rtb_inversion_0[i] * rtb_Sum3_e) + rtb_inversion_0[i + 8] *
                     rtb_MultiportSwitch_h_idx_0) + rtb_inversion_0[i + 12] *
        rtb_Delay;
    }

    /* End of Product: '<S15>/Product3' */

    /* UnitDelay: '<S37>/wm_delay' */
    rtb_wm_delay = Asynchronous_motor_2017a_DW.wm_delay_DSTATE;

    /* Sum: '<S37>/Sum1' incorporates:
     *  Gain: '<S37>/F2'
     *  UnitDelay: '<S37>/wm_predict'
     */
    rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.F2_Gain * rtb_wm_delay -
      Asynchronous_motor_2017a_DW.wm_predict_DSTATE;

    /* Outputs for Enabled SubSystem: '<S18>/sin(thr),cos(thr)1' incorporates:
     *  EnablePort: '<S36>/Enable'
     */
    /* Outputs for Enabled SubSystem: '<S18>/sin(thr),cos(thr)' incorporates:
     *  EnablePort: '<S35>/Enable'
     */
    if (rtsiIsModeUpdateTimeStep(&Asynchronous_motor_2017a_M->solverInfo)) {
      /* Constant: '<S18>/Constant' */
      if (Asynchronous_motor_2017a_P.Constant_Value_c) {
        Asynchronous_motor_2017a_DW.sinthrcosthr_MODE = true;
      } else if (Asynchronous_motor_2017a_DW.sinthrcosthr_MODE) {
        /* Disable for Trigonometry: '<S35>/Trigonometric Function' incorporates:
         *  Outport: '<S35>/sin(thr),cos(thr)'
         */
        Asynchronous_motor_2017a_B.TrigonometricFunction_o1_e =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0;

        /* Disable for Trigonometry: '<S35>/Trigonometric Function' incorporates:
         *  Outport: '<S35>/sin(thr),cos(thr)'
         */
        Asynchronous_motor_2017a_B.TrigonometricFunction_o2_l =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0;

        /* Disable for Outport: '<S35>/sin(thr),cos(thr)' incorporates:
         *  Constant: '<S35>/Constant'
         */
        Asynchronous_motor_2017a_B.Constant_f[0] =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0;
        Asynchronous_motor_2017a_B.Constant_f[1] =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0;
        for (i = 0; i < 16; i++) {
          /* Disable for Assignment: '<S35>/W(2,1)=-wr' incorporates:
           *  Outport: '<S35>/W'
           */
          Asynchronous_motor_2017a_B.W21wr[i] =
            Asynchronous_motor_2017a_P.W_Y0_d;
        }

        Asynchronous_motor_2017a_DW.sinthrcosthr_MODE = false;
      }

      /* End of Constant: '<S18>/Constant' */

      /* Constant: '<S18>/Constant1' */
      if (Asynchronous_motor_2017a_P.Constant1_Value_p) {
        Asynchronous_motor_2017a_DW.sinthrcosthr1_MODE = true;
      } else if (Asynchronous_motor_2017a_DW.sinthrcosthr1_MODE) {
        /* Disable for Trigonometry: '<S36>/Trigonometric Function' incorporates:
         *  Outport: '<S36>/sin(thr),cos(thr)'
         */
        Asynchronous_motor_2017a_B.TrigonometricFunction_o1 =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;

        /* Disable for Trigonometry: '<S36>/Trigonometric Function' incorporates:
         *  Outport: '<S36>/sin(thr),cos(thr)'
         */
        Asynchronous_motor_2017a_B.TrigonometricFunction_o2 =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;

        /* Disable for Outport: '<S36>/sin(thr),cos(thr)' incorporates:
         *  Constant: '<S36>/Constant'
         */
        Asynchronous_motor_2017a_B.Constant[0] =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;
        Asynchronous_motor_2017a_B.Constant[1] =
          Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;
        Asynchronous_motor_2017a_DW.sinthrcosthr1_MODE = false;
      }

      /* End of Constant: '<S18>/Constant1' */
    }

    /* End of Outputs for SubSystem: '<S18>/sin(thr),cos(thr)1' */
    if (Asynchronous_motor_2017a_DW.sinthrcosthr_MODE) {
      /* Constant: '<S35>/Constant' */
      Asynchronous_motor_2017a_B.Constant_f[0] =
        Asynchronous_motor_2017a_P.Constant_Value[0];
      Asynchronous_motor_2017a_B.Constant_f[1] =
        Asynchronous_motor_2017a_P.Constant_Value[1];

      /* Trigonometry: '<S35>/Trigonometric Function' incorporates:
       *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
       */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o1_e = sin
        (Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE);

      /* Trigonometry: '<S35>/Trigonometric Function' incorporates:
       *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
       */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o2_l = cos
        (Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE);

      /* Assignment: '<S35>/W(2,1)=-wr' incorporates:
       *  Assignment: '<S35>/W(1,2)=wr'
       *  Constant: '<S35>/u1'
       */
      memcpy(&Asynchronous_motor_2017a_B.W21wr[0],
             &Asynchronous_motor_2017a_P.u1_Value_a[0], sizeof(real_T) << 4U);

      /* Assignment: '<S35>/W(1,2)=wr' incorporates:
       *  Assignment: '<S35>/W(2,1)=-wr'
       */
      Asynchronous_motor_2017a_B.W21wr[4] = rtb_Switch_idx_1;

      /* Assignment: '<S35>/W(2,1)=-wr' incorporates:
       *  Gain: '<S35>/Gain1'
       */
      Asynchronous_motor_2017a_B.W21wr[1] =
        Asynchronous_motor_2017a_P.Gain1_Gain * rtb_Switch_idx_1;
    }

    /* End of Outputs for SubSystem: '<S18>/sin(thr),cos(thr)' */

    /* Outputs for Enabled SubSystem: '<S18>/sin(thr),cos(thr)1' incorporates:
     *  EnablePort: '<S36>/Enable'
     */
    if (Asynchronous_motor_2017a_DW.sinthrcosthr1_MODE) {
      /* Constant: '<S36>/Constant' */
      Asynchronous_motor_2017a_B.Constant[0] =
        Asynchronous_motor_2017a_P.Constant_Value_a[0];
      Asynchronous_motor_2017a_B.Constant[1] =
        Asynchronous_motor_2017a_P.Constant_Value_a[1];

      /* Trigonometry: '<S36>/Trigonometric Function' incorporates:
       *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
       */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o1 = sin
        (Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE);

      /* Trigonometry: '<S36>/Trigonometric Function' incorporates:
       *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
       */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o2 = cos
        (Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE);

      /* Assignment: '<S36>/W(4,3)=wr' incorporates:
       *  Assignment: '<S36>/W(3,4)=-wr'
       *  Constant: '<S36>/u4'
       */
      memcpy(&Asynchronous_motor_2017a_B.W43wr[0],
             &Asynchronous_motor_2017a_P.u4_Value[0], sizeof(real_T) << 4U);

      /* Assignment: '<S36>/W(3,4)=-wr' incorporates:
       *  Assignment: '<S36>/W(4,3)=wr'
       *  Gain: '<S36>/Gain3'
       */
      Asynchronous_motor_2017a_B.W43wr[14] =
        Asynchronous_motor_2017a_P.Gain3_Gain * rtb_Switch_idx_1;

      /* Assignment: '<S36>/W(4,3)=wr' */
      Asynchronous_motor_2017a_B.W43wr[11] = rtb_Switch_idx_1;
    }

    /* End of Outputs for SubSystem: '<S18>/sin(thr),cos(thr)1' */

    /* Outputs for Enabled SubSystem: '<S18>/sin(beta),cos(beta),sin(th),cos(th)' incorporates:
     *  EnablePort: '<S34>/Enable'
     */
    /* Constant: '<S18>/Constant3' */
    if (Asynchronous_motor_2017a_P.Constant3_Value_k) {
      /* Sum: '<S34>/Sum' incorporates:
       *  Constant: '<S34>/we'
       */
      rtb_Sum3_e = Asynchronous_motor_2017a_P.we_Value - rtb_Switch_idx_1;

      /* Gain: '<S34>/web_psb' incorporates:
       *  DigitalClock: '<S34>/Digital Clock'
       */
      rtb_phimd = Asynchronous_motor_2017a_P.web_psb_Gain *
        (((Asynchronous_motor_2017a_M->Timing.clockTick1+
           Asynchronous_motor_2017a_M->Timing.clockTickH1* 4294967296.0)) *
         0.0001);

      /* Sum: '<S34>/Sum1' incorporates:
       *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
       */
      rtb_Delay = rtb_phimd -
        Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE;

      /* Trigonometry: '<S34>/Trigonometric Function' */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o1_f = sin(rtb_phimd);

      /* Trigonometry: '<S34>/Trigonometric Function' */
      Asynchronous_motor_2017a_B.TrigonometricFunction_o2_k = cos(rtb_phimd);

      /* Trigonometry: '<S34>/Trigonometric Function1' */
      Asynchronous_motor_2017a_B.TrigonometricFunction1_o1 = sin(rtb_Delay);

      /* Trigonometry: '<S34>/Trigonometric Function1' */
      Asynchronous_motor_2017a_B.TrigonometricFunction1_o2 = cos(rtb_Delay);

      /* Assignment: '<S34>/W(4,3)=wr-1' incorporates:
       *  Assignment: '<S34>/W(3,4)=1-wr'
       *  Constant: '<S34>/u3'
       */
      memcpy(&Asynchronous_motor_2017a_B.W43wr1[0],
             &Asynchronous_motor_2017a_P.u3_Value[0], sizeof(real_T) << 4U);

      /* Assignment: '<S34>/W(3,4)=1-wr' incorporates:
       *  Assignment: '<S34>/W(4,3)=wr-1'
       */
      Asynchronous_motor_2017a_B.W43wr1[14] = rtb_Sum3_e;

      /* Assignment: '<S34>/W(4,3)=wr-1' incorporates:
       *  Gain: '<S34>/Gain2'
       */
      Asynchronous_motor_2017a_B.W43wr1[11] =
        Asynchronous_motor_2017a_P.Gain2_Gain * rtb_Sum3_e;
    }

    /* End of Constant: '<S18>/Constant3' */
    /* End of Outputs for SubSystem: '<S18>/sin(beta),cos(beta),sin(th),cos(th)' */

    /* MultiPortSwitch: '<S18>/Multiport Switch' incorporates:
     *  Constant: '<S18>/Constant2'
     */
    switch ((int32_T)Asynchronous_motor_2017a_P.Constant2_Value_n) {
     case 1:
      rtb_MultiportSwitch_idx_0 =
        Asynchronous_motor_2017a_B.TrigonometricFunction_o1_e;
      rtb_MultiportSwitch_idx_1 =
        Asynchronous_motor_2017a_B.TrigonometricFunction_o2_l;
      rtb_MultiportSwitch_idx_2 = Asynchronous_motor_2017a_B.Constant_f[0];
      rtb_Sum3_e = Asynchronous_motor_2017a_B.Constant_f[1];
      break;

     case 2:
      rtb_MultiportSwitch_idx_0 =
        Asynchronous_motor_2017a_B.TrigonometricFunction_o1;
      rtb_MultiportSwitch_idx_1 =
        Asynchronous_motor_2017a_B.TrigonometricFunction_o2;
      rtb_MultiportSwitch_idx_2 = Asynchronous_motor_2017a_B.Constant[0];
      rtb_Sum3_e = Asynchronous_motor_2017a_B.Constant[1];
      break;

     default:
      rtb_MultiportSwitch_idx_0 =
        Asynchronous_motor_2017a_B.TrigonometricFunction1_o1;
      rtb_MultiportSwitch_idx_1 =
        Asynchronous_motor_2017a_B.TrigonometricFunction1_o2;
      rtb_MultiportSwitch_idx_2 =
        Asynchronous_motor_2017a_B.TrigonometricFunction_o1_f;
      rtb_Sum3_e = Asynchronous_motor_2017a_B.TrigonometricFunction_o2_k;
      break;
    }

    /* End of MultiPortSwitch: '<S18>/Multiport Switch' */

    /* Outputs for Enabled SubSystem: '<S17>/Stationary reference frame' incorporates:
     *  EnablePort: '<S32>/Enable'
     */
    /* Outputs for Enabled SubSystem: '<S17>/Rotor reference frame' incorporates:
     *  EnablePort: '<S31>/Enable'
     */
    if (rtsiIsModeUpdateTimeStep(&Asynchronous_motor_2017a_M->solverInfo)) {
      /* Constant: '<S17>/Constant' */
      if (Asynchronous_motor_2017a_P.Constant_Value_c3) {
        Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE = true;
      } else if (Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE) {
        /* Disable for Fcn: '<S31>/ira' incorporates:
         *  Outport: '<S31>/ira,irb'
         */
        Asynchronous_motor_2017a_B.ira_e = Asynchronous_motor_2017a_P.irairb_Y0;

        /* Disable for Fcn: '<S31>/irb' incorporates:
         *  Outport: '<S31>/ira,irb'
         */
        Asynchronous_motor_2017a_B.irb_f = Asynchronous_motor_2017a_P.irairb_Y0;

        /* Disable for Fcn: '<S31>/isa' incorporates:
         *  Outport: '<S31>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isa_i = Asynchronous_motor_2017a_P.isaisb_Y0;

        /* Disable for Fcn: '<S31>/isb' incorporates:
         *  Outport: '<S31>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isb_e = Asynchronous_motor_2017a_P.isaisb_Y0;
        Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE = false;
      }

      /* End of Constant: '<S17>/Constant' */

      /* Constant: '<S17>/Constant1' */
      if (Asynchronous_motor_2017a_P.Constant1_Value_j) {
        Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE = true;
      } else if (Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE) {
        /* Disable for Fcn: '<S32>/ira' incorporates:
         *  Outport: '<S32>/ira,irb'
         */
        Asynchronous_motor_2017a_B.ira_a =
          Asynchronous_motor_2017a_P.irairb_Y0_n;

        /* Disable for Fcn: '<S32>/irb' incorporates:
         *  Outport: '<S32>/ira,irb'
         */
        Asynchronous_motor_2017a_B.irb_b =
          Asynchronous_motor_2017a_P.irairb_Y0_n;

        /* Disable for Fcn: '<S32>/isa' incorporates:
         *  Outport: '<S32>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isa_o =
          Asynchronous_motor_2017a_P.isaisb_Y0_p;

        /* Disable for Fcn: '<S32>/isb' incorporates:
         *  Outport: '<S32>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isb_d =
          Asynchronous_motor_2017a_P.isaisb_Y0_p;
        Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE = false;
      }

      /* End of Constant: '<S17>/Constant1' */
    }

    /* End of Outputs for SubSystem: '<S17>/Stationary reference frame' */
    if (Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE) {
      /* Fcn: '<S31>/ira' */
      Asynchronous_motor_2017a_B.ira_e = rtb_Sum2[2];

      /* Fcn: '<S31>/irb' */
      Asynchronous_motor_2017a_B.irb_f = -(1.7320508075688772 * rtb_Sum2[3] +
        rtb_Sum2[2]) / 2.0;

      /* Fcn: '<S31>/isa' */
      Asynchronous_motor_2017a_B.isa_i = rtb_Sum2[0] * rtb_MultiportSwitch_idx_1
        + rtb_MultiportSwitch_idx_0 * rtb_Sum2[1];

      /* Fcn: '<S31>/isb' */
      Asynchronous_motor_2017a_B.isb_e = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_0 - rtb_MultiportSwitch_idx_1) * rtb_Sum2[0] + (
        -1.7320508075688772 * rtb_MultiportSwitch_idx_1 -
        rtb_MultiportSwitch_idx_0) * rtb_Sum2[1]) / 2.0;
    }

    /* End of Outputs for SubSystem: '<S17>/Rotor reference frame' */

    /* Outputs for Enabled SubSystem: '<S17>/Stationary reference frame' incorporates:
     *  EnablePort: '<S32>/Enable'
     */
    if (Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE) {
      /* Fcn: '<S32>/ira' */
      Asynchronous_motor_2017a_B.ira_a = rtb_MultiportSwitch_idx_1 * rtb_Sum2[2]
        - rtb_MultiportSwitch_idx_0 * rtb_Sum2[3];

      /* Fcn: '<S32>/irb' */
      Asynchronous_motor_2017a_B.irb_b = ((-rtb_MultiportSwitch_idx_1 -
        1.7320508075688772 * rtb_MultiportSwitch_idx_0) * rtb_Sum2[2] +
        (rtb_MultiportSwitch_idx_0 - 1.7320508075688772 *
         rtb_MultiportSwitch_idx_1) * rtb_Sum2[3]) / 2.0;

      /* Fcn: '<S32>/isa' */
      Asynchronous_motor_2017a_B.isa_o = rtb_Sum2[0];

      /* Fcn: '<S32>/isb' */
      Asynchronous_motor_2017a_B.isb_d = -(1.7320508075688772 * rtb_Sum2[1] +
        rtb_Sum2[0]) / 2.0;
    }

    /* End of Outputs for SubSystem: '<S17>/Stationary reference frame' */

    /* Outputs for Enabled SubSystem: '<S17>/Synchronous reference frame' incorporates:
     *  EnablePort: '<S33>/Enable'
     */
    if (rtsiIsModeUpdateTimeStep(&Asynchronous_motor_2017a_M->solverInfo)) {
      /* Constant: '<S17>/Constant2' */
      if (Asynchronous_motor_2017a_P.Constant2_Value_a) {
        Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MODE = true;
      } else if (Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MODE) {
        /* Disable for Fcn: '<S33>/ira' incorporates:
         *  Outport: '<S33>/ira,irb'
         */
        Asynchronous_motor_2017a_B.ira = Asynchronous_motor_2017a_P.irairb_Y0_n1;

        /* Disable for Fcn: '<S33>/irb' incorporates:
         *  Outport: '<S33>/ira,irb'
         */
        Asynchronous_motor_2017a_B.irb = Asynchronous_motor_2017a_P.irairb_Y0_n1;

        /* Disable for Fcn: '<S33>/isa' incorporates:
         *  Outport: '<S33>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isa = Asynchronous_motor_2017a_P.isaisb_Y0_g;

        /* Disable for Fcn: '<S33>/isb' incorporates:
         *  Outport: '<S33>/isa,isb'
         */
        Asynchronous_motor_2017a_B.isb = Asynchronous_motor_2017a_P.isaisb_Y0_g;
        Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MODE = false;
      }

      /* End of Constant: '<S17>/Constant2' */
    }

    if (Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MODE) {
      /* Fcn: '<S33>/ira' */
      Asynchronous_motor_2017a_B.ira = rtb_MultiportSwitch_idx_1 * rtb_Sum2[2] +
        rtb_MultiportSwitch_idx_0 * rtb_Sum2[3];

      /* Fcn: '<S33>/irb' */
      Asynchronous_motor_2017a_B.irb = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_0 - rtb_MultiportSwitch_idx_1) * rtb_Sum2[2] + (
        -1.7320508075688772 * rtb_MultiportSwitch_idx_1 -
        rtb_MultiportSwitch_idx_0) * rtb_Sum2[3]) / 2.0;

      /* Fcn: '<S33>/isa' */
      Asynchronous_motor_2017a_B.isa = rtb_Sum2[0] * rtb_Sum3_e + rtb_Sum2[1] *
        rtb_MultiportSwitch_idx_2;

      /* Fcn: '<S33>/isb' */
      Asynchronous_motor_2017a_B.isb = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_2 - rtb_Sum3_e) * rtb_Sum2[0] +
        (-1.7320508075688772 * rtb_Sum3_e - rtb_MultiportSwitch_idx_2) *
        rtb_Sum2[1]) / 2.0;
    }

    /* End of Outputs for SubSystem: '<S17>/Synchronous reference frame' */

    /* MultiPortSwitch: '<S17>/Multiport Switch' incorporates:
     *  Constant: '<S17>/Constant3'
     */
    switch ((int32_T)Asynchronous_motor_2017a_P.Constant3_Value_a) {
     case 1:
      rtb_MultiportSwitch_h_idx_0 = Asynchronous_motor_2017a_B.ira_e;
      rtb_Delay = Asynchronous_motor_2017a_B.irb_f;
      break;

     case 2:
      rtb_MultiportSwitch_h_idx_0 = Asynchronous_motor_2017a_B.ira_a;
      rtb_Delay = Asynchronous_motor_2017a_B.irb_b;
      break;

     default:
      rtb_MultiportSwitch_h_idx_0 = Asynchronous_motor_2017a_B.ira;
      rtb_Delay = Asynchronous_motor_2017a_B.irb;
      break;
    }

    /* End of MultiPortSwitch: '<S17>/Multiport Switch' */

    /* MultiPortSwitch: '<S17>/Multiport Switch1' incorporates:
     *  Constant: '<S17>/Constant4'
     */
    switch ((int32_T)Asynchronous_motor_2017a_P.Constant4_Value_d) {
     case 1:
      rtb_MultiportSwitch1_l_idx_0 = Asynchronous_motor_2017a_B.isa_i;
      rtb_MultiportSwitch1_l_idx_1 = Asynchronous_motor_2017a_B.isb_e;
      break;

     case 2:
      rtb_MultiportSwitch1_l_idx_0 = Asynchronous_motor_2017a_B.isa_o;
      rtb_MultiportSwitch1_l_idx_1 = Asynchronous_motor_2017a_B.isb_d;
      break;

     default:
      rtb_MultiportSwitch1_l_idx_0 = Asynchronous_motor_2017a_B.isa;
      rtb_MultiportSwitch1_l_idx_1 = Asynchronous_motor_2017a_B.isb;
      break;
    }

    /* End of MultiPortSwitch: '<S17>/Multiport Switch1' */

    /* Gain: '<S12>/Gain Vr_Vs1' incorporates:
     *  Gain: '<S17>/ib'
     */
    Asynchronous_motor_2017a_B.GainVr_Vs1[0] =
      Asynchronous_motor_2017a_P.ib_Gain * rtb_MultiportSwitch1_l_idx_0 *
      Asynchronous_motor_2017a_P.GainVr_Vs1_Gain[0];
    Asynchronous_motor_2017a_B.GainVr_Vs1[1] =
      Asynchronous_motor_2017a_P.ib_Gain * rtb_MultiportSwitch1_l_idx_1 *
      Asynchronous_motor_2017a_P.GainVr_Vs1_Gain[1];

    /* S-Function (sfun_spssw_discc): '<S62>/State-Space' incorporates:
     *  Constant: '<S38>/DC'
     */

    /* S-Function block: <S62>/State-Space */
    {
      real_T accum;

      /* Circuit has switches */
      int_T *switch_status = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *switch_status_init = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;
      int_T *SwitchChange = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.SW_CHG;
      int_T *gState = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.IDX_OUT_SW;
      real_T *DxCol = (real_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.DX_COL;
      real_T *tmp2 = (real_T*)Asynchronous_motor_2017a_DW.StateSpace_PWORK.TMP2;
      real_T *uswlast = (real_T*)
        Asynchronous_motor_2017a_DW.StateSpace_PWORK.USWLAST;
      int_T newState;
      int_T swChanged = 0;
      int loopsToDo = 20;
      real_T temp;

      /* keep an initial copy of switch_status*/
      memcpy(switch_status_init, switch_status, 6 * sizeof(int_T));
      memcpy(uswlast, &Asynchronous_motor_2017a_B.StateSpace_o1[0], 6*sizeof
             (real_T));
      do {
        if (loopsToDo == 1) {          /* Need to reset some variables: */
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
          real_T *Ds = (real_T*)Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS;

          {
            int_T i1;
            real_T *y0 = &Asynchronous_motor_2017a_B.StateSpace_o1[0];
            for (i1=0; i1 < 8; i1++) {
              accum = 0.0;

              {
                int_T i2;
                const real_T *u0;
                for (i2=0; i2 < 6; i2++) {
                  accum += *(Ds++) * 0.0;
                }

                accum += *(Ds++) * Asynchronous_motor_2017a_B.GainVr_Vs1[0];
                accum += *(Ds++) * Asynchronous_motor_2017a_B.GainVr_Vs1[1];
                accum += *(Ds++) * Asynchronous_motor_2017a_P.Vdc;
              }

              y0[i1] = accum;
            }
          }

          swChanged = 0;

          {
            int_T i1;
            real_T *y0 = &Asynchronous_motor_2017a_B.StateSpace_o1[0];
            for (i1=0; i1 < 6; i1++) {
              newState = ((y0[i1] > 0.0) && (gState[i1] > 0)) || (y0[i1] < 0.0) ?
                1 : (((y0[i1] > 0.0) && gState[i1] == 0) ? 0 : switch_status[i1]);
              swChanged = ((SwitchChange[i1] = newState - switch_status[i1]) !=
                           0) ? 1 : swChanged;
              switch_status[i1] = newState;/* Keep new state */
            }
          }
        }

        /*
         * Compute new As, Bs, Cs and Ds matrixes:
         * --------------------------------------
         */
        if (swChanged) {
          real_T *Ds = (real_T*)Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS;
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
        }                              /* if (swChanged) */
      } while (swChanged > 0 && --loopsToDo > 0);

      if (loopsToDo == 0) {
        real_T *Ds = (real_T*)Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &Asynchronous_motor_2017a_B.StateSpace_o1[0];
          for (i1=0; i1 < 8; i1++) {
            accum = 0.0;

            {
              int_T i2;
              const real_T *u0;
              for (i2=0; i2 < 6; i2++) {
                accum += *(Ds++) * 0.0;
              }

              accum += *(Ds++) * Asynchronous_motor_2017a_B.GainVr_Vs1[0];
              accum += *(Ds++) * Asynchronous_motor_2017a_B.GainVr_Vs1[1];
              accum += *(Ds++) * Asynchronous_motor_2017a_P.Vdc;
            }

            y0[i1] = accum;
          }
        }
      }

      /* Output new switches states */
      {
        int_T i1;
        real_T *y1 = &Asynchronous_motor_2017a_B.StateSpace_o2[0];
        for (i1=0; i1 < 6; i1++) {
          y1[i1] = (real_T)switch_status[i1];
        }
      }
    }

    /* Sum: '<S19>/Sum2' incorporates:
     *  Gain: '<S19>/1-1'
     *  Product: '<S19>/Mult1'
     */
    rtb_phimd = Asynchronous_motor_2017a_P.u1_Gain[0] * rtb_Sum2[0] * rtb_xk1[1];

    /* Sum: '<S19>/Sum2' incorporates:
     *  Gain: '<S19>/1-1'
     *  Product: '<S19>/Mult1'
     */
    Asynchronous_motor_2017a_B.Sum2 = Asynchronous_motor_2017a_P.u1_Gain[1] *
      rtb_Sum2[1] * rtb_xk1[0] + rtb_phimd;

    /* Gain: '<S14>/1\p' incorporates:
     *  DiscreteIntegrator: '<S14>/Rotor angle thetam'
     */
    rtb_phimd = Asynchronous_motor_2017a_P.up_Gain *
      Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE;

    /* Gain: '<S14>/1\p1' */
    Asynchronous_motor_2017a_B.wTethr[0] = Asynchronous_motor_2017a_P.up1_Gain[0]
      * rtb_Switch_idx_1;
    Asynchronous_motor_2017a_B.wTethr[1] = Asynchronous_motor_2017a_P.up1_Gain[1]
      * Asynchronous_motor_2017a_B.Sum2;
    Asynchronous_motor_2017a_B.wTethr[2] = Asynchronous_motor_2017a_P.up1_Gain[2]
      * rtb_phimd;

    /* Gain: '<Root>/Gain' */
    rtb_phimd = Asynchronous_motor_2017a_P.Gain_Gain_e *
      Asynchronous_motor_2017a_B.wTethr[0];

    /* Sum: '<S17>/Sum2' */
    rtb_phimd = 0.0 - rtb_MultiportSwitch_h_idx_0;

    /* Gain: '<S12>/Gain Vr_Vs' */
    rtb_Switch_idx_0 = Asynchronous_motor_2017a_P.GainVr_Vs_Gain[0] *
      Asynchronous_motor_2017a_B.StateSpace_o1[6];

    /* Gain: '<S16>/1_Vb' incorporates:
     *  Constant: '<S30>/Constant6'
     */
    rtb_u_Vb_idx_0 = Asynchronous_motor_2017a_P.u_Vb_Gain *
      Asynchronous_motor_2017a_P.Constant6_Value[0];
    rtb_u_Vb_idx_2 = Asynchronous_motor_2017a_P.u_Vb_Gain * rtb_Switch_idx_0;

    /* Sum: '<S17>/Sum2' */
    rtb_phimd -= rtb_Delay;

    /* Gain: '<S16>/1_Vb' incorporates:
     *  Constant: '<S30>/Constant6'
     *  Gain: '<S12>/Gain Vr_Vs'
     */
    rtb_u_Vb_idx_1 = Asynchronous_motor_2017a_P.u_Vb_Gain *
      Asynchronous_motor_2017a_P.Constant6_Value[1];
    rtb_u_Vb_idx_3 = Asynchronous_motor_2017a_P.GainVr_Vs_Gain[1] *
      Asynchronous_motor_2017a_B.StateSpace_o1[7] *
      Asynchronous_motor_2017a_P.u_Vb_Gain;

    /* Outputs for Enabled SubSystem: '<S16>/Stationary reference frame' incorporates:
     *  EnablePort: '<S28>/Enable'
     */
    /* Outputs for Enabled SubSystem: '<S16>/Rotor reference frame' incorporates:
     *  EnablePort: '<S27>/Enable'
     */
    if (rtsiIsModeUpdateTimeStep(&Asynchronous_motor_2017a_M->solverInfo)) {
      /* Constant: '<S16>/Constant' */
      if (Asynchronous_motor_2017a_P.Constant_Value_i) {
        Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE_o = true;
      } else if (Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE_o) {
        /* Disable for Fcn: '<S27>/vqr' incorporates:
         *  Outport: '<S27>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vqr_p = Asynchronous_motor_2017a_P.vqrvdr_Y0;

        /* Disable for Fcn: '<S27>/vdr' incorporates:
         *  Outport: '<S27>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vdr_i = Asynchronous_motor_2017a_P.vqrvdr_Y0;

        /* Disable for Fcn: '<S27>/vqs' incorporates:
         *  Outport: '<S27>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vqs_g = Asynchronous_motor_2017a_P.vqsvds_Y0;

        /* Disable for Fcn: '<S27>/vds' incorporates:
         *  Outport: '<S27>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vds_n = Asynchronous_motor_2017a_P.vqsvds_Y0;
        Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE_o = false;
      }

      /* End of Constant: '<S16>/Constant' */

      /* Constant: '<S16>/Constant1' */
      if (Asynchronous_motor_2017a_P.Constant1_Value_mg) {
        Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE_l = true;
      } else if (Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE_l) {
        /* Disable for Fcn: '<S28>/vqr' incorporates:
         *  Outport: '<S28>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vqr_a =
          Asynchronous_motor_2017a_P.vqrvdr_Y0_b;

        /* Disable for Fcn: '<S28>/vdr' incorporates:
         *  Outport: '<S28>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vdr_d =
          Asynchronous_motor_2017a_P.vqrvdr_Y0_b;

        /* Disable for Fcn: '<S28>/vqs' incorporates:
         *  Outport: '<S28>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vqs_h =
          Asynchronous_motor_2017a_P.vqsvds_Y0_j;

        /* Disable for Fcn: '<S28>/vds' incorporates:
         *  Outport: '<S28>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vds_b =
          Asynchronous_motor_2017a_P.vqsvds_Y0_j;
        Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE_l = false;
      }

      /* End of Constant: '<S16>/Constant1' */
    }

    /* End of Outputs for SubSystem: '<S16>/Stationary reference frame' */
    if (Asynchronous_motor_2017a_DW.Rotorreferenceframe_MODE_o) {
      /* Fcn: '<S27>/vdr' */
      Asynchronous_motor_2017a_B.vdr_i = -0.57735026918962573 * rtb_u_Vb_idx_1;

      /* Fcn: '<S27>/vds' */
      Asynchronous_motor_2017a_B.vds_n = ((rtb_MultiportSwitch_idx_0 -
        1.7320508075688772 * rtb_MultiportSwitch_idx_1) * rtb_u_Vb_idx_3 + 2.0 *
        rtb_MultiportSwitch_idx_0 * rtb_u_Vb_idx_2) * 0.33333333333333331;

      /* Fcn: '<S27>/vqr' */
      Asynchronous_motor_2017a_B.vqr_p = (2.0 * rtb_u_Vb_idx_0 + rtb_u_Vb_idx_1)
        * 0.33333333333333331;

      /* Fcn: '<S27>/vqs' */
      Asynchronous_motor_2017a_B.vqs_g = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_0 + rtb_MultiportSwitch_idx_1) * rtb_u_Vb_idx_3
        + 2.0 * rtb_MultiportSwitch_idx_1 * rtb_u_Vb_idx_2) *
        0.33333333333333331;
    }

    /* End of Outputs for SubSystem: '<S16>/Rotor reference frame' */

    /* Outputs for Enabled SubSystem: '<S16>/Stationary reference frame' incorporates:
     *  EnablePort: '<S28>/Enable'
     */
    if (Asynchronous_motor_2017a_DW.Stationaryreferenceframe_MODE_l) {
      /* Fcn: '<S28>/vdr' */
      Asynchronous_motor_2017a_B.vdr_d = ((-rtb_MultiportSwitch_idx_0 -
        1.7320508075688772 * rtb_MultiportSwitch_idx_1) * rtb_u_Vb_idx_1 + -2.0 *
        rtb_MultiportSwitch_idx_0 * rtb_u_Vb_idx_0) * 0.33333333333333331;

      /* Fcn: '<S28>/vds' */
      Asynchronous_motor_2017a_B.vds_b = -0.57735026918962573 * rtb_u_Vb_idx_3;

      /* Fcn: '<S28>/vqr' */
      Asynchronous_motor_2017a_B.vqr_a = ((rtb_MultiportSwitch_idx_1 -
        1.7320508075688772 * rtb_MultiportSwitch_idx_0) * rtb_u_Vb_idx_1 + 2.0 *
        rtb_MultiportSwitch_idx_1 * rtb_u_Vb_idx_0) * 0.33333333333333331;

      /* Fcn: '<S28>/vqs' */
      Asynchronous_motor_2017a_B.vqs_h = (2.0 * rtb_u_Vb_idx_2 + rtb_u_Vb_idx_3)
        * 0.33333333333333331;
    }

    /* End of Outputs for SubSystem: '<S16>/Stationary reference frame' */

    /* Outputs for Enabled SubSystem: '<S16>/Synchronous reference frame' incorporates:
     *  EnablePort: '<S29>/Enable'
     */
    if (rtsiIsModeUpdateTimeStep(&Asynchronous_motor_2017a_M->solverInfo)) {
      /* Constant: '<S16>/Constant2' */
      if (Asynchronous_motor_2017a_P.Constant2_Value_g) {
        Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MOD_a = true;
      } else if (Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MOD_a) {
        /* Disable for Fcn: '<S29>/vqr' incorporates:
         *  Outport: '<S29>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vqr = Asynchronous_motor_2017a_P.vqrvdr_Y0_g;

        /* Disable for Fcn: '<S29>/vdr' incorporates:
         *  Outport: '<S29>/vqr,vdr'
         */
        Asynchronous_motor_2017a_B.vdr = Asynchronous_motor_2017a_P.vqrvdr_Y0_g;

        /* Disable for Fcn: '<S29>/vqs' incorporates:
         *  Outport: '<S29>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vqs = Asynchronous_motor_2017a_P.vqsvds_Y0_b;

        /* Disable for Fcn: '<S29>/vds' incorporates:
         *  Outport: '<S29>/vqs,vds'
         */
        Asynchronous_motor_2017a_B.vds = Asynchronous_motor_2017a_P.vqsvds_Y0_b;
        Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MOD_a = false;
      }

      /* End of Constant: '<S16>/Constant2' */
    }

    if (Asynchronous_motor_2017a_DW.Synchronousreferenceframe_MOD_a) {
      /* Fcn: '<S29>/vdr' */
      Asynchronous_motor_2017a_B.vdr = ((rtb_MultiportSwitch_idx_0 -
        1.7320508075688772 * rtb_MultiportSwitch_idx_1) * rtb_u_Vb_idx_1 + 2.0 *
        rtb_MultiportSwitch_idx_0 * rtb_u_Vb_idx_0) / 3.0;

      /* Fcn: '<S29>/vds' */
      Asynchronous_motor_2017a_B.vds = ((rtb_MultiportSwitch_idx_2 -
        1.7320508075688772 * rtb_Sum3_e) * rtb_u_Vb_idx_3 + 2.0 *
        rtb_MultiportSwitch_idx_2 * rtb_u_Vb_idx_2) / 3.0;

      /* Fcn: '<S29>/vqr' */
      Asynchronous_motor_2017a_B.vqr = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_0 + rtb_MultiportSwitch_idx_1) * rtb_u_Vb_idx_1
        + 2.0 * rtb_MultiportSwitch_idx_1 * rtb_u_Vb_idx_0) / 3.0;

      /* Fcn: '<S29>/vqs' */
      Asynchronous_motor_2017a_B.vqs = ((1.7320508075688772 *
        rtb_MultiportSwitch_idx_2 + rtb_Sum3_e) * rtb_u_Vb_idx_3 + 2.0 *
        rtb_Sum3_e * rtb_u_Vb_idx_2) / 3.0;
    }

    /* End of Outputs for SubSystem: '<S16>/Synchronous reference frame' */

    /* MultiPortSwitch: '<S16>/Multiport Switch' incorporates:
     *  Constant: '<S16>/Constant3'
     */
    switch ((int32_T)Asynchronous_motor_2017a_P.Constant3_Value_m) {
     case 1:
      /* MultiPortSwitch: '<S16>/Multiport Switch' */
      rtb_MultiportSwitch_p[0] = Asynchronous_motor_2017a_B.vqr_p;
      rtb_MultiportSwitch_p[1] = Asynchronous_motor_2017a_B.vdr_i;
      break;

     case 2:
      /* MultiPortSwitch: '<S16>/Multiport Switch' */
      rtb_MultiportSwitch_p[0] = Asynchronous_motor_2017a_B.vqr_a;
      rtb_MultiportSwitch_p[1] = Asynchronous_motor_2017a_B.vdr_d;
      break;

     default:
      /* MultiPortSwitch: '<S16>/Multiport Switch' */
      rtb_MultiportSwitch_p[0] = Asynchronous_motor_2017a_B.vqr;
      rtb_MultiportSwitch_p[1] = Asynchronous_motor_2017a_B.vdr;
      break;
    }

    /* End of MultiPortSwitch: '<S16>/Multiport Switch' */

    /* MultiPortSwitch: '<S16>/Multiport Switch1' incorporates:
     *  Constant: '<S16>/Constant4'
     */
    switch ((int32_T)Asynchronous_motor_2017a_P.Constant4_Value_h) {
     case 1:
      /* MultiPortSwitch: '<S16>/Multiport Switch1' */
      rtb_MultiportSwitch1[0] = Asynchronous_motor_2017a_B.vqs_g;
      rtb_MultiportSwitch1[1] = Asynchronous_motor_2017a_B.vds_n;
      break;

     case 2:
      /* MultiPortSwitch: '<S16>/Multiport Switch1' */
      rtb_MultiportSwitch1[0] = Asynchronous_motor_2017a_B.vqs_h;
      rtb_MultiportSwitch1[1] = Asynchronous_motor_2017a_B.vds_b;
      break;

     default:
      /* MultiPortSwitch: '<S16>/Multiport Switch1' */
      rtb_MultiportSwitch1[0] = Asynchronous_motor_2017a_B.vqs;
      rtb_MultiportSwitch1[1] = Asynchronous_motor_2017a_B.vds;
      break;
    }

    /* End of MultiPortSwitch: '<S16>/Multiport Switch1' */

    /* Gain: '<S12>/unit conversion' incorporates:
     *  Sum: '<S17>/Sum3'
     */
    rtb_unitconversion[2] = Asynchronous_motor_2017a_P.unitconversion_Gain[2] *
      rtb_phimd;
    rtb_unitconversion[11] = ((0.0 - rtb_MultiportSwitch1_l_idx_0) -
      rtb_MultiportSwitch1_l_idx_1) *
      Asynchronous_motor_2017a_P.unitconversion_Gain[11];
    rtb_unitconversion[0] = Asynchronous_motor_2017a_P.unitconversion_Gain[0] *
      rtb_MultiportSwitch_h_idx_0;
    rtb_unitconversion[3] = rtb_Sum2[2] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[3];
    rtb_unitconversion[5] = rtb_xk1[2] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[5];
    rtb_unitconversion[7] = rtb_MultiportSwitch_p[0] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[7];
    rtb_unitconversion[9] = Asynchronous_motor_2017a_P.unitconversion_Gain[9] *
      rtb_MultiportSwitch1_l_idx_0;
    rtb_unitconversion[12] = rtb_Sum2[0] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[12];
    rtb_unitconversion[14] = rtb_xk1[0] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[14];
    rtb_unitconversion[16] = rtb_MultiportSwitch1[0] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[16];
    rtb_unitconversion[1] = Asynchronous_motor_2017a_P.unitconversion_Gain[1] *
      rtb_Delay;
    rtb_unitconversion[4] = rtb_Sum2[3] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[4];
    rtb_unitconversion[6] = rtb_xk1[3] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[6];
    rtb_unitconversion[8] = rtb_MultiportSwitch_p[1] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[8];
    rtb_unitconversion[10] = Asynchronous_motor_2017a_P.unitconversion_Gain[10] *
      rtb_MultiportSwitch1_l_idx_1;
    rtb_unitconversion[13] = rtb_Sum2[1] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[13];
    rtb_unitconversion[15] = rtb_xk1[1] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[15];
    rtb_unitconversion[17] = rtb_MultiportSwitch1[1] *
      Asynchronous_motor_2017a_P.unitconversion_Gain[17];

    /* Switch: '<S15>/Switch2' incorporates:
     *  Constant: '<S15>/Constant5'
     *  Constant: '<S15>/Lm_nosat'
     */
    if (Asynchronous_motor_2017a_P.Constant5_Value >=
        Asynchronous_motor_2017a_P.Switch2_Threshold) {
      rtb_MultiportSwitch_idx_0 = Asynchronous_motor_2017a_B.Switch_m;
    } else {
      rtb_MultiportSwitch_idx_0 = Asynchronous_motor_2017a_P.Lm_nosat_Value;
    }

    /* Gain: '<S12>/unit conversion' incorporates:
     *  Switch: '<S15>/Switch2'
     */
    rtb_unitconversion[18] = Asynchronous_motor_2017a_P.unitconversion_Gain[18] *
      rtb_MultiportSwitch_idx_0;

    /* Gain: '<S52>/Gain3' */
    rtb_phimd = rtb_unitconversion[10];
    rtb_MultiportSwitch_idx_0 = rtb_unitconversion[9];
    rtb_MultiportSwitch_idx_1 = rtb_unitconversion[11];
    for (i = 0; i < 3; i++) {
      /* Gain: '<S52>/Gain1' incorporates:
       *  Gain: '<S52>/Gain3'
       */
      rtb_MathFunction[i] = ((Asynchronous_motor_2017a_P.Gain3_Gain_o[i + 3] *
        rtb_phimd + Asynchronous_motor_2017a_P.Gain3_Gain_o[i] *
        rtb_MultiportSwitch_idx_0) + Asynchronous_motor_2017a_P.Gain3_Gain_o[i +
        6] * rtb_MultiportSwitch_idx_1) *
        Asynchronous_motor_2017a_P.Gain1_Gain_g;
    }

    /* End of Gain: '<S52>/Gain3' */

    /* ComplexToMagnitudeAngle: '<Root>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<Root>/Real-Imag to Complex'
     */
    Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2 = rt_atan2d_snf
      (rtb_unitconversion[5], rtb_unitconversion[6]);

    /* RelationalOperator: '<S53>/Compare' incorporates:
     *  Constant: '<S51>/Constant'
     *  Constant: '<S53>/Constant'
     */
    rtb_Compare = (uint8_T)
      (Asynchronous_motor_2017a_P.AlphaBetaZerotodq0_Alignment ==
       Asynchronous_motor_2017a_P.CompareToConstant_const);

    /* Outputs for Enabled SubSystem: '<S51>/Subsystem1' incorporates:
     *  EnablePort: '<S56>/Enable'
     */
    if (rtb_Compare > 0) {
      /* Fcn: '<S56>/Fcn' incorporates:
       *  Fcn: '<S56>/Fcn1'
       */
      rtb_Sum3_e = sin(Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);
      rtb_phimd = cos(Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);

      /* Fcn: '<S56>/Fcn' */
      Asynchronous_motor_2017a_B.Fcn_l = rtb_MathFunction[0] * rtb_phimd +
        rtb_MathFunction[1] * rtb_Sum3_e;

      /* Fcn: '<S56>/Fcn1' */
      Asynchronous_motor_2017a_B.Fcn1_n = -rtb_MathFunction[0] * rtb_Sum3_e +
        rtb_MathFunction[1] * rtb_phimd;
    }

    /* End of Outputs for SubSystem: '<S51>/Subsystem1' */

    /* Outputs for Enabled SubSystem: '<S51>/Subsystem - pi//2 delay' incorporates:
     *  EnablePort: '<S55>/Enable'
     */
    /* RelationalOperator: '<S54>/Compare' incorporates:
     *  Constant: '<S51>/Constant'
     *  Constant: '<S54>/Constant'
     */
    if (Asynchronous_motor_2017a_P.AlphaBetaZerotodq0_Alignment ==
        Asynchronous_motor_2017a_P.CompareToConstant1_const) {
      /* Fcn: '<S55>/Fcn' incorporates:
       *  Fcn: '<S55>/Fcn1'
       */
      rtb_phimd = cos(Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);
      rtb_MultiportSwitch_idx_0 = sin
        (Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);

      /* Fcn: '<S55>/Fcn' */
      Asynchronous_motor_2017a_B.Fcn_a = rtb_MathFunction[0] *
        rtb_MultiportSwitch_idx_0 - rtb_MathFunction[1] * rtb_phimd;

      /* Fcn: '<S55>/Fcn1' */
      Asynchronous_motor_2017a_B.Fcn1_h = rtb_MathFunction[0] * rtb_phimd +
        rtb_MathFunction[1] * rtb_MultiportSwitch_idx_0;
    }

    /* End of RelationalOperator: '<S54>/Compare' */
    /* End of Outputs for SubSystem: '<S51>/Subsystem - pi//2 delay' */

    /* Switch: '<S51>/Switch' */
    if (rtb_Compare != 0) {
      /* Switch: '<S51>/Switch' */
      Asynchronous_motor_2017a_B.Switch[0] = Asynchronous_motor_2017a_B.Fcn_l;
      Asynchronous_motor_2017a_B.Switch[1] = Asynchronous_motor_2017a_B.Fcn1_n;
    } else {
      /* Switch: '<S51>/Switch' */
      Asynchronous_motor_2017a_B.Switch[0] = Asynchronous_motor_2017a_B.Fcn_a;
      Asynchronous_motor_2017a_B.Switch[1] = Asynchronous_motor_2017a_B.Fcn1_h;
    }

    /* End of Switch: '<S51>/Switch' */

    /* Sum: '<Root>/Sum2' incorporates:
     *  Constant: '<Root>/Constant'
     */
    Asynchronous_motor_2017a_B.Sum2_g = Asynchronous_motor_2017a_P.Ids_rated -
      Asynchronous_motor_2017a_B.Switch[0];

    /* Gain: '<S1>/Gain' */
    Asynchronous_motor_2017a_B.Gain = Asynchronous_motor_2017a_P.Kpc *
      Asynchronous_motor_2017a_B.Sum2_g;

    /* Product: '<S10>/Divide' */
    Asynchronous_motor_2017a_B.Divide = Asynchronous_motor_2017a_B.Switch[1] /
      Asynchronous_motor_2017a_B.Switch[0];

    /* Gain: '<S10>/Gain1' */
    Asynchronous_motor_2017a_B.Gain1 = Asynchronous_motor_2017a_P.np *
      Asynchronous_motor_2017a_B.wTethr[0];
  }

  /* Clock: '<S10>/Clock' incorporates:
   *  Clock: '<S48>/Clock'
   */
  rtb_MultiportSwitch_idx_0 = Asynchronous_motor_2017a_M->Timing.t[0];

  /* Switch: '<S10>/Switch' incorporates:
   *  Clock: '<S10>/Clock'
   */
  if (rtb_MultiportSwitch_idx_0 > Asynchronous_motor_2017a_P.Tsw) {
    rtb_MultiportSwitch_idx_1 = Asynchronous_motor_2017a_B.Divide;
  } else {
    rtb_MultiportSwitch_idx_1 = 0.0;
  }

  /* Sum: '<S10>/Add' incorporates:
   *  Gain: '<S10>/Gain'
   *  Switch: '<S10>/Switch'
   */
  rtb_Sum3_e = 1.0 / Asynchronous_motor_2017a_P.tau_r *
    rtb_MultiportSwitch_idx_1 + Asynchronous_motor_2017a_B.Gain1;
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* RelationalOperator: '<S58>/Compare' incorporates:
     *  Constant: '<S57>/Constant'
     *  Constant: '<S58>/Constant'
     */
    rtb_Compare = (uint8_T)
      (Asynchronous_motor_2017a_P.AlphaBetaZerotodq0_Alignment_b ==
       Asynchronous_motor_2017a_P.CompareToConstant_const_l);

    /* Outputs for Enabled SubSystem: '<S57>/Subsystem1' incorporates:
     *  EnablePort: '<S61>/Enable'
     */
    if (rtb_Compare > 0) {
      /* Fcn: '<S61>/Fcn' incorporates:
       *  Fcn: '<S61>/Fcn1'
       */
      rtb_phimd = sin(Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);

      /* Fcn: '<S61>/Fcn' */
      Asynchronous_motor_2017a_B.Fcn = rtb_unitconversion[6] * cos
        (Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2) +
        rtb_unitconversion[5] * rtb_phimd;
    }

    /* End of Outputs for SubSystem: '<S57>/Subsystem1' */

    /* Outputs for Enabled SubSystem: '<S57>/Subsystem - pi//2 delay' incorporates:
     *  EnablePort: '<S60>/Enable'
     */
    /* RelationalOperator: '<S59>/Compare' incorporates:
     *  Constant: '<S57>/Constant'
     *  Constant: '<S59>/Constant'
     */
    if (Asynchronous_motor_2017a_P.AlphaBetaZerotodq0_Alignment_b ==
        Asynchronous_motor_2017a_P.CompareToConstant1_const_e) {
      /* Fcn: '<S60>/Fcn' */
      Asynchronous_motor_2017a_B.Fcn_b = rtb_unitconversion[6] * sin
        (Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2) -
        rtb_unitconversion[5] * cos
        (Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);
    }

    /* End of RelationalOperator: '<S59>/Compare' */
    /* End of Outputs for SubSystem: '<S57>/Subsystem - pi//2 delay' */

    /* Switch: '<S57>/Switch' */
    if (rtb_Compare != 0) {
      rtb_Switch_idx_0 = Asynchronous_motor_2017a_B.Fcn;
    } else {
      rtb_Switch_idx_0 = Asynchronous_motor_2017a_B.Fcn_b;
    }

    /* End of Switch: '<S57>/Switch' */

    /* Gain: '<S10>/Gain5' incorporates:
     *  Gain: '<S10>/Gain4'
     */
    Asynchronous_motor_2017a_B.Gain5 = Asynchronous_motor_2017a_P.Rr *
      Asynchronous_motor_2017a_P.Lm * rtb_Switch_idx_0 *
      Asynchronous_motor_2017a_P.Gain5_Gain;
  }

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S10>/Gain2'
   *  Gain: '<S10>/Gain3'
   *  Integrator: '<S1>/Integrator'
   *  Product: '<S10>/Product'
   *  Sum: '<S10>/Add2'
   *  Sum: '<S1>/Sum'
   */
  rtb_MultiportSwitch1_l_idx_0 = ((0.0 - rtb_Sum3_e *
    Asynchronous_motor_2017a_B.Switch[1] * Asynchronous_motor_2017a_P.sigma *
    Asynchronous_motor_2017a_P.Ls) - Asynchronous_motor_2017a_B.Gain5) +
    (Asynchronous_motor_2017a_B.Gain +
     Asynchronous_motor_2017a_X.Integrator_CSTATE);

  /* Saturate: '<S1>/Saturation' incorporates:
   *  Saturate: '<S2>/Saturation'
   */
  rtb_MultiportSwitch_idx_2 = -Asynchronous_motor_2017a_P.V_max * 2.0;
  rtb_MultiportSwitch_h_idx_0 = Asynchronous_motor_2017a_P.V_max * 2.0;
  if (rtb_MultiportSwitch1_l_idx_0 > rtb_MultiportSwitch_h_idx_0) {
    rtb_phimd = rtb_MultiportSwitch_h_idx_0;
  } else if (rtb_MultiportSwitch1_l_idx_0 < rtb_MultiportSwitch_idx_2) {
    rtb_phimd = rtb_MultiportSwitch_idx_2;
  } else {
    rtb_phimd = rtb_MultiportSwitch1_l_idx_0;
  }

  /* End of Saturate: '<S1>/Saturation' */

  /* Gain: '<S1>/Gain1' incorporates:
   *  Gain: '<S1>/Gain2'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S1>/Sum3'
   */
  Asynchronous_motor_2017a_B.Gain1_m = (Asynchronous_motor_2017a_B.Sum2_g -
    (rtb_MultiportSwitch1_l_idx_0 - rtb_phimd) * Asynchronous_motor_2017a_P.Kac)
    * Asynchronous_motor_2017a_P.Kic;

  /* Step: '<Root>/Step' incorporates:
   *  Step: '<Root>/Step1'
   */
  rtb_MultiportSwitch_idx_1 = Asynchronous_motor_2017a_M->Timing.t[0];
  if (rtb_MultiportSwitch_idx_1 < Asynchronous_motor_2017a_P.Step_Time) {
    rtb_Delay = Asynchronous_motor_2017a_P.Step_Y0;
  } else {
    rtb_Delay = Asynchronous_motor_2017a_P.Step_YFinal;
  }

  /* Sum: '<Root>/Sum' incorporates:
   *  Gain: '<Root>/Gain1'
   *  Step: '<Root>/Step'
   */
  rtb_MultiportSwitch1_l_idx_0 = Asynchronous_motor_2017a_P.Gain1_Gain_i *
    rtb_Delay - Asynchronous_motor_2017a_B.wTethr[0];

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/Gain'
   *  Integrator: '<S3>/Integrator'
   */
  rtb_MultiportSwitch1_l_idx_1 = Asynchronous_motor_2017a_P.Kps *
    rtb_MultiportSwitch1_l_idx_0 +
    Asynchronous_motor_2017a_X.Integrator_CSTATE_b;

  /* Saturate: '<S3>/Saturation' */
  if (rtb_MultiportSwitch1_l_idx_1 > Asynchronous_motor_2017a_P.Te_rated) {
    rtb_u_Vb_idx_0 = Asynchronous_motor_2017a_P.Te_rated;
  } else if (rtb_MultiportSwitch1_l_idx_1 < -Asynchronous_motor_2017a_P.Te_rated)
  {
    rtb_u_Vb_idx_0 = -Asynchronous_motor_2017a_P.Te_rated;
  } else {
    rtb_u_Vb_idx_0 = rtb_MultiportSwitch1_l_idx_1;
  }

  /* End of Saturate: '<S3>/Saturation' */

  /* Sum: '<Root>/Sum1' incorporates:
   *  Gain: '<S3>/Gain2'
   */
  rtb_Delay = 1.0 / Asynchronous_motor_2017a_P.Kt * rtb_u_Vb_idx_0 -
    Asynchronous_motor_2017a_B.Switch[1];
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* Gain: '<S10>/Gain8' incorporates:
     *  Product: '<S10>/Product2'
     */
    Asynchronous_motor_2017a_B.Gain8 = Asynchronous_motor_2017a_P.Lm /
      Asynchronous_motor_2017a_P.Lr * (rtb_Switch_idx_0 *
      Asynchronous_motor_2017a_B.Gain1);
  }

  /* Sum: '<S2>/Sum1' incorporates:
   *  Gain: '<S10>/Gain6'
   *  Gain: '<S10>/Gain7'
   *  Gain: '<S2>/Gain'
   *  Integrator: '<S2>/Integrator'
   *  Product: '<S10>/Product1'
   *  Sum: '<S10>/Add1'
   *  Sum: '<S2>/Sum'
   */
  rtb_Sum3_e = (rtb_Sum3_e * Asynchronous_motor_2017a_B.Switch[0] *
                Asynchronous_motor_2017a_P.sigma * Asynchronous_motor_2017a_P.Ls
                + Asynchronous_motor_2017a_B.Gain8) +
    (Asynchronous_motor_2017a_P.Kpc * rtb_Delay +
     Asynchronous_motor_2017a_X.Integrator_CSTATE_c);

  /* Saturate: '<S2>/Saturation' */
  if (rtb_Sum3_e > rtb_MultiportSwitch_h_idx_0) {
    rtb_MultiportSwitch_idx_2 = rtb_MultiportSwitch_h_idx_0;
  } else if (!(rtb_Sum3_e < rtb_MultiportSwitch_idx_2)) {
    rtb_MultiportSwitch_idx_2 = rtb_Sum3_e;
  }

  /* Gain: '<S2>/Gain1' incorporates:
   *  Gain: '<S2>/Gain2'
   *  Sum: '<S2>/Sum2'
   *  Sum: '<S2>/Sum3'
   */
  Asynchronous_motor_2017a_B.Gain1_k = (rtb_Delay - (rtb_Sum3_e -
    rtb_MultiportSwitch_idx_2) * Asynchronous_motor_2017a_P.Kac) *
    Asynchronous_motor_2017a_P.Kic;

  /* Gain: '<S3>/Gain1' incorporates:
   *  Gain: '<S3>/Gain3'
   *  Sum: '<S3>/Sum1'
   *  Sum: '<S3>/Sum2'
   */
  Asynchronous_motor_2017a_B.Gain1_f = (rtb_MultiportSwitch1_l_idx_0 -
    (rtb_MultiportSwitch1_l_idx_1 - rtb_u_Vb_idx_0) *
    Asynchronous_motor_2017a_P.Kas) * Asynchronous_motor_2017a_P.Kis;
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* Switch: '<S22>/IC' incorporates:
     *  DigitalClock: '<S22>/Digital Clock'
     */
    if ((((Asynchronous_motor_2017a_M->Timing.clockTick1+
           Asynchronous_motor_2017a_M->Timing.clockTickH1* 4294967296.0)) *
         0.0001) >= Asynchronous_motor_2017a_P.IC_Threshold) {
      /* MultiPortSwitch: '<S18>/Multiport Switch1' incorporates:
       *  Assignment: '<S34>/W(4,3)=wr-1'
       *  Assignment: '<S35>/W(2,1)=-wr'
       *  Assignment: '<S36>/W(4,3)=wr'
       *  Constant: '<S18>/Constant4'
       *  Sum: '<S26>/Sum5'
       */
      switch ((int32_T)Asynchronous_motor_2017a_P.Constant4_Value_n) {
       case 1:
        memcpy(&rtb_Lminrows24col24[0], &Asynchronous_motor_2017a_B.W21wr[0],
               sizeof(real_T) << 4U);
        break;

       case 2:
        memcpy(&rtb_Lminrows24col24[0], &Asynchronous_motor_2017a_B.W43wr[0],
               sizeof(real_T) << 4U);
        break;

       default:
        memcpy(&rtb_Lminrows24col24[0], &Asynchronous_motor_2017a_B.W43wr1[0],
               sizeof(real_T) << 4U);
        break;
      }

      /* End of MultiPortSwitch: '<S18>/Multiport Switch1' */

      /* Switch: '<S15>/Switch1' incorporates:
       *  Constant: '<S15>/Constant3'
       */
      rtb_LogicalOperator2 = (Asynchronous_motor_2017a_P.Constant3_Value >=
        Asynchronous_motor_2017a_P.Switch1_Threshold);
      for (i = 0; i < 16; i++) {
        /* Switch: '<S15>/Switch1' incorporates:
         *  Constant: '<S15>/Constant4'
         *  Product: '<S21>/Product1'
         *  Sum: '<S15>/Sum1'
         *  Sum: '<S26>/Sum5'
         */
        if (rtb_LogicalOperator2) {
          rtb_Delay = Asynchronous_motor_2017a_B.RLinv[i];
        } else {
          rtb_Delay = Asynchronous_motor_2017a_P.Constant4_Value[i];
        }

        /* Gain: '<S26>/wbase*Ts//2' incorporates:
         *  Sum: '<S15>/Sum1'
         *  Sum: '<S26>/Sum5'
         *  Switch: '<S15>/Switch1'
         */
        rtb_Switch_idx_0 = ((0.0 - rtb_Lminrows24col24[i]) - rtb_Delay) *
          Asynchronous_motor_2017a_P.wbaseTs2_Gain;
        rtb_Lminrows24col24[i] = rtb_Switch_idx_0;

        /* Sum: '<S26>/Sum1' incorporates:
         *  Constant: '<S26>/u5'
         *  Sum: '<S26>/Sum5'
         */
        rtb_Sum3_e = Asynchronous_motor_2017a_P.u5_Value_m[i];
        rtb_inversion_0[i] = rtb_Sum3_e - rtb_Switch_idx_0;

        /* Sum: '<S26>/Sum5' incorporates:
         *  Constant: '<S26>/u5'
         */
        rtb_Lminrows24col24_0[i] = rtb_Sum3_e + rtb_Switch_idx_0;
      }

      /* Product: '<S26>/inversion' incorporates:
       *  Sum: '<S26>/Sum1'
       */
      rt_invd4x4_snf(rtb_inversion_0, rtb_Lminrows24col24);

      /* Product: '<S26>/Product4' incorporates:
       *  Product: '<S26>/inversion'
       */
      for (i = 0; i < 4; i++) {
        i_0 = i << 2;
        rtb_Delay = rtb_Lminrows24col24_0[i_0 + 1];
        rtb_Switch_idx_0 = rtb_Lminrows24col24_0[i_0];
        rtb_Sum3_e = rtb_Lminrows24col24_0[i_0 + 2];
        rtb_MultiportSwitch_h_idx_0 = rtb_Lminrows24col24_0[i_0 + 3];
        for (Linv_tmp = 0; Linv_tmp < 4; Linv_tmp++) {
          rtb_inversion_0[Linv_tmp + i_0] = ((rtb_Lminrows24col24[Linv_tmp + 4] *
            rtb_Delay + rtb_Switch_idx_0 * rtb_Lminrows24col24[Linv_tmp]) +
            rtb_Lminrows24col24[Linv_tmp + 8] * rtb_Sum3_e) +
            rtb_Lminrows24col24[Linv_tmp + 12] * rtb_MultiportSwitch_h_idx_0;
        }
      }

      /* End of Product: '<S26>/Product4' */

      /* Sum: '<S22>/sum' incorporates:
       *  UnitDelay: '<S22>/voltages'
       */
      rtb_MultiportSwitch1_l_idx_0 = rtb_MultiportSwitch1[0] +
        Asynchronous_motor_2017a_DW.voltages_DSTATE[0];
      rtb_MultiportSwitch1_l_idx_1 = rtb_MultiportSwitch_p[0] +
        Asynchronous_motor_2017a_DW.voltages_DSTATE[2];
      rtb_u_Vb_idx_0 = rtb_MultiportSwitch1[1] +
        Asynchronous_motor_2017a_DW.voltages_DSTATE[1];
      rtb_u_Vb_idx_2 = rtb_MultiportSwitch_p[1] +
        Asynchronous_motor_2017a_DW.voltages_DSTATE[3];

      /* Product: '<S22>/Product2' */
      rtb_Switch_idx_0 = rtb_xk1[1];
      rtb_Sum3_e = rtb_xk1[0];
      rtb_MultiportSwitch_h_idx_0 = rtb_xk1[2];
      rtb_Delay = rtb_xk1[3];
      for (i = 0; i < 4; i++) {
        /* Switch: '<S22>/IC' incorporates:
         *  Gain: '<S26>/wbase*Ts//2 '
         *  Product: '<S22>/Product1'
         *  Product: '<S26>/inversion'
         *  Sum: '<S22>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
         */
        rtb_IC[i] = (((rtb_Lminrows24col24[i + 4] *
                       Asynchronous_motor_2017a_P.wbaseTs2_Gain_f *
                       rtb_u_Vb_idx_0 +
                       Asynchronous_motor_2017a_P.wbaseTs2_Gain_f *
                       rtb_Lminrows24col24[i] * rtb_MultiportSwitch1_l_idx_0) +
                      rtb_Lminrows24col24[i + 8] *
                      Asynchronous_motor_2017a_P.wbaseTs2_Gain_f *
                      rtb_MultiportSwitch1_l_idx_1) + rtb_Lminrows24col24[i + 12]
                     * Asynchronous_motor_2017a_P.wbaseTs2_Gain_f *
                     rtb_u_Vb_idx_2) + (((rtb_inversion_0[i + 4] *
          rtb_Switch_idx_0 + rtb_inversion_0[i] * rtb_Sum3_e) +
          rtb_inversion_0[i + 8] * rtb_MultiportSwitch_h_idx_0) +
          rtb_inversion_0[i + 12] * rtb_Delay);
      }

      /* End of Product: '<S22>/Product2' */
    } else {
      /* Switch: '<S22>/IC' */
      rtb_IC[0] = rtb_xk1[0];
      rtb_IC[1] = rtb_xk1[1];
      rtb_IC[2] = rtb_xk1[2];
      rtb_IC[3] = rtb_xk1[3];
    }

    /* End of Switch: '<S22>/IC' */

    /* Gain: '<S14>/F' */
    Asynchronous_motor_2017a_B.F = Asynchronous_motor_2017a_P.F_Gain *
      rtb_Switch_idx_1;
  }

  /* Step: '<Root>/Step1' */
  if (rtb_MultiportSwitch_idx_1 < Asynchronous_motor_2017a_P.Step1_Time) {
    rtb_MultiportSwitch_idx_1 = Asynchronous_motor_2017a_P.Step1_Y0;
  } else {
    rtb_MultiportSwitch_idx_1 = Asynchronous_motor_2017a_P.Step1_YFinal;
  }

  /* Gain: '<S14>/1_2H' incorporates:
   *  Gain: '<S14>/Unit conversion'
   *  Step: '<Root>/Step1'
   *  Sum: '<S14>/Sum'
   */
  Asynchronous_motor_2017a_B.u_2H = ((Asynchronous_motor_2017a_B.Sum2 -
    Asynchronous_motor_2017a_P.Unitconversion_Gain * rtb_MultiportSwitch_idx_1)
    - Asynchronous_motor_2017a_B.F) * Asynchronous_motor_2017a_P.u_2H_Gain;
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
    if (Asynchronous_motor_2017a_DW.Rotorspeedwm_SYSTEM_ENABLE != 0) {
      /* DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
      Asynchronous_motor_2017a_B.Rotorspeedwm =
        Asynchronous_motor_2017a_DW.Rotorspeedwm_DSTATE;
    } else {
      /* DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
      Asynchronous_motor_2017a_B.Rotorspeedwm =
        Asynchronous_motor_2017a_P.Rotorspeedwm_gainval *
        Asynchronous_motor_2017a_B.u_2H +
        Asynchronous_motor_2017a_DW.Rotorspeedwm_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S14>/Rotor speed(wm)' */

    /* Gain: '<S14>/web_psb' */
    rtb_web_psb = Asynchronous_motor_2017a_P.web_psb_Gain_l * rtb_Switch_idx_1;
  }

  /* Fcn: '<S6>/alpha' incorporates:
   *  Fcn: '<S6>/beta'
   */
  rtb_Switch_idx_1 = sin(1.5707963267948966 -
    Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);
  rtb_Sum3_e = cos(1.5707963267948966 -
                   Asynchronous_motor_2017a_B.ComplextoMagnitudeAngle_o2);
  rtb_Delay = rtb_Sum3_e * rtb_phimd + rtb_Switch_idx_1 *
    rtb_MultiportSwitch_idx_2;

  /* Fcn: '<S6>/beta' */
  rtb_phimd = -rtb_Switch_idx_1 * rtb_phimd + rtb_Sum3_e *
    rtb_MultiportSwitch_idx_2;

  /* Gain: '<S39>/Gain' */
  rtb_MultiportSwitch_idx_2 = Asynchronous_motor_2017a_P.Gain_Gain_ep *
    rtb_Delay;

  /* Sum: '<S39>/Subtract2' incorporates:
   *  Constant: '<S45>/Constant'
   *  Constant: '<S46>/Constant'
   *  Constant: '<S47>/Constant'
   *  Gain: '<S39>/Gain1'
   *  Gain: '<S39>/Gain2'
   *  Gain: '<S39>/Gain3'
   *  RelationalOperator: '<S45>/Compare'
   *  RelationalOperator: '<S46>/Compare'
   *  RelationalOperator: '<S47>/Compare'
   *  Sum: '<S39>/Subtract'
   *  Sum: '<S39>/Subtract1'
   */
  rtb_Compare = (uint8_T)(((((uint32_T)(rtb_MultiportSwitch_idx_2 - rtb_phimd >
    Asynchronous_motor_2017a_P.B_const ? (int32_T)
    Asynchronous_motor_2017a_P.Gain2_Gain_j : 0) << 1) + (uint32_T)(rtb_phimd >
    Asynchronous_motor_2017a_P.A_const ? (int32_T)
    Asynchronous_motor_2017a_P.Gain1_Gain_g0 : 0)) + ((uint32_T)((0.0 -
    rtb_MultiportSwitch_idx_2) - rtb_phimd > Asynchronous_motor_2017a_P.C_const ?
    (int32_T)Asynchronous_motor_2017a_P.Gain3_Gain_c : 0) << 2)) >> 2);

  /* Gain: '<S44>/Gain' */
  rtb_MultiportSwitch_idx_2 = Asynchronous_motor_2017a_P.pxyz * rtb_phimd;

  /* Gain: '<S44>/Gain2' */
  rtb_phimd *= Asynchronous_motor_2017a_P.Gain2_Gain_m;

  /* Gain: '<S44>/Gain1' */
  rtb_Delay *= Asynchronous_motor_2017a_P.Gain1_Gain_p;

  /* Gain: '<S44>/Gain3' incorporates:
   *  Sum: '<S44>/Subtract'
   */
  rtb_MultiportSwitch1_l_idx_1 = (rtb_phimd + rtb_Delay) *
    Asynchronous_motor_2017a_P.pxyz;

  /* Gain: '<S44>/Gain4' incorporates:
   *  Sum: '<S44>/Subtract1'
   */
  rtb_Delay = (rtb_phimd - rtb_Delay) * Asynchronous_motor_2017a_P.pxyz;

  /* MultiPortSwitch: '<S40>/Multiport Switch' incorporates:
   *  Gain: '<S40>/Gain'
   *  Gain: '<S40>/Gain1'
   *  Gain: '<S40>/Gain2'
   *  Sum: '<S39>/Subtract2'
   */
  switch ((int32_T)((uint32_T)rtb_Compare >> 5)) {
   case 1:
    rtb_Sum3_e = rtb_Delay;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' */
    rtb_Switch_idx_1 = rtb_MultiportSwitch1_l_idx_1;
    break;

   case 2:
    rtb_Sum3_e = rtb_MultiportSwitch1_l_idx_1;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' incorporates:
     *  Gain: '<S40>/Gain'
     */
    rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.Gain_Gain_o *
      rtb_MultiportSwitch_idx_2;
    break;

   case 3:
    rtb_Sum3_e = Asynchronous_motor_2017a_P.Gain2_Gain_o * rtb_Delay;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' incorporates:
     *  Gain: '<S40>/Gain2'
     */
    rtb_Switch_idx_1 = rtb_MultiportSwitch_idx_2;
    break;

   case 4:
    rtb_Sum3_e = Asynchronous_motor_2017a_P.Gain_Gain_o *
      rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' incorporates:
     *  Gain: '<S40>/Gain'
     */
    rtb_Switch_idx_1 = rtb_Delay;
    break;

   case 5:
    rtb_Sum3_e = rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' incorporates:
     *  Gain: '<S40>/Gain1'
     */
    rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.Gain1_Gain_n *
      rtb_MultiportSwitch1_l_idx_1;
    break;

   default:
    rtb_Sum3_e = Asynchronous_motor_2017a_P.Gain1_Gain_n *
      rtb_MultiportSwitch1_l_idx_1;

    /* MultiPortSwitch: '<S40>/Multiport Switch1' incorporates:
     *  Gain: '<S40>/Gain1'
     *  Gain: '<S40>/Gain2'
     */
    rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.Gain2_Gain_o * rtb_Delay;
    break;
  }

  /* End of MultiPortSwitch: '<S40>/Multiport Switch' */

  /* Gain: '<S42>/Gain' incorporates:
   *  Constant: '<S42>/Constant'
   *  Sum: '<S42>/Subtract'
   */
  rtb_MultiportSwitch_idx_2 = ((Asynchronous_motor_2017a_P.Tsw - rtb_Sum3_e) -
    rtb_Switch_idx_1) * Asynchronous_motor_2017a_P.Gain_Gain_n;

  /* Sum: '<S42>/Subtract1' incorporates:
   *  Gain: '<S42>/Gain1'
   */
  rtb_Delay = Asynchronous_motor_2017a_P.Gain1_Gain_ii * rtb_Sum3_e +
    rtb_MultiportSwitch_idx_2;

  /* Sum: '<S42>/Subtract2' incorporates:
   *  Gain: '<S42>/Gain2'
   */
  rtb_Switch_idx_1 = Asynchronous_motor_2017a_P.Gain2_Gain_d * rtb_Switch_idx_1
    + rtb_Delay;

  /* MultiPortSwitch: '<S41>/Multiport Switch1' incorporates:
   *  Sum: '<S39>/Subtract2'
   */
  switch ((int32_T)((uint32_T)rtb_Compare >> 5)) {
   case 1:
    rtb_u_Vb_idx_0 = rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_Switch_idx_1;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_Delay;
    break;

   case 2:
    rtb_u_Vb_idx_0 = rtb_Switch_idx_1;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_Delay;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_MultiportSwitch_idx_2;
    break;

   case 3:
    rtb_u_Vb_idx_0 = rtb_Delay;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_Switch_idx_1;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_MultiportSwitch_idx_2;
    break;

   case 4:
    rtb_u_Vb_idx_0 = rtb_Delay;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_Switch_idx_1;
    break;

   case 5:
    rtb_u_Vb_idx_0 = rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_Delay;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_Switch_idx_1;
    break;

   default:
    rtb_u_Vb_idx_0 = rtb_Switch_idx_1;

    /* MultiPortSwitch: '<S41>/Multiport Switch2' */
    rtb_MultiportSwitch1_l_idx_1 = rtb_MultiportSwitch_idx_2;

    /* MultiPortSwitch: '<S41>/Multiport Switch3' */
    rtb_phimd = rtb_Delay;
    break;
  }

  /* End of MultiPortSwitch: '<S41>/Multiport Switch1' */

  /* Lookup_n-D: '<S48>/Look-Up Table1' incorporates:
   *  Constant: '<S48>/Constant'
   *  Math: '<S48>/Math Function'
   */
  rtb_Switch_idx_1 = look1_binlxpw(rt_remd_snf(rtb_MultiportSwitch_idx_0,
    Asynchronous_motor_2017a_P.Constant_Value_o),
    Asynchronous_motor_2017a_P.LookUpTable1_bp01Data,
    Asynchronous_motor_2017a_P.triangle_rep_seq_y, 2U);

  /* RelationalOperator: '<S43>/Relational Operator' */
  rtb_LogicalOperator2 = (rtb_Switch_idx_1 > rtb_phimd);

  /* DataTypeConversion: '<S43>/Cast To Double' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[0] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S43>/Cast To Double1' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   *  Logic: '<S43>/Logical Operator'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[1] = !rtb_LogicalOperator2;

  /* RelationalOperator: '<S43>/Relational Operator1' */
  rtb_LogicalOperator2 = (rtb_Switch_idx_1 > rtb_u_Vb_idx_0);

  /* DataTypeConversion: '<S43>/Cast To Double2' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[2] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S43>/Cast To Double3' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   *  Logic: '<S43>/Logical Operator1'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[3] = !rtb_LogicalOperator2;

  /* RelationalOperator: '<S43>/Relational Operator2' */
  rtb_LogicalOperator2 = (rtb_Switch_idx_1 > rtb_MultiportSwitch1_l_idx_1);

  /* DataTypeConversion: '<S43>/Cast To Double4' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[4] = rtb_LogicalOperator2;

  /* DataTypeConversion: '<S43>/Cast To Double5' incorporates:
   *  Concatenate: '<S43>/Vector Concatenate'
   *  Logic: '<S43>/Logical Operator2'
   */
  Asynchronous_motor_2017a_B.VectorConcatenate[5] = !rtb_LogicalOperator2;
  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(Asynchronous_motor_2017a_M->rtwLogInfo,
                        (Asynchronous_motor_2017a_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
      /* Update for UnitDelay: '<S22>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE[0] = rtb_IC[0];

      /* Update for UnitDelay: '<S20>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[0] = rtb_xk1[0];

      /* Update for UnitDelay: '<S22>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE[1] = rtb_IC[1];

      /* Update for UnitDelay: '<S20>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[1] = rtb_xk1[1];

      /* Update for UnitDelay: '<S22>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE[2] = rtb_IC[2];

      /* Update for UnitDelay: '<S20>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[2] = rtb_xk1[2];

      /* Update for UnitDelay: '<S22>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE[3] = rtb_IC[3];

      /* Update for UnitDelay: '<S20>/fluxes' */
      Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[3] = rtb_xk1[3];

      /* Update for DiscreteIntegrator: '<S14>/Rotor angle thetam' */
      Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE +=
        Asynchronous_motor_2017a_P.Rotoranglethetam_gainval * rtb_web_psb;

      /* Update for UnitDelay: '<S37>/wm_delay' */
      Asynchronous_motor_2017a_DW.wm_delay_DSTATE =
        Asynchronous_motor_2017a_B.Rotorspeedwm;

      /* Update for UnitDelay: '<S37>/wm_predict' */
      Asynchronous_motor_2017a_DW.wm_predict_DSTATE = rtb_wm_delay;

      /* Update for S-Function (sfun_spssw_discc): '<S62>/State-Space' incorporates:
       *  Constant: '<S38>/DC'
       */
      {
        int_T *gState = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.G_STATE;

        /* Store switch gates values for next step */
        {
          int_T i1;
          const real_T *u1 = &Asynchronous_motor_2017a_B.VectorConcatenate[0];
          for (i1=0; i1 < 6; i1++) {
            *(gState++) = (int_T) u1[i1];
          }
        }
      }

      /* Update for UnitDelay: '<S22>/voltages' */
      Asynchronous_motor_2017a_DW.voltages_DSTATE[0] = rtb_MultiportSwitch1[0];
      Asynchronous_motor_2017a_DW.voltages_DSTATE[2] = rtb_MultiportSwitch_p[0];
      Asynchronous_motor_2017a_DW.voltages_DSTATE[1] = rtb_MultiportSwitch1[1];
      Asynchronous_motor_2017a_DW.voltages_DSTATE[3] = rtb_MultiportSwitch_p[1];

      /* Update for DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
      Asynchronous_motor_2017a_DW.Rotorspeedwm_SYSTEM_ENABLE = 0U;
      Asynchronous_motor_2017a_DW.Rotorspeedwm_DSTATE =
        Asynchronous_motor_2017a_P.Rotorspeedwm_gainval *
        Asynchronous_motor_2017a_B.u_2H +
        Asynchronous_motor_2017a_B.Rotorspeedwm;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Asynchronous_motor_2017a_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(Asynchronous_motor_2017a_M)!=-1) &&
          !((rtmGetTFinal(Asynchronous_motor_2017a_M)-
             (((Asynchronous_motor_2017a_M->Timing.clockTick1+
                Asynchronous_motor_2017a_M->Timing.clockTickH1* 4294967296.0)) *
              0.0001)) > (((Asynchronous_motor_2017a_M->Timing.clockTick1+
                            Asynchronous_motor_2017a_M->Timing.clockTickH1*
                            4294967296.0)) * 0.0001) * (DBL_EPSILON))) {
        rtmSetErrorStatus(Asynchronous_motor_2017a_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&Asynchronous_motor_2017a_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Asynchronous_motor_2017a_M->Timing.clockTick0)) {
      ++Asynchronous_motor_2017a_M->Timing.clockTickH0;
    }

    Asynchronous_motor_2017a_M->Timing.t[0] = rtsiGetSolverStopTime
      (&Asynchronous_motor_2017a_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.0001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.0001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      Asynchronous_motor_2017a_M->Timing.clockTick1++;
      if (!Asynchronous_motor_2017a_M->Timing.clockTick1) {
        Asynchronous_motor_2017a_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Asynchronous_motor_2017a_derivatives(void)
{
  XDot_Asynchronous_motor_2017a_T *_rtXdot;
  _rtXdot = ((XDot_Asynchronous_motor_2017a_T *)
             Asynchronous_motor_2017a_M->derivs);

  /* Derivatives for Integrator: '<S1>/Integrator' */
  _rtXdot->Integrator_CSTATE = Asynchronous_motor_2017a_B.Gain1_m;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  _rtXdot->Integrator_CSTATE_b = Asynchronous_motor_2017a_B.Gain1_f;

  /* Derivatives for Integrator: '<S2>/Integrator' */
  _rtXdot->Integrator_CSTATE_c = Asynchronous_motor_2017a_B.Gain1_k;
}

/* Model initialize function */
void Asynchronous_motor_2017a_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Asynchronous_motor_2017a_M, 0,
                sizeof(RT_MODEL_Asynchronous_motor_2_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Asynchronous_motor_2017a_M->solverInfo,
                          &Asynchronous_motor_2017a_M->Timing.simTimeStep);
    rtsiSetTPtr(&Asynchronous_motor_2017a_M->solverInfo, &rtmGetTPtr
                (Asynchronous_motor_2017a_M));
    rtsiSetStepSizePtr(&Asynchronous_motor_2017a_M->solverInfo,
                       &Asynchronous_motor_2017a_M->Timing.stepSize0);
    rtsiSetdXPtr(&Asynchronous_motor_2017a_M->solverInfo,
                 &Asynchronous_motor_2017a_M->derivs);
    rtsiSetContStatesPtr(&Asynchronous_motor_2017a_M->solverInfo, (real_T **)
                         &Asynchronous_motor_2017a_M->contStates);
    rtsiSetNumContStatesPtr(&Asynchronous_motor_2017a_M->solverInfo,
      &Asynchronous_motor_2017a_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Asynchronous_motor_2017a_M->solverInfo,
      &Asynchronous_motor_2017a_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Asynchronous_motor_2017a_M->solverInfo,
      &Asynchronous_motor_2017a_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Asynchronous_motor_2017a_M->solverInfo,
      &Asynchronous_motor_2017a_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Asynchronous_motor_2017a_M->solverInfo,
                          (&rtmGetErrorStatus(Asynchronous_motor_2017a_M)));
    rtsiSetRTModelPtr(&Asynchronous_motor_2017a_M->solverInfo,
                      Asynchronous_motor_2017a_M);
  }

  rtsiSetSimTimeStep(&Asynchronous_motor_2017a_M->solverInfo, MAJOR_TIME_STEP);
  Asynchronous_motor_2017a_M->intgData.y = Asynchronous_motor_2017a_M->odeY;
  Asynchronous_motor_2017a_M->intgData.f[0] = Asynchronous_motor_2017a_M->odeF[0];
  Asynchronous_motor_2017a_M->intgData.f[1] = Asynchronous_motor_2017a_M->odeF[1];
  Asynchronous_motor_2017a_M->intgData.f[2] = Asynchronous_motor_2017a_M->odeF[2];
  Asynchronous_motor_2017a_M->contStates = ((X_Asynchronous_motor_2017a_T *)
    &Asynchronous_motor_2017a_X);
  rtsiSetSolverData(&Asynchronous_motor_2017a_M->solverInfo, (void *)
                    &Asynchronous_motor_2017a_M->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&Asynchronous_motor_2017a_M->solverInfo,
    false);
  rtsiSetSolverName(&Asynchronous_motor_2017a_M->solverInfo,"ode3");
  rtmSetTPtr(Asynchronous_motor_2017a_M,
             &Asynchronous_motor_2017a_M->Timing.tArray[0]);
  rtmSetTFinal(Asynchronous_motor_2017a_M, 1.0);
  Asynchronous_motor_2017a_M->Timing.stepSize0 = 0.0001;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    Asynchronous_motor_2017a_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(Asynchronous_motor_2017a_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(Asynchronous_motor_2017a_M->rtwLogInfo, (NULL));
    rtliSetLogT(Asynchronous_motor_2017a_M->rtwLogInfo, "tout");
    rtliSetLogX(Asynchronous_motor_2017a_M->rtwLogInfo, "");
    rtliSetLogXFinal(Asynchronous_motor_2017a_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(Asynchronous_motor_2017a_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(Asynchronous_motor_2017a_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(Asynchronous_motor_2017a_M->rtwLogInfo, 0);
    rtliSetLogDecimation(Asynchronous_motor_2017a_M->rtwLogInfo, 1);
    rtliSetLogY(Asynchronous_motor_2017a_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(Asynchronous_motor_2017a_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(Asynchronous_motor_2017a_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &Asynchronous_motor_2017a_B), 0,
                sizeof(B_Asynchronous_motor_2017a_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Asynchronous_motor_2017a_X, 0,
                  sizeof(X_Asynchronous_motor_2017a_T));
  }

  /* states (dwork) */
  (void) memset((void *)&Asynchronous_motor_2017a_DW, 0,
                sizeof(DW_Asynchronous_motor_2017a_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(Asynchronous_motor_2017a_M->rtwLogInfo, 0.0,
    rtmGetTFinal(Asynchronous_motor_2017a_M),
    Asynchronous_motor_2017a_M->Timing.stepSize0, (&rtmGetErrorStatus
    (Asynchronous_motor_2017a_M)));

  /* Start for S-Function (sfun_spssw_discc): '<S62>/State-Space' incorporates:
   *  Constant: '<S38>/DC'
   */

  /* S-Function block: <S62>/State-Space */
  {
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS = (real_T*)calloc(8 * 9,
      sizeof(real_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(8,
      sizeof(real_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(9,
      sizeof(real_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc
      (6, sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(6,
      sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.G_STATE = (int_T*)calloc(6,
      sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(6,
      sizeof(real_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(6,
      sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(6,
      sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)
      calloc(6, sizeof(int_T));
    Asynchronous_motor_2017a_DW.StateSpace_PWORK.USWLAST = (real_T*)calloc(6,
      sizeof(real_T));
  }

  {
    int32_T i;

    /* InitializeConditions for UnitDelay: '<S22>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE[0] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition[0];

    /* InitializeConditions for UnitDelay: '<S20>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[0] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition_p[0];

    /* InitializeConditions for UnitDelay: '<S22>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE[1] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition[1];

    /* InitializeConditions for UnitDelay: '<S20>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[1] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition_p[1];

    /* InitializeConditions for UnitDelay: '<S22>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE[2] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition[2];

    /* InitializeConditions for UnitDelay: '<S20>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[2] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition_p[2];

    /* InitializeConditions for UnitDelay: '<S22>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE[3] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition[3];

    /* InitializeConditions for UnitDelay: '<S20>/fluxes' */
    Asynchronous_motor_2017a_DW.fluxes_DSTATE_d[3] =
      Asynchronous_motor_2017a_P.fluxes_InitialCondition_p[3];

    /* InitializeConditions for DiscreteIntegrator: '<S14>/Rotor angle thetam' */
    Asynchronous_motor_2017a_DW.Rotoranglethetam_DSTATE =
      Asynchronous_motor_2017a_P.Rotoranglethetam_IC;

    /* InitializeConditions for UnitDelay: '<S37>/wm_delay' */
    Asynchronous_motor_2017a_DW.wm_delay_DSTATE =
      Asynchronous_motor_2017a_P.wm_delay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S37>/wm_predict' */
    Asynchronous_motor_2017a_DW.wm_predict_DSTATE =
      Asynchronous_motor_2017a_P.wm_predict_InitialCondition;

    /* InitializeConditions for S-Function (sfun_spssw_discc): '<S62>/State-Space' incorporates:
     *  Constant: '<S38>/DC'
     */
    {
      int32_T i, j;
      real_T *Ds = (real_T*)Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS;

      /* Copy and transpose D */
      for (i=0; i<8; i++) {
        for (j=0; j<9; j++)
          Ds[i*9 + j] = (Asynchronous_motor_2017a_P.StateSpace_DS_param[i + j*8]);
      }

      {
        /* Switches work vectors */
        int_T *switch_status = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS;
        int_T *gState = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.G_STATE;
        real_T *yswitch = (real_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.Y_SWITCH;
        int_T *switchTypes = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_TYPES;
        int_T *idxOutSw = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.IDX_OUT_SW;
        int_T *switch_status_init = (int_T*)
          Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;

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

    /* InitializeConditions for Integrator: '<S1>/Integrator' */
    Asynchronous_motor_2017a_X.Integrator_CSTATE =
      Asynchronous_motor_2017a_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S3>/Integrator' */
    Asynchronous_motor_2017a_X.Integrator_CSTATE_b =
      Asynchronous_motor_2017a_P.Integrator_IC_f;

    /* InitializeConditions for Integrator: '<S2>/Integrator' */
    Asynchronous_motor_2017a_X.Integrator_CSTATE_c =
      Asynchronous_motor_2017a_P.Integrator_IC_b;

    /* InitializeConditions for UnitDelay: '<S22>/voltages' */
    Asynchronous_motor_2017a_DW.voltages_DSTATE[0] =
      Asynchronous_motor_2017a_P.voltages_InitialCondition;
    Asynchronous_motor_2017a_DW.voltages_DSTATE[1] =
      Asynchronous_motor_2017a_P.voltages_InitialCondition;
    Asynchronous_motor_2017a_DW.voltages_DSTATE[2] =
      Asynchronous_motor_2017a_P.voltages_InitialCondition;
    Asynchronous_motor_2017a_DW.voltages_DSTATE[3] =
      Asynchronous_motor_2017a_P.voltages_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
    Asynchronous_motor_2017a_DW.Rotorspeedwm_DSTATE =
      Asynchronous_motor_2017a_P.Rotorspeedwm_IC;

    /* SystemInitialize for Enabled SubSystem: '<S15>/Saturation' */
    /* SystemInitialize for Switch: '<S21>/Switch' incorporates:
     *  Outport: '<S21>/Lm'
     */
    Asynchronous_motor_2017a_B.Switch_m = Asynchronous_motor_2017a_P.Lm_Y0;

    /* End of SystemInitialize for SubSystem: '<S15>/Saturation' */

    /* SystemInitialize for Enabled SubSystem: '<S18>/sin(thr),cos(thr)' */
    /* SystemInitialize for Trigonometry: '<S35>/Trigonometric Function' incorporates:
     *  Outport: '<S35>/sin(thr),cos(thr)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o1_e =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0;

    /* SystemInitialize for Trigonometry: '<S35>/Trigonometric Function' incorporates:
     *  Outport: '<S35>/sin(thr),cos(thr)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o2_l =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0;

    /* SystemInitialize for Outport: '<S35>/sin(thr),cos(thr)' incorporates:
     *  Constant: '<S35>/Constant'
     */
    Asynchronous_motor_2017a_B.Constant_f[0] =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0;
    Asynchronous_motor_2017a_B.Constant_f[1] =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0;

    /* End of SystemInitialize for SubSystem: '<S18>/sin(thr),cos(thr)' */

    /* SystemInitialize for Enabled SubSystem: '<S18>/sin(thr),cos(thr)1' */
    /* SystemInitialize for Trigonometry: '<S36>/Trigonometric Function' incorporates:
     *  Outport: '<S36>/sin(thr),cos(thr)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o1 =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;

    /* SystemInitialize for Trigonometry: '<S36>/Trigonometric Function' incorporates:
     *  Outport: '<S36>/sin(thr),cos(thr)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o2 =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;

    /* SystemInitialize for Outport: '<S36>/sin(thr),cos(thr)' incorporates:
     *  Constant: '<S36>/Constant'
     */
    Asynchronous_motor_2017a_B.Constant[0] =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;
    Asynchronous_motor_2017a_B.Constant[1] =
      Asynchronous_motor_2017a_P.sinthrcosthr_Y0_n;

    /* End of SystemInitialize for SubSystem: '<S18>/sin(thr),cos(thr)1' */

    /* SystemInitialize for Enabled SubSystem: '<S18>/sin(beta),cos(beta),sin(th),cos(th)' */
    /* SystemInitialize for Trigonometry: '<S34>/Trigonometric Function1' incorporates:
     *  Outport: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction1_o1 =
      Asynchronous_motor_2017a_P.sinbetacosbetasinthcosth_Y0;

    /* SystemInitialize for Trigonometry: '<S34>/Trigonometric Function1' incorporates:
     *  Outport: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction1_o2 =
      Asynchronous_motor_2017a_P.sinbetacosbetasinthcosth_Y0;

    /* SystemInitialize for Trigonometry: '<S34>/Trigonometric Function' incorporates:
     *  Outport: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o1_f =
      Asynchronous_motor_2017a_P.sinbetacosbetasinthcosth_Y0;

    /* SystemInitialize for Trigonometry: '<S34>/Trigonometric Function' incorporates:
     *  Outport: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
     */
    Asynchronous_motor_2017a_B.TrigonometricFunction_o2_k =
      Asynchronous_motor_2017a_P.sinbetacosbetasinthcosth_Y0;

    /* SystemInitialize for Enabled SubSystem: '<S18>/sin(thr),cos(thr)1' */
    /* SystemInitialize for Enabled SubSystem: '<S18>/sin(thr),cos(thr)' */
    /* SystemInitialize for Enabled SubSystem: '<S15>/Saturation' */
    for (i = 0; i < 16; i++) {
      /* SystemInitialize for Product: '<S21>/inversion' incorporates:
       *  Outport: '<S21>/Linv'
       */
      Asynchronous_motor_2017a_B.Linv[i] = Asynchronous_motor_2017a_P.Linv_Y0;

      /* SystemInitialize for Product: '<S21>/Product1' incorporates:
       *  Outport: '<S21>/R*Linv'
       */
      Asynchronous_motor_2017a_B.RLinv[i] = Asynchronous_motor_2017a_P.RLinv_Y0;

      /* SystemInitialize for Assignment: '<S35>/W(2,1)=-wr' incorporates:
       *  Outport: '<S35>/W'
       */
      Asynchronous_motor_2017a_B.W21wr[i] = Asynchronous_motor_2017a_P.W_Y0_d;

      /* SystemInitialize for Assignment: '<S36>/W(4,3)=wr' incorporates:
       *  Outport: '<S36>/W'
       */
      Asynchronous_motor_2017a_B.W43wr[i] = Asynchronous_motor_2017a_P.W_Y0_o;

      /* SystemInitialize for Assignment: '<S34>/W(4,3)=wr-1' incorporates:
       *  Outport: '<S34>/W'
       */
      Asynchronous_motor_2017a_B.W43wr1[i] = Asynchronous_motor_2017a_P.W_Y0;
    }

    /* End of SystemInitialize for SubSystem: '<S15>/Saturation' */
    /* End of SystemInitialize for SubSystem: '<S18>/sin(thr),cos(thr)' */
    /* End of SystemInitialize for SubSystem: '<S18>/sin(thr),cos(thr)1' */
    /* End of SystemInitialize for SubSystem: '<S18>/sin(beta),cos(beta),sin(th),cos(th)' */

    /* SystemInitialize for Enabled SubSystem: '<S17>/Rotor reference frame' */
    /* SystemInitialize for Fcn: '<S31>/ira' incorporates:
     *  Outport: '<S31>/ira,irb'
     */
    Asynchronous_motor_2017a_B.ira_e = Asynchronous_motor_2017a_P.irairb_Y0;

    /* SystemInitialize for Fcn: '<S31>/irb' incorporates:
     *  Outport: '<S31>/ira,irb'
     */
    Asynchronous_motor_2017a_B.irb_f = Asynchronous_motor_2017a_P.irairb_Y0;

    /* SystemInitialize for Fcn: '<S31>/isa' incorporates:
     *  Outport: '<S31>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isa_i = Asynchronous_motor_2017a_P.isaisb_Y0;

    /* SystemInitialize for Fcn: '<S31>/isb' incorporates:
     *  Outport: '<S31>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isb_e = Asynchronous_motor_2017a_P.isaisb_Y0;

    /* End of SystemInitialize for SubSystem: '<S17>/Rotor reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S17>/Stationary reference frame' */
    /* SystemInitialize for Fcn: '<S32>/ira' incorporates:
     *  Outport: '<S32>/ira,irb'
     */
    Asynchronous_motor_2017a_B.ira_a = Asynchronous_motor_2017a_P.irairb_Y0_n;

    /* SystemInitialize for Fcn: '<S32>/irb' incorporates:
     *  Outport: '<S32>/ira,irb'
     */
    Asynchronous_motor_2017a_B.irb_b = Asynchronous_motor_2017a_P.irairb_Y0_n;

    /* SystemInitialize for Fcn: '<S32>/isa' incorporates:
     *  Outport: '<S32>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isa_o = Asynchronous_motor_2017a_P.isaisb_Y0_p;

    /* SystemInitialize for Fcn: '<S32>/isb' incorporates:
     *  Outport: '<S32>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isb_d = Asynchronous_motor_2017a_P.isaisb_Y0_p;

    /* End of SystemInitialize for SubSystem: '<S17>/Stationary reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S17>/Synchronous reference frame' */
    /* SystemInitialize for Fcn: '<S33>/ira' incorporates:
     *  Outport: '<S33>/ira,irb'
     */
    Asynchronous_motor_2017a_B.ira = Asynchronous_motor_2017a_P.irairb_Y0_n1;

    /* SystemInitialize for Fcn: '<S33>/irb' incorporates:
     *  Outport: '<S33>/ira,irb'
     */
    Asynchronous_motor_2017a_B.irb = Asynchronous_motor_2017a_P.irairb_Y0_n1;

    /* SystemInitialize for Fcn: '<S33>/isa' incorporates:
     *  Outport: '<S33>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isa = Asynchronous_motor_2017a_P.isaisb_Y0_g;

    /* SystemInitialize for Fcn: '<S33>/isb' incorporates:
     *  Outport: '<S33>/isa,isb'
     */
    Asynchronous_motor_2017a_B.isb = Asynchronous_motor_2017a_P.isaisb_Y0_g;

    /* End of SystemInitialize for SubSystem: '<S17>/Synchronous reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S16>/Rotor reference frame' */
    /* SystemInitialize for Fcn: '<S27>/vqr' incorporates:
     *  Outport: '<S27>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vqr_p = Asynchronous_motor_2017a_P.vqrvdr_Y0;

    /* SystemInitialize for Fcn: '<S27>/vdr' incorporates:
     *  Outport: '<S27>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vdr_i = Asynchronous_motor_2017a_P.vqrvdr_Y0;

    /* SystemInitialize for Fcn: '<S27>/vqs' incorporates:
     *  Outport: '<S27>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vqs_g = Asynchronous_motor_2017a_P.vqsvds_Y0;

    /* SystemInitialize for Fcn: '<S27>/vds' incorporates:
     *  Outport: '<S27>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vds_n = Asynchronous_motor_2017a_P.vqsvds_Y0;

    /* End of SystemInitialize for SubSystem: '<S16>/Rotor reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S16>/Stationary reference frame' */
    /* SystemInitialize for Fcn: '<S28>/vqr' incorporates:
     *  Outport: '<S28>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vqr_a = Asynchronous_motor_2017a_P.vqrvdr_Y0_b;

    /* SystemInitialize for Fcn: '<S28>/vdr' incorporates:
     *  Outport: '<S28>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vdr_d = Asynchronous_motor_2017a_P.vqrvdr_Y0_b;

    /* SystemInitialize for Fcn: '<S28>/vqs' incorporates:
     *  Outport: '<S28>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vqs_h = Asynchronous_motor_2017a_P.vqsvds_Y0_j;

    /* SystemInitialize for Fcn: '<S28>/vds' incorporates:
     *  Outport: '<S28>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vds_b = Asynchronous_motor_2017a_P.vqsvds_Y0_j;

    /* End of SystemInitialize for SubSystem: '<S16>/Stationary reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S16>/Synchronous reference frame' */
    /* SystemInitialize for Fcn: '<S29>/vqr' incorporates:
     *  Outport: '<S29>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vqr = Asynchronous_motor_2017a_P.vqrvdr_Y0_g;

    /* SystemInitialize for Fcn: '<S29>/vdr' incorporates:
     *  Outport: '<S29>/vqr,vdr'
     */
    Asynchronous_motor_2017a_B.vdr = Asynchronous_motor_2017a_P.vqrvdr_Y0_g;

    /* SystemInitialize for Fcn: '<S29>/vqs' incorporates:
     *  Outport: '<S29>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vqs = Asynchronous_motor_2017a_P.vqsvds_Y0_b;

    /* SystemInitialize for Fcn: '<S29>/vds' incorporates:
     *  Outport: '<S29>/vqs,vds'
     */
    Asynchronous_motor_2017a_B.vds = Asynchronous_motor_2017a_P.vqsvds_Y0_b;

    /* End of SystemInitialize for SubSystem: '<S16>/Synchronous reference frame' */

    /* SystemInitialize for Enabled SubSystem: '<S51>/Subsystem1' */
    /* SystemInitialize for Fcn: '<S56>/Fcn' incorporates:
     *  Outport: '<S56>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn_l = Asynchronous_motor_2017a_P.dq_Y0_k[0];

    /* SystemInitialize for Fcn: '<S56>/Fcn1' incorporates:
     *  Outport: '<S56>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn1_n = Asynchronous_motor_2017a_P.dq_Y0_k[1];

    /* End of SystemInitialize for SubSystem: '<S51>/Subsystem1' */

    /* SystemInitialize for Enabled SubSystem: '<S51>/Subsystem - pi//2 delay' */
    /* SystemInitialize for Fcn: '<S55>/Fcn' incorporates:
     *  Outport: '<S55>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn_a = Asynchronous_motor_2017a_P.dq_Y0[0];

    /* SystemInitialize for Fcn: '<S55>/Fcn1' incorporates:
     *  Outport: '<S55>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn1_h = Asynchronous_motor_2017a_P.dq_Y0[1];

    /* End of SystemInitialize for SubSystem: '<S51>/Subsystem - pi//2 delay' */

    /* SystemInitialize for Enabled SubSystem: '<S57>/Subsystem1' */
    /* SystemInitialize for Fcn: '<S61>/Fcn' incorporates:
     *  Outport: '<S61>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn = Asynchronous_motor_2017a_P.dq_Y0_m[0];

    /* End of SystemInitialize for SubSystem: '<S57>/Subsystem1' */

    /* SystemInitialize for Enabled SubSystem: '<S57>/Subsystem - pi//2 delay' */
    /* SystemInitialize for Fcn: '<S60>/Fcn' incorporates:
     *  Outport: '<S60>/dq'
     */
    Asynchronous_motor_2017a_B.Fcn_b = Asynchronous_motor_2017a_P.dq_Y0_h[0];

    /* End of SystemInitialize for SubSystem: '<S57>/Subsystem - pi//2 delay' */
  }

  /* Enable for DiscreteIntegrator: '<S14>/Rotor speed(wm)' */
  Asynchronous_motor_2017a_DW.Rotorspeedwm_SYSTEM_ENABLE = 1U;
}

/* Model terminate function */
void Asynchronous_motor_2017a_terminate(void)
{
  /* Terminate for S-Function (sfun_spssw_discc): '<S62>/State-Space' incorporates:
   *  Constant: '<S38>/DC'
   */

  /* S-Function block: <S62>/State-Space */
  {
    /* Free memory */
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.DS);
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.DX_COL);
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.TMP2);

    /*
     * Circuit has switches*/
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.G_STATE);
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS);
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.SW_CHG);
    free(Asynchronous_motor_2017a_DW.StateSpace_PWORK.SWITCH_STATUS_INIT);
  }
}
