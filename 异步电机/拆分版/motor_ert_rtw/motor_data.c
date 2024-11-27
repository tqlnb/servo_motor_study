/*
 * File: motor_data.c
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

/* Constant parameters (default storage) */
const ConstP rtConstP = {
  /* Expression: SM.Linv
   * Referenced by: '<S8>/Constant2'
   */
  { 11.778689645591433, 0.0, -11.431106508699781, -0.0, 0.0, 11.778689645591433,
    0.0, -11.431106508699781, -11.431106508699781, 0.0, 11.778689645591433, -0.0,
    -0.0, -11.431106508699781, -0.0, 11.778689645591433 },

  /* Expression: SM.RLinv
   * Referenced by: '<S8>/Constant4'
   */
  { 0.26562892041659192, 0.0, -0.33119505489978723, 0.0, 0.0,
    0.26562892041659192, 0.0, -0.33119505489978723, -0.25779034616210356, 0.0,
    0.34126562995894349, 0.0, 0.0, -0.25779034616210356, 0.0,
    0.34126562995894349 },

  /* Expression: eye(4,4)
   * Referenced by: '<S19>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  /* Expression: S.D
   * Referenced by: '<S34>/State-Space'
   */
  { -50000.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, -50000.0, 0.0,
    0.0, 0.0, 0.0, -50000.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, 0.0,
    -50000.0, 50000.0, 0.0, 0.0, 50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0,
    0.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    50000.0, -50000.0, 0.0, 50000.0, 50000.0, -50000.0, 0.0, 0.0, -50000.0,
    50000.0, -50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0, -50000.0, 50000.0,
    50000.0, -100000.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0 }
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
