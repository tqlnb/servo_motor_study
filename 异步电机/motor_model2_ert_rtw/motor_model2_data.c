/*
 * File: motor_model2_data.c
 *
 * Code generated for Simulink model 'motor_model2'.
 *
 * Model version                  : 13.19
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Mon Oct 28 15:28:25 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex-A
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "motor_model2.h"

/* Constant parameters (default storage) */
const ConstP rtConstP = {
  /* Expression: SM.Linv
   * Referenced by: '<S16>/Constant2'
   */
  { 11.778689645591433, 0.0, -11.431106508699781, -0.0, 0.0, 11.778689645591433,
    0.0, -11.431106508699781, -11.431106508699781, 0.0, 11.778689645591433, -0.0,
    -0.0, -11.431106508699781, -0.0, 11.778689645591433 },

  /* Expression: SM.RLinv
   * Referenced by: '<S16>/Constant4'
   */
  { 0.26562892041659192, 0.0, -0.33119505489978723, 0.0, 0.0,
    0.26562892041659192, 0.0, -0.33119505489978723, -0.25779034616210356, 0.0,
    0.34126562995894349, 0.0, 0.0, -0.25779034616210356, 0.0,
    0.34126562995894349 },

  /* Expression: eye(4,4)
   * Referenced by: '<S27>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  /* Expression: S.D
   * Referenced by: '<S63>/State-Space'
   */
  { -50000.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, -50000.0, 0.0,
    0.0, 0.0, 0.0, -50000.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, 0.0,
    -50000.0, 50000.0, 0.0, 0.0, 50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0,
    0.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    50000.0, -50000.0, 0.0, 50000.0, 50000.0, -50000.0, 0.0, 0.0, -50000.0,
    50000.0, -50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0, -50000.0, 50000.0,
    50000.0, -100000.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0 },

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S53>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: rep_seq_y
   * Referenced by: '<S49>/Look-Up Table1'
   */
  { 0.0, 0.0005, 0.0 },

  /* Expression: rep_seq_t - min(rep_seq_t)
   * Referenced by: '<S49>/Look-Up Table1'
   */
  { 0.0, 0.0005, 0.001 }
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
