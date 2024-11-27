/*
 * File: control_data.c
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

/* Constant parameters (default storage) */
const ConstP rtConstP = {
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S20>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: rep_seq_t
   * Referenced by: '<S18>/Look-Up Table'
   */
  { 0.0, 5.0E-5, 0.0001 },

  /* Expression: rep_seq_y
   * Referenced by: '<S18>/Look-Up Table'
   */
  { 0.0, 5.0E-5, 0.0 }
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
