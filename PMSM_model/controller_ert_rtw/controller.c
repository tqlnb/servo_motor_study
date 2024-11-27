/*
 * File: controller.c
 *
 * Code generated for Simulink model 'controller'.
 *
 * Model version                  : 1.1
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Mon Nov 18 18:11:07 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex-M
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "controller.h"
#include "rtwtypes.h"

/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Model step function */
void controller_step(void)
{
  real_T rtb_FilterCoefficient;
  real_T rtb_Sum;

  /* Sum: '<Root>/Sum' incorporates:
   *  DiscreteIntegrator: '<Root>/Discrete-Time Integrator'
   *  Inport: '<Root>/In1'
   */
  rtb_Sum = rtU.In1 - rtDW.DiscreteTimeIntegrator_DSTATE;

  /* Gain: '<S36>/Filter Coefficient' incorporates:
   *  DiscreteIntegrator: '<S28>/Filter'
   *  Gain: '<S27>/Derivative Gain'
   *  Sum: '<S28>/SumD'
   */
  rtb_FilterCoefficient = (0.0 * rtb_Sum - rtDW.Filter_DSTATE) * 100.0;

  /* Outport: '<Root>/Out1' incorporates:
   *  DiscreteIntegrator: '<S33>/Integrator'
   *  Gain: '<S38>/Proportional Gain'
   *  Sum: '<S42>/Sum'
   */
  rtY.Out1 = (100.0 * rtb_Sum + rtDW.Integrator_DSTATE) + rtb_FilterCoefficient;

  /* Update for DiscreteIntegrator: '<Root>/Discrete-Time Integrator' incorporates:
   *  Inport: '<Root>/In2'
   */
  rtDW.DiscreteTimeIntegrator_DSTATE += 0.2 * rtU.In2;

  /* Update for DiscreteIntegrator: '<S33>/Integrator' */
  rtDW.Integrator_DSTATE += 0.2 * rtb_Sum;

  /* Update for DiscreteIntegrator: '<S28>/Filter' */
  rtDW.Filter_DSTATE += 0.2 * rtb_FilterCoefficient;
}

/* Model initialize function */
void controller_initialize(void)
{
  /* (no initialization code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
