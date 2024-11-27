/*
 * File: controller.h
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

#ifndef RTW_HEADER_controller_h_
#define RTW_HEADER_controller_h_
#ifndef controller_COMMON_INCLUDES_
#define controller_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* controller_COMMON_INCLUDES_ */

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteTimeIntegrator_DSTATE;/* '<Root>/Discrete-Time Integrator' */
  real_T Integrator_DSTATE;            /* '<S33>/Integrator' */
  real_T Filter_DSTATE;                /* '<S28>/Filter' */
} DW;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T In1;                          /* '<Root>/In1' */
  real_T In2;                          /* '<Root>/In2' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Out1;                         /* '<Root>/Out1' */
} ExtY;

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Model entry point functions */
extern void controller_initialize(void);
extern void controller_step(void);

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S30>/Integral Gain' : Eliminated nontunable gain of 1
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'controller'
 * '<S1>'   : 'controller/Discrete PID Controller'
 * '<S2>'   : 'controller/Discrete PID Controller/Anti-windup'
 * '<S3>'   : 'controller/Discrete PID Controller/D Gain'
 * '<S4>'   : 'controller/Discrete PID Controller/Filter'
 * '<S5>'   : 'controller/Discrete PID Controller/Filter ICs'
 * '<S6>'   : 'controller/Discrete PID Controller/I Gain'
 * '<S7>'   : 'controller/Discrete PID Controller/Ideal P Gain'
 * '<S8>'   : 'controller/Discrete PID Controller/Ideal P Gain Fdbk'
 * '<S9>'   : 'controller/Discrete PID Controller/Integrator'
 * '<S10>'  : 'controller/Discrete PID Controller/Integrator ICs'
 * '<S11>'  : 'controller/Discrete PID Controller/N Copy'
 * '<S12>'  : 'controller/Discrete PID Controller/N Gain'
 * '<S13>'  : 'controller/Discrete PID Controller/P Copy'
 * '<S14>'  : 'controller/Discrete PID Controller/Parallel P Gain'
 * '<S15>'  : 'controller/Discrete PID Controller/Reset Signal'
 * '<S16>'  : 'controller/Discrete PID Controller/Saturation'
 * '<S17>'  : 'controller/Discrete PID Controller/Saturation Fdbk'
 * '<S18>'  : 'controller/Discrete PID Controller/Sum'
 * '<S19>'  : 'controller/Discrete PID Controller/Sum Fdbk'
 * '<S20>'  : 'controller/Discrete PID Controller/Tracking Mode'
 * '<S21>'  : 'controller/Discrete PID Controller/Tracking Mode Sum'
 * '<S22>'  : 'controller/Discrete PID Controller/Tsamp - Integral'
 * '<S23>'  : 'controller/Discrete PID Controller/Tsamp - Ngain'
 * '<S24>'  : 'controller/Discrete PID Controller/postSat Signal'
 * '<S25>'  : 'controller/Discrete PID Controller/preSat Signal'
 * '<S26>'  : 'controller/Discrete PID Controller/Anti-windup/Passthrough'
 * '<S27>'  : 'controller/Discrete PID Controller/D Gain/Internal Parameters'
 * '<S28>'  : 'controller/Discrete PID Controller/Filter/Disc. Forward Euler Filter'
 * '<S29>'  : 'controller/Discrete PID Controller/Filter ICs/Internal IC - Filter'
 * '<S30>'  : 'controller/Discrete PID Controller/I Gain/Internal Parameters'
 * '<S31>'  : 'controller/Discrete PID Controller/Ideal P Gain/Passthrough'
 * '<S32>'  : 'controller/Discrete PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S33>'  : 'controller/Discrete PID Controller/Integrator/Discrete'
 * '<S34>'  : 'controller/Discrete PID Controller/Integrator ICs/Internal IC'
 * '<S35>'  : 'controller/Discrete PID Controller/N Copy/Disabled'
 * '<S36>'  : 'controller/Discrete PID Controller/N Gain/Internal Parameters'
 * '<S37>'  : 'controller/Discrete PID Controller/P Copy/Disabled'
 * '<S38>'  : 'controller/Discrete PID Controller/Parallel P Gain/Internal Parameters'
 * '<S39>'  : 'controller/Discrete PID Controller/Reset Signal/Disabled'
 * '<S40>'  : 'controller/Discrete PID Controller/Saturation/Passthrough'
 * '<S41>'  : 'controller/Discrete PID Controller/Saturation Fdbk/Disabled'
 * '<S42>'  : 'controller/Discrete PID Controller/Sum/Sum_PID'
 * '<S43>'  : 'controller/Discrete PID Controller/Sum Fdbk/Disabled'
 * '<S44>'  : 'controller/Discrete PID Controller/Tracking Mode/Disabled'
 * '<S45>'  : 'controller/Discrete PID Controller/Tracking Mode Sum/Passthrough'
 * '<S46>'  : 'controller/Discrete PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S47>'  : 'controller/Discrete PID Controller/Tsamp - Ngain/Passthrough'
 * '<S48>'  : 'controller/Discrete PID Controller/postSat Signal/Forward_Path'
 * '<S49>'  : 'controller/Discrete PID Controller/preSat Signal/Forward_Path'
 */
#endif                                 /* RTW_HEADER_controller_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
