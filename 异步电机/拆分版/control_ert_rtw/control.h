/*
 * File: control.h
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

#ifndef RTW_HEADER_control_h_
#define RTW_HEADER_control_h_
#ifndef control_COMMON_INCLUDES_
#define control_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* control_COMMON_INCLUDES_ */

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S1>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_g;/* '<S3>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_p;/* '<S2>/Discrete-Time Integrator' */
} DW;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S20>/Gain3'
   */
  real_T Gain3_Gain[9];

  /* Expression: rep_seq_t
   * Referenced by: '<S18>/Look-Up Table'
   */
  real_T LookUpTable_XData[3];

  /* Expression: rep_seq_y
   * Referenced by: '<S18>/Look-Up Table'
   */
  real_T LookUpTable_YData[3];
} ConstP;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T ia;                           /* '<Root>/ia' */
  real_T ib;                           /* '<Root>/ib' */
  real_T ic;                           /* '<Root>/ic' */
  real_T vd;                           /* '<Root>/vd' */
  real_T vq;                           /* '<Root>/vq' */
  real_T wm;                           /* '<Root>/wm' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T gate[6];                      /* '<Root>/gate' */
} ExtY;

/* Real-time Model Data Structure */
struct tag_RTM {
  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
  } Timing;
};

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Constant parameters (default storage) */
extern const ConstP rtConstP;

/* Model entry point functions */
extern void control_initialize(void);
extern void control_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S4>/0' : Unused code path elimination
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
 * '<Root>' : 'control'
 * '<S1>'   : 'control/ACR_d'
 * '<S2>'   : 'control/ACR_q'
 * '<S3>'   : 'control/ASR'
 * '<S4>'   : 'control/Park to Clarke Angle Transform'
 * '<S5>'   : 'control/SVPWM'
 * '<S6>'   : 'control/abc to dq0'
 * '<S7>'   : 'control/feedforward'
 * '<S8>'   : 'control/powergui'
 * '<S9>'   : 'control/SVPWM/Sector'
 * '<S10>'  : 'control/SVPWM/T1T2'
 * '<S11>'  : 'control/SVPWM/allocate'
 * '<S12>'  : 'control/SVPWM/fuse'
 * '<S13>'  : 'control/SVPWM/gating'
 * '<S14>'  : 'control/SVPWM/xyz'
 * '<S15>'  : 'control/SVPWM/Sector/准则A'
 * '<S16>'  : 'control/SVPWM/Sector/准则B'
 * '<S17>'  : 'control/SVPWM/Sector/准则C'
 * '<S18>'  : 'control/SVPWM/gating/triangle'
 * '<S19>'  : 'control/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S20>'  : 'control/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S21>'  : 'control/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S22>'  : 'control/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S23>'  : 'control/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S24>'  : 'control/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S25>'  : 'control/feedforward/Alpha-Beta-Zero to dq0'
 * '<S26>'  : 'control/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S27>'  : 'control/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S28>'  : 'control/feedforward/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S29>'  : 'control/feedforward/Alpha-Beta-Zero to dq0/Subsystem1'
 */
#endif                                 /* RTW_HEADER_control_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
