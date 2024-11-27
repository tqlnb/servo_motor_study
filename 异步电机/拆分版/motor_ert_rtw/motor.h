/*
 * File: motor.h
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

#ifndef RTW_HEADER_motor_h_
#define RTW_HEADER_motor_h_
#ifndef motor_COMMON_INCLUDES_
#define motor_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* motor_COMMON_INCLUDES_ */

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real_T ib[4];                        /* '<S10>/ib' */
  real_T StateSpace_o1[8];             /* '<S34>/State-Space' */
  real_T StateSpace_o2[6];             /* '<S34>/State-Space' */
  real_T W21wr[16];                    /* '<S28>/W(2,1)=-wr' */
  real_T Constant_e[2];                /* '<S28>/Constant' */
  real_T fluxes_DSTATE[4];             /* '<S15>/fluxes' */
  real_T fluxes_DSTATE_l[4];           /* '<S13>/fluxes' */
  real_T voltages_DSTATE[4];           /* '<S15>/voltages' */
  real_T inversion_DWORK4[16];         /* '<S19>/inversion' */
  real_T inversion_DWORK4_k[16];       /* '<S14>/inversion' */
  real_T TrigonometricFunction_o1_d;   /* '<S28>/Trigonometric Function' */
  real_T TrigonometricFunction_o2_j;   /* '<S28>/Trigonometric Function' */
  real_T ira;                          /* '<S26>/ira' */
  real_T irb;                          /* '<S26>/irb' */
  real_T isa;                          /* '<S26>/isa' */
  real_T isb;                          /* '<S26>/isb' */
  real_T ira_p;                        /* '<S24>/ira' */
  real_T irb_n;                        /* '<S24>/irb' */
  real_T isa_m;                        /* '<S24>/isa' */
  real_T isb_p;                        /* '<S24>/isb' */
  real_T vdr;                          /* '<S22>/vdr' */
  real_T vds;                          /* '<S22>/vds' */
  real_T vqr;                          /* '<S22>/vqr' */
  real_T vqs;                          /* '<S22>/vqs' */
  real_T vdr_l;                        /* '<S20>/vdr' */
  real_T vds_f;                        /* '<S20>/vds' */
  real_T vqr_j;                        /* '<S20>/vqr' */
  real_T vqs_l;                        /* '<S20>/vqs' */
  real_T Rotoranglethetam_DSTATE;      /* '<S7>/Rotor angle thetam' */
  real_T wm_delay_DSTATE;              /* '<S30>/wm_delay' */
  real_T wm_predict_DSTATE;            /* '<S30>/wm_predict' */
  real_T Rotorspeedwm_DSTATE;          /* '<S7>/Rotor speed(wm)' */
  struct {
    void *AS;
    void *BS;
    void *CS;
    void *DS;
    void *DX_COL;
    void *BD_COL;
    void *TMP1;
    void *TMP2;
    void *XTMP;
    void *SWITCH_STATUS;
    void *SWITCH_STATUS_INIT;
    void *SW_CHG;
    void *G_STATE;
    void *USWLAST;
    void *XKM12;
    void *XKP12;
    void *XLAST;
    void *ULAST;
    void *IDX_SW_CHG;
    void *Y_SWITCH;
    void *SWITCH_TYPES;
    void *IDX_OUT_SW;
    void *SWITCH_TOPO_SAVED_IDX;
    void *SWITCH_MAP;
  } StateSpace_PWORK;                  /* '<S34>/State-Space' */

  int_T StateSpace_IWORK[11];          /* '<S34>/State-Space' */
  uint8_T Rotorspeedwm_SYSTEM_ENABLE;  /* '<S7>/Rotor speed(wm)' */
  boolean_T sinthrcosthr_MODE;         /* '<S11>/sin(thr),cos(thr)' */
  boolean_T Synchronousreferenceframe_MODE;/* '<S10>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE;  /* '<S10>/Rotor reference frame' */
  boolean_T Synchronousreferenceframe_MOD_e;/* '<S9>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE_p;/* '<S9>/Rotor reference frame' */
} DW;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: SM.Linv
   * Referenced by: '<S8>/Constant2'
   */
  real_T Constant2_Value[16];

  /* Expression: SM.RLinv
   * Referenced by: '<S8>/Constant4'
   */
  real_T Constant4_Value[16];

  /* Expression: eye(4,4)
   * Referenced by: '<S19>/u5'
   */
  real_T u5_Value_l[16];

  /* Expression: S.D
   * Referenced by: '<S34>/State-Space'
   */
  real_T StateSpace_DS_param[72];
} ConstP;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T Tm;                           /* '<Root>/Tm' */
  real_T gate[6];                      /* '<Root>/gate' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T ia;                           /* '<Root>/ia' */
  real_T ib;                           /* '<Root>/ib' */
  real_T ic;                           /* '<Root>/ic' */
  real_T vd;                           /* '<Root>/vd' */
  real_T vq;                           /* '<Root>/vq' */
  real_T wm;                           /* '<Root>/wm' */
} ExtY;

/* Real-time Model Data Structure */
struct tag_RTM {
  const char_T * volatile errorStatus;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
  } Timing;
};

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* External data declarations for dependent source files */
extern const real_T motor_RGND;        /* real_T ground */

/* Constant parameters (default storage) */
extern const ConstP rtConstP;

/* Model entry point functions */
extern void motor_initialize(void);
extern void motor_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S32>/0 4' : Unused code path elimination
 * Block '<S32>/1//Ron' : Unused code path elimination
 * Block '<S32>/Saturation' : Unused code path elimination
 * Block '<S32>/Switch' : Unused code path elimination
 * Block '<S32>/Unit Delay' : Unused code path elimination
 * Block '<S33>/Switch' : Unused code path elimination
 * Block '<S33>/Vf Devices & Clamping Diodes' : Unused code path elimination
 * Block '<S33>/Vf Diodes' : Unused code path elimination
 * Block '<S5>/Gain Vr_Vs' : Eliminated nontunable gain of 1
 * Block '<S5>/Gain Vr_Vs1' : Eliminated nontunable gain of 1
 * Block '<S32>/Data Type Conversion' : Eliminate redundant data type conversion
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
 * '<Root>' : 'motor'
 * '<S1>'   : 'motor/Asynchronous Machine SI Units'
 * '<S2>'   : 'motor/DC Voltage Source2'
 * '<S3>'   : 'motor/Universal Bridge1'
 * '<S4>'   : 'motor/powergui'
 * '<S5>'   : 'motor/Asynchronous Machine SI Units/Electrical model'
 * '<S6>'   : 'motor/Asynchronous Machine SI Units/Measurements'
 * '<S7>'   : 'motor/Asynchronous Machine SI Units/Mechanical model'
 * '<S8>'   : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model'
 * '<S9>'   : 'motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation'
 * '<S10>'  : 'motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation'
 * '<S11>'  : 'motor/Asynchronous Machine SI Units/Electrical model/sin,cos'
 * '<S12>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Electromagnetic Torque'
 * '<S13>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Flux Prediction'
 * '<S14>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation'
 * '<S15>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR'
 * '<S16>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Laq=Lad'
 * '<S17>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Matrix L'
 * '<S18>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/phimqd'
 * '<S19>'  : 'motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR/Subsystem'
 * '<S20>'  : 'motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Rotor reference frame'
 * '<S21>'  : 'motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Stationary reference frame'
 * '<S22>'  : 'motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Synchronous reference frame'
 * '<S23>'  : 'motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/transit'
 * '<S24>'  : 'motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Rotor reference frame'
 * '<S25>'  : 'motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Stationary reference frame'
 * '<S26>'  : 'motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Synchronous reference frame'
 * '<S27>'  : 'motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(beta),cos(beta),sin(th),cos(th)'
 * '<S28>'  : 'motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)'
 * '<S29>'  : 'motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)1'
 * '<S30>'  : 'motor/Asynchronous Machine SI Units/Mechanical model/Delay Prediction'
 * '<S31>'  : 'motor/DC Voltage Source2/Model'
 * '<S32>'  : 'motor/Universal Bridge1/Model'
 * '<S33>'  : 'motor/Universal Bridge1/Model/Vf 1'
 * '<S34>'  : 'motor/powergui/EquivalentModel1'
 * '<S35>'  : 'motor/powergui/EquivalentModel1/Gates'
 * '<S36>'  : 'motor/powergui/EquivalentModel1/Sources'
 * '<S37>'  : 'motor/powergui/EquivalentModel1/Status'
 * '<S38>'  : 'motor/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_motor_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
