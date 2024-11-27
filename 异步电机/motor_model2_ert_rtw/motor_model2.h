/*
 * File: motor_model2.h
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

#ifndef RTW_HEADER_motor_model2_h_
#define RTW_HEADER_motor_model2_h_
#ifndef motor_model2_COMMON_INCLUDES_
#define motor_model2_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* motor_model2_COMMON_INCLUDES_ */

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real_T ib[4];                        /* '<S18>/ib' */
  real_T StateSpace_o1[8];             /* '<S63>/State-Space' */
  real_T StateSpace_o2[6];             /* '<S63>/State-Space' */
  real_T VectorConcatenate[6];         /* '<S44>/Vector Concatenate' */
  real_T W21wr[16];                    /* '<S36>/W(2,1)=-wr' */
  real_T Constant_f[2];                /* '<S36>/Constant' */
  real_T fluxes_DSTATE[4];             /* '<S23>/fluxes' */
  real_T fluxes_DSTATE_d[4];           /* '<S21>/fluxes' */
  real_T voltages_DSTATE[4];           /* '<S23>/voltages' */
  real_T inversion_DWORK4[16];         /* '<S27>/inversion' */
  real_T inversion_DWORK4_e[16];       /* '<S22>/inversion' */
  real_T Gain1_m;                      /* '<S2>/Gain1' */
  real_T Gain1_k;                      /* '<S3>/Gain1' */
  real_T Gain1_f;                      /* '<S4>/Gain1' */
  real_T Rotorspeedwm;                 /* '<S15>/Rotor speed(wm)' */
  real_T TrigonometricFunction_o1_e;   /* '<S36>/Trigonometric Function' */
  real_T TrigonometricFunction_o2_l;   /* '<S36>/Trigonometric Function' */
  real_T ira;                          /* '<S34>/ira' */
  real_T irb;                          /* '<S34>/irb' */
  real_T isa;                          /* '<S34>/isa' */
  real_T isb;                          /* '<S34>/isb' */
  real_T ira_e;                        /* '<S32>/ira' */
  real_T irb_f;                        /* '<S32>/irb' */
  real_T isa_i;                        /* '<S32>/isa' */
  real_T isb_e;                        /* '<S32>/isb' */
  real_T vdr;                          /* '<S30>/vdr' */
  real_T vds;                          /* '<S30>/vds' */
  real_T vqr;                          /* '<S30>/vqr' */
  real_T vqs;                          /* '<S30>/vqs' */
  real_T vdr_i;                        /* '<S28>/vdr' */
  real_T vds_n;                        /* '<S28>/vds' */
  real_T vqr_p;                        /* '<S28>/vqr' */
  real_T vqs_g;                        /* '<S28>/vqs' */
  real_T Rotoranglethetam_DSTATE;      /* '<S15>/Rotor angle thetam' */
  real_T wm_delay_DSTATE;              /* '<S38>/wm_delay' */
  real_T wm_predict_DSTATE;            /* '<S38>/wm_predict' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S2>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_p;/* '<S3>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_f;/* '<S4>/Discrete-Time Integrator' */
  real_T Rotorspeedwm_DSTATE;          /* '<S15>/Rotor speed(wm)' */
  int_T StateSpace_IWORK[11];          /* '<S63>/State-Space' */
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
  } StateSpace_PWORK;                  /* '<S63>/State-Space' */

  uint8_T Rotorspeedwm_SYSTEM_ENABLE;  /* '<S15>/Rotor speed(wm)' */
  boolean_T sinthrcosthr_MODE;         /* '<S19>/sin(thr),cos(thr)' */
  boolean_T Synchronousreferenceframe_MODE;/* '<S18>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE;  /* '<S18>/Rotor reference frame' */
  boolean_T Synchronousreferenceframe_MOD_a;/* '<S17>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE_o;/* '<S17>/Rotor reference frame' */
} DW;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: SM.Linv
   * Referenced by: '<S16>/Constant2'
   */
  real_T Constant2_Value[16];

  /* Expression: SM.RLinv
   * Referenced by: '<S16>/Constant4'
   */
  real_T Constant4_Value[16];

  /* Expression: eye(4,4)
   * Referenced by: '<S27>/u5'
   */
  real_T u5_Value_m[16];

  /* Expression: S.D
   * Referenced by: '<S63>/State-Space'
   */
  real_T StateSpace_DS_param[72];

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S53>/Gain3'
   */
  real_T Gain3_Gain[9];

  /* Expression: rep_seq_y
   * Referenced by: '<S49>/Look-Up Table1'
   */
  real_T LookUpTable1_tableData[3];

  /* Expression: rep_seq_t - min(rep_seq_t)
   * Referenced by: '<S49>/Look-Up Table1'
   */
  real_T LookUpTable1_bp01Data[3];
} ConstP;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T wish_speed;                   /* '<Root>/wish_speed' */
  real_T load_torque;                  /* '<Root>/load_torque' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T is_a;                         /* '<Root>/is_a' */
  real_T is_b;                         /* '<Root>/is_b' */
  real_T is_c;                         /* '<Root>/is_c' */
  real_T v_d;                          /* '<Root>/v_d' */
  real_T v_q;                          /* '<Root>/v_q' */
  real_T wm;                           /* '<Root>/wm' */
} ExtY;

/* Real-time Model Data Structure */
struct tag_RTM {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    SimTimeStep simTimeStep;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* External data declarations for dependent source files */
extern const real_T motor_model2_RGND; /* real_T ground */

/* Constant parameters (default storage) */
extern const ConstP rtConstP;

/* Model entry point functions */
extern void motor_model2_initialize(void);
extern void motor_model2_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/Constant' : Unused code path elimination
 * Block '<Root>/Constant1' : Unused code path elimination
 * Block '<Root>/Scope' : Unused code path elimination
 * Block '<Root>/Scope1' : Unused code path elimination
 * Block '<S1>/Gain' : Unused code path elimination
 * Block '<S7>/0' : Unused code path elimination
 * Block '<S1>/Scope' : Unused code path elimination
 * Block '<S1>/Scope1' : Unused code path elimination
 * Block '<S50>/0 4' : Unused code path elimination
 * Block '<S50>/1//Ron' : Unused code path elimination
 * Block '<S50>/Saturation' : Unused code path elimination
 * Block '<S50>/Switch' : Unused code path elimination
 * Block '<S50>/Unit Delay' : Unused code path elimination
 * Block '<S51>/Switch' : Unused code path elimination
 * Block '<S51>/Vf Devices & Clamping Diodes' : Unused code path elimination
 * Block '<S51>/Vf Diodes' : Unused code path elimination
 * Block '<S13>/Gain Vr_Vs' : Eliminated nontunable gain of 1
 * Block '<S13>/Gain Vr_Vs1' : Eliminated nontunable gain of 1
 * Block '<S49>/Output' : Eliminate redundant signal conversion block
 * Block '<S50>/Data Type Conversion' : Eliminate redundant data type conversion
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
 * '<Root>' : 'motor_model2'
 * '<S1>'   : 'motor_model2/Subsystem'
 * '<S2>'   : 'motor_model2/Subsystem/ACR_d'
 * '<S3>'   : 'motor_model2/Subsystem/ACR_q'
 * '<S4>'   : 'motor_model2/Subsystem/ASR'
 * '<S5>'   : 'motor_model2/Subsystem/Asynchronous Machine SI Units'
 * '<S6>'   : 'motor_model2/Subsystem/DC Voltage Source2'
 * '<S7>'   : 'motor_model2/Subsystem/Park to Clarke Angle Transform'
 * '<S8>'   : 'motor_model2/Subsystem/SVPWM'
 * '<S9>'   : 'motor_model2/Subsystem/Universal Bridge1'
 * '<S10>'  : 'motor_model2/Subsystem/abc to dq0'
 * '<S11>'  : 'motor_model2/Subsystem/feedforward'
 * '<S12>'  : 'motor_model2/Subsystem/powergui'
 * '<S13>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model'
 * '<S14>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Measurements'
 * '<S15>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Mechanical model'
 * '<S16>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model'
 * '<S17>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation'
 * '<S18>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/dq to abc transformation'
 * '<S19>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/sin,cos'
 * '<S20>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Electromagnetic Torque'
 * '<S21>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Flux Prediction'
 * '<S22>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation'
 * '<S23>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR'
 * '<S24>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Laq=Lad'
 * '<S25>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Matrix L'
 * '<S26>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/phimqd'
 * '<S27>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR/Subsystem'
 * '<S28>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Rotor reference frame'
 * '<S29>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Stationary reference frame'
 * '<S30>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Synchronous reference frame'
 * '<S31>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/transit'
 * '<S32>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Rotor reference frame'
 * '<S33>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Stationary reference frame'
 * '<S34>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Synchronous reference frame'
 * '<S35>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(beta),cos(beta),sin(th),cos(th)'
 * '<S36>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)'
 * '<S37>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)1'
 * '<S38>'  : 'motor_model2/Subsystem/Asynchronous Machine SI Units/Mechanical model/Delay Prediction'
 * '<S39>'  : 'motor_model2/Subsystem/DC Voltage Source2/Model'
 * '<S40>'  : 'motor_model2/Subsystem/SVPWM/Sector'
 * '<S41>'  : 'motor_model2/Subsystem/SVPWM/T1T2'
 * '<S42>'  : 'motor_model2/Subsystem/SVPWM/allocate'
 * '<S43>'  : 'motor_model2/Subsystem/SVPWM/fuse'
 * '<S44>'  : 'motor_model2/Subsystem/SVPWM/gating'
 * '<S45>'  : 'motor_model2/Subsystem/SVPWM/xyz'
 * '<S46>'  : 'motor_model2/Subsystem/SVPWM/Sector/准则A'
 * '<S47>'  : 'motor_model2/Subsystem/SVPWM/Sector/准则B'
 * '<S48>'  : 'motor_model2/Subsystem/SVPWM/Sector/准则C'
 * '<S49>'  : 'motor_model2/Subsystem/SVPWM/gating/triangle'
 * '<S50>'  : 'motor_model2/Subsystem/Universal Bridge1/Model'
 * '<S51>'  : 'motor_model2/Subsystem/Universal Bridge1/Model/Vf 1'
 * '<S52>'  : 'motor_model2/Subsystem/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S53>'  : 'motor_model2/Subsystem/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S54>'  : 'motor_model2/Subsystem/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S55>'  : 'motor_model2/Subsystem/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S56>'  : 'motor_model2/Subsystem/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S57>'  : 'motor_model2/Subsystem/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S58>'  : 'motor_model2/Subsystem/feedforward/Alpha-Beta-Zero to dq0'
 * '<S59>'  : 'motor_model2/Subsystem/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S60>'  : 'motor_model2/Subsystem/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S61>'  : 'motor_model2/Subsystem/feedforward/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S62>'  : 'motor_model2/Subsystem/feedforward/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S63>'  : 'motor_model2/Subsystem/powergui/EquivalentModel1'
 * '<S64>'  : 'motor_model2/Subsystem/powergui/EquivalentModel1/Gates'
 * '<S65>'  : 'motor_model2/Subsystem/powergui/EquivalentModel1/Sources'
 * '<S66>'  : 'motor_model2/Subsystem/powergui/EquivalentModel1/Status'
 * '<S67>'  : 'motor_model2/Subsystem/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_motor_model2_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
