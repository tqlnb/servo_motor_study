/*
 * File: sumAll.h
 *
 * Code generated for Simulink model 'sumAll'.
 *
 * Model version                  : 1.6
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Thu Oct 10 21:22:43 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_sumAll_h_
#define RTW_HEADER_sumAll_h_
#ifndef sumAll_COMMON_INCLUDES_
#define sumAll_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* sumAll_COMMON_INCLUDES_ */

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
  real_T ib[4];                        /* '<S40>/ib' */
  real_T StateSpace_o1[8];             /* '<S64>/State-Space' */
  real_T StateSpace_o2[6];             /* '<S64>/State-Space' */
  real_T VectorConcatenate[6];         /* '<S15>/Vector Concatenate' */
  real_T W21wr[16];                    /* '<S58>/W(2,1)=-wr' */
  real_T Constant_n[2];                /* '<S58>/Constant' */
  real_T fluxes_DSTATE[4];             /* '<S45>/fluxes' */
  real_T fluxes_DSTATE_m[4];           /* '<S43>/fluxes' */
  real_T voltages_DSTATE[4];           /* '<S45>/voltages' */
  real_T inversion_DWORK4[16];         /* '<S49>/inversion' */
  real_T inversion_DWORK4_c[16];       /* '<S44>/inversion' */
  real_T Clock;                        /* '<S10>/Clock' */
  real_T LookUpTable;                  /* '<S20>/Look-Up Table' */
  real_T TrigonometricFunction_o1_a;   /* '<S58>/Trigonometric Function' */
  real_T TrigonometricFunction_o2_a;   /* '<S58>/Trigonometric Function' */
  real_T ira;                          /* '<S56>/ira' */
  real_T irb;                          /* '<S56>/irb' */
  real_T isa;                          /* '<S56>/isa' */
  real_T isb;                          /* '<S56>/isb' */
  real_T ira_d;                        /* '<S54>/ira' */
  real_T irb_f;                        /* '<S54>/irb' */
  real_T isa_p;                        /* '<S54>/isa' */
  real_T isb_j;                        /* '<S54>/isb' */
  real_T vdr;                          /* '<S52>/vdr' */
  real_T vds;                          /* '<S52>/vds' */
  real_T vqr;                          /* '<S52>/vqr' */
  real_T vqs;                          /* '<S52>/vqs' */
  real_T vdr_g;                        /* '<S50>/vdr' */
  real_T vds_k;                        /* '<S50>/vds' */
  real_T vqr_o;                        /* '<S50>/vqr' */
  real_T vqs_m;                        /* '<S50>/vqs' */
  real_T Rotoranglethetam_DSTATE;      /* '<S37>/Rotor angle thetam' */
  real_T wm_delay_DSTATE;              /* '<S60>/wm_delay' */
  real_T wm_predict_DSTATE;            /* '<S60>/wm_predict' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S4>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_l;/* '<S5>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_b;/* '<S6>/Discrete-Time Integrator' */
  real_T Rotorspeedwm_DSTATE;          /* '<S37>/Rotor speed(wm)' */
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
  } StateSpace_PWORK;                  /* '<S64>/State-Space' */

  int_T StateSpace_IWORK[11];          /* '<S64>/State-Space' */
  uint8_T Rotorspeedwm_SYSTEM_ENABLE;  /* '<S37>/Rotor speed(wm)' */
  boolean_T sinthrcosthr_MODE;         /* '<S41>/sin(thr),cos(thr)' */
  boolean_T Synchronousreferenceframe_MODE;/* '<S40>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE;  /* '<S40>/Rotor reference frame' */
  boolean_T Synchronousreferenceframe_MOD_g;/* '<S39>/Synchronous reference frame' */
  boolean_T Rotorreferenceframe_MODE_b;/* '<S39>/Rotor reference frame' */
} DW;

/* Constant parameters (default storage) */
typedef struct {
  /* Expression: SM.Linv
   * Referenced by: '<S38>/Constant2'
   */
  real_T Constant2_Value[16];

  /* Expression: SM.RLinv
   * Referenced by: '<S38>/Constant4'
   */
  real_T Constant4_Value[16];

  /* Expression: eye(4,4)
   * Referenced by: '<S49>/u5'
   */
  real_T u5_Value_p[16];

  /* Expression: S.D
   * Referenced by: '<S64>/State-Space'
   */
  real_T StateSpace_DS_param[72];

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S22>/Gain3'
   */
  real_T Gain3_Gain[9];

  /* Expression: rep_seq_t
   * Referenced by: '<S20>/Look-Up Table'
   */
  real_T LookUpTable_XData[3];

  /* Expression: rep_seq_y
   * Referenced by: '<S20>/Look-Up Table'
   */
  real_T LookUpTable_YData[3];
} ConstP;

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
    uint32_T clockTick1;
    struct {
      uint32_T TID[2];
    } TaskCounters;
  } Timing;
};

/* Block signals and states (default storage) */
extern DW rtDW;

/* External data declarations for dependent source files */
extern const real_T sumAll_RGND;       /* real_T ground */

/* Constant parameters (default storage) */
extern const ConstP rtConstP;

/* Model entry point functions */
extern void sumAll_initialize(void);
extern void sumAll_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S1>/Gain' : Unused code path elimination
 * Block '<S7>/0' : Unused code path elimination
 * Block '<S1>/Scope' : Unused code path elimination
 * Block '<S1>/Scope1' : Unused code path elimination
 * Block '<S1>/Scope2' : Unused code path elimination
 * Block '<S62>/0 4' : Unused code path elimination
 * Block '<S62>/1//Ron' : Unused code path elimination
 * Block '<S62>/Saturation' : Unused code path elimination
 * Block '<S62>/Switch' : Unused code path elimination
 * Block '<S62>/Unit Delay' : Unused code path elimination
 * Block '<S63>/Switch' : Unused code path elimination
 * Block '<S63>/Vf Devices & Clamping Diodes' : Unused code path elimination
 * Block '<S63>/Vf Diodes' : Unused code path elimination
 * Block '<S35>/Gain Vr_Vs' : Eliminated nontunable gain of 1
 * Block '<S35>/Gain Vr_Vs1' : Eliminated nontunable gain of 1
 * Block '<S62>/Data Type Conversion' : Eliminate redundant data type conversion
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
 * '<Root>' : 'sumAll'
 * '<S1>'   : 'sumAll/ctrl'
 * '<S2>'   : 'sumAll/motor'
 * '<S3>'   : 'sumAll/powergui'
 * '<S4>'   : 'sumAll/ctrl/ACR_d'
 * '<S5>'   : 'sumAll/ctrl/ACR_q'
 * '<S6>'   : 'sumAll/ctrl/ASR'
 * '<S7>'   : 'sumAll/ctrl/Park to Clarke Angle Transform'
 * '<S8>'   : 'sumAll/ctrl/SVPWM'
 * '<S9>'   : 'sumAll/ctrl/abc to dq0'
 * '<S10>'  : 'sumAll/ctrl/feedforward'
 * '<S11>'  : 'sumAll/ctrl/SVPWM/Sector'
 * '<S12>'  : 'sumAll/ctrl/SVPWM/T1T2'
 * '<S13>'  : 'sumAll/ctrl/SVPWM/allocate'
 * '<S14>'  : 'sumAll/ctrl/SVPWM/fuse'
 * '<S15>'  : 'sumAll/ctrl/SVPWM/gating'
 * '<S16>'  : 'sumAll/ctrl/SVPWM/xyz'
 * '<S17>'  : 'sumAll/ctrl/SVPWM/Sector/准则A'
 * '<S18>'  : 'sumAll/ctrl/SVPWM/Sector/准则B'
 * '<S19>'  : 'sumAll/ctrl/SVPWM/Sector/准则C'
 * '<S20>'  : 'sumAll/ctrl/SVPWM/gating/triangle'
 * '<S21>'  : 'sumAll/ctrl/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S22>'  : 'sumAll/ctrl/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S23>'  : 'sumAll/ctrl/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S24>'  : 'sumAll/ctrl/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S25>'  : 'sumAll/ctrl/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S26>'  : 'sumAll/ctrl/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S27>'  : 'sumAll/ctrl/feedforward/Alpha-Beta-Zero to dq0'
 * '<S28>'  : 'sumAll/ctrl/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S29>'  : 'sumAll/ctrl/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S30>'  : 'sumAll/ctrl/feedforward/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S31>'  : 'sumAll/ctrl/feedforward/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S32>'  : 'sumAll/motor/Asynchronous Machine SI Units'
 * '<S33>'  : 'sumAll/motor/DC Voltage Source2'
 * '<S34>'  : 'sumAll/motor/Universal Bridge1'
 * '<S35>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model'
 * '<S36>'  : 'sumAll/motor/Asynchronous Machine SI Units/Measurements'
 * '<S37>'  : 'sumAll/motor/Asynchronous Machine SI Units/Mechanical model'
 * '<S38>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model'
 * '<S39>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation'
 * '<S40>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation'
 * '<S41>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/sin,cos'
 * '<S42>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Electromagnetic Torque'
 * '<S43>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Flux Prediction'
 * '<S44>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation'
 * '<S45>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR'
 * '<S46>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Laq=Lad'
 * '<S47>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Matrix L'
 * '<S48>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/phimqd'
 * '<S49>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR/Subsystem'
 * '<S50>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Rotor reference frame'
 * '<S51>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Stationary reference frame'
 * '<S52>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Synchronous reference frame'
 * '<S53>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/transit'
 * '<S54>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Rotor reference frame'
 * '<S55>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Stationary reference frame'
 * '<S56>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Synchronous reference frame'
 * '<S57>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(beta),cos(beta),sin(th),cos(th)'
 * '<S58>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)'
 * '<S59>'  : 'sumAll/motor/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)1'
 * '<S60>'  : 'sumAll/motor/Asynchronous Machine SI Units/Mechanical model/Delay Prediction'
 * '<S61>'  : 'sumAll/motor/DC Voltage Source2/Model'
 * '<S62>'  : 'sumAll/motor/Universal Bridge1/Model'
 * '<S63>'  : 'sumAll/motor/Universal Bridge1/Model/Vf 1'
 * '<S64>'  : 'sumAll/powergui/EquivalentModel1'
 * '<S65>'  : 'sumAll/powergui/EquivalentModel1/Gates'
 * '<S66>'  : 'sumAll/powergui/EquivalentModel1/Sources'
 * '<S67>'  : 'sumAll/powergui/EquivalentModel1/Status'
 * '<S68>'  : 'sumAll/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_sumAll_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
