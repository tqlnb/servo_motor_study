/*
 * Asynchronous_motor_2017a.h
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

#ifndef RTW_HEADER_Asynchronous_motor_2017a_h_
#define RTW_HEADER_Asynchronous_motor_2017a_h_
#ifndef Asynchronous_motor_2017a_COMMON_INCLUDES_
#define Asynchronous_motor_2017a_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#endif                           /* Asynchronous_motor_2017a_COMMON_INCLUDES_ */

#include "Asynchronous_motor_2017a_types.h"
#include "rtGetNaN.h"
#include <float.h>
#include <string.h>
#include <stddef.h>
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T GainVr_Vs1[2];                /* '<S12>/Gain Vr_Vs1' */
  real_T StateSpace_o1[8];             /* '<S62>/State-Space' */
  real_T StateSpace_o2[6];             /* '<S62>/State-Space' */
  real_T Sum2;                         /* '<S19>/Sum2' */
  real_T wTethr[3];                    /* '<S14>/1\p1' */
  real_T ComplextoMagnitudeAngle_o2;   /* '<Root>/Complex to Magnitude-Angle' */
  real_T Switch[2];                    /* '<S51>/Switch' */
  real_T Sum2_g;                       /* '<Root>/Sum2' */
  real_T Gain;                         /* '<S1>/Gain' */
  real_T Divide;                       /* '<S10>/Divide' */
  real_T Gain1;                        /* '<S10>/Gain1' */
  real_T Gain5;                        /* '<S10>/Gain5' */
  real_T Gain1_m;                      /* '<S1>/Gain1' */
  real_T Gain8;                        /* '<S10>/Gain8' */
  real_T Gain1_k;                      /* '<S2>/Gain1' */
  real_T Gain1_f;                      /* '<S3>/Gain1' */
  real_T F;                            /* '<S14>/F' */
  real_T u_2H;                         /* '<S14>/1_2H' */
  real_T Rotorspeedwm;                 /* '<S14>/Rotor speed(wm)' */
  real_T VectorConcatenate[6];         /* '<S43>/Vector Concatenate' */
  real_T Fcn;                          /* '<S61>/Fcn' */
  real_T Fcn_b;                        /* '<S60>/Fcn' */
  real_T Fcn_l;                        /* '<S56>/Fcn' */
  real_T Fcn1_n;                       /* '<S56>/Fcn1' */
  real_T Fcn_a;                        /* '<S55>/Fcn' */
  real_T Fcn1_h;                       /* '<S55>/Fcn1' */
  real_T Constant[2];                  /* '<S36>/Constant' */
  real_T TrigonometricFunction_o1;     /* '<S36>/Trigonometric Function' */
  real_T TrigonometricFunction_o2;     /* '<S36>/Trigonometric Function' */
  real_T W43wr[16];                    /* '<S36>/W(4,3)=wr' */
  real_T Constant_f[2];                /* '<S35>/Constant' */
  real_T TrigonometricFunction_o1_e;   /* '<S35>/Trigonometric Function' */
  real_T TrigonometricFunction_o2_l;   /* '<S35>/Trigonometric Function' */
  real_T W21wr[16];                    /* '<S35>/W(2,1)=-wr' */
  real_T TrigonometricFunction_o1_f;   /* '<S34>/Trigonometric Function' */
  real_T TrigonometricFunction_o2_k;   /* '<S34>/Trigonometric Function' */
  real_T TrigonometricFunction1_o1;    /* '<S34>/Trigonometric Function1' */
  real_T TrigonometricFunction1_o2;    /* '<S34>/Trigonometric Function1' */
  real_T W43wr1[16];                   /* '<S34>/W(4,3)=wr-1' */
  real_T ira;                          /* '<S33>/ira' */
  real_T irb;                          /* '<S33>/irb' */
  real_T isa;                          /* '<S33>/isa' */
  real_T isb;                          /* '<S33>/isb' */
  real_T ira_a;                        /* '<S32>/ira' */
  real_T irb_b;                        /* '<S32>/irb' */
  real_T isa_o;                        /* '<S32>/isa' */
  real_T isb_d;                        /* '<S32>/isb' */
  real_T ira_e;                        /* '<S31>/ira' */
  real_T irb_f;                        /* '<S31>/irb' */
  real_T isa_i;                        /* '<S31>/isa' */
  real_T isb_e;                        /* '<S31>/isb' */
  real_T vdr;                          /* '<S29>/vdr' */
  real_T vds;                          /* '<S29>/vds' */
  real_T vqr;                          /* '<S29>/vqr' */
  real_T vqs;                          /* '<S29>/vqs' */
  real_T vdr_d;                        /* '<S28>/vdr' */
  real_T vds_b;                        /* '<S28>/vds' */
  real_T vqr_a;                        /* '<S28>/vqr' */
  real_T vqs_h;                        /* '<S28>/vqs' */
  real_T vdr_i;                        /* '<S27>/vdr' */
  real_T vds_n;                        /* '<S27>/vds' */
  real_T vqr_p;                        /* '<S27>/vqr' */
  real_T vqs_g;                        /* '<S27>/vqs' */
  real_T Switch_m;                     /* '<S21>/Switch' */
  real_T Linv[16];                     /* '<S21>/inversion' */
  real_T RLinv[16];                    /* '<S21>/Product1' */
} B_Asynchronous_motor_2017a_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T fluxes_DSTATE[4];             /* '<S22>/fluxes' */
  real_T fluxes_DSTATE_d[4];           /* '<S20>/fluxes' */
  real_T Rotoranglethetam_DSTATE;      /* '<S14>/Rotor angle thetam' */
  real_T wm_delay_DSTATE;              /* '<S37>/wm_delay' */
  real_T wm_predict_DSTATE;            /* '<S37>/wm_predict' */
  real_T voltages_DSTATE[4];           /* '<S22>/voltages' */
  real_T Rotorspeedwm_DSTATE;          /* '<S14>/Rotor speed(wm)' */
  real_T inversion_DWORK4[16];         /* '<S26>/inversion' */
  real_T inversion_DWORK4_e[16];       /* '<S21>/inversion' */
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
  } StateSpace_PWORK;                  /* '<S62>/State-Space' */

  uint32_T m_bpIndex;                  /* '<S21>/1-D Lookup Table' */
  int_T StateSpace_IWORK[11];          /* '<S62>/State-Space' */
  uint8_T Rotorspeedwm_SYSTEM_ENABLE;  /* '<S14>/Rotor speed(wm)' */
  boolean_T sinthrcosthr1_MODE;        /* '<S18>/sin(thr),cos(thr)1' */
  boolean_T sinthrcosthr_MODE;         /* '<S18>/sin(thr),cos(thr)' */
  boolean_T Synchronousreferenceframe_MODE;/* '<S17>/Synchronous reference frame' */
  boolean_T Stationaryreferenceframe_MODE;/* '<S17>/Stationary reference frame' */
  boolean_T Rotorreferenceframe_MODE;  /* '<S17>/Rotor reference frame' */
  boolean_T Synchronousreferenceframe_MOD_a;/* '<S16>/Synchronous reference frame' */
  boolean_T Stationaryreferenceframe_MODE_l;/* '<S16>/Stationary reference frame' */
  boolean_T Rotorreferenceframe_MODE_o;/* '<S16>/Rotor reference frame' */
} DW_Asynchronous_motor_2017a_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S1>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S3>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S2>/Integrator' */
} X_Asynchronous_motor_2017a_T;

/* State derivatives (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S1>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S3>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S2>/Integrator' */
} XDot_Asynchronous_motor_2017a_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE;         /* '<S1>/Integrator' */
  boolean_T Integrator_CSTATE_b;       /* '<S3>/Integrator' */
  boolean_T Integrator_CSTATE_c;       /* '<S2>/Integrator' */
} XDis_Asynchronous_motor_2017a_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Parameters (default storage) */
struct P_Asynchronous_motor_2017a_T_ {
  real_T Ids_rated;                    /* Variable: Ids_rated
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Kac;                          /* Variable: Kac
                                        * Referenced by:
                                        *   '<S1>/Gain2'
                                        *   '<S2>/Gain2'
                                        */
  real_T Kas;                          /* Variable: Kas
                                        * Referenced by: '<S3>/Gain3'
                                        */
  real_T Kic;                          /* Variable: Kic
                                        * Referenced by:
                                        *   '<S1>/Gain1'
                                        *   '<S2>/Gain1'
                                        */
  real_T Kis;                          /* Variable: Kis
                                        * Referenced by: '<S3>/Gain1'
                                        */
  real_T Kpc;                          /* Variable: Kpc
                                        * Referenced by:
                                        *   '<S1>/Gain'
                                        *   '<S2>/Gain'
                                        */
  real_T Kps;                          /* Variable: Kps
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T Kt;                           /* Variable: Kt
                                        * Referenced by: '<S3>/Gain2'
                                        */
  real_T Lm;                           /* Variable: Lm
                                        * Referenced by:
                                        *   '<S10>/Gain4'
                                        *   '<S10>/Gain8'
                                        */
  real_T Lr;                           /* Variable: Lr
                                        * Referenced by: '<S10>/Gain8'
                                        */
  real_T Ls;                           /* Variable: Ls
                                        * Referenced by:
                                        *   '<S10>/Gain3'
                                        *   '<S10>/Gain7'
                                        */
  real_T Rr;                           /* Variable: Rr
                                        * Referenced by: '<S10>/Gain4'
                                        */
  real_T Te_rated;                     /* Variable: Te_rated
                                        * Referenced by: '<S3>/Saturation'
                                        */
  real_T Tsw;                          /* Variable: Tsw
                                        * Referenced by:
                                        *   '<S10>/Switch'
                                        *   '<S42>/Constant'
                                        */
  real_T V_max;                        /* Variable: V_max
                                        * Referenced by:
                                        *   '<S1>/Saturation'
                                        *   '<S2>/Saturation'
                                        */
  real_T Vdc;                          /* Variable: Vdc
                                        * Referenced by: '<S38>/DC'
                                        */
  real_T np;                           /* Variable: np
                                        * Referenced by: '<S10>/Gain1'
                                        */
  real_T pxyz;                         /* Variable: pxyz
                                        * Referenced by:
                                        *   '<S44>/Gain'
                                        *   '<S44>/Gain3'
                                        *   '<S44>/Gain4'
                                        */
  real_T sigma;                        /* Variable: sigma
                                        * Referenced by:
                                        *   '<S10>/Gain2'
                                        *   '<S10>/Gain6'
                                        */
  real_T tau_r;                        /* Variable: tau_r
                                        * Referenced by: '<S10>/Gain'
                                        */
  real_T AlphaBetaZerotodq0_Alignment;
                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                  * Referenced by: '<S51>/Constant'
                                  */
  real_T AlphaBetaZerotodq0_Alignment_b;
                               /* Mask Parameter: AlphaBetaZerotodq0_Alignment_b
                                * Referenced by: '<S57>/Constant'
                                */
  real_T CompareToConstant_const;     /* Mask Parameter: CompareToConstant_const
                                       * Referenced by: '<S53>/Constant'
                                       */
  real_T CompareToConstant1_const;   /* Mask Parameter: CompareToConstant1_const
                                      * Referenced by: '<S54>/Constant'
                                      */
  real_T CompareToConstant_const_l; /* Mask Parameter: CompareToConstant_const_l
                                     * Referenced by: '<S58>/Constant'
                                     */
  real_T CompareToConstant1_const_e;
                                   /* Mask Parameter: CompareToConstant1_const_e
                                    * Referenced by: '<S59>/Constant'
                                    */
  real_T A_const;                      /* Mask Parameter: A_const
                                        * Referenced by: '<S45>/Constant'
                                        */
  real_T B_const;                      /* Mask Parameter: B_const
                                        * Referenced by: '<S46>/Constant'
                                        */
  real_T C_const;                      /* Mask Parameter: C_const
                                        * Referenced by: '<S47>/Constant'
                                        */
  real_T triangle_rep_seq_y[3];        /* Mask Parameter: triangle_rep_seq_y
                                        * Referenced by: '<S48>/Look-Up Table1'
                                        */
  real_T Constant1_Value;              /* Expression: SM.Lsat(1)
                                        * Referenced by: '<S21>/Constant1'
                                        */
  real_T Linv_Y0;                      /* Computed Parameter: Linv_Y0
                                        * Referenced by: '<S21>/Linv'
                                        */
  real_T RLinv_Y0;                     /* Computed Parameter: RLinv_Y0
                                        * Referenced by: '<S21>/R*Linv'
                                        */
  real_T Lm_Y0;                        /* Computed Parameter: Lm_Y0
                                        * Referenced by: '<S21>/Lm'
                                        */
  real_T u1_Value[2];                  /* Expression: [1/SM.Lls 1/SM.Llr]
                                        * Referenced by: '<S25>/u1'
                                        */
  real_T u2_Value[2];                  /* Expression: [ SM.Lls SM.Llr ]
                                        * Referenced by: '<S23>/u2'
                                        */
  real_T Delay_InitialCondition;       /* Expression: SM.Lsat(1)
                                        * Referenced by: '<S21>/Delay'
                                        */
  real_T uDLookupTable_tableData[2];
                           /* Expression: [ 0 SM.Phisat(2:end)./SM.Lsat(2:end) ]
                            * Referenced by: '<S21>/1-D Lookup Table'
                            */
  real_T uDLookupTable_bp01Data[2];    /* Expression: SM.Phisat
                                        * Referenced by: '<S21>/1-D Lookup Table'
                                        */
  real_T u1_Value_j[16];               /* Expression: zeros(4,4)
                                        * Referenced by: '<S24>/u1'
                                        */
  real_T u5_Value[16];                 /* Expression: SM.Ll
                                        * Referenced by: '<S24>/u5'
                                        */
  real_T u1_Value_e[16];               /* Expression: SM.R
                                        * Referenced by: '<S21>/u1'
                                        */
  real_T Lm_nosat_Value;               /* Expression: SM.Lm
                                        * Referenced by: '<S15>/Lm_nosat'
                                        */
  real_T Constant2_Value[16];          /* Expression: SM.Linv
                                        * Referenced by: '<S15>/Constant2'
                                        */
  real_T Constant4_Value[16];          /* Expression: SM.RLinv
                                        * Referenced by: '<S15>/Constant4'
                                        */
  real_T Constant3_Value;              /* Expression: SM.ensat
                                        * Referenced by: '<S15>/Constant3'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S15>/Switch1'
                                        */
  real_T Constant4_Value_n;            /* Expression: SM.ctrl
                                        * Referenced by: '<S18>/Constant4'
                                        */
  real_T wbaseTs2_Gain;                /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S26>/wbase*Ts//2'
                                        */
  real_T u5_Value_m[16];               /* Expression: eye(4,4)
                                        * Referenced by: '<S26>/u5'
                                        */
  real_T wbaseTs2_Gain_f;              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S26>/wbase*Ts//2 '
                                        */
  real_T vqrvdr_Y0;                    /* Expression: 0
                                        * Referenced by: '<S27>/vqr,vdr'
                                        */
  real_T vqsvds_Y0;                    /* Expression: 0
                                        * Referenced by: '<S27>/vqs,vds'
                                        */
  real_T vqrvdr_Y0_b;                  /* Expression: 0
                                        * Referenced by: '<S28>/vqr,vdr'
                                        */
  real_T vqsvds_Y0_j;                  /* Expression: 0
                                        * Referenced by: '<S28>/vqs,vds'
                                        */
  real_T vqrvdr_Y0_g;                  /* Expression: 0
                                        * Referenced by: '<S29>/vqr,vdr'
                                        */
  real_T vqsvds_Y0_b;                  /* Expression: 0
                                        * Referenced by: '<S29>/vqs,vds'
                                        */
  real_T irairb_Y0;                    /* Expression: 0
                                        * Referenced by: '<S31>/ira,irb'
                                        */
  real_T isaisb_Y0;                    /* Expression: 0
                                        * Referenced by: '<S31>/isa,isb'
                                        */
  real_T irairb_Y0_n;                  /* Expression: 0
                                        * Referenced by: '<S32>/ira,irb'
                                        */
  real_T isaisb_Y0_p;                  /* Expression: 0
                                        * Referenced by: '<S32>/isa,isb'
                                        */
  real_T irairb_Y0_n1;                 /* Expression: 0
                                        * Referenced by: '<S33>/ira,irb'
                                        */
  real_T isaisb_Y0_g;                  /* Expression: 0
                                        * Referenced by: '<S33>/isa,isb'
                                        */
  real_T sinbetacosbetasinthcosth_Y0;  /* Expression: 0
                                        * Referenced by: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
                                        */
  real_T W_Y0;                         /* Expression: 0
                                        * Referenced by: '<S34>/W'
                                        */
  real_T we_Value;                     /* Expression: 1
                                        * Referenced by: '<S34>/we'
                                        */
  real_T Gain2_Gain;                   /* Expression: -1
                                        * Referenced by: '<S34>/Gain2'
                                        */
  real_T web_psb_Gain;                 /* Expression: SM.web
                                        * Referenced by: '<S34>/web_psb'
                                        */
  real_T u3_Value[16];
              /* Expression: [ 0 1  0  0; -1  0  0  0;  0  0  0  0;  0  0  0  0]
               * Referenced by: '<S34>/u3'
               */
  real_T sinthrcosthr_Y0;              /* Expression: 0
                                        * Referenced by: '<S35>/sin(thr),cos(thr)'
                                        */
  real_T W_Y0_d;                       /* Expression: 0
                                        * Referenced by: '<S35>/W'
                                        */
  real_T Constant_Value[2];            /* Expression: [0; 0]
                                        * Referenced by: '<S35>/Constant'
                                        */
  real_T Gain1_Gain;                   /* Expression: -1
                                        * Referenced by: '<S35>/Gain1'
                                        */
  real_T u1_Value_a[16];               /* Expression: zeros(4,4)
                                        * Referenced by: '<S35>/u1'
                                        */
  real_T sinthrcosthr_Y0_n;            /* Expression: 0
                                        * Referenced by: '<S36>/sin(thr),cos(thr)'
                                        */
  real_T W_Y0_o;                       /* Computed Parameter: W_Y0_o
                                        * Referenced by: '<S36>/W'
                                        */
  real_T Constant_Value_a[2];          /* Expression: [0; 0]
                                        * Referenced by: '<S36>/Constant'
                                        */
  real_T Gain3_Gain;                   /* Expression: -1
                                        * Referenced by: '<S36>/Gain3'
                                        */
  real_T u4_Value[16];                 /* Expression: zeros(4,4)
                                        * Referenced by: '<S36>/u4'
                                        */
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S55>/dq'
                                        */
  real_T dq_Y0_k[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S56>/dq'
                                        */
  real_T dq_Y0_h[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S60>/dq'
                                        */
  real_T dq_Y0_m[2];                   /* Expression: [0,0]
                                        * Referenced by: '<S61>/dq'
                                        */
  real_T Constant3_Value_a;            /* Expression: SM.ctrl
                                        * Referenced by: '<S17>/Constant3'
                                        */
  real_T fluxes_InitialCondition[4];   /* Expression: SM.phiqd0
                                        * Referenced by: '<S22>/fluxes'
                                        */
  real_T Gain_Gain;                    /* Expression: 2
                                        * Referenced by: '<S20>/Gain'
                                        */
  real_T fluxes_InitialCondition_p[4]; /* Expression: SM.phiqd0
                                        * Referenced by: '<S20>/fluxes'
                                        */
  real_T Constant_Value_a2;            /* Expression: SM.ensat
                                        * Referenced by: '<S15>/Constant'
                                        */
  real_T Constant1_Value_m;            /* Expression: SM.ensat
                                        * Referenced by: '<S15>/Constant1'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S15>/Switch'
                                        */
  real_T Constant2_Value_n;            /* Expression: SM.ctrl
                                        * Referenced by: '<S18>/Constant2'
                                        */
  real_T Rotoranglethetam_gainval;
                                 /* Computed Parameter: Rotoranglethetam_gainval
                                  * Referenced by: '<S14>/Rotor angle thetam'
                                  */
  real_T Rotoranglethetam_IC;          /* Expression: SM.tho
                                        * Referenced by: '<S14>/Rotor angle thetam'
                                        */
  real_T wm_delay_InitialCondition;    /* Expression: SM.wmo
                                        * Referenced by: '<S37>/wm_delay'
                                        */
  real_T F2_Gain;                      /* Expression: 2
                                        * Referenced by: '<S37>/F2'
                                        */
  real_T wm_predict_InitialCondition;  /* Expression: SM.wmo
                                        * Referenced by: '<S37>/wm_predict'
                                        */
  real_T Constant4_Value_d;            /* Expression: SM.ctrl
                                        * Referenced by: '<S17>/Constant4'
                                        */
  real_T ib_Gain;                      /* Expression: SM.ib
                                        * Referenced by: '<S17>/ib'
                                        */
  real_T GainVr_Vs1_Gain[2];           /* Expression: SM.Gain_Vr_Vs
                                        * Referenced by: '<S12>/Gain Vr_Vs1'
                                        */
  real_T StateSpace_DS_param[72];      /* Expression: S.D
                                        * Referenced by: '<S62>/State-Space'
                                        */
  real_T u1_Gain[2];                   /* Expression: [1 -1]
                                        * Referenced by: '<S19>/1-1'
                                        */
  real_T up_Gain;                      /* Expression: 1/SM.p
                                        * Referenced by: '<S14>/1\p'
                                        */
  real_T up1_Gain[3];                  /* Expression: [SM.Nb2;SM.Tb2;1]
                                        * Referenced by: '<S14>/1\p1'
                                        */
  real_T Gain_Gain_e;                  /* Expression: 30/pi
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Constant3_Value_m;            /* Expression: SM.ctrl
                                        * Referenced by: '<S16>/Constant3'
                                        */
  real_T Constant6_Value[2];           /* Expression: [0;0]
                                        * Referenced by: '<S30>/Constant6'
                                        */
  real_T GainVr_Vs_Gain[2];            /* Expression: SM.Gain_Vr_Vs
                                        * Referenced by: '<S12>/Gain Vr_Vs'
                                        */
  real_T u_Vb_Gain;                    /* Expression: 1/SM.Vb
                                        * Referenced by: '<S16>/1_Vb'
                                        */
  real_T Constant4_Value_h;            /* Expression: SM.ctrl
                                        * Referenced by: '<S16>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: SM.ensat
                                        * Referenced by: '<S15>/Constant5'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S15>/Switch2'
                                        */
  real_T unitconversion_Gain[19];
  /* Expression: [SM.ib2*SM.kIr*ones(5,1); SM.phib2*SM.kVr; SM.phib2*SM.kVr; SM.Vb2*SM.kVr; SM.Vb2*SM.kVr ; SM.ib2*ones(5,1); SM.phib2; SM.phib2; SM.Vb2; SM.Vb2; SM.phib2/SM.ib2]
   * Referenced by: '<S12>/unit conversion'
   */
  real_T Gain3_Gain_o[9];
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S52>/Gain3'
   */
  real_T Gain1_Gain_g;                 /* Expression: 2/3
                                        * Referenced by: '<S52>/Gain1'
                                        */
  real_T Gain5_Gain;                   /* Expression: 1 / (Lr ^ 2)
                                        * Referenced by: '<S10>/Gain5'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S1>/Integrator'
                                        */
  real_T Step_Time;                    /* Expression: 0.4
                                        * Referenced by: '<Root>/Step'
                                        */
  real_T Step_Y0;                      /* Expression: 500
                                        * Referenced by: '<Root>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 1400
                                        * Referenced by: '<Root>/Step'
                                        */
  real_T Gain1_Gain_i;                 /* Expression: pi/30
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T Integrator_IC_f;              /* Expression: 0
                                        * Referenced by: '<S3>/Integrator'
                                        */
  real_T Integrator_IC_b;              /* Expression: 0
                                        * Referenced by: '<S2>/Integrator'
                                        */
  real_T voltages_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S22>/voltages'
                                        */
  real_T IC_Threshold;                 /* Expression: Ts
                                        * Referenced by: '<S22>/IC'
                                        */
  real_T Step1_Time;                   /* Expression: 0.7
                                        * Referenced by: '<Root>/Step1'
                                        */
  real_T Step1_Y0;                     /* Expression: 1
                                        * Referenced by: '<Root>/Step1'
                                        */
  real_T Step1_YFinal;                 /* Expression: 40
                                        * Referenced by: '<Root>/Step1'
                                        */
  real_T Unitconversion_Gain;          /* Expression: 1/SM.Tb2
                                        * Referenced by: '<S14>/Unit conversion'
                                        */
  real_T F_Gain;                       /* Expression: SM.F
                                        * Referenced by: '<S14>/F'
                                        */
  real_T u_2H_Gain;                    /* Expression: 1/(2*SM.H)
                                        * Referenced by: '<S14>/1_2H'
                                        */
  real_T Rotorspeedwm_gainval;       /* Computed Parameter: Rotorspeedwm_gainval
                                      * Referenced by: '<S14>/Rotor speed(wm)'
                                      */
  real_T Rotorspeedwm_IC;              /* Expression: SM.wmo
                                        * Referenced by: '<S14>/Rotor speed(wm)'
                                        */
  real_T web_psb_Gain_l;               /* Expression: SM.web
                                        * Referenced by: '<S14>/web_psb'
                                        */
  real_T Gain_Gain_ep;                 /* Expression: sqrt(3)
                                        * Referenced by: '<S39>/Gain'
                                        */
  real_T Gain_Gain_o;                  /* Expression: -1
                                        * Referenced by: '<S40>/Gain'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: 0.5
                                        * Referenced by: '<S44>/Gain2'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: sqrt(3) / 2
                                        * Referenced by: '<S44>/Gain1'
                                        */
  real_T Gain1_Gain_n;                 /* Expression: -1
                                        * Referenced by: '<S40>/Gain1'
                                        */
  real_T Gain2_Gain_o;                 /* Expression: -1
                                        * Referenced by: '<S40>/Gain2'
                                        */
  real_T Gain_Gain_n;                  /* Expression: 1 / 4
                                        * Referenced by: '<S42>/Gain'
                                        */
  real_T Gain1_Gain_ii;                /* Expression: 1 / 2
                                        * Referenced by: '<S42>/Gain1'
                                        */
  real_T Gain2_Gain_d;                 /* Expression: 1 / 2
                                        * Referenced by: '<S42>/Gain2'
                                        */
  real_T Constant_Value_o;             /* Expression: period
                                        * Referenced by: '<S48>/Constant'
                                        */
  real_T LookUpTable1_bp01Data[3];     /* Expression: rep_seq_t - min(rep_seq_t)
                                        * Referenced by: '<S48>/Look-Up Table1'
                                        */
  boolean_T Constant_Value_c;          /* Expression: SM.ctrl==1
                                        * Referenced by: '<S18>/Constant'
                                        */
  boolean_T Constant1_Value_p;         /* Expression: SM.ctrl==2
                                        * Referenced by: '<S18>/Constant1'
                                        */
  boolean_T Constant3_Value_k;         /* Expression: SM.ctrl==3
                                        * Referenced by: '<S18>/Constant3'
                                        */
  boolean_T Constant_Value_c3;         /* Expression: SM.ctrl==1
                                        * Referenced by: '<S17>/Constant'
                                        */
  boolean_T Constant1_Value_j;         /* Expression: SM.ctrl==2
                                        * Referenced by: '<S17>/Constant1'
                                        */
  boolean_T Constant2_Value_a;         /* Expression: SM.ctrl==3
                                        * Referenced by: '<S17>/Constant2'
                                        */
  boolean_T Constant_Value_i;          /* Expression: SM.ctrl==1
                                        * Referenced by: '<S16>/Constant'
                                        */
  boolean_T Constant1_Value_mg;        /* Expression: SM.ctrl==2
                                        * Referenced by: '<S16>/Constant1'
                                        */
  boolean_T Constant2_Value_g;         /* Expression: SM.ctrl==3
                                        * Referenced by: '<S16>/Constant2'
                                        */
  uint8_T Gain3_Gain_c;                /* Computed Parameter: Gain3_Gain_c
                                        * Referenced by: '<S39>/Gain3'
                                        */
  uint8_T Gain2_Gain_j;                /* Computed Parameter: Gain2_Gain_j
                                        * Referenced by: '<S39>/Gain2'
                                        */
  uint8_T Gain1_Gain_g0;               /* Computed Parameter: Gain1_Gain_g0
                                        * Referenced by: '<S39>/Gain1'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Asynchronous_motor_20_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_Asynchronous_motor_2017a_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_Asynchronous_motor_2017a_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[3];
  real_T odeF[3][3];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_P;

/* Block signals (default storage) */
extern B_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_B;

/* Continuous states (default storage) */
extern X_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_X;

/* Block states (default storage) */
extern DW_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_DW;

/* External data declarations for dependent source files */
extern const real_T Asynchronous_motor_2017a_RGND;/* real_T ground */

/* Model entry point functions */
extern void Asynchronous_motor_2017a_initialize(void);
extern void Asynchronous_motor_2017a_step(void);
extern void Asynchronous_motor_2017a_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Asynchronous_motor_2_T *const Asynchronous_motor_2017a_M;

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
 * '<Root>' : 'Asynchronous_motor_2017a'
 * '<S1>'   : 'Asynchronous_motor_2017a/ACR_d'
 * '<S2>'   : 'Asynchronous_motor_2017a/ACR_q'
 * '<S3>'   : 'Asynchronous_motor_2017a/ASR'
 * '<S4>'   : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units'
 * '<S5>'   : 'Asynchronous_motor_2017a/DC Voltage Source2'
 * '<S6>'   : 'Asynchronous_motor_2017a/Park to Clarke Angle Transform'
 * '<S7>'   : 'Asynchronous_motor_2017a/SVPWM'
 * '<S8>'   : 'Asynchronous_motor_2017a/Universal Bridge1'
 * '<S9>'   : 'Asynchronous_motor_2017a/abc to dq0'
 * '<S10>'  : 'Asynchronous_motor_2017a/feedforward'
 * '<S11>'  : 'Asynchronous_motor_2017a/powergui'
 * '<S12>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model'
 * '<S13>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Measurements'
 * '<S14>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Mechanical model'
 * '<S15>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model'
 * '<S16>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation'
 * '<S17>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/dq to abc transformation'
 * '<S18>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/sin,cos'
 * '<S19>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Electromagnetic Torque'
 * '<S20>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Flux Prediction'
 * '<S21>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation'
 * '<S22>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR'
 * '<S23>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Laq=Lad'
 * '<S24>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/Matrix L'
 * '<S25>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/Saturation/phimqd'
 * '<S26>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/Asynchronous Machine State-space model/phiqd_SR/Subsystem'
 * '<S27>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Rotor reference frame'
 * '<S28>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Stationary reference frame'
 * '<S29>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/Synchronous reference frame'
 * '<S30>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/abc to dq  transformation/transit'
 * '<S31>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Rotor reference frame'
 * '<S32>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Stationary reference frame'
 * '<S33>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/dq to abc transformation/Synchronous reference frame'
 * '<S34>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(beta),cos(beta),sin(th),cos(th)'
 * '<S35>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)'
 * '<S36>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Electrical model/sin,cos/sin(thr),cos(thr)1'
 * '<S37>'  : 'Asynchronous_motor_2017a/Asynchronous Machine SI Units/Mechanical model/Delay Prediction'
 * '<S38>'  : 'Asynchronous_motor_2017a/DC Voltage Source2/Model'
 * '<S39>'  : 'Asynchronous_motor_2017a/SVPWM/Sector'
 * '<S40>'  : 'Asynchronous_motor_2017a/SVPWM/T1T2'
 * '<S41>'  : 'Asynchronous_motor_2017a/SVPWM/allocate'
 * '<S42>'  : 'Asynchronous_motor_2017a/SVPWM/fuse'
 * '<S43>'  : 'Asynchronous_motor_2017a/SVPWM/gating'
 * '<S44>'  : 'Asynchronous_motor_2017a/SVPWM/xyz'
 * '<S45>'  : 'Asynchronous_motor_2017a/SVPWM/Sector/准则A'
 * '<S46>'  : 'Asynchronous_motor_2017a/SVPWM/Sector/准则B'
 * '<S47>'  : 'Asynchronous_motor_2017a/SVPWM/Sector/准则C'
 * '<S48>'  : 'Asynchronous_motor_2017a/SVPWM/gating/triangle'
 * '<S49>'  : 'Asynchronous_motor_2017a/Universal Bridge1/Model'
 * '<S50>'  : 'Asynchronous_motor_2017a/Universal Bridge1/Model/Vf 1'
 * '<S51>'  : 'Asynchronous_motor_2017a/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S52>'  : 'Asynchronous_motor_2017a/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S53>'  : 'Asynchronous_motor_2017a/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S54>'  : 'Asynchronous_motor_2017a/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S55>'  : 'Asynchronous_motor_2017a/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S56>'  : 'Asynchronous_motor_2017a/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S57>'  : 'Asynchronous_motor_2017a/feedforward/Alpha-Beta-Zero to dq0'
 * '<S58>'  : 'Asynchronous_motor_2017a/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S59>'  : 'Asynchronous_motor_2017a/feedforward/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S60>'  : 'Asynchronous_motor_2017a/feedforward/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S61>'  : 'Asynchronous_motor_2017a/feedforward/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S62>'  : 'Asynchronous_motor_2017a/powergui/EquivalentModel1'
 * '<S63>'  : 'Asynchronous_motor_2017a/powergui/EquivalentModel1/Gates'
 * '<S64>'  : 'Asynchronous_motor_2017a/powergui/EquivalentModel1/Sources'
 * '<S65>'  : 'Asynchronous_motor_2017a/powergui/EquivalentModel1/Status'
 * '<S66>'  : 'Asynchronous_motor_2017a/powergui/EquivalentModel1/Yout'
 */
#endif                              /* RTW_HEADER_Asynchronous_motor_2017a_h_ */
