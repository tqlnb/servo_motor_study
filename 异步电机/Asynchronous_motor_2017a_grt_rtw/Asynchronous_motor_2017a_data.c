/*
 * Asynchronous_motor_2017a_data.c
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

#include "Asynchronous_motor_2017a.h"

/* Block parameters (default storage) */
P_Asynchronous_motor_2017a_T Asynchronous_motor_2017a_P = {
  /* Variable: Ids_rated
   * Referenced by: '<Root>/Constant'
   */
  9.4040441629117311,

  /* Variable: Kac
   * Referenced by:
   *   '<S1>/Gain2'
   *   '<S2>/Gain2'
   */
  0.636801285969825,

  /* Variable: Kas
   * Referenced by: '<S3>/Gain3'
   */
  0.72152170507471169,

  /* Variable: Kic
   * Referenced by:
   *   '<S1>/Gain1'
   *   '<S2>/Gain1'
   */
  289.61527294964543,

  /* Variable: Kis
   * Referenced by: '<S3>/Gain1'
   */
  24.626907165410849,

  /* Variable: Kpc
   * Referenced by:
   *   '<S1>/Gain'
   *   '<S2>/Gain'
   */
  1.5703485875928103,

  /* Variable: Kps
   * Referenced by: '<S3>/Gain'
   */
  1.3859596918105916,

  /* Variable: Kt
   * Referenced by: '<S3>/Gain2'
   */
  1.6153968022055993,

  /* Variable: Lm
   * Referenced by:
   *   '<S10>/Gain4'
   *   '<S10>/Gain8'
   */
  0.059,

  /* Variable: Lr
   * Referenced by: '<S10>/Gain8'
   */
  0.060793999999999994,

  /* Variable: Ls
   * Referenced by:
   *   '<S10>/Gain3'
   *   '<S10>/Gain7'
   */
  0.060793999999999994,

  /* Variable: Rr
   * Referenced by: '<S10>/Gain4'
   */
  0.379,

  /* Variable: Te_rated
   * Referenced by: '<S3>/Saturation'
   */
  48.400544337535294,

  /* Variable: Tsw
   * Referenced by:
   *   '<S10>/Switch'
   *   '<S42>/Constant'
   */
  0.001,

  /* Variable: V_max
   * Referenced by:
   *   '<S1>/Saturation'
   *   '<S2>/Saturation'
   */
  179.62924780409975,

  /* Variable: Vdc
   * Referenced by: '<S38>/DC'
   */
  500.0,

  /* Variable: np
   * Referenced by: '<S10>/Gain1'
   */
  2.0,

  /* Variable: pxyz
   * Referenced by:
   *   '<S44>/Gain'
   *   '<S44>/Gain3'
   *   '<S44>/Gain4'
   */
  3.4641016151377543E-6,

  /* Variable: sigma
   * Referenced by:
   *   '<S10>/Gain2'
   *   '<S10>/Gain6'
   */
  0.058148172073290927,

  /* Variable: tau_r
   * Referenced by: '<S10>/Gain'
   */
  0.16040633245382585,

  /* Mask Parameter: AlphaBetaZerotodq0_Alignment
   * Referenced by: '<S51>/Constant'
   */
  2.0,

  /* Mask Parameter: AlphaBetaZerotodq0_Alignment_b
   * Referenced by: '<S57>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const
   * Referenced by: '<S53>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const
   * Referenced by: '<S54>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const_l
   * Referenced by: '<S58>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const_e
   * Referenced by: '<S59>/Constant'
   */
  2.0,

  /* Mask Parameter: A_const
   * Referenced by: '<S45>/Constant'
   */
  0.0,

  /* Mask Parameter: B_const
   * Referenced by: '<S46>/Constant'
   */
  0.0,

  /* Mask Parameter: C_const
   * Referenced by: '<S47>/Constant'
   */
  0.0,

  /* Mask Parameter: triangle_rep_seq_y
   * Referenced by: '<S48>/Look-Up Table1'
   */
  { 0.0, 0.0005, 0.0 },

  /* Expression: SM.Lsat(1)
   * Referenced by: '<S21>/Constant1'
   */
  0.0,

  /* Computed Parameter: Linv_Y0
   * Referenced by: '<S21>/Linv'
   */
  0.0,

  /* Computed Parameter: RLinv_Y0
   * Referenced by: '<S21>/R*Linv'
   */
  0.0,

  /* Computed Parameter: Lm_Y0
   * Referenced by: '<S21>/Lm'
   */
  0.0,

  /* Expression: [1/SM.Lls 1/SM.Llr]
   * Referenced by: '<S25>/u1'
   */
  { 23.209796154291286, 23.209796154291286 },

  /* Expression: [ SM.Lls SM.Llr ]
   * Referenced by: '<S23>/u2'
   */
  { 0.043085255611566793, 0.043085255611566793 },

  /* Expression: SM.Lsat(1)
   * Referenced by: '<S21>/Delay'
   */
  0.0,

  /* Expression: [ 0 SM.Phisat(2:end)./SM.Lsat(2:end) ]
   * Referenced by: '<S21>/1-D Lookup Table'
   */
  { 0.0, 1.0 },

  /* Expression: SM.Phisat
   * Referenced by: '<S21>/1-D Lookup Table'
   */
  { 0.0, 1.0 },

  /* Expression: zeros(4,4)
   * Referenced by: '<S24>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  /* Expression: SM.Ll
   * Referenced by: '<S24>/u5'
   */
  { 0.043085255611566793, 0.0, 0.0, 0.0, 0.0, 0.043085255611566793, 0.0, 0.0,
    0.0, 0.0, 0.043085255611566793, 0.0, 0.0, 0.0, 0.0, 0.043085255611566793 },

  /* Expression: SM.R
   * Referenced by: '<S21>/u1'
   */
  { 0.022551652892561985, 0.0, 0.0, 0.0, 0.0, 0.022551652892561985, 0.0, 0.0,
    0.0, 0.0, 0.028973140495867768, 0.0, 0.0, 0.0, 0.0, 0.028973140495867768 },

  /* Expression: SM.Lm
   * Referenced by: '<S15>/Lm_nosat'
   */
  1.4169621410715947,

  /* Expression: SM.Linv
   * Referenced by: '<S15>/Constant2'
   */
  { 11.778689645591433, 0.0, -11.431106508699781, -0.0, 0.0, 11.778689645591433,
    0.0, -11.431106508699781, -11.431106508699781, 0.0, 11.778689645591433, -0.0,
    -0.0, -11.431106508699781, -0.0, 11.778689645591433 },

  /* Expression: SM.RLinv
   * Referenced by: '<S15>/Constant4'
   */
  { 0.26562892041659192, 0.0, -0.33119505489978723, 0.0, 0.0,
    0.26562892041659192, 0.0, -0.33119505489978723, -0.25779034616210356, 0.0,
    0.34126562995894349, 0.0, 0.0, -0.25779034616210356, 0.0,
    0.34126562995894349 },

  /* Expression: SM.ensat
   * Referenced by: '<S15>/Constant3'
   */
  0.0,

  /* Expression: 0.5
   * Referenced by: '<S15>/Switch1'
   */
  0.5,

  /* Expression: SM.ctrl
   * Referenced by: '<S18>/Constant4'
   */
  2.0,

  /* Expression: SM.web*Ts/2
   * Referenced by: '<S26>/wbase*Ts//2'
   */
  0.015707963267948967,

  /* Expression: eye(4,4)
   * Referenced by: '<S26>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  /* Expression: SM.web*Ts/2
   * Referenced by: '<S26>/wbase*Ts//2 '
   */
  0.015707963267948967,

  /* Expression: 0
   * Referenced by: '<S27>/vqr,vdr'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S27>/vqs,vds'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S28>/vqr,vdr'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S28>/vqs,vds'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S29>/vqr,vdr'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S29>/vqs,vds'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S31>/ira,irb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S31>/isa,isb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S32>/ira,irb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S32>/isa,isb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S33>/ira,irb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S33>/isa,isb'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S34>/sin(beta),cos(beta), sin(th),cos(th)'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S34>/W'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S34>/we'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S34>/Gain2'
   */
  -1.0,

  /* Expression: SM.web
   * Referenced by: '<S34>/web_psb'
   */
  314.15926535897933,

  /* Expression: [ 0 1  0  0; -1  0  0  0;  0  0  0  0;  0  0  0  0]
   * Referenced by: '<S34>/u3'
   */
  { 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  /* Expression: 0
   * Referenced by: '<S35>/sin(thr),cos(thr)'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S35>/W'
   */
  0.0,

  /* Expression: [0; 0]
   * Referenced by: '<S35>/Constant'
   */
  { 0.0, 0.0 },

  /* Expression: -1
   * Referenced by: '<S35>/Gain1'
   */
  -1.0,

  /* Expression: zeros(4,4)
   * Referenced by: '<S35>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  /* Expression: 0
   * Referenced by: '<S36>/sin(thr),cos(thr)'
   */
  0.0,

  /* Computed Parameter: W_Y0_o
   * Referenced by: '<S36>/W'
   */
  0.0,

  /* Expression: [0; 0]
   * Referenced by: '<S36>/Constant'
   */
  { 0.0, 0.0 },

  /* Expression: -1
   * Referenced by: '<S36>/Gain3'
   */
  -1.0,

  /* Expression: zeros(4,4)
   * Referenced by: '<S36>/u4'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  /* Expression: [0,0]
   * Referenced by: '<S55>/dq'
   */
  { 0.0, 0.0 },

  /* Expression: [0,0]
   * Referenced by: '<S56>/dq'
   */
  { 0.0, 0.0 },

  /* Expression: [0,0]
   * Referenced by: '<S60>/dq'
   */
  { 0.0, 0.0 },

  /* Expression: [0,0]
   * Referenced by: '<S61>/dq'
   */
  { 0.0, 0.0 },

  /* Expression: SM.ctrl
   * Referenced by: '<S17>/Constant3'
   */
  2.0,

  /* Expression: SM.phiqd0
   * Referenced by: '<S22>/fluxes'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: 2
   * Referenced by: '<S20>/Gain'
   */
  2.0,

  /* Expression: SM.phiqd0
   * Referenced by: '<S20>/fluxes'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: SM.ensat
   * Referenced by: '<S15>/Constant'
   */
  0.0,

  /* Expression: SM.ensat
   * Referenced by: '<S15>/Constant1'
   */
  0.0,

  /* Expression: 0.5
   * Referenced by: '<S15>/Switch'
   */
  0.5,

  /* Expression: SM.ctrl
   * Referenced by: '<S18>/Constant2'
   */
  2.0,

  /* Computed Parameter: Rotoranglethetam_gainval
   * Referenced by: '<S14>/Rotor angle thetam'
   */
  0.0001,

  /* Expression: SM.tho
   * Referenced by: '<S14>/Rotor angle thetam'
   */
  0.0,

  /* Expression: SM.wmo
   * Referenced by: '<S37>/wm_delay'
   */
  0.0,

  /* Expression: 2
   * Referenced by: '<S37>/F2'
   */
  2.0,

  /* Expression: SM.wmo
   * Referenced by: '<S37>/wm_predict'
   */
  0.0,

  /* Expression: SM.ctrl
   * Referenced by: '<S17>/Constant4'
   */
  2.0,

  /* Expression: SM.ib
   * Referenced by: '<S17>/ib'
   */
  13.731987951966302,

  /* Expression: SM.Gain_Vr_Vs
   * Referenced by: '<S12>/Gain Vr_Vs1'
   */
  { 1.0, 1.0 },

  /* Expression: S.D
   * Referenced by: '<S62>/State-Space'
   */
  { -50000.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, -50000.0, 0.0,
    0.0, 0.0, 0.0, -50000.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, 0.0,
    -50000.0, 50000.0, 0.0, 0.0, 50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0,
    0.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    50000.0, -50000.0, 0.0, 50000.0, 50000.0, -50000.0, 0.0, 0.0, -50000.0,
    50000.0, -50000.0, -50000.0, 0.0, 0.0, 50000.0, -50000.0, -50000.0, 50000.0,
    50000.0, -100000.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0 },

  /* Expression: [1 -1]
   * Referenced by: '<S19>/1-1'
   */
  { 1.0, -1.0 },

  /* Expression: 1/SM.p
   * Referenced by: '<S14>/1\p'
   */
  0.5,

  /* Expression: [SM.Nb2;SM.Tb2;1]
   * Referenced by: '<S14>/1\p1'
   */
  { 157.07963267948966, 23.554931577600509, 1.0 },

  /* Expression: 30/pi
   * Referenced by: '<Root>/Gain'
   */
  9.5492965855137211,

  /* Expression: SM.ctrl
   * Referenced by: '<S16>/Constant3'
   */
  2.0,

  /* Expression: [0;0]
   * Referenced by: '<S30>/Constant6'
   */
  { 0.0, 0.0 },

  /* Expression: SM.Gain_Vr_Vs
   * Referenced by: '<S12>/Gain Vr_Vs'
   */
  { 1.0, 1.0 },

  /* Expression: 1/SM.Vb
   * Referenced by: '<S16>/1_Vb'
   */
  0.0055670221426890416,

  /* Expression: SM.ctrl
   * Referenced by: '<S16>/Constant4'
   */
  2.0,

  /* Expression: SM.ensat
   * Referenced by: '<S15>/Constant5'
   */
  0.0,

  /* Expression: 0.5
   * Referenced by: '<S15>/Switch2'
   */
  0.5,

  /* Expression: [SM.ib2*SM.kIr*ones(5,1); SM.phib2*SM.kVr; SM.phib2*SM.kVr; SM.Vb2*SM.kVr; SM.Vb2*SM.kVr ; SM.ib2*ones(5,1); SM.phib2; SM.phib2; SM.Vb2; SM.Vb2; SM.phib2/SM.ib2]
   * Referenced by: '<S12>/unit conversion'
   */
  { 13.731987951966302, 13.731987951966302, 13.731987951966302,
    13.731987951966302, 13.731987951966302, 0.57177765423802906,
    0.57177765423802906, 179.62924780409972, 179.62924780409972,
    13.731987951966302, 13.731987951966302, 13.731987951966302,
    13.731987951966302, 13.731987951966302, 0.57177765423802906,
    0.57177765423802906, 179.62924780409972, 179.62924780409972,
    0.041638374300798559 },

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S52>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: 2/3
   * Referenced by: '<S52>/Gain1'
   */
  0.66666666666666663,

  /* Expression: 1 / (Lr ^ 2)
   * Referenced by: '<S10>/Gain5'
   */
  270.56932718377169,

  /* Expression: 0
   * Referenced by: '<S1>/Integrator'
   */
  0.0,

  /* Expression: 0.4
   * Referenced by: '<Root>/Step'
   */
  0.4,

  /* Expression: 500
   * Referenced by: '<Root>/Step'
   */
  500.0,

  /* Expression: 1400
   * Referenced by: '<Root>/Step'
   */
  1400.0,

  /* Expression: pi/30
   * Referenced by: '<Root>/Gain1'
   */
  0.10471975511965977,

  /* Expression: 0
   * Referenced by: '<S3>/Integrator'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S2>/Integrator'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S22>/voltages'
   */
  0.0,

  /* Expression: Ts
   * Referenced by: '<S22>/IC'
   */
  0.0001,

  /* Expression: 0.7
   * Referenced by: '<Root>/Step1'
   */
  0.7,

  /* Expression: 1
   * Referenced by: '<Root>/Step1'
   */
  1.0,

  /* Expression: 40
   * Referenced by: '<Root>/Step1'
   */
  40.0,

  /* Expression: 1/SM.Tb2
   * Referenced by: '<S14>/Unit conversion'
   */
  0.042453954778240446,

  /* Expression: SM.F
   * Referenced by: '<S14>/F'
   */
  0.0,

  /* Expression: 1/(2*SM.H)
   * Referenced by: '<S14>/1_2H'
   */
  5.9506091980420592,

  /* Computed Parameter: Rotorspeedwm_gainval
   * Referenced by: '<S14>/Rotor speed(wm)'
   */
  5.0E-5,

  /* Expression: SM.wmo
   * Referenced by: '<S14>/Rotor speed(wm)'
   */
  0.0,

  /* Expression: SM.web
   * Referenced by: '<S14>/web_psb'
   */
  314.15926535897933,

  /* Expression: sqrt(3)
   * Referenced by: '<S39>/Gain'
   */
  1.7320508075688772,

  /* Expression: -1
   * Referenced by: '<S40>/Gain'
   */
  -1.0,

  /* Expression: 0.5
   * Referenced by: '<S44>/Gain2'
   */
  0.5,

  /* Expression: sqrt(3) / 2
   * Referenced by: '<S44>/Gain1'
   */
  0.8660254037844386,

  /* Expression: -1
   * Referenced by: '<S40>/Gain1'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<S40>/Gain2'
   */
  -1.0,

  /* Expression: 1 / 4
   * Referenced by: '<S42>/Gain'
   */
  0.25,

  /* Expression: 1 / 2
   * Referenced by: '<S42>/Gain1'
   */
  0.5,

  /* Expression: 1 / 2
   * Referenced by: '<S42>/Gain2'
   */
  0.5,

  /* Expression: period
   * Referenced by: '<S48>/Constant'
   */
  0.001,

  /* Expression: rep_seq_t - min(rep_seq_t)
   * Referenced by: '<S48>/Look-Up Table1'
   */
  { 0.0, 0.0005, 0.001 },

  /* Expression: SM.ctrl==1
   * Referenced by: '<S18>/Constant'
   */
  false,

  /* Expression: SM.ctrl==2
   * Referenced by: '<S18>/Constant1'
   */
  true,

  /* Expression: SM.ctrl==3
   * Referenced by: '<S18>/Constant3'
   */
  false,

  /* Expression: SM.ctrl==1
   * Referenced by: '<S17>/Constant'
   */
  false,

  /* Expression: SM.ctrl==2
   * Referenced by: '<S17>/Constant1'
   */
  true,

  /* Expression: SM.ctrl==3
   * Referenced by: '<S17>/Constant2'
   */
  false,

  /* Expression: SM.ctrl==1
   * Referenced by: '<S16>/Constant'
   */
  false,

  /* Expression: SM.ctrl==2
   * Referenced by: '<S16>/Constant1'
   */
  true,

  /* Expression: SM.ctrl==3
   * Referenced by: '<S16>/Constant2'
   */
  false,

  /* Computed Parameter: Gain3_Gain_c
   * Referenced by: '<S39>/Gain3'
   */
  128U,

  /* Computed Parameter: Gain2_Gain_j
   * Referenced by: '<S39>/Gain2'
   */
  128U,

  /* Computed Parameter: Gain1_Gain_g0
   * Referenced by: '<S39>/Gain1'
   */
  128U
};
