/*
 * filterbankFIR_private.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "filterbankFIR".
 *
 * Model version              : 1.11
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C source code generated on : Fri Oct 19 11:23:20 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_filterbankFIR_private_h_
#define RTW_HEADER_filterbankFIR_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "filterbankFIR.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#include "dsp_rt.h"                    /* DSP System Toolbox general run time support functions */

extern real_T rt_powd_snf(real_T u0, real_T u1);
extern void RandSrcInitState_GZ(const uint32_T seed[], uint32_T state[], int32_T
  nChans);
extern void RandSrc_GZ_D(real_T y[], const real_T mean[], int32_T meanLen, const
  real_T xstd[], int32_T xstdLen, uint32_T state[], int32_T nChans, int32_T
  nSamps);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern void MWDSPCG_FFT_Interleave_R2BR_D(const real_T x[], creal_T y[], int32_T
  nChans, int32_T nRows);
extern void MWDSPCG_R2DIT_TBLS_Z(creal_T y[], int32_T nChans, int32_T nRows,
  int32_T fftLen, int32_T offset, const real_T tablePtr[], int32_T twiddleStep,
  boolean_T isInverse);
extern void MWDSPCG_FFT_DblLen_Z(creal_T y[], int32_T nChans, int32_T nRows,
  const real_T twiddleTable[], int32_T twiddleStep);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern void filterbankFIR_Compressor2_Init(DW_Compressor2_filterbankFIR_T
  *localDW);
extern void filterbankFIR_Compressor2_Start(DW_Compressor2_filterbankFIR_T
  *localDW, P_Compressor2_filterbankFIR_T *localP);
extern void filterbankFIR_Compressor2(const real_T rtu_0[512],
  B_Compressor2_filterbankFIR_T *localB, DW_Compressor2_filterbankFIR_T *localDW,
  P_Compressor2_filterbankFIR_T *localP);
extern void filterbankFIR_Compressor2_Term(DW_Compressor2_filterbankFIR_T
  *localDW);

#endif                                 /* RTW_HEADER_filterbankFIR_private_h_ */
