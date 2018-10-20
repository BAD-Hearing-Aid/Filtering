/*
 * filterbankFIR_types.h
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

#ifndef RTW_HEADER_filterbankFIR_types_h_
#define RTW_HEADER_filterbankFIR_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#ifndef typedef_audio_simulink_DynamicRangeCo_T
#define typedef_audio_simulink_DynamicRangeCo_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  real_T pSampleRateInherit;
  real_T Threshold;
  real_T AttackTime;
  real_T ReleaseTime;
  real_T pNumChannels;
  real_T pAlphaA;
  real_T pAlphaR;
  real_T pLevelDetectionState;
  real_T MakeUpGain;
  real_T KneeWidth;
  real_T pMakeUpGain;
  real_T Ratio;
} audio_simulink_DynamicRangeCo_T;

#endif                                 /*typedef_audio_simulink_DynamicRangeCo_T*/

#ifndef struct_mdfd2164f8fdb796fffdf8db2fd22180f
#define struct_mdfd2164f8fdb796fffdf8db2fd22180f

struct mdfd2164f8fdb796fffdf8db2fd22180f
{
  int32_T S0_isInitialized;
  real_T W0_states[30];
  real_T P0_InitialStates;
  real_T P1_Coefficients[31];
};

#endif                                 /*struct_mdfd2164f8fdb796fffdf8db2fd22180f*/

#ifndef typedef_dsp_FIRFilter_0_filterbankFIR_T
#define typedef_dsp_FIRFilter_0_filterbankFIR_T

typedef struct mdfd2164f8fdb796fffdf8db2fd22180f dsp_FIRFilter_0_filterbankFIR_T;

#endif                                 /*typedef_dsp_FIRFilter_0_filterbankFIR_T*/

#ifndef struct_mdbdff4a3c1f76c9fac00de1992c764099
#define struct_mdbdff4a3c1f76c9fac00de1992c764099

struct mdbdff4a3c1f76c9fac00de1992c764099
{
  int32_T S0_isInitialized;
  real_T W0_states[3];
  real_T P0_InitialStates;
  real_T P1_Coefficients[4];
};

#endif                                 /*struct_mdbdff4a3c1f76c9fac00de1992c764099*/

#ifndef typedef_dsp_FIRFilter_0_filterbankF_n_T
#define typedef_dsp_FIRFilter_0_filterbankF_n_T

typedef struct mdbdff4a3c1f76c9fac00de1992c764099
  dsp_FIRFilter_0_filterbankF_n_T;

#endif                                 /*typedef_dsp_FIRFilter_0_filterbankF_n_T*/

#ifndef struct_md5lPSJr8Lv2wDv1jvJuDcwF
#define struct_md5lPSJr8Lv2wDv1jvJuDcwF

struct md5lPSJr8Lv2wDv1jvJuDcwF
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  dsp_FIRFilter_0_filterbankF_n_T cSFunObject;
};

#endif                                 /*struct_md5lPSJr8Lv2wDv1jvJuDcwF*/

#ifndef typedef_dspcodegen_FIRFilter_filterba_T
#define typedef_dspcodegen_FIRFilter_filterba_T

typedef struct md5lPSJr8Lv2wDv1jvJuDcwF dspcodegen_FIRFilter_filterba_T;

#endif                                 /*typedef_dspcodegen_FIRFilter_filterba_T*/

#ifndef typedef_cell_wrap_filterbankFIR_T
#define typedef_cell_wrap_filterbankFIR_T

typedef struct {
  uint32_T f1[8];
} cell_wrap_filterbankFIR_T;

#endif                                 /*typedef_cell_wrap_filterbankFIR_T*/

#ifndef struct_mdklfPOfagD0w0TY11WZaESF
#define struct_mdklfPOfagD0w0TY11WZaESF

struct mdklfPOfagD0w0TY11WZaESF
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  dsp_FIRFilter_0_filterbankFIR_T cSFunObject;
  real_T Numerator[31];
};

#endif                                 /*struct_mdklfPOfagD0w0TY11WZaESF*/

#ifndef typedef_dspcodegen_FIRFilter_filter_n_T
#define typedef_dspcodegen_FIRFilter_filter_n_T

typedef struct mdklfPOfagD0w0TY11WZaESF dspcodegen_FIRFilter_filter_n_T;

#endif                                 /*typedef_dspcodegen_FIRFilter_filter_n_T*/

#ifndef typedef_dsp_SpectrumEstimator_filterb_T
#define typedef_dsp_SpectrumEstimator_filterb_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  real_T pPeriodogramMatrix[5120];
  real_T pW[257];
  real_T pWindowData[512];
  real_T pWindowPower;
  real_T pNumAvgsCounter;
  real_T pNewPeriodogramIdx;
  real_T ReferenceLoad;
  real_T pReferenceLoad;
} dsp_SpectrumEstimator_filterb_T;

#endif                                 /*typedef_dsp_SpectrumEstimator_filterb_T*/

#ifndef struct_tag_syFXq48HXjkgxE7hLpgsSCB
#define struct_tag_syFXq48HXjkgxE7hLpgsSCB

struct tag_syFXq48HXjkgxE7hLpgsSCB
{
  real_T fgrid[32];
  real_T des[32];
  real_T wt[32];
  creal_T H[32];
  real_T error[32];
  real_T iextr[3];
  real_T fextr[3];
};

#endif                                 /*struct_tag_syFXq48HXjkgxE7hLpgsSCB*/

#ifndef typedef_syFXq48HXjkgxE7hLpgsSCB_filte_T
#define typedef_syFXq48HXjkgxE7hLpgsSCB_filte_T

typedef struct tag_syFXq48HXjkgxE7hLpgsSCB syFXq48HXjkgxE7hLpgsSCB_filte_T;

#endif                                 /*typedef_syFXq48HXjkgxE7hLpgsSCB_filte_T*/

#ifndef typedef_dsp_CrossSpectrumEstimator_fi_T
#define typedef_dsp_CrossSpectrumEstimator_fi_T

typedef struct {
  int32_T isInitialized;
  boolean_T isSetupComplete;
  creal_T pPeriodogramMatrix[5120];
  real_T pWindowData[512];
  real_T pWindowPower;
  real_T pNumAvgsCounter;
  real_T pNewPeriodogramIdx;
} dsp_CrossSpectrumEstimator_fi_T;

#endif                                 /*typedef_dsp_CrossSpectrumEstimator_fi_T*/

#ifndef typedef_dsp_simulink_TransferFunction_T
#define typedef_dsp_simulink_TransferFunction_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  dsp_SpectrumEstimator_filterb_T pSpectrumEstimator;
  dsp_CrossSpectrumEstimator_fi_T pCrossSpectrumEstimator;
  real_T pFrameCounter;
  real_T pFrameDelay;
} dsp_simulink_TransferFunction_T;

#endif                                 /*typedef_dsp_simulink_TransferFunction_T*/

#ifndef struct_mdgWHtPYTNyNSQTANgBr3RzE
#define struct_mdgWHtPYTNyNSQTANgBr3RzE

struct mdgWHtPYTNyNSQTANgBr3RzE
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  cell_wrap_filterbankFIR_T inputVarSize;
  dspcodegen_FIRFilter_filter_n_T pfilter;
  int8_T pNumChannels;
  real_T CutoffFrequency;
};

#endif                                 /*struct_mdgWHtPYTNyNSQTANgBr3RzE*/

#ifndef typedef_dsp_simulink_VariableBandwidt_T
#define typedef_dsp_simulink_VariableBandwidt_T

typedef struct mdgWHtPYTNyNSQTANgBr3RzE dsp_simulink_VariableBandwidt_T;

#endif                                 /*typedef_dsp_simulink_VariableBandwidt_T*/

#ifndef struct_md59aPrGTCjfgMH5P3jG0gsH
#define struct_md59aPrGTCjfgMH5P3jG0gsH

struct md59aPrGTCjfgMH5P3jG0gsH
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  cell_wrap_filterbankFIR_T inputVarSize;
  int32_T NumChannels;
  dspcodegen_FIRFilter_filterba_T *FilterObj;
};

#endif                                 /*struct_md59aPrGTCjfgMH5P3jG0gsH*/

#ifndef typedef_dsp_Differentiator_filterbank_T
#define typedef_dsp_Differentiator_filterbank_T

typedef struct md59aPrGTCjfgMH5P3jG0gsH dsp_Differentiator_filterbank_T;

#endif                                 /*typedef_dsp_Differentiator_filterbank_T*/

#ifndef typedef_dsp_private_PhaseDifferentiat_T
#define typedef_dsp_private_PhaseDifferentiat_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  cell_wrap_filterbankFIR_T inputVarSize;
  creal_T pLastSample;
  int32_T pNumChans;
} dsp_private_PhaseDifferentiat_T;

#endif                                 /*typedef_dsp_private_PhaseDifferentiat_T*/

#ifndef typedef_dsp_PhaseExtractor_filterbank_T
#define typedef_dsp_PhaseExtractor_filterbank_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  cell_wrap_filterbankFIR_T inputVarSize;
  dsp_private_PhaseDifferentiat_T pPhaseDifferentiator;
} dsp_PhaseExtractor_filterbank_T;

#endif                                 /*typedef_dsp_PhaseExtractor_filterbank_T*/

/* Parameters for system: '<S4>/Compressor2' */
typedef struct P_Compressor2_filterbankFIR_T_ P_Compressor2_filterbankFIR_T;

/* Parameters (default storage) */
typedef struct P_filterbankFIR_T_ P_filterbankFIR_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_filterbankFIR_T RT_MODEL_filterbankFIR_T;

#endif                                 /* RTW_HEADER_filterbankFIR_types_h_ */
