/*
 * filterbankFIR.c
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

#include "filterbankFIR.h"
#include "filterbankFIR_private.h"
#define filterbankFIR_SampleRate       (44100.0)

/* Block signals (default storage) */
B_filterbankFIR_T filterbankFIR_B;

/* Block states (default storage) */
DW_filterbankFIR_T filterbankFIR_DW;

/* Real-time model */
RT_MODEL_filterbankFIR_T filterbankFIR_M_;
RT_MODEL_filterbankFIR_T *const filterbankFIR_M = &filterbankFIR_M_;

/* Forward declaration for local functions */
static void filterbankFIR_abs(const real_T x[512], real_T y[512]);
static void filterbankFIR_log10(const real_T x[512], real_T b_x[512]);
static void filterbankFIR_power(const real_T a_data[], const int32_T *a_size,
  real_T y_data[], int32_T *y_size);
static void filterba_compressor_computeGain(const
  audio_simulink_DynamicRangeCo_T *obj, const real_T xG[512], real_T G[512]);
static void filt_CompressorBase_detectLevel(audio_simulink_DynamicRangeCo_T *obj,
  const real_T x[512], real_T yout[512], B_Compressor2_filterbankFIR_T *localB);
static void filterbankFIR_power_n(const real_T b[512], real_T y[512]);
static void filterbankFIR_SystemCore_step(audio_simulink_DynamicRangeCo_T *obj,
  const real_T varargin_1[512], real_T varargout_1[512],
  B_Compressor2_filterbankFIR_T *localB);
static void SystemProp_matlabCodegenSetAnyP(audio_simulink_DynamicRangeCo_T *obj);
static void filterbankFI_SystemCore_release(audio_simulink_DynamicRangeCo_T *obj);
static void filterbankFIR_SystemCore_delete(audio_simulink_DynamicRangeCo_T *obj);
static void matlabCodegenHandle_matlabCodeg(audio_simulink_DynamicRangeCo_T *obj);

/* Forward declaration for local functions */
static real_T filterbankFIR_sum(const real_T x[31]);
static void filterbankFIR_bsxfun(const real_T a[512], const real_T b[512],
  real_T c[512]);
static void SpectrumEstimatorBase_computeFF(const real_T xin[512], creal_T xout
  [512]);
static void SpectrumEstimator_computeWindow(const
  dsp_SpectrumEstimator_filterb_T *obj, const real_T x[512], real_T Pxx[512]);
static real_T filterbankFIR_mod(real_T x);
static void SpectrumEstimatorBase_updatePer(dsp_SpectrumEstimator_filterb_T *obj,
  const real_T P[512]);
static void SpectrumEstimatorBase_getPeriod(const
  dsp_SpectrumEstimator_filterb_T *obj, real_T P[512]);
static void SpectrumEstimator_convertAndSca(const real_T P[512], real_T Pout[257]);
static void SpectrumEstimatorBase_compute_n(const real_T xin[1024], creal_T
  xout[1024]);
static void CrossSpectrumEstimator_computeW(const
  dsp_CrossSpectrumEstimator_fi_T *obj, const real_T x[512], const real_T y[512],
  creal_T Pxy[512]);
static void SpectrumEstimatorBase_updateP_n(dsp_CrossSpectrumEstimator_fi_T *obj,
  const creal_T P[512]);
static void SpectrumEstimatorBase_getPeri_n(const
  dsp_CrossSpectrumEstimator_fi_T *obj, creal_T P[512]);
static void CrossSpectrumEstimator_convertA(const creal_T P[512], creal_T Pout
  [257]);
static void SystemProp_matlabCodege_nzusynb(dsp_simulink_TransferFunction_T *obj,
  boolean_T value);
static void fi_SystemCore_release_nzusynbuf(dsp_simulink_TransferFunction_T *obj);
static void filte_SystemCore_delete_nzusynb(dsp_simulink_TransferFunction_T *obj);
static void matlabCodegenHandle_mat_nzusynb(dsp_simulink_TransferFunction_T *obj);
static void SystemProp_matlabCodeg_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj,
  boolean_T value);
static void f_SystemCore_release_nzusynbufg(dsp_SpectrumEstimator_filterb_T *obj);
static void filt_SystemCore_delete_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj);
static void matlabCodegenHandle_ma_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj);
static void SystemProp_matlabCodegen_nzusyn(dsp_PhaseExtractor_filterbank_T *obj,
  boolean_T value);
static void filt_SystemCore_release_nzusynb(dsp_PhaseExtractor_filterbank_T *obj);
static void filter_SystemCore_delete_nzusyn(dsp_PhaseExtractor_filterbank_T *obj);
static void matlabCodegenHandle_matl_nzusyn(dsp_PhaseExtractor_filterbank_T *obj);
static void SystemProp_matlabCodegenS_nzusy(dsp_private_PhaseDifferentiat_T *obj,
  boolean_T value);
static void fil_SystemCore_release_nzusynbu(dsp_private_PhaseDifferentiat_T *obj);
static void filterb_SystemCore_delete_nzusy(dsp_private_PhaseDifferentiat_T *obj);
static void matlabCodegenHandle_matla_nzusy(dsp_private_PhaseDifferentiat_T *obj);
static void SystemProp_matlabCodegenSe_nzus(dsp_Differentiator_filterbank_T *obj,
  boolean_T value);
static void filterb_SystemCore_release_nzus(dspcodegen_FIRFilter_filterba_T *obj);
static void filt_Differentiator_releaseImpl(dsp_Differentiator_filterbank_T *obj);
static void SystemCore_releaseWrapper_nzus(dsp_Differentiator_filterbank_T *obj);
static void filter_SystemCore_release_nzusy(dsp_Differentiator_filterbank_T *obj);
static void filterba_SystemCore_delete_nzus(dsp_Differentiator_filterbank_T *obj);
static void matlabCodegenHandle_matlab_nzus(dsp_Differentiator_filterbank_T *obj);
static void SystemProp_matlabCodegenSet_nzu(dspcodegen_FIRFilter_filterba_T *obj,
  boolean_T value);
static void filterban_SystemCore_delete_nzu(dspcodegen_FIRFilter_filterba_T *obj);
static void matlabCodegenHandle_matlabC_nzu(dspcodegen_FIRFilter_filterba_T *obj);
static void SystemProp_matlabCodegenSetA_nz(dsp_simulink_VariableBandwidt_T *obj,
  boolean_T value);
static void filterbank_SystemCore_release_n(dspcodegen_FIRFilter_filter_n_T *obj);
static void VariableBandwidthFilterBase_rel(dsp_simulink_VariableBandwidt_T *obj);
static void filte_SystemCore_releaseWrapper(dsp_simulink_VariableBandwidt_T *obj);
static void filterban_SystemCore_release_nz(dsp_simulink_VariableBandwidt_T *obj);
static void filterbank_SystemCore_delete_nz(dsp_simulink_VariableBandwidt_T *obj);
static void matlabCodegenHandle_matlabCo_nz(dsp_simulink_VariableBandwidt_T *obj);
static void SystemProp_matlabCodegenSetAn_n(dspcodegen_FIRFilter_filter_n_T *obj,
  boolean_T value);
static void filterbankF_SystemCore_delete_n(dspcodegen_FIRFilter_filter_n_T *obj);
static void matlabCodegenHandle_matlabCod_n(dspcodegen_FIRFilter_filter_n_T *obj);
static void filterbankFIR_abs(const real_T x[512], real_T y[512])
{
  int32_T k;
  for (k = 0; k < 512; k++) {
    y[k] = fabs(x[k]);
  }
}

static void filterbankFIR_log10(const real_T x[512], real_T b_x[512])
{
  int32_T k;
  for (k = 0; k < 512; k++) {
    b_x[k] = log10(x[k]);
  }
}

static void filterbankFIR_power(const real_T a_data[], const int32_T *a_size,
  real_T y_data[], int32_T *y_size)
{
  real_T b_z1_data[512];
  int32_T loop_ub;
  if (0 <= *a_size - 1) {
    memcpy(&b_z1_data[0], &y_data[0], *a_size * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub < *a_size; loop_ub++) {
    b_z1_data[loop_ub] = a_data[loop_ub] * a_data[loop_ub];
  }

  *y_size = *a_size;
  if (0 <= *a_size - 1) {
    memcpy(&y_data[0], &b_z1_data[0], *a_size * sizeof(real_T));
  }
}

static void filterba_compressor_computeGain(const
  audio_simulink_DynamicRangeCo_T *obj, const real_T xG[512], real_T G[512])
{
  boolean_T ind2[512];
  real_T c_data[512];
  int16_T d_data[512];
  int16_T e_data[512];
  int32_T trueCount;
  real_T y;
  real_T y_0;
  real_T y_1;
  int32_T i;
  real_T tmp_data[512];
  real_T xG_0[512];
  int32_T d_size_idx_0;
  boolean_T ind2_0;
  trueCount = 0;
  for (i = 0; i < 512; i++) {
    G[i] = xG[i];
    ind2_0 = ((xG[i] - obj->Threshold) * 2.0 > obj->KneeWidth);
    if (ind2_0) {
      trueCount++;
    }

    ind2[i] = ind2_0;
  }

  d_size_idx_0 = trueCount;
  trueCount = 0;
  for (i = 0; i < 512; i++) {
    if (ind2[i]) {
      d_data[trueCount] = (int16_T)(i + 1);
      trueCount++;
    }
  }

  for (trueCount = 0; trueCount < d_size_idx_0; trueCount++) {
    c_data[trueCount] = (xG[d_data[trueCount] - 1] - obj->Threshold) /
      obj->Ratio;
  }

  trueCount = 0;
  for (i = 0; i < 512; i++) {
    if (ind2[i]) {
      G[i] = obj->Threshold + c_data[trueCount];
      trueCount++;
    }
  }

  if (obj->KneeWidth != 0.0) {
    for (trueCount = 0; trueCount < 512; trueCount++) {
      xG_0[trueCount] = xG[trueCount] - obj->Threshold;
    }

    filterbankFIR_abs(xG_0, tmp_data);
    trueCount = 0;
    for (i = 0; i < 512; i++) {
      ind2_0 = (2.0 * tmp_data[i] <= obj->KneeWidth);
      if (ind2_0) {
        trueCount++;
      }

      ind2[i] = ind2_0;
    }

    d_size_idx_0 = trueCount;
    trueCount = 0;
    for (i = 0; i < 512; i++) {
      if (ind2[i]) {
        e_data[trueCount] = (int16_T)(i + 1);
        trueCount++;
      }
    }

    y = 1.0 / obj->Ratio;
    y_0 = obj->KneeWidth / 2.0;
    y_1 = 2.0 * obj->KneeWidth;
    for (trueCount = 0; trueCount < d_size_idx_0; trueCount++) {
      xG_0[trueCount] = (xG[e_data[trueCount] - 1] - obj->Threshold) + y_0;
    }

    filterbankFIR_power(xG_0, &d_size_idx_0, tmp_data, &i);
    for (trueCount = 0; trueCount < i; trueCount++) {
      c_data[trueCount] = (y - 1.0) * tmp_data[trueCount] / y_1;
    }

    trueCount = 0;
    for (i = 0; i < 512; i++) {
      if (ind2[i]) {
        G[i] = xG[i] + c_data[trueCount];
        trueCount++;
      }
    }
  }

  for (trueCount = 0; trueCount < 512; trueCount++) {
    G[trueCount] -= xG[trueCount];
  }
}

static void filt_CompressorBase_detectLevel(audio_simulink_DynamicRangeCo_T *obj,
  const real_T x[512], real_T yout[512], B_Compressor2_filterbankFIR_T *localB)
{
  real_T alphaA;
  real_T alphaR;
  int32_T b_i;
  memset(&localB->y[0], 0, 513U * sizeof(real_T));
  localB->y[0] = obj->pLevelDetectionState;
  alphaA = obj->pAlphaA;
  alphaR = obj->pAlphaR;
  for (b_i = 0; b_i < 512; b_i++) {
    if (x[b_i] <= localB->y[b_i]) {
      localB->y[b_i + 1] = (1.0 - alphaA) * x[b_i] + alphaA * localB->y[b_i];
    } else {
      localB->y[b_i + 1] = (1.0 - alphaR) * x[b_i] + alphaR * localB->y[b_i];
    }

    yout[b_i] = 0.0;
  }

  for (b_i = 0; b_i < 512; b_i++) {
    yout[(b_i + 1) - 1] = localB->y[(b_i + 2) - 1];
  }

  obj->pLevelDetectionState = localB->y[512];
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

static void filterbankFIR_power_n(const real_T b[512], real_T y[512])
{
  int32_T k;
  for (k = 0; k < 512; k++) {
    y[k] = rt_powd_snf(10.0, b[k]);
  }
}

static void filterbankFIR_SystemCore_step(audio_simulink_DynamicRangeCo_T *obj,
  const real_T varargin_1[512], real_T varargout_1[512],
  B_Compressor2_filterbankFIR_T *localB)
{
  real_T Fs;
  real_T varargin_1_0[512];
  real_T b_z1[512];
  int32_T k;
  if (obj->TunablePropsChanged) {
    obj->TunablePropsChanged = false;
    Fs = obj->pSampleRateInherit;
    if (obj->AttackTime != 0.0) {
      obj->pAlphaA = exp(-2.1972245773362196 / (obj->AttackTime *
        obj->pSampleRateInherit));
    } else {
      obj->pAlphaA = 0.0;
    }

    if (obj->ReleaseTime != 0.0) {
      obj->pAlphaR = exp(-2.1972245773362196 / (obj->ReleaseTime * Fs));
    } else {
      obj->pAlphaR = 0.0;
    }

    obj->pMakeUpGain = obj->MakeUpGain;
  }

  filterbankFIR_abs(varargin_1, varargin_1_0);
  for (k = 0; k < 512; k++) {
    if (varargin_1_0[k] > 2.2204460492503131E-16) {
      b_z1[k] = varargin_1_0[k];
    } else {
      b_z1[k] = 2.2204460492503131E-16;
    }
  }

  filterbankFIR_log10(b_z1, varargin_1_0);
  for (k = 0; k < 512; k++) {
    b_z1[k] = 20.0 * varargin_1_0[k];
  }

  filterba_compressor_computeGain(obj, b_z1, varargin_1_0);
  memcpy(&b_z1[0], &varargin_1_0[0], sizeof(real_T) << 9U);
  filt_CompressorBase_detectLevel(obj, b_z1, varargin_1_0, localB);
  for (k = 0; k < 512; k++) {
    b_z1[k] = (varargin_1_0[k] + obj->pMakeUpGain) / 20.0;
  }

  filterbankFIR_power_n(b_z1, varargin_1_0);
  for (k = 0; k < 512; k++) {
    varargout_1[k] = varargin_1[k] * varargin_1_0[k];
  }
}

static void SystemProp_matlabCodegenSetAnyP(audio_simulink_DynamicRangeCo_T *obj)
{
  obj->matlabCodegenIsDeleted = true;
}

static void filterbankFI_SystemCore_release(audio_simulink_DynamicRangeCo_T *obj)
{
  if ((obj->isInitialized == 1) && obj->isSetupComplete) {
    obj->pNumChannels = -1.0;
  }
}

static void filterbankFIR_SystemCore_delete(audio_simulink_DynamicRangeCo_T *obj)
{
  filterbankFI_SystemCore_release(obj);
}

static void matlabCodegenHandle_matlabCodeg(audio_simulink_DynamicRangeCo_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenSetAnyP(obj);
    filterbankFIR_SystemCore_delete(obj);
  }
}

/*
 * System initialize for atomic system:
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 */
void filterbankFIR_Compressor2_Init(DW_Compressor2_filterbankFIR_T *localDW)
{
  /* InitializeConditions for MATLABSystem: '<S4>/Compressor2' */
  localDW->obj.pLevelDetectionState = 0.0;
  if (localDW->obj.AttackTime != 0.0) {
    localDW->obj.pAlphaA = exp(-2.1972245773362196 / (localDW->obj.AttackTime *
      localDW->obj.pSampleRateInherit));
  } else {
    localDW->obj.pAlphaA = 0.0;
  }

  if (localDW->obj.ReleaseTime != 0.0) {
    localDW->obj.pAlphaR = exp(-2.1972245773362196 / (localDW->obj.ReleaseTime *
      localDW->obj.pSampleRateInherit));
  } else {
    localDW->obj.pAlphaR = 0.0;
  }

  /* End of InitializeConditions for MATLABSystem: '<S4>/Compressor2' */
}

/*
 * Start for atomic system:
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 */
void filterbankFIR_Compressor2_Start(DW_Compressor2_filterbankFIR_T *localDW,
  P_Compressor2_filterbankFIR_T *localP)
{
  /* Start for MATLABSystem: '<S4>/Compressor2' */
  localDW->obj.matlabCodegenIsDeleted = true;
  localDW->obj.pNumChannels = -1.0;
  localDW->obj.pSampleRateInherit = -1.0;
  localDW->obj.isInitialized = 0;
  localDW->obj.matlabCodegenIsDeleted = false;
  localDW->objisempty = true;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.Ratio = localP->Compressor2_Ratio;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.Threshold = localP->Compressor2_Threshold;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.KneeWidth = localP->Compressor2_KneeWidth;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.AttackTime = localP->Compressor2_AttackTime;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.ReleaseTime = localP->Compressor2_ReleaseTime;
  if (localDW->obj.isInitialized == 1) {
    localDW->obj.TunablePropsChanged = true;
  }

  localDW->obj.MakeUpGain = localP->Compressor2_MakeUpGain;
  localDW->obj.isSetupComplete = false;
  localDW->obj.isInitialized = 1;
  localDW->obj.pSampleRateInherit = 44100.0;
  localDW->obj.pLevelDetectionState = 0.0;
  if (localDW->obj.AttackTime != 0.0) {
    localDW->obj.pAlphaA = exp(-2.1972245773362196 / (localDW->obj.AttackTime *
      localDW->obj.pSampleRateInherit));
  } else {
    localDW->obj.pAlphaA = 0.0;
  }

  if (localDW->obj.ReleaseTime != 0.0) {
    localDW->obj.pAlphaR = exp(-2.1972245773362196 / (localDW->obj.ReleaseTime *
      localDW->obj.pSampleRateInherit));
  } else {
    localDW->obj.pAlphaR = 0.0;
  }

  localDW->obj.pNumChannels = 1.0;
  localDW->obj.pMakeUpGain = localDW->obj.MakeUpGain;
  localDW->obj.isSetupComplete = true;
  localDW->obj.TunablePropsChanged = false;

  /* End of Start for MATLABSystem: '<S4>/Compressor2' */
}

/*
 * Output and update for atomic system:
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 */
void filterbankFIR_Compressor2(const real_T rtu_0[512],
  B_Compressor2_filterbankFIR_T *localB, DW_Compressor2_filterbankFIR_T *localDW,
  P_Compressor2_filterbankFIR_T *localP)
{
  boolean_T p;
  boolean_T p_0;

  /* MATLABSystem: '<S4>/Compressor2' */
  p = false;
  p_0 = true;
  if (!(localDW->obj.Ratio == localP->Compressor2_Ratio)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.Ratio = localP->Compressor2_Ratio;
  }

  p = false;
  p_0 = true;
  if (!(localDW->obj.Threshold == localP->Compressor2_Threshold)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.Threshold = localP->Compressor2_Threshold;
  }

  p = false;
  p_0 = true;
  if (!(localDW->obj.KneeWidth == localP->Compressor2_KneeWidth)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.KneeWidth = localP->Compressor2_KneeWidth;
  }

  p = false;
  p_0 = true;
  if (!(localDW->obj.AttackTime == localP->Compressor2_AttackTime)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.AttackTime = localP->Compressor2_AttackTime;
  }

  p = false;
  p_0 = true;
  if (!(localDW->obj.ReleaseTime == localP->Compressor2_ReleaseTime)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.ReleaseTime = localP->Compressor2_ReleaseTime;
  }

  p = false;
  p_0 = true;
  if (!(localDW->obj.MakeUpGain == localP->Compressor2_MakeUpGain)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    if (localDW->obj.isInitialized == 1) {
      localDW->obj.TunablePropsChanged = true;
    }

    localDW->obj.MakeUpGain = localP->Compressor2_MakeUpGain;
  }

  filterbankFIR_SystemCore_step(&localDW->obj, rtu_0, localB->Compressor2,
    localB);

  /* End of MATLABSystem: '<S4>/Compressor2' */
}

/*
 * Termination for atomic system:
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 *    synthesized block
 */
void filterbankFIR_Compressor2_Term(DW_Compressor2_filterbankFIR_T *localDW)
{
  /* Terminate for MATLABSystem: '<S4>/Compressor2' */
  matlabCodegenHandle_matlabCodeg(&localDW->obj);
}

static real_T filterbankFIR_sum(const real_T x[31])
{
  real_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 30; k++) {
    y += x[k + 1];
  }

  return y;
}

void RandSrcInitState_GZ(const uint32_T seed[], uint32_T state[], int32_T nChans)
{
  int32_T i;

  /* InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source1' */
  /* RandSrcInitState_GZ */
  for (i = 0; i < nChans; i++) {
    state[i << 1] = 362436069U;
    state[(i << 1) + 1] = seed[i] == 0U ? 521288629U : seed[i];
  }

  /* End of InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source1' */
}

void RandSrc_GZ_D(real_T y[], const real_T mean[], int32_T meanLen, const real_T
                  xstd[], int32_T xstdLen, uint32_T state[], int32_T nChans,
                  int32_T nSamps)
{
  int32_T i;
  int32_T j;
  real_T r;
  real_T x;
  real_T s;
  real_T y_0;
  int32_T chan;
  uint32_T icng;
  uint32_T jsr;
  int32_T samp;
  static const real_T vt[65] = { 0.340945, 0.4573146, 0.5397793, 0.6062427,
    0.6631691, 0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,
    0.9490791, 0.9820005, 1.0138492, 1.044781, 1.0749254, 1.1043917, 1.1332738,
    1.161653, 1.189601, 1.2171815, 1.2444516, 1.2714635, 1.298265, 1.3249008,
    1.3514125, 1.3778399, 1.4042211, 1.4305929, 1.4569915, 1.4834527, 1.5100122,
    1.5367061, 1.5635712, 1.5906454, 1.617968, 1.6455802, 1.6735255, 1.7018503,
    1.7306045, 1.7598422, 1.7896223, 1.8200099, 1.851077, 1.8829044, 1.9155831,
    1.9492166, 1.9839239, 2.0198431, 2.0571356, 2.095993, 2.136645, 2.1793713,
    2.2245175, 2.2725186, 2.3239338, 2.3795008, 2.4402218, 2.5075117, 2.5834658,
    2.6713916, 2.7769942, 2.7769942, 2.7769942, 2.7769942 };

  /* S-Function (sdsprandsrc2): '<Root>/Random Source1' */
  /* RandSrc_GZ_D */
  for (chan = 0; chan < nChans; chan++) {
    icng = state[chan << 1];
    jsr = state[(chan << 1) + 1];
    for (samp = 0; samp < nSamps; samp++) {
      icng = 69069U * icng + 1234567U;
      jsr ^= jsr << 13;
      jsr ^= jsr >> 17;
      jsr ^= jsr << 5;
      i = (int32_T)(icng + jsr);
      j = (i & 63) + 1;
      r = (real_T)i * 4.6566128730773926E-10 * vt[j];
      if (!(fabs(r) <= vt[j - 1])) {
        x = (fabs(r) - vt[j - 1]) / (vt[j] - vt[j - 1]);
        icng = 69069U * icng + 1234567U;
        jsr ^= jsr << 13;
        jsr ^= jsr >> 17;
        jsr ^= jsr << 5;
        y_0 = (real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10 + 0.5;
        s = x + y_0;
        if (s > 1.301198) {
          r = r < 0.0 ? 0.4878992 * x - 0.4878992 : 0.4878992 - 0.4878992 * x;
        } else {
          if (!(s <= 0.9689279)) {
            x = 0.4878992 - 0.4878992 * x;
            if (y_0 > 12.67706 - exp(-0.5 * x * x) * 12.37586) {
              r = r < 0.0 ? -x : x;
            } else {
              if (!(exp(-0.5 * vt[j] * vt[j]) + y_0 * 0.01958303 / vt[j] <= exp(
                    -0.5 * r * r))) {
                do {
                  icng = 69069U * icng + 1234567U;
                  jsr ^= jsr << 13;
                  jsr ^= jsr >> 17;
                  jsr ^= jsr << 5;
                  x = log((real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10
                          + 0.5) / 2.776994;
                  icng = 69069U * icng + 1234567U;
                  jsr ^= jsr << 13;
                  jsr ^= jsr >> 17;
                  jsr ^= jsr << 5;
                } while (log((real_T)(int32_T)(icng + jsr) *
                             2.328306436538696E-10 + 0.5) * -2.0 <= x * x);

                r = r < 0.0 ? x - 2.776994 : 2.776994 - x;
              }
            }
          }
        }
      }

      y[chan * nSamps + samp] = xstd[xstdLen > 1 ? chan : 0] * r + mean[meanLen >
        1 ? chan : 0];
    }

    state[chan << 1] = icng;
    state[(chan << 1) + 1] = jsr;
  }

  /* End of S-Function (sdsprandsrc2): '<Root>/Random Source1' */
}

static void filterbankFIR_bsxfun(const real_T a[512], const real_T b[512],
  real_T c[512])
{
  int32_T k;
  for (k = 0; k < 512; k++) {
    c[k] = a[k] * b[k];
  }
}

static void SpectrumEstimatorBase_computeFF(const real_T xin[512], creal_T xout
  [512])
{
  int32_T nRowsD4;
  int32_T ix;
  int32_T ju;
  int32_T iy;
  int32_T iheight;
  int32_T istart;
  int32_T j;
  int32_T ihi;
  boolean_T tst;
  static const real_T tmp[257] = { 1.0, 0.9999247018391445, 0.99969881869620425,
    0.99932238458834954, 0.99879545620517241, 0.99811811290014918,
    0.99729045667869021, 0.996312612182778, 0.99518472667219693,
    0.99390697000235606, 0.99247953459871, 0.99090263542778, 0.989176509964781,
    0.98730141815785843, 0.98527764238894122, 0.98310548743121629,
    0.98078528040323043, 0.97831737071962765, 0.97570213003852857,
    0.97293995220556018, 0.970031253194544, 0.96697647104485207,
    0.96377606579543984, 0.96043051941556579, 0.95694033573220882,
    0.95330604035419386, 0.94952818059303667, 0.94560732538052128,
    0.94154406518302081, 0.937339011912575, 0.932992798834739,
    0.92850608047321559, 0.92387953251128674, 0.91911385169005777,
    0.91420975570353069, 0.90916798309052238, 0.90398929312344334,
    0.89867446569395382, 0.89322430119551532, 0.88763962040285393,
    0.881921264348355, 0.8760700941954066, 0.87008699110871146,
    0.8639728561215867, 0.85772861000027212, 0.8513551931052652,
    0.84485356524970712, 0.83822470555483808, 0.83146961230254524,
    0.82458930278502529, 0.81758481315158371, 0.81045719825259477,
    0.80320753148064494, 0.79583690460888357, 0.78834642762660634,
    0.78073722857209449, 0.773010453362737, 0.765167265622459,
    0.75720884650648457, 0.74913639452345937, 0.74095112535495922,
    0.73265427167241282, 0.724247082951467, 0.71573082528381859,
    0.70710678118654757, 0.69837624940897292, 0.68954054473706683,
    0.680600997795453, 0.67155895484701833, 0.66241577759017178,
    0.65317284295377676, 0.64383154288979139, 0.63439328416364549,
    0.62485948814238634, 0.61523159058062682, 0.60551104140432555,
    0.59569930449243336, 0.58579785745643886, 0.57580819141784534,
    0.56573181078361312, 0.55557023301960218, 0.54532498842204646,
    0.53499761988709715, 0.524589682678469, 0.51410274419322166,
    0.50353838372571758, 0.49289819222978404, 0.48218377207912272,
    0.47139673682599764, 0.46053871095824, 0.44961132965460654,
    0.43861623853852766, 0.42755509343028208, 0.41642956009763715,
    0.40524131400498986, 0.3939920400610481, 0.38268343236508978,
    0.37131719395183749, 0.35989503653498811, 0.34841868024943456,
    0.33688985339222005, 0.32531029216226293, 0.31368174039889152,
    0.30200594931922808, 0.29028467725446233, 0.27851968938505306,
    0.26671275747489837, 0.25486565960451457, 0.24298017990326387,
    0.23105810828067111, 0.2191012401568698, 0.20711137619221856,
    0.19509032201612825, 0.18303988795514095, 0.17096188876030122,
    0.15885814333386145, 0.14673047445536175, 0.13458070850712617,
    0.1224106751992162, 0.11022220729388306, 0.0980171403295606,
    0.0857973123444399, 0.073564563599667426, 0.061320736302208578,
    0.049067674327418015, 0.036807222941358832, 0.024541228522912288,
    0.012271538285719925, 0.0, -0.012271538285719925, -0.024541228522912288,
    -0.036807222941358832, -0.049067674327418015, -0.061320736302208578,
    -0.073564563599667426, -0.0857973123444399, -0.0980171403295606,
    -0.11022220729388306, -0.1224106751992162, -0.13458070850712617,
    -0.14673047445536175, -0.15885814333386145, -0.17096188876030122,
    -0.18303988795514095, -0.19509032201612825, -0.20711137619221856,
    -0.2191012401568698, -0.23105810828067111, -0.24298017990326387,
    -0.25486565960451457, -0.26671275747489837, -0.27851968938505306,
    -0.29028467725446233, -0.30200594931922808, -0.31368174039889152,
    -0.32531029216226293, -0.33688985339222005, -0.34841868024943456,
    -0.35989503653498811, -0.37131719395183749, -0.38268343236508978,
    -0.3939920400610481, -0.40524131400498986, -0.41642956009763715,
    -0.42755509343028208, -0.43861623853852766, -0.44961132965460654,
    -0.46053871095824, -0.47139673682599764, -0.48218377207912272,
    -0.49289819222978404, -0.50353838372571758, -0.51410274419322166,
    -0.524589682678469, -0.53499761988709715, -0.54532498842204646,
    -0.55557023301960218, -0.56573181078361312, -0.57580819141784534,
    -0.58579785745643886, -0.59569930449243336, -0.60551104140432555,
    -0.61523159058062682, -0.62485948814238634, -0.63439328416364549,
    -0.64383154288979139, -0.65317284295377676, -0.66241577759017178,
    -0.67155895484701833, -0.680600997795453, -0.68954054473706683,
    -0.69837624940897292, -0.70710678118654757, -0.71573082528381859,
    -0.724247082951467, -0.73265427167241282, -0.74095112535495922,
    -0.74913639452345937, -0.75720884650648457, -0.765167265622459,
    -0.773010453362737, -0.78073722857209449, -0.78834642762660634,
    -0.79583690460888357, -0.80320753148064494, -0.81045719825259477,
    -0.81758481315158371, -0.82458930278502529, -0.83146961230254524,
    -0.83822470555483808, -0.84485356524970712, -0.8513551931052652,
    -0.85772861000027212, -0.8639728561215867, -0.87008699110871146,
    -0.8760700941954066, -0.881921264348355, -0.88763962040285393,
    -0.89322430119551532, -0.89867446569395382, -0.90398929312344334,
    -0.90916798309052238, -0.91420975570353069, -0.91911385169005777,
    -0.92387953251128674, -0.92850608047321559, -0.932992798834739,
    -0.937339011912575, -0.94154406518302081, -0.94560732538052128,
    -0.94952818059303667, -0.95330604035419386, -0.95694033573220882,
    -0.96043051941556579, -0.96377606579543984, -0.96697647104485207,
    -0.970031253194544, -0.97293995220556018, -0.97570213003852857,
    -0.97831737071962765, -0.98078528040323043, -0.98310548743121629,
    -0.98527764238894122, -0.98730141815785843, -0.989176509964781,
    -0.99090263542778, -0.99247953459871, -0.99390697000235606,
    -0.99518472667219693, -0.996312612182778, -0.99729045667869021,
    -0.99811811290014918, -0.99879545620517241, -0.99932238458834954,
    -0.99969881869620425, -0.9999247018391445, -1.0 };

  static const real_T tmp_0[257] = { 0.0, -0.012271538285719925,
    -0.024541228522912288, -0.036807222941358832, -0.049067674327418015,
    -0.061320736302208578, -0.073564563599667426, -0.0857973123444399,
    -0.0980171403295606, -0.11022220729388306, -0.1224106751992162,
    -0.13458070850712617, -0.14673047445536175, -0.15885814333386145,
    -0.17096188876030122, -0.18303988795514095, -0.19509032201612825,
    -0.20711137619221856, -0.2191012401568698, -0.23105810828067111,
    -0.24298017990326387, -0.25486565960451457, -0.26671275747489837,
    -0.27851968938505306, -0.29028467725446233, -0.30200594931922808,
    -0.31368174039889152, -0.32531029216226293, -0.33688985339222005,
    -0.34841868024943456, -0.35989503653498811, -0.37131719395183749,
    -0.38268343236508978, -0.3939920400610481, -0.40524131400498986,
    -0.41642956009763715, -0.42755509343028208, -0.43861623853852766,
    -0.44961132965460654, -0.46053871095824, -0.47139673682599764,
    -0.48218377207912272, -0.49289819222978404, -0.50353838372571758,
    -0.51410274419322166, -0.524589682678469, -0.53499761988709715,
    -0.54532498842204646, -0.55557023301960218, -0.56573181078361312,
    -0.57580819141784534, -0.58579785745643886, -0.59569930449243336,
    -0.60551104140432555, -0.61523159058062682, -0.62485948814238634,
    -0.63439328416364549, -0.64383154288979139, -0.65317284295377676,
    -0.66241577759017178, -0.67155895484701833, -0.680600997795453,
    -0.68954054473706683, -0.69837624940897292, -0.70710678118654757,
    -0.71573082528381859, -0.724247082951467, -0.73265427167241282,
    -0.74095112535495922, -0.74913639452345937, -0.75720884650648457,
    -0.765167265622459, -0.773010453362737, -0.78073722857209449,
    -0.78834642762660634, -0.79583690460888357, -0.80320753148064494,
    -0.81045719825259477, -0.81758481315158371, -0.82458930278502529,
    -0.83146961230254524, -0.83822470555483808, -0.84485356524970712,
    -0.8513551931052652, -0.85772861000027212, -0.8639728561215867,
    -0.87008699110871146, -0.8760700941954066, -0.881921264348355,
    -0.88763962040285393, -0.89322430119551532, -0.89867446569395382,
    -0.90398929312344334, -0.90916798309052238, -0.91420975570353069,
    -0.91911385169005777, -0.92387953251128674, -0.92850608047321559,
    -0.932992798834739, -0.937339011912575, -0.94154406518302081,
    -0.94560732538052128, -0.94952818059303667, -0.95330604035419386,
    -0.95694033573220882, -0.96043051941556579, -0.96377606579543984,
    -0.96697647104485207, -0.970031253194544, -0.97293995220556018,
    -0.97570213003852857, -0.97831737071962765, -0.98078528040323043,
    -0.98310548743121629, -0.98527764238894122, -0.98730141815785843,
    -0.989176509964781, -0.99090263542778, -0.99247953459871,
    -0.99390697000235606, -0.99518472667219693, -0.996312612182778,
    -0.99729045667869021, -0.99811811290014918, -0.99879545620517241,
    -0.99932238458834954, -0.99969881869620425, -0.9999247018391445, -1.0,
    -0.9999247018391445, -0.99969881869620425, -0.99932238458834954,
    -0.99879545620517241, -0.99811811290014918, -0.99729045667869021,
    -0.996312612182778, -0.99518472667219693, -0.99390697000235606,
    -0.99247953459871, -0.99090263542778, -0.989176509964781,
    -0.98730141815785843, -0.98527764238894122, -0.98310548743121629,
    -0.98078528040323043, -0.97831737071962765, -0.97570213003852857,
    -0.97293995220556018, -0.970031253194544, -0.96697647104485207,
    -0.96377606579543984, -0.96043051941556579, -0.95694033573220882,
    -0.95330604035419386, -0.94952818059303667, -0.94560732538052128,
    -0.94154406518302081, -0.937339011912575, -0.932992798834739,
    -0.92850608047321559, -0.92387953251128674, -0.91911385169005777,
    -0.91420975570353069, -0.90916798309052238, -0.90398929312344334,
    -0.89867446569395382, -0.89322430119551532, -0.88763962040285393,
    -0.881921264348355, -0.8760700941954066, -0.87008699110871146,
    -0.8639728561215867, -0.85772861000027212, -0.8513551931052652,
    -0.84485356524970712, -0.83822470555483808, -0.83146961230254524,
    -0.82458930278502529, -0.81758481315158371, -0.81045719825259477,
    -0.80320753148064494, -0.79583690460888357, -0.78834642762660634,
    -0.78073722857209449, -0.773010453362737, -0.765167265622459,
    -0.75720884650648457, -0.74913639452345937, -0.74095112535495922,
    -0.73265427167241282, -0.724247082951467, -0.71573082528381859,
    -0.70710678118654757, -0.69837624940897292, -0.68954054473706683,
    -0.680600997795453, -0.67155895484701833, -0.66241577759017178,
    -0.65317284295377676, -0.64383154288979139, -0.63439328416364549,
    -0.62485948814238634, -0.61523159058062682, -0.60551104140432555,
    -0.59569930449243336, -0.58579785745643886, -0.57580819141784534,
    -0.56573181078361312, -0.55557023301960218, -0.54532498842204646,
    -0.53499761988709715, -0.524589682678469, -0.51410274419322166,
    -0.50353838372571758, -0.49289819222978404, -0.48218377207912272,
    -0.47139673682599764, -0.46053871095824, -0.44961132965460654,
    -0.43861623853852766, -0.42755509343028208, -0.41642956009763715,
    -0.40524131400498986, -0.3939920400610481, -0.38268343236508978,
    -0.37131719395183749, -0.35989503653498811, -0.34841868024943456,
    -0.33688985339222005, -0.32531029216226293, -0.31368174039889152,
    -0.30200594931922808, -0.29028467725446233, -0.27851968938505306,
    -0.26671275747489837, -0.25486565960451457, -0.24298017990326387,
    -0.23105810828067111, -0.2191012401568698, -0.20711137619221856,
    -0.19509032201612825, -0.18303988795514095, -0.17096188876030122,
    -0.15885814333386145, -0.14673047445536175, -0.13458070850712617,
    -0.1224106751992162, -0.11022220729388306, -0.0980171403295606,
    -0.0857973123444399, -0.073564563599667426, -0.061320736302208578,
    -0.049067674327418015, -0.036807222941358832, -0.024541228522912288,
    -0.012271538285719925, -0.0 };

  real_T temp_re;
  real_T temp_im;
  real_T twid_re;
  real_T twid_im;
  int32_T temp_re_tmp;
  const real_T *costab;
  const real_T *sintab;
  costab = &tmp[0];
  sintab = &tmp_0[0];
  nRowsD4 = 128;
  ix = 0;
  ju = 0;
  iy = 0;
  for (iheight = 0; iheight < 511; iheight++) {
    xout[iy].re = xin[ix];
    xout[iy].im = 0.0;
    iy = 512;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  xout[iy].re = xin[ix];
  xout[iy].im = 0.0;
  for (ix = 0; ix <= 511; ix += 2) {
    temp_re = xout[ix + 1].re;
    temp_im = xout[ix + 1].im;
    xout[ix + 1].re = xout[ix].re - xout[ix + 1].re;
    xout[ix + 1].im = xout[ix].im - xout[ix + 1].im;
    xout[ix].re += temp_re;
    xout[ix].im += temp_im;
  }

  ix = 2;
  ju = 4;
  iheight = 509;
  while (nRowsD4 > 0) {
    for (iy = 0; iy < iheight; iy += ju) {
      temp_re_tmp = iy + ix;
      temp_re = xout[temp_re_tmp].re;
      temp_im = xout[iy + ix].im;
      xout[temp_re_tmp].re = xout[iy].re - xout[temp_re_tmp].re;
      xout[temp_re_tmp].im = xout[iy].im - temp_im;
      xout[iy].re += temp_re;
      xout[iy].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < 256; j += nRowsD4) {
      twid_re = costab[j];
      twid_im = sintab[j];
      iy = istart;
      ihi = istart + iheight;
      while (iy < ihi) {
        temp_re_tmp = iy + ix;
        temp_re = xout[temp_re_tmp].re * twid_re - xout[iy + ix].im * twid_im;
        temp_im = xout[iy + ix].im * twid_re + xout[iy + ix].re * twid_im;
        xout[temp_re_tmp].re = xout[iy].re - temp_re;
        xout[temp_re_tmp].im = xout[iy].im - temp_im;
        xout[iy].re += temp_re;
        xout[iy].im += temp_im;
        iy += ju;
      }

      istart++;
    }

    nRowsD4 /= 2;
    ix = ju;
    ju += ju;
    iheight -= ix;
  }
}

static void SpectrumEstimator_computeWindow(const
  dsp_SpectrumEstimator_filterb_T *obj, const real_T x[512], real_T Pxx[512])
{
  int32_T i;
  real_T tmp[512];
  filterbankFIR_bsxfun(x, obj->pWindowData, tmp);
  SpectrumEstimatorBase_computeFF(tmp, filterbankFIR_B.X);
  for (i = 0; i < 512; i++) {
    Pxx[i] = filterbankFIR_B.X[i].re * filterbankFIR_B.X[i].re -
      filterbankFIR_B.X[i].im * -filterbankFIR_B.X[i].im;
  }
}

static real_T filterbankFIR_mod(real_T x)
{
  real_T r;
  if ((!rtIsInf(x)) && (!rtIsNaN(x))) {
    if (x == 0.0) {
      r = 0.0;
    } else {
      r = fmod(x, 11.0);
      if (r == 0.0) {
        r = 0.0;
      } else {
        if (x < 0.0) {
          r += 11.0;
        }
      }
    }
  } else {
    r = (rtNaN);
  }

  return r;
}

static void SpectrumEstimatorBase_updatePer(dsp_SpectrumEstimator_filterb_T *obj,
  const real_T P[512])
{
  real_T u1;
  if (obj->pNumAvgsCounter + 1.0 < 10.0) {
    obj->pNumAvgsCounter++;
  } else {
    obj->pNumAvgsCounter = 10.0;
  }

  u1 = filterbankFIR_mod(obj->pNewPeriodogramIdx + 1.0);
  if ((1.0 > u1) || rtIsNaN(u1)) {
    obj->pNewPeriodogramIdx = 1.0;
  } else {
    obj->pNewPeriodogramIdx = u1;
  }

  memcpy(&obj->pPeriodogramMatrix[((int32_T)obj->pNewPeriodogramIdx << 9) + -512],
         &P[0], sizeof(real_T) << 9U);
}

static void SpectrumEstimatorBase_getPeriod(const
  dsp_SpectrumEstimator_filterb_T *obj, real_T P[512])
{
  int32_T xoffset;
  int32_T j;
  int32_T b_j;
  memcpy(&P[0], &obj->pPeriodogramMatrix[0], sizeof(real_T) << 9U);
  for (j = 0; j < 9; j++) {
    xoffset = (j + 1) << 9;
    for (b_j = 0; b_j < 512; b_j++) {
      P[b_j] += obj->pPeriodogramMatrix[xoffset + b_j];
    }
  }

  for (j = 0; j < 512; j++) {
    P[j] /= obj->pNumAvgsCounter;
  }
}

static void SpectrumEstimator_convertAndSca(const real_T P[512], real_T Pout[257])
{
  int32_T i;
  Pout[0] = P[0];
  for (i = 0; i < 255; i++) {
    Pout[i + 1] = P[i + 1] * 2.0;
  }

  Pout[256] = P[256];
}

static void SpectrumEstimatorBase_compute_n(const real_T xin[1024], creal_T
  xout[1024])
{
  int32_T yoff;
  int32_T ix;
  int32_T ju;
  int32_T iy;
  int32_T i;
  int32_T k;
  int32_T istart;
  int32_T j;
  int32_T ihi;
  int32_T chan;
  boolean_T tst;
  static const real_T tmp[257] = { 1.0, 0.9999247018391445, 0.99969881869620425,
    0.99932238458834954, 0.99879545620517241, 0.99811811290014918,
    0.99729045667869021, 0.996312612182778, 0.99518472667219693,
    0.99390697000235606, 0.99247953459871, 0.99090263542778, 0.989176509964781,
    0.98730141815785843, 0.98527764238894122, 0.98310548743121629,
    0.98078528040323043, 0.97831737071962765, 0.97570213003852857,
    0.97293995220556018, 0.970031253194544, 0.96697647104485207,
    0.96377606579543984, 0.96043051941556579, 0.95694033573220882,
    0.95330604035419386, 0.94952818059303667, 0.94560732538052128,
    0.94154406518302081, 0.937339011912575, 0.932992798834739,
    0.92850608047321559, 0.92387953251128674, 0.91911385169005777,
    0.91420975570353069, 0.90916798309052238, 0.90398929312344334,
    0.89867446569395382, 0.89322430119551532, 0.88763962040285393,
    0.881921264348355, 0.8760700941954066, 0.87008699110871146,
    0.8639728561215867, 0.85772861000027212, 0.8513551931052652,
    0.84485356524970712, 0.83822470555483808, 0.83146961230254524,
    0.82458930278502529, 0.81758481315158371, 0.81045719825259477,
    0.80320753148064494, 0.79583690460888357, 0.78834642762660634,
    0.78073722857209449, 0.773010453362737, 0.765167265622459,
    0.75720884650648457, 0.74913639452345937, 0.74095112535495922,
    0.73265427167241282, 0.724247082951467, 0.71573082528381859,
    0.70710678118654757, 0.69837624940897292, 0.68954054473706683,
    0.680600997795453, 0.67155895484701833, 0.66241577759017178,
    0.65317284295377676, 0.64383154288979139, 0.63439328416364549,
    0.62485948814238634, 0.61523159058062682, 0.60551104140432555,
    0.59569930449243336, 0.58579785745643886, 0.57580819141784534,
    0.56573181078361312, 0.55557023301960218, 0.54532498842204646,
    0.53499761988709715, 0.524589682678469, 0.51410274419322166,
    0.50353838372571758, 0.49289819222978404, 0.48218377207912272,
    0.47139673682599764, 0.46053871095824, 0.44961132965460654,
    0.43861623853852766, 0.42755509343028208, 0.41642956009763715,
    0.40524131400498986, 0.3939920400610481, 0.38268343236508978,
    0.37131719395183749, 0.35989503653498811, 0.34841868024943456,
    0.33688985339222005, 0.32531029216226293, 0.31368174039889152,
    0.30200594931922808, 0.29028467725446233, 0.27851968938505306,
    0.26671275747489837, 0.25486565960451457, 0.24298017990326387,
    0.23105810828067111, 0.2191012401568698, 0.20711137619221856,
    0.19509032201612825, 0.18303988795514095, 0.17096188876030122,
    0.15885814333386145, 0.14673047445536175, 0.13458070850712617,
    0.1224106751992162, 0.11022220729388306, 0.0980171403295606,
    0.0857973123444399, 0.073564563599667426, 0.061320736302208578,
    0.049067674327418015, 0.036807222941358832, 0.024541228522912288,
    0.012271538285719925, 0.0, -0.012271538285719925, -0.024541228522912288,
    -0.036807222941358832, -0.049067674327418015, -0.061320736302208578,
    -0.073564563599667426, -0.0857973123444399, -0.0980171403295606,
    -0.11022220729388306, -0.1224106751992162, -0.13458070850712617,
    -0.14673047445536175, -0.15885814333386145, -0.17096188876030122,
    -0.18303988795514095, -0.19509032201612825, -0.20711137619221856,
    -0.2191012401568698, -0.23105810828067111, -0.24298017990326387,
    -0.25486565960451457, -0.26671275747489837, -0.27851968938505306,
    -0.29028467725446233, -0.30200594931922808, -0.31368174039889152,
    -0.32531029216226293, -0.33688985339222005, -0.34841868024943456,
    -0.35989503653498811, -0.37131719395183749, -0.38268343236508978,
    -0.3939920400610481, -0.40524131400498986, -0.41642956009763715,
    -0.42755509343028208, -0.43861623853852766, -0.44961132965460654,
    -0.46053871095824, -0.47139673682599764, -0.48218377207912272,
    -0.49289819222978404, -0.50353838372571758, -0.51410274419322166,
    -0.524589682678469, -0.53499761988709715, -0.54532498842204646,
    -0.55557023301960218, -0.56573181078361312, -0.57580819141784534,
    -0.58579785745643886, -0.59569930449243336, -0.60551104140432555,
    -0.61523159058062682, -0.62485948814238634, -0.63439328416364549,
    -0.64383154288979139, -0.65317284295377676, -0.66241577759017178,
    -0.67155895484701833, -0.680600997795453, -0.68954054473706683,
    -0.69837624940897292, -0.70710678118654757, -0.71573082528381859,
    -0.724247082951467, -0.73265427167241282, -0.74095112535495922,
    -0.74913639452345937, -0.75720884650648457, -0.765167265622459,
    -0.773010453362737, -0.78073722857209449, -0.78834642762660634,
    -0.79583690460888357, -0.80320753148064494, -0.81045719825259477,
    -0.81758481315158371, -0.82458930278502529, -0.83146961230254524,
    -0.83822470555483808, -0.84485356524970712, -0.8513551931052652,
    -0.85772861000027212, -0.8639728561215867, -0.87008699110871146,
    -0.8760700941954066, -0.881921264348355, -0.88763962040285393,
    -0.89322430119551532, -0.89867446569395382, -0.90398929312344334,
    -0.90916798309052238, -0.91420975570353069, -0.91911385169005777,
    -0.92387953251128674, -0.92850608047321559, -0.932992798834739,
    -0.937339011912575, -0.94154406518302081, -0.94560732538052128,
    -0.94952818059303667, -0.95330604035419386, -0.95694033573220882,
    -0.96043051941556579, -0.96377606579543984, -0.96697647104485207,
    -0.970031253194544, -0.97293995220556018, -0.97570213003852857,
    -0.97831737071962765, -0.98078528040323043, -0.98310548743121629,
    -0.98527764238894122, -0.98730141815785843, -0.989176509964781,
    -0.99090263542778, -0.99247953459871, -0.99390697000235606,
    -0.99518472667219693, -0.996312612182778, -0.99729045667869021,
    -0.99811811290014918, -0.99879545620517241, -0.99932238458834954,
    -0.99969881869620425, -0.9999247018391445, -1.0 };

  static const real_T tmp_0[257] = { 0.0, -0.012271538285719925,
    -0.024541228522912288, -0.036807222941358832, -0.049067674327418015,
    -0.061320736302208578, -0.073564563599667426, -0.0857973123444399,
    -0.0980171403295606, -0.11022220729388306, -0.1224106751992162,
    -0.13458070850712617, -0.14673047445536175, -0.15885814333386145,
    -0.17096188876030122, -0.18303988795514095, -0.19509032201612825,
    -0.20711137619221856, -0.2191012401568698, -0.23105810828067111,
    -0.24298017990326387, -0.25486565960451457, -0.26671275747489837,
    -0.27851968938505306, -0.29028467725446233, -0.30200594931922808,
    -0.31368174039889152, -0.32531029216226293, -0.33688985339222005,
    -0.34841868024943456, -0.35989503653498811, -0.37131719395183749,
    -0.38268343236508978, -0.3939920400610481, -0.40524131400498986,
    -0.41642956009763715, -0.42755509343028208, -0.43861623853852766,
    -0.44961132965460654, -0.46053871095824, -0.47139673682599764,
    -0.48218377207912272, -0.49289819222978404, -0.50353838372571758,
    -0.51410274419322166, -0.524589682678469, -0.53499761988709715,
    -0.54532498842204646, -0.55557023301960218, -0.56573181078361312,
    -0.57580819141784534, -0.58579785745643886, -0.59569930449243336,
    -0.60551104140432555, -0.61523159058062682, -0.62485948814238634,
    -0.63439328416364549, -0.64383154288979139, -0.65317284295377676,
    -0.66241577759017178, -0.67155895484701833, -0.680600997795453,
    -0.68954054473706683, -0.69837624940897292, -0.70710678118654757,
    -0.71573082528381859, -0.724247082951467, -0.73265427167241282,
    -0.74095112535495922, -0.74913639452345937, -0.75720884650648457,
    -0.765167265622459, -0.773010453362737, -0.78073722857209449,
    -0.78834642762660634, -0.79583690460888357, -0.80320753148064494,
    -0.81045719825259477, -0.81758481315158371, -0.82458930278502529,
    -0.83146961230254524, -0.83822470555483808, -0.84485356524970712,
    -0.8513551931052652, -0.85772861000027212, -0.8639728561215867,
    -0.87008699110871146, -0.8760700941954066, -0.881921264348355,
    -0.88763962040285393, -0.89322430119551532, -0.89867446569395382,
    -0.90398929312344334, -0.90916798309052238, -0.91420975570353069,
    -0.91911385169005777, -0.92387953251128674, -0.92850608047321559,
    -0.932992798834739, -0.937339011912575, -0.94154406518302081,
    -0.94560732538052128, -0.94952818059303667, -0.95330604035419386,
    -0.95694033573220882, -0.96043051941556579, -0.96377606579543984,
    -0.96697647104485207, -0.970031253194544, -0.97293995220556018,
    -0.97570213003852857, -0.97831737071962765, -0.98078528040323043,
    -0.98310548743121629, -0.98527764238894122, -0.98730141815785843,
    -0.989176509964781, -0.99090263542778, -0.99247953459871,
    -0.99390697000235606, -0.99518472667219693, -0.996312612182778,
    -0.99729045667869021, -0.99811811290014918, -0.99879545620517241,
    -0.99932238458834954, -0.99969881869620425, -0.9999247018391445, -1.0,
    -0.9999247018391445, -0.99969881869620425, -0.99932238458834954,
    -0.99879545620517241, -0.99811811290014918, -0.99729045667869021,
    -0.996312612182778, -0.99518472667219693, -0.99390697000235606,
    -0.99247953459871, -0.99090263542778, -0.989176509964781,
    -0.98730141815785843, -0.98527764238894122, -0.98310548743121629,
    -0.98078528040323043, -0.97831737071962765, -0.97570213003852857,
    -0.97293995220556018, -0.970031253194544, -0.96697647104485207,
    -0.96377606579543984, -0.96043051941556579, -0.95694033573220882,
    -0.95330604035419386, -0.94952818059303667, -0.94560732538052128,
    -0.94154406518302081, -0.937339011912575, -0.932992798834739,
    -0.92850608047321559, -0.92387953251128674, -0.91911385169005777,
    -0.91420975570353069, -0.90916798309052238, -0.90398929312344334,
    -0.89867446569395382, -0.89322430119551532, -0.88763962040285393,
    -0.881921264348355, -0.8760700941954066, -0.87008699110871146,
    -0.8639728561215867, -0.85772861000027212, -0.8513551931052652,
    -0.84485356524970712, -0.83822470555483808, -0.83146961230254524,
    -0.82458930278502529, -0.81758481315158371, -0.81045719825259477,
    -0.80320753148064494, -0.79583690460888357, -0.78834642762660634,
    -0.78073722857209449, -0.773010453362737, -0.765167265622459,
    -0.75720884650648457, -0.74913639452345937, -0.74095112535495922,
    -0.73265427167241282, -0.724247082951467, -0.71573082528381859,
    -0.70710678118654757, -0.69837624940897292, -0.68954054473706683,
    -0.680600997795453, -0.67155895484701833, -0.66241577759017178,
    -0.65317284295377676, -0.64383154288979139, -0.63439328416364549,
    -0.62485948814238634, -0.61523159058062682, -0.60551104140432555,
    -0.59569930449243336, -0.58579785745643886, -0.57580819141784534,
    -0.56573181078361312, -0.55557023301960218, -0.54532498842204646,
    -0.53499761988709715, -0.524589682678469, -0.51410274419322166,
    -0.50353838372571758, -0.49289819222978404, -0.48218377207912272,
    -0.47139673682599764, -0.46053871095824, -0.44961132965460654,
    -0.43861623853852766, -0.42755509343028208, -0.41642956009763715,
    -0.40524131400498986, -0.3939920400610481, -0.38268343236508978,
    -0.37131719395183749, -0.35989503653498811, -0.34841868024943456,
    -0.33688985339222005, -0.32531029216226293, -0.31368174039889152,
    -0.30200594931922808, -0.29028467725446233, -0.27851968938505306,
    -0.26671275747489837, -0.25486565960451457, -0.24298017990326387,
    -0.23105810828067111, -0.2191012401568698, -0.20711137619221856,
    -0.19509032201612825, -0.18303988795514095, -0.17096188876030122,
    -0.15885814333386145, -0.14673047445536175, -0.13458070850712617,
    -0.1224106751992162, -0.11022220729388306, -0.0980171403295606,
    -0.0857973123444399, -0.073564563599667426, -0.061320736302208578,
    -0.049067674327418015, -0.036807222941358832, -0.024541228522912288,
    -0.012271538285719925, -0.0 };

  real_T temp_re;
  real_T temp_im;
  real_T twid_re;
  real_T twid_im;
  int32_T temp_re_tmp;
  const real_T *costab;
  const real_T *sintab;
  costab = &tmp[0];
  sintab = &tmp_0[0];
  for (chan = 0; chan < 2; chan++) {
    yoff = chan << 9;
    ix = chan << 9;
    ju = 0;
    iy = yoff;
    for (k = 0; k < 511; k++) {
      xout[iy].re = xin[ix];
      xout[iy].im = 0.0;
      iy = 512;
      tst = true;
      while (tst) {
        iy >>= 1;
        ju ^= iy;
        tst = ((ju & iy) == 0);
      }

      iy = yoff + ju;
      ix++;
    }

    xout[iy].re = xin[ix];
    xout[iy].im = 0.0;
    for (ix = yoff; ix <= yoff + 510; ix += 2) {
      temp_re = xout[ix + 1].re;
      temp_im = xout[ix + 1].im;
      xout[ix + 1].re = xout[ix].re - xout[ix + 1].re;
      xout[ix + 1].im = xout[ix].im - xout[ix + 1].im;
      xout[ix].re += temp_re;
      xout[ix].im += temp_im;
    }

    ix = 2;
    ju = 4;
    k = 128;
    iy = 509;
    while (k > 0) {
      i = yoff;
      ihi = yoff + iy;
      while (i < ihi) {
        temp_re_tmp = i + ix;
        temp_re = xout[temp_re_tmp].re;
        temp_im = xout[i + ix].im;
        xout[temp_re_tmp].re = xout[i].re - xout[temp_re_tmp].re;
        xout[temp_re_tmp].im = xout[i].im - temp_im;
        xout[i].re += temp_re;
        xout[i].im += temp_im;
        i += ju;
      }

      istart = yoff + 1;
      for (j = k; j < 256; j += k) {
        twid_re = costab[j];
        twid_im = sintab[j];
        i = istart;
        ihi = istart + iy;
        while (i < ihi) {
          temp_re_tmp = i + ix;
          temp_re = xout[temp_re_tmp].re * twid_re - xout[i + ix].im * twid_im;
          temp_im = xout[i + ix].im * twid_re + xout[i + ix].re * twid_im;
          xout[temp_re_tmp].re = xout[i].re - temp_re;
          xout[temp_re_tmp].im = xout[i].im - temp_im;
          xout[i].re += temp_re;
          xout[i].im += temp_im;
          i += ju;
        }

        istart++;
      }

      k /= 2;
      ix = ju;
      ju += ju;
      iy -= ix;
    }
  }
}

static void CrossSpectrumEstimator_computeW(const
  dsp_CrossSpectrumEstimator_fi_T *obj, const real_T x[512], const real_T y[512],
  creal_T Pxy[512])
{
  real_T tmp[512];
  real_T tmp_0[512];
  int32_T i;
  real_T Z_re;
  real_T Z_im;
  filterbankFIR_bsxfun(x, obj->pWindowData, tmp);
  filterbankFIR_bsxfun(y, obj->pWindowData, tmp_0);
  for (i = 0; i < 512; i++) {
    filterbankFIR_B.dv0[i] = tmp[i];
    filterbankFIR_B.dv0[i + 512] = tmp_0[i];
  }

  SpectrumEstimatorBase_compute_n(filterbankFIR_B.dv0, filterbankFIR_B.Z);
  for (i = 0; i < 512; i++) {
    Z_re = filterbankFIR_B.Z[512 + i].re;
    Z_im = -filterbankFIR_B.Z[512 + i].im;
    Pxy[i].re = filterbankFIR_B.Z[i].re * Z_re - filterbankFIR_B.Z[i].im * Z_im;
    Pxy[i].im = filterbankFIR_B.Z[i].re * Z_im + filterbankFIR_B.Z[i].im * Z_re;
  }
}

static void SpectrumEstimatorBase_updateP_n(dsp_CrossSpectrumEstimator_fi_T *obj,
  const creal_T P[512])
{
  real_T u1;
  if (obj->pNumAvgsCounter + 1.0 < 10.0) {
    obj->pNumAvgsCounter++;
  } else {
    obj->pNumAvgsCounter = 10.0;
  }

  u1 = filterbankFIR_mod(obj->pNewPeriodogramIdx + 1.0);
  if ((1.0 > u1) || rtIsNaN(u1)) {
    obj->pNewPeriodogramIdx = 1.0;
  } else {
    obj->pNewPeriodogramIdx = u1;
  }

  memcpy(&obj->pPeriodogramMatrix[((int32_T)obj->pNewPeriodogramIdx << 9) + -512],
         &P[0], sizeof(creal_T) << 9U);
}

static void SpectrumEstimatorBase_getPeri_n(const
  dsp_CrossSpectrumEstimator_fi_T *obj, creal_T P[512])
{
  int32_T xoffset;
  int32_T j;
  int32_T b_j;
  creal_T P_0;
  creal_T P_1;
  memcpy(&P[0], &obj->pPeriodogramMatrix[0], sizeof(creal_T) << 9U);
  for (j = 0; j < 9; j++) {
    xoffset = (j + 1) << 9;
    for (b_j = 0; b_j < 512; b_j++) {
      P_1.re = obj->pPeriodogramMatrix[xoffset + b_j].re + P[b_j].re;
      P_1.im = obj->pPeriodogramMatrix[xoffset + b_j].im + P[b_j].im;
      P[b_j] = P_1;
    }
  }

  for (j = 0; j < 512; j++) {
    if (P[j].im == 0.0) {
      P_0.re = P[j].re / obj->pNumAvgsCounter;
      P_0.im = 0.0;
    } else if (P[j].re == 0.0) {
      P_0.re = 0.0;
      P_0.im = P[j].im / obj->pNumAvgsCounter;
    } else {
      P_0.re = P[j].re / obj->pNumAvgsCounter;
      P_0.im = P[j].im / obj->pNumAvgsCounter;
    }

    P[j] = P_0;
  }
}

static void CrossSpectrumEstimator_convertA(const creal_T P[512], creal_T Pout
  [257])
{
  int32_T i;
  real_T ar;
  real_T ai;
  if (P[0].im == 0.0) {
    Pout[0].re = P[0].re;
    Pout[0].im = 0.0;
  } else if (P[0].re == 0.0) {
    Pout[0].re = 0.0;
    Pout[0].im = P[0].im;
  } else {
    Pout[0] = P[0];
  }

  for (i = 0; i < 255; i++) {
    ar = P[i + 1].re * 2.0;
    ai = P[i + 1].im * 2.0;
    if (ai == 0.0) {
      Pout[i + 1].re = ar;
      Pout[i + 1].im = 0.0;
    } else if (ar == 0.0) {
      Pout[i + 1].re = 0.0;
      Pout[i + 1].im = ai;
    } else {
      Pout[i + 1].re = ar;
      Pout[i + 1].im = ai;
    }
  }

  if (P[256].im == 0.0) {
    Pout[256].re = P[256].re;
    Pout[256].im = 0.0;
  } else if (P[256].re == 0.0) {
    Pout[256].re = 0.0;
    Pout[256].im = P[256].im;
  } else {
    Pout[256] = P[256];
  }
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

void MWDSPCG_FFT_Interleave_R2BR_D(const real_T x[], creal_T y[], int32_T nChans,
  int32_T nRows)
{
  int32_T br_j;
  int32_T yIdx;
  int32_T uIdx;
  int32_T j;
  int32_T nChansBy2;
  int32_T bit_fftLen;

  /* S-Function (sdspfft2): '<Root>/FFT1' */
  /* Bit-reverses the input data simultaneously with the interleaving operation,
     obviating the need for explicit data reordering later.  This requires an
     FFT with bit-reversed inputs.
   */
  br_j = 0;
  yIdx = 0;
  uIdx = 0;
  for (nChansBy2 = nChans >> 1; nChansBy2 != 0; nChansBy2--) {
    for (j = nRows; j - 1 > 0; j--) {
      y[yIdx + br_j].re = x[uIdx];
      y[yIdx + br_j].im = x[uIdx + nRows];
      uIdx++;

      /* Compute next bit-reversed destination index */
      bit_fftLen = nRows;
      do {
        bit_fftLen = (int32_T)((uint32_T)bit_fftLen >> 1);
        br_j ^= bit_fftLen;
      } while (!((br_j & bit_fftLen) != 0));
    }

    br_j += yIdx;
    y[br_j].re = x[uIdx];
    uIdx += nRows;
    y[br_j].im = x[uIdx];
    uIdx++;
    yIdx += nRows << 1;
    br_j = 0;
  }

  /* For an odd number of channels, prepare the last channel
     for a double-length real signal algorithm.  No actual
     interleaving is required, just a copy of the last column
     of real data, but now placed in bit-reversed order.
     We need to cast the real u pointer to a cDType_T pointer,
     in order to fake the interleaving, and cut the number
     of elements in half (half as many complex interleaved
     elements as compared to real non-interleaved elements).
   */
  if ((nChans & 1U) != 0U) {
    for (j = nRows >> 1; j - 1 > 0; j--) {
      y[yIdx + br_j].re = x[uIdx];
      y[yIdx + br_j].im = x[uIdx + 1];
      uIdx += 2;

      /* Compute next bit-reversed destination index */
      bit_fftLen = nRows >> 1;
      do {
        bit_fftLen = (int32_T)((uint32_T)bit_fftLen >> 1);
        br_j ^= bit_fftLen;
      } while (!((br_j & bit_fftLen) != 0));
    }

    br_j += yIdx;
    y[br_j].re = x[uIdx];
    y[br_j].im = x[uIdx + 1];
  }

  /* End of S-Function (sdspfft2): '<Root>/FFT1' */
}

void MWDSPCG_R2DIT_TBLS_Z(creal_T y[], int32_T nChans, int32_T nRows, int32_T
  fftLen, int32_T offset, const real_T tablePtr[], int32_T twiddleStep,
  boolean_T isInverse)
{
  int32_T nHalf;
  real_T twidRe;
  real_T twidIm;
  int32_T nQtr;
  int32_T fwdInvFactor;
  int32_T iCh;
  int32_T offsetCh;
  int32_T idelta;
  int32_T ix;
  int32_T k;
  int32_T kratio;
  int32_T istart;
  int32_T i1;
  int32_T j;
  int32_T i2;
  real_T tmp_re;
  real_T tmp_im;

  /* S-Function (sdspfft2): '<Root>/FFT1' */
  /* DSP System Toolbox Decimation in Time FFT  */
  /* Computation performed using table lookup  */
  /* Output type: complex real_T */
  nHalf = (fftLen >> 1) * twiddleStep;
  nQtr = nHalf >> 1;
  fwdInvFactor = isInverse ? -1 : 1;

  /* For each channel */
  offsetCh = offset;
  for (iCh = 0; iCh < nChans; iCh++) {
    /* Perform butterflies for the first stage, where no multiply is required. */
    for (ix = offsetCh; ix < (fftLen + offsetCh) - 1; ix += 2) {
      tmp_re = y[ix + 1].re;
      tmp_im = y[ix + 1].im;
      y[ix + 1].re = y[ix].re - tmp_re;
      y[ix + 1].im = y[ix].im - tmp_im;
      y[ix].re += tmp_re;
      y[ix].im += tmp_im;
    }

    idelta = 2;
    k = fftLen >> 2;
    kratio = k * twiddleStep;
    while (k > 0) {
      i1 = offsetCh;

      /* Perform the first butterfly in each remaining stage, where no multiply is required. */
      for (ix = 0; ix < k; ix++) {
        i2 = i1 + idelta;
        tmp_re = y[i2].re;
        tmp_im = y[i2].im;
        y[i2].re = y[i1].re - y[i2].re;
        y[i2].im = y[i1].im - tmp_im;
        y[i1].re += tmp_re;
        y[i1].im += tmp_im;
        i1 += idelta << 1;
      }

      istart = offsetCh;

      /* Perform remaining butterflies */
      for (j = kratio; j < nHalf; j += kratio) {
        i1 = istart + 1;
        twidRe = tablePtr[j];
        twidIm = tablePtr[j + nQtr] * (real_T)fwdInvFactor;
        for (ix = 0; ix < k; ix++) {
          i2 = i1 + idelta;
          tmp_re = y[i2].re * twidRe - y[i2].im * twidIm;
          tmp_im = y[i2].re * twidIm + y[i2].im * twidRe;
          y[i2].re = y[i1].re - tmp_re;
          y[i2].im = y[i1].im - tmp_im;
          y[i1].re += tmp_re;
          y[i1].im += tmp_im;
          i1 += idelta << 1;
        }

        istart++;
      }

      idelta <<= 1;
      k >>= 1;
      kratio >>= 1;
    }

    /* Point to next channel */
    offsetCh += nRows;
  }

  /* End of S-Function (sdspfft2): '<Root>/FFT1' */
}

void MWDSPCG_FFT_DblLen_Z(creal_T y[], int32_T nChans, int32_T nRows, const
  real_T twiddleTable[], int32_T twiddleStep)
{
  real_T accRe;
  real_T tempOut0Re;
  real_T tempOut0Im;
  int32_T N2;
  int32_T N4;
  int32_T W4;
  int32_T yIdx;
  int32_T ix;
  int32_T k;
  real_T cTemp_re;
  real_T cTemp_im;
  int32_T tempOut0Re_tmp;

  /* S-Function (sdspfft2): '<Root>/FFT1' */
  /* In-place "double-length" data recovery
     Table-based mem-optimized twiddle computation

     Used to recover linear-ordered length-N point complex FFT result
     from a linear-ordered complex length-N/2 point FFT, performed
     on N interleaved real values.
   */
  N2 = nRows >> 1;
  N4 = N2 >> 1;
  W4 = N4 * twiddleStep;
  yIdx = (nChans - 1) * nRows;
  if (nRows > 2) {
    tempOut0Re = y[N4 + yIdx].re;
    tempOut0Im = y[N4 + yIdx].im;
    y[(N2 + N4) + yIdx].re = tempOut0Re;
    y[(N2 + N4) + yIdx].im = tempOut0Im;
    y[N4 + yIdx].re = tempOut0Re;
    y[N4 + yIdx].im = -tempOut0Im;
  }

  if (nRows > 1) {
    y[N2 + yIdx].re = y[yIdx].re - y[yIdx].im;
    y[N2 + yIdx].im = 0.0;
  }

  y[yIdx].re += y[yIdx].im;
  y[yIdx].im = 0.0;
  k = twiddleStep;
  for (ix = 1; ix < N4; ix++) {
    tempOut0Re_tmp = ix + yIdx;
    tempOut0Re = (y[(N2 - ix) + yIdx].re + y[tempOut0Re_tmp].re) / 2.0;
    tempOut0Im = (y[ix + yIdx].im - y[(N2 - ix) + yIdx].im) / 2.0;
    accRe = (y[(N2 - ix) + yIdx].re - y[ix + yIdx].re) / 2.0;
    y[tempOut0Re_tmp].re = (y[(N2 - ix) + yIdx].im + y[ix + yIdx].im) / 2.0;
    y[tempOut0Re_tmp].im = accRe;
    cTemp_re = y[ix + yIdx].re * twiddleTable[k] - -twiddleTable[W4 - k] * y[ix
      + yIdx].im;
    cTemp_im = y[ix + yIdx].im * twiddleTable[k] + -twiddleTable[W4 - k] * y[ix
      + yIdx].re;
    accRe = cTemp_im;
    y[tempOut0Re_tmp].re = tempOut0Re + cTemp_re;
    y[tempOut0Re_tmp].im = tempOut0Im + cTemp_im;
    cTemp_im = -y[ix + yIdx].im;
    tempOut0Re_tmp = (nRows - ix) + yIdx;
    y[tempOut0Re_tmp].re = y[ix + yIdx].re;
    y[tempOut0Re_tmp].im = cTemp_im;
    y[(N2 + ix) + yIdx].re = tempOut0Re - cTemp_re;
    y[(N2 + ix) + yIdx].im = tempOut0Im - accRe;
    tempOut0Re = y[(N2 + ix) + yIdx].im;
    y[(N2 - ix) + yIdx].re = y[(N2 + ix) + yIdx].re;
    y[(N2 - ix) + yIdx].im = -tempOut0Re;
    k += twiddleStep;
  }

  /* End of S-Function (sdspfft2): '<Root>/FFT1' */
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

static void SystemProp_matlabCodege_nzusynb(dsp_simulink_TransferFunction_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void fi_SystemCore_release_nzusynbuf(dsp_simulink_TransferFunction_T *obj)
{
  if ((obj->isInitialized == 1) && obj->isSetupComplete) {
    if (obj->pSpectrumEstimator.isInitialized == 1) {
      obj->pSpectrumEstimator.isInitialized = 2;
    }

    if (obj->pCrossSpectrumEstimator.isInitialized == 1) {
      obj->pCrossSpectrumEstimator.isInitialized = 2;
    }
  }
}

static void filte_SystemCore_delete_nzusynb(dsp_simulink_TransferFunction_T *obj)
{
  fi_SystemCore_release_nzusynbuf(obj);
}

static void matlabCodegenHandle_mat_nzusynb(dsp_simulink_TransferFunction_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodege_nzusynb(obj, true);
    filte_SystemCore_delete_nzusynb(obj);
  }
}

static void SystemProp_matlabCodeg_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void f_SystemCore_release_nzusynbufg(dsp_SpectrumEstimator_filterb_T *obj)
{
  if (obj->isInitialized == 1) {
    obj->isInitialized = 2;
  }
}

static void filt_SystemCore_delete_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj)
{
  f_SystemCore_release_nzusynbufg(obj);
}

static void matlabCodegenHandle_ma_nzusynbu(dsp_SpectrumEstimator_filterb_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodeg_nzusynbu(obj, true);
    filt_SystemCore_delete_nzusynbu(obj);
  }
}

static void SystemProp_matlabCodegen_nzusyn(dsp_PhaseExtractor_filterbank_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void filt_SystemCore_release_nzusynb(dsp_PhaseExtractor_filterbank_T *obj)
{
  if ((obj->isInitialized == 1) && obj->isSetupComplete &&
      (obj->pPhaseDifferentiator.isInitialized == 1)) {
    obj->pPhaseDifferentiator.isInitialized = 2;
    if (obj->pPhaseDifferentiator.isSetupComplete) {
      obj->pPhaseDifferentiator.pNumChans = -1;
    }
  }
}

static void filter_SystemCore_delete_nzusyn(dsp_PhaseExtractor_filterbank_T *obj)
{
  filt_SystemCore_release_nzusynb(obj);
}

static void matlabCodegenHandle_matl_nzusyn(dsp_PhaseExtractor_filterbank_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegen_nzusyn(obj, true);
    filter_SystemCore_delete_nzusyn(obj);
  }
}

static void SystemProp_matlabCodegenS_nzusy(dsp_private_PhaseDifferentiat_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void fil_SystemCore_release_nzusynbu(dsp_private_PhaseDifferentiat_T *obj)
{
  if (obj->isInitialized == 1) {
    obj->isInitialized = 2;
    if (obj->isSetupComplete) {
      obj->pNumChans = -1;
    }
  }
}

static void filterb_SystemCore_delete_nzusy(dsp_private_PhaseDifferentiat_T *obj)
{
  fil_SystemCore_release_nzusynbu(obj);
}

static void matlabCodegenHandle_matla_nzusy(dsp_private_PhaseDifferentiat_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenS_nzusy(obj, true);
    filterb_SystemCore_delete_nzusy(obj);
  }
}

static void SystemProp_matlabCodegenSe_nzus(dsp_Differentiator_filterbank_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void filterb_SystemCore_release_nzus(dspcodegen_FIRFilter_filterba_T *obj)
{
  if (obj->isInitialized == 1) {
    obj->isInitialized = 2;
  }
}

static void filt_Differentiator_releaseImpl(dsp_Differentiator_filterbank_T *obj)
{
  filterb_SystemCore_release_nzus(obj->FilterObj);
  obj->NumChannels = -1;
}

static void SystemCore_releaseWrapper_nzus(dsp_Differentiator_filterbank_T *obj)
{
  if (obj->isSetupComplete) {
    filt_Differentiator_releaseImpl(obj);
  }
}

static void filter_SystemCore_release_nzusy(dsp_Differentiator_filterbank_T *obj)
{
  if (obj->isInitialized == 1) {
    SystemCore_releaseWrapper_nzus(obj);
  }
}

static void filterba_SystemCore_delete_nzus(dsp_Differentiator_filterbank_T *obj)
{
  filter_SystemCore_release_nzusy(obj);
}

static void matlabCodegenHandle_matlab_nzus(dsp_Differentiator_filterbank_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenSe_nzus(obj, true);
    filterba_SystemCore_delete_nzus(obj);
  }
}

static void SystemProp_matlabCodegenSet_nzu(dspcodegen_FIRFilter_filterba_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void filterban_SystemCore_delete_nzu(dspcodegen_FIRFilter_filterba_T *obj)
{
  filterb_SystemCore_release_nzus(obj);
}

static void matlabCodegenHandle_matlabC_nzu(dspcodegen_FIRFilter_filterba_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenSet_nzu(obj, true);
    filterban_SystemCore_delete_nzu(obj);
  }
}

static void SystemProp_matlabCodegenSetA_nz(dsp_simulink_VariableBandwidt_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void filterbank_SystemCore_release_n(dspcodegen_FIRFilter_filter_n_T *obj)
{
  if (obj->isInitialized == 1) {
    obj->isInitialized = 2;
  }
}

static void VariableBandwidthFilterBase_rel(dsp_simulink_VariableBandwidt_T *obj)
{
  filterbank_SystemCore_release_n(&obj->pfilter);
  obj->pNumChannels = -1;
}

static void filte_SystemCore_releaseWrapper(dsp_simulink_VariableBandwidt_T *obj)
{
  if (obj->isSetupComplete) {
    VariableBandwidthFilterBase_rel(obj);
  }
}

static void filterban_SystemCore_release_nz(dsp_simulink_VariableBandwidt_T *obj)
{
  if (obj->isInitialized == 1) {
    filte_SystemCore_releaseWrapper(obj);
  }
}

static void filterbank_SystemCore_delete_nz(dsp_simulink_VariableBandwidt_T *obj)
{
  filterban_SystemCore_release_nz(obj);
}

static void matlabCodegenHandle_matlabCo_nz(dsp_simulink_VariableBandwidt_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenSetA_nz(obj, true);
    filterbank_SystemCore_delete_nz(obj);
  }
}

static void SystemProp_matlabCodegenSetAn_n(dspcodegen_FIRFilter_filter_n_T *obj,
  boolean_T value)
{
  obj->matlabCodegenIsDeleted = value;
}

static void filterbankF_SystemCore_delete_n(dspcodegen_FIRFilter_filter_n_T *obj)
{
  filterbank_SystemCore_release_n(obj);
}

static void matlabCodegenHandle_matlabCod_n(dspcodegen_FIRFilter_filter_n_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    SystemProp_matlabCodegenSetAn_n(obj, true);
    filterbankF_SystemCore_delete_n(obj);
  }
}

/* Model step function */
void filterbankFIR_step(void)
{
  /* local block i/o variables */
  real_T rtb_OriginalSignal[512];
  real_T rtb_Gain1_f[512];
  real_T rtb_Gain2_g[512];
  real_T rtb_Gain3_m[512];
  real_T rtb_Gain4_e[512];
  real_T rtb_Gain5_p[512];
  real_T rtb_Gain6_m[512];
  real_T rtb_Gain7_j[512];
  real_T rtb_Gain8_a[512];
  real_T rtb_Sum2[512];
  int32_T i;
  int32_T j;
  real_T rtb_DifferentiatorFilter_0[257];
  real_T tmp[257];
  dsp_FIRFilter_0_filterbankF_n_T *obj;
  boolean_T p;
  boolean_T p_0;
  dsp_simulink_VariableBandwidt_T *obj_0;
  real_T newCutoff;
  real_T pCoefficients[31];
  dspcodegen_FIRFilter_filter_n_T *obj_1;
  real_T val[31];
  dsp_FIRFilter_0_filterbankFIR_T *obj_2;
  real_T acc2;
  real_T pxx[257];
  real_T rtb_Gain1[512];
  real_T rtb_GeneratedFilterBlock[512];
  real_T rtb_Gain2[512];
  real_T rtb_Gain3[512];
  real_T rtb_Gain4[512];
  real_T rtb_Gain5[512];
  real_T rtb_Gain6[512];
  real_T rtb_Gain7[512];
  real_T rtb_Gain8[512];
  real_T rtb_Gain9[512];
  real_T rtb_Gain10[512];
  real_T rtb_Gain11[512];
  real_T rtb_Gain12[512];
  real_T rtb_Gain13[512];
  real_T rtb_Gain14[512];
  real_T rtb_Gain15[512];
  static const real_T tmp_0[31] = { 0.0, 0.010926199633097156,
    0.043227271178699567, 0.095491502812526274, 0.16543469682057088,
    0.24999999999999994, 0.34549150281252627, 0.44773576836617329,
    0.55226423163382665, 0.65450849718747373, 0.74999999999999989,
    0.8345653031794289, 0.90450849718747373, 0.95677272882130049,
    0.98907380036690284, 1.0, 0.98907380036690284, 0.95677272882130049,
    0.90450849718747373, 0.8345653031794289, 0.74999999999999989,
    0.65450849718747373, 0.55226423163382665, 0.44773576836617329,
    0.34549150281252627, 0.24999999999999994, 0.16543469682057088,
    0.095491502812526274, 0.043227271178699567, 0.010926199633097156, 0.0 };

  static const real_T tmp_1[512] = { 0.0, 3.7649080427748505E-5,
    0.00015059065189787502, 0.00033880770582522812, 0.0006022718974137975,
    0.00094094354992541041, 0.0013547716606548965, 0.0018436939086109994,
    0.0024076366639015356, 0.0030465149988219697, 0.0037602327006450165,
    0.0045486822861099951, 0.0054117450176094928, 0.0063492909210707826,
    0.0073611788055293892, 0.0084472562843918575, 0.0096073597983847847,
    0.010841314640186173, 0.012148934980735715, 0.013530023897219912,
    0.014984373402728013, 0.016511764477573965, 0.01811196710228008,
    0.019784740292217107, 0.021529832133895588, 0.023346979822903069,
    0.025235909703481663, 0.02719633730973936, 0.029227967408489597,
    0.03133049404371252, 0.033503600582630522, 0.035746959763392205,
    0.038060233744356631, 0.040443074154971115, 0.042895122148234655,
    0.04541600845473881, 0.048005353438278331, 0.050662767153023092,
    0.053387849402242338, 0.056180189798573033, 0.059039367825822475,
    0.0619649529022967, 0.06495650444564427, 0.068013571939206652,
    0.071135694999863941, 0.0743224034473674, 0.077573217375146442,
    0.080887647222580961, 0.084265193848727382, 0.087705348607487355,
    0.091207593424208144, 0.094771400873702616, 0.098396234259677529,
    0.10208154769555822, 0.10582678618669683, 0.10963138571395276,
    0.1134947733186315, 0.11741636718877052, 0.12139557674675772,
    0.12543180273827031, 0.12952443732252039, 0.13367286416379359,
    0.1378764585242665, 0.1421345873580907, 0.14644660940672621,
    0.15081187529551354, 0.15522972763146653, 0.15969950110227343,
    0.16422052257649078, 0.16879211120491411, 0.17341357852311157,
    0.17808422855510425, 0.18280335791817726, 0.18757025592880677,
    0.19238420470968659, 0.19724447929783723, 0.20215034775378327,
    0.20710107127178057, 0.21209590429107733, 0.21713409460819338,
    0.22221488349019886, 0.22733750578897677, 0.23250119005645137,
    0.23770515866076558, 0.24294862790338917, 0.24823080813714121,
    0.25355090388510793, 0.25890811396043856, 0.2643016315870011,
    0.26973064452088, 0.2751943351726967, 0.28069188073073614,
    0.2862224532848589, 0.29178521995118134, 0.29737934299750507,
    0.30300397996947592, 0.30865828381745508, 0.3143414030240812,
    0.32005248173250589, 0.32579065987528277, 0.33155507330389,
    0.33734485391886848, 0.34315912980055419, 0.3489970253403859,
    0.35485766137276886, 0.3607401553074735, 0.36664362126255079,
    0.37256717019774266, 0.378509910048368, 0.38447094585966435,
    0.39044937992156514, 0.39644431190389073, 0.40245483899193585,
    0.40848005602242948, 0.41451905561984931, 0.4205709283330693,
    0.42663476277231915, 0.43270964574643689, 0.43879466240039189,
    0.44488889635305839, 0.45099142983521961, 0.45710134382778006,
    0.46321771820016627, 0.46933963184889566, 0.47546616283629095,
    0.48159638852932052, 0.48772938573854385, 0.49386423085714004,
    0.49999999999999994, 0.50613576914286, 0.512270614261456,
    0.51840361147067948, 0.524533837163709, 0.53066036815110429,
    0.53678228179983367, 0.54289865617221988, 0.54900857016478033,
    0.55511110364694149, 0.56120533759960811, 0.567290354253563,
    0.57336523722768085, 0.57942907166693058, 0.58548094438015064,
    0.59151994397757046, 0.5975451610080641, 0.60355568809610927,
    0.60955062007843486, 0.61552905414033554, 0.62149008995163191,
    0.62743282980225723, 0.63335637873744921, 0.6392598446925265,
    0.645142338627231, 0.651002974659614, 0.6568408701994457,
    0.66265514608113141, 0.66844492669611, 0.67420934012471723,
    0.67994751826749411, 0.6856585969759188, 0.69134171618254481,
    0.696996020030524, 0.70262065700249487, 0.70821478004881855,
    0.71377754671514093, 0.71930811926926363, 0.72480566482730335,
    0.73026935547912009, 0.73569836841299885, 0.74109188603956133,
    0.746449096114892, 0.75176919186285873, 0.75705137209661078,
    0.76229484133923431, 0.76749880994354847, 0.77266249421102318,
    0.777785116509801, 0.78286590539180656, 0.78790409570892272,
    0.79289892872821943, 0.79784965224621662, 0.80275552070216272,
    0.80761579529031335, 0.81242974407119317, 0.81719664208182263,
    0.82191577144489569, 0.82658642147688832, 0.831207888795086,
    0.83577947742350922, 0.84030049889772651, 0.84477027236853353,
    0.84918812470448635, 0.85355339059327373, 0.8578654126419093,
    0.86212354147573333, 0.86632713583620635, 0.87047556267747939,
    0.87456819726172963, 0.87860442325324239, 0.88258363281122953,
    0.88650522668136844, 0.89036861428604719, 0.89417321381330317,
    0.89791845230444167, 0.90160376574032242, 0.90522859912629738,
    0.90879240657579174, 0.91229465139251253, 0.91573480615127267,
    0.919112352777419, 0.92242678262485356, 0.9256775965526326,
    0.928864305000136, 0.93198642806079335, 0.93504349555435562,
    0.93803504709770325, 0.94096063217417747, 0.94381981020142691,
    0.94661215059775761, 0.94933723284697691, 0.95199464656172172,
    0.95458399154526119, 0.95710487785176535, 0.95955692584502894,
    0.96193976625564337, 0.96425304023660774, 0.96649639941736942,
    0.96866950595628742, 0.97077203259151035, 0.97280366269026053,
    0.97476409029651834, 0.97665302017709688, 0.97847016786610441,
    0.98021525970778289, 0.98188803289772, 0.983488235522426, 0.985015626597272,
    0.98646997610278, 0.98785106501926423, 0.98915868535981377,
    0.99039264020161522, 0.9915527437156082, 0.99263882119447056,
    0.99365070907892916, 0.99458825498239056, 0.99545131771389,
    0.996239767299355, 0.99695348500117809, 0.99759236333609835,
    0.998156306091389, 0.9986452283393451, 0.99905905645007453,
    0.9993977281025862, 0.99966119229417472, 0.99984940934810207,
    0.99996235091957231, 1.0, 0.99996235091957231, 0.99984940934810207,
    0.99966119229417472, 0.9993977281025862, 0.99905905645007453,
    0.9986452283393451, 0.998156306091389, 0.99759236333609835,
    0.99695348500117809, 0.996239767299355, 0.99545131771389,
    0.99458825498239056, 0.99365070907892916, 0.99263882119447056,
    0.9915527437156082, 0.99039264020161522, 0.98915868535981377,
    0.98785106501926423, 0.98646997610278, 0.985015626597272, 0.983488235522426,
    0.98188803289772, 0.98021525970778289, 0.97847016786610441,
    0.97665302017709688, 0.97476409029651834, 0.97280366269026053,
    0.97077203259151035, 0.96866950595628742, 0.96649639941736942,
    0.96425304023660774, 0.96193976625564337, 0.95955692584502894,
    0.95710487785176535, 0.95458399154526119, 0.95199464656172172,
    0.94933723284697691, 0.94661215059775761, 0.94381981020142691,
    0.94096063217417747, 0.93803504709770325, 0.93504349555435562,
    0.93198642806079335, 0.928864305000136, 0.9256775965526326,
    0.92242678262485356, 0.919112352777419, 0.91573480615127267,
    0.91229465139251253, 0.90879240657579174, 0.90522859912629738,
    0.90160376574032242, 0.89791845230444167, 0.89417321381330317,
    0.89036861428604719, 0.88650522668136844, 0.88258363281122953,
    0.87860442325324239, 0.87456819726172963, 0.87047556267747939,
    0.86632713583620635, 0.86212354147573333, 0.8578654126419093,
    0.85355339059327373, 0.84918812470448635, 0.84477027236853353,
    0.84030049889772651, 0.83577947742350922, 0.831207888795086,
    0.82658642147688832, 0.82191577144489569, 0.81719664208182263,
    0.81242974407119317, 0.80761579529031335, 0.80275552070216272,
    0.79784965224621662, 0.79289892872821943, 0.78790409570892272,
    0.78286590539180656, 0.777785116509801, 0.77266249421102318,
    0.76749880994354847, 0.76229484133923431, 0.75705137209661078,
    0.75176919186285873, 0.746449096114892, 0.74109188603956133,
    0.73569836841299885, 0.73026935547912009, 0.72480566482730335,
    0.71930811926926363, 0.71377754671514093, 0.70821478004881855,
    0.70262065700249487, 0.696996020030524, 0.69134171618254481,
    0.6856585969759188, 0.67994751826749411, 0.67420934012471723,
    0.66844492669611, 0.66265514608113141, 0.6568408701994457, 0.651002974659614,
    0.645142338627231, 0.6392598446925265, 0.63335637873744921,
    0.62743282980225723, 0.62149008995163191, 0.61552905414033554,
    0.60955062007843486, 0.60355568809610927, 0.5975451610080641,
    0.59151994397757046, 0.58548094438015064, 0.57942907166693058,
    0.57336523722768085, 0.567290354253563, 0.56120533759960811,
    0.55511110364694149, 0.54900857016478033, 0.54289865617221988,
    0.53678228179983367, 0.53066036815110429, 0.524533837163709,
    0.51840361147067948, 0.512270614261456, 0.50613576914286,
    0.49999999999999994, 0.49386423085714004, 0.48772938573854385,
    0.48159638852932052, 0.47546616283629095, 0.46933963184889566,
    0.46321771820016627, 0.45710134382778006, 0.45099142983521961,
    0.44488889635305839, 0.43879466240039189, 0.43270964574643689,
    0.42663476277231915, 0.4205709283330693, 0.41451905561984931,
    0.40848005602242948, 0.40245483899193585, 0.39644431190389073,
    0.39044937992156514, 0.38447094585966435, 0.378509910048368,
    0.37256717019774266, 0.36664362126255079, 0.3607401553074735,
    0.35485766137276886, 0.3489970253403859, 0.34315912980055419,
    0.33734485391886848, 0.33155507330389, 0.32579065987528277,
    0.32005248173250589, 0.3143414030240812, 0.30865828381745508,
    0.30300397996947592, 0.29737934299750507, 0.29178521995118134,
    0.2862224532848589, 0.28069188073073614, 0.2751943351726967,
    0.26973064452088, 0.2643016315870011, 0.25890811396043856,
    0.25355090388510793, 0.24823080813714121, 0.24294862790338917,
    0.23770515866076558, 0.23250119005645137, 0.22733750578897677,
    0.22221488349019886, 0.21713409460819338, 0.21209590429107733,
    0.20710107127178057, 0.20215034775378327, 0.19724447929783723,
    0.19238420470968659, 0.18757025592880677, 0.18280335791817726,
    0.17808422855510425, 0.17341357852311157, 0.16879211120491411,
    0.16422052257649078, 0.15969950110227343, 0.15522972763146653,
    0.15081187529551354, 0.14644660940672621, 0.1421345873580907,
    0.1378764585242665, 0.13367286416379359, 0.12952443732252039,
    0.12543180273827031, 0.12139557674675772, 0.11741636718877052,
    0.1134947733186315, 0.10963138571395276, 0.10582678618669683,
    0.10208154769555822, 0.098396234259677529, 0.094771400873702616,
    0.091207593424208144, 0.087705348607487355, 0.084265193848727382,
    0.080887647222580961, 0.077573217375146442, 0.0743224034473674,
    0.071135694999863941, 0.068013571939206652, 0.06495650444564427,
    0.0619649529022967, 0.059039367825822475, 0.056180189798573033,
    0.053387849402242338, 0.050662767153023092, 0.048005353438278331,
    0.04541600845473881, 0.042895122148234655, 0.040443074154971115,
    0.038060233744356631, 0.035746959763392205, 0.033503600582630522,
    0.03133049404371252, 0.029227967408489597, 0.02719633730973936,
    0.025235909703481663, 0.023346979822903069, 0.021529832133895588,
    0.019784740292217107, 0.01811196710228008, 0.016511764477573965,
    0.014984373402728013, 0.013530023897219912, 0.012148934980735715,
    0.010841314640186173, 0.0096073597983847847, 0.0084472562843918575,
    0.0073611788055293892, 0.0063492909210707826, 0.0054117450176094928,
    0.0045486822861099951, 0.0037602327006450165, 0.0030465149988219697,
    0.0024076366639015356, 0.0018436939086109994, 0.0013547716606548965,
    0.00094094354992541041, 0.0006022718974137975, 0.00033880770582522812,
    0.00015059065189787502, 3.7649080427748505E-5 };

  creal_T pyx;

  /* S-Function (sdsprandsrc2): '<Root>/Random Source1' */
  RandSrc_GZ_D(filterbankFIR_B.RandomSource1,
               &filterbankFIR_P.RandomSource1_MeanVal, 1,
               &filterbankFIR_P.RandomSource1_VarianceRTP, 1,
               filterbankFIR_DW.RandomSource1_STATE_DWORK, 1, 512);

  /* DiscreteFir: '<S32>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coefficien[j];
    }

    for (j = 0; j < 702 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coefficien[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 189; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states[i];
  }

  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S26>/Gain1' */
    rtb_Gain1[i] = filterbankFIR_P.gainConsts[0] * rtb_GeneratedFilterBlock[i];

    /* DiscreteFir: '<S33>/Generated Filter Block' */
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_b[j];
    }

    for (j = 0; j < 1583 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_b[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_f[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* DiscreteFir: '<S33>/Generated Filter Block' */
  /* Update delay line for next frame */
  for (i = 1070; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_f[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_f[i];
  }

  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_f[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S26>/Gain2' */
    rtb_Gain2[i] = filterbankFIR_P.gainConsts[1] * rtb_GeneratedFilterBlock[i];

    /* DiscreteFir: '<S34>/Generated Filter Block' */
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_f[j];
    }

    for (j = 0; j < 1262 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_f[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_c[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* DiscreteFir: '<S34>/Generated Filter Block' */
  /* Update delay line for next frame */
  for (i = 749; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_c[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_c[i];
  }

  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_c[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S26>/Gain3' */
    rtb_Gain3[i] = filterbankFIR_P.gainConsts[2] * rtb_GeneratedFilterBlock[i];

    /* DiscreteFir: '<S36>/Generated Filter Block' */
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_d[j];
    }

    for (j = 0; j < 989 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_d[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_l[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* DiscreteFir: '<S36>/Generated Filter Block' */
  /* Update delay line for next frame */
  for (i = 476; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_l[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_l[i];
  }

  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_l[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S27>/Gain4' */
    rtb_Gain4[i] = filterbankFIR_P.gainConsts[3] * rtb_GeneratedFilterBlock[i];

    /* DiscreteFir: '<S35>/Generated Filter Block' */
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_p[j];
    }

    for (j = 0; j < 785 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_p[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_j[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* DiscreteFir: '<S35>/Generated Filter Block' */
  /* Update delay line for next frame */
  for (i = 272; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_j[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_j[i];
  }

  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_j[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S27>/Gain5' */
    rtb_Gain5[i] = filterbankFIR_P.gainConsts[4] * rtb_GeneratedFilterBlock[i];

    /* DiscreteFir: '<S37>/Generated Filter Block' */
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_n[j];
    }

    for (j = 0; j < 627 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_n[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_e[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* DiscreteFir: '<S37>/Generated Filter Block' */
  /* Update delay line for next frame */
  for (i = 114; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_e[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_e[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_e[511 - i] =
      filterbankFIR_B.RandomSource1[i];

    /* Gain: '<S27>/Gain6' */
    rtb_Gain6[i] = filterbankFIR_P.gainConsts[5] * rtb_GeneratedFilterBlock[i];
  }

  /* DiscreteFir: '<S40>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 498; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_l[j];
    }

    for (j = 0; j < 498 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_l[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_eh[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 498; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 499; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_l[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 498; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_eh[497 - i] =
      filterbankFIR_B.RandomSource1[i + 14];
  }

  /* Gain: '<S28>/Gain7' */
  for (i = 0; i < 512; i++) {
    rtb_Gain7[i] = filterbankFIR_P.gainConsts[6] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S28>/Gain7' */

  /* DiscreteFir: '<S38>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 394; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_p1[j];
    }

    for (j = 0; j < 394 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_p1[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_jw[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 394; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 395; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_p1[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 394; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_jw[393 - i] =
      filterbankFIR_B.RandomSource1[i + 118];
  }

  /* Gain: '<S28>/Gain8' */
  for (i = 0; i < 512; i++) {
    rtb_Gain8[i] = filterbankFIR_P.gainConsts[7] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S28>/Gain8' */

  /* DiscreteFir: '<S39>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 313; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_a[j];
    }

    for (j = 0; j < 313 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_a[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_i[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 313; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 314; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_a[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 313; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_i[312 - i] =
      filterbankFIR_B.RandomSource1[i + 199];
  }

  /* Gain: '<S28>/Gain9' */
  for (i = 0; i < 512; i++) {
    rtb_Gain9[i] = filterbankFIR_P.gainConsts[8] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S28>/Gain9' */

  /* DiscreteFir: '<S42>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 249; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_pa[j];
    }

    for (j = 0; j < 249 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_pa[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_fb[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 249; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 250; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_pa[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 249; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_fb[248 - i] =
      filterbankFIR_B.RandomSource1[i + 263];
  }

  /* Gain: '<S29>/Gain10' */
  for (i = 0; i < 512; i++) {
    rtb_Gain10[i] = filterbankFIR_P.gainConsts[9] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S29>/Gain10' */

  /* DiscreteFir: '<S41>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 198; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_e[j];
    }

    for (j = 0; j < 198 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_e[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_d[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 198; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 199; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_e[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 198; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_d[197 - i] =
      filterbankFIR_B.RandomSource1[i + 314];
  }

  /* Gain: '<S29>/Gain11' */
  for (i = 0; i < 512; i++) {
    rtb_Gain11[i] = filterbankFIR_P.gainConsts[10] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S29>/Gain11' */

  /* DiscreteFir: '<S43>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 157; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_a4[j];
    }

    for (j = 0; j < 157 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_a4[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_g[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 157; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 158; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_a4[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 157; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_g[156 - i] =
      filterbankFIR_B.RandomSource1[i + 355];
  }

  /* Gain: '<S29>/Gain12' */
  for (i = 0; i < 512; i++) {
    rtb_Gain12[i] = filterbankFIR_P.gainConsts[11] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S29>/Gain12' */

  /* DiscreteFir: '<S46>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 153; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_k[j];
    }

    for (j = 0; j < 153 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_k[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_b[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 153; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 154; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_k[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 153; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_b[152 - i] =
      filterbankFIR_B.RandomSource1[i + 359];
  }

  /* Gain: '<S30>/Gain13 ' */
  for (i = 0; i < 512; i++) {
    rtb_Gain13[i] = filterbankFIR_P.gainConsts[12] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S30>/Gain13 ' */

  /* DiscreteFir: '<S44>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 122; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_j[j];
    }

    for (j = 0; j < 122 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_j[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_fi[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 122; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 123; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_j[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 122; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_fi[121 - i] =
      filterbankFIR_B.RandomSource1[i + 390];
  }

  /* Gain: '<S30>/Gain14' */
  for (i = 0; i < 512; i++) {
    rtb_Gain14[i] = filterbankFIR_P.gainConsts[13] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S30>/Gain14' */

  /* DiscreteFir: '<S45>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 97; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_bz[j];
    }

    for (j = 0; j < 97 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_bz[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_o[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 97; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 98; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_bz[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 97; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_o[96 - i] =
      filterbankFIR_B.RandomSource1[i + 415];
  }

  /* Gain: '<S30>/Gain15' */
  for (i = 0; i < 512; i++) {
    rtb_Gain15[i] = filterbankFIR_P.gainConsts[14] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S30>/Gain15' */

  /* DiscreteFir: '<S47>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 77; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_o[j];
    }

    for (j = 0; j < 77 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_o[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_bf[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 77; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 78; j++) {
      newCutoff += filterbankFIR_B.RandomSource1[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_o[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 77; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_bf[76 - i] =
      filterbankFIR_B.RandomSource1[i + 435];
  }

  /* Sum: '<S2>/Sum2' incorporates:
   *  Gain: '<S31>/Gain16'
   *  Sum: '<S26>/Add'
   *  Sum: '<S27>/Add'
   *  Sum: '<S28>/Add'
   *  Sum: '<S29>/Add'
   *  Sum: '<S30>/Add'
   */
  for (i = 0; i < 512; i++) {
    filterbankFIR_B.Sum2[i] = ((((((rtb_Gain1[i] + rtb_Gain2[i]) + rtb_Gain3[i])
      + ((rtb_Gain4[i] + rtb_Gain5[i]) + rtb_Gain6[i])) + ((rtb_Gain7[i] +
      rtb_Gain8[i]) + rtb_Gain9[i])) + ((rtb_Gain10[i] + rtb_Gain11[i]) +
      rtb_Gain12[i])) + ((rtb_Gain13[i] + rtb_Gain14[i]) + rtb_Gain15[i])) +
      filterbankFIR_P.gainConsts[15] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Sum: '<S2>/Sum2' */

  /* Outputs for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
  /* MATLABSystem: '<S48>/TFE' incorporates:
   *  Buffer: '<S48>/Buffer_x'
   *  Buffer: '<S48>/Buffer_y'
   */
  if (filterbankFIR_DW.obj.pFrameCounter >= filterbankFIR_DW.obj.pFrameDelay) {
    if (filterbankFIR_DW.obj.pSpectrumEstimator.isInitialized != 1) {
      filterbankFIR_DW.obj.pSpectrumEstimator.isSetupComplete = false;
      filterbankFIR_DW.obj.pSpectrumEstimator.isInitialized = 1;
      newCutoff = 0.0;
      for (i = 0; i < 512; i++) {
        filterbankFIR_DW.obj.pSpectrumEstimator.pWindowData[i] = tmp_1[i];
        newCutoff += filterbankFIR_DW.obj.pSpectrumEstimator.pWindowData[i] *
          filterbankFIR_DW.obj.pSpectrumEstimator.pWindowData[i];
      }

      filterbankFIR_DW.obj.pSpectrumEstimator.pWindowPower = newCutoff;
      for (i = 0; i < 257; i++) {
        filterbankFIR_DW.obj.pSpectrumEstimator.pW[i] = 0.001953125 * (real_T)i;
      }

      filterbankFIR_DW.obj.pSpectrumEstimator.pReferenceLoad =
        filterbankFIR_DW.obj.pSpectrumEstimator.ReferenceLoad;
      filterbankFIR_DW.obj.pSpectrumEstimator.isSetupComplete = true;
      filterbankFIR_DW.obj.pSpectrumEstimator.TunablePropsChanged = false;
      filterbankFIR_DW.obj.pSpectrumEstimator.pNumAvgsCounter = 0.0;
      filterbankFIR_DW.obj.pSpectrumEstimator.pNewPeriodogramIdx = 0.0;
      memset(&filterbankFIR_DW.obj.pSpectrumEstimator.pPeriodogramMatrix[0], 0,
             5120U * sizeof(real_T));
    }

    if (filterbankFIR_DW.obj.pSpectrumEstimator.TunablePropsChanged) {
      filterbankFIR_DW.obj.pSpectrumEstimator.TunablePropsChanged = false;
      filterbankFIR_DW.obj.pSpectrumEstimator.pReferenceLoad =
        filterbankFIR_DW.obj.pSpectrumEstimator.ReferenceLoad;
    }

    SpectrumEstimator_computeWindow(&filterbankFIR_DW.obj.pSpectrumEstimator,
      filterbankFIR_B.RandomSource1, rtb_Gain1);
    for (i = 0; i < 512; i++) {
      rtb_GeneratedFilterBlock[i] = rtb_Gain1[i] /
        filterbankFIR_DW.obj.pSpectrumEstimator.pWindowPower;
    }

    SpectrumEstimatorBase_updatePer(&filterbankFIR_DW.obj.pSpectrumEstimator,
      rtb_GeneratedFilterBlock);
    SpectrumEstimatorBase_getPeriod(&filterbankFIR_DW.obj.pSpectrumEstimator,
      rtb_Gain1);
    SpectrumEstimator_convertAndSca(rtb_Gain1, pxx);
    for (i = 0; i < 257; i++) {
      pxx[i] /= filterbankFIR_DW.obj.pSpectrumEstimator.pReferenceLoad;
    }

    if (filterbankFIR_DW.obj.pCrossSpectrumEstimator.isInitialized != 1) {
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.isSetupComplete = false;
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.isInitialized = 1;
      newCutoff = 0.0;
      for (i = 0; i < 512; i++) {
        filterbankFIR_DW.obj.pCrossSpectrumEstimator.pWindowData[i] = tmp_1[i];
        newCutoff += filterbankFIR_DW.obj.pCrossSpectrumEstimator.pWindowData[i]
          * filterbankFIR_DW.obj.pCrossSpectrumEstimator.pWindowData[i];
      }

      filterbankFIR_DW.obj.pCrossSpectrumEstimator.pWindowPower = newCutoff;
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.isSetupComplete = true;
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.pNumAvgsCounter = 0.0;
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.pNewPeriodogramIdx = 0.0;
      memset(&filterbankFIR_DW.obj.pCrossSpectrumEstimator.pPeriodogramMatrix[0],
             0, 5120U * sizeof(creal_T));
    }

    CrossSpectrumEstimator_computeW
      (&filterbankFIR_DW.obj.pCrossSpectrumEstimator, filterbankFIR_B.Sum2,
       filterbankFIR_B.RandomSource1, filterbankFIR_B.b_P);
    newCutoff = filterbankFIR_DW.obj.pCrossSpectrumEstimator.pWindowPower;
    for (i = 0; i < 512; i++) {
      if (filterbankFIR_B.b_P[i].im == 0.0) {
        filterbankFIR_B.b_P_m[i].re = filterbankFIR_B.b_P[i].re / newCutoff;
        filterbankFIR_B.b_P_m[i].im = 0.0;
      } else if (filterbankFIR_B.b_P[i].re == 0.0) {
        filterbankFIR_B.b_P_m[i].re = 0.0;
        filterbankFIR_B.b_P_m[i].im = filterbankFIR_B.b_P[i].im / newCutoff;
      } else {
        filterbankFIR_B.b_P_m[i].re = filterbankFIR_B.b_P[i].re / newCutoff;
        filterbankFIR_B.b_P_m[i].im = filterbankFIR_B.b_P[i].im / newCutoff;
      }
    }

    SpectrumEstimatorBase_updateP_n
      (&filterbankFIR_DW.obj.pCrossSpectrumEstimator, filterbankFIR_B.b_P_m);
    SpectrumEstimatorBase_getPeri_n
      (&filterbankFIR_DW.obj.pCrossSpectrumEstimator, filterbankFIR_B.b_P);
    CrossSpectrumEstimator_convertA(filterbankFIR_B.b_P, filterbankFIR_B.pyx);
    for (i = 0; i < 257; i++) {
      if (filterbankFIR_B.pyx[i].im == 0.0) {
        pyx.re = filterbankFIR_B.pyx[i].re / pxx[i];
        pyx.im = 0.0;
      } else if (filterbankFIR_B.pyx[i].re == 0.0) {
        pyx.re = 0.0;
        pyx.im = filterbankFIR_B.pyx[i].im / pxx[i];
      } else {
        pyx.re = filterbankFIR_B.pyx[i].re / pxx[i];
        pyx.im = filterbankFIR_B.pyx[i].im / pxx[i];
      }

      filterbankFIR_B.pyx[i] = pyx;
    }
  } else {
    filterbankFIR_DW.obj.pFrameCounter++;
    memset(&filterbankFIR_B.pyx[0], 0, 257U * sizeof(creal_T));
  }

  /* End of Outputs for SubSystem: '<S3>/Discrete Transfer Function Estimator' */

  /* MATLABSystem: '<S3>/Phase Extractor' incorporates:
   *  MATLABSystem: '<S48>/TFE'
   */
  if (filterbankFIR_DW.obj_g.pPhaseDifferentiator.isInitialized != 1) {
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.isSetupComplete = false;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.isInitialized = 1;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pNumChans = 1;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.re = 1.0;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.im = 0.0;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.isSetupComplete = true;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.re = 1.0;
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.im = 0.0;
  }

  filterbankFIR_B.a[0].re =
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.re;
  filterbankFIR_B.a[0].im =
    -filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.im;
  for (i = 0; i < 256; i++) {
    /* Outputs for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
    filterbankFIR_B.a[i + 1].re = filterbankFIR_B.pyx[(i + 1) - 1].re;
    filterbankFIR_B.a[i + 1].im = -filterbankFIR_B.pyx[(i + 1) - 1].im;

    /* End of Outputs for SubSystem: '<S3>/Discrete Transfer Function Estimator' */
  }

  for (i = 0; i < 257; i++) {
    /* Outputs for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
    pxx[i] = rt_atan2d_snf(filterbankFIR_B.pyx[i].re * filterbankFIR_B.a[i].im +
      filterbankFIR_B.pyx[i].im * filterbankFIR_B.a[i].re, filterbankFIR_B.pyx[i]
      .re * filterbankFIR_B.a[i].re - filterbankFIR_B.pyx[i].im *
      filterbankFIR_B.a[i].im);

    /* End of Outputs for SubSystem: '<S3>/Discrete Transfer Function Estimator' */
  }

  for (i = 0; i < 256; i++) {
    pxx[i + 1] += pxx[i];
  }

  /* Gain: '<S3>/Gain1' incorporates:
   *  MATLABSystem: '<S3>/Phase Extractor'
   */
  for (i = 0; i < 257; i++) {
    pxx[i] *= filterbankFIR_P.Gain1_Gain;
  }

  /* End of Gain: '<S3>/Gain1' */

  /* MATLABSystem: '<S3>/Differentiator Filter' */
  if (filterbankFIR_DW.obj_m.FilterObj->isInitialized != 1) {
    filterbankFIR_DW.obj_m.FilterObj->isSetupComplete = false;
    filterbankFIR_DW.obj_m.FilterObj->isInitialized = 1;
    filterbankFIR_DW.obj_m.FilterObj->isSetupComplete = true;

    /* System object Initialization function: dsp.FIRFilter */
    filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[0] =
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
    filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[1] =
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
    filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[2] =
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
  }

  obj = &filterbankFIR_DW.obj_m.FilterObj->cSFunObject;

  /* System object Outputs function: dsp.FIRFilter */
  /* Consume delay line and beginning of input samples */
  acc2 = pxx[0] * filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P1_Coefficients
    [0];
  newCutoff = acc2;
  for (j = 0; j < 3; j++) {
    acc2 = obj->P1_Coefficients[j + 1] * obj->W0_states[j];
    newCutoff += acc2;
  }

  rtb_DifferentiatorFilter_0[0] = newCutoff;
  newCutoff = 0.0;
  for (j = 0; j < 2; j++) {
    acc2 = pxx[1 - j] * obj->P1_Coefficients[j];
    newCutoff += acc2;
  }

  for (j = 0; j < 2; j++) {
    acc2 = obj->P1_Coefficients[j + 2] * obj->W0_states[j];
    newCutoff += acc2;
  }

  rtb_DifferentiatorFilter_0[1] = newCutoff;
  newCutoff = 0.0;
  for (j = 0; j < 3; j++) {
    acc2 = pxx[2 - j] * obj->P1_Coefficients[j];
    newCutoff += acc2;
  }

  acc2 = filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[0] *
    filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P1_Coefficients[3];
  rtb_DifferentiatorFilter_0[2] = newCutoff + acc2;

  /* DiscreteFir: '<S25>/Generated Filter Block' incorporates:
   *  DiscreteFir: '<S10>/Generated Filter Block'
   *  DiscreteFir: '<S11>/Generated Filter Block'
   *  DiscreteFir: '<S12>/Generated Filter Block'
   *  DiscreteFir: '<S13>/Generated Filter Block'
   *  DiscreteFir: '<S14>/Generated Filter Block'
   *  DiscreteFir: '<S15>/Generated Filter Block'
   *  DiscreteFir: '<S16>/Generated Filter Block'
   *  DiscreteFir: '<S17>/Generated Filter Block'
   *  DiscreteFir: '<S18>/Generated Filter Block'
   *  DiscreteFir: '<S19>/Generated Filter Block'
   *  DiscreteFir: '<S20>/Generated Filter Block'
   *  DiscreteFir: '<S21>/Generated Filter Block'
   *  DiscreteFir: '<S22>/Generated Filter Block'
   *  DiscreteFir: '<S23>/Generated Filter Block'
   *  DiscreteFir: '<S24>/Generated Filter Block'
   *  DiscreteFir: '<S32>/Generated Filter Block'
   *  DiscreteFir: '<S33>/Generated Filter Block'
   *  DiscreteFir: '<S34>/Generated Filter Block'
   *  DiscreteFir: '<S35>/Generated Filter Block'
   *  DiscreteFir: '<S36>/Generated Filter Block'
   *  DiscreteFir: '<S37>/Generated Filter Block'
   *  DiscreteFir: '<S38>/Generated Filter Block'
   *  DiscreteFir: '<S39>/Generated Filter Block'
   *  DiscreteFir: '<S40>/Generated Filter Block'
   *  DiscreteFir: '<S41>/Generated Filter Block'
   *  DiscreteFir: '<S42>/Generated Filter Block'
   *  DiscreteFir: '<S43>/Generated Filter Block'
   *  DiscreteFir: '<S44>/Generated Filter Block'
   *  DiscreteFir: '<S45>/Generated Filter Block'
   *  DiscreteFir: '<S46>/Generated Filter Block'
   *  DiscreteFir: '<S47>/Generated Filter Block'
   */
  /* MATLABSystem: '<S3>/Differentiator Filter' */
  /* Consume remaining input samples */
  for (i = 3; i < 257; i++) {
    newCutoff = pxx[i] * obj->P1_Coefficients[0];
    acc2 = pxx[i - 1] * obj->P1_Coefficients[1];
    newCutoff += acc2;
    acc2 = pxx[i - 2] * obj->P1_Coefficients[2];
    newCutoff += acc2;
    acc2 = pxx[i - 3] * obj->P1_Coefficients[3];
    rtb_DifferentiatorFilter_0[i] = newCutoff + acc2;
  }

  /* Update delay line for next frame */
  filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[2] = pxx[254];
  filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[1] = pxx[255];
  filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[0] = pxx[256];

  /* Start for MATLABSystem: '<Root>/Lowpass Filter' */
  newCutoff = filterbankFIR_DW.obj_k.CutoffFrequency;

  /* MATLABSystem: '<Root>/Lowpass Filter' incorporates:
   *  MATLABSystem: '<S3>/Differentiator Filter'
   */
  p = false;
  p_0 = true;
  if (!(newCutoff == filterbankFIR_P.LowpassFilter_fc)) {
    p_0 = false;
  }

  if (p_0) {
    p = true;
  }

  if (!p) {
    p = (filterbankFIR_DW.obj_k.isInitialized == 1);
    if (p) {
      filterbankFIR_DW.obj_k.TunablePropsChanged = true;
    }

    filterbankFIR_DW.obj_k.CutoffFrequency = filterbankFIR_P.LowpassFilter_fc;
  }

  obj_0 = &filterbankFIR_DW.obj_k;
  if (filterbankFIR_DW.obj_k.TunablePropsChanged) {
    filterbankFIR_DW.obj_k.TunablePropsChanged = false;
    newCutoff = filterbankFIR_DW.obj_k.CutoffFrequency * 2.0;
    newCutoff /= filterbankFIR_SampleRate;
    for (i = 0; i < 31; i++) {
      pCoefficients[i] = 0.0;
      if (-15 + i == 0) {
        pCoefficients[i] = newCutoff;
      } else {
        pCoefficients[i] = sin((-15.0 + (real_T)i) * newCutoff *
          3.1415926535897931) / ((-15.0 + (real_T)i) * 3.1415926535897931);
      }
    }

    obj_1 = &filterbankFIR_DW.obj_k.pfilter;
    for (i = 0; i < 31; i++) {
      val[i] = tmp_0[i] * pCoefficients[i];
    }

    for (i = 0; i < 31; i++) {
      obj_1->cSFunObject.P1_Coefficients[i] = val[i];
    }

    for (i = 0; i < 31; i++) {
      obj_0->pfilter.Numerator[i] = val[i];
    }

    memcpy(&val[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0], 31U * sizeof
           (real_T));
    if (filterbankFIR_sum(val) > 0.0) {
      memcpy(&pCoefficients[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0],
             31U * sizeof(real_T));
      memcpy(&val[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0], 31U * sizeof
             (real_T));
      newCutoff = filterbankFIR_sum(val);
      for (i = 0; i < 31; i++) {
        pCoefficients[i] /= newCutoff;
      }

      obj_1 = &filterbankFIR_DW.obj_k.pfilter;
      for (i = 0; i < 31; i++) {
        obj_1->cSFunObject.P1_Coefficients[i] = pCoefficients[i];
      }

      for (i = 0; i < 31; i++) {
        obj_0->pfilter.Numerator[i] = pCoefficients[i];
      }
    }
  }

  obj_1 = &filterbankFIR_DW.obj_k.pfilter;
  if (obj_0->pfilter.isInitialized != 1) {
    obj_0->pfilter.isSetupComplete = false;
    obj_0->pfilter.isInitialized = 1;
    obj_0->pfilter.isSetupComplete = true;

    /* System object Initialization function: dsp.FIRFilter */
    for (i = 0; i < 30; i++) {
      obj_1->cSFunObject.W0_states[i] = obj_1->cSFunObject.P0_InitialStates;
    }
  }

  obj_2 = &obj_0->pfilter.cSFunObject;

  /* System object Outputs function: dsp.FIRFilter */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 30; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      acc2 = rtb_DifferentiatorFilter_0[i - j] * obj_2->P1_Coefficients[j];
      newCutoff += acc2;
    }

    for (j = 0; j < 30 - i; j++) {
      acc2 = obj_2->P1_Coefficients[(i + j) + 1] * obj_2->W0_states[j];
      newCutoff += acc2;
    }

    tmp[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 30; i < 257; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 31; j++) {
      acc2 = rtb_DifferentiatorFilter_0[i - j] * obj_2->P1_Coefficients[j];
      newCutoff += acc2;
    }

    tmp[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 30; i++) {
    obj_2->W0_states[29 - i] = rtb_DifferentiatorFilter_0[i + 227];
  }

  memcpy(&filterbankFIR_B.LowpassFilter[0], &tmp[0], 257U * sizeof(real_T));

  /* S-Function (sdspsubmtrx): '<Root>/Submatrix' */
  memcpy(&filterbankFIR_B.Submatrix[0], &filterbankFIR_B.LowpassFilter[30], 227U
         * sizeof(real_T));

  /* S-Function (sdspimpgen2): '<Root>/Discrete  Impulse' */
  /* DSP System Toolbox Discrete Impulse Generator (sdspimpgen2) - '<Root>/Discrete  Impulse' */
  {
    real_T *yPtr = &rtb_OriginalSignal[0];
    int_T sampCount = 512;
    while (sampCount--) {
      if (filterbankFIR_DW.DiscreteImpulse_COUNT > 0) {
        if (filterbankFIR_DW.DiscreteImpulse_COUNT == 1) {
          *yPtr++ = (1.0);
        } else {
          *yPtr++ = (0.0);
        }

        filterbankFIR_DW.DiscreteImpulse_COUNT--;
      } else {
        *yPtr++ = (0.0);
      }
    }
  }

  /* DiscreteFir: '<S10>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_m[j];
    }

    for (j = 0; j < 989 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_m[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_n[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 476; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_n[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_n[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_n[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S4>/Gain1' */
    rtb_Gain1_f[i] = filterbankFIR_P.gainConsts[0] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S4>/Compressor' */
  filterbankFIR_Compressor2(rtb_Gain1_f, &filterbankFIR_B.Compressor,
    &filterbankFIR_DW.Compressor, &filterbankFIR_P.Compressor);

  /* DiscreteFir: '<S11>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_bs[j];
    }

    for (j = 0; j < 1583 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_bs[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_cf[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 1070; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_cf[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_cf[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_cf[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S4>/Gain2' */
    rtb_Gain2_g[i] = filterbankFIR_P.gainConsts[1] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S4>/Compressor1' */
  filterbankFIR_Compressor2(rtb_Gain2_g, &filterbankFIR_B.Compressor1,
    &filterbankFIR_DW.Compressor1, &filterbankFIR_P.Compressor1);

  /* DiscreteFir: '<S12>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_i[j];
    }

    for (j = 0; j < 1262 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_i[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_a[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 749; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_a[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_a[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_a[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S4>/Gain3' */
    rtb_Gain3_m[i] = filterbankFIR_P.gainConsts[2] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S4>/Compressor2' */
  filterbankFIR_Compressor2(rtb_Gain3_m, &filterbankFIR_B.Compressor2,
    &filterbankFIR_DW.Compressor2, &filterbankFIR_P.Compressor2);

  /* DiscreteFir: '<S14>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_oh[j];
    }

    for (j = 0; j < 989 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_oh[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_bs[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 476; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_bs[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_bs[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_bs[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S5>/Gain4' */
    rtb_Gain4_e[i] = filterbankFIR_P.gainConsts[3] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S5>/Compressor' */
  filterbankFIR_Compressor2(rtb_Gain4_e, &filterbankFIR_B.Compressor_p,
    &filterbankFIR_DW.Compressor_p, &filterbankFIR_P.Compressor_p);

  /* DiscreteFir: '<S13>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_h[j];
    }

    for (j = 0; j < 785 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_h[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_et[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 272; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_et[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_et[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_et[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S5>/Gain5' */
    rtb_Gain5_p[i] = filterbankFIR_P.gainConsts[4] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S5>/Compressor1' */
  filterbankFIR_Compressor2(rtb_Gain5_p, &filterbankFIR_B.Compressor1_p,
    &filterbankFIR_DW.Compressor1_p, &filterbankFIR_P.Compressor1_p);

  /* DiscreteFir: '<S15>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_n0[j];
    }

    for (j = 0; j < 627 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_n0[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_lf[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 114; i >= 0; i--) {
    filterbankFIR_DW.GeneratedFilterBlock_states_lf[512 + i] =
      filterbankFIR_DW.GeneratedFilterBlock_states_lf[i];
  }

  for (i = 0; i < 512; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_lf[511 - i] =
      rtb_OriginalSignal[i];

    /* Gain: '<S5>/Gain6' */
    rtb_Gain6_m[i] = filterbankFIR_P.gainConsts[5] * rtb_GeneratedFilterBlock[i];
  }

  /* MATLABSystem: '<S5>/Compressor2' */
  filterbankFIR_Compressor2(rtb_Gain6_m, &filterbankFIR_B.Compressor2_p,
    &filterbankFIR_DW.Compressor2_p, &filterbankFIR_P.Compressor2_p);

  /* DiscreteFir: '<S18>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 498; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_c[j];
    }

    for (j = 0; j < 498 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_c[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_ed[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 498; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 499; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_c[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 498; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_ed[497 - i] =
      rtb_OriginalSignal[i + 14];
  }

  /* Gain: '<S6>/Gain7' */
  for (i = 0; i < 512; i++) {
    rtb_Gain7_j[i] = filterbankFIR_P.gainConsts[6] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S6>/Gain7' */

  /* MATLABSystem: '<S6>/Compressor' */
  filterbankFIR_Compressor2(rtb_Gain7_j, &filterbankFIR_B.Compressor_pn,
    &filterbankFIR_DW.Compressor_pn, &filterbankFIR_P.Compressor_pn);

  /* DiscreteFir: '<S16>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 394; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_g[j];
    }

    for (j = 0; j < 394 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffici_g[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_nq[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 394; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 395; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffici_g[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 394; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_nq[393 - i] =
      rtb_OriginalSignal[i + 118];
  }

  /* Gain: '<S6>/Gain8' */
  for (i = 0; i < 512; i++) {
    rtb_Gain8_a[i] = filterbankFIR_P.gainConsts[7] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S6>/Gain8' */

  /* MATLABSystem: '<S6>/Compressor1' */
  filterbankFIR_Compressor2(rtb_Gain8_a, &filterbankFIR_B.Compressor1_pn,
    &filterbankFIR_DW.Compressor1_pn, &filterbankFIR_P.Compressor1_pn);

  /* DiscreteFir: '<S17>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 313; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_nt[j];
    }

    for (j = 0; j < 313 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_nt[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_ec[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 313; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 314; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_nt[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 313; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_ec[312 - i] =
      rtb_OriginalSignal[i + 199];
  }

  /* Gain: '<S6>/Gain9' */
  for (i = 0; i < 512; i++) {
    rtb_Gain1[i] = filterbankFIR_P.gainConsts[8] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S6>/Gain9' */

  /* DiscreteFir: '<S20>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 249; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_f2[j];
    }

    for (j = 0; j < 249 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_f2[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_i2[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 249; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 250; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_f2[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 249; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_i2[248 - i] =
      rtb_OriginalSignal[i + 263];
  }

  /* Gain: '<S7>/Gain10' */
  for (i = 0; i < 512; i++) {
    rtb_Gain2[i] = filterbankFIR_P.gainConsts[9] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S7>/Gain10' */

  /* DiscreteFir: '<S19>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 198; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_li[j];
    }

    for (j = 0; j < 198 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_li[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_ou[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 198; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 199; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_li[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 198; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_ou[197 - i] =
      rtb_OriginalSignal[i + 314];
  }

  /* Gain: '<S7>/Gain11' */
  for (i = 0; i < 512; i++) {
    rtb_Gain3[i] = filterbankFIR_P.gainConsts[10] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S7>/Gain11' */

  /* DiscreteFir: '<S21>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 157; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_by[j];
    }

    for (j = 0; j < 157 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_by[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_fu[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 157; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 158; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_by[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 157; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_fu[156 - i] =
      rtb_OriginalSignal[i + 355];
  }

  /* Gain: '<S7>/Gain12' */
  for (i = 0; i < 512; i++) {
    rtb_Gain4[i] = filterbankFIR_P.gainConsts[11] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S7>/Gain12' */

  /* DiscreteFir: '<S24>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 153; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_bn[j];
    }

    for (j = 0; j < 153 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_bn[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_iu[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 153; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 154; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_bn[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 153; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_iu[152 - i] =
      rtb_OriginalSignal[i + 359];
  }

  /* Gain: '<S8>/Gain13 ' */
  for (i = 0; i < 512; i++) {
    rtb_Gain5[i] = filterbankFIR_P.gainConsts[12] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S8>/Gain13 ' */

  /* DiscreteFir: '<S22>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 122; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_f1[j];
    }

    for (j = 0; j < 122 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_f1[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_ck[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 122; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 123; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_f1[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 122; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_ck[121 - i] =
      rtb_OriginalSignal[i + 390];
  }

  /* Gain: '<S8>/Gain14' */
  for (i = 0; i < 512; i++) {
    rtb_Gain6[i] = filterbankFIR_P.gainConsts[13] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S8>/Gain14' */

  /* DiscreteFir: '<S23>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 97; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_km[j];
    }

    for (j = 0; j < 97 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_km[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_bt[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 97; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 98; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_km[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 97; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_bt[96 - i] =
      rtb_OriginalSignal[i + 415];
  }

  /* Gain: '<S8>/Gain15' */
  for (i = 0; i < 512; i++) {
    rtb_Gain7[i] = filterbankFIR_P.gainConsts[14] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Gain: '<S8>/Gain15' */

  /* DiscreteFir: '<S25>/Generated Filter Block' */
  /* Consume delay line and beginning of input samples */
  for (i = 0; i < 77; i++) {
    newCutoff = 0.0;
    for (j = 0; j < i + 1; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_ns[j];
    }

    for (j = 0; j < 77 - i; j++) {
      newCutoff += filterbankFIR_P.GeneratedFilterBlock_Coeffic_ns[(i + j) + 1] *
        filterbankFIR_DW.GeneratedFilterBlock_states_av[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Consume remaining input samples */
  for (i = 77; i < 512; i++) {
    newCutoff = 0.0;
    for (j = 0; j < 78; j++) {
      newCutoff += rtb_OriginalSignal[i - j] *
        filterbankFIR_P.GeneratedFilterBlock_Coeffic_ns[j];
    }

    rtb_GeneratedFilterBlock[i] = newCutoff;
  }

  /* Update delay line for next frame */
  for (i = 0; i < 77; i++) {
    filterbankFIR_DW.GeneratedFilterBlock_states_av[76 - i] =
      rtb_OriginalSignal[i + 435];
  }

  /* Sum: '<S1>/Sum2' incorporates:
   *  Gain: '<S9>/Gain16'
   *  Sum: '<S4>/Add'
   *  Sum: '<S5>/Add'
   *  Sum: '<S6>/Add'
   *  Sum: '<S7>/Add'
   *  Sum: '<S8>/Add'
   */
  for (i = 0; i < 512; i++) {
    rtb_Sum2[i] = ((((((filterbankFIR_B.Compressor.Compressor2[i] +
                        filterbankFIR_B.Compressor1.Compressor2[i]) +
                       filterbankFIR_B.Compressor2.Compressor2[i]) +
                      ((filterbankFIR_B.Compressor_p.Compressor2[i] +
                        filterbankFIR_B.Compressor1_p.Compressor2[i]) +
                       filterbankFIR_B.Compressor2_p.Compressor2[i])) +
                     ((filterbankFIR_B.Compressor_pn.Compressor2[i] +
                       filterbankFIR_B.Compressor1_pn.Compressor2[i]) +
                      rtb_Gain1[i])) + ((rtb_Gain2[i] + rtb_Gain3[i]) +
      rtb_Gain4[i])) + ((rtb_Gain5[i] + rtb_Gain6[i]) + rtb_Gain7[i])) +
      filterbankFIR_P.gainConsts[15] * rtb_GeneratedFilterBlock[i];
  }

  /* End of Sum: '<S1>/Sum2' */
  /* S-Function (sdspfft2): '<Root>/FFT1' */
  MWDSPCG_FFT_Interleave_R2BR_D(&rtb_Sum2[0U], &filterbankFIR_B.FFT1[0U], 1, 512);
  MWDSPCG_R2DIT_TBLS_Z(&filterbankFIR_B.FFT1[0U], 1, 512, 256, 0,
                       &filterbankFIR_ConstP.FFT1_TwiddleTable[0U], 2, false);
  MWDSPCG_FFT_DblLen_Z(&filterbankFIR_B.FFT1[0U], 1, 512,
                       &filterbankFIR_ConstP.FFT1_TwiddleTable[0U], 1);

  /* ComplexToMagnitudeAngle: '<Root>/Complex to Magnitude-Angle1' */
  for (i = 0; i < 512; i++) {
    filterbankFIR_B.ComplextoMagnitudeAngle1[i] = rt_hypotd_snf
      (filterbankFIR_B.FFT1[i].re, filterbankFIR_B.FFT1[i].im);
  }

  /* End of ComplexToMagnitudeAngle: '<Root>/Complex to Magnitude-Angle1' */

  /* Matfile logging */
  rt_UpdateTXYLogVars(filterbankFIR_M->rtwLogInfo,
                      (&filterbankFIR_M->Timing.taskTime0));

  /* signal main to stop simulation */
  {                                    /* Sample time: [0.011609977324263039s, 0.0s] */
    if ((rtmGetTFinal(filterbankFIR_M)!=-1) &&
        !((rtmGetTFinal(filterbankFIR_M)-filterbankFIR_M->Timing.taskTime0) >
          filterbankFIR_M->Timing.taskTime0 * (DBL_EPSILON))) {
      rtmSetErrorStatus(filterbankFIR_M, "Simulation finished");
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++filterbankFIR_M->Timing.clockTick0)) {
    ++filterbankFIR_M->Timing.clockTickH0;
  }

  filterbankFIR_M->Timing.taskTime0 = filterbankFIR_M->Timing.clockTick0 *
    filterbankFIR_M->Timing.stepSize0 + filterbankFIR_M->Timing.clockTickH0 *
    filterbankFIR_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void filterbankFIR_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)filterbankFIR_M, 0,
                sizeof(RT_MODEL_filterbankFIR_T));
  rtmSetTFinal(filterbankFIR_M, 0.0);
  filterbankFIR_M->Timing.stepSize0 = 0.011609977324263039;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    filterbankFIR_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(filterbankFIR_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(filterbankFIR_M->rtwLogInfo, (NULL));
    rtliSetLogT(filterbankFIR_M->rtwLogInfo, "");
    rtliSetLogX(filterbankFIR_M->rtwLogInfo, "");
    rtliSetLogXFinal(filterbankFIR_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(filterbankFIR_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(filterbankFIR_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(filterbankFIR_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(filterbankFIR_M->rtwLogInfo, 1);
    rtliSetLogY(filterbankFIR_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(filterbankFIR_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(filterbankFIR_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &filterbankFIR_B), 0,
                sizeof(B_filterbankFIR_T));

  /* states (dwork) */
  (void) memset((void *)&filterbankFIR_DW, 0,
                sizeof(DW_filterbankFIR_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(filterbankFIR_M->rtwLogInfo, 0.0,
    rtmGetTFinal(filterbankFIR_M), filterbankFIR_M->Timing.stepSize0,
    (&rtmGetErrorStatus(filterbankFIR_M)));

  {
    dsp_PhaseExtractor_filterbank_T *obj;
    dspcodegen_FIRFilter_filterba_T *iobj_0;
    dsp_simulink_VariableBandwidt_T *obj_0;
    dspcodegen_FIRFilter_filter_n_T *obj_1;
    real_T val[31];
    real_T newCutoff;
    real_T pCoefficients[31];
    dsp_simulink_TransferFunction_T *obj_2;
    boolean_T flag;
    real_T a[512];
    real_T b[512];
    int32_T i;
    static const real_T tmp[31] = { 0.0, 0.010926199633097156,
      0.043227271178699567, 0.095491502812526274, 0.16543469682057088,
      0.24999999999999994, 0.34549150281252627, 0.44773576836617329,
      0.55226423163382665, 0.65450849718747373, 0.74999999999999989,
      0.8345653031794289, 0.90450849718747373, 0.95677272882130049,
      0.98907380036690284, 1.0, 0.98907380036690284, 0.95677272882130049,
      0.90450849718747373, 0.8345653031794289, 0.74999999999999989,
      0.65450849718747373, 0.55226423163382665, 0.44773576836617329,
      0.34549150281252627, 0.24999999999999994, 0.16543469682057088,
      0.095491502812526274, 0.043227271178699567, 0.010926199633097156, 0.0 };

    static const real_T tmp_0[512] = { 0.0, 3.7649080427748505E-5,
      0.00015059065189787502, 0.00033880770582522812, 0.0006022718974137975,
      0.00094094354992541041, 0.0013547716606548965, 0.0018436939086109994,
      0.0024076366639015356, 0.0030465149988219697, 0.0037602327006450165,
      0.0045486822861099951, 0.0054117450176094928, 0.0063492909210707826,
      0.0073611788055293892, 0.0084472562843918575, 0.0096073597983847847,
      0.010841314640186173, 0.012148934980735715, 0.013530023897219912,
      0.014984373402728013, 0.016511764477573965, 0.01811196710228008,
      0.019784740292217107, 0.021529832133895588, 0.023346979822903069,
      0.025235909703481663, 0.02719633730973936, 0.029227967408489597,
      0.03133049404371252, 0.033503600582630522, 0.035746959763392205,
      0.038060233744356631, 0.040443074154971115, 0.042895122148234655,
      0.04541600845473881, 0.048005353438278331, 0.050662767153023092,
      0.053387849402242338, 0.056180189798573033, 0.059039367825822475,
      0.0619649529022967, 0.06495650444564427, 0.068013571939206652,
      0.071135694999863941, 0.0743224034473674, 0.077573217375146442,
      0.080887647222580961, 0.084265193848727382, 0.087705348607487355,
      0.091207593424208144, 0.094771400873702616, 0.098396234259677529,
      0.10208154769555822, 0.10582678618669683, 0.10963138571395276,
      0.1134947733186315, 0.11741636718877052, 0.12139557674675772,
      0.12543180273827031, 0.12952443732252039, 0.13367286416379359,
      0.1378764585242665, 0.1421345873580907, 0.14644660940672621,
      0.15081187529551354, 0.15522972763146653, 0.15969950110227343,
      0.16422052257649078, 0.16879211120491411, 0.17341357852311157,
      0.17808422855510425, 0.18280335791817726, 0.18757025592880677,
      0.19238420470968659, 0.19724447929783723, 0.20215034775378327,
      0.20710107127178057, 0.21209590429107733, 0.21713409460819338,
      0.22221488349019886, 0.22733750578897677, 0.23250119005645137,
      0.23770515866076558, 0.24294862790338917, 0.24823080813714121,
      0.25355090388510793, 0.25890811396043856, 0.2643016315870011,
      0.26973064452088, 0.2751943351726967, 0.28069188073073614,
      0.2862224532848589, 0.29178521995118134, 0.29737934299750507,
      0.30300397996947592, 0.30865828381745508, 0.3143414030240812,
      0.32005248173250589, 0.32579065987528277, 0.33155507330389,
      0.33734485391886848, 0.34315912980055419, 0.3489970253403859,
      0.35485766137276886, 0.3607401553074735, 0.36664362126255079,
      0.37256717019774266, 0.378509910048368, 0.38447094585966435,
      0.39044937992156514, 0.39644431190389073, 0.40245483899193585,
      0.40848005602242948, 0.41451905561984931, 0.4205709283330693,
      0.42663476277231915, 0.43270964574643689, 0.43879466240039189,
      0.44488889635305839, 0.45099142983521961, 0.45710134382778006,
      0.46321771820016627, 0.46933963184889566, 0.47546616283629095,
      0.48159638852932052, 0.48772938573854385, 0.49386423085714004,
      0.49999999999999994, 0.50613576914286, 0.512270614261456,
      0.51840361147067948, 0.524533837163709, 0.53066036815110429,
      0.53678228179983367, 0.54289865617221988, 0.54900857016478033,
      0.55511110364694149, 0.56120533759960811, 0.567290354253563,
      0.57336523722768085, 0.57942907166693058, 0.58548094438015064,
      0.59151994397757046, 0.5975451610080641, 0.60355568809610927,
      0.60955062007843486, 0.61552905414033554, 0.62149008995163191,
      0.62743282980225723, 0.63335637873744921, 0.6392598446925265,
      0.645142338627231, 0.651002974659614, 0.6568408701994457,
      0.66265514608113141, 0.66844492669611, 0.67420934012471723,
      0.67994751826749411, 0.6856585969759188, 0.69134171618254481,
      0.696996020030524, 0.70262065700249487, 0.70821478004881855,
      0.71377754671514093, 0.71930811926926363, 0.72480566482730335,
      0.73026935547912009, 0.73569836841299885, 0.74109188603956133,
      0.746449096114892, 0.75176919186285873, 0.75705137209661078,
      0.76229484133923431, 0.76749880994354847, 0.77266249421102318,
      0.777785116509801, 0.78286590539180656, 0.78790409570892272,
      0.79289892872821943, 0.79784965224621662, 0.80275552070216272,
      0.80761579529031335, 0.81242974407119317, 0.81719664208182263,
      0.82191577144489569, 0.82658642147688832, 0.831207888795086,
      0.83577947742350922, 0.84030049889772651, 0.84477027236853353,
      0.84918812470448635, 0.85355339059327373, 0.8578654126419093,
      0.86212354147573333, 0.86632713583620635, 0.87047556267747939,
      0.87456819726172963, 0.87860442325324239, 0.88258363281122953,
      0.88650522668136844, 0.89036861428604719, 0.89417321381330317,
      0.89791845230444167, 0.90160376574032242, 0.90522859912629738,
      0.90879240657579174, 0.91229465139251253, 0.91573480615127267,
      0.919112352777419, 0.92242678262485356, 0.9256775965526326,
      0.928864305000136, 0.93198642806079335, 0.93504349555435562,
      0.93803504709770325, 0.94096063217417747, 0.94381981020142691,
      0.94661215059775761, 0.94933723284697691, 0.95199464656172172,
      0.95458399154526119, 0.95710487785176535, 0.95955692584502894,
      0.96193976625564337, 0.96425304023660774, 0.96649639941736942,
      0.96866950595628742, 0.97077203259151035, 0.97280366269026053,
      0.97476409029651834, 0.97665302017709688, 0.97847016786610441,
      0.98021525970778289, 0.98188803289772, 0.983488235522426,
      0.985015626597272, 0.98646997610278, 0.98785106501926423,
      0.98915868535981377, 0.99039264020161522, 0.9915527437156082,
      0.99263882119447056, 0.99365070907892916, 0.99458825498239056,
      0.99545131771389, 0.996239767299355, 0.99695348500117809,
      0.99759236333609835, 0.998156306091389, 0.9986452283393451,
      0.99905905645007453, 0.9993977281025862, 0.99966119229417472,
      0.99984940934810207, 0.99996235091957231, 1.0, 0.99996235091957231,
      0.99984940934810207, 0.99966119229417472, 0.9993977281025862,
      0.99905905645007453, 0.9986452283393451, 0.998156306091389,
      0.99759236333609835, 0.99695348500117809, 0.996239767299355,
      0.99545131771389, 0.99458825498239056, 0.99365070907892916,
      0.99263882119447056, 0.9915527437156082, 0.99039264020161522,
      0.98915868535981377, 0.98785106501926423, 0.98646997610278,
      0.985015626597272, 0.983488235522426, 0.98188803289772,
      0.98021525970778289, 0.97847016786610441, 0.97665302017709688,
      0.97476409029651834, 0.97280366269026053, 0.97077203259151035,
      0.96866950595628742, 0.96649639941736942, 0.96425304023660774,
      0.96193976625564337, 0.95955692584502894, 0.95710487785176535,
      0.95458399154526119, 0.95199464656172172, 0.94933723284697691,
      0.94661215059775761, 0.94381981020142691, 0.94096063217417747,
      0.93803504709770325, 0.93504349555435562, 0.93198642806079335,
      0.928864305000136, 0.9256775965526326, 0.92242678262485356,
      0.919112352777419, 0.91573480615127267, 0.91229465139251253,
      0.90879240657579174, 0.90522859912629738, 0.90160376574032242,
      0.89791845230444167, 0.89417321381330317, 0.89036861428604719,
      0.88650522668136844, 0.88258363281122953, 0.87860442325324239,
      0.87456819726172963, 0.87047556267747939, 0.86632713583620635,
      0.86212354147573333, 0.8578654126419093, 0.85355339059327373,
      0.84918812470448635, 0.84477027236853353, 0.84030049889772651,
      0.83577947742350922, 0.831207888795086, 0.82658642147688832,
      0.82191577144489569, 0.81719664208182263, 0.81242974407119317,
      0.80761579529031335, 0.80275552070216272, 0.79784965224621662,
      0.79289892872821943, 0.78790409570892272, 0.78286590539180656,
      0.777785116509801, 0.77266249421102318, 0.76749880994354847,
      0.76229484133923431, 0.75705137209661078, 0.75176919186285873,
      0.746449096114892, 0.74109188603956133, 0.73569836841299885,
      0.73026935547912009, 0.72480566482730335, 0.71930811926926363,
      0.71377754671514093, 0.70821478004881855, 0.70262065700249487,
      0.696996020030524, 0.69134171618254481, 0.6856585969759188,
      0.67994751826749411, 0.67420934012471723, 0.66844492669611,
      0.66265514608113141, 0.6568408701994457, 0.651002974659614,
      0.645142338627231, 0.6392598446925265, 0.63335637873744921,
      0.62743282980225723, 0.62149008995163191, 0.61552905414033554,
      0.60955062007843486, 0.60355568809610927, 0.5975451610080641,
      0.59151994397757046, 0.58548094438015064, 0.57942907166693058,
      0.57336523722768085, 0.567290354253563, 0.56120533759960811,
      0.55511110364694149, 0.54900857016478033, 0.54289865617221988,
      0.53678228179983367, 0.53066036815110429, 0.524533837163709,
      0.51840361147067948, 0.512270614261456, 0.50613576914286,
      0.49999999999999994, 0.49386423085714004, 0.48772938573854385,
      0.48159638852932052, 0.47546616283629095, 0.46933963184889566,
      0.46321771820016627, 0.45710134382778006, 0.45099142983521961,
      0.44488889635305839, 0.43879466240039189, 0.43270964574643689,
      0.42663476277231915, 0.4205709283330693, 0.41451905561984931,
      0.40848005602242948, 0.40245483899193585, 0.39644431190389073,
      0.39044937992156514, 0.38447094585966435, 0.378509910048368,
      0.37256717019774266, 0.36664362126255079, 0.3607401553074735,
      0.35485766137276886, 0.3489970253403859, 0.34315912980055419,
      0.33734485391886848, 0.33155507330389, 0.32579065987528277,
      0.32005248173250589, 0.3143414030240812, 0.30865828381745508,
      0.30300397996947592, 0.29737934299750507, 0.29178521995118134,
      0.2862224532848589, 0.28069188073073614, 0.2751943351726967,
      0.26973064452088, 0.2643016315870011, 0.25890811396043856,
      0.25355090388510793, 0.24823080813714121, 0.24294862790338917,
      0.23770515866076558, 0.23250119005645137, 0.22733750578897677,
      0.22221488349019886, 0.21713409460819338, 0.21209590429107733,
      0.20710107127178057, 0.20215034775378327, 0.19724447929783723,
      0.19238420470968659, 0.18757025592880677, 0.18280335791817726,
      0.17808422855510425, 0.17341357852311157, 0.16879211120491411,
      0.16422052257649078, 0.15969950110227343, 0.15522972763146653,
      0.15081187529551354, 0.14644660940672621, 0.1421345873580907,
      0.1378764585242665, 0.13367286416379359, 0.12952443732252039,
      0.12543180273827031, 0.12139557674675772, 0.11741636718877052,
      0.1134947733186315, 0.10963138571395276, 0.10582678618669683,
      0.10208154769555822, 0.098396234259677529, 0.094771400873702616,
      0.091207593424208144, 0.087705348607487355, 0.084265193848727382,
      0.080887647222580961, 0.077573217375146442, 0.0743224034473674,
      0.071135694999863941, 0.068013571939206652, 0.06495650444564427,
      0.0619649529022967, 0.059039367825822475, 0.056180189798573033,
      0.053387849402242338, 0.050662767153023092, 0.048005353438278331,
      0.04541600845473881, 0.042895122148234655, 0.040443074154971115,
      0.038060233744356631, 0.035746959763392205, 0.033503600582630522,
      0.03133049404371252, 0.029227967408489597, 0.02719633730973936,
      0.025235909703481663, 0.023346979822903069, 0.021529832133895588,
      0.019784740292217107, 0.01811196710228008, 0.016511764477573965,
      0.014984373402728013, 0.013530023897219912, 0.012148934980735715,
      0.010841314640186173, 0.0096073597983847847, 0.0084472562843918575,
      0.0073611788055293892, 0.0063492909210707826, 0.0054117450176094928,
      0.0045486822861099951, 0.0037602327006450165, 0.0030465149988219697,
      0.0024076366639015356, 0.0018436939086109994, 0.0013547716606548965,
      0.00094094354992541041, 0.0006022718974137975, 0.00033880770582522812,
      0.00015059065189787502, 3.7649080427748505E-5 };

    /* Start for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
    /* Start for MATLABSystem: '<S48>/TFE' */
    filterbankFIR_DW.obj.pSpectrumEstimator.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj.matlabCodegenIsDeleted = true;
    obj_2 = &filterbankFIR_DW.obj;
    filterbankFIR_DW.obj.isInitialized = 0;
    obj_2->pSpectrumEstimator.isInitialized = 0;
    obj_2->pSpectrumEstimator.matlabCodegenIsDeleted = false;
    obj_2->pCrossSpectrumEstimator.isInitialized = 0;
    filterbankFIR_DW.obj.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.objisempty = true;
    obj_2 = &filterbankFIR_DW.obj;
    filterbankFIR_DW.obj.isSetupComplete = false;
    filterbankFIR_DW.obj.isInitialized = 1;
    flag = (obj_2->pSpectrumEstimator.isInitialized == 1);
    if (flag) {
      obj_2->pSpectrumEstimator.TunablePropsChanged = true;
    }

    obj_2->pSpectrumEstimator.ReferenceLoad = 1.0;
    obj_2->pSpectrumEstimator.isSetupComplete = false;
    obj_2->pSpectrumEstimator.isInitialized = 1;
    for (i = 0; i < 512; i++) {
      obj_2->pSpectrumEstimator.pWindowData[i] = tmp_0[i];
    }

    for (i = 0; i < 512; i++) {
      a[i] = obj_2->pSpectrumEstimator.pWindowData[i];
    }

    for (i = 0; i < 512; i++) {
      b[i] = obj_2->pSpectrumEstimator.pWindowData[i];
    }

    newCutoff = 0.0;
    for (i = 0; i < 512; i++) {
      newCutoff += a[i] * b[i];
    }

    obj_2->pSpectrumEstimator.pWindowPower = newCutoff;
    for (i = 0; i < 257; i++) {
      obj_2->pSpectrumEstimator.pW[i] = 0.001953125 * (real_T)i;
    }

    obj_2->pSpectrumEstimator.pReferenceLoad =
      obj_2->pSpectrumEstimator.ReferenceLoad;
    obj_2->pSpectrumEstimator.isSetupComplete = true;
    obj_2->pSpectrumEstimator.TunablePropsChanged = false;
    obj_2->pCrossSpectrumEstimator.isSetupComplete = false;
    obj_2->pCrossSpectrumEstimator.isInitialized = 1;
    for (i = 0; i < 512; i++) {
      obj_2->pCrossSpectrumEstimator.pWindowData[i] = tmp_0[i];
    }

    for (i = 0; i < 512; i++) {
      a[i] = obj_2->pCrossSpectrumEstimator.pWindowData[i];
    }

    for (i = 0; i < 512; i++) {
      b[i] = obj_2->pCrossSpectrumEstimator.pWindowData[i];
    }

    newCutoff = 0.0;
    for (i = 0; i < 512; i++) {
      newCutoff += a[i] * b[i];
    }

    obj_2->pCrossSpectrumEstimator.pWindowPower = newCutoff;
    obj_2->pCrossSpectrumEstimator.isSetupComplete = true;
    filterbankFIR_DW.obj.pFrameCounter = 0.0;
    filterbankFIR_DW.obj.pFrameDelay = 0.0;
    filterbankFIR_DW.obj.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<S48>/TFE' */
    /* End of Start for SubSystem: '<S3>/Discrete Transfer Function Estimator' */

    /* Start for MATLABSystem: '<S3>/Phase Extractor' */
    filterbankFIR_DW.obj_g.pPhaseDifferentiator.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_g.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_g.isInitialized = 0;
    filterbankFIR_DW.obj_g.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.objisempty_o = true;
    obj = &filterbankFIR_DW.obj_g;
    filterbankFIR_DW.obj_g.isSetupComplete = false;
    filterbankFIR_DW.obj_g.isInitialized = 1;
    obj->pPhaseDifferentiator.isInitialized = 0;
    obj->pPhaseDifferentiator.pNumChans = -1;
    obj->pPhaseDifferentiator.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.obj_g.isSetupComplete = true;

    /* Start for MATLABSystem: '<S3>/Differentiator Filter' */
    filterbankFIR_DW.gobj_1.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.gobj_0.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_m.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_m.isInitialized = 0;
    filterbankFIR_DW.obj_m.NumChannels = -1;
    filterbankFIR_DW.obj_m.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.objisempty_g = true;
    iobj_0 = &filterbankFIR_DW.gobj_0;
    filterbankFIR_DW.obj_m.isSetupComplete = false;
    filterbankFIR_DW.obj_m.isInitialized = 1;
    filterbankFIR_DW.gobj_0.isInitialized = 0;

    /* System object Constructor function: dsp.FIRFilter */
    iobj_0->cSFunObject.P0_InitialStates = 0.0;
    iobj_0->cSFunObject.P1_Coefficients[0] = -0.13106369956607014;
    iobj_0->cSFunObject.P1_Coefficients[1] = 1.3091932765180649;
    iobj_0->cSFunObject.P1_Coefficients[2] = -1.3091932765180649;
    iobj_0->cSFunObject.P1_Coefficients[3] = 0.13106369956607014;
    filterbankFIR_DW.gobj_0.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.obj_m.FilterObj = &filterbankFIR_DW.gobj_0;
    filterbankFIR_DW.obj_m.NumChannels = 1;
    filterbankFIR_DW.obj_m.isSetupComplete = true;

    /* Start for MATLABSystem: '<Root>/Lowpass Filter' */
    filterbankFIR_DW.obj_k.pfilter.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_k.matlabCodegenIsDeleted = true;
    filterbankFIR_DW.obj_k.isInitialized = 0;
    filterbankFIR_DW.obj_k.pNumChannels = -1;
    filterbankFIR_DW.obj_k.matlabCodegenIsDeleted = false;
    filterbankFIR_DW.objisempty_l = true;
    flag = (filterbankFIR_DW.obj_k.isInitialized == 1);
    if (flag) {
      filterbankFIR_DW.obj_k.TunablePropsChanged = true;
    }

    filterbankFIR_DW.obj_k.CutoffFrequency = filterbankFIR_P.LowpassFilter_fc;
    obj_0 = &filterbankFIR_DW.obj_k;
    filterbankFIR_DW.obj_k.isSetupComplete = false;
    filterbankFIR_DW.obj_k.isInitialized = 1;
    filterbankFIR_DW.obj_k.pNumChannels = 1;
    obj_1 = &filterbankFIR_DW.obj_k.pfilter;
    obj_0->pfilter.isInitialized = 0;

    /* System object Constructor function: dsp.FIRFilter */
    obj_1->cSFunObject.P0_InitialStates = 0.0;
    memcpy(&val[0], &tmp[0], 31U * sizeof(real_T));
    for (i = 0; i < 31; i++) {
      obj_1->cSFunObject.P1_Coefficients[i] = val[i];
    }

    for (i = 0; i < 31; i++) {
      obj_0->pfilter.Numerator[i] = val[i];
    }

    obj_0->pfilter.matlabCodegenIsDeleted = false;
    obj_0->pfilter.isSetupComplete = false;
    obj_0->pfilter.isInitialized = 1;
    obj_0->pfilter.isSetupComplete = true;
    newCutoff = filterbankFIR_DW.obj_k.CutoffFrequency * 2.0;
    newCutoff /= filterbankFIR_SampleRate;
    for (i = 0; i < 31; i++) {
      pCoefficients[i] = 0.0;
      if (-15 + i == 0) {
        pCoefficients[i] = newCutoff;
      } else {
        pCoefficients[i] = sin((-15.0 + (real_T)i) * newCutoff *
          3.1415926535897931) / ((-15.0 + (real_T)i) * 3.1415926535897931);
      }
    }

    obj_1 = &filterbankFIR_DW.obj_k.pfilter;
    for (i = 0; i < 31; i++) {
      val[i] = tmp[i] * pCoefficients[i];
    }

    for (i = 0; i < 31; i++) {
      obj_1->cSFunObject.P1_Coefficients[i] = val[i];
    }

    for (i = 0; i < 31; i++) {
      obj_0->pfilter.Numerator[i] = val[i];
    }

    memcpy(&pCoefficients[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0], 31U *
           sizeof(real_T));
    if (filterbankFIR_sum(pCoefficients) > 0.0) {
      memcpy(&val[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0], 31U * sizeof
             (real_T));
      memcpy(&pCoefficients[0], &filterbankFIR_DW.obj_k.pfilter.Numerator[0],
             31U * sizeof(real_T));
      newCutoff = filterbankFIR_sum(pCoefficients);
      for (i = 0; i < 31; i++) {
        val[i] /= newCutoff;
      }

      obj_1 = &filterbankFIR_DW.obj_k.pfilter;
      for (i = 0; i < 31; i++) {
        obj_1->cSFunObject.P1_Coefficients[i] = val[i];
      }

      for (i = 0; i < 31; i++) {
        obj_0->pfilter.Numerator[i] = val[i];
      }
    }

    filterbankFIR_DW.obj_k.isSetupComplete = true;
    filterbankFIR_DW.obj_k.TunablePropsChanged = false;

    /* End of Start for MATLABSystem: '<Root>/Lowpass Filter' */
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor,
      &filterbankFIR_P.Compressor);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor1,
      &filterbankFIR_P.Compressor1);

    /* Start for MATLABSystem: '<S4>/Compressor2' */
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor2,
      &filterbankFIR_P.Compressor2);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor_p,
      &filterbankFIR_P.Compressor_p);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor1_p,
      &filterbankFIR_P.Compressor1_p);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor2_p,
      &filterbankFIR_P.Compressor2_p);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor_pn,
      &filterbankFIR_P.Compressor_pn);
    filterbankFIR_Compressor2_Start(&filterbankFIR_DW.Compressor1_pn,
      &filterbankFIR_P.Compressor1_pn);
  }

  {
    dsp_simulink_VariableBandwidt_T *obj;
    dspcodegen_FIRFilter_filter_n_T *obj_0;
    int32_T i;

    /* InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source1' */
    filterbankFIR_DW.RandomSource1_SEED_DWORK = 81472U;
    RandSrcInitState_GZ(&filterbankFIR_DW.RandomSource1_SEED_DWORK,
                        filterbankFIR_DW.RandomSource1_STATE_DWORK, 1);

    /* InitializeConditions for DiscreteFir: '<S32>/Generated Filter Block' */
    for (i = 0; i < 702; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialSta;
    }

    /* End of InitializeConditions for DiscreteFir: '<S32>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S33>/Generated Filter Block' */
    for (i = 0; i < 1583; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_f[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_n;
    }

    /* End of InitializeConditions for DiscreteFir: '<S33>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S34>/Generated Filter Block' */
    for (i = 0; i < 1262; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_c[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_e;
    }

    /* End of InitializeConditions for DiscreteFir: '<S34>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S36>/Generated Filter Block' */
    for (i = 0; i < 989; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_l[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_f;
    }

    /* End of InitializeConditions for DiscreteFir: '<S36>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S35>/Generated Filter Block' */
    for (i = 0; i < 785; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_j[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_m;
    }

    /* End of InitializeConditions for DiscreteFir: '<S35>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S37>/Generated Filter Block' */
    for (i = 0; i < 627; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_e[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_i;
    }

    /* End of InitializeConditions for DiscreteFir: '<S37>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S40>/Generated Filter Block' */
    for (i = 0; i < 498; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_eh[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_me;
    }

    /* End of InitializeConditions for DiscreteFir: '<S40>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S38>/Generated Filter Block' */
    for (i = 0; i < 394; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_jw[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_p;
    }

    /* End of InitializeConditions for DiscreteFir: '<S38>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S39>/Generated Filter Block' */
    for (i = 0; i < 313; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_i[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_na;
    }

    /* End of InitializeConditions for DiscreteFir: '<S39>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S42>/Generated Filter Block' */
    for (i = 0; i < 249; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_fb[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_h;
    }

    /* End of InitializeConditions for DiscreteFir: '<S42>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S41>/Generated Filter Block' */
    for (i = 0; i < 198; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_d[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_k;
    }

    /* End of InitializeConditions for DiscreteFir: '<S41>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S43>/Generated Filter Block' */
    for (i = 0; i < 157; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_g[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_fa;
    }

    /* End of InitializeConditions for DiscreteFir: '<S43>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S46>/Generated Filter Block' */
    for (i = 0; i < 153; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_b[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_h4;
    }

    /* End of InitializeConditions for DiscreteFir: '<S46>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S44>/Generated Filter Block' */
    for (i = 0; i < 122; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_fi[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_d;
    }

    /* End of InitializeConditions for DiscreteFir: '<S44>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S45>/Generated Filter Block' */
    for (i = 0; i < 97; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_o[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_g;
    }

    /* End of InitializeConditions for DiscreteFir: '<S45>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S47>/Generated Filter Block' */
    for (i = 0; i < 77; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_bf[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_c;
    }

    /* End of InitializeConditions for DiscreteFir: '<S47>/Generated Filter Block' */

    /* InitializeConditions for S-Function (sdspimpgen2): '<Root>/Discrete  Impulse' */

    /* DSP System Toolbox Discrete Impulse Generator (sdspimpgen2) - '<Root>/Discrete  Impulse' */
    filterbankFIR_DW.DiscreteImpulse_COUNT = ((int_T)0) + 1;

    /* InitializeConditions for DiscreteFir: '<S10>/Generated Filter Block' */
    for (i = 0; i < 989; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_n[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_cv;
    }

    /* End of InitializeConditions for DiscreteFir: '<S10>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S11>/Generated Filter Block' */
    for (i = 0; i < 1583; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_cf[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_j;
    }

    /* End of InitializeConditions for DiscreteFir: '<S11>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S12>/Generated Filter Block' */
    for (i = 0; i < 1262; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_a[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_nu;
    }

    /* End of InitializeConditions for DiscreteFir: '<S12>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S14>/Generated Filter Block' */
    for (i = 0; i < 989; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_bs[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_ci;
    }

    /* End of InitializeConditions for DiscreteFir: '<S14>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S13>/Generated Filter Block' */
    for (i = 0; i < 785; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_et[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_pu;
    }

    /* End of InitializeConditions for DiscreteFir: '<S13>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S15>/Generated Filter Block' */
    for (i = 0; i < 627; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_lf[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_kq;
    }

    /* End of InitializeConditions for DiscreteFir: '<S15>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S18>/Generated Filter Block' */
    for (i = 0; i < 498; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_ed[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_hv;
    }

    /* End of InitializeConditions for DiscreteFir: '<S18>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S16>/Generated Filter Block' */
    for (i = 0; i < 394; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_nq[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_ip;
    }

    /* End of InitializeConditions for DiscreteFir: '<S16>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S17>/Generated Filter Block' */
    for (i = 0; i < 313; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_ec[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_b;
    }

    /* End of InitializeConditions for DiscreteFir: '<S17>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S20>/Generated Filter Block' */
    for (i = 0; i < 249; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_i2[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_fu;
    }

    /* End of InitializeConditions for DiscreteFir: '<S20>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S19>/Generated Filter Block' */
    for (i = 0; i < 198; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_ou[i] =
        filterbankFIR_P.GeneratedFilterBlock_InitialS_l;
    }

    /* End of InitializeConditions for DiscreteFir: '<S19>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S21>/Generated Filter Block' */
    for (i = 0; i < 157; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_fu[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_m0;
    }

    /* End of InitializeConditions for DiscreteFir: '<S21>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S24>/Generated Filter Block' */
    for (i = 0; i < 153; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_iu[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_h1;
    }

    /* End of InitializeConditions for DiscreteFir: '<S24>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S22>/Generated Filter Block' */
    for (i = 0; i < 122; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_ck[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_nz;
    }

    /* End of InitializeConditions for DiscreteFir: '<S22>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S23>/Generated Filter Block' */
    for (i = 0; i < 97; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_bt[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_ns;
    }

    /* End of InitializeConditions for DiscreteFir: '<S23>/Generated Filter Block' */

    /* InitializeConditions for DiscreteFir: '<S25>/Generated Filter Block' */
    for (i = 0; i < 77; i++) {
      filterbankFIR_DW.GeneratedFilterBlock_states_av[i] =
        filterbankFIR_P.GeneratedFilterBlock_Initial_dg;
    }

    /* End of InitializeConditions for DiscreteFir: '<S25>/Generated Filter Block' */
    /* SystemInitialize for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
    /* InitializeConditions for MATLABSystem: '<S48>/TFE' */
    if (filterbankFIR_DW.obj.pSpectrumEstimator.isInitialized == 1) {
      filterbankFIR_DW.obj.pSpectrumEstimator.pNumAvgsCounter = 0.0;
      filterbankFIR_DW.obj.pSpectrumEstimator.pNewPeriodogramIdx = 0.0;
      memset(&filterbankFIR_DW.obj.pSpectrumEstimator.pPeriodogramMatrix[0], 0,
             5120U * sizeof(real_T));
    }

    if (filterbankFIR_DW.obj.pCrossSpectrumEstimator.isInitialized == 1) {
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.pNumAvgsCounter = 0.0;
      filterbankFIR_DW.obj.pCrossSpectrumEstimator.pNewPeriodogramIdx = 0.0;
      memset(&filterbankFIR_DW.obj.pCrossSpectrumEstimator.pPeriodogramMatrix[0],
             0, 5120U * sizeof(creal_T));
    }

    /* End of InitializeConditions for MATLABSystem: '<S48>/TFE' */
    /* End of SystemInitialize for SubSystem: '<S3>/Discrete Transfer Function Estimator' */

    /* InitializeConditions for MATLABSystem: '<S3>/Phase Extractor' */
    if (filterbankFIR_DW.obj_g.pPhaseDifferentiator.isInitialized == 1) {
      filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.re = 1.0;
      filterbankFIR_DW.obj_g.pPhaseDifferentiator.pLastSample.im = 0.0;
    }

    /* End of InitializeConditions for MATLABSystem: '<S3>/Phase Extractor' */

    /* InitializeConditions for MATLABSystem: '<S3>/Differentiator Filter' */
    if (filterbankFIR_DW.obj_m.FilterObj->isInitialized == 1) {
      /* System object Initialization function: dsp.FIRFilter */
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[0] =
        filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[1] =
        filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
      filterbankFIR_DW.obj_m.FilterObj->cSFunObject.W0_states[2] =
        filterbankFIR_DW.obj_m.FilterObj->cSFunObject.P0_InitialStates;
    }

    /* End of InitializeConditions for MATLABSystem: '<S3>/Differentiator Filter' */

    /* InitializeConditions for MATLABSystem: '<Root>/Lowpass Filter' */
    obj = &filterbankFIR_DW.obj_k;
    obj_0 = &filterbankFIR_DW.obj_k.pfilter;
    if (obj->pfilter.isInitialized == 1) {
      /* System object Initialization function: dsp.FIRFilter */
      for (i = 0; i < 30; i++) {
        obj_0->cSFunObject.W0_states[i] = obj_0->cSFunObject.P0_InitialStates;
      }
    }

    /* End of InitializeConditions for MATLABSystem: '<Root>/Lowpass Filter' */
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor1);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor2);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor_p);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor1_p);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor2_p);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor_pn);
    filterbankFIR_Compressor2_Init(&filterbankFIR_DW.Compressor1_pn);
  }
}

/* Model terminate function */
void filterbankFIR_terminate(void)
{
  /* Terminate for Atomic SubSystem: '<S3>/Discrete Transfer Function Estimator' */
  /* Terminate for MATLABSystem: '<S48>/TFE' */
  matlabCodegenHandle_mat_nzusynb(&filterbankFIR_DW.obj);
  matlabCodegenHandle_ma_nzusynbu(&filterbankFIR_DW.obj.pSpectrumEstimator);

  /* End of Terminate for SubSystem: '<S3>/Discrete Transfer Function Estimator' */

  /* Terminate for MATLABSystem: '<S3>/Phase Extractor' */
  matlabCodegenHandle_matl_nzusyn(&filterbankFIR_DW.obj_g);
  matlabCodegenHandle_matla_nzusy(&filterbankFIR_DW.obj_g.pPhaseDifferentiator);

  /* Terminate for MATLABSystem: '<S3>/Differentiator Filter' */
  matlabCodegenHandle_matlab_nzus(&filterbankFIR_DW.obj_m);
  matlabCodegenHandle_matlabC_nzu(&filterbankFIR_DW.gobj_0);
  matlabCodegenHandle_matlabC_nzu(&filterbankFIR_DW.gobj_1);

  /* Terminate for MATLABSystem: '<Root>/Lowpass Filter' */
  matlabCodegenHandle_matlabCo_nz(&filterbankFIR_DW.obj_k);
  matlabCodegenHandle_matlabCod_n(&filterbankFIR_DW.obj_k.pfilter);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor1);

  /* Terminate for MATLABSystem: '<S4>/Compressor2' */
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor2);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor_p);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor1_p);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor2_p);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor_pn);
  filterbankFIR_Compressor2_Term(&filterbankFIR_DW.Compressor1_pn);
}
