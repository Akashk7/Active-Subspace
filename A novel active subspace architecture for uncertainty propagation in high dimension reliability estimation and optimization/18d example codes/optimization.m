 % clc;clear;close all;
 
%% Optimization
t0 = 4*ones(1,n_design_var);

func_PRS = @(t_PRS)(srgtsPRSEvaluate((normalizeuv(t_PRS,lbx_obj,ubx_obj)*W1_obj_llm), srgt_PRS_obj_1));
[xmin_PRS,fval_PRS,exitflag_PRS,output_PRS]=fmincon(func_PRS,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,W1_ls_04_llm,W1_ls_05_llm,W1_ls_06_llm,W1_ls_07_llm,W1_ls_08_llm,W1_ls_09_llm,W1_ls_10_llm,W1_ls_11_llm,W1_ls_12_llm,...
srgt_PRS_c01,srgt_PRS_c02,srgt_PRS_c03,srgt_PRS_c04,srgt_PRS_c05,srgt_PRS_c06,srgt_PRS_c07,srgt_PRS_c08,srgt_PRS_c09,srgt_PRS_c10,srgt_PRS_c11,srgt_PRS_c12,@srgtsPRSEvaluate));

%% Validation

%bounds
lbx_c01=lbx;
lbx_c02=lbx;
lbx_c03=[lbx(1:16) lbx(18)];
lbx_c04=lbx(1:17);
lbx_c05=lbx;
lbx_c06=lbx;
lbx_c07=lbx;
lbx_c08=lbx;
lbx_c09=lbx;
lbx_c10=lbx;
lbx_c11=lbx(1:15);
lbx_c12=lbx(1:15);

ubx_c01=ubx;
ubx_c02=ubx;
ubx_c03=[ubx(1:16) ubx(18)];
ubx_c04=ubx(1:17);
ubx_c05=ubx;
ubx_c06=ubx;
ubx_c07=ubx;
ubx_c08=ubx;
ubx_c09=ubx;
ubx_c10=ubx;
ubx_c11=ubx(1:15);
ubx_c12=ubx(1:15);

%PRS
xmin_PRS_obj = normalizeuv(xmin_PRS,lbx_obj, ubx_obj);
xmin_PRS_c01=normalizeuv(xmin_PRS, lbx_c01, ubx_c01);                       
xmin_PRS_c02=normalizeuv(xmin_PRS, lbx_c02, ubx_c02);
xmin_PRS_c03=normalizeuv([xmin_PRS(1:16) xmin_PRS(18)], lbx_c03, ubx_c03);
xmin_PRS_c04=normalizeuv(xmin_PRS(1:17), lbx_c04, ubx_c04);
xmin_PRS_c05=normalizeuv(xmin_PRS, lbx_c05, ubx_c05);
xmin_PRS_c06=normalizeuv(xmin_PRS, lbx_c06, ubx_c06);
xmin_PRS_c07=normalizeuv(xmin_PRS, lbx_c07, ubx_c07);
xmin_PRS_c08=normalizeuv(xmin_PRS, lbx_c08, ubx_c08);
xmin_PRS_c09=normalizeuv(xmin_PRS, lbx_c09, ubx_c09);
xmin_PRS_c10=normalizeuv(xmin_PRS, lbx_c10, ubx_c10);
xmin_PRS_c11=normalizeuv(xmin_PRS(1:15), lbx_c11, ubx_c11);
xmin_PRS_c12=normalizeuv(xmin_PRS(1:15), lbx_c12, ubx_c12);

mm_obj_PRS = srgtsPRSEvaluate(xmin_PRS_obj*W1_obj_llm, srgt_PRS_obj_1);
mm_psf_PRS_c01 = srgtsPRSEvaluate(xmin_PRS_c01*W1_ls_01_llm, srgt_PRS_c01);
mm_psf_PRS_c02 = srgtsPRSEvaluate(xmin_PRS_c02*W1_ls_02_llm, srgt_PRS_c02);
mm_psf_PRS_c03 = srgtsPRSEvaluate(xmin_PRS_c03*W1_ls_03_llm, srgt_PRS_c03);
mm_psf_PRS_c04 = srgtsPRSEvaluate(xmin_PRS_c04*W1_ls_04_llm, srgt_PRS_c04);
mm_psf_PRS_c05 = srgtsPRSEvaluate(xmin_PRS_c05*W1_ls_05_llm, srgt_PRS_c05);
mm_psf_PRS_c06 = srgtsPRSEvaluate(xmin_PRS_c06*W1_ls_06_llm, srgt_PRS_c06);
mm_psf_PRS_c07 = srgtsPRSEvaluate(xmin_PRS_c07*W1_ls_07_llm, srgt_PRS_c07);
mm_psf_PRS_c08 = srgtsPRSEvaluate(xmin_PRS_c08*W1_ls_08_llm, srgt_PRS_c08);
mm_psf_PRS_c09 = srgtsPRSEvaluate(xmin_PRS_c09*W1_ls_09_llm, srgt_PRS_c09);
mm_psf_PRS_c10 = srgtsPRSEvaluate(xmin_PRS_c10*W1_ls_10_llm, srgt_PRS_c10);
mm_psf_PRS_c11 = srgtsPRSEvaluate(xmin_PRS_c11*W1_ls_11_llm, srgt_PRS_c11);
mm_psf_PRS_c12 = srgtsPRSEvaluate(xmin_PRS_c12*W1_ls_12_llm, srgt_PRS_c12);

act_obj_PRS = Weightc(xmin_PRS);
act_psf_PRS_c01 = mcspsfconstraint1c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c02 = mcspsfconstraint2c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c03 = mcspsfconstraint3c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c04 = mcspsfconstraint4c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c05 = mcspsfconstraint5c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c06 = mcspsfconstraint6c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c07 = mcspsfconstraint7c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c08 = mcspsfconstraint8c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c09 = mcspsfconstraint9c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c10 = mcspsfconstraint10c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c11 = mcspsfconstraint11c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c12 = mcspsfconstraint12c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);

%% Variance of prediction error [stderr_c0i_PRS_1, i=1,2,3,...,12]
[predvar_c01_PRS,stderr_c01_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c01*W1_ls_01_llm, t1_c01, srgt_PRS_c01));
[predvar_c02_PRS,stderr_c02_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c02*W1_ls_02_llm, t1_c02, srgt_PRS_c02));
[predvar_c03_PRS,stderr_c03_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c03*W1_ls_03_llm, t1_c03, srgt_PRS_c03));
[predvar_c04_PRS,stderr_c04_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c04*W1_ls_04_llm, t1_c04, srgt_PRS_c04));
[predvar_c05_PRS,stderr_c05_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c05*W1_ls_05_llm, t1_c05, srgt_PRS_c05));
[predvar_c06_PRS,stderr_c06_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c06*W1_ls_06_llm, t1_c06, srgt_PRS_c06));
[predvar_c07_PRS,stderr_c07_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c07*W1_ls_07_llm, t1_c07, srgt_PRS_c07));
[predvar_c08_PRS,stderr_c08_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c08*W1_ls_08_llm, t1_c08, srgt_PRS_c08));
[predvar_c09_PRS,stderr_c09_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c09*W1_ls_09_llm, t1_c09, srgt_PRS_c09));
[predvar_c10_PRS,stderr_c10_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c10*W1_ls_10_llm, t1_c10, srgt_PRS_c10));
[predvar_c11_PRS,stderr_c11_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c11*W1_ls_11_llm, t1_c11, srgt_PRS_c11));
[predvar_c12_PRS,stderr_c12_PRS] = (srgtsPRSPredictionVariance(xmin_PRS_c12*W1_ls_12_llm, t1_c12, srgt_PRS_c12));

stderr_c01_PRS_1 = sqrt(predvar_c01_PRS+stderr_c01_PRS);
stderr_c02_PRS_1 = sqrt(predvar_c02_PRS+stderr_c02_PRS);
stderr_c03_PRS_1 = sqrt(predvar_c03_PRS+stderr_c03_PRS);
stderr_c04_PRS_1 = sqrt(predvar_c04_PRS+stderr_c04_PRS);
stderr_c05_PRS_1 = sqrt(predvar_c05_PRS+stderr_c05_PRS);
stderr_c06_PRS_1 = sqrt(predvar_c06_PRS+stderr_c06_PRS);
stderr_c07_PRS_1 = sqrt(predvar_c07_PRS+stderr_c07_PRS);
stderr_c08_PRS_1 = sqrt(predvar_c08_PRS+stderr_c08_PRS);
stderr_c09_PRS_1 = sqrt(predvar_c09_PRS+stderr_c09_PRS);
stderr_c10_PRS_1 = sqrt(predvar_c10_PRS+stderr_c10_PRS);
stderr_c11_PRS_1 = sqrt(predvar_c11_PRS+stderr_c11_PRS);
stderr_c12_PRS_1 = sqrt(predvar_c12_PRS+stderr_c12_PRS);
stderr_PRS = [stderr_c01_PRS_1;stderr_c02_PRS_1;stderr_c03_PRS_1;stderr_c04_PRS_1;stderr_c05_PRS_1;stderr_c06_PRS_1;stderr_c07_PRS_1;stderr_c08_PRS_1;stderr_c09_PRS_1;stderr_c10_PRS_1;stderr_c11_PRS_1;stderr_c12_PRS_1];

