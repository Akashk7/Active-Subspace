 % clc;clear;close all;

%% Optimization
t0 = 6*ones(1,n_var);

func_RBF = @(t_RBF)(srgtsRBFEvaluate((normalizeuv(t_RBF,lbx_obj,ubx_obj)*W1_obj_llm), srgt_RBF_obj_1));
[xmin_RBF,fval_RBF,exitflag_RBF,output_RBF]=fmincon(func_RBF,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,...
srgt_RBF_c01,srgt_RBF_c02,srgt_RBF_c03,@srgtsRBFEvaluate));

func_KRG = @(t_KRG)(srgtsKRGEvaluate((normalizeuv(t_KRG,lbx_obj,ubx_obj)*W1_obj_llm), srgt_KRG_obj_1));
[xmin_KRG,fval_KRG,exitflag_KRG,output_KRG]=fmincon(func_KRG,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,...
srgt_KRG_c01,srgt_KRG_c02,srgt_KRG_c03,@srgtsKRGEvaluate));

func_PRS = @(t_PRS)(srgtsPRSEvaluate((normalizeuv(t_PRS,lbx_obj,ubx_obj)*W1_obj_llm), srgt_PRS_obj_1));
[xmin_PRS,fval_PRS,exitflag_PRS,output_PRS]=fmincon(func_PRS,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,...
srgt_PRS_c01,srgt_PRS_c02,srgt_PRS_c03,@srgtsPRSEvaluate));

func_WAS = @(t_WAS)(srgtsWASEvaluate((normalizeuv(t_WAS,lbx_obj,ubx_obj)*W1_obj_llm), srgt_WAS_obj_1));
[xmin_WAS,fval_WAS,exitflag_WAS,output_WAS]=fmincon(func_WAS,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,...
srgt_WAS_c01,srgt_WAS_c02,srgt_WAS_c03,@srgtsWASEvaluate));

%% Validation

%bounds
lbx_c01 = lbx;
lbx_c02 = lbx;
lbx_c03 = lbx;

ubx_c01 = ubx;
ubx_c02 = ubx;
ubx_c03 = ubx;

%RBF
xmin_RBF_obj = normalizeuv(xmin_RBF, lbx_obj, ubx_obj);
xmin_RBF_c01 = normalizeuv(xmin_RBF, lbx_c01, ubx_c01);
xmin_RBF_c02 = normalizeuv(xmin_RBF, lbx_c02, ubx_c02);
xmin_RBF_c03 = normalizeuv(xmin_RBF, lbx_c03, ubx_c03);

mm_obj_RBF = srgtsRBFEvaluate(xmin_RBF_obj*W1_obj_llm, srgt_RBF_obj_1);
mm_psf_RBF_c01 = srgtsRBFEvaluate(xmin_RBF_c01*W1_ls_01_llm, srgt_RBF_c01);
mm_psf_RBF_c02 = srgtsRBFEvaluate(xmin_RBF_c02*W1_ls_02_llm, srgt_RBF_c02);
mm_psf_RBF_c03 = srgtsRBFEvaluate(xmin_RBF_c03*W1_ls_03_llm, srgt_RBF_c03);

act_obj_RBF = Weightc(xmin_RBF);
act_psf_RBF_c01 = mcspsfconstraint1c(xmin_RBF,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_RBF_c02 = mcspsfconstraint2c(xmin_RBF,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_RBF_c03 = mcspsfconstraint3c(xmin_RBF,lbx,ubx,sd,Nmcs,Pftarget);

%KRG
xmin_KRG_obj = normalizeuv(xmin_KRG, lbx_obj, ubx_obj);
xmin_KRG_c01 = normalizeuv(xmin_KRG, lbx_c01, ubx_c01);
xmin_KRG_c02 = normalizeuv(xmin_KRG, lbx_c02, ubx_c02);
xmin_KRG_c03 = normalizeuv(xmin_KRG, lbx_c03, ubx_c03);

mm_obj_KRG = srgtsKRGEvaluate(xmin_KRG_obj*W1_obj_llm, srgt_KRG_obj_1);
mm_psf_KRG_c01 = srgtsKRGEvaluate(xmin_KRG_c01*W1_ls_01_llm, srgt_KRG_c01);
mm_psf_KRG_c02 = srgtsKRGEvaluate(xmin_KRG_c02*W1_ls_02_llm, srgt_KRG_c02);
mm_psf_KRG_c03 = srgtsKRGEvaluate(xmin_KRG_c03*W1_ls_03_llm, srgt_KRG_c03);

act_obj_KRG = Weightc(xmin_KRG);
act_psf_KRG_c01 = mcspsfconstraint1c(xmin_KRG,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_KRG_c02 = mcspsfconstraint2c(xmin_KRG,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_KRG_c03 = mcspsfconstraint3c(xmin_KRG,lbx,ubx,sd,Nmcs,Pftarget);
                                
%PRS
%PRS
xmin_PRS_obj = normalizeuv(xmin_PRS, lbx_obj, ubx_obj);
xmin_PRS_c01 = normalizeuv(xmin_PRS, lbx_c01, ubx_c01);
xmin_PRS_c02 = normalizeuv(xmin_PRS, lbx_c02, ubx_c02);
xmin_PRS_c03 = normalizeuv(xmin_PRS, lbx_c03, ubx_c03);

mm_obj_PRS = srgtsPRSEvaluate(xmin_PRS_obj*W1_obj_llm, srgt_PRS_obj_1);
mm_psf_PRS_c01 = srgtsPRSEvaluate(xmin_PRS_c01*W1_ls_01_llm, srgt_PRS_c01);
mm_psf_PRS_c02 = srgtsPRSEvaluate(xmin_PRS_c02*W1_ls_02_llm, srgt_PRS_c02);
mm_psf_PRS_c03 = srgtsPRSEvaluate(xmin_PRS_c03*W1_ls_03_llm, srgt_PRS_c03);

act_obj_PRS = Weightc(xmin_PRS);
act_psf_PRS_c01 = mcspsfconstraint1c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c02 = mcspsfconstraint2c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_PRS_c03 = mcspsfconstraint3c(xmin_PRS,lbx,ubx,sd,Nmcs,Pftarget);

%WAS
%WAS
xmin_WAS_obj = normalizeuv(xmin_WAS, lbx_obj, ubx_obj);
xmin_WAS_c01 = normalizeuv(xmin_WAS, lbx_c01, ubx_c01);
xmin_WAS_c02 = normalizeuv(xmin_WAS, lbx_c02, ubx_c02);
xmin_WAS_c03 = normalizeuv(xmin_WAS, lbx_c03, ubx_c03);

mm_obj_WAS = srgtsWASEvaluate(xmin_WAS_obj*W1_obj_llm, srgt_WAS_obj_1);
mm_psf_WAS_c01 = srgtsWASEvaluate(xmin_WAS_c01*W1_ls_01_llm, srgt_WAS_c01);
mm_psf_WAS_c02 = srgtsWASEvaluate(xmin_WAS_c02*W1_ls_02_llm, srgt_WAS_c02);
mm_psf_WAS_c03 = srgtsWASEvaluate(xmin_WAS_c03*W1_ls_03_llm, srgt_WAS_c03);

act_obj_WAS = Weightc(xmin_WAS);
act_psf_WAS_c01 = mcspsfconstraint1c(xmin_WAS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_WAS_c02 = mcspsfconstraint2c(xmin_WAS,lbx,ubx,sd,Nmcs,Pftarget);
act_psf_WAS_c03 = mcspsfconstraint3c(xmin_WAS,lbx,ubx,sd,Nmcs,Pftarget);

%% Tabulate Results

xmin = [xmin_RBF;xmin_KRG;xmin_PRS;xmin_WAS];
fval = [fval_RBF;fval_KRG;fval_PRS;fval_WAS];
exitflag = [exitflag_RBF;exitflag_KRG;exitflag_PRS;exitflag_WAS];
output = [output_RBF, output_KRG, output_PRS, output_WAS];
metamodel_values = [mm_obj_RBF, mm_obj_KRG, mm_obj_PRS, mm_obj_WAS; mm_psf_RBF_c01, mm_psf_KRG_c01, mm_psf_PRS_c01, mm_psf_WAS_c01; mm_psf_RBF_c02, mm_psf_KRG_c02, mm_psf_PRS_c02, mm_psf_WAS_c02; mm_psf_RBF_c03, mm_psf_KRG_c03, mm_psf_PRS_c03, mm_psf_WAS_c03];
actual_values = [exitflag_RBF,exitflag_KRG,exitflag_PRS,exitflag_WAS;act_obj_RBF, act_obj_KRG, act_obj_PRS, act_obj_WAS; act_psf_RBF_c01, act_psf_KRG_c01, act_psf_PRS_c01, act_psf_WAS_c01; act_psf_RBF_c02, act_psf_KRG_c02, act_psf_PRS_c02, act_psf_WAS_c02; act_psf_RBF_c03, act_psf_KRG_c03, act_psf_PRS_c03, act_psf_WAS_c03];
r2_predicted = [R2_pred_RBF_obj_1, R2_pred_KRG_obj_1, R2_pred_PRS_obj_1, R2_pred_WAS_obj_1; R2_pred_RBF_c01, R2_pred_KRG_c01, R2_pred_PRS_c01, R2_pred_WAS_c01; R2_pred_RBF_c02, R2_pred_KRG_c02, R2_pred_PRS_c02, R2_pred_WAS_c02; R2_pred_RBF_c03, R2_pred_KRG_c03, R2_pred_PRS_c03, R2_pred_WAS_c03];
press_rmse = [PRESSRMS_RBF_obj_1, PRESSRMS_KRG_obj_1, PRESSRMS_PRS_obj_1, PRESSRMS_WAS_obj_1; PRESSRMS_RBF_c01, PRESSRMS_KRG_c01, PRESSRMS_PRS_c01, PRESSRMS_WAS_c01; PRESSRMS_RBF_c02, PRESSRMS_KRG_c02, PRESSRMS_PRS_c02, PRESSRMS_WAS_c02; PRESSRMS_RBF_c03, PRESSRMS_KRG_c03, PRESSRMS_PRS_c03, PRESSRMS_WAS_c03];
excel_table=[mm_obj_RBF,act_obj_RBF,mm_obj_KRG,act_obj_KRG,mm_obj_PRS,act_obj_PRS,mm_obj_WAS,act_obj_WAS; mm_psf_RBF_c01, act_psf_RBF_c01, mm_psf_KRG_c01, act_psf_KRG_c01, mm_psf_PRS_c01, act_psf_PRS_c01, mm_psf_WAS_c01, act_psf_WAS_c01; mm_psf_RBF_c02, act_psf_RBF_c02, mm_psf_KRG_c02, act_psf_KRG_c02, mm_psf_PRS_c02, act_psf_PRS_c02, mm_psf_WAS_c02, act_psf_WAS_c02; mm_psf_RBF_c03, act_psf_RBF_c03, mm_psf_KRG_c03, act_psf_KRG_c03, mm_psf_PRS_c03, act_psf_PRS_c03, mm_psf_WAS_c03, act_psf_WAS_c03];
std__inactive = [std_obj_RBF std_obj_KRG std_obj_PRS std_obj_WAS; std_c01_RBF std_c01_KRG std_c01_PRS std_c01_WAS;std_c02_RBF std_c02_KRG std_c02_PRS std_c02_WAS; std_c03_RBF std_c03_KRG std_c03_PRS std_c03_WAS];    
percentage_variation = [pv_ls_01(1);pv_ls_02(1);pv_ls_03(1)];

excel_table_2_RBF = [act_obj_RBF, mm_obj_RBF, R2_pred_RBF_obj_1, PRESSRMS_RBF_obj_1, std_obj_RBF;act_psf_RBF_c01, mm_psf_RBF_c01, R2_pred_RBF_c01, PRESSRMS_RBF_c01, std_c01_RBF;act_psf_RBF_c02, mm_psf_RBF_c02, R2_pred_RBF_c02, PRESSRMS_RBF_c02, std_c02_RBF;act_psf_RBF_c03, mm_psf_RBF_c03, R2_pred_RBF_c03, PRESSRMS_RBF_c03, std_c03_RBF];
excel_table_2_KRG = [act_obj_KRG, mm_obj_KRG, R2_pred_KRG_obj_1, PRESSRMS_KRG_obj_1, std_obj_KRG;act_psf_KRG_c01, mm_psf_KRG_c01, R2_pred_KRG_c01, PRESSRMS_KRG_c01, std_c01_KRG;act_psf_KRG_c02, mm_psf_KRG_c02, R2_pred_KRG_c02, PRESSRMS_KRG_c02, std_c02_KRG;act_psf_KRG_c03, mm_psf_KRG_c03, R2_pred_KRG_c03, PRESSRMS_KRG_c03, std_c03_KRG];
excel_table_2_PRS = [act_obj_PRS, mm_obj_PRS, R2_pred_PRS_obj_1, PRESSRMS_PRS_obj_1, std_obj_PRS;act_psf_PRS_c01, mm_psf_PRS_c01, R2_pred_PRS_c01, PRESSRMS_PRS_c01, std_c01_PRS;act_psf_PRS_c02, mm_psf_PRS_c02, R2_pred_PRS_c02, PRESSRMS_PRS_c02, std_c02_PRS;act_psf_PRS_c03, mm_psf_PRS_c03, R2_pred_PRS_c03, PRESSRMS_PRS_c03, std_c03_PRS];
excel_table_2_WAS = [act_obj_WAS, mm_obj_WAS, R2_pred_WAS_obj_1, PRESSRMS_WAS_obj_1, std_obj_WAS;act_psf_WAS_c01, mm_psf_WAS_c01, R2_pred_WAS_c01, PRESSRMS_WAS_c01, std_c01_WAS;act_psf_WAS_c02, mm_psf_WAS_c02, R2_pred_WAS_c02, PRESSRMS_WAS_c02, std_c02_WAS;act_psf_WAS_c03, mm_psf_WAS_c03, R2_pred_WAS_c03, PRESSRMS_WAS_c03, std_c03_WAS];


%% Somplots
% h_ls_01 = isomplots([scalenewDoe1_c01(:,1:4) y_ls_01],'Limit State 1');
% psf_01_hat = srgtsRBFEvaluate(t1_ls_01, srgt_RBF_c01); 
% h_psf_01 = isomplots([scalenewDoe1_c01(:,1:4) psf_01_hat],'PSF 1');
% x_psf_01 = [zeros(npoints_c01,1) scalenewDoe1_c01(:,1:3) zeros(npoints_c01,1) scalenewDoe1_c01(:,4) zeros(npoints_c01,1)];
% psf_01_actual = mcspsfconstraint1c1(x_psf_01,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
% h_psf_01_actual = isomplots([scalenewDoe1_c01(:,1:4) psf_01_actual],'PSF 1 Actual');
