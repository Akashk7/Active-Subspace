 % clc;clear;close all;

%% Optimization
t0 = ubx;

func_RBF = @(t_RBF)(srgtsRBFEvaluate((normalizeuv([t_RBF(1:5) t_RBF(7)],lbx_obj,ubx_obj)*W1_obj_llm), srgt_RBF_obj_1));
[xmin_RBF,fval_RBF,exitflag_RBF,output_RBF]=fmincon(func_RBF,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,W1_ls_04_llm,W1_ls_05_llm,W1_ls_06_llm,W1_ls_07_llm,W1_ls_08_llm,W1_ls_09_llm,W1_ls_10_llm,...
srgt_RBF_c01,srgt_RBF_c02,srgt_RBF_c03,srgt_RBF_c04,srgt_RBF_c05,srgt_RBF_c06,srgt_RBF_c07,srgt_RBF_c08,srgt_RBF_c09,srgt_RBF_c10,@srgtsRBFEvaluate));

func_KRG = @(t_KRG)(srgtsKRGEvaluate((normalizeuv([t_KRG(1:5) t_KRG(7)],lbx_obj,ubx_obj)*W1_obj_llm), srgt_KRG_obj_1));
[xmin_KRG,fval_KRG,exitflag_KRG,output_KRG]=fmincon(func_KRG,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,W1_ls_04_llm,W1_ls_05_llm,W1_ls_06_llm,W1_ls_07_llm,W1_ls_08_llm,W1_ls_09_llm,W1_ls_10_llm,...
srgt_KRG_c01,srgt_KRG_c02,srgt_KRG_c03,srgt_KRG_c04,srgt_KRG_c05,srgt_KRG_c06,srgt_KRG_c07,srgt_KRG_c08,srgt_KRG_c09,srgt_KRG_c10,@srgtsKRGEvaluate));

func_PRS = @(t_PRS)(srgtsPRSEvaluate((normalizeuv([t_PRS(1:5) t_PRS(7)],lbx_obj,ubx_obj)*W1_obj_llm), srgt_PRS_obj_1));
[xmin_PRS,fval_PRS,exitflag_PRS,output_PRS]=fmincon(func_PRS,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,W1_ls_04_llm,W1_ls_05_llm,W1_ls_06_llm,W1_ls_07_llm,W1_ls_08_llm,W1_ls_09_llm,W1_ls_10_llm,...
srgt_PRS_c01,srgt_PRS_c02,srgt_PRS_c03,srgt_PRS_c04,srgt_PRS_c05,srgt_PRS_c06,srgt_PRS_c07,srgt_PRS_c08,srgt_PRS_c09,srgt_PRS_c10,@srgtsPRSEvaluate));

func_WAS = @(t_WAS)(srgtsWASEvaluate((normalizeuv([t_WAS(1:5) t_WAS(7)],lbx_obj,ubx_obj)*W1_obj_llm), srgt_WAS_obj_1));
[xmin_WAS,fval_WAS,exitflag_WAS,output_WAS]=fmincon(func_WAS,t0,[],[],[],[],lbx,ubx,@(x)nonlcon_active_subspace_c(x,lbx,ubx,...
W1_ls_01_llm,W1_ls_02_llm,W1_ls_03_llm,W1_ls_04_llm,W1_ls_05_llm,W1_ls_06_llm,W1_ls_07_llm,W1_ls_08_llm,W1_ls_09_llm,W1_ls_10_llm,...
srgt_WAS_c01,srgt_WAS_c02,srgt_WAS_c03,srgt_WAS_c04,srgt_WAS_c05,srgt_WAS_c06,srgt_WAS_c07,srgt_WAS_c08,srgt_WAS_c09,srgt_WAS_c10,@srgtsWASEvaluate));

%% Validation

%bounds
lbx_c01=[lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c02=[lbx(1) lbx(2) lbx(3)];
lbx_c03=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(7)];
lbx_c04=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c05=[lbx(2) lbx(3) lbx(7)];
lbx_c06=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c07=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c08=[lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c09=[lbx(1) lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c10=[lbx(3) lbx(5) lbx(6) lbx(7)];

ubx_c01=[ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c02=[ubx(1) ubx(2) ubx(3)];
ubx_c03=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(7)];
ubx_c04=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c05=[ubx(2) ubx(3) ubx(7)];
ubx_c06=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c07=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c08=[ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c09=[ubx(1) ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c10=[ubx(3) ubx(5) ubx(6) ubx(7)];

%RBF
xmin_RBF_obj = normalizeuv([xmin_RBF(1:5) xmin_RBF(7)],lbx_obj, ubx_obj);
xmin_RBF_c01 = normalizeuv([xmin_RBF(2:4) xmin_RBF(6)], lbx_c01, ubx_c01);
xmin_RBF_c02 = normalizeuv(xmin_RBF(1:3), lbx_c02, ubx_c02);
xmin_RBF_c03 = normalizeuv([xmin_RBF(1:3) xmin_RBF(5) xmin_RBF(7)], lbx_c03, ubx_c03);
xmin_RBF_c04 = normalizeuv([xmin_RBF(1:3) xmin_RBF(5:7)], lbx_c04, ubx_c04);
xmin_RBF_c05 = normalizeuv([xmin_RBF(2:3) xmin_RBF(7)], lbx_c05, ubx_c05);
xmin_RBF_c06 = normalizeuv([xmin_RBF(1:3) xmin_RBF(5:7)], lbx_c06, ubx_c06);
xmin_RBF_c07 = normalizeuv([xmin_RBF(1:3) xmin_RBF(5:7)], lbx_c07, ubx_c07);
xmin_RBF_c08 = normalizeuv([xmin_RBF(2:4) xmin_RBF(6)], lbx_c08, ubx_c08);
xmin_RBF_c09 = normalizeuv([xmin_RBF(1:4) xmin_RBF(6)], lbx_c09, ubx_c09);
xmin_RBF_c10 = normalizeuv([xmin_RBF(3) xmin_RBF(5:7)], lbx_c10, ubx_c10);

mm_obj_RBF = srgtsRBFEvaluate(xmin_RBF_obj*W1_obj_llm, srgt_RBF_obj_1);
mm_psf_RBF_c01 = srgtsRBFEvaluate(xmin_RBF_c01*W1_ls_01_llm, srgt_RBF_c01);
mm_psf_RBF_c02 = srgtsRBFEvaluate(xmin_RBF_c02*W1_ls_02_llm, srgt_RBF_c02);
mm_psf_RBF_c03 = srgtsRBFEvaluate(xmin_RBF_c03*W1_ls_03_llm, srgt_RBF_c03);
mm_psf_RBF_c04 = srgtsRBFEvaluate(xmin_RBF_c04*W1_ls_04_llm, srgt_RBF_c04);
mm_psf_RBF_c05 = srgtsRBFEvaluate(xmin_RBF_c05*W1_ls_05_llm, srgt_RBF_c05);
mm_psf_RBF_c06 = srgtsRBFEvaluate(xmin_RBF_c06*W1_ls_06_llm, srgt_RBF_c06);
mm_psf_RBF_c07 = srgtsRBFEvaluate(xmin_RBF_c07*W1_ls_07_llm, srgt_RBF_c07);
mm_psf_RBF_c08 = srgtsRBFEvaluate(xmin_RBF_c08*W1_ls_08_llm, srgt_RBF_c08);
mm_psf_RBF_c09 = srgtsRBFEvaluate(xmin_RBF_c09*W1_ls_09_llm, srgt_RBF_c09);
mm_psf_RBF_c10 = srgtsRBFEvaluate(xmin_RBF_c10*W1_ls_10_llm, srgt_RBF_c10);

act_obj_RBF = Weightc([xmin_RBF(1:5) xmin_RBF(7)]);
act_psf_RBF_c01 = mcspsfconstraint1c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c02 = mcspsfconstraint2c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c03 = mcspsfconstraint3c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c04 = mcspsfconstraint4c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c05 = mcspsfconstraint5c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c06 = mcspsfconstraint6c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c07 = mcspsfconstraint7c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c08 = mcspsfconstraint8c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c09 = mcspsfconstraint9c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_RBF_c10 = mcspsfconstraint10c(xmin_RBF,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%KRG
xmin_KRG_obj = normalizeuv([xmin_KRG(1:5) xmin_KRG(7)],lbx_obj, ubx_obj);
xmin_KRG_c01 = normalizeuv([xmin_KRG(2:4) xmin_KRG(6)], lbx_c01, ubx_c01);
xmin_KRG_c02 = normalizeuv(xmin_KRG(1:3), lbx_c02, ubx_c02);
xmin_KRG_c03 = normalizeuv([xmin_KRG(1:3) xmin_KRG(5) xmin_KRG(7)], lbx_c03, ubx_c03);
xmin_KRG_c04 = normalizeuv([xmin_KRG(1:3) xmin_KRG(5:7)], lbx_c04, ubx_c04);
xmin_KRG_c05 = normalizeuv([xmin_KRG(2:3) xmin_KRG(7)], lbx_c05, ubx_c05);
xmin_KRG_c06 = normalizeuv([xmin_KRG(1:3) xmin_KRG(5:7)], lbx_c06, ubx_c06);
xmin_KRG_c07 = normalizeuv([xmin_KRG(1:3) xmin_KRG(5:7)], lbx_c07, ubx_c07);
xmin_KRG_c08 = normalizeuv([xmin_KRG(2:4) xmin_KRG(6)], lbx_c08, ubx_c08);
xmin_KRG_c09 = normalizeuv([xmin_KRG(1:4) xmin_KRG(6)], lbx_c09, ubx_c09);
xmin_KRG_c10 = normalizeuv([xmin_KRG(3) xmin_KRG(5:7)], lbx_c10, ubx_c10);

mm_obj_KRG = srgtsKRGEvaluate(xmin_KRG_obj*W1_obj_llm, srgt_KRG_obj_1);
mm_psf_KRG_c01 = srgtsKRGEvaluate(xmin_KRG_c01*W1_ls_01_llm, srgt_KRG_c01);
mm_psf_KRG_c02 = srgtsKRGEvaluate(xmin_KRG_c02*W1_ls_02_llm, srgt_KRG_c02);
mm_psf_KRG_c03 = srgtsKRGEvaluate(xmin_KRG_c03*W1_ls_03_llm, srgt_KRG_c03);
mm_psf_KRG_c04 = srgtsKRGEvaluate(xmin_KRG_c04*W1_ls_04_llm, srgt_KRG_c04);
mm_psf_KRG_c05 = srgtsKRGEvaluate(xmin_KRG_c05*W1_ls_05_llm, srgt_KRG_c05);
mm_psf_KRG_c06 = srgtsKRGEvaluate(xmin_KRG_c06*W1_ls_06_llm, srgt_KRG_c06);
mm_psf_KRG_c07 = srgtsKRGEvaluate(xmin_KRG_c07*W1_ls_07_llm, srgt_KRG_c07);
mm_psf_KRG_c08 = srgtsKRGEvaluate(xmin_KRG_c08*W1_ls_08_llm, srgt_KRG_c08);
mm_psf_KRG_c09 = srgtsKRGEvaluate(xmin_KRG_c09*W1_ls_09_llm, srgt_KRG_c09);
mm_psf_KRG_c10 = srgtsKRGEvaluate(xmin_KRG_c10*W1_ls_10_llm, srgt_KRG_c10);

act_obj_KRG = Weightc([xmin_KRG(1:5) xmin_KRG(7)]);
act_psf_KRG_c01 = mcspsfconstraint1c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c02 = mcspsfconstraint2c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c03 = mcspsfconstraint3c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c04 = mcspsfconstraint4c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c05 = mcspsfconstraint5c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c06 = mcspsfconstraint6c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c07 = mcspsfconstraint7c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c08 = mcspsfconstraint8c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c09 = mcspsfconstraint9c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_KRG_c10 = mcspsfconstraint10c(xmin_KRG,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
                                
%PRS
%PRS
xmin_PRS_obj = normalizeuv([xmin_PRS(1:5) xmin_PRS(7)],lbx_obj, ubx_obj);
xmin_PRS_c01 = normalizeuv([xmin_PRS(2:4) xmin_PRS(6)], lbx_c01, ubx_c01);
xmin_PRS_c02 = normalizeuv(xmin_PRS(1:3), lbx_c02, ubx_c02);
xmin_PRS_c03 = normalizeuv([xmin_PRS(1:3) xmin_PRS(5) xmin_PRS(7)], lbx_c03, ubx_c03);
xmin_PRS_c04 = normalizeuv([xmin_PRS(1:3) xmin_PRS(5:7)], lbx_c04, ubx_c04);
xmin_PRS_c05 = normalizeuv([xmin_PRS(2:3) xmin_PRS(7)], lbx_c05, ubx_c05);
xmin_PRS_c06 = normalizeuv([xmin_PRS(1:3) xmin_PRS(5:7)], lbx_c06, ubx_c06);
xmin_PRS_c07 = normalizeuv([xmin_PRS(1:3) xmin_PRS(5:7)], lbx_c07, ubx_c07);
xmin_PRS_c08 = normalizeuv([xmin_PRS(2:4) xmin_PRS(6)], lbx_c08, ubx_c08);
xmin_PRS_c09 = normalizeuv([xmin_PRS(1:4) xmin_PRS(6)], lbx_c09, ubx_c09);
xmin_PRS_c10 = normalizeuv([xmin_PRS(3) xmin_PRS(5:7)], lbx_c10, ubx_c10);

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

act_obj_PRS = Weightc([xmin_PRS(1:5) xmin_PRS(7)]);
act_psf_PRS_c01 = mcspsfconstraint1c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c02 = mcspsfconstraint2c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c03 = mcspsfconstraint3c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c04 = mcspsfconstraint4c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c05 = mcspsfconstraint5c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c06 = mcspsfconstraint6c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c07 = mcspsfconstraint7c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c08 = mcspsfconstraint8c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c09 = mcspsfconstraint9c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_PRS_c10 = mcspsfconstraint10c(xmin_PRS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%WAS
%WAS
xmin_WAS_obj = normalizeuv([xmin_WAS(1:5) xmin_WAS(7)],lbx_obj, ubx_obj);
xmin_WAS_c01 = normalizeuv([xmin_WAS(2:4) xmin_WAS(6)], lbx_c01, ubx_c01);
xmin_WAS_c02 = normalizeuv(xmin_WAS(1:3), lbx_c02, ubx_c02);
xmin_WAS_c03 = normalizeuv([xmin_WAS(1:3) xmin_WAS(5) xmin_WAS(7)], lbx_c03, ubx_c03);
xmin_WAS_c04 = normalizeuv([xmin_WAS(1:3) xmin_WAS(5:7)], lbx_c04, ubx_c04);
xmin_WAS_c05 = normalizeuv([xmin_WAS(2:3) xmin_WAS(7)], lbx_c05, ubx_c05);
xmin_WAS_c06 = normalizeuv([xmin_WAS(1:3) xmin_WAS(5:7)], lbx_c06, ubx_c06);
xmin_WAS_c07 = normalizeuv([xmin_WAS(1:3) xmin_WAS(5:7)], lbx_c07, ubx_c07);
xmin_WAS_c08 = normalizeuv([xmin_WAS(2:4) xmin_WAS(6)], lbx_c08, ubx_c08);
xmin_WAS_c09 = normalizeuv([xmin_WAS(1:4) xmin_WAS(6)], lbx_c09, ubx_c09);
xmin_WAS_c10 = normalizeuv([xmin_WAS(3) xmin_WAS(5:7)], lbx_c10, ubx_c10);

mm_obj_WAS = srgtsWASEvaluate(xmin_WAS_obj*W1_obj_llm, srgt_WAS_obj_1);
mm_psf_WAS_c01 = srgtsWASEvaluate(xmin_WAS_c01*W1_ls_01_llm, srgt_WAS_c01);
mm_psf_WAS_c02 = srgtsWASEvaluate(xmin_WAS_c02*W1_ls_02_llm, srgt_WAS_c02);
mm_psf_WAS_c03 = srgtsWASEvaluate(xmin_WAS_c03*W1_ls_03_llm, srgt_WAS_c03);
mm_psf_WAS_c04 = srgtsWASEvaluate(xmin_WAS_c04*W1_ls_04_llm, srgt_WAS_c04);
mm_psf_WAS_c05 = srgtsWASEvaluate(xmin_WAS_c05*W1_ls_05_llm, srgt_WAS_c05);
mm_psf_WAS_c06 = srgtsWASEvaluate(xmin_WAS_c06*W1_ls_06_llm, srgt_WAS_c06);
mm_psf_WAS_c07 = srgtsWASEvaluate(xmin_WAS_c07*W1_ls_07_llm, srgt_WAS_c07);
mm_psf_WAS_c08 = srgtsWASEvaluate(xmin_WAS_c08*W1_ls_08_llm, srgt_WAS_c08);
mm_psf_WAS_c09 = srgtsWASEvaluate(xmin_WAS_c09*W1_ls_09_llm, srgt_WAS_c09);
mm_psf_WAS_c10 = srgtsWASEvaluate(xmin_WAS_c10*W1_ls_10_llm, srgt_WAS_c10);

act_obj_WAS = Weightc([xmin_WAS(1:5) xmin_WAS(7)]);
act_psf_WAS_c01 = mcspsfconstraint1c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c02 = mcspsfconstraint2c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c03 = mcspsfconstraint3c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c04 = mcspsfconstraint4c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c05 = mcspsfconstraint5c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c06 = mcspsfconstraint6c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c07 = mcspsfconstraint7c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c08 = mcspsfconstraint8c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c09 = mcspsfconstraint9c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
act_psf_WAS_c10 = mcspsfconstraint10c(xmin_WAS,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%% Prediction Variance
predvar_c01_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c01*W1_ls_01_llm, srgt_KRG_c01));
predvar_c02_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c02*W1_ls_02_llm, srgt_KRG_c02));
predvar_c03_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c03*W1_ls_03_llm, srgt_KRG_c03));
predvar_c04_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c04*W1_ls_04_llm, srgt_KRG_c04));
predvar_c05_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c05*W1_ls_05_llm, srgt_KRG_c05));
predvar_c06_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c06*W1_ls_06_llm, srgt_KRG_c06));
predvar_c07_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c07*W1_ls_07_llm, srgt_KRG_c07));
predvar_c08_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c08*W1_ls_08_llm, srgt_KRG_c08));
predvar_c09_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c09*W1_ls_09_llm, srgt_KRG_c09));
predvar_c10_KRG = (srgtsKRGPredictionVariance(xmin_KRG_c10*W1_ls_10_llm, srgt_KRG_c10));

stderr_c01_KRG_1 = sqrt(predvar_c01_KRG);
stderr_c02_KRG_1 = sqrt(predvar_c02_KRG);
stderr_c03_KRG_1 = sqrt(predvar_c03_KRG);
stderr_c04_KRG_1 = sqrt(predvar_c04_KRG);
stderr_c05_KRG_1 = sqrt(predvar_c05_KRG);
stderr_c06_KRG_1 = sqrt(predvar_c06_KRG);
stderr_c07_KRG_1 = sqrt(predvar_c07_KRG);
stderr_c08_KRG_1 = sqrt(predvar_c08_KRG);
stderr_c09_KRG_1 = sqrt(predvar_c09_KRG);
stderr_c10_KRG_1 = sqrt(predvar_c10_KRG);

stderr_KRG = [stderr_c01_KRG_1;stderr_c02_KRG_1;stderr_c03_KRG_1;stderr_c04_KRG_1;stderr_c05_KRG_1;stderr_c06_KRG_1;stderr_c07_KRG_1;stderr_c08_KRG_1;stderr_c09_KRG_1;stderr_c10_KRG_1];

