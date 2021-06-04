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
%std__inactive = [std_c01_RBF_enrico std_c01_KRG_enrico std_c01_PRS_enrico std_c01_WAS_enrico;std_c02_RBF_enrico std_c02_KRG_enrico std_c02_PRS_enrico std_c02_WAS_enrico; std_c03_RBF_enrico std_c03_KRG_enrico std_c03_PRS_enrico std_c03_WAS_enrico];    
percentage_variation = [pv_ls_01(1);pv_ls_02(1);pv_ls_03(1)];

%excel_table_2_RBF = [act_psf_RBF_c01, mm_psf_RBF_c01, R2_pred_RBF_c01, PRESSRMS_RBF_c01, std_c01_RBF_enrico;act_psf_RBF_c02, mm_psf_RBF_c02, R2_pred_RBF_c02, PRESSRMS_RBF_c02, std_c02_RBF_enrico;act_psf_RBF_c03, mm_psf_RBF_c03, R2_pred_RBF_c03, PRESSRMS_RBF_c03, std_c03_RBF_enrico];
excel_table_2_KRG = [act_psf_KRG_c01, mm_psf_KRG_c01, R2_pred_KRG_c01, PRESSRMS_KRG_c01, std_c01_KRG_enrico;act_psf_KRG_c02, mm_psf_KRG_c02, R2_pred_KRG_c02, PRESSRMS_KRG_c02, std_c02_KRG_enrico;act_psf_KRG_c03, mm_psf_KRG_c03, R2_pred_KRG_c03, PRESSRMS_KRG_c03, std_c03_KRG_enrico];
% excel_table_2_PRS = [act_psf_PRS_c01, mm_psf_PRS_c01, R2_pred_PRS_c01, PRESSRMS_PRS_c01, std_c01_PRS_enrico;act_psf_PRS_c02, mm_psf_PRS_c02, R2_pred_PRS_c02, PRESSRMS_PRS_c02, std_c02_PRS_enrico;act_psf_PRS_c03, mm_psf_PRS_c03, R2_pred_PRS_c03, PRESSRMS_PRS_c03, std_c03_PRS_enrico];
% excel_table_2_WAS = [act_psf_WAS_c01, mm_psf_WAS_c01, R2_pred_WAS_c01, PRESSRMS_WAS_c01, std_c01_WAS_enrico;act_psf_WAS_c02, mm_psf_WAS_c02, R2_pred_WAS_c02, PRESSRMS_WAS_c02, std_c02_WAS_enrico;act_psf_WAS_c03, mm_psf_WAS_c03, R2_pred_WAS_c03, PRESSRMS_WAS_c03, std_c03_WAS_enrico];
