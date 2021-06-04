%% Tabulate Results

xmin = [xmin_RBF;xmin_KRG;xmin_PRS;xmin_WAS];
fval = [fval_RBF;fval_KRG;fval_PRS;fval_WAS];
exitflag = [exitflag_RBF;exitflag_KRG;exitflag_PRS;exitflag_WAS];
output = [output_RBF, output_KRG, output_PRS, output_WAS];
metamodel_values = [mm_obj_RBF, mm_obj_KRG, mm_obj_PRS, mm_obj_WAS; mm_psf_RBF_c01, mm_psf_KRG_c01, mm_psf_PRS_c01, mm_psf_WAS_c01; mm_psf_RBF_c02, mm_psf_KRG_c02, mm_psf_PRS_c02, mm_psf_WAS_c02; mm_psf_RBF_c03, mm_psf_KRG_c03, mm_psf_PRS_c03, mm_psf_WAS_c03; mm_psf_RBF_c04, mm_psf_KRG_c04, mm_psf_PRS_c04, mm_psf_WAS_c04; mm_psf_RBF_c05, mm_psf_KRG_c05, mm_psf_PRS_c05, mm_psf_WAS_c05; mm_psf_RBF_c06, mm_psf_KRG_c06, mm_psf_PRS_c06, mm_psf_WAS_c06; mm_psf_RBF_c07, mm_psf_KRG_c07, mm_psf_PRS_c07, mm_psf_WAS_c07; mm_psf_RBF_c08, mm_psf_KRG_c08, mm_psf_PRS_c08, mm_psf_WAS_c08; mm_psf_RBF_c09, mm_psf_KRG_c09, mm_psf_PRS_c09, mm_psf_WAS_c09; mm_psf_RBF_c10, mm_psf_KRG_c10, mm_psf_PRS_c10, mm_psf_WAS_c10];
actual_values = [exitflag_RBF,exitflag_KRG,exitflag_PRS,exitflag_WAS;act_obj_RBF, act_obj_KRG, act_obj_PRS, act_obj_WAS; act_psf_RBF_c01, act_psf_KRG_c01, act_psf_PRS_c01, act_psf_WAS_c01; act_psf_RBF_c02, act_psf_KRG_c02, act_psf_PRS_c02, act_psf_WAS_c02; act_psf_RBF_c03, act_psf_KRG_c03, act_psf_PRS_c03, act_psf_WAS_c03; act_psf_RBF_c04, act_psf_KRG_c04, act_psf_PRS_c04, act_psf_WAS_c04; act_psf_RBF_c05, act_psf_KRG_c05, act_psf_PRS_c05, act_psf_WAS_c05; act_psf_RBF_c06, act_psf_KRG_c06, act_psf_PRS_c06, act_psf_WAS_c06; act_psf_RBF_c07, act_psf_KRG_c07, act_psf_PRS_c07, act_psf_WAS_c07; act_psf_RBF_c08, act_psf_KRG_c08, act_psf_PRS_c08, act_psf_WAS_c08; act_psf_RBF_c09, act_psf_KRG_c09, act_psf_PRS_c09, act_psf_WAS_c09; act_psf_RBF_c10, act_psf_KRG_c10, act_psf_PRS_c10, act_psf_WAS_c10];
r2_predicted = [R2_pred_RBF_obj_1, R2_pred_KRG_obj_1, R2_pred_PRS_obj_1, R2_pred_WAS_obj_1; R2_pred_RBF_c01, R2_pred_KRG_c01, R2_pred_PRS_c01, R2_pred_WAS_c01; R2_pred_RBF_c02, R2_pred_KRG_c02, R2_pred_PRS_c02, R2_pred_WAS_c02; R2_pred_RBF_c03, R2_pred_KRG_c03, R2_pred_PRS_c03, R2_pred_WAS_c03; R2_pred_RBF_c04, R2_pred_KRG_c04, R2_pred_PRS_c04, R2_pred_WAS_c04; R2_pred_RBF_c05, R2_pred_KRG_c05, R2_pred_PRS_c05, R2_pred_WAS_c05; R2_pred_RBF_c06, R2_pred_KRG_c06, R2_pred_PRS_c06, R2_pred_WAS_c06; R2_pred_RBF_c07, R2_pred_KRG_c07, R2_pred_PRS_c07, R2_pred_WAS_c07; R2_pred_RBF_c08, R2_pred_KRG_c08, R2_pred_PRS_c08, R2_pred_WAS_c08; R2_pred_RBF_c09, R2_pred_KRG_c09, R2_pred_PRS_c09, R2_pred_WAS_c09; R2_pred_RBF_c10, R2_pred_KRG_c10, R2_pred_PRS_c10, R2_pred_WAS_c10];
press_rmse = [PRESSRMS_RBF_obj_1, PRESSRMS_KRG_obj_1, PRESSRMS_PRS_obj_1, PRESSRMS_WAS_obj_1; PRESSRMS_RBF_c01, PRESSRMS_KRG_c01, PRESSRMS_PRS_c01, PRESSRMS_WAS_c01; PRESSRMS_RBF_c02, PRESSRMS_KRG_c02, PRESSRMS_PRS_c02, PRESSRMS_WAS_c02; PRESSRMS_RBF_c03, PRESSRMS_KRG_c03, PRESSRMS_PRS_c03, PRESSRMS_WAS_c03; PRESSRMS_RBF_c04, PRESSRMS_KRG_c04, PRESSRMS_PRS_c04, PRESSRMS_WAS_c04; PRESSRMS_RBF_c05, PRESSRMS_KRG_c05, PRESSRMS_PRS_c05, PRESSRMS_WAS_c05; PRESSRMS_RBF_c06, PRESSRMS_KRG_c06, PRESSRMS_PRS_c06, PRESSRMS_WAS_c06; PRESSRMS_RBF_c07, PRESSRMS_KRG_c07, PRESSRMS_PRS_c07, PRESSRMS_WAS_c07; PRESSRMS_RBF_c08, PRESSRMS_KRG_c08, PRESSRMS_PRS_c08, PRESSRMS_WAS_c08; PRESSRMS_RBF_c09, PRESSRMS_KRG_c09, PRESSRMS_PRS_c09, PRESSRMS_WAS_c09; PRESSRMS_RBF_c10, PRESSRMS_KRG_c10, PRESSRMS_PRS_c10, PRESSRMS_WAS_c10];
excel_table=[mm_obj_RBF,act_obj_RBF,mm_obj_KRG,act_obj_KRG,mm_obj_PRS,act_obj_PRS,mm_obj_WAS,act_obj_WAS; mm_psf_RBF_c01, act_psf_RBF_c01, mm_psf_KRG_c01, act_psf_KRG_c01, mm_psf_PRS_c01, act_psf_PRS_c01, mm_psf_WAS_c01, act_psf_WAS_c01; mm_psf_RBF_c02, act_psf_RBF_c02, mm_psf_KRG_c02, act_psf_KRG_c02, mm_psf_PRS_c02, act_psf_PRS_c02, mm_psf_WAS_c02, act_psf_WAS_c02; mm_psf_RBF_c03, act_psf_RBF_c03, mm_psf_KRG_c03, act_psf_KRG_c03, mm_psf_PRS_c03, act_psf_PRS_c03, mm_psf_WAS_c03, act_psf_WAS_c03; mm_psf_RBF_c04, act_psf_RBF_c04, mm_psf_KRG_c04, act_psf_KRG_c04, mm_psf_PRS_c04, act_psf_PRS_c04, mm_psf_WAS_c04, act_psf_WAS_c04; mm_psf_RBF_c05, act_psf_RBF_c05, mm_psf_KRG_c05, act_psf_KRG_c05, mm_psf_PRS_c05, act_psf_PRS_c05, mm_psf_WAS_c05, act_psf_WAS_c05; mm_psf_RBF_c06, act_psf_RBF_c06, mm_psf_KRG_c06, act_psf_KRG_c06, mm_psf_PRS_c06, act_psf_PRS_c06, mm_psf_WAS_c06, act_psf_WAS_c06; mm_psf_RBF_c07, act_psf_RBF_c07, mm_psf_KRG_c07, act_psf_KRG_c07, mm_psf_PRS_c07, act_psf_PRS_c07, mm_psf_WAS_c07, act_psf_WAS_c07; mm_psf_RBF_c08, act_psf_RBF_c08,  mm_psf_KRG_c08, act_psf_KRG_c08,  mm_psf_PRS_c08, act_psf_PRS_c08,  mm_psf_WAS_c08, act_psf_WAS_c08; mm_psf_RBF_c09, act_psf_RBF_c09, mm_psf_KRG_c09, act_psf_KRG_c09, mm_psf_PRS_c09, act_psf_PRS_c09, mm_psf_WAS_c09, act_psf_WAS_c09; mm_psf_RBF_c10, act_psf_RBF_c10, mm_psf_KRG_c10, act_psf_KRG_c10, mm_psf_PRS_c10, act_psf_PRS_c10, mm_psf_WAS_c10, act_psf_WAS_c10];
%std__inactive = [std_obj_RBF std_obj_KRG std_obj_PRS std_obj_WAS; std_c01_RBF std_c01_KRG std_c01_PRS std_c01_WAS;std_c02_RBF std_c02_KRG std_c02_PRS std_c02_WAS; std_c03_RBF std_c03_KRG std_c03_PRS std_c03_WAS; std_c04_RBF std_c04_KRG std_c04_PRS std_c04_WAS; std_c05_RBF std_c05_KRG std_c05_PRS std_c05_WAS; std_c06_RBF std_c06_KRG std_c06_PRS std_c06_WAS; std_c07_RBF std_c07_KRG std_c07_PRS std_c07_WAS; std_c08_RBF std_c08_KRG std_c08_PRS std_c08_WAS; std_c09_RBF std_c09_KRG std_c09_PRS std_c09_WAS;  std_c10_RBF std_c10_KRG std_c10_PRS std_c10_WAS ];     
percentage_variation = [pv_ls_01(1);pv_ls_02(1);pv_ls_03(1);pv_ls_04(1);pv_ls_05(1);pv_ls_06(1);pv_ls_07(1);pv_ls_08(1);pv_ls_09(1);pv_ls_10(1)];
mean_all_variables = [mean(y_obj_1);mean(y_c01);mean(y_c02);mean(y_c03);mean(y_c04);mean(y_c05);mean(y_c06);mean(y_c07);mean(y_c08);mean(y_c09);mean(y_c10)];

% excel_table_2_RBF = [act_obj_RBF, mm_obj_RBF, R2_pred_RBF_obj_1, PRESSRMS_RBF_obj_1, std_obj_RBF;act_psf_RBF_c01, mm_psf_RBF_c01, R2_pred_RBF_c01, PRESSRMS_RBF_c01, std_c01_RBF;act_psf_RBF_c02, mm_psf_RBF_c02, R2_pred_RBF_c02, PRESSRMS_RBF_c02, std_c02_RBF;act_psf_RBF_c03, mm_psf_RBF_c03, R2_pred_RBF_c03, PRESSRMS_RBF_c03, std_c03_RBF;act_psf_RBF_c04, mm_psf_RBF_c04, R2_pred_RBF_c04, PRESSRMS_RBF_c04, std_c04_RBF; act_psf_RBF_c05, mm_psf_RBF_c05, R2_pred_RBF_c05, PRESSRMS_RBF_c05, std_c05_RBF;act_psf_RBF_c06, mm_psf_RBF_c06, R2_pred_RBF_c06, PRESSRMS_RBF_c06, std_c06_RBF;act_psf_RBF_c07, mm_psf_RBF_c07, R2_pred_RBF_c07, PRESSRMS_RBF_c07, std_c07_RBF;act_psf_RBF_c08, mm_psf_RBF_c08, R2_pred_RBF_c08, PRESSRMS_RBF_c08, std_c08_RBF;act_psf_RBF_c09, mm_psf_RBF_c09, R2_pred_RBF_c09, PRESSRMS_RBF_c09, std_c09_RBF;act_psf_RBF_c10, mm_psf_RBF_c10, R2_pred_RBF_c10, PRESSRMS_RBF_c10, std_c10_RBF];
excel_table_2_KRG = [act_obj_KRG, mm_obj_KRG, R2_pred_KRG_obj_1, PRESSRMS_KRG_obj_1, std_obj_KRG;act_psf_KRG_c01, mm_psf_KRG_c01, R2_pred_KRG_c01, PRESSRMS_KRG_c01, std_c01_KRG;act_psf_KRG_c02, mm_psf_KRG_c02, R2_pred_KRG_c02, PRESSRMS_KRG_c02, std_c02_KRG;act_psf_KRG_c03, mm_psf_KRG_c03, R2_pred_KRG_c03, PRESSRMS_KRG_c03, std_c03_KRG;act_psf_KRG_c04, mm_psf_KRG_c04, R2_pred_KRG_c04, PRESSRMS_KRG_c04, std_c04_KRG; act_psf_KRG_c05, mm_psf_KRG_c05, R2_pred_KRG_c05, PRESSRMS_KRG_c05, std_c05_KRG;act_psf_KRG_c06, mm_psf_KRG_c06, R2_pred_KRG_c06, PRESSRMS_KRG_c06, std_c06_KRG;act_psf_KRG_c07, mm_psf_KRG_c07, R2_pred_KRG_c07, PRESSRMS_KRG_c07, std_c07_KRG;act_psf_KRG_c08, mm_psf_KRG_c08, R2_pred_KRG_c08, PRESSRMS_KRG_c08, std_c08_KRG;act_psf_KRG_c09, mm_psf_KRG_c09, R2_pred_KRG_c09, PRESSRMS_KRG_c09, std_c09_KRG;act_psf_KRG_c10, mm_psf_KRG_c10, R2_pred_KRG_c10, PRESSRMS_KRG_c10, std_c10_KRG];
excel_table_2_PRS = [act_obj_PRS, mm_obj_PRS, R2_pred_PRS_obj_1, PRESSRMS_PRS_obj_1, std_obj_PRS;act_psf_PRS_c01, mm_psf_PRS_c01, R2_pred_PRS_c01, PRESSRMS_PRS_c01, std_c01_PRS;act_psf_PRS_c02, mm_psf_PRS_c02, R2_pred_PRS_c02, PRESSRMS_PRS_c02, std_c02_PRS;act_psf_PRS_c03, mm_psf_PRS_c03, R2_pred_PRS_c03, PRESSRMS_PRS_c03, std_c03_PRS;act_psf_PRS_c04, mm_psf_PRS_c04, R2_pred_PRS_c04, PRESSRMS_PRS_c04, std_c04_PRS; act_psf_PRS_c05, mm_psf_PRS_c05, R2_pred_PRS_c05, PRESSRMS_PRS_c05, std_c05_PRS;act_psf_PRS_c06, mm_psf_PRS_c06, R2_pred_PRS_c06, PRESSRMS_PRS_c06, std_c06_PRS;act_psf_PRS_c07, mm_psf_PRS_c07, R2_pred_PRS_c07, PRESSRMS_PRS_c07, std_c07_PRS;act_psf_PRS_c08, mm_psf_PRS_c08, R2_pred_PRS_c08, PRESSRMS_PRS_c08, std_c08_PRS;act_psf_PRS_c09, mm_psf_PRS_c09, R2_pred_PRS_c09, PRESSRMS_PRS_c09, std_c09_PRS;act_psf_PRS_c10, mm_psf_PRS_c10, R2_pred_PRS_c10, PRESSRMS_PRS_c10, std_c10_PRS];
% excel_table_2_WAS = [act_obj_WAS, mm_obj_WAS, R2_pred_WAS_obj_1, PRESSRMS_WAS_obj_1, std_obj_WAS;act_psf_WAS_c01, mm_psf_WAS_c01, R2_pred_WAS_c01, PRESSRMS_WAS_c01, std_c01_WAS;act_psf_WAS_c02, mm_psf_WAS_c02, R2_pred_WAS_c02, PRESSRMS_WAS_c02, std_c02_WAS;act_psf_WAS_c03, mm_psf_WAS_c03, R2_pred_WAS_c03, PRESSRMS_WAS_c03, std_c03_WAS;act_psf_WAS_c04, mm_psf_WAS_c04, R2_pred_WAS_c04, PRESSRMS_WAS_c04, std_c04_WAS; act_psf_WAS_c05, mm_psf_WAS_c05, R2_pred_WAS_c05, PRESSRMS_WAS_c05, std_c05_WAS;act_psf_WAS_c06, mm_psf_WAS_c06, R2_pred_WAS_c06, PRESSRMS_WAS_c06, std_c06_WAS;act_psf_WAS_c07, mm_psf_WAS_c07, R2_pred_WAS_c07, PRESSRMS_WAS_c07, std_c07_WAS;act_psf_WAS_c08, mm_psf_WAS_c08, R2_pred_WAS_c08, PRESSRMS_WAS_c08, std_c08_WAS;act_psf_WAS_c09, mm_psf_WAS_c09, R2_pred_WAS_c09, PRESSRMS_WAS_c09, std_c09_WAS;act_psf_WAS_c10, mm_psf_WAS_c10, R2_pred_WAS_c10, PRESSRMS_WAS_c10, std_c10_WAS];
