%% Tabulate Results

xmin = xmin_PRS;
fval = fval_PRS;
exitflag = exitflag_PRS;
output = output_PRS;
metamodel_values = [mm_obj_PRS; mm_psf_PRS_c01; mm_psf_PRS_c02; mm_psf_PRS_c03; mm_psf_PRS_c04; mm_psf_PRS_c05; mm_psf_PRS_c06; mm_psf_PRS_c07; mm_psf_PRS_c08; mm_psf_PRS_c09; mm_psf_PRS_c10; mm_psf_PRS_c11; mm_psf_PRS_c12];
percentage_variation = [pv_ls_01(1);pv_ls_02(1);pv_ls_03(1);pv_ls_04(1);pv_ls_05(1);pv_ls_06(1);pv_ls_07(1);pv_ls_08(1);pv_ls_09(1);pv_ls_10(1);pv_ls_11(1);pv_ls_12(1)];
mean_all_variables = [mean(y_obj_1);mean(y_c01);mean(y_c02);mean(y_c03);mean(y_c04);mean(y_c05);mean(y_c06);mean(y_c07);mean(y_c08);mean(y_c09);mean(y_c10);mean(y_c11);mean(y_c12)];
excel_table_2_PRS = [act_obj_PRS, mm_obj_PRS, R2_pred_PRS_obj_1, PRESSRMS_PRS_obj_1, std_obj_PRS;act_psf_PRS_c01, mm_psf_PRS_c01, R2_pred_PRS_c01, PRESSRMS_PRS_c01, std_c01_PRS;act_psf_PRS_c02, mm_psf_PRS_c02, R2_pred_PRS_c02, PRESSRMS_PRS_c02, std_c02_PRS;act_psf_PRS_c03, mm_psf_PRS_c03, R2_pred_PRS_c03, PRESSRMS_PRS_c03, std_c03_PRS;act_psf_PRS_c04, mm_psf_PRS_c04, R2_pred_PRS_c04, PRESSRMS_PRS_c04, std_c04_PRS; act_psf_PRS_c05, mm_psf_PRS_c05, R2_pred_PRS_c05, PRESSRMS_PRS_c05, std_c05_PRS;act_psf_PRS_c06, mm_psf_PRS_c06, R2_pred_PRS_c06, PRESSRMS_PRS_c06, std_c06_PRS;act_psf_PRS_c07, mm_psf_PRS_c07, R2_pred_PRS_c07, PRESSRMS_PRS_c07, std_c07_PRS;act_psf_PRS_c08, mm_psf_PRS_c08, R2_pred_PRS_c08, PRESSRMS_PRS_c08, std_c08_PRS;act_psf_PRS_c09, mm_psf_PRS_c09, R2_pred_PRS_c09, PRESSRMS_PRS_c09, std_c09_PRS;act_psf_PRS_c10, mm_psf_PRS_c10, R2_pred_PRS_c10, PRESSRMS_PRS_c10, std_c10_PRS; act_psf_PRS_c11, mm_psf_PRS_c11, R2_pred_PRS_c11, PRESSRMS_PRS_c11, std_c11_PRS;act_psf_PRS_c12, mm_psf_PRS_c12, R2_pred_PRS_c10, PRESSRMS_PRS_c12, std_c12_PRS];
