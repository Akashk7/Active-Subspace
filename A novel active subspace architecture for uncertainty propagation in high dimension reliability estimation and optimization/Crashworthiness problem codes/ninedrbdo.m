clc;clear;close all;

%% Input variable bounds
n_design_var = 7;
n_random_var = 4;
n_var=n_design_var+n_random_var;

lbx = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
ubx = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
sd=[0.03 0.03 0.03 0.03 0.05 0.03 0.03];
mup1 = 0.345;
mup2 = 0.192;
mup3 = 0;
mup4 = 0;
sdp1 = 0.006;
sdp2 = 0.006;
sdp3=10;
sdp4=10;
mup = [mup1 mup2 mup3 mup4];
sdp = [sdp1 sdp2 sdp3 sdp4];

%% Sampling from the input variable space
npoints=3*n_var;
scalenewDoe1 =zeros(npoints,n_var);
for i=1:n_design_var
scalenewDoe1(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints,1));
end
% 
%% Test DoE
n_test=50;
x_test=zeros(n_test,n_design_var);

for i=1:n_design_var
x_test(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(n_test,1));
end

%% Constructing Metamodels:
Nmcs=1e6;
Pftarget=0.00135;

% Objective function:
x_obj = [scalenewDoe1(:,1:5) scalenewDoe1(:,7)];
y_obj=Weightc(x_obj);
[srgt_RBF_obj, PRESSRMS_RBF_obj, eXV_RBF_obj, srgtOPT_RBF_obj, Y_hat_RBF_obj, R2_pred_RBF_obj] = metamodel_RBF(x_obj,y_obj);
[srgt_KRG_obj, PRESSRMS_KRG_obj, eXV_KRG_obj, srgtOPT_KRG_obj, Y_hat_KRG_obj, pred_var_KRG_obj, R2_pred_KRG_obj] = metamodel_KRG(x_obj,y_obj);
[srgt_PRS_obj, PRESSRMS_PRS_obj, eXV_PRS_obj, srgtOPT_PRS_obj, Y_hat_PRS_obj, pred_var_PRS_obj, R2_pred_PRS_obj] = metamodel_PRS(x_obj,y_obj);
[srgt_WAS_obj, PRESSRMS_WAS_obj, eXV_WAS_obj, srgtOPT_WAS_obj, Y_hat_WAS_obj, pred_var_WAS_obj, R2_pred_WAS_obj] = metamodel_WAS(x_obj,y_obj);

x_test_1 = [x_test(:,1:5) x_test(:,7)];
y_test_obj= Weightc(x_test_1);

y_predicted_RBF_obj = srgtsRBFEvaluate(x_test_1, srgt_RBF_obj);
y_predicted_KRG_obj = srgtsKRGEvaluate(x_test_1, srgt_KRG_obj);
y_predicted_PRS_obj = srgtsPRSEvaluate(x_test_1, srgt_PRS_obj);
y_predicted_WAS_obj = srgtsWASEvaluate(x_test_1, srgt_WAS_obj);

[R2_RBF_obj,RMSE_RBF_obj]=metrics(y_test_obj,y_predicted_RBF_obj);
[R2_KRG_obj,RMSE_KRG_obj]=metrics(y_test_obj,y_predicted_KRG_obj);
[R2_PRS_obj,RMSE_PRS_obj]=metrics(y_test_obj,y_predicted_PRS_obj);
[R2_WAS_obj,RMSE_WAS_obj]=metrics(y_test_obj,y_predicted_WAS_obj);

%Constraint1
[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(scalenewDoe1,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
x_c01 = [scalenewDoe1(:,2) scalenewDoe1(:,3) scalenewDoe1(:,4) scalenewDoe1(:,6)];
y_c01 = psf_c01;
[srgt_RBF_c01, PRESSRMS_RBF_c01, eXV_RBF_c01, srgtOPT_RBF_c01, Y_hat_RBF_c01, R2_pred_RBF_c01] = metamodel_RBF(x_c01,y_c01);
[srgt_KRG_c01, PRESSRMS_KRG_c01, eXV_KRG_c01, srgtOPT_KRG_c01, Y_hat_KRG_c01, pred_var_KRG_c01, R2_pred_KRG_c01] = metamodel_KRG(x_c01,y_c01);
[srgt_PRS_c01, PRESSRMS_PRS_c01, eXV_PRS_c01, srgtOPT_PRS_c01, Y_hat_PRS_c01, pred_var_PRS_c01, R2_pred_PRS_c01] = metamodel_PRS(x_c01,y_c01);
[srgt_WAS_c01, PRESSRMS_WAS_c01, eXV_WAS_c01, srgtOPT_WAS_c01, Y_hat_WAS_c01, pred_var_WAS_c01, R2_pred_WAS_c01] = metamodel_WAS(x_c01,y_c01);

%Compute the active subspace
k_c01=1;
alpha=10;
m_c01=4;
npoints_active_susbpace_c01 = 10*m_c01;
M_c01 = min(ceil(alpha*k_c01*log(m_c01)),npoints_active_susbpace_c01-1);
p_c01 = 39;
[b_c01_llm,C_c01_llm,W_c01_llm,Ev_c01_llm]=llm(x_c01,y_c01,M_c01,p_c01);

%Define the dimension of active subspace
pv_c01_llm=cumsum(diag(Ev_c01_llm))/sum(diag(Ev_c01_llm));
n_c01_llm = 1;

%Building the metamodel in the approximated active space:
lbx_c01 = [lbx(2) lbx(3) lbx(4) lbx(6)];
ubx_c01 = [ubx(2) ubx(3) ubx(4) ubx(6)];
x_norm_c01 = ((x_c01 - lbx_c01)./((ubx_c01-lbx_c01)./2))-1;
W1_c01_llm = W_c01_llm(:,1:n_c01_llm);
t_c01_llm = x_norm_c01*W1_c01_llm;
[srgt_RBF_c01_as_llm, PRESSRMS_RBF_c01_as_llm, eXV_RBF_c01_as_llm, srgtOPT_RBF_c01_as_llm,Y_hat_RBF_c01_as_llm, R2_pred_RBF_c01_as_llm] = metamodel_RBF(t_c01_llm,y_c01);
[srgt_KRG_c01_as_llm, PRESSRMS_KRG_c01_as_llm, eXV_KRG_c01_as_llm, srgtOPT_KRG_c01_as_llm,Y_hat_KRG_c01_as_llm, pred_Var_KRG_c01_as_llm, R2_pred_KRG_c01_as_llm] = metamodel_KRG(t_c01_llm,y_c01);
[srgt_PRS_c01_as_llm, PRESSRMS_PRS_c01_as_llm, eXV_PRS_c01_as_llm, srgtOPT_PRS_c01_as_llm,Y_hat_PRS_c01_as_llm, pred_Var_PRS_c01_as_llm, R2_pred_PRS_c01_as_llm] = metamodel_PRS(t_c01_llm,y_c01);
[srgt_WAS_c01_as_llm, PRESSRMS_WAS_c01_as_llm, eXV_WAS_c01_as_llm, srgtOPT_WAS_c01_as_llm,Y_hat_WAS_c01_as_llm, pred_Var_WAS_c01_as_llm, R2_pred_WAS_c01_as_llm] = metamodel_WAS(t_c01_llm,y_c01);

%Compute the active subspace
k_ls_01=1;
alpha=10;
m_ls_01=4;
npoints_active_susbpace_ls_01 = 10*m_ls_01;
M_ls_01 = min(ceil(alpha*k_ls_01*log(m_ls_01)),npoints_active_susbpace_ls_01-1);
p_ls_01 = 39;
y_ls_01 = constraint1c([x_c01 repmat([mup2, mup3, mup4],npoints,1)]); 
[b_ls_01_llm,C_ls_01_llm,W_ls_01_llm,Ev_ls_01_llm]=llm(x_c01,y_ls_01,M_ls_01,p_ls_01);

%Define the dimension of active subspace
pv_ls_01_llm=cumsum(diag(Ev_ls_01_llm))/sum(diag(Ev_ls_01_llm));
n_ls_01_llm = 1;

%Mapping to the original space
N_ls_01 = 400;
n_design_var_ls_01 = 4;
[x_new_ls_01_llm, t_new_ls_01_llm]= maptooriginalspace(n_ls_01_llm,W_ls_01_llm,N_ls_01,n_design_var_ls_01,lbx_c01,ubx_c01);
x_new_ls_01_llm_1 = [zeros(size(x_new_ls_01_llm,1),1) x_new_ls_01_llm(:,1:3) zeros(size(x_new_ls_01_llm,1),1) x_new_ls_01_llm(:,4) zeros(size(x_new_ls_01_llm,1),1)]; 
[psf_new_ls_01_llm,pf_new_ls_01_llm,re_new_ls_01_llm,beta_new_ls_01_llm] = mcspsfconstraint1c(x_new_ls_01_llm_1,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

y_ls_01_new_llm = psf_new_ls_01_llm;
t_ls_01_llm_new = t_new_ls_01_llm;

[srgt_RBF_ls_01_as_llm_new, PRESSRMS_RBF_ls_01_as_llm_new, eXV_RBF_ls_01_as_llm_new, srgtOPT_RBF_ls_01_as_llm_new,Y_hat_RBF_ls_01_as_llm_new,R2_pred_RBF_ls_01_as_llm_new] = metamodel_RBF(t_ls_01_llm_new,y_ls_01_new_llm);
[srgt_KRG_ls_01_as_llm_new, PRESSRMS_KRG_ls_01_as_llm_new, eXV_KRG_ls_01_as_llm_new, srgtOPT_KRG_ls_01_as_llm_new,Y_hat_KRG_ls_01_as_llm_new, pred_Var_KRG_ls_01_as_llm_new, R2_pred_KRG_ls_01_as_llm_new] = metamodel_KRG(t_ls_01_llm_new,y_ls_01_new_llm);
[srgt_PRS_ls_01_as_llm_new, PRESSRMS_PRS_ls_01_as_llm_new, eXV_PRS_ls_01_as_llm_new, srgtOPT_PRS_ls_01_as_llm_new,Y_hat_PRS_ls_01_as_llm_new, pred_Var_PRS_ls_01_as_llm_new, R2_pred_PRS_ls_01_as_llm_new] = metamodel_PRS(t_ls_01_llm_new,y_ls_01_new_llm);
[srgt_WAS_ls_01_as_llm_new, PRESSRMS_WAS_ls_01_as_llm_new, eXV_WAS_ls_01_as_llm_new, srgtOPT_WAS_ls_01_as_llm_new,Y_hat_WAS_ls_01_as_llm_new, pred_Var_WAS_ls_01_as_llm_new, R2_pred_WAS_ls_01_as_llm_new] = metamodel_WAS(t_ls_01_llm_new,y_ls_01_new_llm);

%Mapping to the original space
N_c01_llm = 400;
n_design_var_c01 = 4;
[x_new_c01_llm, t_new_c01_llm]= maptooriginalspace(n_c01_llm,W_c01_llm,N_c01_llm,n_design_var_c01,lbx_c01,ubx_c01);
x_new_c01_llm_1 = [zeros(size(x_new_c01_llm,1),1) x_new_c01_llm(:,1:3) zeros(size(x_new_c01_llm,1),1) x_new_c01_llm(:,4) zeros(size(x_new_c01_llm,1),1)]; 
[psf_new_c01_llm,pf_new_c01_llm,re_new_c01_llm,beta_new_c01_llm] = mcspsfconstraint1c(x_new_c01_llm_1,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

y_c01_new_llm = psf_new_c01_llm;
t_c01_llm_new = t_new_c01_llm;

[srgt_RBF_c01_as_llm_new, PRESSRMS_RBF_c01_as_llm_new, eXV_RBF_c01_as_llm_new, srgtOPT_RBF_c01_as_llm_new,Y_hat_RBF_c01_as_llm_new,R2_pred_RBF_c01_as_llm_new] = metamodel_RBF(t_c01_llm_new,y_c01_new_llm);
[srgt_KRG_c01_as_llm_new, PRESSRMS_KRG_c01_as_llm_new, eXV_KRG_c01_as_llm_new, srgtOPT_KRG_c01_as_llm_new,Y_hat_KRG_c01_as_llm_new, pred_Var_KRG_c01_as_llm_new, R2_pred_KRG_c01_as_llm_new] = metamodel_KRG(t_c01_llm_new,y_c01_new_llm);
[srgt_PRS_c01_as_llm_new, PRESSRMS_PRS_c01_as_llm_new, eXV_PRS_c01_as_llm_new, srgtOPT_PRS_c01_as_llm_new,Y_hat_PRS_c01_as_llm_new, pred_Var_PRS_c01_as_llm_new, R2_pred_PRS_c01_as_llm_new] = metamodel_PRS(t_c01_llm_new,y_c01_new_llm);
[srgt_WAS_c01_as_llm_new, PRESSRMS_WAS_c01_as_llm_new, eXV_WAS_c01_as_llm_new, srgtOPT_WAS_c01_as_llm_new,Y_hat_WAS_c01_as_llm_new, pred_Var_WAS_c01_as_llm_new, R2_pred_WAS_c01_as_llm_new] = metamodel_WAS(t_c01_llm_new,y_c01_new_llm);

%Building the metamodel in the actual active space:
[b_c01_actual,C_c01_actual,W_c01_actual,Ev_c01_actual]= asmc1(x_c01,mup);

pv_c01_actual=cumsum(diag(Ev_c01_actual))/sum(diag(Ev_c01_actual));
n_c01_actual = 1;

W1_c01_actual = W_c01_actual(:,1:n_c01_actual);
t_c01_actual = x_norm_c01*W1_c01_actual;

[srgt_RBF_c01_as_actual, PRESSRMS_RBF_c01_as_actual, eXV_RBF_c01_as_actual, srgtOPT_RBF_c01_as_actual,Y_hat_RBF_c01_as_actual,R2_pred_RBF_c01_as_actual] = metamodel_RBF(t_c01_actual,y_c01);
[srgt_KRG_c01_as_actual, PRESSRMS_KRG_c01_as_actual, eXV_KRG_c01_as_actual, srgtOPT_KRG_c01_as_actual,Y_hat_KRG_c01_as_actual, pred_Var_KRG_c01_as_actual, R2_pred_KRG_c01_as_actual] = metamodel_KRG(t_c01_actual,y_c01);
[srgt_PRS_c01_as_actual, PRESSRMS_PRS_c01_as_actual, eXV_PRS_c01_as_actual, srgtOPT_PRS_c01_as_actual,Y_hat_PRS_c01_as_actual, pred_Var_PRS_c01_as_actual, R2_pred_PRS_c01_as_actual] = metamodel_PRS(t_c01_actual,y_c01);
[srgt_WAS_c01_as_actual, PRESSRMS_WAS_c01_as_actual, eXV_WAS_c01_as_actual, srgtOPT_WAS_c01_as_actual,Y_hat_WAS_c01_as_actual, pred_Var_WAS_c01_as_actual, R2_pred_WAS_c01_as_actual] = metamodel_WAS(t_c01_actual,y_c01);

%Mapping to original space:
N_c01 = 150;
n_design_var_c01 = 4;
[x_new_c01_actual, t_new_c01_actual]= maptooriginalspace(n_c01_actual,W_c01_actual,N_c01,n_design_var_c01,lbx_c01,ubx_c01);
x_new_c01_actual_1 = [zeros(size(x_new_c01_actual,1),1) x_new_c01_actual(:,1:3) zeros(size(x_new_c01_actual,1),1) x_new_c01_actual(:,4) zeros(size(x_new_c01_actual,1),1)]; 
[psf_new_c01_actual,pf_new_c01_actual,re_new_c01_actual,beta_new_c01_actual] = mcspsfconstraint1c(x_new_c01_actual_1,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

y_c01_new_actual = psf_new_c01_actual;
t_c01_actual_new = t_new_c01_actual;

[srgt_RBF_c01_as_actual_new, PRESSRMS_RBF_c01_as_actual_new, eXV_RBF_c01_as_actual_new, srgtOPT_RBF_c01_as_actual_new,Y_hat_RBF_c01_as_actual_new,R2_pred_RBF_c01_as_actual_new] = metamodel_RBF(t_c01_actual_new,y_c01_new_actual);
[srgt_KRG_c01_as_actual_new, PRESSRMS_KRG_c01_as_actual_new, eXV_KRG_c01_as_actual_new, srgtOPT_KRG_c01_as_actual_new,Y_hat_KRG_c01_as_actual_new, pred_Var_KRG_c01_as_actual_new, R2_pred_KRG_c01_as_actual_new] = metamodel_KRG(t_c01_actual_new,y_c01_new_actual);
[srgt_PRS_c01_as_actual_new, PRESSRMS_PRS_c01_as_actual_new, eXV_PRS_c01_as_actual_new, srgtOPT_PRS_c01_as_actual_new,Y_hat_PRS_c01_as_actual_new, pred_Var_PRS_c01_as_actual_new, R2_pred_PRS_c01_as_actual_new] = metamodel_PRS(t_c01_actual_new,y_c01_new_actual);
[srgt_WAS_c01_as_actual_new, PRESSRMS_WAS_c01_as_actual_new, eXV_WAS_c01_as_actual_new, srgtOPT_WAS_c01_as_actual_new,Y_hat_WAS_c01_as_actual_new, pred_Var_WAS_c01_as_actual_new, R2_pred_WAS_c01_as_actual_new] = metamodel_WAS(t_c01_actual_new,y_c01_new_actual);

r2_predicted = [R2_pred_KRG_c01, R2_pred_RBF_c01, R2_pred_PRS_c01, R2_pred_WAS_c01; R2_pred_KRG_c01_as_llm, R2_pred_RBF_c01_as_llm, R2_pred_PRS_c01_as_llm, R2_pred_WAS_c01_as_llm; R2_pred_KRG_c01_as_actual, R2_pred_RBF_c01_as_actual, R2_pred_PRS_c01_as_actual, R2_pred_WAS_c01_as_actual];
r2_predicted_new = [R2_pred_KRG_c01, R2_pred_RBF_c01, R2_pred_PRS_c01, R2_pred_WAS_c01; R2_pred_KRG_c01_as_actual_new, R2_pred_RBF_c01_as_actual_new, R2_pred_PRS_c01_as_actual_new, R2_pred_WAS_c01_as_actual_new; R2_pred_KRG_c01_as_llm_new, R2_pred_RBF_c01_as_llm_new, R2_pred_PRS_c01_as_llm_new, R2_pred_WAS_c01_as_llm_new; R2_pred_KRG_ls_01_as_llm_new, R2_pred_RBF_ls_01_as_llm_new, R2_pred_PRS_ls_01_as_llm_new, R2_pred_WAS_ls_01_as_llm_new];

%Constraint2
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2(scalenewDoe1,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c02 = psf_c02;
[srgt_RBF_c02, PRESSRMS_RBF_c02, eXV_RBF_c02, srgtOPT_RBF_c02,Y_hat_RBF_c02,R2_pred_RBF_c02] = metamodel_RBF(x,y_c2);

%Constraint3
psf_c3 = mcspsfconstraint3(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c3 = psf_c3;
[srgt_RBF_c3, PRESSRMS_RBF_c3, eXV_RBF_c3, srgtOPT_RBF_c3,Y_hat_RBF_c3,R2_pred_c3] = metamodel_RBF(x,y_c3);

%Constraint4
psf_c4 = mcspsfconstraint4(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c4 = psf_c4;
[srgt_RBF_c4, PRESSRMS_RBF_c4, eXV_RBF_c4, srgtOPT_RBF_c4,Y_hat_RBF_c4,R2_pred_c4] = metamodel_RBF(x,y_c4);

%Constraint5
psf_c5 = mcspsfconstraint5(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c5 = psf_c5;
[srgt_RBF_c5, PRESSRMS_RBF_c5, eXV_RBF_c5, srgtOPT_RBF_c5,Y_hat_RBF_c5,R2_pred_c5] = metamodel_RBF(x,y_c5);

%Constraint6
psf_c6 = mcspsfconstraint6(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c6 = psf_c6;
[srgt_RBF_c6, PRESSRMS_RBF_c6, eXV_RBF_c6, srgtOPT_RBF_c6,Y_hat_RBF_c6,R2_pred_c6] = metamodel_RBF(x,y_c6);

%Constraint7
psf_c7 = mcspsfconstraint7(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c7 = psf_c7;
[srgt_RBF_c7, PRESSRMS_RBF_c7, eXV_RBF_c7, srgtOPT_RBF_c7,Y_hat_RBF_c7,R2_pred_c7] = metamodel_RBF(x,y_c7);

%Constraint8
psf_c8 = mcspsfconstraint8(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c8 = psf_c8;
[srgt_RBF_c8, PRESSRMS_RBF_c8, eXV_RBF_c8, srgtOPT_RBF_c8,Y_hat_RBF_c8,R2_pred_c8] = metamodel_RBF(x,y_c8);

%Constraint9
psf_c9 = mcspsfconstraint9(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c9 = psf_c9;
[srgt_RBF_c9, PRESSRMS_RBF_c9, eXV_RBF_c9, srgtOPT_RBF_c9,Y_hat_RBF_c9,R2_pred_c9] = metamodel_RBF(x,y_c9);

%Constraint10
psf_c10 = mcspsfconstraint10(x,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
y_c10 = psf_c10;
[srgt_RBF_c10, PRESSRMS_RBF_c10, eXV_RBF_c10, srgtOPT_RBF_c10,Y_hat_RBF_c10,R2_pred_c10] = metamodel_RBF(x,y_c10);

%% Optimization:
func = @(t)(srgtsRBFEvaluate(t, srgt_RBF_obj));
t0=[0.511 1.417 0.500 1.352 0.658 1.473 0.500 0.345 0.192];
%t0=lbx;
[xmin,fval,exitflag,output]=fmincon(func,t0,[],[],[],[],lbx,ubx,@(t)nonlcon(t,srgt_RBF_c1,srgt_RBF_c2,srgt_RBF_c3,...
                                    srgt_RBF_c4,srgt_RBF_c5,srgt_RBF_c6,srgt_RBF_c7,srgt_RBF_c8,srgt_RBF_c9,srgt_RBF_c10,@srgtsRBFEvaluate));
%% Constraint Surrogate Values at the optima
surr_obj = srgtsRBFEvaluate(xmin, srgt_RBF_obj);
surr_psf_c1= srgtsRBFEvaluate(xmin, srgt_RBF_c1);
surr_psf_c2= srgtsRBFEvaluate(xmin, srgt_RBF_c2);
surr_psf_c3= srgtsRBFEvaluate(xmin, srgt_RBF_c3);
surr_psf_c4= srgtsRBFEvaluate(xmin, srgt_RBF_c4);
surr_psf_c5= srgtsRBFEvaluate(xmin, srgt_RBF_c5);
surr_psf_c6= srgtsRBFEvaluate(xmin, srgt_RBF_c6);
surr_psf_c7= srgtsRBFEvaluate(xmin, srgt_RBF_c7);
surr_psf_c8= srgtsRBFEvaluate(xmin, srgt_RBF_c8);
surr_psf_c9= srgtsRBFEvaluate(xmin, srgt_RBF_c9);
surr_psf_c10= srgtsRBFEvaluate(xmin, srgt_RBF_c10);

%% Validation
act_obj = Weight(xmin);
act_psf_c1 = mcspsfconstraint1(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c2 = mcspsfconstraint2(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);                                
act_psf_c3 = mcspsfconstraint3(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);                                
act_psf_c4 = mcspsfconstraint4(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c5 = mcspsfconstraint5(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c6 = mcspsfconstraint6(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c7 = mcspsfconstraint7(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c8 = mcspsfconstraint8(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget8);
act_psf_c9 = mcspsfconstraint9(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
act_psf_c10 = mcspsfconstraint10(xmin,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);


%% Standard Deviation:

error=abs([act_obj-surr_obj,act_psf_c1-surr_psf_c1,act_psf_c2-surr_psf_c2,act_psf_c3-surr_psf_c3,act_psf_c4-surr_psf_c4,act_psf_c5-surr_psf_c5,act_psf_c6-surr_psf_c6,act_psf_c7-surr_psf_c7,act_psf_c8-surr_psf_c8,act_psf_c9-surr_psf_c9,act_psf_c10-surr_psf_c10]);
press_rmse=[PRESSRMS_RBF_obj,PRESSRMS_RBF_c1,PRESSRMS_RBF_c2,PRESSRMS_RBF_c3,PRESSRMS_RBF_c4,PRESSRMS_RBF_c5,PRESSRMS_RBF_c6,PRESSRMS_RBF_c7,PRESSRMS_RBF_c8,PRESSRMS_RBF_c9,PRESSRMS_RBF_c10];
std_dev = error./press_rmse;
%% Normalized PRESSRMSE
y_obj_range=max(y_obj)-min(y_obj);
psf_c1_range=max(y_c1)-min(y_c1);
psf_c2_range=max(y_c2)-min(y_c2);
psf_c3_range=max(y_c3)-min(y_c3);
psf_c4_range=max(y_c4)-min(y_c4);
psf_c5_range=max(y_c5)-min(y_c5);
psf_c6_range=max(y_c6)-min(y_c6);
psf_c7_range=max(y_c7)-min(y_c7);
psf_c8_range=max(y_c8)-min(y_c8);
psf_c9_range=max(y_c9)-min(y_c9);
psf_c10_range=max(y_c10)-min(y_c10);
psf_range=[y_obj_range,psf_c1_range,psf_c2_range,psf_c3_range,psf_c4_range,psf_c5_range,psf_c6_range,psf_c7_range,psf_c8_range,psf_c9_range,psf_c10_range];
normalized_PRESS_RMSE=press_rmse./psf_range;

%% Lower Bound and Upper Bound Calculation
surr_values=[surr_obj,surr_psf_c1,surr_psf_c2,surr_psf_c3,surr_psf_c4,surr_psf_c5,surr_psf_c6,surr_psf_c7,surr_psf_c8,surr_psf_c9,surr_psf_c10];
upper_bound = surr_values+3*normalized_PRESS_RMSE;
lower_bound = surr_values-3*normalized_PRESS_RMSE;
std_dev_1=error./normalized_PRESS_RMSE;