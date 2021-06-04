clc;clear;close all;

%% Input variable bounds
n_design_var = 2;
n_random_var = 0;
n_var=n_design_var+n_random_var;

lbx = [0.1 0.1];
ubx = [10 10];
sd=[0.3 0.3];
Nmcs=1e6;
Pftarget = 0.00135;

%% Load Data
data = load('Initial DoE.mat');
scalenewDoe1_c01 = data.scalenewDoe1_c01;
scalenewDoe1_c02 = data.scalenewDoe1_c02;
scalenewDoe1_c03 = data.scalenewDoe1_c03;
scalenewDoe1_obj = data.scalenewDoe1_obj;

data_1 = load('AS DoE.mat');
t1_ls_01_new = data_1.t1_ls_01_new;
t1_ls_02_new = data_1.t1_ls_02_new;
t1_ls_03_new = data_1.t1_ls_03_new;
t1_obj_new = data_1.t1_obj_new;

%% Constraint 1

% Bounds
lbx_c01 = lbx(1:2);
ubx_c01 = ubx(1:2);

% Sample from the input space
n_var_c01 = 2;
npoints_c01=10*n_var_c01;
% scalenewDoe1_c01 =zeros(npoints_c01,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c01(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c01,1));
% end

% Sampling points from sample mean and standard deviation
doe1_ls_01 = scalenewDoe1_c01;
doe1_norm_ls_01 = normalizeuv(doe1_ls_01(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_01 = doe1_ls_01;
x_c01 = doe1_ls_01;
x_norm_c01 = doe1_norm_ls_01;

% Evaluating the Limit State 1
y_ls_01  = constraint1c(x_ls_01);

% Computing active subspace
k_ls_01 = 1;
m_c01 = 2;
npoints_active_susbpace_ls_01 = 10*m_c01;
M_ls_01 = npoints_active_susbpace_ls_01-1;
p_ls_01 = floor(10*m_c01)-1;
[b_ls_01_llm,C_ls_01_llm,W_ls_01_llm,Ev_ls_01_llm,avg_r2_pred_ls_01_llm]=llm(x_norm_c01,y_ls_01,M_ls_01,p_ls_01);

% Define active Subspace dimension
pv_ls_01 = cumsum(diag(Ev_ls_01_llm))/sum(diag(Ev_ls_01_llm));
n_ls_01  = active_subspace_dimension(diag(Ev_ls_01_llm));

%Compute active subspace and active variable
W1_ls_01_llm = W_ls_01_llm(:,1:n_ls_01);
t1_ls_01 = x_norm_c01*W1_ls_01_llm;

% Doe in the active subspace
t1_ls_01_min = -diag(sign(W1_ls_01_llm)'*W1_ls_01_llm);
t1_ls_01_max = diag(sign(W1_ls_01_llm)'*W1_ls_01_llm);
t1_ls_01_range = t1_ls_01_max'-t1_ls_01_min';
N_ls_01 = 10*n_ls_01+1;
% t1_ls_01_new = t1_ls_01_min'+((t1_ls_01_range).*lhsdesign(N_ls_01,n_ls_01));

% Mapping to original space
sp_c01 = [0 0];
lbx_norm_c01 = normalizeuv(lbx_c01,lbx_c01,ubx_c01);
ubx_norm_c01 = normalizeuv(ubx_c01,lbx_c01,ubx_c01);
xmin_c01 = zeros(N_ls_01,m_c01);
fval_c01 = zeros(N_ls_01,1);
exitflag_c01 = zeros(N_ls_01,1);

for i=1:N_ls_01
func_c01 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xmin_c01(i,:),fval_c01(i,1),exitflag_c01(i,1)] = fmincon(func_c01,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t1_ls_01_new(i,:)));
end

npoints_as=10*n_ls_01+1;
npoints_as_c01 = 10*n_ls_01+1;
x_norm_c01_org = xmin_c01(exitflag_c01==1,:);
x_norm_c01_org_chose = x_norm_c01_org(1:npoints_as_c01,:); 
x_c01_org_new = denormalizeuv(x_norm_c01_org_chose,lbx_c01,ubx_c01);

% Evaluate PSF
x_c01_new = x_c01_org_new;
[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(x_c01_new,lbx,ubx,sd,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
mean_psf_c01 = mean(psf_c01);
t1_c01 = x_norm_c01_org_chose*W1_ls_01_llm;
y_c01 = psf_c01;
[srgt_RBF_c01, PRESSRMS_RBF_c01, eXV_RBF_c01, srgtOPT_RBF_c01, Y_hat_RBF_c01, R2_pred_RBF_c01] = metamodel_RBF(t1_c01,y_c01);
[srgt_KRG_c01, PRESSRMS_KRG_c01, eXV_KRG_c01, srgtOPT_KRG_c01, Y_hat_KRG_c01, pred_var_KRG_c01, R2_pred_KRG_c01] = metamodel_KRG(t1_c01,y_c01);
[srgt_PRS_c01, PRESSRMS_PRS_c01, eXV_PRS_c01, srgtOPT_PRS_c01, Y_hat_PRS_c01, pred_var_PRS_c01, R2_pred_PRS_c01] = metamodel_PRS(t1_c01,y_c01);
[srgt_WAS_c01, PRESSRMS_WAS_c01, eXV_WAS_c01, srgtOPT_WAS_c01, Y_hat_WAS_c01, pred_var_WAS_c01, R2_pred_WAS_c01] = metamodel_WAS(t1_c01,y_c01);

%% Constraint 2

% Bounds
lbx_c02 = lbx(1:2);
ubx_c02 = ubx(1:2);

% Sample from the input space
n_var_c02 = 2;
npoints_c02=10*n_var_c02;
% scalenewDoe1_c02 = scalenewDoe1_c01;

% Sampling points from sample mean and standard deviation
doe1_ls_02 = scalenewDoe1_c02;
doe1_norm_ls_02 = normalizeuv(doe1_ls_02(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 2
x_ls_02 = doe1_ls_02;
x_c02 = doe1_ls_02;
x_norm_c02 = doe1_norm_ls_02;

% Evaluating the Limit State 2
y_ls_02  = constraint2c(x_ls_02);

% Computing active subspace
k_ls_02=1;
m_c02=2;
npoints_active_susbpace_ls_02 = 10*m_c02;
M_ls_02 = npoints_active_susbpace_ls_02-1; 
p_ls_02 = floor(10*m_c02)-1;
[b_ls_02_llm,C_ls_02_llm,W_ls_02_llm,Ev_ls_02_llm]=llm(x_norm_c02,y_ls_02,M_ls_02,p_ls_02);

% Define active Subspace dimension
pv_ls_02 = cumsum(diag(Ev_ls_02_llm))/sum(diag(Ev_ls_02_llm));
n_ls_02  = active_subspace_dimension(diag(Ev_ls_01_llm));

% Compute active subspace and active variable
W1_ls_02_llm = W_ls_02_llm(:,1:n_ls_02);
t1_ls_02 = x_norm_c02*W1_ls_02_llm;

% Doe in the active subspace
t1_ls_02_min = -diag(sign(W1_ls_02_llm)'*W1_ls_02_llm);
t1_ls_02_max = diag(sign(W1_ls_02_llm)'*W1_ls_02_llm);
t1_ls_02_range = t1_ls_02_max'-t1_ls_02_min';
N_ls_02 = 10*n_ls_02+1;
% t1_ls_02_new = t1_ls_02_min'+((t1_ls_02_range).*lhsdesign(N_ls_02,n_ls_02));

% Mapping to original space
sp_c02 = [0 0];
lbx_norm_c02 = normalizeuv(lbx_c02,lbx_c02,ubx_c02);
ubx_norm_c02 = normalizeuv(ubx_c02,lbx_c02,ubx_c02);

xmin_c02 = zeros(N_ls_02,m_c02);
fval_c02 = zeros(N_ls_02,1);
exitflag_c02 = zeros(N_ls_02,1);

for i=1:N_ls_02
func_c02 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xmin_c02(i,:),fval_c02(i,1),exitflag_c02(i,1)] = fmincon(func_c02,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t1_ls_02_new(i,:)));
end

npoints_as_c02=10*n_ls_02+1;
x_norm_c02_org = xmin_c02(exitflag_c02==1,:);
x_norm_c02_org_chose = x_norm_c02_org(1:npoints_as_c02,:);
x_c02_org_new = denormalizeuv(x_norm_c02_org_chose,lbx_c02,ubx_c02);

% Evaluate PSF
x_c02_new = x_c02_org_new;
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2c(x_c02_new,lbx,ubx,sd,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
mean_psf_c02 = mean(psf_c02);
t1_c02 = x_norm_c02_org_chose*W1_ls_02_llm;
y_c02 = psf_c02;
[srgt_RBF_c02, PRESSRMS_RBF_c02, eXV_RBF_c02, srgtOPT_RBF_c02, Y_hat_RBF_c02, R2_pred_RBF_c02] = metamodel_RBF(t1_c02,y_c02);
[srgt_KRG_c02, PRESSRMS_KRG_c02, eXV_KRG_c02, srgtOPT_KRG_c02, Y_hat_KRG_c02, pred_var_KRG_c02, R2_pred_KRG_c02] = metamodel_KRG(t1_c02,y_c02);
[srgt_PRS_c02, PRESSRMS_PRS_c02, eXV_PRS_c02, srgtOPT_PRS_c02, Y_hat_PRS_c02, pred_var_PRS_c02, R2_pred_PRS_c02] = metamodel_PRS(t1_c02,y_c02);
[srgt_WAS_c02, PRESSRMS_WAS_c02, eXV_WAS_c02, srgtOPT_WAS_c02, Y_hat_WAS_c02, pred_var_WAS_c02, R2_pred_WAS_c02] = metamodel_WAS(t1_c02,y_c02);

%% Constraint 3

% Bounds
lbx_c03 = lbx;
ubx_c03 = ubx;

% Sample from the input space
n_var_c03 = 2;
npoints_c03=10*n_var_c03;
% scalenewDoe1_c03 = scalenewDoe1_c01;

% Sampling points from sample mean and standard deviation
doe1_ls_03 = scalenewDoe1_c03;
doe1_norm_ls_03 = normalizeuv(doe1_ls_03(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 3
x_ls_03 = doe1_ls_03;
x_c03 = doe1_ls_03;
x_norm_c03 = doe1_norm_ls_03;

% Evaluating the Limit State 3
y_ls_03  = constraint3c(x_ls_03);

% Computing active subspace
k_ls_03=1;
m_c03=2;
npoints_active_susbpace_ls_03 = 10*m_c03;
M_ls_03 = npoints_active_susbpace_ls_03-1;
p_ls_03 = floor(10*m_c03)-1;
[b_ls_03_llm,C_ls_03_llm,W_ls_03_llm,Ev_ls_03_llm]=llm(x_norm_c03,y_ls_03,M_ls_03,p_ls_03);

% Define active Subspace dimension
pv_ls_03 = cumsum(diag(Ev_ls_03_llm))/sum(diag(Ev_ls_03_llm));
n_ls_03  = active_subspace_dimension(diag(Ev_ls_03_llm));

% Compute active subspace and active variable
W1_ls_03_llm = W_ls_03_llm(:,1:n_ls_03);
t1_ls_03 = x_norm_c03*W1_ls_03_llm;

% Doe in the active subspace
t1_ls_03_min = -diag(sign(W1_ls_03_llm)'*W1_ls_03_llm);
t1_ls_03_max = diag(sign(W1_ls_03_llm)'*W1_ls_03_llm);
t1_ls_03_range = t1_ls_03_max'-t1_ls_03_min';
N_ls_03 = 10*n_ls_03+1;
% t1_ls_03_new = t1_ls_03_min'+((t1_ls_03_range).*lhsdesign(N_ls_03,n_ls_03));

% Mapping to original space
sp_c03 = [0 0];
lbx_norm_c03 = normalizeuv(lbx_c03,lbx_c03,ubx_c03);
ubx_norm_c03 = normalizeuv(ubx_c03,lbx_c03,ubx_c03);
xmin_c03 = zeros(N_ls_03,m_c03);
fval_c03 = zeros(N_ls_03,1);
exitflag_c03 = zeros(N_ls_03,1);

for i=1:N_ls_03
func_c03 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xmin_c03(i,:),fval_c03(i,1),exitflag_c03(i,1)] = fmincon(func_c03,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t1_ls_03_new(i,:)));
end

npoints_as_c03=10*n_ls_03+1;
x_norm_c03_org = xmin_c03(exitflag_c03==1,:);
x_norm_c03_org_chose= x_norm_c03_org(1:npoints_as_c03,:);
x_c03_org_new = denormalizeuv(x_norm_c03_org_chose,lbx_c03,ubx_c03);

% Evaluate PSF
x_c03_new =  x_c03_org_new;
[psf_c03,pf_c03,re_c03,beta_c03] = mcspsfconstraint3c(x_c03_new,lbx,ubx,sd,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
mean_psf_c03 = mean(psf_c03);
t1_c03 = x_norm_c03_org_chose*W1_ls_03_llm;
y_c03 = psf_c03;
[srgt_RBF_c03, PRESSRMS_RBF_c03, eXV_RBF_c03, srgtOPT_RBF_c03, Y_hat_RBF_c03, R2_pred_RBF_c03] = metamodel_RBF(t1_c03,y_c03);
[srgt_KRG_c03, PRESSRMS_KRG_c03, eXV_KRG_c03, srgtOPT_KRG_c03, Y_hat_KRG_c03, pred_var_KRG_c03, R2_pred_KRG_c03] = metamodel_KRG(t1_c03,y_c03);
[srgt_PRS_c03, PRESSRMS_PRS_c03, eXV_PRS_c03, srgtOPT_PRS_c03, Y_hat_PRS_c03, pred_var_PRS_c03, R2_pred_PRS_c03] = metamodel_PRS(t1_c03,y_c03);
[srgt_WAS_c03, PRESSRMS_WAS_c03, eXV_WAS_c03, srgtOPT_WAS_c03, Y_hat_WAS_c03, pred_var_WAS_c03, R2_pred_WAS_c03] = metamodel_WAS(t1_c03,y_c03);

%% Objective Function

% Bounds
lbx_obj = lbx;
ubx_obj = ubx;

% Sample from the input space
n_var_obj = 2;
npoints_obj=10*n_var_obj;
% scalenewDoe1_obj = scalenewDoe1_c01;

% Sampling points from sample mean and standard deviation
doe1_obj = scalenewDoe1_obj;
doe1_norm_obj = normalizeuv(doe1_obj,lbx,ubx);

% Input Doe for Constraint 1
x_obj = doe1_obj;
x_norm_obj = doe1_norm_obj;

% Evaluating the Limit State 1
y_obj  = Weightc(x_obj);

% Computing active subspace
k_obj=1;
m_obj=2;
npoints_active_susbpace_obj = 10*m_obj;
M_obj = npoints_active_susbpace_obj-1;
p_obj = floor(10*m_obj)-1;
[b_obj_llm,C_obj_llm,W_obj_llm,Ev_obj_llm]=llm(x_norm_obj,y_obj,M_obj,p_obj);

% Define active Subspace dimension
pv_obj = cumsum(diag(Ev_obj_llm))/sum(diag(Ev_obj_llm));
n_obj  = active_subspace_dimension(diag(Ev_obj_llm));

% Compute active subspace and active variable
W1_obj_llm = W_obj_llm(:,1:n_obj);
t1_obj = x_norm_obj*W1_obj_llm;

% Doe in the active subspace
t1_obj_min = -diag(sign(W1_obj_llm)'*W1_obj_llm);
t1_obj_max = diag(sign(W1_obj_llm)'*W1_obj_llm);
t1_obj_range = t1_obj_max'-t1_obj_min';
N_obj = 10*n_obj+1;
% t1_obj_new = t1_obj_min'+((t1_obj_range).*lhsdesign(N_obj,n_obj));

% Mapping to original space
sp_obj = [0 0];
lbx_norm_obj = normalizeuv(lbx_obj,lbx_obj,ubx_obj);
ubx_norm_obj = normalizeuv(ubx_obj,lbx_obj,ubx_obj);
xmin_obj = zeros(N_obj,m_obj);
fval_obj = zeros(N_obj,1);
exitflag_obj = zeros(N_obj,1);

for i=1:N_obj
func_obj =  @(x_opt) (x_opt*zeros(m_obj,1));
[xmin_obj(i,:),fval_obj(i,1),exitflag_obj(i,1)] = fmincon(func_obj,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t1_obj_new(i,:)));
end

npoints_as_obj=10*n_obj+1;
x_norm_obj_org = xmin_obj(exitflag_obj==1,:);
x_norm_obj_org_chose = x_norm_obj_org(1:npoints_as_obj,:);
x_obj_org_new = denormalizeuv(x_norm_obj_org_chose,lbx_obj,ubx_obj);

% Evaluate PSF
x_obj_new =  x_obj_org_new;
y_obj_new = Weightc(x_obj_new);

%Build the metamodel in active subspace of Limit State
mean_obj_new = mean(y_obj_new);
obj_new_range = max(y_obj_new)-min(y_obj_new);
t1_obj_1 = x_norm_obj_org_chose*W1_obj_llm;
y_obj_1 = y_obj_new;
[srgt_RBF_obj_1, PRESSRMS_RBF_obj_1, eXV_RBF_obj_1, srgtOPT_RBF_obj_1, Y_hat_RBF_obj_1, R2_pred_RBF_obj_1] = metamodel_RBF(t1_obj_1,y_obj_1);
[srgt_KRG_obj_1, PRESSRMS_KRG_obj_1, eXV_KRG_obj_1, srgtOPT_KRG_obj_1, Y_hat_KRG_obj_1, pred_var_KRG_obj_1, R2_pred_KRG_obj_1] = metamodel_KRG(t1_obj_1,y_obj_1);
[srgt_PRS_obj_1, PRESSRMS_PRS_obj_1, eXV_PRS_obj_1, srgtOPT_PRS_obj_1, Y_hat_PRS_obj_1, pred_var_PRS_obj_1, R2_pred_PRS_obj_1] = metamodel_PRS(t1_obj_1,y_obj_1);
[srgt_WAS_obj_1, PRESSRMS_WAS_obj_1, eXV_WAS_obj_1, srgtOPT_WAS_obj_1, Y_hat_WAS_obj_1, pred_var_WAS_obj_1, R2_pred_WAS_obj_1] = metamodel_WAS(t1_obj_1,y_obj_1);

