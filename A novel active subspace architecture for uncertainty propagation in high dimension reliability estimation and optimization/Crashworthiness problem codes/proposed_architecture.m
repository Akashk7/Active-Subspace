clc;clear;close all;

%% Input variable bounds
n_design_var = 7;
n_random_var = 4;
n_var=n_design_var+n_random_var;

lbx = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
ubx = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
sd = [0.03 0.03 0.03 0.03 0.05 0.03 0.03];
mup1 = 0.345;
mup2 = 0.192;
mup3 = 0;
mup4 = 0;
sdp1 = 0.006;
sdp2 = 0.006;
sdp3 = 7;
sdp4 = 7;
mup = [mup1 mup2 mup3 mup4];
sdp = [sdp1 sdp2 sdp3 sdp4];

lbp = mup-3*sdp;
ubp = mup+3*sdp;

Nmcs=1e6;
Pftarget=0.0013;

load('Experiment Journal.mat');

%% Constraint 1

% Bounds
lbx_c01 = [lbx(2:4) lbx(6)];
ubx_c01 = [ubx(2:4) ubx(6)];

% Sample from the input space
n_var_c01 = 4;
npoints_c01=10*n_var_c01;
% scalenewDoe1_c01 =zeros(npoints_c01,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c01(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c01,1));
% end
% 
% for i=1:n_random_var
%     scalenewDoe1_c01(:,i+n_design_var) = repmat(mup(i),npoints_c01,1);
% end

% Sampling points from sample mean and standard deviation
doe1_ls_01 = scalenewDoe1_c01;
doe1_norm_ls_01 = normalizeuv(doe1_ls_01(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_01 = [doe1_ls_01(:,2:4) doe1_ls_01(:,6) doe1_ls_01(:,9:11)];
x_c01 = [doe1_ls_01(:,2:4) doe1_ls_01(:,6)];
x_norm_c01 = [doe1_norm_ls_01(:,2:4) doe1_norm_ls_01(:,6)];

% Evaluating the Limit State 1
y_ls_01  = constraint1c(x_ls_01);

% Computing active subspace
k_ls_01=2;
m_c01=4;
npoints_active_susbpace_ls_01 = 10*m_c01;
M_ls_01 = min(ceil(10*m_c01*log(m_c01)),npoints_active_susbpace_ls_01-1);
p_ls_01 =floor(10*m_c01)-1;
[b_ls_01_llm,C_ls_01_llm,W_ls_01_llm,Ev_ls_01_llm]=llm(x_norm_c01,y_ls_01,M_ls_01,p_ls_01);

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
sp_c01 = [0 0 0 0];
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
x_c01_new = [zeros(npoints_as,1) x_c01_org_new(:,1:3) zeros(npoints_as,1) x_c01_org_new(:,4) zeros(npoints_as,1)];
[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(x_c01_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c01 = x_norm_c01_org_chose*W1_ls_01_llm;
y_c01 = psf_c01;
[srgt_RBF_c01, PRESSRMS_RBF_c01, eXV_RBF_c01, srgtOPT_RBF_c01, Y_hat_RBF_c01, R2_pred_RBF_c01] = metamodel_RBF(t1_c01,y_c01);
[srgt_KRG_c01, PRESSRMS_KRG_c01, eXV_KRG_c01, srgtOPT_KRG_c01, Y_hat_KRG_c01, pred_var_KRG_c01, R2_pred_KRG_c01] = metamodel_KRG(t1_c01,y_c01);
[srgt_PRS_c01, PRESSRMS_PRS_c01, eXV_PRS_c01, srgtOPT_PRS_c01, Y_hat_PRS_c01, pred_var_PRS_c01, R2_pred_PRS_c01] = metamodel_PRS(t1_c01,y_c01);
[srgt_WAS_c01, PRESSRMS_WAS_c01, eXV_WAS_c01, srgtOPT_WAS_c01, Y_hat_WAS_c01, pred_var_WAS_c01, R2_pred_WAS_c01] = metamodel_WAS(t1_c01,y_c01);

%% Constraint 2

% Bounds
lbx_c02 = lbx(1:3);
ubx_c02 = ubx(1:3);

% Sample from the input space
n_var_c02 = 3;
npoints_c02=10*n_var_c02;
% scalenewDoe1_c02 =zeros(npoints_c02,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c02(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c02,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c02(:,i+n_design_var) = repmat(mup(i),npoints_c02,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_02 = scalenewDoe1_c02;
doe1_norm_ls_02 = normalizeuv(doe1_ls_02(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_02 = [doe1_ls_02(:,1:3) doe1_ls_02(:,8) doe1_ls_02(:,10)];
x_c02 = doe1_ls_02(:,1:3);
x_norm_c02 = doe1_norm_ls_02(:,1:3);

% Evaluating the Limit State 1
y_ls_02  = constraint2c(x_ls_02);

% Computing active subspace
k_ls_02=2;
m_c02=3;
npoints_active_susbpace_ls_02 = 10*m_c02;
M_ls_02 = min(ceil(10*m_c02*log(m_c02)),npoints_active_susbpace_ls_02-1);
p_ls_02 =floor(10*m_c02)-1;
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
sp_c02 = [0 0 0];
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
x_c02_new =  [x_c02_org_new(:,1:3) zeros(npoints_as_c02,4)];
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2c(x_c02_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c02 = x_norm_c02_org_chose*W1_ls_02_llm;
y_c02 = psf_c02;
[srgt_RBF_c02, PRESSRMS_RBF_c02, eXV_RBF_c02, srgtOPT_RBF_c02, Y_hat_RBF_c02, R2_pred_RBF_c02] = metamodel_RBF(t1_c02,y_c02);
[srgt_KRG_c02, PRESSRMS_KRG_c02, eXV_KRG_c02, srgtOPT_KRG_c02, Y_hat_KRG_c02, pred_var_KRG_c02, R2_pred_KRG_c02] = metamodel_KRG(t1_c02,y_c02);
[srgt_PRS_c02, PRESSRMS_PRS_c02, eXV_PRS_c02, srgtOPT_PRS_c02, Y_hat_PRS_c02, pred_var_PRS_c02, R2_pred_PRS_c02] = metamodel_PRS(t1_c02,y_c02);
[srgt_WAS_c02, PRESSRMS_WAS_c02, eXV_WAS_c02, srgtOPT_WAS_c02, Y_hat_WAS_c02, pred_var_WAS_c02, R2_pred_WAS_c02] = metamodel_WAS(t1_c02,y_c02);

%% Constraint 3

% Bounds
lbx_c03 = [lbx(1:3) lbx(5) lbx(7)];
ubx_c03 = [ubx(1:3) ubx(5) ubx(7)];

% Sample from the input space
n_var_c03 = 5;
npoints_c03=10*n_var_c03;
% scalenewDoe1_c03 =zeros(npoints_c03,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c03(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c03,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c03(:,i+n_design_var) = repmat(mup(i),npoints_c03,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_03 = scalenewDoe1_c03;
doe1_norm_ls_03 = normalizeuv(doe1_ls_03(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_03 = [doe1_ls_03(:,1:3) doe1_ls_03(:,5) doe1_ls_03(:,7) doe1_ls_03(:,8:10)];
x_c03 = [doe1_ls_03(:,1:3) doe1_ls_03(:,5) doe1_ls_03(:,7)];
x_norm_c03 = [doe1_norm_ls_03(:,1:3) doe1_norm_ls_03(:,5) doe1_norm_ls_03(:,7)];

% Evaluating the Limit State 1
y_ls_03  = constraint3c(x_ls_03);

% Computing active subspace
k_ls_03=2;
m_c03=5;
npoints_active_susbpace_ls_03 = 10*m_c03;
M_ls_03 = min(ceil(10*m_c03*log(m_c03)),npoints_active_susbpace_ls_03-1);
p_ls_03 =floor(10*m_c03)-1;
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
sp_c03 = [0 0 0 0 0];
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
x_c03_new =  [x_c03_org_new(:,1:3) zeros(npoints_as,1) x_c03_org_new(:,4) zeros(npoints_as,1) x_c03_org_new(:,5)];
[psf_c03,pf_c03,re_c03,beta_c03] = mcspsfconstraint3c(x_c03_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c03 = x_norm_c03_org_chose*W1_ls_03_llm;
y_c03 = psf_c03;
[srgt_RBF_c03, PRESSRMS_RBF_c03, eXV_RBF_c03, srgtOPT_RBF_c03, Y_hat_RBF_c03, R2_pred_RBF_c03] = metamodel_RBF(t1_c03,y_c03);
[srgt_KRG_c03, PRESSRMS_KRG_c03, eXV_KRG_c03, srgtOPT_KRG_c03, Y_hat_KRG_c03, pred_var_KRG_c03, R2_pred_KRG_c03] = metamodel_KRG(t1_c03,y_c03);
[srgt_PRS_c03, PRESSRMS_PRS_c03, eXV_PRS_c03, srgtOPT_PRS_c03, Y_hat_PRS_c03, pred_var_PRS_c03, R2_pred_PRS_c03] = metamodel_PRS(t1_c03,y_c03);
[srgt_WAS_c03, PRESSRMS_WAS_c03, eXV_WAS_c03, srgtOPT_WAS_c03, Y_hat_WAS_c03, pred_var_WAS_c03, R2_pred_WAS_c03] = metamodel_WAS(t1_c03,y_c03);

%% Constraint 4

% Bounds
lbx_c04 = [lbx(1:3) lbx(5:7)];
ubx_c04 = [ubx(1:3) ubx(5:7)];

% Sample from the input space
n_var_c04 = 6;
npoints_c04=10*n_var_c04;
% scalenewDoe1_c04 =zeros(npoints_c04,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c04(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c04,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c04(:,i+n_design_var) = repmat(mup(i),npoints_c04,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_04 = scalenewDoe1_c04;
doe1_norm_ls_04 = normalizeuv(doe1_ls_04(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_04 = [doe1_ls_04(:,1:3) doe1_ls_04(:,5:7) doe1_ls_04(:,8:10)];
x_c04 = [doe1_ls_04(:,1:3) doe1_ls_04(:,5:7)];
x_norm_c04 = [doe1_norm_ls_04(:,1:3) doe1_norm_ls_04(:,5:7)];

% Evaluating the Limit State 1
y_ls_04  = constraint4c(x_ls_04);

% Computing active subspace
k_ls_04=2;
m_c04=6;
npoints_active_susbpace_ls_04 = 10*m_c04;
M_ls_04 = min(ceil(10*m_c04*log(m_c04)),npoints_active_susbpace_ls_04-1);
p_ls_04 =floor(10*m_c04)-1;
[b_ls_04_llm,C_ls_04_llm,W_ls_04_llm,Ev_ls_04_llm]=llm(x_norm_c04,y_ls_04,M_ls_04,p_ls_04);

% Define active Subspace dimension
pv_ls_04 = cumsum(diag(Ev_ls_04_llm))/sum(diag(Ev_ls_04_llm));
n_ls_04  = active_subspace_dimension(diag(Ev_ls_04_llm));

% Compute active subspace and active variable
W1_ls_04_llm = W_ls_04_llm(:,1:n_ls_04);
t1_ls_04 = x_norm_c04*W1_ls_04_llm;

% Doe in the active subspace
t1_ls_04_min = -diag(sign(W1_ls_04_llm)'*W1_ls_04_llm);
t1_ls_04_max = diag(sign(W1_ls_04_llm)'*W1_ls_04_llm);
t1_ls_04_range = t1_ls_04_max'-t1_ls_04_min';
N_ls_04 = 10*n_ls_04+1;
% t1_ls_04_new = t1_ls_04_min'+((t1_ls_04_range).*lhsdesign(N_ls_04,n_ls_04));

% Mapping to original space
sp_c04 = [0 0 0 0 0 0];
lbx_norm_c04 = normalizeuv(lbx_c04,lbx_c04,ubx_c04);
ubx_norm_c04 = normalizeuv(ubx_c04,lbx_c04,ubx_c04);
xmin_c04 = zeros(N_ls_04,m_c04);
fval_c04 = zeros(N_ls_04,1);
exitflag_c04 = zeros(N_ls_04,1);

for i=1:N_ls_04
func_c04 =  @(x_opt) (x_opt*zeros(m_c04,1));
[xmin_c04(i,:),fval_c04(i,1),exitflag_c04(i,1)] = fmincon(func_c04,sp_c04,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x_opt)nonlconstage1(x_opt,W1_ls_04_llm,t1_ls_04_new(i,:)));
end

npoints_as_c04=10*n_ls_04+1;
x_norm_c04_org = xmin_c04(exitflag_c04==1,:);
x_norm_c04_org_chose = x_norm_c04_org(1:npoints_as_c04,:);
x_c04_org_new = denormalizeuv(x_norm_c04_org_chose,lbx_c04,ubx_c04);

% Evaluate PSF
x_c04_new =  [x_c04_org_new(:,1:3) zeros(npoints_as,1) x_c04_org_new(:,4:6)];
[psf_c04,pf_c04,re_c04,beta_c04] = mcspsfconstraint4c(x_c04_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c04 = x_norm_c04_org_chose*W1_ls_04_llm;
y_c04 = psf_c04;
[srgt_RBF_c04, PRESSRMS_RBF_c04, eXV_RBF_c04, srgtOPT_RBF_c04, Y_hat_RBF_c04, R2_pred_RBF_c04] = metamodel_RBF(t1_c04,y_c04);
[srgt_KRG_c04, PRESSRMS_KRG_c04, eXV_KRG_c04, srgtOPT_KRG_c04, Y_hat_KRG_c04, pred_var_KRG_c04, R2_pred_KRG_c04] = metamodel_KRG(t1_c04,y_c04);
[srgt_PRS_c04, PRESSRMS_PRS_c04, eXV_PRS_c04, srgtOPT_PRS_c04, Y_hat_PRS_c04, pred_var_PRS_c04, R2_pred_PRS_c04] = metamodel_PRS(t1_c04,y_c04);
[srgt_WAS_c04, PRESSRMS_WAS_c04, eXV_WAS_c04, srgtOPT_WAS_c04, Y_hat_WAS_c04, pred_var_WAS_c04, R2_pred_WAS_c04] = metamodel_WAS(t1_c04,y_c04);

%% Constraint 5

% Bounds
lbx_c05 = [lbx(2:3) lbx(7)];
ubx_c05 = [ubx(2:3) ubx(7)];

% Sample from the input space
n_var_c05 = 3;
npoints_c05=10*n_var_c05;
% scalenewDoe1_c05 =zeros(npoints_c05,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c05(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c05,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c05(:,i+n_design_var) = repmat(mup(i),npoints_c05,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_05 = scalenewDoe1_c05;
doe1_norm_ls_05 = normalizeuv(doe1_ls_05(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_05 = [doe1_ls_05(:,2:3) doe1_ls_05(:,7) doe1_ls_05(:,8:10)];
x_c05 = [doe1_ls_05(:,2:3) doe1_ls_05(:,7)];
x_norm_c05 = [doe1_norm_ls_05(:,2:3) doe1_norm_ls_05(:,7)];

% Evaluating the Limit State 1
y_ls_05  = constraint5c(x_ls_05);

% Computing active subspace
k_ls_05=2;
m_c05=3;
npoints_active_susbpace_ls_05 = 10*m_c05;
M_ls_05 = min(ceil(10*m_c05*log(m_c05)),npoints_active_susbpace_ls_05-1);
p_ls_05 =floor(10*m_c05)-1;
[b_ls_05_llm,C_ls_05_llm,W_ls_05_llm,Ev_ls_05_llm]=llm(x_norm_c05,y_ls_05,M_ls_05,p_ls_05);

% Define active Subspace dimension
pv_ls_05 = cumsum(diag(Ev_ls_05_llm))/sum(diag(Ev_ls_05_llm));
n_ls_05  = active_subspace_dimension(diag(Ev_ls_05_llm));

% Compute active subspace and active variable
W1_ls_05_llm = W_ls_05_llm(:,1:n_ls_05);
t1_ls_05 = x_norm_c05*W1_ls_05_llm;

% Doe in the active subspace
t1_ls_05_min = -diag(sign(W1_ls_05_llm)'*W1_ls_05_llm);
t1_ls_05_max = diag(sign(W1_ls_05_llm)'*W1_ls_05_llm);
t1_ls_05_range = t1_ls_05_max'-t1_ls_05_min';
N_ls_05 = 10*n_ls_05+1;
% t1_ls_05_new = t1_ls_05_min'+((t1_ls_05_range).*lhsdesign(N_ls_05,n_ls_05));

% Mapping to original space
sp_c05 = [0 0 0];
lbx_norm_c05 = normalizeuv(lbx_c05,lbx_c05,ubx_c05);
ubx_norm_c05 = normalizeuv(ubx_c05,lbx_c05,ubx_c05);
xmin_c05 = zeros(N_ls_05,m_c05);
fval_c05 = zeros(N_ls_05,1);
exitflag_c05 = zeros(N_ls_05,1);

for i=1:N_ls_05
func_c05 =  @(x_opt) (x_opt*zeros(m_c05,1));
[xmin_c05(i,:),fval_c05(i,1),exitflag_c05(i,1)] = fmincon(func_c05,sp_c05,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x_opt)nonlconstage1(x_opt,W1_ls_05_llm,t1_ls_05_new(i,:)));
end

npoints_as_c05 = 10*n_ls_05+1;
x_norm_c05_org = xmin_c05(exitflag_c05==1,:);
x_norm_c05_org_chose = x_norm_c05_org(1:npoints_as_c05,:);
x_c05_org_new = denormalizeuv(x_norm_c05_org_chose,lbx_c05,ubx_c05);

% Evaluate PSF
x_c05_new =  [zeros(npoints_as,1) x_c05_org_new(:,1:2) zeros(npoints_as,3) x_c05_org_new(:,3)];
[psf_c05,pf_c05,re_c05,beta_c05] = mcspsfconstraint5c(x_c05_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c05 = x_norm_c05_org_chose*W1_ls_05_llm;
y_c05 = psf_c05;
[srgt_RBF_c05, PRESSRMS_RBF_c05, eXV_RBF_c05, srgtOPT_RBF_c05, Y_hat_RBF_c05, R2_pred_RBF_c05] = metamodel_RBF(t1_c05,y_c05);
[srgt_KRG_c05, PRESSRMS_KRG_c05, eXV_KRG_c05, srgtOPT_KRG_c05, Y_hat_KRG_c05, pred_var_KRG_c05, R2_pred_KRG_c05] = metamodel_KRG(t1_c05,y_c05);
[srgt_PRS_c05, PRESSRMS_PRS_c05, eXV_PRS_c05, srgtOPT_PRS_c05, Y_hat_PRS_c05, pred_var_PRS_c05, R2_pred_PRS_c05] = metamodel_PRS(t1_c05,y_c05);
[srgt_WAS_c05, PRESSRMS_WAS_c05, eXV_WAS_c05, srgtOPT_WAS_c05, Y_hat_WAS_c05, pred_var_WAS_c05, R2_pred_WAS_c05] = metamodel_WAS(t1_c05,y_c05);

%% Constraint 6

% Bounds
lbx_c06 = [lbx(1:3) lbx(5:7)];
ubx_c06 = [ubx(1:3) ubx(5:7)];

% Sample from the input space
n_var_c06 = 6;
npoints_c06=10*n_var_c06;
% scalenewDoe1_c06 =zeros(npoints_c06,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c06(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c06,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c06(:,i+n_design_var) = repmat(mup(i),npoints_c06,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_06 = scalenewDoe1_c06;
doe1_norm_ls_06 = normalizeuv(doe1_ls_06(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_06 = [doe1_ls_06(:,1:3) doe1_ls_06(:,5:7) doe1_ls_06(:,8:11)];
x_c06 = [doe1_ls_06(:,1:3) doe1_ls_06(:,5:7)];
x_norm_c06 = [doe1_norm_ls_06(:,1:3) doe1_norm_ls_06(:,5:7)];

% Evaluating the Limit State 1
y_ls_06  = constraint6c(x_ls_06);

% Computing active subspace
k_ls_06=2;
m_c06=6;
npoints_active_susbpace_ls_06 = 10*m_c06;
M_ls_06 = min(ceil(10*m_c06*log(m_c06)),npoints_active_susbpace_ls_06-1);
p_ls_06 =floor(10*m_c06)-1;
[b_ls_06_llm,C_ls_06_llm,W_ls_06_llm,Ev_ls_06_llm]=llm(x_norm_c06,y_ls_06,M_ls_06,p_ls_06);

% Define active Subspace dimension
pv_ls_06 = cumsum(diag(Ev_ls_06_llm))/sum(diag(Ev_ls_06_llm));
n_ls_06  = active_subspace_dimension(diag(Ev_ls_05_llm));

% Compute active subspace and active variable
W1_ls_06_llm = W_ls_06_llm(:,1:n_ls_06);
t1_ls_06 = x_norm_c06*W1_ls_06_llm;

% Doe in the active subspace
t1_ls_06_min = -diag(sign(W1_ls_06_llm)'*W1_ls_06_llm);
t1_ls_06_max = diag(sign(W1_ls_06_llm)'*W1_ls_06_llm);
t1_ls_06_range = t1_ls_06_max'-t1_ls_06_min';
N_ls_06 = 10*n_ls_06+1;
% t1_ls_06_new = t1_ls_06_min'+((t1_ls_06_range).*lhsdesign(N_ls_06,n_ls_06));

% Mapping to original space
sp_c06 = [0 0 0 0 0 0];
lbx_norm_c06 = normalizeuv(lbx_c06,lbx_c06,ubx_c06);
ubx_norm_c06 = normalizeuv(ubx_c06,lbx_c06,ubx_c06);
xmin_c06 = zeros(N_ls_06,m_c06);
fval_c06 = zeros(N_ls_06,1);
exitflag_c06 = zeros(N_ls_06,1);

for i=1:N_ls_06
func_c06 =  @(x_opt) (x_opt*zeros(m_c06,1));
[xmin_c06(i,:),fval_c06(i,1),exitflag_c06(i,1)] = fmincon(func_c06,sp_c06,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x_opt)nonlconstage1(x_opt,W1_ls_06_llm,t1_ls_06_new(i,:)));
end

npoints_as_c06 = 10*n_ls_06+1;
x_norm_c06_org = xmin_c06(exitflag_c06==1,:);
x_norm_c06_org_chose = x_norm_c06_org(1:npoints_as_c06,:);
% x_c06_org_new = denormalizeuv(x_norm_c06_org_chose,lbx_c06,ubx_c06);

% Evaluate PSF
x_c06_new =  [x_c06_org_new(:,1:3) zeros(npoints_as,1) x_c06_org_new(:,4:6)];
[psf_c06,pf_c06,re_c06,beta_c06] = mcspsfconstraint6c(x_c06_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c06 = x_norm_c06_org_chose*W1_ls_06_llm;
y_c06 = psf_c06;
[srgt_RBF_c06, PRESSRMS_RBF_c06, eXV_RBF_c06, srgtOPT_RBF_c06, Y_hat_RBF_c06, R2_pred_RBF_c06] = metamodel_RBF(t1_c06,y_c06);
[srgt_KRG_c06, PRESSRMS_KRG_c06, eXV_KRG_c06, srgtOPT_KRG_c06, Y_hat_KRG_c06, pred_var_KRG_c06, R2_pred_KRG_c06] = metamodel_KRG(t1_c06,y_c06);
[srgt_PRS_c06, PRESSRMS_PRS_c06, eXV_PRS_c06, srgtOPT_PRS_c06, Y_hat_PRS_c06, pred_var_PRS_c06, R2_pred_PRS_c06] = metamodel_PRS(t1_c06,y_c06);
[srgt_WAS_c06, PRESSRMS_WAS_c06, eXV_WAS_c06, srgtOPT_WAS_c06, Y_hat_WAS_c06, pred_var_WAS_c06, R2_pred_WAS_c06] = metamodel_WAS(t1_c06,y_c06);

%% Constraint 7

% Bounds
lbx_c07 = [lbx(1:3) lbx(5:7)];
ubx_c07 = [ubx(1:3) ubx(5:7)];

% Sample from the input space
n_var_c07 = 6;
npoints_c07=10*n_var_c07;
% scalenewDoe1_c07 =zeros(npoints_c07,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c07(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c07,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c07(:,i+n_design_var) = repmat(mup(i),npoints_c07,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_07 = scalenewDoe1_c07;
doe1_norm_ls_07 = normalizeuv(doe1_ls_07(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_07 = [doe1_ls_07(:,1:3) doe1_ls_07(:,5:7) doe1_ls_07(:,8:11)];
x_c07 = [doe1_ls_07(:,1:3) doe1_ls_07(:,5:7)];
x_norm_c07 = [doe1_norm_ls_07(:,1:3) doe1_norm_ls_07(:,5:7)];

% Evaluating the Limit State 1
y_ls_07  = constraint7c(x_ls_07);

% Computing active subspace
k_ls_07=2;
m_c07=6;
npoints_active_susbpace_ls_07 = 10*m_c07;
M_ls_07 = min(ceil(10*m_c07*log(m_c07)),npoints_active_susbpace_ls_07-1);
p_ls_07 = floor(10*m_c07)-1;
[b_ls_07_llm,C_ls_07_llm,W_ls_07_llm,Ev_ls_07_llm]=llm(x_norm_c07,y_ls_07,M_ls_07,p_ls_07);

% Define active Subspace dimension
pv_ls_07 = cumsum(diag(Ev_ls_07_llm))/sum(diag(Ev_ls_07_llm));
n_ls_07  = active_subspace_dimension(diag(Ev_ls_07_llm));

% Compute active subspace and active variable
W1_ls_07_llm = W_ls_07_llm(:,1:n_ls_07);
t1_ls_07 = x_norm_c07*W1_ls_07_llm;

% Doe in the active subspace
t1_ls_07_min = -diag(sign(W1_ls_07_llm)'*W1_ls_07_llm);
t1_ls_07_max = diag(sign(W1_ls_07_llm)'*W1_ls_07_llm);
t1_ls_07_range = t1_ls_07_max'-t1_ls_07_min';
N_ls_07 = 10*n_ls_07+1;
% t1_ls_07_new = t1_ls_07_min'+((t1_ls_07_range).*lhsdesign(N_ls_07,n_ls_07));

% Mapping to original space
sp_c07 = [0 0 0 0 0 0];
lbx_norm_c07 = normalizeuv(lbx_c07,lbx_c07,ubx_c07);
ubx_norm_c07 = normalizeuv(ubx_c07,lbx_c07,ubx_c07);
xmin_c07 = zeros(N_ls_07,m_c07);
fval_c07 = zeros(N_ls_07,1);
exitflag_c07 = zeros(N_ls_07,1);

for i=1:N_ls_07
func_c07 =  @(x_opt) (x_opt*zeros(m_c07,1));
[xmin_c07(i,:),fval_c07(i,1),exitflag_c07(i,1)] = fmincon(func_c07,sp_c07,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x_opt)nonlconstage1(x_opt,W1_ls_07_llm,t1_ls_07_new(i,:)));
end

npoints_as_c07=10*n_ls_07+1;
x_norm_c07_org = xmin_c07(exitflag_c07==1,:);
x_norm_c07_org_chose = x_norm_c07_org(1:npoints_as_c07,:);
% x_c07_org_new = denormalizeuv(x_norm_c07_org_chose,lbx_c07,ubx_c07);

% Evaluate PSF
x_c07_new =  [x_c07_org_new(:,1:3) zeros(npoints_as,1) x_c07_org_new(:,4:6)];
[psf_c07,pf_c07,re_c07,beta_c07] = mcspsfconstraint7c(x_c07_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c07 = x_norm_c07_org_chose*W1_ls_07_llm;
y_c07 = psf_c07;
[srgt_RBF_c07, PRESSRMS_RBF_c07, eXV_RBF_c07, srgtOPT_RBF_c07, Y_hat_RBF_c07, R2_pred_RBF_c07] = metamodel_RBF(t1_c07,y_c07);
[srgt_KRG_c07, PRESSRMS_KRG_c07, eXV_KRG_c07, srgtOPT_KRG_c07, Y_hat_KRG_c07, pred_var_KRG_c07, R2_pred_KRG_c07] = metamodel_KRG(t1_c07,y_c07);
[srgt_PRS_c07, PRESSRMS_PRS_c07, eXV_PRS_c07, srgtOPT_PRS_c07, Y_hat_PRS_c07, pred_var_PRS_c07, R2_pred_PRS_c07] = metamodel_PRS(t1_c07,y_c07);
[srgt_WAS_c07, PRESSRMS_WAS_c07, eXV_WAS_c07, srgtOPT_WAS_c07, Y_hat_WAS_c07, pred_var_WAS_c07, R2_pred_WAS_c07] = metamodel_WAS(t1_c07,y_c07);

%% Constraint 8

% Bounds
lbx_c08 = [lbx(2:4) lbx(6)];
ubx_c08 = [ubx(2:4) ubx(6)];

% Sample from the input space
n_var_c08 = 4;
npoints_c08=10*n_var_c08;
% scalenewDoe1_c08 =zeros(npoints_c08,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c08(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c08,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c08(:,i+n_design_var) = repmat(mup(i),npoints_c08,1);
% end
%Normalizing the sample
doe1_ls_08 = scalenewDoe1_c08;
doe1_norm_ls_08 = normalizeuv(scalenewDoe1_c08(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_08 = [doe1_ls_08(:,2:4) doe1_ls_08(:,6) doe1_ls_08(:,10:11)];
x_c08 = [doe1_ls_08(:,2:4) doe1_ls_08(:,6)];
x_norm_c08 = [doe1_norm_ls_08(:,2:4) doe1_norm_ls_08(:,6)];

% Evaluating the Limit State 1
y_ls_08  = constraint8c(x_ls_08);

% Computing active subspace
k_ls_08=2;
m_c08=4;
npoints_active_susbpace_ls_08 = 10*m_c08;
M_ls_08 = min(ceil(10*m_c08*log(m_c08)),npoints_active_susbpace_ls_08-1);
p_ls_08 =floor(10*m_c08)-1;
[b_ls_08_llm,C_ls_08_llm,W_ls_08_llm,Ev_ls_08_llm]=llm(x_norm_c08,y_ls_08,M_ls_08,p_ls_08);

% Define active Subspace dimension
pv_ls_08 = cumsum(diag(Ev_ls_08_llm))/sum(diag(Ev_ls_08_llm));
n_ls_08  = active_subspace_dimension(diag(Ev_ls_08_llm));

% Compute active subspace and active variable
W1_ls_08_llm = W_ls_08_llm(:,1:n_ls_08);
t1_ls_08 = x_norm_c08*W1_ls_08_llm;

% Doe in the active subspace
t1_ls_08_min = -diag(sign(W1_ls_08_llm)'*W1_ls_08_llm);
t1_ls_08_max = diag(sign(W1_ls_08_llm)'*W1_ls_08_llm);
t1_ls_08_range = t1_ls_08_max'-t1_ls_08_min';
N_ls_08 = 10*n_ls_08+1;
% t1_ls_08_new = t1_ls_08_min'+((t1_ls_08_range).*lhsdesign(N_ls_08,n_ls_08));

% Mapping to original space
sp_c08 = [0 0 0 0];
lbx_norm_c08 = normalizeuv(lbx_c08,lbx_c08,ubx_c08);
ubx_norm_c08 = normalizeuv(ubx_c08,lbx_c08,ubx_c08);
xmin_c08 = zeros(N_ls_08,m_c08);
fval_c08 = zeros(N_ls_08,1);
exitflag_c08 = zeros(N_ls_08,1);

for i=1:N_ls_08
func_c08 =  @(x_opt) (x_opt*zeros(m_c08,1));
[xmin_c08(i,:),fval_c08(i,1),exitflag_c08(i,1)] = fmincon(func_c08,sp_c08,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x_opt)nonlconstage1(x_opt,W1_ls_08_llm,t1_ls_08_new(i,:)));
end

npoints_as_c08=10*n_ls_08+1;
x_norm_c08_org = xmin_c08(exitflag_c08==1,:);
x_norm_c08_org_chose_1 = x_norm_c08_org(1:npoints_as_c08,:);
[t_bnd_c08,x_bnd_c08] = zonotope_vertices(W1_ls_08_llm,1e4,1e5);
x_norm_c08_org_chose = x_norm_c08_org_chose_1;
% x_c08_org_new = denormalizeuv(x_norm_c08_org_chose,lbx_c08,ubx_c08);

% Evaluate PSF
x_c08_new = [zeros(npoints_as,1) x_c08_org_new(:,1:3) zeros(npoints_as,1) x_c08_org_new(:,4) zeros(npoints_as,1)];
[psf_c08,pf_c08,re_c08,beta_c08] = mcspsfconstraint8c(x_c08_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

% Build the metamodel in active subspace of Limit State
t1_c08 = x_norm_c08_org_chose*W1_ls_08_llm;
y_c08 = psf_c08;
[srgt_RBF_c08, PRESSRMS_RBF_c08, eXV_RBF_c08, srgtOPT_RBF_c08, Y_hat_RBF_c08, R2_pred_RBF_c08] = metamodel_RBF(t1_c08,y_c08);
[srgt_KRG_c08, PRESSRMS_KRG_c08, eXV_KRG_c08, srgtOPT_KRG_c08, Y_hat_KRG_c08, pred_var_KRG_c08, R2_pred_KRG_c08] = metamodel_KRG(t1_c08,y_c08);
[srgt_PRS_c08, PRESSRMS_PRS_c08, eXV_PRS_c08, srgtOPT_PRS_c08, Y_hat_PRS_c08, pred_var_PRS_c08, R2_pred_PRS_c08] = metamodel_PRS(t1_c08,y_c08);
[srgt_WAS_c08, PRESSRMS_WAS_c08, eXV_WAS_c08, srgtOPT_WAS_c08, Y_hat_WAS_c08, pred_var_WAS_c08, R2_pred_WAS_c08] = metamodel_WAS(t1_c08,y_c08);

%% Constraint 9

% Bounds
lbx_c09 = [lbx(1:4) lbx(6)];
ubx_c09 = [ubx(1:4) ubx(6)];

% Sample from the input space
n_var_c09 = 5;
npoints_c09=10*n_var_c09;
% scalenewDoe1_c09 =zeros(npoints_c09,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c09(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c09,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c09(:,i+n_design_var) = repmat(mup(i),npoints_c09,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_09 = scalenewDoe1_c09;
doe1_norm_ls_09 = normalizeuv(doe1_ls_09(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 1
x_ls_09 = [doe1_ls_09(:,1:4) doe1_ls_09(:,6) doe1_ls_09(:,8) doe1_ls_09(:,10)];
x_c09 = [doe1_ls_09(:,1:4) doe1_ls_09(:,6)];
x_norm_c09 = [doe1_norm_ls_09(:,1:4) doe1_norm_ls_09(:,6)];

% Evaluating the Limit State 1
y_ls_09  = constraint9c(x_ls_09);

% Computing active subspace
k_ls_09=2;
m_c09=5;
npoints_active_susbpace_ls_09 = 10*m_c09;
M_ls_09 = min(ceil(10*m_c09*log(m_c09)),npoints_active_susbpace_ls_09-1);
p_ls_09 = floor(10*m_c09)-1;
[b_ls_09_llm,C_ls_09_llm,W_ls_09_llm,Ev_ls_09_llm]=llm(x_norm_c09,y_ls_09,M_ls_09,p_ls_09);

% Define active Subspace dimension
pv_ls_09 = cumsum(diag(Ev_ls_09_llm))/sum(diag(Ev_ls_09_llm));
n_ls_09  = active_subspace_dimension(diag(Ev_ls_09_llm));

% Compute active subspace and active variable
W1_ls_09_llm = W_ls_09_llm(:,1:n_ls_09);
t1_ls_09 = x_norm_c09*W1_ls_09_llm;

% Doe in the active subspace
t1_ls_09_min = -diag(sign(W1_ls_09_llm)'*W1_ls_09_llm);
t1_ls_09_max = diag(sign(W1_ls_09_llm)'*W1_ls_09_llm);
t1_ls_09_range = t1_ls_09_max'-t1_ls_09_min';
N_ls_09 = 10*n_ls_09+1;
% t1_ls_09_new = t1_ls_09_min'+((t1_ls_09_range).*lhsdesign(N_ls_09,n_ls_09));

% Mapping to original space
sp_c09 = [0 0 0 0 0];
lbx_norm_c09 = normalizeuv(lbx_c09,lbx_c09,ubx_c09);
ubx_norm_c09 = normalizeuv(ubx_c09,lbx_c09,ubx_c09);
xmin_c09 = zeros(N_ls_09,m_c09);
fval_c09 = zeros(N_ls_09,1);
exitflag_c09 = zeros(N_ls_09,1);

for i=1:N_ls_09
func_c09 =  @(x_opt) (x_opt*zeros(m_c09,1));
[xmin_c09(i,:),fval_c09(i,1),exitflag_c09(i,1)] = fmincon(func_c09,sp_c09,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x_opt)nonlconstage1(x_opt,W1_ls_09_llm,t1_ls_09_new(i,:)));
end

npoints_as_c09=10*n_ls_09+1;
x_norm_c09_org = xmin_c09(exitflag_c09==1,:);
x_norm_c09_org_chose = x_norm_c09_org(1:npoints_as_c09,:);
% x_c09_org_new = denormalizeuv(x_norm_c09_org_chose,lbx_c09,ubx_c09);

% Evaluate PSF
x_c09_new = [x_c09_org_new(:,1:4) zeros(npoints_as,1) x_c09_org_new(:,5) zeros(npoints_as,1)];
[psf_c09,pf_c09,re_c09,beta_c09] = mcspsfconstraint9c(x_c09_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c09 = x_norm_c09_org_chose*W1_ls_09_llm;
y_c09 = psf_c09;
[srgt_RBF_c09, PRESSRMS_RBF_c09, eXV_RBF_c09, srgtOPT_RBF_c09, Y_hat_RBF_c09, R2_pred_RBF_c09] = metamodel_RBF(t1_c09,y_c09);
[srgt_KRG_c09, PRESSRMS_KRG_c09, eXV_KRG_c09, srgtOPT_KRG_c09, Y_hat_KRG_c09, pred_var_KRG_c09, R2_pred_KRG_c09] = metamodel_KRG(t1_c09,y_c09);
[srgt_PRS_c09, PRESSRMS_PRS_c09, eXV_PRS_c09, srgtOPT_PRS_c09, Y_hat_PRS_c09, pred_var_PRS_c09, R2_pred_PRS_c09] = metamodel_PRS(t1_c09,y_c09);
[srgt_WAS_c09, PRESSRMS_WAS_c09, eXV_WAS_c09, srgtOPT_WAS_c09, Y_hat_WAS_c09, pred_var_WAS_c09, R2_pred_WAS_c09] = metamodel_WAS(t1_c09,y_c09);

%% Constraint 10

% Bounds
lbx_c10 = [lbx(3) lbx(5:7)];
ubx_c10 = [ubx(3) ubx(5:7)];

% Sample from the input space
n_var_c10 = 4;
npoints_c10=10*n_var_c10;
% scalenewDoe1_c10 =zeros(npoints_c10,n_var);
% 
% for i=1:n_design_var
% scalenewDoe1_c10(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_c10,1));
% end
% for i=1:n_random_var
%     scalenewDoe1_c10(:,i+n_design_var) = repmat(mup(i),npoints_c10,1);
% end
% Sampling points from sample mean and standard deviation
doe1_ls_10 = scalenewDoe1_c10;
doe1_norm_ls_10 = normalizeuv(doe1_ls_10(:,1:n_design_var),lbx,ubx);

% Input Doe for Constraint 10
x_ls_10 = [doe1_ls_10(:,3) doe1_ls_10(:,5:7) doe1_ls_10(:,9:11)];
x_c10 = [doe1_ls_10(:,3) doe1_ls_10(:,5:7)];
x_norm_c10 = [doe1_norm_ls_10(:,3) doe1_norm_ls_10(:,5:7)];

% Evaluating the Limit State 10
y_ls_10  = constraint10c(x_ls_10);

% Computing active subspace
k_ls_10=2;
m_c10=4;
npoints_active_susbpace_ls_10 = 10*m_c10;
M_ls_10 = min(ceil(10*m_c10*log(m_c10)),npoints_active_susbpace_ls_10-1);
p_ls_10 = floor(10*m_c10)-1;
[b_ls_10_llm,C_ls_10_llm,W_ls_10_llm,Ev_ls_10_llm]=llm(x_norm_c10,y_ls_10,M_ls_10,p_ls_10);

% Define active Subspace dimension
pv_ls_10 = cumsum(diag(Ev_ls_10_llm))/sum(diag(Ev_ls_10_llm));
n_ls_10  = active_subspace_dimension(diag(Ev_ls_10_llm));

% Compute active subspace and active variable
W1_ls_10_llm = W_ls_10_llm(:,1:n_ls_10);
t1_ls_10 = x_norm_c10*W1_ls_10_llm;

% Doe in the active subspace
t1_ls_10_min = -diag(sign(W1_ls_10_llm)'*W1_ls_10_llm);
t1_ls_10_max = diag(sign(W1_ls_10_llm)'*W1_ls_10_llm);
t1_ls_10_range = t1_ls_10_max'-t1_ls_10_min';
N_ls_10 = 10*n_ls_10+1;
% t1_ls_10_new = t1_ls_10_min'+((t1_ls_10_range).*lhsdesign(N_ls_10,n_ls_10));

% Mapping to original space
sp_c10 = [0 0 0 0];
lbx_norm_c10 = normalizeuv(lbx_c10,lbx_c10,ubx_c10);
ubx_norm_c10 = normalizeuv(ubx_c10,lbx_c10,ubx_c10);
xmin_c10 = zeros(N_ls_10,m_c10);
fval_c10 = zeros(N_ls_10,1);
exitflag_c10 = zeros(N_ls_10,1);

for i=1:N_ls_10
func_c10 =  @(x_opt) (x_opt*zeros(m_c10,1));
[xmin_c10(i,:),fval_c10(i,1),exitflag_c10(i,1)] = fmincon(func_c10,sp_c10,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x_opt)nonlconstage1(x_opt,W1_ls_10_llm,t1_ls_10_new(i,:)));
end

npoints_as_c10=10*n_ls_10+1;
x_norm_c10_org = xmin_c10(exitflag_c10==1,:);
x_norm_c10_org_chose = x_norm_c10_org(1:npoints_as_c10,:);
% x_c10_org_new = denormalizeuv(x_norm_c10_org_chose,lbx_c10,ubx_c10);

% Evaluate PSF
x_c10_new = [zeros(npoints_as,2) x_c10_org_new(:,1) zeros(npoints_as,1) x_c10_org_new(:,2:4)];
[psf_c10,pf_c10,re_c10,beta_c10] = mcspsfconstraint10c(x_c10_new,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

%Build the metamodel in active subspace of Limit State
t1_c10 = x_norm_c10_org_chose*W1_ls_10_llm;
y_c10 = psf_c10;
[srgt_RBF_c10, PRESSRMS_RBF_c10, eXV_RBF_c10, srgtOPT_RBF_c10, Y_hat_RBF_c10, R2_pred_RBF_c10] = metamodel_RBF(t1_c10,y_c10);
[srgt_KRG_c10, PRESSRMS_KRG_c10, eXV_KRG_c10, srgtOPT_KRG_c10, Y_hat_KRG_c10, pred_var_KRG_c10, R2_pred_KRG_c10] = metamodel_KRG(t1_c10,y_c10);
[srgt_PRS_c10, PRESSRMS_PRS_c10, eXV_PRS_c10, srgtOPT_PRS_c10, Y_hat_PRS_c10, pred_var_PRS_c10, R2_pred_PRS_c10] = metamodel_PRS(t1_c10,y_c10);
[srgt_WAS_c10, PRESSRMS_WAS_c10, eXV_WAS_c10, srgtOPT_WAS_c10, Y_hat_WAS_c10, pred_var_WAS_c10, R2_pred_WAS_c10] = metamodel_WAS(t1_c10,y_c10);

%% Objective Function

% Bounds
lbx_obj = [lbx(1:5) lbx(7)];
ubx_obj = [ubx(1:5) ubx(7)];

% Sample from the input space
n_var_obj = 6;
npoints_obj=10*n_var_obj;
% scalenewDoe1_obj =zeros(npoints_obj,n_design_var);
% 
% for i=1:n_design_var
% scalenewDoe1_obj(:,i) = lbx(i)+((ubx(i)-lbx(i)).*lhsdesign(npoints_obj,1));
% end

% Sampling points from sample mean and standard deviation
doe1_obj = scalenewDoe1_obj;
doe1_norm_obj = normalizeuv(doe1_obj,lbx,ubx);

% Input Doe for Constraint 1
x_obj = [doe1_obj(:,1:5) doe1_obj(:,7)];
x_norm_obj = [doe1_norm_obj(:,1:5) doe1_norm_obj(:,7)];

% Evaluating the Limit State 1
y_obj  = Weightc(x_obj);

% Computing active subspace
k_obj=2;
m_obj=6;
npoints_active_susbpace_obj = 10*m_obj;
M_obj = min(ceil(10*m_obj*log(m_obj)),npoints_active_susbpace_obj-1);
p_obj =floor(10*m_obj)-1;
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
sp_obj = [0 0 0 0 0 0];
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
% x_obj_org_new = denormalizeuv(x_norm_obj_org_chose,lbx_obj,ubx_obj);

% Evaluate PSF
x_obj_new =  x_obj_org_new;
y_obj_new = Weightc(x_obj_new);

%Build the metamodel in active subspace of Limit State
t1_obj_1 = x_norm_obj_org_chose*W1_obj_llm;
y_obj_1 = y_obj_new;
[srgt_RBF_obj_1, PRESSRMS_RBF_obj_1, eXV_RBF_obj_1, srgtOPT_RBF_obj_1, Y_hat_RBF_obj_1, R2_pred_RBF_obj_1] = metamodel_RBF(t1_obj_1,y_obj_1);
[srgt_KRG_obj_1, PRESSRMS_KRG_obj_1, eXV_KRG_obj_1, srgtOPT_KRG_obj_1, Y_hat_KRG_obj_1, pred_var_KRG_obj_1, R2_pred_KRG_obj_1] = metamodel_KRG(t1_obj_1,y_obj_1);
[srgt_PRS_obj_1, PRESSRMS_PRS_obj_1, eXV_PRS_obj_1, srgtOPT_PRS_obj_1, Y_hat_PRS_obj_1, pred_var_PRS_obj_1, R2_pred_PRS_obj_1] = metamodel_PRS(t1_obj_1,y_obj_1);
[srgt_WAS_obj_1, PRESSRMS_WAS_obj_1, eXV_WAS_obj_1, srgtOPT_WAS_obj_1, Y_hat_WAS_obj_1, pred_var_WAS_obj_1, R2_pred_WAS_obj_1] = metamodel_WAS(t1_obj_1,y_obj_1);

%% Range of Response

psf_c01_range = max(psf_c01)-min(psf_c01);
psf_c02_range = max(psf_c02)-min(psf_c02);
psf_c03_range = max(psf_c03)-min(psf_c03);
psf_c04_range = max(psf_c04)-min(psf_c04);
psf_c05_range = max(psf_c05)-min(psf_c05);
psf_c06_range = max(psf_c06)-min(psf_c06);
psf_c07_range = max(psf_c07)-min(psf_c07);
psf_c08_range = max(psf_c08)-min(psf_c08);
psf_c09_range = max(psf_c09)-min(psf_c09);
psf_c10_range = max(psf_c10)-min(psf_c10);
obj_range = max(y_obj_new)-min(y_obj_new);
psf_range = [psf_c01_range;psf_c02_range;psf_c03_range;psf_c04_range;psf_c05_range;psf_c06_range;psf_c07_range;psf_c08_range;psf_c09_range;psf_c10_range];

