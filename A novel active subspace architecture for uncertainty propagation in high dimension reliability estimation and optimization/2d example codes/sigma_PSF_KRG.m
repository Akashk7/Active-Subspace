clc;close all;
%% Constraint 1

%Finding t at which PSF is one
sp_t_c01_KRG = zeros(1,n_ls_01);
func_t_c01_KRG =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_KRG, fval_t_c01_KRG, exitflag_t_c01_KRG] = fmincon(func_t_c01_KRG,sp_t_c01_KRG,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c01,@srgtsKRGEvaluate,1));

%Validation
func_c01_KRG_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_KRG_1, fval_c01_KRG_1, exitflag_c01_KRG_1] = fmincon(func_c01_KRG_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_KRG));
xopt_c01_KRG_1_denorm = denormalizeuv(xopt_c01_KRG_1,lbx_c01,ubx_c01);
xopt_c01_KRG_1_org = xopt_c01_KRG_1_denorm;
[psf_c01_KRG_1,pf_c01_KRG_1,re_c01_KRG_1,beta_c01_KRG_1] = mcspsfconstraint1c(xopt_c01_KRG_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c01_KRG_1_mm = srgtsKRGEvaluate(t_c01_KRG,srgt_KRG_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);
t2_ls_01_llm = xopt_c01_KRG_1*W2_ls_01_llm;

sp_c01_KRG = zeros(1,m_c01);
fmin_c01_KRG_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_KRG, fvalmin_1_c01_KRG, exitflagmin_1_c01_KRG]=fmincon(fmin_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm));

fmax_c01_KRG_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_KRG, fvalmax_1_c01_KRG, exitflagmax_1_c01_KRG]=fmincon(fmax_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm));

x_bnd_c01_KRG_1 = [xmin_1_c01_KRG;xmax_1_c01_KRG];

ind_c01_KRG = 2*(m_c01-n_ls_01);
x_bnd_c01_KRG_1_org = denormalizeuv(x_bnd_c01_KRG_1,lbx,ubx);
x_bnd_c01_KRG_2_1 = [t_c01_KRG (t2_ls_01_llm+fvalmin_1_c01_KRG)/2]*(W_ls_01_llm');
x_bnd_c01_KRG_2_2 = [t_c01_KRG (t2_ls_01_llm-fvalmax_1_c01_KRG)/2]*(W_ls_01_llm');
x_bnd_c01_KRG_2 = [x_bnd_c01_KRG_2_1;x_bnd_c01_KRG_2_2];
x_bnd_c01_KRG_2_denorm = denormalizeuv(x_bnd_c01_KRG_2,lbx,ubx);
x_bnd_c01_KRG_3_1 = (x_bnd_c01_KRG_2_1+xopt_c01_KRG_1)/2; 
x_bnd_c01_KRG_3_2 = (x_bnd_c01_KRG_2_2+xopt_c01_KRG_1)/2; 
x_bnd_c01_KRG_3 = [x_bnd_c01_KRG_3_1;x_bnd_c01_KRG_3_2];
x_bnd_c01_KRG_3_denorm = denormalizeuv(x_bnd_c01_KRG_3,lbx,ubx);
x_bnd_c01_KRG_org = [x_bnd_c01_KRG_1_org;x_bnd_c01_KRG_2_denorm;x_bnd_c01_KRG_3_denorm];

[psf_c01_KRG_2, pf_c01_KRG_2, re_c01_KRG_2, beta_c01_KRG_2] = mcspsfconstraint1c(x_bnd_c01_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);

[mean_c01_KRG_enrico, std_c01_KRG_enrico] = momentsusingderrico(fvalmin_1_c01_KRG,-fvalmin_1_c01_KRG,W_ls_01_llm,sd,t_c01_KRG,lbx,ubx,@mcspsfconstraint1c);
% mean_c01_KRG = 1/6*(psf_c01_KRG_2(1)+psf_c01_KRG_1)+(4/6*(psf_c01_KRG_2(2)));
% std_c01_KRG = (1/6*((psf_c01_KRG_2(1)-mean_c01_KRG).^2))+(1/6*((psf_c01_KRG_1-mean_c01_KRG).^2));
y_c01_KRG = constraint1c(x_bnd_c01_KRG_org);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_KRG = zeros(1,n_ls_02);
func_t_c02_KRG =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_KRG, fval_t_c02_KRG, exitflag_t_c02_KRG] = fmincon(func_t_c02_KRG,sp_t_c02_KRG,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c02,@srgtsKRGEvaluate,1));

%Validation
func_c02_KRG_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_KRG_1, fval_c02_KRG_1, exitflag_c02_KRG_1] = fmincon(func_c02_KRG_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_KRG));
xopt_c02_KRG_1_denorm = denormalizeuv(xopt_c02_KRG_1,lbx_c02,ubx_c02);
xopt_c02_KRG_1_org = xopt_c02_KRG_1_denorm;
[psf_c02_KRG_1,pf_c02_KRG_1,re_c02_KRG_1,beta_c02_KRG_1] = mcspsfconstraint2c(xopt_c02_KRG_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_KRG_1_mm = srgtsKRGEvaluate(t_c02_KRG,srgt_KRG_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);
t2_ls_02_llm = xopt_c02_KRG_1*W2_ls_02_llm;

sp_c02_KRG = zeros(1,m_c02);
fmin_c02_KRG_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_KRG,fvalmin_1_c02_KRG,exitflagmin_1_c02_KRG]=fmincon(fmin_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm));

fmax_c02_KRG_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_KRG,fvalmax_1_c02_KRG,exitflagmax_1_c02_KRG]=fmincon(fmax_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm));

x_bnd_c02_KRG_1 = [xmin_1_c02_KRG;xmax_1_c02_KRG];

ind_c02_KRG = 2*(m_c02-n_ls_02);
x_bnd_c02_KRG_1_org = denormalizeuv(x_bnd_c02_KRG_1,lbx,ubx);
x_bnd_c02_KRG_2_1 = [t_c02_KRG (t2_ls_02_llm+fvalmin_1_c02_KRG)/2]*(W_ls_02_llm');
x_bnd_c02_KRG_2_2 = [t_c02_KRG (t2_ls_02_llm-fvalmax_1_c02_KRG)/2]*(W_ls_02_llm');
x_bnd_c02_KRG_2 = [x_bnd_c02_KRG_2_1;x_bnd_c02_KRG_2_2];
x_bnd_c02_KRG_2_denorm = denormalizeuv(x_bnd_c02_KRG_2,lbx,ubx);
x_bnd_c02_KRG_3_1 = (x_bnd_c02_KRG_2_1+xopt_c02_KRG_1)/2; 
x_bnd_c02_KRG_3_2 = (x_bnd_c02_KRG_2_2+xopt_c02_KRG_1)/2; 
x_bnd_c02_KRG_3 = [x_bnd_c02_KRG_3_1;x_bnd_c02_KRG_3_2];
x_bnd_c02_KRG_3_denorm = denormalizeuv(x_bnd_c02_KRG_3,lbx,ubx);
x_bnd_c02_KRG_org = [x_bnd_c02_KRG_1_org;x_bnd_c02_KRG_2_denorm;x_bnd_c02_KRG_3_denorm];

[psf_c02_KRG_2,pf_c02_KRG_2,re_c02_KRG_2,beta_c02_KRG_2] = mcspsfconstraint2c(x_bnd_c02_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c02_KRG_enrico, std_c02_KRG_enrico] = momentsusingderrico(fvalmin_1_c02_KRG,-fvalmin_1_c02_KRG,W_ls_02_llm,sd,t_c02_KRG,lbx,ubx,@mcspsfconstraint2c);
% mean_c02_KRG = 1/6*(psf_c02_KRG_1+psf_c02_KRG_2(2))+(4/6*(psf_c02_KRG_2(1)));
% std_c02_KRG = (1/6*((psf_c02_KRG_1-mean_c02_KRG).^2))+(1/6*((psf_c02_KRG_2(2)-mean_c02_KRG).^2));
y_c02_KRG = constraint2c(x_bnd_c02_KRG_org);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_KRG = zeros(1,n_ls_03);
func_t_c03_KRG =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_KRG, fval_t_c03_KRG, exitflag_t_c03_KRG] = fmincon(func_t_c03_KRG,sp_t_c03_KRG,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c03,@srgtsKRGEvaluate,1));

%Validation
func_c03_KRG_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_KRG_1, fval_c03_KRG_1, exitflag_c03_KRG_1] = fmincon(func_c03_KRG_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_KRG));
xopt_c03_KRG_1_denorm = denormalizeuv(xopt_c03_KRG_1,lbx_c03,ubx_c03); 
xopt_c03_KRG_1_org = xopt_c03_KRG_1_denorm;
[psf_c03_KRG_1,pf_c03_KRG_1,re_c03_KRG_1,beta_c03_KRG_1] = mcspsfconstraint3c(xopt_c03_KRG_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c03_KRG_1_mm = srgtsKRGEvaluate(t_c03_KRG,srgt_KRG_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);
t2_ls_03_llm = xopt_c03_KRG_1*W2_ls_03_llm;

sp_c03_KRG = zeros(1,m_c03);
fmin_c03_KRG_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_KRG,fvalmin_1_c03_KRG,exitflagmin_1_c03_KRG]=fmincon(fmin_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm));

fmax_c03_KRG_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_KRG,fvalmax_1_c03_KRG,exitflagmax_1_c03_KRG]=fmincon(fmax_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm));

x_bnd_c03_KRG_1 = [xmin_1_c03_KRG;xmax_1_c03_KRG];

ind_c03_KRG = 2*(m_c03-n_ls_03);
x_bnd_c03_KRG_1_org = denormalizeuv(x_bnd_c03_KRG_1,lbx,ubx);
x_bnd_c03_KRG_2_1 = [t_c03_KRG (t2_ls_03_llm+fvalmin_1_c03_KRG)/2]*(W_ls_03_llm');
x_bnd_c03_KRG_2_2 = [t_c03_KRG (t2_ls_03_llm-fvalmax_1_c03_KRG)/2]*(W_ls_03_llm');
x_bnd_c03_KRG_2 = [x_bnd_c03_KRG_2_1;x_bnd_c03_KRG_2_2];
x_bnd_c03_KRG_2_denorm = denormalizeuv(x_bnd_c03_KRG_2,lbx,ubx);
x_bnd_c03_KRG_3_1 = (x_bnd_c03_KRG_2_1+xopt_c03_KRG_1)/2; 
x_bnd_c03_KRG_3_2 = (x_bnd_c03_KRG_2_2+xopt_c03_KRG_1)/2; 
x_bnd_c03_KRG_3 = [x_bnd_c03_KRG_3_1;x_bnd_c03_KRG_3_2];
x_bnd_c03_KRG_3_denorm = denormalizeuv(x_bnd_c03_KRG_3,lbx,ubx);
x_bnd_c03_KRG_org = [x_bnd_c03_KRG_1_org;x_bnd_c03_KRG_2_denorm];

[psf_c03_KRG_2,pf_c03_KRG_2,re_c03_KRG_2,beta_c03_KRG_2] = mcspsfconstraint3c(x_bnd_c03_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c03_KRG_enrico, std_c03_KRG_enrico] = momentsusingderrico(fvalmin_1_c03_KRG,-fvalmin_1_c03_KRG,W_ls_03_llm,sd,t_c03_KRG,lbx,ubx,@mcspsfconstraint3c);
y_c03_KRG = constraint3c(x_bnd_c03_KRG_org);
