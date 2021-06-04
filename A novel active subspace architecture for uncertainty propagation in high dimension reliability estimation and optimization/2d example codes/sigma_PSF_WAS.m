
%% WAS
clc;close all

%% Constraint 1

%Finding active variable at which PSF from metamodel is one:
sp_t_c01_WAS = zeros(1,n_ls_01);
func_t_c01_WAS =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_WAS, fval_t_c01_WAS, exitflag_t_c01_WAS] = fmincon(func_t_c01_WAS,sp_t_c01_WAS,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c01,@srgtsWASEvaluate,mm_psf_WAS_c01));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c01_WAS_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_WAS_1, fval_c01_WAS_1, exitflag_c01_WAS_1] = fmincon(func_c01_WAS_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_WAS));
xopt_c01_WAS_1_denorm = denormalizeuv(xopt_c01_WAS_1,lbx_c01,ubx_c01);
xopt_c01_WAS_1_org = xopt_c01_WAS_1_denorm;
[psf_c01_WAS_1,pf_c01_WAS_1,re_c01_WAS_1,beta_c01_WAS_1] = mcspsfconstraint1c(xopt_c01_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Finding Bounds along inactive subspace of the limit states
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);

sp_c01_WAS = zeros(1,m_c01);
fmin_c01_WAS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_WAS, fvalmin_1_c01_WAS, exitflagmin_1_c01_WAS]=fmincon(fmin_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm));

fmax_c01_WAS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_WAS, fvalmax_1_c01_WAS, exitflagmax_1_c01_WAS]=fmincon(fmax_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states
[mean_c01_WAS_enrico, std_c01_WAS_enrico] = momentsusingderrico(fvalmin_1_c01_WAS,-fvalmin_1_c01_WAS,W_ls_01_llm,sd,t_c01_WAS,lbx,ubx,@mcspsfconstraint1c);

%% Constraint 2

%Finding active variable at which PSF from metamodel is one:
sp_t_c02_WAS = zeros(1,n_ls_02);
func_t_c02_WAS =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_WAS, fval_t_c02_WAS, exitflag_t_c02_WAS] = fmincon(func_t_c02_WAS,sp_t_c02_WAS,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c02,@srgtsWASEvaluate,mm_psf_WAS_c02));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c02_WAS_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_WAS_1, fval_c02_WAS_1, exitflag_c02_WAS_1] = fmincon(func_c02_WAS_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_WAS));
xopt_c02_WAS_1_denorm = denormalizeuv(xopt_c02_WAS_1,lbx_c02,ubx_c02);
xopt_c02_WAS_1_org = xopt_c02_WAS_1_denorm;
[psf_c02_WAS_1,pf_c02_WAS_1,re_c02_WAS_1,beta_c02_WAS_1] = mcspsfconstraint2c(xopt_c02_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_WAS_1_mm = srgtsWASEvaluate(t_c02_WAS,srgt_WAS_c02);

%Finding Bounds along inactive subspace of the limit states:
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);

sp_c02_WAS = zeros(1,m_c02);
fmin_c02_WAS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_WAS,fvalmin_1_c02_WAS,exitflagmin_1_c02_WAS]=fmincon(fmin_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm));

fmax_c02_WAS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_WAS,fvalmax_1_c02_WAS,exitflagmax_1_c02_WAS]=fmincon(fmax_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c02_WAS_enrico, std_c02_WAS_enrico] = momentsusingderrico(fvalmin_1_c02_WAS,-fvalmin_1_c02_WAS,W_ls_02_llm,sd,t_c02_WAS,lbx,ubx,@mcspsfconstraint2c);

%% Constraint 3

%Finding active variable at which PSF from metamodel is one:
sp_t_c03_WAS = zeros(1,n_ls_03);
func_t_c03_WAS =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_WAS, fval_t_c03_WAS, exitflag_t_c03_WAS] = fmincon(func_t_c03_WAS,sp_t_c03_WAS,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c03,@srgtsWASEvaluate,mm_psf_WAS_c03));

%Finding active variable at which PSF from metamodel is one:
func_c03_WAS_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_WAS_1, fval_c03_WAS_1, exitflag_c03_WAS_1] = fmincon(func_c03_WAS_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_WAS));
xopt_c03_WAS_1_denorm = denormalizeuv(xopt_c03_WAS_1,lbx_c03,ubx_c03); 
xopt_c03_WAS_1_org = xopt_c03_WAS_1_denorm;
[psf_c03_WAS_1,pf_c03_WAS_1,re_c03_WAS_1,beta_c03_WAS_1] = mcspsfconstraint3c(xopt_c03_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);

sp_c03_WAS = zeros(1,m_c03);
fmin_c03_WAS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_WAS,fvalmin_1_c03_WAS,exitflagmin_1_c03_WAS]=fmincon(fmin_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm));

fmax_c03_WAS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_WAS,fvalmax_1_c03_WAS,exitflagmax_1_c03_WAS]=fmincon(fmax_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c03_WAS_enrico, std_c03_WAS_enrico] = momentsusingderrico(fvalmin_1_c03_WAS,-fvalmin_1_c03_WAS,W_ls_03_llm,sd,t_c03_WAS,lbx,ubx,@mcspsfconstraint3c);
