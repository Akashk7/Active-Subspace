 clc;close all

%% Constraint 1

%Finding active variable at which PSF from metamodel is one:
sp_t_c01_RBF = zeros(1,n_ls_01);
func_t_c01_RBF =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_RBF, fval_t_c01_RBF, exitflag_t_c01_RBF] = fmincon(func_t_c01_RBF,sp_t_c01_RBF,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c01,@srgtsRBFEvaluate,1));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c01_RBF_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_RBF_1, fval_c01_RBF_1, exitflag_c01_RBF_1] = fmincon(func_c01_RBF_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_RBF));
xopt_c01_RBF_1_denorm = denormalizeuv(xopt_c01_RBF_1,lbx_c01,ubx_c01);
xopt_c01_RBF_1_org = xopt_c01_RBF_1_denorm;
[psf_c01_RBF_1,pf_c01_RBF_1,re_c01_RBF_1,beta_c01_RBF_1] = mcspsfconstraint1c(xopt_c01_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Finding Bounds along inactive subspace of the limit states
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);

sp_c01_RBF = zeros(1,m_c01);
fmin_c01_RBF_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_RBF, fvalmin_1_c01_RBF, exitflagmin_1_c01_RBF]=fmincon(fmin_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

fmax_c01_RBF_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_RBF, fvalmax_1_c01_RBF, exitflagmax_1_c01_RBF]=fmincon(fmax_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states
[mean_c01_RBF_enrico, std_c01_RBF_enrico] = momentsusingderrico(fvalmin_1_c01_RBF,-fvalmin_1_c01_RBF,W_ls_01_llm,sd,t_c01_RBF,lbx,ubx,@mcspsfconstraint1c);

%% Constraint 2

%Finding active variable at which PSF from metamodel is one:
sp_t_c02_RBF = zeros(1,n_ls_02);
func_t_c02_RBF =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_RBF, fval_t_c02_RBF, exitflag_t_c02_RBF] = fmincon(func_t_c02_RBF,sp_t_c02_RBF,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c02,@srgtsRBFEvaluate,1));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c02_RBF_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_RBF_1, fval_c02_RBF_1, exitflag_c02_RBF_1] = fmincon(func_c02_RBF_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_RBF));
xopt_c02_RBF_1_denorm = denormalizeuv(xopt_c02_RBF_1,lbx_c02,ubx_c02);
xopt_c02_RBF_1_org = xopt_c02_RBF_1_denorm;
[psf_c02_RBF_1,pf_c02_RBF_1,re_c02_RBF_1,beta_c02_RBF_1] = mcspsfconstraint2c(xopt_c02_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_RBF_1_mm = srgtsRBFEvaluate(t_c02_RBF,srgt_RBF_c02);

%Finding Bounds along inactive subspace of the limit states:
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);

sp_c02_RBF = zeros(1,m_c02);
fmin_c02_RBF_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_RBF,fvalmin_1_c02_RBF,exitflagmin_1_c02_RBF]=fmincon(fmin_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

fmax_c02_RBF_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_RBF,fvalmax_1_c02_RBF,exitflagmax_1_c02_RBF]=fmincon(fmax_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c02_RBF_enrico, std_c02_RBF_enrico] = momentsusingderrico(fvalmin_1_c02_RBF,-fvalmin_1_c02_RBF,W_ls_02_llm,sd,t_c02_RBF,lbx,ubx,@mcspsfconstraint2c);

%% Constraint 3

%Finding active variable at which PSF from metamodel is one:
sp_t_c03_RBF = zeros(1,n_ls_03);
func_t_c03_RBF =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_RBF, fval_t_c03_RBF, exitflag_t_c03_RBF] = fmincon(func_t_c03_RBF,sp_t_c03_RBF,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c03,@srgtsRBFEvaluate,1));

%Finding active variable at which PSF from metamodel is one:
func_c03_RBF_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_RBF_1, fval_c03_RBF_1, exitflag_c03_RBF_1] = fmincon(func_c03_RBF_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_RBF));
xopt_c03_RBF_1_denorm = denormalizeuv(xopt_c03_RBF_1,lbx_c03,ubx_c03); 
xopt_c03_RBF_1_org = xopt_c03_RBF_1_denorm;
[psf_c03_RBF_1,pf_c03_RBF_1,re_c03_RBF_1,beta_c03_RBF_1] = mcspsfconstraint3c(xopt_c03_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);

sp_c03_RBF = zeros(1,m_c03);
fmin_c03_RBF_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_RBF,fvalmin_1_c03_RBF,exitflagmin_1_c03_RBF]=fmincon(fmin_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

fmax_c03_RBF_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_RBF,fvalmax_1_c03_RBF,exitflagmax_1_c03_RBF]=fmincon(fmax_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c03_RBF_enrico, std_c03_RBF_enrico] = momentsusingderrico(fvalmin_1_c03_RBF,-fvalmin_1_c03_RBF,W_ls_03_llm,sd,t_c03_RBF,lbx,ubx,@mcspsfconstraint3c);
 clc;close all

%% Constraint 1

%Finding active variable at which PSF from metamodel is one:
sp_t_c01_RBF = zeros(1,n_ls_01);
func_t_c01_RBF =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_RBF, fval_t_c01_RBF, exitflag_t_c01_RBF] = fmincon(func_t_c01_RBF,sp_t_c01_RBF,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c01,@srgtsRBFEvaluate,1));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c01_RBF_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_RBF_1, fval_c01_RBF_1, exitflag_c01_RBF_1] = fmincon(func_c01_RBF_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_RBF));
xopt_c01_RBF_1_denorm = denormalizeuv(xopt_c01_RBF_1,lbx_c01,ubx_c01);
xopt_c01_RBF_1_org = xopt_c01_RBF_1_denorm;
[psf_c01_RBF_1,pf_c01_RBF_1,re_c01_RBF_1,beta_c01_RBF_1] = mcspsfconstraint1c(xopt_c01_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Finding Bounds along inactive subspace of the limit states
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);

sp_c01_RBF = zeros(1,m_c01);
fmin_c01_RBF_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_RBF, fvalmin_1_c01_RBF, exitflagmin_1_c01_RBF]=fmincon(fmin_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

fmax_c01_RBF_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_RBF, fvalmax_1_c01_RBF, exitflagmax_1_c01_RBF]=fmincon(fmax_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states
[mean_c01_RBF_enrico, std_c01_RBF_enrico] = momentsusingderrico(fvalmin_1_c01_RBF,-fvalmin_1_c01_RBF,W_ls_01_llm,sd,t_c01_RBF,lbx,ubx,@mcspsfconstraint1c);

%% Constraint 2

%Finding active variable at which PSF from metamodel is one:
sp_t_c02_RBF = zeros(1,n_ls_02);
func_t_c02_RBF =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_RBF, fval_t_c02_RBF, exitflag_t_c02_RBF] = fmincon(func_t_c02_RBF,sp_t_c02_RBF,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c02,@srgtsRBFEvaluate,1));

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
func_c02_RBF_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_RBF_1, fval_c02_RBF_1, exitflag_c02_RBF_1] = fmincon(func_c02_RBF_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_RBF));
xopt_c02_RBF_1_denorm = denormalizeuv(xopt_c02_RBF_1,lbx_c02,ubx_c02);
xopt_c02_RBF_1_org = xopt_c02_RBF_1_denorm;
[psf_c02_RBF_1,pf_c02_RBF_1,re_c02_RBF_1,beta_c02_RBF_1] = mcspsfconstraint2c(xopt_c02_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_RBF_1_mm = srgtsRBFEvaluate(t_c02_RBF,srgt_RBF_c02);

%Finding Bounds along inactive subspace of the limit states:
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);

sp_c02_RBF = zeros(1,m_c02);
fmin_c02_RBF_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_RBF,fvalmin_1_c02_RBF,exitflagmin_1_c02_RBF]=fmincon(fmin_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

fmax_c02_RBF_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_RBF,fvalmax_1_c02_RBF,exitflagmax_1_c02_RBF]=fmincon(fmax_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c02_RBF_enrico, std_c02_RBF_enrico] = momentsusingderrico(fvalmin_1_c02_RBF,-fvalmin_1_c02_RBF,W_ls_02_llm,sd,t_c02_RBF,lbx,ubx,@mcspsfconstraint2c);

%% Constraint 3

%Finding active variable at which PSF from metamodel is one:
sp_t_c03_RBF = zeros(1,n_ls_03);
func_t_c03_RBF =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_RBF, fval_t_c03_RBF, exitflag_t_c03_RBF] = fmincon(func_t_c03_RBF,sp_t_c03_RBF,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c03,@srgtsRBFEvaluate,1));

%Finding active variable at which PSF from metamodel is one:
func_c03_RBF_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_RBF_1, fval_c03_RBF_1, exitflag_c03_RBF_1] = fmincon(func_c03_RBF_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_RBF));
xopt_c03_RBF_1_denorm = denormalizeuv(xopt_c03_RBF_1,lbx_c03,ubx_c03); 
xopt_c03_RBF_1_org = xopt_c03_RBF_1_denorm;
[psf_c03_RBF_1,pf_c03_RBF_1,re_c03_RBF_1,beta_c03_RBF_1] = mcspsfconstraint3c(xopt_c03_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);

%Evaluating actual PSF value at active variable at which PSF from metamodel is one: 
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);

sp_c03_RBF = zeros(1,m_c03);
fmin_c03_RBF_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_RBF,fvalmin_1_c03_RBF,exitflagmin_1_c03_RBF]=fmincon(fmin_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

fmax_c03_RBF_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_RBF,fvalmax_1_c03_RBF,exitflagmax_1_c03_RBF]=fmincon(fmax_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

%Computing standard deviation of PSF along inactive subspace of the limit states:
[mean_c03_RBF_enrico, std_c03_RBF_enrico] = momentsusingderrico(fvalmin_1_c03_RBF,-fvalmin_1_c03_RBF,W_ls_03_llm,sd,t_c03_RBF,lbx,ubx,@mcspsfconstraint3c);
