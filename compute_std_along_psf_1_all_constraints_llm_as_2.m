 clc;close all

%% Constraint 1

%Finding t at which PSF is one
sp_t_c01_RBF = zeros(1,n_ls_01);
func_t_c01_RBF =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_RBF, fval_t_c01_RBF, exitflag_t_c01_RBF] = fmincon(func_t_c01_RBF,sp_t_c01_RBF,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c01,@srgtsRBFEvaluate,1));

%Validation
func_c01_RBF_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_RBF_1, fval_c01_RBF_1, exitflag_c01_RBF_1] = fmincon(func_c01_RBF_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_RBF));
xopt_c01_RBF_1_denorm = denormalizeuv(xopt_c01_RBF_1,lbx_c01,ubx_c01);
xopt_c01_RBF_1_org = xopt_c01_RBF_1_denorm;
[psf_c01_RBF_1,pf_c01_RBF_1,re_c01_RBF_1,beta_c01_RBF_1] = mcspsfconstraint1c(xopt_c01_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c01_RBF_1_mm = srgtsRBFEvaluate(t_c01_RBF,srgt_RBF_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);
t2_ls_01_llm = xopt_c01_RBF_1*W2_ls_01_llm;

sp_c01_RBF = zeros(1,m_c01);
fmin_c01_RBF_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_RBF, fvalmin_1_c01_RBF, exitflagmin_1_c01_RBF]=fmincon(fmin_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

fmax_c01_RBF_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_RBF, fvalmax_1_c01_RBF, exitflagmax_1_c01_RBF]=fmincon(fmax_c01_RBF_1,sp_c01_RBF,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_RBF,W1_ls_01_llm));

x_bnd_c01_RBF_1 = [xmin_1_c01_RBF;xmax_1_c01_RBF];

ind_c01_RBF = 2*(m_c01-n_ls_01);
x_bnd_c01_RBF_1_org = denormalizeuv(x_bnd_c01_RBF_1,lbx,ubx);
x_bnd_c01_RBF_2_1 = [t_c01_RBF (t2_ls_01_llm+fvalmin_1_c01_RBF)/2]*(W_ls_01_llm');
x_bnd_c01_RBF_2_2 = [t_c01_RBF (t2_ls_01_llm-fvalmax_1_c01_RBF)/2]*(W_ls_01_llm');
x_bnd_c01_RBF_2 = [x_bnd_c01_RBF_2_1;x_bnd_c01_RBF_2_2];
x_bnd_c01_RBF_2_denorm = denormalizeuv(x_bnd_c01_RBF_2,lbx,ubx);
x_bnd_c01_RBF_3_1 = (x_bnd_c01_RBF_2_1+xopt_c01_RBF_1)/2; 
x_bnd_c01_RBF_3_2 = (x_bnd_c01_RBF_2_2+xopt_c01_RBF_1)/2; 
x_bnd_c01_RBF_3 = [x_bnd_c01_RBF_3_1;x_bnd_c01_RBF_3_2];
x_bnd_c01_RBF_3_denorm = denormalizeuv(x_bnd_c01_RBF_3,lbx,ubx);
x_bnd_c01_RBF_org = [x_bnd_c01_RBF_1_org;x_bnd_c01_RBF_2_denorm;x_bnd_c01_RBF_3_denorm];

[psf_c01_RBF_2, pf_c01_RBF_2, re_c01_RBF_2, beta_c01_RBF_2] = mcspsfconstraint1c(x_bnd_c01_RBF_org,lbx,ubx,sd,Nmcs,Pftarget);

[mean_c01_RBF, std_c01_RBF] = momentsusingyounetal([psf_c01_RBF_2;psf_c01_RBF_1]);
[mean_c01_RBF_enrico, std_c01_RBF_enrico] = momentsusingderrico(fvalmin_1_c01_RBF,-fvalmin_1_c01_RBF,W_ls_01_llm,sd,t_c01_RBF,lbx,ubx,@mcspsfconstraint1c);
% mean_c01_RBF = 1/6*(psf_c01_RBF_2(1)+psf_c01_RBF_1)+(4/6*(psf_c01_RBF_2(2)));
% std_c01_RBF = (1/6*((psf_c01_RBF_2(1)-mean_c01_RBF).^2))+(1/6*((psf_c01_RBF_1-mean_c01_RBF).^2));
y_c01_RBF = constraint1c(x_bnd_c01_RBF_org);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_RBF = zeros(1,n_ls_02);
func_t_c02_RBF =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_RBF, fval_t_c02_RBF, exitflag_t_c02_RBF] = fmincon(func_t_c02_RBF,sp_t_c02_RBF,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c02,@srgtsRBFEvaluate,1));

%Validation
func_c02_RBF_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_RBF_1, fval_c02_RBF_1, exitflag_c02_RBF_1] = fmincon(func_c02_RBF_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_RBF));
xopt_c02_RBF_1_denorm = denormalizeuv(xopt_c02_RBF_1,lbx_c02,ubx_c02);
xopt_c02_RBF_1_org = xopt_c02_RBF_1_denorm;
[psf_c02_RBF_1,pf_c02_RBF_1,re_c02_RBF_1,beta_c02_RBF_1] = mcspsfconstraint2c(xopt_c02_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_RBF_1_mm = srgtsRBFEvaluate(t_c02_RBF,srgt_RBF_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);
t2_ls_02_llm = xopt_c02_RBF_1*W2_ls_02_llm;

sp_c02_RBF = zeros(1,m_c02);
fmin_c02_RBF_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_RBF,fvalmin_1_c02_RBF,exitflagmin_1_c02_RBF]=fmincon(fmin_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

fmax_c02_RBF_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_RBF,fvalmax_1_c02_RBF,exitflagmax_1_c02_RBF]=fmincon(fmax_c02_RBF_1,sp_c02_RBF,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_RBF,W1_ls_02_llm));

x_bnd_c02_RBF_1 = [xmin_1_c02_RBF;xmax_1_c02_RBF];

ind_c02_RBF = 2*(m_c02-n_ls_02);
x_bnd_c02_RBF_1_org = denormalizeuv(x_bnd_c02_RBF_1,lbx,ubx);
x_bnd_c02_RBF_2_1 = [t_c02_RBF (t2_ls_02_llm+fvalmin_1_c02_RBF)/2]*(W_ls_02_llm');
x_bnd_c02_RBF_2_2 = [t_c02_RBF (t2_ls_02_llm-fvalmax_1_c02_RBF)/2]*(W_ls_02_llm');
x_bnd_c02_RBF_2 = [x_bnd_c02_RBF_2_1;x_bnd_c02_RBF_2_2];
x_bnd_c02_RBF_2_denorm = denormalizeuv(x_bnd_c02_RBF_2,lbx,ubx);
x_bnd_c02_RBF_3_1 = (x_bnd_c02_RBF_2_1+xopt_c02_RBF_1)/2; 
x_bnd_c02_RBF_3_2 = (x_bnd_c02_RBF_2_2+xopt_c02_RBF_1)/2; 
x_bnd_c02_RBF_3 = [x_bnd_c02_RBF_3_1;x_bnd_c02_RBF_3_2];
x_bnd_c02_RBF_3_denorm = denormalizeuv(x_bnd_c02_RBF_3,lbx,ubx);
x_bnd_c02_RBF_org = [x_bnd_c02_RBF_1_org;x_bnd_c02_RBF_2_denorm;x_bnd_c02_RBF_3_denorm];

[psf_c02_RBF_2,pf_c02_RBF_2,re_c02_RBF_2,beta_c02_RBF_2] = mcspsfconstraint2c(x_bnd_c02_RBF_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c02_RBF, std_c02_RBF] = momentsusingyounetal([psf_c02_RBF_2;psf_c02_RBF_1]);
[mean_c02_RBF_enrico, std_c02_RBF_enrico] = momentsusingderrico(fvalmin_1_c02_RBF,-fvalmin_1_c02_RBF,W_ls_02_llm,sd,t_c02_RBF,lbx,ubx,@mcspsfconstraint2c);
% mean_c02_RBF = 1/6*(psf_c02_RBF_1+psf_c02_RBF_2(2))+(4/6*(psf_c02_RBF_2(1)));
% std_c02_RBF = (1/6*((psf_c02_RBF_1-mean_c02_RBF).^2))+(1/6*((psf_c02_RBF_2(2)-mean_c02_RBF).^2));
y_c02_RBF = constraint2c(x_bnd_c02_RBF_org);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_RBF = zeros(1,n_ls_03);
func_t_c03_RBF =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_RBF, fval_t_c03_RBF, exitflag_t_c03_RBF] = fmincon(func_t_c03_RBF,sp_t_c03_RBF,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_c03,@srgtsRBFEvaluate,1));

%Validation
func_c03_RBF_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_RBF_1, fval_c03_RBF_1, exitflag_c03_RBF_1] = fmincon(func_c03_RBF_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_RBF));
xopt_c03_RBF_1_denorm = denormalizeuv(xopt_c03_RBF_1,lbx_c03,ubx_c03); 
xopt_c03_RBF_1_org = xopt_c03_RBF_1_denorm;
[psf_c03_RBF_1,pf_c03_RBF_1,re_c03_RBF_1,beta_c03_RBF_1] = mcspsfconstraint3c(xopt_c03_RBF_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c03_RBF_1_mm = srgtsRBFEvaluate(t_c03_RBF,srgt_RBF_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);
t2_ls_03_llm = xopt_c03_RBF_1*W2_ls_03_llm;

sp_c03_RBF = zeros(1,m_c03);
fmin_c03_RBF_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_RBF,fvalmin_1_c03_RBF,exitflagmin_1_c03_RBF]=fmincon(fmin_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

fmax_c03_RBF_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_RBF,fvalmax_1_c03_RBF,exitflagmax_1_c03_RBF]=fmincon(fmax_c03_RBF_1,sp_c03_RBF,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_RBF,W1_ls_03_llm));

x_bnd_c03_RBF_1 = [xmin_1_c03_RBF;xmax_1_c03_RBF];

ind_c03_RBF = 2*(m_c03-n_ls_03);
x_bnd_c03_RBF_1_org = denormalizeuv(x_bnd_c03_RBF_1,lbx,ubx);
x_bnd_c03_RBF_2_1 = [t_c03_RBF (t2_ls_03_llm+fvalmin_1_c03_RBF)/2]*(W_ls_03_llm');
x_bnd_c03_RBF_2_2 = [t_c03_RBF (t2_ls_03_llm-fvalmax_1_c03_RBF)/2]*(W_ls_03_llm');
x_bnd_c03_RBF_2 = [x_bnd_c03_RBF_2_1;x_bnd_c03_RBF_2_2];
x_bnd_c03_RBF_2_denorm = denormalizeuv(x_bnd_c03_RBF_2,lbx,ubx);
x_bnd_c03_RBF_3_1 = (x_bnd_c03_RBF_2_1+xopt_c03_RBF_1)/2; 
x_bnd_c03_RBF_3_2 = (x_bnd_c03_RBF_2_2+xopt_c03_RBF_1)/2; 
x_bnd_c03_RBF_3 = [x_bnd_c03_RBF_3_1;x_bnd_c03_RBF_3_2];
x_bnd_c03_RBF_3_denorm = denormalizeuv(x_bnd_c03_RBF_3,lbx,ubx);
x_bnd_c03_RBF_org = [x_bnd_c03_RBF_1_org;x_bnd_c03_RBF_2_denorm];

[psf_c03_RBF_2,pf_c03_RBF_2,re_c03_RBF_2,beta_c03_RBF_2] = mcspsfconstraint3c(x_bnd_c03_RBF_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c03_RBF, std_c03_RBF] = momentsusingyounetal([psf_c03_RBF_2;psf_c03_RBF_1]);
[mean_c03_RBF_enrico, std_c03_RBF_enrico] = momentsusingderrico(fvalmin_1_c03_RBF,-fvalmin_1_c03_RBF,W_ls_03_llm,sd,t_c03_RBF,lbx,ubx,@mcspsfconstraint3c);

% mean_c03_RBF = 1/6*(psf_c03_RBF_2(1)+psf_c03_RBF_1)+(4/6*(psf_c03_RBF_2(2)));
% std_c03_RBF = (1/6*((psf_c03_RBF_2(1)-mean_c03_RBF).^2))+(1/6*((psf_c03_RBF_1-mean_c03_RBF).^2));
% std_c03_RBF = std([psf_c03_RBF_2;psf_c03_RBF_1]);

y_c03_RBF = constraint3c(x_bnd_c03_RBF_org);

%% Weight

%Finding t at which PSF is one
sp_t_obj_RBF = zeros(1,n_obj);
func_t_obj_RBF =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_RBF, fval_t_obj_RBF, exitflag_t_obj_RBF] = fmincon(func_t_obj_RBF,sp_t_obj_RBF,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_RBF_obj_1,@srgtsRBFEvaluate,mean(y_obj_new)));

%Validation
func_obj_RBF_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_RBF_1, fval_obj_RBF_1, exitflag_obj_RBF_1] = fmincon(func_obj_RBF_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_RBF));
xopt_obj_RBF_1_denorm = denormalizeuv(xopt_obj_RBF_1,lbx_obj,ubx_obj); 
xopt_obj_RBF_1_org = xopt_obj_RBF_1_denorm;
obj_RBF_1 = Weightc(xopt_obj_RBF_1_org);
obj_RBF_1_mm = srgtsRBFEvaluate(t_obj_RBF,srgt_RBF_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);

sp_obj_RBF = zeros(1,m_obj);
fmin_obj_RBF_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_RBF,fvalmin_1_obj_RBF,exitflagmin_1_obj_RBF]=fmincon(fmin_obj_RBF_1,sp_obj_RBF,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_RBF,W1_obj_llm));     

fmax_obj_RBF_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_RBF,fvalmax_1_obj_RBF,exitflagmax_1_obj_RBF]=fmincon(fmax_obj_RBF_1,sp_obj_RBF,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_RBF,W1_obj_llm));

x_bnd_obj_RBF_1 = [xmin_1_obj_RBF;xmax_1_obj_RBF];

ind_obj_RBF = 2*(m_obj-n_obj);
x_bnd_obj_RBF_org = denormalizeuv(x_bnd_obj_RBF_1,lbx_obj,ubx_obj);
obj_RBF_2 = Weightc(x_bnd_obj_RBF_org);

std_obj_RBF = std([obj_RBF_2;obj_RBF_1]);
y_obj_RBF = Weightc(x_bnd_obj_RBF_org);

%% KRG
clc;close all

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

sp_c01_KRG = zeros(1,m_c01);
fmin_c01_KRG_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_KRG,fvalmin_1_c01_KRG,exitflagmin_1_c01_KRG]=fmincon(fmin_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm));

fmax_c01_KRG_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_KRG,fvalmax_1_c01_KRG,exitflagmax_1_c01_KRG]=fmincon(fmax_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm));

x_bnd_c01_KRG_1 = [xmin_1_c01_KRG;xmax_1_c01_KRG];

ind_c01_KRG = 2*(m_c01-n_ls_01);
x_bnd_c01_KRG_org = denormalizeuv(x_bnd_c01_KRG_1,lbx,ubx);
[psf_c01_KRG_2,pf_c01_KRG_2,re_c01_KRG_2,beta_c01_KRG_2] = mcspsfconstraint1c(x_bnd_c01_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c01_KRG = mean([psf_c01_KRG_2;psf_c01_KRG_1]);
std_c01_KRG = std([psf_c01_KRG_2;psf_c01_KRG_1]);
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

sp_c02_KRG = zeros(1,m_c02);
fmin_c02_KRG_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_KRG,fvalmin_1_c02_KRG,exitflagmin_1_c02_KRG]=fmincon(fmin_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm));

fmax_c02_KRG_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_KRG,fvalmax_1_c02_KRG,exitflagmax_1_c02_KRG]=fmincon(fmax_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm));

x_bnd_c02_KRG_1 = [xmin_1_c02_KRG;xmax_1_c02_KRG];

ind_c02_KRG = 2*(m_c02-n_ls_02);
x_bnd_c02_KRG_org = denormalizeuv(x_bnd_c02_KRG_1,lbx,ubx);
[psf_c02_KRG_2,pf_c02_KRG_2,re_c02_KRG_2,beta_c02_KRG_2] = mcspsfconstraint2c(x_bnd_c02_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c02_KRG = mean([psf_c02_KRG_2;psf_c02_KRG_1]);
std_c02_KRG = std([psf_c02_KRG_2;psf_c02_KRG_1]);
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


sp_c03_KRG = zeros(1,m_c03);
fmin_c03_KRG_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_KRG,fvalmin_1_c03_KRG,exitflagmin_1_c03_KRG]=fmincon(fmin_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm));

fmax_c03_KRG_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_KRG,fvalmax_1_c03_KRG,exitflagmax_1_c03_KRG]=fmincon(fmax_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm));

x_bnd_c03_KRG_1 = [xmin_1_c03_KRG;xmax_1_c03_KRG];

ind_c03_KRG = 2*(m_c03-n_ls_03);
x_bnd_c03_KRG_org = denormalizeuv(x_bnd_c03_KRG_1,lbx,ubx);
[psf_c03_KRG_2,pf_c03_KRG_2,re_c03_KRG_2,beta_c03_KRG_2] = mcspsfconstraint3c(x_bnd_c03_KRG_org,lbx,ubx,sd,Nmcs,Pftarget);

std_c03_KRG = std([psf_c03_KRG_2;psf_c03_KRG_1]);
y_c03_KRG = constraint3c(x_bnd_c03_KRG_org);

%% Weight

%Finding t at which PSF is one
sp_t_obj_KRG = zeros(1,n_obj);
func_t_obj_KRG =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_KRG, fval_t_obj_KRG, exitflag_t_obj_KRG] = fmincon(func_t_obj_KRG,sp_t_obj_KRG,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_obj_1,@srgtsKRGEvaluate,mean(y_obj_new)));

%Validation
func_obj_KRG_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_KRG_1, fval_obj_KRG_1, exitflag_obj_KRG_1] = fmincon(func_obj_KRG_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_KRG));
xopt_obj_KRG_1_denorm = denormalizeuv(xopt_obj_KRG_1,lbx_obj,ubx_obj); 
xopt_obj_KRG_1_org = xopt_obj_KRG_1_denorm;
obj_KRG_1 = Weightc(xopt_obj_KRG_1_org);
obj_KRG_1_mm = srgtsKRGEvaluate(t_obj_KRG,srgt_KRG_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);

sp_obj_KRG = zeros(1,m_obj);
fmin_obj_KRG_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_KRG,fvalmin_1_obj_KRG,exitflagmin_1_obj_KRG]=fmincon(fmin_obj_KRG_1,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm));     

fmax_obj_KRG_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_KRG,fvalmax_1_obj_KRG,exitflagmax_1_obj_KRG]=fmincon(fmax_obj_KRG_1,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm));

x_bnd_obj_KRG_1 = [xmin_1_obj_KRG;xmax_1_obj_KRG];

ind_obj_KRG = 2*(m_obj-n_obj);
x_bnd_obj_KRG_org = denormalizeuv(x_bnd_obj_KRG_1,lbx_obj,ubx_obj);
obj_KRG_2 = Weightc(x_bnd_obj_KRG_org);

std_obj_KRG = std([obj_KRG_2;obj_KRG_1]);
y_obj_KRG = Weightc(x_bnd_obj_KRG_org);

%% PRS
clc;close all

%% Constraint 1

%Finding t at which PSF is one
sp_t_c01_PRS = zeros(1,n_ls_01);
func_t_c01_PRS =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_PRS, fval_t_c01_PRS, exitflag_t_c01_PRS] = fmincon(func_t_c01_PRS,sp_t_c01_PRS,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c01,@srgtsPRSEvaluate,1));

%Validation
func_c01_PRS_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_PRS_1, fval_c01_PRS_1, exitflag_c01_PRS_1] = fmincon(func_c01_PRS_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_PRS));
xopt_c01_PRS_1_denorm = denormalizeuv(xopt_c01_PRS_1,lbx_c01,ubx_c01);
xopt_c01_PRS_1_org = xopt_c01_PRS_1_denorm;
[psf_c01_PRS_1,pf_c01_PRS_1,re_c01_PRS_1,beta_c01_PRS_1] = mcspsfconstraint1c(xopt_c01_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c01_PRS_1_mm = srgtsPRSEvaluate(t_c01_PRS,srgt_PRS_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);

sp_c01_PRS = zeros(1,m_c01);
fmin_c01_PRS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_PRS,fvalmin_1_c01_PRS,exitflagmin_1_c01_PRS]=fmincon(fmin_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_PRS,W1_ls_01_llm));

fmax_c01_PRS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_PRS,fvalmax_1_c01_PRS,exitflagmax_1_c01_PRS]=fmincon(fmax_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_PRS,W1_ls_01_llm));

x_bnd_c01_PRS_1 = [xmin_1_c01_PRS;xmax_1_c01_PRS];

ind_c01_PRS = 2*(m_c01-n_ls_01);
x_bnd_c01_PRS_org = denormalizeuv(x_bnd_c01_PRS_1,lbx,ubx);
[psf_c01_PRS_2,pf_c01_PRS_2,re_c01_PRS_2,beta_c01_PRS_2] = mcspsfconstraint1c(x_bnd_c01_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c01_PRS = mean([psf_c01_PRS_2;psf_c01_PRS_1]);
std_c01_PRS = std([psf_c01_PRS_2;psf_c01_PRS_1]);
y_c01_PRS = constraint1c(x_bnd_c01_PRS_org);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_PRS = zeros(1,n_ls_02);
func_t_c02_PRS =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_PRS, fval_t_c02_PRS, exitflag_t_c02_PRS] = fmincon(func_t_c02_PRS,sp_t_c02_PRS,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c02,@srgtsPRSEvaluate,1));

%Validation
func_c02_PRS_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_PRS_1, fval_c02_PRS_1, exitflag_c02_PRS_1] = fmincon(func_c02_PRS_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_PRS));
xopt_c02_PRS_1_denorm = denormalizeuv(xopt_c02_PRS_1,lbx_c02,ubx_c02);
xopt_c02_PRS_1_org = xopt_c02_PRS_1_denorm;
[psf_c02_PRS_1,pf_c02_PRS_1,re_c02_PRS_1,beta_c02_PRS_1] = mcspsfconstraint2c(xopt_c02_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_PRS_1_mm = srgtsPRSEvaluate(t_c02_PRS,srgt_PRS_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);

sp_c02_PRS = zeros(1,m_c02);
fmin_c02_PRS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_PRS,fvalmin_1_c02_PRS,exitflagmin_1_c02_PRS]=fmincon(fmin_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_PRS,W1_ls_02_llm));

fmax_c02_PRS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_PRS,fvalmax_1_c02_PRS,exitflagmax_1_c02_PRS]=fmincon(fmax_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_PRS,W1_ls_02_llm));

x_bnd_c02_PRS_1 = [xmin_1_c02_PRS;xmax_1_c02_PRS];

ind_c02_PRS = 2*(m_c02-n_ls_02);
x_bnd_c02_PRS_org = denormalizeuv(x_bnd_c02_PRS_1,lbx,ubx);
[psf_c02_PRS_2,pf_c02_PRS_2,re_c02_PRS_2,beta_c02_PRS_2] = mcspsfconstraint2c(x_bnd_c02_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c02_PRS = mean([psf_c02_PRS_2;psf_c02_PRS_1]);
std_c02_PRS = std([psf_c02_PRS_2;psf_c02_PRS_1]);
y_c02_PRS = constraint2c(x_bnd_c02_PRS_org);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_PRS = zeros(1,n_ls_03);
func_t_c03_PRS =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_PRS, fval_t_c03_PRS, exitflag_t_c03_PRS] = fmincon(func_t_c03_PRS,sp_t_c03_PRS,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c03,@srgtsPRSEvaluate,1));

%Validation
func_c03_PRS_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_PRS_1, fval_c03_PRS_1, exitflag_c03_PRS_1] = fmincon(func_c03_PRS_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_PRS));
xopt_c03_PRS_1_denorm = denormalizeuv(xopt_c03_PRS_1,lbx_c03,ubx_c03); 
xopt_c03_PRS_1_org = xopt_c03_PRS_1_denorm;
[psf_c03_PRS_1,pf_c03_PRS_1,re_c03_PRS_1,beta_c03_PRS_1] = mcspsfconstraint3c(xopt_c03_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c03_PRS_1_mm = srgtsPRSEvaluate(t_c03_PRS,srgt_PRS_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);


sp_c03_PRS = zeros(1,m_c03);
fmin_c03_PRS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_PRS,fvalmin_1_c03_PRS,exitflagmin_1_c03_PRS]=fmincon(fmin_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_PRS,W1_ls_03_llm));

fmax_c03_PRS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_PRS,fvalmax_1_c03_PRS,exitflagmax_1_c03_PRS]=fmincon(fmax_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_PRS,W1_ls_03_llm));

x_bnd_c03_PRS_1 = [xmin_1_c03_PRS;xmax_1_c03_PRS];

ind_c03_PRS = 2*(m_c03-n_ls_03);
x_bnd_c03_PRS_org = denormalizeuv(x_bnd_c03_PRS_1,lbx,ubx);
[psf_c03_PRS_2,pf_c03_PRS_2,re_c03_PRS_2,beta_c03_PRS_2] = mcspsfconstraint3c(x_bnd_c03_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);

std_c03_PRS = std([psf_c03_PRS_2;psf_c03_PRS_1]);
y_c03_PRS = constraint3c(x_bnd_c03_PRS_org);

%% Weight

%Finding t at which PSF is one
sp_t_obj_PRS = zeros(1,n_obj);
func_t_obj_PRS =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_PRS, fval_t_obj_PRS, exitflag_t_obj_PRS] = fmincon(func_t_obj_PRS,sp_t_obj_PRS,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_obj_1,@srgtsPRSEvaluate,mean(y_obj_new)));

%Validation
func_obj_PRS_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_PRS_1, fval_obj_PRS_1, exitflag_obj_PRS_1] = fmincon(func_obj_PRS_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_PRS));
xopt_obj_PRS_1_denorm = denormalizeuv(xopt_obj_PRS_1,lbx_obj,ubx_obj); 
xopt_obj_PRS_1_org = xopt_obj_PRS_1_denorm;
obj_PRS_1 = Weightc(xopt_obj_PRS_1_org);
obj_PRS_1_mm = srgtsPRSEvaluate(t_obj_PRS,srgt_PRS_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);

sp_obj_PRS = zeros(1,m_obj);
fmin_obj_PRS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_PRS,fvalmin_1_obj_PRS,exitflagmin_1_obj_PRS]=fmincon(fmin_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_PRS,W1_obj_llm));     

fmax_obj_PRS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_PRS,fvalmax_1_obj_PRS,exitflagmax_1_obj_PRS]=fmincon(fmax_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_PRS,W1_obj_llm));

x_bnd_obj_PRS_1 = [xmin_1_obj_PRS;xmax_1_obj_PRS];

ind_obj_PRS = 2*(m_obj-n_obj);
x_bnd_obj_PRS_org = denormalizeuv(x_bnd_obj_PRS_1,lbx_obj,ubx_obj);
obj_PRS_2 = Weightc(x_bnd_obj_PRS_org);

std_obj_PRS = std([obj_PRS_2;obj_PRS_1]);
y_obj_PRS = Weightc(x_bnd_obj_PRS_org);

%% WAS

clc;close all

%% Constraint 1

%Finding t at which PSF is one
sp_t_c01_WAS = zeros(1,n_ls_01);
func_t_c01_WAS =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_WAS, fval_t_c01_WAS, exitflag_t_c01_WAS] = fmincon(func_t_c01_WAS,sp_t_c01_WAS,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c01,@srgtsWASEvaluate,1));

%Validation
func_c01_WAS_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_WAS_1, fval_c01_WAS_1, exitflag_c01_WAS_1] = fmincon(func_c01_WAS_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_WAS));
xopt_c01_WAS_1_denorm = denormalizeuv(xopt_c01_WAS_1,lbx_c01,ubx_c01);
xopt_c01_WAS_1_org = xopt_c01_WAS_1_denorm;
[psf_c01_WAS_1,pf_c01_WAS_1,re_c01_WAS_1,beta_c01_WAS_1] = mcspsfconstraint1c(xopt_c01_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c01_WAS_1_mm = srgtsWASEvaluate(t_c01_WAS,srgt_WAS_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);

sp_c01_WAS = zeros(1,m_c01);
fmin_c01_WAS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_WAS,fvalmin_1_c01_WAS,exitflagmin_1_c01_WAS]=fmincon(fmin_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm));

fmax_c01_WAS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_WAS,fvalmax_1_c01_WAS,exitflagmax_1_c01_WAS]=fmincon(fmax_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm));

x_bnd_c01_WAS_1 = [xmin_1_c01_WAS;xmax_1_c01_WAS];

ind_c01_WAS = 2*(m_c01-n_ls_01);
x_bnd_c01_WAS_org = denormalizeuv(x_bnd_c01_WAS_1,lbx,ubx);
[psf_c01_WAS_2,pf_c01_WAS_2,re_c01_WAS_2,beta_c01_WAS_2] = mcspsfconstraint1c(x_bnd_c01_WAS_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c01_WAS = mean([psf_c01_WAS_2;psf_c01_WAS_1]);
std_c01_WAS = std([psf_c01_WAS_2;psf_c01_WAS_1]);
y_c01_WAS = constraint1c(x_bnd_c01_WAS_org);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_WAS = zeros(1,n_ls_02);
func_t_c02_WAS =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_WAS, fval_t_c02_WAS, exitflag_t_c02_WAS] = fmincon(func_t_c02_WAS,sp_t_c02_WAS,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c02,@srgtsWASEvaluate,1));

%Validation
func_c02_WAS_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_WAS_1, fval_c02_WAS_1, exitflag_c02_WAS_1] = fmincon(func_c02_WAS_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_WAS));
xopt_c02_WAS_1_denorm = denormalizeuv(xopt_c02_WAS_1,lbx_c02,ubx_c02);
xopt_c02_WAS_1_org = xopt_c02_WAS_1_denorm;
[psf_c02_WAS_1,pf_c02_WAS_1,re_c02_WAS_1,beta_c02_WAS_1] = mcspsfconstraint2c(xopt_c02_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_WAS_1_mm = srgtsWASEvaluate(t_c02_WAS,srgt_WAS_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);

sp_c02_WAS = zeros(1,m_c02);
fmin_c02_WAS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_WAS,fvalmin_1_c02_WAS,exitflagmin_1_c02_WAS]=fmincon(fmin_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm));

fmax_c02_WAS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_WAS,fvalmax_1_c02_WAS,exitflagmax_1_c02_WAS]=fmincon(fmax_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm));

x_bnd_c02_WAS_1 = [xmin_1_c02_WAS;xmax_1_c02_WAS];

ind_c02_WAS = 2*(m_c02-n_ls_02);
x_bnd_c02_WAS_org = denormalizeuv(x_bnd_c02_WAS_1,lbx,ubx);
[psf_c02_WAS_2,pf_c02_WAS_2,re_c02_WAS_2,beta_c02_WAS_2] = mcspsfconstraint2c(x_bnd_c02_WAS_org,lbx,ubx,sd,Nmcs,Pftarget);

mean_c02_WAS = mean([psf_c02_WAS_2;psf_c02_WAS_1]);
std_c02_WAS = std([psf_c02_WAS_2;psf_c02_WAS_1]);
y_c02_WAS = constraint2c(x_bnd_c02_WAS_org);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_WAS = zeros(1,n_ls_03);
func_t_c03_WAS =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_WAS, fval_t_c03_WAS, exitflag_t_c03_WAS] = fmincon(func_t_c03_WAS,sp_t_c03_WAS,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c03,@srgtsWASEvaluate,1));

%Validation
func_c03_WAS_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_WAS_1, fval_c03_WAS_1, exitflag_c03_WAS_1] = fmincon(func_c03_WAS_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_WAS));
xopt_c03_WAS_1_denorm = denormalizeuv(xopt_c03_WAS_1,lbx_c03,ubx_c03); 
xopt_c03_WAS_1_org = xopt_c03_WAS_1_denorm;
[psf_c03_WAS_1,pf_c03_WAS_1,re_c03_WAS_1,beta_c03_WAS_1] = mcspsfconstraint3c(xopt_c03_WAS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c03_WAS_1_mm = srgtsWASEvaluate(t_c03_WAS,srgt_WAS_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);


sp_c03_WAS = zeros(1,m_c03);
fmin_c03_WAS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_WAS,fvalmin_1_c03_WAS,exitflagmin_1_c03_WAS]=fmincon(fmin_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm));

fmax_c03_WAS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_WAS,fvalmax_1_c03_WAS,exitflagmax_1_c03_WAS]=fmincon(fmax_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm));

x_bnd_c03_WAS_1 = [xmin_1_c03_WAS;xmax_1_c03_WAS];

ind_c03_WAS = 2*(m_c03-n_ls_03);
x_bnd_c03_WAS_org = denormalizeuv(x_bnd_c03_WAS_1,lbx,ubx);
[psf_c03_WAS_2,pf_c03_WAS_2,re_c03_WAS_2,beta_c03_WAS_2] = mcspsfconstraint3c(x_bnd_c03_WAS_org,lbx,ubx,sd,Nmcs,Pftarget);

std_c03_WAS = std([psf_c03_WAS_2;psf_c03_WAS_1]);
y_c03_WAS = constraint3c(x_bnd_c03_WAS_org);

%% Weight

%Finding t at which PSF is one
sp_t_obj_WAS = zeros(1,n_obj);
func_t_obj_WAS =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_WAS, fval_t_obj_WAS, exitflag_t_obj_WAS] = fmincon(func_t_obj_WAS,sp_t_obj_WAS,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_obj_1,@srgtsWASEvaluate,mean(y_obj_new)));

%Validation
func_obj_WAS_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_WAS_1, fval_obj_WAS_1, exitflag_obj_WAS_1] = fmincon(func_obj_WAS_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_WAS));
xopt_obj_WAS_1_denorm = denormalizeuv(xopt_obj_WAS_1,lbx_obj,ubx_obj); 
xopt_obj_WAS_1_org = xopt_obj_WAS_1_denorm;
obj_WAS_1 = Weightc(xopt_obj_WAS_1_org);
obj_WAS_1_mm = srgtsWASEvaluate(t_obj_WAS,srgt_WAS_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);

sp_obj_WAS = zeros(1,m_obj);
fmin_obj_WAS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_WAS,fvalmin_1_obj_WAS,exitflagmin_1_obj_WAS]=fmincon(fmin_obj_WAS_1,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm));     

fmax_obj_WAS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_WAS,fvalmax_1_obj_WAS,exitflagmax_1_obj_WAS]=fmincon(fmax_obj_WAS_1,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm));

x_bnd_obj_WAS_1 = [xmin_1_obj_WAS;xmax_1_obj_WAS];

ind_obj_WAS = 2*(m_obj-n_obj);
x_bnd_obj_WAS_org = denormalizeuv(x_bnd_obj_WAS_1,lbx_obj,ubx_obj);
obj_WAS_2 = Weightc(x_bnd_obj_WAS_org);

std_obj_WAS = std([obj_WAS_2;obj_WAS_1]);
y_obj_WAS = Weightc(x_bnd_obj_WAS_org);


