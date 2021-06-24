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
xopt_c01_KRG_1_org = [0 xopt_c01_KRG_1_denorm(:,1:3) 0 xopt_c01_KRG_1_denorm(:,4) 0];
[psf_c01_KRG_1,pf_c01_KRG_1,re_c01_KRG_1,beta_c01_KRG_1] = mcspsfconstraint1c(xopt_c01_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c01_KRG_1_mm = srgtsKRGEvaluate(t_c01_KRG,srgt_KRG_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);
t2_ls_01_KRG = xopt_c01_KRG_1*W2_ls_01_llm;

sp_c01_KRG = zeros(1,m_c01);
fmin_c01_KRG_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_KRG,fvalmin_1_c01_KRG,exitflagmin_1_c01_KRG]=fmincon(fmin_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,2),W2_ls_01_llm(:,3),t2_ls_01_KRG(2),t2_ls_01_KRG(3)));

fmax_c01_KRG_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_KRG,fvalmax_1_c01_KRG,exitflagmax_1_c01_KRG]=fmincon(fmax_c01_KRG_1,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,2),W2_ls_01_llm(:,3),t2_ls_01_KRG(2),t2_ls_01_KRG(3)));   

fmin_c01_KRG_2 = @(xmin_c01_2)(xmin_c01_2*W2_ls_01_llm(:,2));
[xmin_2_c01_KRG,fvalmin_2_c01_KRG,exitflagmin_2_c01_KRG]=fmincon(fmin_c01_KRG_2,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,3),t2_ls_01_KRG(1),t2_ls_01_KRG(3)));    

fmax_c01_KRG_2 = @(xmax_c01_2)(-(xmax_c01_2*W2_ls_01_llm(:,2)));
[xmax_2_c01_KRG,fvalmax_2_c01_KRG,exitflagmax_2_c01_KRG]=fmincon(fmax_c01_KRG_2,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,3),t2_ls_01_KRG(1),t2_ls_01_KRG(3)));    

fmin_c01_KRG_3 = @(xmin_c01_3)(xmin_c01_3*W2_ls_01_llm(:,3));
[xmin_3_c01_KRG,fvalmin_3_c01_KRG,exitflagmin_3_c01_KRG]=fmincon(fmin_c01_KRG_3,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,2),t2_ls_01_KRG(1),t2_ls_01_KRG(2)));       

fmax_c01_KRG_3 = @(xmax_c01_3)(-(xmax_c01_3*W2_ls_01_llm(:,3)));
[xmax_3_c01_KRG,fvalmax_3_c01_KRG,exitflagmax_3_c01_KRG]=fmincon(fmax_c01_KRG_3,sp_c01_KRG,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_KRG,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,2),t2_ls_01_KRG(1),t2_ls_01_KRG(2)));    

x_bnd_c01_KRG_1 = [xmin_1_c01_KRG;xmax_1_c01_KRG;xmin_2_c01_KRG;xmax_2_c01_KRG;xmin_3_c01_KRG;xmax_3_c01_KRG];

ind_c01_KRG = 2*(m_c01-n_ls_01);
x_bnd_c01_KRG_org = denormalizeuv([zeros(ind_c01_KRG,1) x_bnd_c01_KRG_1(:,1:3) zeros(ind_c01_KRG,1) x_bnd_c01_KRG_1(:,4) zeros(ind_c01_KRG,1)],lbx,ubx);
[psf_c01_KRG_2,pf_c01_KRG_2,re_c01_KRG_2,beta_c01_KRG_2] = mcspsfconstraint1c(x_bnd_c01_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

mean_c01_KRG = mean([psf_c01_KRG_2;psf_c01_KRG_1]);
std_c01_KRG = std([psf_c01_KRG_2;psf_c01_KRG_1]);
y_c01_KRG = constraint1c([x_bnd_c01_KRG_org(:,2:4) x_bnd_c01_KRG_org(:,6) repmat([mup2 mup3 mup4],ind_c01_KRG,1)]);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_KRG = zeros(1,n_ls_02);
func_t_c02_KRG =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_KRG, fval_t_c02_KRG, exitflag_t_c02_KRG] = fmincon(func_t_c02_KRG,sp_t_c02_KRG,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c02,@srgtsKRGEvaluate,1));

%Validation
func_c02_KRG_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_KRG_1, fval_c02_KRG_1, exitflag_c02_KRG_1] = fmincon(func_c02_KRG_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_KRG));
xopt_c02_KRG_1_denorm = denormalizeuv(xopt_c02_KRG_1,lbx_c02,ubx_c02);
xopt_c02_KRG_1_org = [xopt_c02_KRG_1_denorm(:,1:3) zeros(1,4)];
[psf_c02_KRG_1,pf_c02_KRG_1,re_c02_KRG_1,beta_c02_KRG_1] = mcspsfconstraint2c(xopt_c02_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c02_KRG_1_mm = srgtsKRGEvaluate(t_c02_KRG,srgt_KRG_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);
t2_ls_02_llm = xopt_c02_KRG_1*W2_ls_02_llm;

sp_c02_KRG = zeros(1,m_c02);
fmin_c02_KRG_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_KRG,fvalmin_1_c02_KRG,exitflagmin_1_c02_KRG]=fmincon(fmin_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm,W2_ls_02_llm(:,2),t2_ls_02_llm(2)));

fmax_c02_KRG_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_KRG,fvalmax_1_c02_KRG,exitflagmax_1_c02_KRG]=fmincon(fmax_c02_KRG_1,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm,W2_ls_02_llm(:,2),t2_ls_02_llm(2)));

fmin_c02_KRG_2 = @(xmin_c02_2)(xmin_c02_2*W2_ls_02_llm(:,2));
[xmin_2_c02_KRG,fvalmin_2_c02_KRG,exitflagmin_2_c02_KRG]=fmincon(fmin_c02_KRG_2,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm,W2_ls_02_llm(:,1),t2_ls_02_llm(1)));

fmax_c02_KRG_2 = @(xmax_c02_2)(-(xmax_c02_2*W2_ls_02_llm(:,2)));
[xmax_2_c02_KRG,fvalmax_2_c02_KRG,exitflagmax_2_c02_KRG]=fmincon(fmax_c02_KRG_2,sp_c02_KRG,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_KRG,W1_ls_02_llm,W2_ls_02_llm(:,1),t2_ls_02_llm(1)));


x_bnd_c02_KRG_1 = [xmin_1_c02_KRG;xmax_1_c02_KRG;xmin_2_c02_KRG;xmax_2_c02_KRG];

ind_c02_KRG = 2*(m_c02-n_ls_02);
x_bnd_c02_KRG_org = denormalizeuv([x_bnd_c02_KRG_1(:,1:3) zeros(ind_c02_KRG,4)],lbx,ubx);
[psf_c02_KRG_2,pf_c02_KRG_2,re_c02_KRG_2,beta_c02_KRG_2] = mcspsfconstraint2c(x_bnd_c02_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

mean_c02_KRG = mean([psf_c02_KRG_2;psf_c02_KRG_1]);
std_c02_KRG = std([psf_c02_KRG_2;psf_c02_KRG_1]);
y_c02_KRG = constraint2c([x_bnd_c02_KRG_org(:,1:3) repmat([mup1 mup3],ind_c02_KRG,1)]);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_KRG = zeros(1,n_ls_03);
func_t_c03_KRG =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_KRG, fval_t_c03_KRG, exitflag_t_c03_KRG] = fmincon(func_t_c03_KRG,sp_t_c03_KRG,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c03,@srgtsKRGEvaluate,1));

%Validation
func_c03_KRG_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_KRG_1, fval_c03_KRG_1, exitflag_c03_KRG_1] = fmincon(func_c03_KRG_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_KRG));
xopt_c03_KRG_1_denorm = denormalizeuv(xopt_c03_KRG_1,lbx_c03,ubx_c03); 
xopt_c03_KRG_1_org = [xopt_c03_KRG_1_denorm(:,1:3) 0 xopt_c03_KRG_1_denorm(:,4) 0 xopt_c03_KRG_1_denorm(:,5)];
[psf_c03_KRG_1,pf_c03_KRG_1,re_c03_KRG_1,beta_c03_KRG_1] = mcspsfconstraint3c(xopt_c03_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c03_KRG_1_mm = srgtsKRGEvaluate(t_c03_KRG,srgt_KRG_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);
t2_ls_03_llm = xopt_c03_KRG_1*W2_ls_03_llm;

sp_c03_KRG = zeros(1,m_c03);
fmin_c03_KRG_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_KRG,fvalmin_1_c03_KRG,exitflagmin_1_c03_KRG]=fmincon(fmin_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(2),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmax_c03_KRG_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_KRG,fvalmax_1_c03_KRG,exitflagmax_1_c03_KRG]=fmincon(fmax_c03_KRG_1,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(2),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmin_c03_KRG_2 = @(xmin_c03_2)(xmin_c03_2*W2_ls_03_llm(:,2));
[xmin_2_c03_KRG,fvalmin_2_c03_KRG,exitflagmin_2_c03_KRG]=fmincon(fmin_c03_KRG_2,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmax_c03_KRG_2 = @(xmax_c03_2)(-(xmax_c03_2*W2_ls_03_llm(:,2)));
[xmax_2_c03_KRG,fvalmax_2_c03_KRG,exitflagmax_2_c03_KRG]=fmincon(fmax_c03_KRG_2,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmin_c03_KRG_3 = @(xmin_c03_3)(xmin_c03_3*W2_ls_03_llm(:,3));
[xmin_3_c03_KRG,fvalmin_3_c03_KRG,exitflagmin_3_c03_KRG]=fmincon(fmin_c03_KRG_3,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(4)));

fmax_c03_KRG_3 = @(xmax_c03_3)(-(xmax_c03_3*W2_ls_03_llm(:,3)));
[xmax_3_c03_KRG,fvalmax_3_c03_KRG,exitflagmax_3_c03_KRG]=fmincon(fmax_c03_KRG_3,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(4)));

fmin_c03_KRG_4 = @(xmin_c03_4)(xmin_c03_4*W2_ls_03_llm(:,4));
[xmin_4_c03_KRG,fvalmin_4_c03_KRG,exitflagmin_4_c03_KRG]=fmincon(fmin_c03_KRG_4,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(3)));

fmax_c03_KRG_4 = @(xmax_c03_4)(-(xmax_c03_4*W2_ls_03_llm(:,4)));
[xmax_4_c03_KRG,fvalmax_4_c03_KRG,exitflagmax_4_c03_KRG]=fmincon(fmax_c03_KRG_4,sp_c03_KRG,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_KRG,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(3)));


x_bnd_c03_KRG_1 = [xmin_1_c03_KRG;xmax_1_c03_KRG;xmin_2_c03_KRG;xmax_2_c03_KRG;xmin_3_c03_KRG;xmax_3_c03_KRG;xmin_4_c03_KRG;xmax_4_c03_KRG];

ind_c03_KRG = 2*(m_c03-n_ls_03);
x_bnd_c03_KRG_org = denormalizeuv([x_bnd_c03_KRG_1(:,1:3) zeros(ind_c03_KRG,1) x_bnd_c03_KRG_1(:,4) zeros(ind_c03_KRG,1) x_bnd_c03_KRG_1(:,5)],lbx,ubx);
[psf_c03_KRG_2,pf_c03_KRG_2,re_c03_KRG_2,beta_c03_KRG_2] = mcspsfconstraint3c(x_bnd_c03_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c03_KRG = std([psf_c03_KRG_2;psf_c03_KRG_1]);
y_c03_KRG = constraint3c([x_bnd_c03_KRG_org repmat([mup1 mup2 mup3],ind_c03_KRG,1)]);

%% Constraint 4

%Finding t at which PSF is one
sp_t_c04_KRG = zeros(1,n_ls_04);
func_t_c04_KRG =  @(t_opt) (t_opt*zeros(n_ls_04,1));
[t_c04_KRG, fval_t_c04_KRG, exitflag_t_c04_KRG] = fmincon(func_t_c04_KRG,sp_t_c04_KRG,[],[],[],[],t1_ls_04_min ,t1_ls_04_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c04,@srgtsKRGEvaluate,1));

%Validation
func_c04_KRG_1 =  @(x_opt) (x_opt*zeros(m_c04,1));
[xopt_c04_KRG_1, fval_c04_KRG_1, exitflag_c04_KRG_1] = fmincon(func_c04_KRG_1,sp_c04,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x_opt)nonlconstage1(x_opt,W1_ls_04_llm,t_c04_KRG));
xopt_c04_KRG_1_denorm = denormalizeuv(xopt_c04_KRG_1,lbx_c04,ubx_c04); 
xopt_c04_KRG_1_org = [xopt_c04_KRG_1_denorm(:,1:3) 0 xopt_c04_KRG_1_denorm(:,4:6)];
[psf_c04_KRG_1,pf_c04_KRG_1,re_c04_KRG_1,beta_c04_KRG_1] = mcspsfconstraint4c(xopt_c04_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c04_KRG_1_mm = srgtsKRGEvaluate(t_c04_KRG,srgt_KRG_c04);

%Bounds in inactive subspace
W2_ls_04_llm = W_ls_04_llm(:,n_ls_04+1:m_c04);
t2_ls_04_llm = xopt_c04_KRG_1*W2_ls_04_llm;

sp_c04_KRG = zeros(1,m_c04);
fmin_c04_KRG_1 = @(xmin_c04_1)(xmin_c04_1*W2_ls_04_llm(:,1));
[xmin_1_c04_KRG,fvalmin_1_c04_KRG,exitflagmin_1_c04_KRG]=fmincon(fmin_c04_KRG_1,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_KRG_1 = @(xmax_c04_1)(-(xmax_c04_1*W2_ls_04_llm(:,1)));
[xmax_1_c04_KRG,fvalmax_1_c04_KRG,exitflagmax_1_c04_KRG]=fmincon(fmax_c04_KRG_1,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_KRG_2 = @(xmin_c04_2)(xmin_c04_2*W2_ls_04_llm(:,2));
[xmin_2_c04_KRG,fvalmin_2_c04_KRG,exitflagmin_2_c04_KRG]=fmincon(fmin_c04_KRG_2,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_KRG_2 = @(xmax_c04_2)(-(xmax_c04_2*W2_ls_04_llm(:,2)));
[xmax_2_c04_KRG,fvalmax_2_c04_KRG,exitflagmax_2_c04_KRG]=fmincon(fmax_c04_KRG_2,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_KRG_3 = @(xmin_c04_3)(xmin_c04_3*W2_ls_04_llm(:,3));
[xmin_3_c04_KRG,fvalmin_3_c04_KRG,exitflagmin_3_c04_KRG]=fmincon(fmin_c04_KRG_3,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_KRG_3 = @(xmax_c04_3)(-(xmax_c04_3*W2_ls_04_llm(:,3)));
[xmax_3_c04_KRG,fvalmax_3_c04_KRG,exitflagmax_3_c04_KRG]=fmincon(fmax_c04_KRG_3,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_KRG_4 = @(xmin_c04_4)(xmin_c04_4*W2_ls_04_llm(:,4));
[xmin_4_c04_KRG,fvalmin_4_c04_KRG,exitflagmin_4_c04_KRG]=fmincon(fmin_c04_KRG_4,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(5)));

fmax_c04_KRG_4 = @(xmax_c04_4)(-(xmax_c04_4*W2_ls_04_llm(:,4)));
[xmax_4_c04_KRG,fvalmax_4_c04_KRG,exitflagmax_4_c04_KRG]=fmincon(fmax_c04_KRG_4,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(5)));

fmin_c04_KRG_5 = @(xmin_c04_5)(xmin_c04_5*W2_ls_04_llm(:,5));
[xmin_5_c04_KRG,fvalmin_5_c04_KRG,exitflagmin_5_c04_KRG]=fmincon(fmin_c04_KRG_5,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4)));

fmax_c04_KRG_5 = @(xmax_c04_5)(-(xmax_c04_5*W2_ls_04_llm(:,5)));
[xmax_5_c04_KRG,fvalmax_5_c04_KRG,exitflagmax_5_c04_KRG]=fmincon(fmax_c04_KRG_5,sp_c04_KRG,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_KRG,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4)));

x_bnd_c04_KRG_1 = [xmin_1_c04_KRG;xmax_1_c04_KRG;xmin_2_c04_KRG;xmax_2_c04_KRG;xmin_3_c04_KRG;xmax_3_c04_KRG;xmin_4_c04_KRG;xmax_4_c04_KRG;xmin_5_c04_KRG;xmax_5_c04_KRG];

ind_c04_KRG = 2*(m_c04-n_ls_04);
x_bnd_c04_KRG_org = denormalizeuv([x_bnd_c04_KRG_1(:,1:3) zeros(ind_c04_KRG,1) x_bnd_c04_KRG_1(:,4:6)],lbx,ubx);
[psf_c04_KRG_2,pf_c04_KRG_2,re_c04_KRG_2,beta_c04_KRG_2] = mcspsfconstraint4c(x_bnd_c04_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c04_KRG = std([psf_c04_KRG_2;psf_c04_KRG_1]);
y_c04_KRG = constraint4c([x_bnd_c04_KRG_org repmat([mup1 mup2 mup3],ind_c04_KRG,1)]);

%% Constraint 5
%Finding t at which PSF is one
sp_t_c05_KRG = zeros(1,n_ls_05);
func_t_c05_KRG =  @(t_opt) (t_opt*zeros(n_ls_05,1));
[t_c05_KRG, fval_t_c05_KRG, exitflag_t_c05_KRG] = fmincon(func_t_c05_KRG,sp_t_c05_KRG,[],[],[],[],t1_ls_05_min ,t1_ls_05_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c05,@srgtsKRGEvaluate,1));

%Validation
func_c05_KRG_1 =  @(x_opt) (x_opt*zeros(m_c05,1));
[xopt_c05_KRG_1, fval_c05_KRG_1, exitflag_c05_KRG_1] = fmincon(func_c05_KRG_1,sp_c05,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x_opt)nonlconstage1(x_opt,W1_ls_05_llm,t_c05_KRG));
xopt_c05_KRG_1_denorm = denormalizeuv(xopt_c05_KRG_1,lbx_c05,ubx_c05); 
xopt_c05_KRG_1_org = [0 xopt_c05_KRG_1_denorm(:,1:2) 0 0 0 xopt_c05_KRG_1_denorm(:,3)];
[psf_c05_KRG_1,pf_c05_KRG_1,re_c05_KRG_1,beta_c05_KRG_1] = mcspsfconstraint5c(xopt_c05_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c05_KRG_1_mm = srgtsKRGEvaluate(t_c05_KRG,srgt_KRG_c05);

%Bounds in inactive subspace
W2_ls_05_llm = W_ls_05_llm(:,n_ls_05+1:m_c05);
t2_ls_05_llm = xopt_c05_KRG_1*W2_ls_05_llm;

sp_c05_KRG = zeros(1,m_c05);
fmin_c05_KRG_1 = @(xmin_c05_1)(xmin_c05_1*W2_ls_05_llm(:,1));
[xmin_1_c05_KRG,fvalmin_1_c05_KRG,exitflagmin_1_c05_KRG]=fmincon(fmin_c05_KRG_1,sp_c05_KRG,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_KRG,W1_ls_05_llm,W2_ls_05_llm(:,2),t2_ls_05_llm(:,2)));

fmax_c05_KRG_1 = @(xmax_c05_1)(-(xmax_c05_1*W2_ls_05_llm(:,1)));
[xmax_1_c05_KRG,fvalmax_1_c05_KRG,exitflagmax_1_c05_KRG]=fmincon(fmax_c05_KRG_1,sp_c05_KRG,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_KRG,W1_ls_05_llm,W2_ls_05_llm(:,2),t2_ls_05_llm(:,2)));

fmin_c05_KRG_2 = @(xmin_c05_2)(xmin_c05_2*W2_ls_05_llm(:,2));
[xmin_2_c05_KRG,fvalmin_2_c05_KRG,exitflagmin_2_c05_KRG]=fmincon(fmin_c05_KRG_2,sp_c05_KRG,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_KRG,W1_ls_05_llm,W2_ls_05_llm(:,1),t2_ls_05_llm(:,1)));

fmax_c05_KRG_2 = @(xmax_c05_2)(-(xmax_c05_2*W2_ls_05_llm(:,2)));
[xmax_2_c05_KRG,fvalmax_2_c05_KRG,exitflagmax_2_c05_KRG]=fmincon(fmax_c05_KRG_2,sp_c05_KRG,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_KRG,W1_ls_05_llm,W2_ls_05_llm(:,1),t2_ls_05_llm(:,1)));


x_bnd_c05_KRG_1 = [xmin_1_c05_KRG;xmax_1_c05_KRG;xmin_2_c05_KRG;xmax_2_c05_KRG];

ind_c05_KRG = 2*(m_c05-n_ls_05);

x_bnd_c05_KRG_org = denormalizeuv([zeros(ind_c05_KRG,1) x_bnd_c05_KRG_1(:,1:2) zeros(ind_c05_KRG,3) x_bnd_c05_KRG_1(:,3)],lbx,ubx);
[psf_c05_KRG_2,pf_c05_KRG_2,re_c05_KRG_2,beta_c05_KRG_2] = mcspsfconstraint5c(x_bnd_c05_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c05_KRG = std([psf_c05_KRG_2;psf_c05_KRG_1]);
y_c05_KRG = constraint5c([x_bnd_c05_KRG_org repmat([mup1 mup2 mup3],ind_c05_KRG,1)]);

%% Constraint 6
%Finding t at which PSF is one
sp_t_c06_KRG = zeros(1,n_ls_06);
func_t_c06_KRG =  @(t_opt) (t_opt*zeros(n_ls_06,1));
[t_c06_KRG, fval_t_c06_KRG, exitflag_t_c06_KRG] = fmincon(func_t_c06_KRG,sp_t_c06_KRG,[],[],[],[],t1_ls_06_min ,t1_ls_06_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c06,@srgtsKRGEvaluate,mean(psf_c06)));

%Validation
func_c06_KRG_1 =  @(x_opt) (x_opt*zeros(m_c06,1));
[xopt_c06_KRG_1, fval_c06_KRG_1, exitflag_c06_KRG_1] = fmincon(func_c06_KRG_1,sp_c06,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x_opt)nonlconstage1(x_opt,W1_ls_06_llm,t_c06_KRG));
xopt_c06_KRG_1_denorm = denormalizeuv(xopt_c06_KRG_1,lbx_c06,ubx_c06); 
xopt_c06_KRG_1_org = [xopt_c06_KRG_1_denorm(:,1:3) 0 xopt_c06_KRG_1_denorm(:,4:6)];
[psf_c06_KRG_1,pf_c06_KRG_1,re_c06_KRG_1,beta_c06_KRG_1] = mcspsfconstraint6c(xopt_c06_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c06_KRG_1_mm = srgtsKRGEvaluate(t_c06_KRG,srgt_KRG_c06);

%Bounds in inactive subspace
W2_ls_06_llm = W_ls_06_llm(:,n_ls_06+1:m_c06);
t2_ls_06_llm = xopt_c06_KRG_1*W2_ls_06_llm;

sp_c06_KRG = zeros(1,m_c06);
fmin_c06_KRG_1 = @(xmin_c06_1)(xmin_c06_1*W2_ls_06_llm(:,1));
[xmin_1_c06_KRG,fvalmin_1_c06_KRG,exitflagmin_1_c06_KRG]=fmincon(fmin_c06_KRG_1,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_KRG_1 = @(xmax_c06_1)(-(xmax_c06_1*W2_ls_06_llm(:,1)));
[xmax_1_c06_KRG,fvalmax_1_c06_KRG,exitflagmax_1_c06_KRG]=fmincon(fmax_c06_KRG_1,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_KRG_2 = @(xmin_c06_2)(xmin_c06_2*W2_ls_06_llm(:,2));
[xmin_2_c06_KRG,fvalmin_2_c06_KRG,exitflagmin_2_c06_KRG]=fmincon(fmin_c06_KRG_2,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_KRG_2 = @(xmax_c06_2)(-(xmax_c06_2*W2_ls_06_llm(:,2)));
[xmax_2_c06_KRG,fvalmax_2_c06_KRG,exitflagmax_2_c06_KRG]=fmincon(fmax_c06_KRG_2,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_KRG_3 = @(xmin_c06_3)(xmin_c06_3*W2_ls_06_llm(:,3));
[xmin_3_c06_KRG,fvalmin_3_c06_KRG,exitflagmin_3_c06_KRG]=fmincon(fmin_c06_KRG_3,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_KRG_3 = @(xmax_c06_3)(-(xmax_c06_3*W2_ls_06_llm(:,3)));
[xmax_3_c06_KRG,fvalmax_3_c06_KRG,exitflagmax_3_c06_KRG]=fmincon(fmax_c06_KRG_3,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_KRG_4 = @(xmin_c06_4)(xmin_c06_4*W2_ls_06_llm(:,4));
[xmin_4_c06_KRG,fvalmin_4_c06_KRG,exitflagmin_4_c06_KRG]=fmincon(fmin_c06_KRG_4,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(5)));

fmax_c06_KRG_4 = @(xmax_c06_4)(-(xmax_c06_4*W2_ls_06_llm(:,4)));
[xmax_4_c06_KRG,fvalmax_4_c06_KRG,exitflagmax_4_c06_KRG]=fmincon(fmax_c06_KRG_4,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(5)));

fmin_c06_KRG_5 = @(xmin_c06_5)(xmin_c06_5*W2_ls_06_llm(:,5));
[xmin_5_c06_KRG,fvalmin_5_c06_KRG,exitflagmin_5_c06_KRG]=fmincon(fmin_c06_KRG_5,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4)));

fmax_c06_KRG_5 = @(xmax_c06_5)(-(xmax_c06_5*W2_ls_06_llm(:,5)));
[xmax_5_c06_KRG,fvalmax_5_c06_KRG,exitflagmax_5_c06_KRG]=fmincon(fmax_c06_KRG_5,sp_c06_KRG,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_KRG,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4)));

x_bnd_c06_KRG_1 = [xmin_1_c06_KRG;xmax_1_c06_KRG;xmin_2_c06_KRG;xmax_2_c06_KRG;xmin_3_c06_KRG;xmax_3_c06_KRG;xmin_4_c06_KRG;xmax_4_c06_KRG;xmin_5_c06_KRG;xmax_5_c06_KRG];

ind_c06_KRG = 2*(m_c06-n_ls_06);
x_bnd_c06_KRG_org = denormalizeuv([x_bnd_c06_KRG_1(:,1:3) zeros(ind_c06_KRG,1) x_bnd_c06_KRG_1(:,4:6)],lbx,ubx);
[psf_c06_KRG_2,pf_c06_KRG_2,re_c06_KRG_2,beta_c06_KRG_2] = mcspsfconstraint6c(x_bnd_c06_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c06_KRG = std([psf_c06_KRG_2;psf_c06_KRG_1]);
y_c06_KRG = constraint6c([x_bnd_c06_KRG_org repmat([mup1 mup2 mup3 mup4],ind_c06_KRG,1)]);

%% Constraint 7
%Finding t at which PSF is one
sp_t_c07_KRG = zeros(1,n_ls_07);
func_t_c07_KRG =  @(t_opt) (t_opt*zeros(n_ls_07,1));
[t_c07_KRG, fval_t_c07_KRG, exitflag_t_c07_KRG] = fmincon(func_t_c07_KRG,sp_t_c07_KRG,[],[],[],[],t1_ls_07_min ,t1_ls_07_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c07,@srgtsKRGEvaluate,mean(psf_c07)));

%Validation
func_c07_KRG_1 =  @(x_opt) (x_opt*zeros(m_c07,1));
[xopt_c07_KRG_1, fval_c07_KRG_1, exitflag_c07_KRG_1] = fmincon(func_c07_KRG_1,sp_c07,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x_opt)nonlconstage1(x_opt,W1_ls_07_llm,t_c07_KRG));
xopt_c07_KRG_1_denorm = denormalizeuv(xopt_c07_KRG_1,lbx_c07,ubx_c07); 
xopt_c07_KRG_1_org = [xopt_c07_KRG_1_denorm(:,1:3) 0 xopt_c07_KRG_1_denorm(:,4:6)];
[psf_c07_KRG_1,pf_c07_KRG_1,re_c07_KRG_1,beta_c07_KRG_1] = mcspsfconstraint7c(xopt_c07_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c07_KRG_1_mm = srgtsKRGEvaluate(t_c07_KRG,srgt_KRG_c07);

%Bounds in inactive subspace
W2_ls_07_llm = W_ls_07_llm(:,n_ls_07+1:m_c07);
t2_ls_07_llm = xopt_c07_KRG_1*W2_ls_07_llm;

sp_c07_KRG = zeros(1,m_c07);
fmin_c07_KRG_1 = @(xmin_c07_1)(xmin_c07_1*W2_ls_07_llm(:,1));
[xmin_1_c07_KRG,fvalmin_1_c07_KRG,exitflagmin_1_c07_KRG]=fmincon(fmin_c07_KRG_1,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_KRG_1 = @(xmax_c07_1)(-(xmax_c07_1*W2_ls_07_llm(:,1)));
[xmax_1_c07_KRG,fvalmax_1_c07_KRG,exitflagmax_1_c07_KRG]=fmincon(fmax_c07_KRG_1,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_KRG_2 = @(xmin_c07_2)(xmin_c07_2*W2_ls_07_llm(:,2));
[xmin_2_c07_KRG,fvalmin_2_c07_KRG,exitflagmin_2_c07_KRG]=fmincon(fmin_c07_KRG_2,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_KRG_2 = @(xmax_c07_2)(-(xmax_c07_2*W2_ls_07_llm(:,2)));
[xmax_2_c07_KRG,fvalmax_2_c07_KRG,exitflagmax_2_c07_KRG]=fmincon(fmax_c07_KRG_2,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_KRG_3 = @(xmin_c07_3)(xmin_c07_3*W2_ls_07_llm(:,3));
[xmin_3_c07_KRG,fvalmin_3_c07_KRG,exitflagmin_3_c07_KRG]=fmincon(fmin_c07_KRG_3,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_KRG_3 = @(xmax_c07_3)(-(xmax_c07_3*W2_ls_07_llm(:,3)));
[xmax_3_c07_KRG,fvalmax_3_c07_KRG,exitflagmax_3_c07_KRG]=fmincon(fmax_c07_KRG_3,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_KRG_4 = @(xmin_c07_4)(xmin_c07_4*W2_ls_07_llm(:,4));
[xmin_4_c07_KRG,fvalmin_4_c07_KRG,exitflagmin_4_c07_KRG]=fmincon(fmin_c07_KRG_4,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(5)));

fmax_c07_KRG_4 = @(xmax_c07_4)(-(xmax_c07_4*W2_ls_07_llm(:,4)));
[xmax_4_c07_KRG,fvalmax_4_c07_KRG,exitflagmax_4_c07_KRG]=fmincon(fmax_c07_KRG_4,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(5)));

fmin_c07_KRG_5 = @(xmin_c07_5)(xmin_c07_5*W2_ls_07_llm(:,5));
[xmin_5_c07_KRG,fvalmin_5_c07_KRG,exitflagmin_5_c07_KRG]=fmincon(fmin_c07_KRG_5,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4)));

fmax_c07_KRG_5 = @(xmax_c07_5)(-(xmax_c07_5*W2_ls_07_llm(:,5)));
[xmax_5_c07_KRG,fvalmax_5_c07_KRG,exitflagmax_5_c07_KRG]=fmincon(fmax_c07_KRG_5,sp_c07_KRG,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_KRG,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4)));

x_bnd_c07_KRG_1 = [xmin_1_c07_KRG;xmax_1_c07_KRG;xmin_2_c07_KRG;xmax_2_c07_KRG;xmin_3_c07_KRG;xmax_3_c07_KRG;xmin_4_c07_KRG;xmax_4_c07_KRG;xmin_5_c07_KRG;xmax_5_c07_KRG];

ind_c07_KRG = 2*(m_c07-n_ls_07);
x_bnd_c07_KRG_org = denormalizeuv([x_bnd_c07_KRG_1(:,1:3) zeros(ind_c07_KRG,1) x_bnd_c07_KRG_1(:,4:6)],lbx,ubx);
[psf_c07_KRG_2,pf_c07_KRG_2,re_c07_KRG_2,beta_c07_KRG_2] = mcspsfconstraint7c(x_bnd_c07_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c07_KRG = std([psf_c07_KRG_2;psf_c07_KRG_1]);
y_c07_KRG = constraint7c([x_bnd_c07_KRG_org repmat([mup1 mup2 mup3 mup4],ind_c07_KRG,1)]);

%% Constraint 8
%Finding t at which PSF is one
sp_t_c08_KRG = zeros(1,n_ls_08);
func_t_c08_KRG =  @(t_opt) (t_opt*zeros(n_ls_08,1));
[t_c08_KRG, fval_t_c08_KRG, exitflag_t_c08_KRG] = fmincon(func_t_c08_KRG,sp_t_c08_KRG,[],[],[],[],t1_ls_08_min ,t1_ls_08_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c08,@srgtsKRGEvaluate,0.9822));

%Validation
func_c08_KRG_1 =  @(x_opt) (x_opt*zeros(m_c08,1));
[xopt_c08_KRG_1, fval_c08_KRG_1, exitflag_c08_KRG_1] = fmincon(func_c08_KRG_1,sp_c08,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x_opt)nonlconstage1(x_opt,W1_ls_08_llm,t_c08_KRG));
xopt_c08_KRG_1_denorm = denormalizeuv(xopt_c08_KRG_1,lbx_c08,ubx_c08); 
xopt_c08_KRG_1_org = [0 xopt_c08_KRG_1_denorm(:,1:3) 0 xopt_c08_KRG_1_denorm(:,4) 0];
[psf_c08_KRG_1,pf_c08_KRG_1,re_c08_KRG_1,beta_c08_KRG_1] = mcspsfconstraint8c(xopt_c08_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c08_KRG_1_mm = srgtsKRGEvaluate(t_c08_KRG,srgt_KRG_c08);

%Bounds in inactive subspace
W2_ls_08_llm = W_ls_08_llm(:,n_ls_08+1:m_c08);
t2_ls_08_KRG = xopt_c08_KRG_1*W2_ls_08_llm;

sp_c08_KRG = zeros(1,m_c08);
fmin_c08_KRG_1 = @(xmin_c08_1)(xmin_c08_1*W2_ls_08_llm(:,1));
[xmin_1_c08_KRG,fvalmin_1_c08_KRG,exitflagmin_1_c08_KRG]=fmincon(fmin_c08_KRG_1,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,2),W2_ls_08_llm(:,3),t2_ls_08_KRG(2),t2_ls_08_KRG(3)));

fmax_c08_KRG_1 = @(xmax_c08_1)(-(xmax_c08_1*W2_ls_08_llm(:,1)));
[xmax_1_c08_KRG,fvalmax_1_c08_KRG,exitflagmax_1_c08_KRG]=fmincon(fmax_c08_KRG_1,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,2),W2_ls_08_llm(:,3),t2_ls_08_KRG(2),t2_ls_08_KRG(3)));   

fmin_c08_KRG_2 = @(xmin_c08_2)(xmin_c08_2*W2_ls_08_llm(:,2));
[xmin_2_c08_KRG,fvalmin_2_c08_KRG,exitflagmin_2_c08_KRG]=fmincon(fmin_c08_KRG_2,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,3),t2_ls_08_KRG(1),t2_ls_08_KRG(3)));    

fmax_c08_KRG_2 = @(xmax_c08_2)(-(xmax_c08_2*W2_ls_08_llm(:,2)));
[xmax_2_c08_KRG,fvalmax_2_c08_KRG,exitflagmax_2_c08_KRG]=fmincon(fmax_c08_KRG_2,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,3),t2_ls_08_KRG(1),t2_ls_08_KRG(3)));    

fmin_c08_KRG_3 = @(xmin_c08_3)(xmin_c08_3*W2_ls_08_llm(:,3));
[xmin_3_c08_KRG,fvalmin_3_c08_KRG,exitflagmin_3_c08_KRG]=fmincon(fmin_c08_KRG_3,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,2),t2_ls_08_KRG(1),t2_ls_08_KRG(2)));       

fmax_c08_KRG_3 = @(xmax_c08_3)(-(xmax_c08_3*W2_ls_08_llm(:,3)));
[xmax_3_c08_KRG,fvalmax_3_c08_KRG,exitflagmax_3_c08_KRG]=fmincon(fmax_c08_KRG_3,sp_c08_KRG,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_KRG,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,2),t2_ls_08_KRG(1),t2_ls_08_KRG(2)));    

x_bnd_c08_KRG_1 = [xmin_1_c08_KRG;xmax_1_c08_KRG;xmin_2_c08_KRG;xmax_2_c08_KRG;xmin_3_c08_KRG;xmax_3_c08_KRG];

ind_c08_KRG = 2*(m_c08-n_ls_08);

x_bnd_c08_KRG_org = denormalizeuv([zeros(ind_c08_KRG,1) x_bnd_c08_KRG_1(:,1:3) zeros(ind_c08_KRG,1) x_bnd_c08_KRG_1(:,4) zeros(ind_c08_KRG,1)],lbx,ubx);
[psf_c08_KRG_2,pf_c08_KRG_2,re_c08_KRG_2,beta_c08_KRG_2] = mcspsfconstraint8c(x_bnd_c08_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c08_KRG = std([psf_c08_KRG_2;psf_c08_KRG_1]);
y_c08_KRG = constraint8c([x_bnd_c08_KRG_org repmat([mup3 mup4],ind_c08_KRG,1)]);

%% Constraint 9
%Finding t at which PSF is one
sp_t_c09_KRG = zeros(1,n_ls_09);
func_t_c09_KRG =  @(t_opt) (t_opt*zeros(n_ls_09,1));
[t_c09_KRG, fval_t_c09_KRG, exitflag_t_c09_KRG] = fmincon(func_t_c09_KRG,sp_t_c09_KRG,[],[],[],[],t1_ls_09_min ,t1_ls_09_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c09,@srgtsKRGEvaluate,1));

%Validation
func_c09_KRG_1 =  @(x_opt) (x_opt*zeros(m_c09,1));
[xopt_c09_KRG_1, fval_c09_KRG_1, exitflag_c09_KRG_1] = fmincon(func_c09_KRG_1,sp_c09,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x_opt)nonlconstage1(x_opt,W1_ls_09_llm,t_c09_KRG));
xopt_c09_KRG_1_denorm = denormalizeuv(xopt_c09_KRG_1,lbx_c09,ubx_c09); 
xopt_c09_KRG_1_org = [xopt_c09_KRG_1_denorm(:,1:4) 0 xopt_c09_KRG_1_denorm(:,5) 0];
[psf_c09_KRG_1,pf_c09_KRG_1,re_c09_KRG_1,beta_c09_KRG_1] = mcspsfconstraint9c(xopt_c09_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c09_KRG_1_mm = srgtsKRGEvaluate(t_c09_KRG,srgt_KRG_c09);

%Bounds in inactive subspace
W2_ls_09_llm = W_ls_09_llm(:,n_ls_09+1:m_c09);
% t2_ls_09_min = -diag(sign(W2_ls_09_llm)'*W2_ls_09_llm);
% t2_ls_09_max = diag(sign(W2_ls_09_llm)'*W2_ls_09_llm);

sp_c09_KRG = zeros(1,m_c09);
fmin_c09_KRG_1 = @(xmin_c09_1)(xmin_c09_1*W2_ls_09_llm(:,1));
[xmin_1_c09_KRG,fvalmin_1_c09_KRG,exitflagmin_1_c09_KRG]=fmincon(fmin_c09_KRG_1,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(2),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmax_c09_KRG_1 = @(xmax_c09_1)(-(xmax_c09_1*W2_ls_09_llm(:,1)));
[xmax_1_c09_KRG,fvalmax_1_c09_KRG,exitflagmax_1_c09_KRG]=fmincon(fmax_c09_KRG_1,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(2),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmin_c09_KRG_2 = @(xmin_c09_2)(xmin_c09_2*W2_ls_09_llm(:,2));
[xmin_2_c09_KRG,fvalmin_2_c09_KRG,exitflagmin_2_c09_KRG]=fmincon(fmin_c09_KRG_2,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmax_c09_KRG_2 = @(xmax_c09_2)(-(xmax_c09_2*W2_ls_09_llm(:,2)));
[xmax_2_c09_KRG,fvalmax_2_c09_KRG,exitflagmax_2_c09_KRG]=fmincon(fmax_c09_KRG_2,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmin_c09_KRG_3 = @(xmin_c09_3)(xmin_c09_3*W2_ls_09_llm(:,3));
[xmin_3_c09_KRG,fvalmin_3_c09_KRG,exitflagmin_3_c09_KRG]=fmincon(fmin_c09_KRG_3,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(4)));

fmax_c09_KRG_3 = @(xmax_c09_3)(-(xmax_c09_3*W2_ls_09_llm(:,3)));
[xmax_3_c09_KRG,fvalmax_3_c09_KRG,exitflagmax_3_c09_KRG]=fmincon(fmax_c09_KRG_3,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(4)));

fmin_c09_KRG_4 = @(xmin_c09_4)(xmin_c09_4*W2_ls_09_llm(:,4));
[xmin_4_c09_KRG,fvalmin_4_c09_KRG,exitflagmin_4_c09_KRG]=fmincon(fmin_c09_KRG_4,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(3)));

fmax_c09_KRG_4 = @(xmax_c09_4)(-(xmax_c09_4*W2_ls_09_llm(:,4)));
[xmax_4_c09_KRG,fvalmax_4_c09_KRG,exitflagmax_4_c09_KRG]=fmincon(fmax_c09_KRG_4,sp_c09_KRG,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_KRG,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(3)));

x_bnd_c09_KRG_1 = [xmin_1_c09_KRG;xmax_1_c09_KRG;xmin_2_c09_KRG;xmax_2_c09_KRG;xmin_3_c09_KRG;xmax_3_c09_KRG;xmin_4_c09_KRG;xmax_4_c09_KRG];

ind_c09_KRG = 2*(m_c09-n_ls_09);
x_bnd_c09_KRG_org = denormalizeuv([x_bnd_c09_KRG_1(:,1:4) zeros(ind_c09_KRG,1) x_bnd_c09_KRG_1(:,5) zeros(ind_c09_KRG,1)],lbx,ubx);
[psf_c09_KRG_2,pf_c09_KRG_2,re_c09_KRG_2,beta_c09_KRG_2] = mcspsfconstraint9c(x_bnd_c09_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

std_c09_KRG = std([psf_c09_KRG_2;psf_c09_KRG_1]);
y_c09_KRG = constraint9c([x_bnd_c09_KRG_org repmat([mup1 mup3],ind_c09_KRG,1)]);

%% Constraint 10
%Finding t at which PSF is one
sp_t_c10_KRG = zeros(1,n_ls_10);
func_t_c10_KRG =  @(t_opt) (t_opt*zeros(n_ls_10,1));
[t_c10_KRG, fval_t_c10_KRG, exitflag_t_c10_KRG] = fmincon(func_t_c10_KRG,sp_t_c10_KRG,[],[],[],[],t1_ls_10_min ,t1_ls_10_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_c10,@srgtsKRGEvaluate,1));

%Validation
func_c10_KRG_1 =  @(x_opt) (x_opt*zeros(m_c10,1));
[xopt_c10_KRG_1, fval_c10_KRG_1, exitflag_c10_KRG_1] = fmincon(func_c10_KRG_1,sp_c10,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x_opt)nonlconstage1(x_opt,W1_ls_10_llm,t_c10_KRG));
xopt_c10_KRG_1_denorm = denormalizeuv(xopt_c10_KRG_1,lbx_c10,ubx_c10);
xopt_c10_KRG_1_org = [0 0 xopt_c10_KRG_1_denorm(:,1) 0 xopt_c10_KRG_1_denorm(:,2:4)];
[psf_c10_KRG_1,pf_c10_KRG_1,re_c10_KRG_1,beta_c10_KRG_1] = mcspsfconstraint10c(xopt_c10_KRG_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c10_KRG_1_mm = srgtsKRGEvaluate(t_c10_KRG,srgt_KRG_c10);

%Bounds in inactive subspace
W2_ls_10_llm = W_ls_10_llm(:,n_ls_10+1:m_c10);
t2_ls_10_KRG = xopt_c10_KRG_1*W2_ls_10_llm;

sp_c10_KRG = zeros(1,m_c10);
fmin_c10_KRG_1 = @(xmin_c10_1)(xmin_c10_1*W2_ls_10_llm(:,1));
[xmin_1_c10_KRG,fvalmin_1_c10_KRG,exitflagmin_1_c10_KRG]=fmincon(fmin_c10_KRG_1,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,2),W2_ls_10_llm(:,3),t2_ls_10_KRG(2),t2_ls_10_KRG(3)));

fmax_c10_KRG_1 = @(xmax_c10_1)(-(xmax_c10_1*W2_ls_10_llm(:,1)));
[xmax_1_c10_KRG,fvalmax_1_c10_KRG,exitflagmax_1_c10_KRG]=fmincon(fmax_c10_KRG_1,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,2),W2_ls_10_llm(:,3),t2_ls_10_KRG(2),t2_ls_10_KRG(3)));   

fmin_c10_KRG_2 = @(xmin_c10_2)(xmin_c10_2*W2_ls_10_llm(:,2));
[xmin_2_c10_KRG,fvalmin_2_c10_KRG,exitflagmin_2_c10_KRG]=fmincon(fmin_c10_KRG_2,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,3),t2_ls_10_KRG(1),t2_ls_10_KRG(3)));    

fmax_c10_KRG_2 = @(xmax_c10_2)(-(xmax_c10_2*W2_ls_10_llm(:,2)));
[xmax_2_c10_KRG,fvalmax_2_c10_KRG,exitflagmax_2_c10_KRG]=fmincon(fmax_c10_KRG_2,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,3),t2_ls_10_KRG(1),t2_ls_10_KRG(3)));    

fmin_c10_KRG_3 = @(xmin_c10_3)(xmin_c10_3*W2_ls_10_llm(:,3));
[xmin_3_c10_KRG,fvalmin_3_c10_KRG,exitflagmin_3_c10_KRG]=fmincon(fmin_c10_KRG_3,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,2),t2_ls_10_KRG(1),t2_ls_10_KRG(2)));       

fmax_c10_KRG_3 = @(xmax_c10_3)(-(xmax_c10_3*W2_ls_10_llm(:,3)));
[xmax_3_c10_KRG,fvalmax_3_c10_KRG,exitflagmax_3_c10_KRG]=fmincon(fmax_c10_KRG_3,sp_c10_KRG,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_KRG,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,2),t2_ls_10_KRG(1),t2_ls_10_KRG(2)));    

x_bnd_c10_KRG_1 = [xmin_1_c10_KRG;xmax_1_c10_KRG;xmin_2_c10_KRG;xmax_2_c10_KRG;xmin_3_c10_KRG;xmax_3_c10_KRG];

ind_c10_KRG = 2*(m_c10-n_ls_10);
x_bnd_c10_KRG_org = denormalizeuv([zeros(ind_c10_KRG,2) x_bnd_c10_KRG_1(:,1) zeros(ind_c10_KRG,1) x_bnd_c10_KRG_1(:,2:4)],lbx,ubx);
[psf_c10_KRG_2,pf_c10_KRG_2,re_c10_KRG_2,beta_c10_KRG_2] = mcspsfconstraint10c(x_bnd_c10_KRG_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);

mean_c10_KRG = mean([psf_c10_KRG_2;psf_c10_KRG_1]);
std_c10_KRG = std([psf_c10_KRG_2;psf_c10_KRG_1]);
y_c10_KRG = constraint10c([x_bnd_c10_KRG_org(:,2:4) x_bnd_c10_KRG_org(:,6) repmat([mup2 mup3 mup4],ind_c10_KRG,1)]);

%% Weight

%Finding t at which PSF is one
sp_t_obj_KRG = zeros(1,n_obj);
func_t_obj_KRG =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_KRG, fval_t_obj_KRG, exitflag_t_obj_KRG] = fmincon(func_t_obj_KRG,sp_t_obj_KRG,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_KRG_obj_1,@srgtsKRGEvaluate,mean(y_obj_new)));

%Validation
func_obj_KRG_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_KRG_1, fval_obj_KRG_1, exitflag_obj_KRG_1] = fmincon(func_obj_KRG_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_KRG));
xopt_obj_KRG_1_denorm = denormalizeuv(xopt_obj_KRG_1,lbx_obj,ubx_obj); 
xopt_obj_KRG_1_org = [xopt_obj_KRG_1_denorm(:,1:5) xopt_obj_KRG_1_denorm(:,6)];
obj_KRG_1 = Weightc(xopt_obj_KRG_1_org);
obj_KRG_1_mm = srgtsKRGEvaluate(t_obj_KRG,srgt_KRG_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);
t2_obj_llm = xopt_obj_KRG_1*W2_obj_llm;

sp_obj_KRG = zeros(1,m_obj);
fmin_obj_KRG_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_KRG,fvalmin_1_obj_KRG,exitflagmin_1_obj_KRG]=fmincon(fmin_obj_KRG_1,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));     

fmax_obj_KRG_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_KRG,fvalmax_1_obj_KRG,exitflagmax_1_obj_KRG]=fmincon(fmax_obj_KRG_1,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));

fmin_obj_KRG_2 = @(xmin_obj_2)(xmin_obj_2*W2_obj_llm(:,2));
[xmin_2_obj_KRG,fvalmin_2_obj_KRG,exitflagmin_2_obj_KRG]=fmincon(fmin_obj_KRG_2,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));    

fmax_obj_KRG_2 = @(xmax_obj_2)(-(xmax_obj_2*W2_obj_llm(:,2)));
[xmax_2_obj_KRG,fvalmax_2_obj_KRG,exitflagmax_2_obj_KRG]=fmincon(fmax_obj_KRG_2,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));          

fmin_obj_KRG_3 = @(xmin_obj_3)(xmin_obj_3*W2_obj_llm(:,3));
[xmin_3_obj_KRG,fvalmin_3_obj_KRG,exitflagmin_3_obj_KRG]=fmincon(fmin_obj_KRG_3,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(4),t2_obj_llm(5))); 

fmax_obj_KRG_3 = @(xmax_obj_3)(-(xmax_obj_3*W2_obj_llm(:,3)));
[xmax_3_obj_KRG,fvalmax_3_obj_KRG,exitflagmax_3_obj_KRG]=fmincon(fmax_obj_KRG_3,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(4),t2_obj_llm(5)));

fmin_obj_KRG_4 = @(xmin_obj_4)(xmin_obj_4*W2_obj_llm(:,4));
[xmin_4_obj_KRG,fvalmin_4_obj_KRG,exitflagmin_4_obj_KRG]=fmincon(fmin_obj_KRG_4,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(5))); 

fmax_obj_KRG_4 = @(xmax_obj_4)(-(xmax_obj_4*W2_obj_llm(:,4)));
[xmax_4_obj_KRG,fvalmax_4_obj_KRG,exitflagmax_4_obj_KRG]=fmincon(fmax_obj_KRG_4,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(5)));    

fmin_obj_KRG_5 = @(xmin_obj_5)(xmin_obj_5*W2_obj_llm(:,5));
[xmin_5_obj_KRG,fvalmin_5_obj_KRG,exitflagmin_5_obj_KRG]=fmincon(fmin_obj_KRG_5,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4)));   

fmax_obj_KRG_5 = @(xmax_obj_5)(-(xmax_obj_5*W2_obj_llm(:,5)));
[xmax_5_obj_KRG,fvalmax_5_obj_KRG,exitflagmax_5_obj_KRG]=fmincon(fmax_obj_KRG_5,sp_obj_KRG,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_KRG,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4))); 


x_bnd_obj_KRG_1 = [xmin_1_obj_KRG;xmax_1_obj_KRG;xmin_2_obj_KRG;xmax_2_obj_KRG;xmin_3_obj_KRG;xmax_3_obj_KRG;xmin_4_obj_KRG;xmax_4_obj_KRG;xmin_5_obj_KRG;xmax_5_obj_KRG];

ind_obj_KRG = 2*(m_obj-n_obj);
x_bnd_obj_KRG_org = denormalizeuv([x_bnd_obj_KRG_1(:,1:5) x_bnd_obj_KRG_1(:,6)],lbx_obj,ubx_obj);
obj_KRG_2 = Weightc(x_bnd_obj_KRG_org);

std_obj_KRG = std([obj_KRG_2;obj_KRG_1]);
y_obj_KRG = Weightc(x_bnd_obj_KRG_org);

