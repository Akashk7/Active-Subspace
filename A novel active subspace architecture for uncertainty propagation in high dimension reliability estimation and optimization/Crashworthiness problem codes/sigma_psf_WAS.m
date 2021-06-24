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
xopt_c01_WAS_1_org = [0 xopt_c01_WAS_1_denorm(:,1:3) 0 xopt_c01_WAS_1_denorm(:,4) 0];
[psf_c01_WAS_1,pf_c01_WAS_1,re_c01_WAS_1,beta_c01_WAS_1] = mcspsfconstraint1c(xopt_c01_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c01_WAS_1_mm = srgtsWASEvaluate(t_c01_WAS,srgt_WAS_c01);

%Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);
t2_ls_01_WAS = xopt_c01_WAS_1*W2_ls_01_llm;

sp_c01_WAS = zeros(1,m_c01);
fmin_c01_WAS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,1));
[xmin_1_c01_WAS,fvalmin_1_c01_WAS,exitflagmin_1_c01_WAS]=fmincon(fmin_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,2),W2_ls_01_llm(:,3),t2_ls_01_WAS(2),t2_ls_01_WAS(3)));

fmax_c01_WAS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,1)));
[xmax_1_c01_WAS,fvalmax_1_c01_WAS,exitflagmax_1_c01_WAS]=fmincon(fmax_c01_WAS_1,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,2),W2_ls_01_llm(:,3),t2_ls_01_WAS(2),t2_ls_01_WAS(3)));   

fmin_c01_WAS_2 = @(xmin_c01_2)(xmin_c01_2*W2_ls_01_llm(:,2));
[xmin_2_c01_WAS,fvalmin_2_c01_WAS,exitflagmin_2_c01_WAS]=fmincon(fmin_c01_WAS_2,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,3),t2_ls_01_WAS(1),t2_ls_01_WAS(3)));    

fmax_c01_WAS_2 = @(xmax_c01_2)(-(xmax_c01_2*W2_ls_01_llm(:,2)));
[xmax_2_c01_WAS,fvalmax_2_c01_WAS,exitflagmax_2_c01_WAS]=fmincon(fmax_c01_WAS_2,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,3),t2_ls_01_WAS(1),t2_ls_01_WAS(3)));    

fmin_c01_WAS_3 = @(xmin_c01_3)(xmin_c01_3*W2_ls_01_llm(:,3));
[xmin_3_c01_WAS,fvalmin_3_c01_WAS,exitflagmin_3_c01_WAS]=fmincon(fmin_c01_WAS_3,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,2),t2_ls_01_WAS(1),t2_ls_01_WAS(2)));       

fmax_c01_WAS_3 = @(xmax_c01_3)(-(xmax_c01_3*W2_ls_01_llm(:,3)));
[xmax_3_c01_WAS,fvalmax_3_c01_WAS,exitflagmax_3_c01_WAS]=fmincon(fmax_c01_WAS_3,sp_c01_WAS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive_c01_1(x,t_c01_WAS,W1_ls_01_llm,W2_ls_01_llm(:,1),W2_ls_01_llm(:,2),t2_ls_01_WAS(1),t2_ls_01_WAS(2)));    

x_bnd_c01_WAS_1 = [xmin_1_c01_WAS;xmax_1_c01_WAS;xmin_2_c01_WAS;xmax_2_c01_WAS;xmin_3_c01_WAS;xmax_3_c01_WAS];

ind_c01_WAS = 2*(m_c01-n_ls_01);
x_bnd_c01_WAS_org = denormalizeuv([zeros(ind_c01_WAS,1) x_bnd_c01_WAS_1(:,1:3) zeros(ind_c01_WAS,1) x_bnd_c01_WAS_1(:,4) zeros(ind_c01_WAS,1)],lbx,ubx);
[psf_c01_WAS_2,pf_c01_WAS_2,re_c01_WAS_2,beta_c01_WAS_2] = mcspsfconstraint1c(x_bnd_c01_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c01_WAS,std_c01_WAS] = momentsusingyounetal([psf_c01_WAS_2;psf_c01_WAS_1]);
% mean_c01_WAS = mean([psf_c01_WAS_2;psf_c01_WAS_1]);
% std_c01_WAS = std([psf_c01_WAS_2;psf_c01_WAS_1]);
y_c01_WAS = constraint1c([x_bnd_c01_WAS_org(:,2:4) x_bnd_c01_WAS_org(:,6) repmat([mup2 mup3 mup4],ind_c01_WAS,1)]);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_WAS = zeros(1,n_ls_02);
func_t_c02_WAS =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_WAS, fval_t_c02_WAS, exitflag_t_c02_WAS] = fmincon(func_t_c02_WAS,sp_t_c02_WAS,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c02,@srgtsWASEvaluate,1));

%Validation
func_c02_WAS_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_WAS_1, fval_c02_WAS_1, exitflag_c02_WAS_1] = fmincon(func_c02_WAS_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_WAS));
xopt_c02_WAS_1_denorm = denormalizeuv(xopt_c02_WAS_1,lbx_c02,ubx_c02);
xopt_c02_WAS_1_org = [xopt_c02_WAS_1_denorm(:,1:3) zeros(1,4)];
[psf_c02_WAS_1,pf_c02_WAS_1,re_c02_WAS_1,beta_c02_WAS_1] = mcspsfconstraint2c(xopt_c02_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c02_WAS_1_mm = srgtsWASEvaluate(t_c02_WAS,srgt_WAS_c02);

%Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);
t2_ls_02_llm = xopt_c02_WAS_1*W2_ls_02_llm;

sp_c02_WAS = zeros(1,m_c02);
fmin_c02_WAS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,1));
[xmin_1_c02_WAS,fvalmin_1_c02_WAS,exitflagmin_1_c02_WAS]=fmincon(fmin_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm,W2_ls_02_llm(:,2),t2_ls_02_llm(2)));

fmax_c02_WAS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,1)));
[xmax_1_c02_WAS,fvalmax_1_c02_WAS,exitflagmax_1_c02_WAS]=fmincon(fmax_c02_WAS_1,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm,W2_ls_02_llm(:,2),t2_ls_02_llm(2)));

fmin_c02_WAS_2 = @(xmin_c02_2)(xmin_c02_2*W2_ls_02_llm(:,2));
[xmin_2_c02_WAS,fvalmin_2_c02_WAS,exitflagmin_2_c02_WAS]=fmincon(fmin_c02_WAS_2,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm,W2_ls_02_llm(:,1),t2_ls_02_llm(1)));

fmax_c02_WAS_2 = @(xmax_c02_2)(-(xmax_c02_2*W2_ls_02_llm(:,2)));
[xmax_2_c02_WAS,fvalmax_2_c02_WAS,exitflagmax_2_c02_WAS]=fmincon(fmax_c02_WAS_2,sp_c02_WAS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive_c02_1(x,t_c02_WAS,W1_ls_02_llm,W2_ls_02_llm(:,1),t2_ls_02_llm(1)));


x_bnd_c02_WAS_1 = [xmin_1_c02_WAS;xmax_1_c02_WAS;xmin_2_c02_WAS;xmax_2_c02_WAS];

ind_c02_WAS = 2*(m_c02-n_ls_02);
x_bnd_c02_WAS_org = denormalizeuv([x_bnd_c02_WAS_1(:,1:3) zeros(ind_c02_WAS,4)],lbx,ubx);
[psf_c02_WAS_2,pf_c02_WAS_2,re_c02_WAS_2,beta_c02_WAS_2] = mcspsfconstraint2c(x_bnd_c02_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c02_WAS,std_c02_WAS] = momentsusingyounetal([psf_c02_WAS_2;psf_c02_WAS_1]);

% mean_c02_WAS = mean([psf_c02_WAS_2;psf_c02_WAS_1]);
% std_c02_WAS = std([psf_c02_WAS_2;psf_c02_WAS_1]);
y_c02_WAS = constraint2c([x_bnd_c02_WAS_org(:,1:3) repmat([mup1 mup3],ind_c02_WAS,1)]);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_WAS = zeros(1,n_ls_03);
func_t_c03_WAS =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_WAS, fval_t_c03_WAS, exitflag_t_c03_WAS] = fmincon(func_t_c03_WAS,sp_t_c03_WAS,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c03,@srgtsWASEvaluate,1));

%Validation
func_c03_WAS_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_WAS_1, fval_c03_WAS_1, exitflag_c03_WAS_1] = fmincon(func_c03_WAS_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_WAS));
xopt_c03_WAS_1_denorm = denormalizeuv(xopt_c03_WAS_1,lbx_c03,ubx_c03); 
xopt_c03_WAS_1_org = [xopt_c03_WAS_1_denorm(:,1:3) 0 xopt_c03_WAS_1_denorm(:,4) 0 xopt_c03_WAS_1_denorm(:,5)];
[psf_c03_WAS_1,pf_c03_WAS_1,re_c03_WAS_1,beta_c03_WAS_1] = mcspsfconstraint3c(xopt_c03_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c03_WAS_1_mm = srgtsWASEvaluate(t_c03_WAS,srgt_WAS_c03);

%Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);
t2_ls_03_llm = xopt_c03_WAS_1*W2_ls_03_llm;

sp_c03_WAS = zeros(1,m_c03);
fmin_c03_WAS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,1));
[xmin_1_c03_WAS,fvalmin_1_c03_WAS,exitflagmin_1_c03_WAS]=fmincon(fmin_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(2),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmax_c03_WAS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,1)));
[xmax_1_c03_WAS,fvalmax_1_c03_WAS,exitflagmax_1_c03_WAS]=fmincon(fmax_c03_WAS_1,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(2),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmin_c03_WAS_2 = @(xmin_c03_2)(xmin_c03_2*W2_ls_03_llm(:,2));
[xmin_2_c03_WAS,fvalmin_2_c03_WAS,exitflagmin_2_c03_WAS]=fmincon(fmin_c03_WAS_2,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmax_c03_WAS_2 = @(xmax_c03_2)(-(xmax_c03_2*W2_ls_03_llm(:,2)));
[xmax_2_c03_WAS,fvalmax_2_c03_WAS,exitflagmax_2_c03_WAS]=fmincon(fmax_c03_WAS_2,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,3),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(3),t2_ls_03_llm(4)));

fmin_c03_WAS_3 = @(xmin_c03_3)(xmin_c03_3*W2_ls_03_llm(:,3));
[xmin_3_c03_WAS,fvalmin_3_c03_WAS,exitflagmin_3_c03_WAS]=fmincon(fmin_c03_WAS_3,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(4)));

fmax_c03_WAS_3 = @(xmax_c03_3)(-(xmax_c03_3*W2_ls_03_llm(:,3)));
[xmax_3_c03_WAS,fvalmax_3_c03_WAS,exitflagmax_3_c03_WAS]=fmincon(fmax_c03_WAS_3,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,4),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(4)));

fmin_c03_WAS_4 = @(xmin_c03_4)(xmin_c03_4*W2_ls_03_llm(:,4));
[xmin_4_c03_WAS,fvalmin_4_c03_WAS,exitflagmin_4_c03_WAS]=fmincon(fmin_c03_WAS_4,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(3)));

fmax_c03_WAS_4 = @(xmax_c03_4)(-(xmax_c03_4*W2_ls_03_llm(:,4)));
[xmax_4_c03_WAS,fvalmax_4_c03_WAS,exitflagmax_4_c03_WAS]=fmincon(fmax_c03_WAS_4,sp_c03_WAS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive_c03_4(x,t_c03_WAS,W1_ls_03_llm,W2_ls_03_llm(:,1),W2_ls_03_llm(:,2),W2_ls_03_llm(:,3),t2_ls_03_llm(1),t2_ls_03_llm(2),t2_ls_03_llm(3)));


x_bnd_c03_WAS_1 = [xmin_1_c03_WAS;xmax_1_c03_WAS;xmin_2_c03_WAS;xmax_2_c03_WAS;xmin_3_c03_WAS;xmax_3_c03_WAS;xmin_4_c03_WAS;xmax_4_c03_WAS];

ind_c03_WAS = 2*(m_c03-n_ls_03);
x_bnd_c03_WAS_org = denormalizeuv([x_bnd_c03_WAS_1(:,1:3) zeros(ind_c03_WAS,1) x_bnd_c03_WAS_1(:,4) zeros(ind_c03_WAS,1) x_bnd_c03_WAS_1(:,5)],lbx,ubx);
[psf_c03_WAS_2,pf_c03_WAS_2,re_c03_WAS_2,beta_c03_WAS_2] = mcspsfconstraint3c(x_bnd_c03_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c03_WAS,std_c03_WAS] = momentsusingyounetal([psf_c03_WAS_2;psf_c03_WAS_1]);

% std_c03_WAS = std([psf_c03_WAS_2;psf_c03_WAS_1]);
y_c03_WAS = constraint3c([x_bnd_c03_WAS_org repmat([mup1 mup2 mup3],ind_c03_WAS,1)]);

%% Constraint 4

%Finding t at which PSF is one
sp_t_c04_WAS = zeros(1,n_ls_04);
func_t_c04_WAS =  @(t_opt) (t_opt*zeros(n_ls_04,1));
[t_c04_WAS, fval_t_c04_WAS, exitflag_t_c04_WAS] = fmincon(func_t_c04_WAS,sp_t_c04_WAS,[],[],[],[],t1_ls_04_min ,t1_ls_04_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c04,@srgtsWASEvaluate,1));

%Validation
func_c04_WAS_1 =  @(x_opt) (x_opt*zeros(m_c04,1));
[xopt_c04_WAS_1, fval_c04_WAS_1, exitflag_c04_WAS_1] = fmincon(func_c04_WAS_1,sp_c04,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x_opt)nonlconstage1(x_opt,W1_ls_04_llm,t_c04_WAS));
xopt_c04_WAS_1_denorm = denormalizeuv(xopt_c04_WAS_1,lbx_c04,ubx_c04); 
xopt_c04_WAS_1_org = [xopt_c04_WAS_1_denorm(:,1:3) 0 xopt_c04_WAS_1_denorm(:,4:6)];
[psf_c04_WAS_1,pf_c04_WAS_1,re_c04_WAS_1,beta_c04_WAS_1] = mcspsfconstraint4c(xopt_c04_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c04_WAS_1_mm = srgtsWASEvaluate(t_c04_WAS,srgt_WAS_c04);

%Bounds in inactive subspace
W2_ls_04_llm = W_ls_04_llm(:,n_ls_04+1:m_c04);
t2_ls_04_llm = xopt_c04_WAS_1*W2_ls_04_llm;

sp_c04_WAS = zeros(1,m_c04);
fmin_c04_WAS_1 = @(xmin_c04_1)(xmin_c04_1*W2_ls_04_llm(:,1));
[xmin_1_c04_WAS,fvalmin_1_c04_WAS,exitflagmin_1_c04_WAS]=fmincon(fmin_c04_WAS_1,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_WAS_1 = @(xmax_c04_1)(-(xmax_c04_1*W2_ls_04_llm(:,1)));
[xmax_1_c04_WAS,fvalmax_1_c04_WAS,exitflagmax_1_c04_WAS]=fmincon(fmax_c04_WAS_1,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_WAS_2 = @(xmin_c04_2)(xmin_c04_2*W2_ls_04_llm(:,2));
[xmin_2_c04_WAS,fvalmin_2_c04_WAS,exitflagmin_2_c04_WAS]=fmincon(fmin_c04_WAS_2,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_WAS_2 = @(xmax_c04_2)(-(xmax_c04_2*W2_ls_04_llm(:,2)));
[xmax_2_c04_WAS,fvalmax_2_c04_WAS,exitflagmax_2_c04_WAS]=fmincon(fmax_c04_WAS_2,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(3),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_WAS_3 = @(xmin_c04_3)(xmin_c04_3*W2_ls_04_llm(:,3));
[xmin_3_c04_WAS,fvalmin_3_c04_WAS,exitflagmin_3_c04_WAS]=fmincon(fmin_c04_WAS_3,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmax_c04_WAS_3 = @(xmax_c04_3)(-(xmax_c04_3*W2_ls_04_llm(:,3)));
[xmax_3_c04_WAS,fvalmax_3_c04_WAS,exitflagmax_3_c04_WAS]=fmincon(fmax_c04_WAS_3,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,4),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(4),t2_ls_04_llm(5)));

fmin_c04_WAS_4 = @(xmin_c04_4)(xmin_c04_4*W2_ls_04_llm(:,4));
[xmin_4_c04_WAS,fvalmin_4_c04_WAS,exitflagmin_4_c04_WAS]=fmincon(fmin_c04_WAS_4,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(5)));

fmax_c04_WAS_4 = @(xmax_c04_4)(-(xmax_c04_4*W2_ls_04_llm(:,4)));
[xmax_4_c04_WAS,fvalmax_4_c04_WAS,exitflagmax_4_c04_WAS]=fmincon(fmax_c04_WAS_4,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,5),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(5)));

fmin_c04_WAS_5 = @(xmin_c04_5)(xmin_c04_5*W2_ls_04_llm(:,5));
[xmin_5_c04_WAS,fvalmin_5_c04_WAS,exitflagmin_5_c04_WAS]=fmincon(fmin_c04_WAS_5,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4)));

fmax_c04_WAS_5 = @(xmax_c04_5)(-(xmax_c04_5*W2_ls_04_llm(:,5)));
[xmax_5_c04_WAS,fvalmax_5_c04_WAS,exitflagmax_5_c04_WAS]=fmincon(fmax_c04_WAS_5,sp_c04_WAS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive_c04_3(x,t_c04_WAS,W1_ls_04_llm,W2_ls_04_llm(:,1),W2_ls_04_llm(:,2),W2_ls_04_llm(:,3),W2_ls_04_llm(:,4),t2_ls_04_llm(1),t2_ls_04_llm(2),t2_ls_04_llm(3),t2_ls_04_llm(4)));

x_bnd_c04_WAS_1 = [xmin_1_c04_WAS;xmax_1_c04_WAS;xmin_2_c04_WAS;xmax_2_c04_WAS;xmin_3_c04_WAS;xmax_3_c04_WAS;xmin_4_c04_WAS;xmax_4_c04_WAS;xmin_5_c04_WAS;xmax_5_c04_WAS];

ind_c04_WAS = 2*(m_c04-n_ls_04);
x_bnd_c04_WAS_org = denormalizeuv([x_bnd_c04_WAS_1(:,1:3) zeros(ind_c04_WAS,1) x_bnd_c04_WAS_1(:,4:6)],lbx,ubx);
[psf_c04_WAS_2,pf_c04_WAS_2,re_c04_WAS_2,beta_c04_WAS_2] = mcspsfconstraint4c(x_bnd_c04_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c04_WAS,std_c04_WAS] = momentsusingyounetal([psf_c04_WAS_2;psf_c04_WAS_1]);
% std_c04_WAS = std([psf_c04_WAS_2;psf_c04_WAS_1]);
y_c04_WAS = constraint4c([x_bnd_c04_WAS_org repmat([mup1 mup2 mup3],ind_c04_WAS,1)]);

%% Constraint 5
%Finding t at which PSF is one
sp_t_c05_WAS = zeros(1,n_ls_05);
func_t_c05_WAS =  @(t_opt) (t_opt*zeros(n_ls_05,1));
[t_c05_WAS, fval_t_c05_WAS, exitflag_t_c05_WAS] = fmincon(func_t_c05_WAS,sp_t_c05_WAS,[],[],[],[],t1_ls_05_min ,t1_ls_05_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c05,@srgtsWASEvaluate,1));

%Validation
func_c05_WAS_1 =  @(x_opt) (x_opt*zeros(m_c05,1));
[xopt_c05_WAS_1, fval_c05_WAS_1, exitflag_c05_WAS_1] = fmincon(func_c05_WAS_1,sp_c05,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x_opt)nonlconstage1(x_opt,W1_ls_05_llm,t_c05_WAS));
xopt_c05_WAS_1_denorm = denormalizeuv(xopt_c05_WAS_1,lbx_c05,ubx_c05); 
xopt_c05_WAS_1_org = [0 xopt_c05_WAS_1_denorm(:,1:2) 0 0 0 xopt_c05_WAS_1_denorm(:,3)];
[psf_c05_WAS_1,pf_c05_WAS_1,re_c05_WAS_1,beta_c05_WAS_1] = mcspsfconstraint5c(xopt_c05_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c05_WAS_1_mm = srgtsWASEvaluate(t_c05_WAS,srgt_WAS_c05);

%Bounds in inactive subspace
W2_ls_05_llm = W_ls_05_llm(:,n_ls_05+1:m_c05);
t2_ls_05_llm = xopt_c05_WAS_1*W2_ls_05_llm;

sp_c05_WAS = zeros(1,m_c05);
fmin_c05_WAS_1 = @(xmin_c05_1)(xmin_c05_1*W2_ls_05_llm(:,1));
[xmin_1_c05_WAS,fvalmin_1_c05_WAS,exitflagmin_1_c05_WAS]=fmincon(fmin_c05_WAS_1,sp_c05_WAS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_WAS,W1_ls_05_llm,W2_ls_05_llm(:,2),t2_ls_05_llm(:,2)));

fmax_c05_WAS_1 = @(xmax_c05_1)(-(xmax_c05_1*W2_ls_05_llm(:,1)));
[xmax_1_c05_WAS,fvalmax_1_c05_WAS,exitflagmax_1_c05_WAS]=fmincon(fmax_c05_WAS_1,sp_c05_WAS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_WAS,W1_ls_05_llm,W2_ls_05_llm(:,2),t2_ls_05_llm(:,2)));

fmin_c05_WAS_2 = @(xmin_c05_2)(xmin_c05_2*W2_ls_05_llm(:,2));
[xmin_2_c05_WAS,fvalmin_2_c05_WAS,exitflagmin_2_c05_WAS]=fmincon(fmin_c05_WAS_2,sp_c05_WAS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_WAS,W1_ls_05_llm,W2_ls_05_llm(:,1),t2_ls_05_llm(:,1)));

fmax_c05_WAS_2 = @(xmax_c05_2)(-(xmax_c05_2*W2_ls_05_llm(:,2)));
[xmax_2_c05_WAS,fvalmax_2_c05_WAS,exitflagmax_2_c05_WAS]=fmincon(fmax_c05_WAS_2,sp_c05_WAS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive_c05_2(x,t_c05_WAS,W1_ls_05_llm,W2_ls_05_llm(:,1),t2_ls_05_llm(:,1)));


x_bnd_c05_WAS_1 = [xmin_1_c05_WAS;xmax_1_c05_WAS;xmin_2_c05_WAS;xmax_2_c05_WAS];

ind_c05_WAS = 2*(m_c05-n_ls_05);

x_bnd_c05_WAS_org = denormalizeuv([zeros(ind_c05_WAS,1) x_bnd_c05_WAS_1(:,1:2) zeros(ind_c05_WAS,3) x_bnd_c05_WAS_1(:,3)],lbx,ubx);
[psf_c05_WAS_2,pf_c05_WAS_2,re_c05_WAS_2,beta_c05_WAS_2] = mcspsfconstraint5c(x_bnd_c05_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c05_WAS,std_c05_WAS] = momentsusingyounetal([psf_c05_WAS_2;psf_c05_WAS_1]);

% std_c05_WAS = std([psf_c05_WAS_2;psf_c05_WAS_1]);
y_c05_WAS = constraint5c([x_bnd_c05_WAS_org repmat([mup1 mup2 mup3],ind_c05_WAS,1)]);

%% Constraint 6
%Finding t at which PSF is one
sp_t_c06_WAS = zeros(1,n_ls_06);
func_t_c06_WAS =  @(t_opt) (t_opt*zeros(n_ls_06,1));
[t_c06_WAS, fval_t_c06_WAS, exitflag_t_c06_WAS] = fmincon(func_t_c06_WAS,sp_t_c06_WAS,[],[],[],[],t1_ls_06_min ,t1_ls_06_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c06,@srgtsWASEvaluate,mean(psf_c06)));

%Validation
func_c06_WAS_1 =  @(x_opt) (x_opt*zeros(m_c06,1));
[xopt_c06_WAS_1, fval_c06_WAS_1, exitflag_c06_WAS_1] = fmincon(func_c06_WAS_1,sp_c06,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x_opt)nonlconstage1(x_opt,W1_ls_06_llm,t_c06_WAS));
xopt_c06_WAS_1_denorm = denormalizeuv(xopt_c06_WAS_1,lbx_c06,ubx_c06); 
xopt_c06_WAS_1_org = [xopt_c06_WAS_1_denorm(:,1:3) 0 xopt_c06_WAS_1_denorm(:,4:6)];
[psf_c06_WAS_1,pf_c06_WAS_1,re_c06_WAS_1,beta_c06_WAS_1] = mcspsfconstraint6c(xopt_c06_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c06_WAS_1_mm = srgtsWASEvaluate(t_c06_WAS,srgt_WAS_c06);

%Bounds in inactive subspace
W2_ls_06_llm = W_ls_06_llm(:,n_ls_06+1:m_c06);
t2_ls_06_llm = xopt_c06_WAS_1*W2_ls_06_llm;

sp_c06_WAS = zeros(1,m_c06);
fmin_c06_WAS_1 = @(xmin_c06_1)(xmin_c06_1*W2_ls_06_llm(:,1));
[xmin_1_c06_WAS,fvalmin_1_c06_WAS,exitflagmin_1_c06_WAS]=fmincon(fmin_c06_WAS_1,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_WAS_1 = @(xmax_c06_1)(-(xmax_c06_1*W2_ls_06_llm(:,1)));
[xmax_1_c06_WAS,fvalmax_1_c06_WAS,exitflagmax_1_c06_WAS]=fmincon(fmax_c06_WAS_1,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_WAS_2 = @(xmin_c06_2)(xmin_c06_2*W2_ls_06_llm(:,2));
[xmin_2_c06_WAS,fvalmin_2_c06_WAS,exitflagmin_2_c06_WAS]=fmincon(fmin_c06_WAS_2,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_WAS_2 = @(xmax_c06_2)(-(xmax_c06_2*W2_ls_06_llm(:,2)));
[xmax_2_c06_WAS,fvalmax_2_c06_WAS,exitflagmax_2_c06_WAS]=fmincon(fmax_c06_WAS_2,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(3),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_WAS_3 = @(xmin_c06_3)(xmin_c06_3*W2_ls_06_llm(:,3));
[xmin_3_c06_WAS,fvalmin_3_c06_WAS,exitflagmin_3_c06_WAS]=fmincon(fmin_c06_WAS_3,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmax_c06_WAS_3 = @(xmax_c06_3)(-(xmax_c06_3*W2_ls_06_llm(:,3)));
[xmax_3_c06_WAS,fvalmax_3_c06_WAS,exitflagmax_3_c06_WAS]=fmincon(fmax_c06_WAS_3,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,4),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(4),t2_ls_06_llm(5)));

fmin_c06_WAS_4 = @(xmin_c06_4)(xmin_c06_4*W2_ls_06_llm(:,4));
[xmin_4_c06_WAS,fvalmin_4_c06_WAS,exitflagmin_4_c06_WAS]=fmincon(fmin_c06_WAS_4,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(5)));

fmax_c06_WAS_4 = @(xmax_c06_4)(-(xmax_c06_4*W2_ls_06_llm(:,4)));
[xmax_4_c06_WAS,fvalmax_4_c06_WAS,exitflagmax_4_c06_WAS]=fmincon(fmax_c06_WAS_4,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,5),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(5)));

fmin_c06_WAS_5 = @(xmin_c06_5)(xmin_c06_5*W2_ls_06_llm(:,5));
[xmin_5_c06_WAS,fvalmin_5_c06_WAS,exitflagmin_5_c06_WAS]=fmincon(fmin_c06_WAS_5,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4)));

fmax_c06_WAS_5 = @(xmax_c06_5)(-(xmax_c06_5*W2_ls_06_llm(:,5)));
[xmax_5_c06_WAS,fvalmax_5_c06_WAS,exitflagmax_5_c06_WAS]=fmincon(fmax_c06_WAS_5,sp_c06_WAS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive_c04_3(x,t_c06_WAS,W1_ls_06_llm,W2_ls_06_llm(:,1),W2_ls_06_llm(:,2),W2_ls_06_llm(:,3),W2_ls_06_llm(:,4),t2_ls_06_llm(1),t2_ls_06_llm(2),t2_ls_06_llm(3),t2_ls_06_llm(4)));

x_bnd_c06_WAS_1 = [xmin_1_c06_WAS;xmax_1_c06_WAS;xmin_2_c06_WAS;xmax_2_c06_WAS;xmin_3_c06_WAS;xmax_3_c06_WAS;xmin_4_c06_WAS;xmax_4_c06_WAS;xmin_5_c06_WAS;xmax_5_c06_WAS];

ind_c06_WAS = 2*(m_c06-n_ls_06);
x_bnd_c06_WAS_org = denormalizeuv([x_bnd_c06_WAS_1(:,1:3) zeros(ind_c06_WAS,1) x_bnd_c06_WAS_1(:,4:6)],lbx,ubx);
[psf_c06_WAS_2,pf_c06_WAS_2,re_c06_WAS_2,beta_c06_WAS_2] = mcspsfconstraint6c(x_bnd_c06_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c06_WAS,std_c06_WAS] = momentsusingyounetal([psf_c06_WAS_2;psf_c06_WAS_1]);

% std_c06_WAS = std([psf_c06_WAS_2;psf_c06_WAS_1]);
y_c06_WAS = constraint6c([x_bnd_c06_WAS_org repmat([mup1 mup2 mup3 mup4],ind_c06_WAS,1)]);

%% Constraint 7
%Finding t at which PSF is one
sp_t_c07_WAS = zeros(1,n_ls_07);
func_t_c07_WAS =  @(t_opt) (t_opt*zeros(n_ls_07,1));
[t_c07_WAS, fval_t_c07_WAS, exitflag_t_c07_WAS] = fmincon(func_t_c07_WAS,sp_t_c07_WAS,[],[],[],[],t1_ls_07_min ,t1_ls_07_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c07,@srgtsWASEvaluate,mean(psf_c07)));

%Validation
func_c07_WAS_1 =  @(x_opt) (x_opt*zeros(m_c07,1));
[xopt_c07_WAS_1, fval_c07_WAS_1, exitflag_c07_WAS_1] = fmincon(func_c07_WAS_1,sp_c07,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x_opt)nonlconstage1(x_opt,W1_ls_07_llm,t_c07_WAS));
xopt_c07_WAS_1_denorm = denormalizeuv(xopt_c07_WAS_1,lbx_c07,ubx_c07); 
xopt_c07_WAS_1_org = [xopt_c07_WAS_1_denorm(:,1:3) 0 xopt_c07_WAS_1_denorm(:,4:6)];
[psf_c07_WAS_1,pf_c07_WAS_1,re_c07_WAS_1,beta_c07_WAS_1] = mcspsfconstraint7c(xopt_c07_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c07_WAS_1_mm = srgtsWASEvaluate(t_c07_WAS,srgt_WAS_c07);

%Bounds in inactive subspace
W2_ls_07_llm = W_ls_07_llm(:,n_ls_07+1:m_c07);
t2_ls_07_llm = xopt_c07_WAS_1*W2_ls_07_llm;

sp_c07_WAS = zeros(1,m_c07);
fmin_c07_WAS_1 = @(xmin_c07_1)(xmin_c07_1*W2_ls_07_llm(:,1));
[xmin_1_c07_WAS,fvalmin_1_c07_WAS,exitflagmin_1_c07_WAS]=fmincon(fmin_c07_WAS_1,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_WAS_1 = @(xmax_c07_1)(-(xmax_c07_1*W2_ls_07_llm(:,1)));
[xmax_1_c07_WAS,fvalmax_1_c07_WAS,exitflagmax_1_c07_WAS]=fmincon(fmax_c07_WAS_1,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_WAS_2 = @(xmin_c07_2)(xmin_c07_2*W2_ls_07_llm(:,2));
[xmin_2_c07_WAS,fvalmin_2_c07_WAS,exitflagmin_2_c07_WAS]=fmincon(fmin_c07_WAS_2,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_WAS_2 = @(xmax_c07_2)(-(xmax_c07_2*W2_ls_07_llm(:,2)));
[xmax_2_c07_WAS,fvalmax_2_c07_WAS,exitflagmax_2_c07_WAS]=fmincon(fmax_c07_WAS_2,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(3),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_WAS_3 = @(xmin_c07_3)(xmin_c07_3*W2_ls_07_llm(:,3));
[xmin_3_c07_WAS,fvalmin_3_c07_WAS,exitflagmin_3_c07_WAS]=fmincon(fmin_c07_WAS_3,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmax_c07_WAS_3 = @(xmax_c07_3)(-(xmax_c07_3*W2_ls_07_llm(:,3)));
[xmax_3_c07_WAS,fvalmax_3_c07_WAS,exitflagmax_3_c07_WAS]=fmincon(fmax_c07_WAS_3,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,4),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(4),t2_ls_07_llm(5)));

fmin_c07_WAS_4 = @(xmin_c07_4)(xmin_c07_4*W2_ls_07_llm(:,4));
[xmin_4_c07_WAS,fvalmin_4_c07_WAS,exitflagmin_4_c07_WAS]=fmincon(fmin_c07_WAS_4,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(5)));

fmax_c07_WAS_4 = @(xmax_c07_4)(-(xmax_c07_4*W2_ls_07_llm(:,4)));
[xmax_4_c07_WAS,fvalmax_4_c07_WAS,exitflagmax_4_c07_WAS]=fmincon(fmax_c07_WAS_4,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,5),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(5)));

fmin_c07_WAS_5 = @(xmin_c07_5)(xmin_c07_5*W2_ls_07_llm(:,5));
[xmin_5_c07_WAS,fvalmin_5_c07_WAS,exitflagmin_5_c07_WAS]=fmincon(fmin_c07_WAS_5,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4)));

fmax_c07_WAS_5 = @(xmax_c07_5)(-(xmax_c07_5*W2_ls_07_llm(:,5)));
[xmax_5_c07_WAS,fvalmax_5_c07_WAS,exitflagmax_5_c07_WAS]=fmincon(fmax_c07_WAS_5,sp_c07_WAS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive_c04_3(x,t_c07_WAS,W1_ls_07_llm,W2_ls_07_llm(:,1),W2_ls_07_llm(:,2),W2_ls_07_llm(:,3),W2_ls_07_llm(:,4),t2_ls_07_llm(1),t2_ls_07_llm(2),t2_ls_07_llm(3),t2_ls_07_llm(4)));

x_bnd_c07_WAS_1 = [xmin_1_c07_WAS;xmax_1_c07_WAS;xmin_2_c07_WAS;xmax_2_c07_WAS;xmin_3_c07_WAS;xmax_3_c07_WAS;xmin_4_c07_WAS;xmax_4_c07_WAS;xmin_5_c07_WAS;xmax_5_c07_WAS];

ind_c07_WAS = 2*(m_c07-n_ls_07);
x_bnd_c07_WAS_org = denormalizeuv([x_bnd_c07_WAS_1(:,1:3) zeros(ind_c07_WAS,1) x_bnd_c07_WAS_1(:,4:6)],lbx,ubx);
[psf_c07_WAS_2,pf_c07_WAS_2,re_c07_WAS_2,beta_c07_WAS_2] = mcspsfconstraint7c(x_bnd_c07_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c07_WAS,std_c07_WAS] = momentsusingyounetal([psf_c07_WAS_2;psf_c07_WAS_1]);

% std_c07_WAS = std([psf_c07_WAS_2;psf_c07_WAS_1]);
y_c07_WAS = constraint7c([x_bnd_c07_WAS_org repmat([mup1 mup2 mup3 mup4],ind_c07_WAS,1)]);

%% Constraint 8
%Finding t at which PSF is one
sp_t_c08_WAS = zeros(1,n_ls_08);
func_t_c08_WAS =  @(t_opt) (t_opt*zeros(n_ls_08,1));
[t_c08_WAS, fval_t_c08_WAS, exitflag_t_c08_WAS] = fmincon(func_t_c08_WAS,sp_t_c08_WAS,[],[],[],[],t1_ls_08_min ,t1_ls_08_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c08,@srgtsWASEvaluate,0.9822));

%Validation
func_c08_WAS_1 =  @(x_opt) (x_opt*zeros(m_c08,1));
[xopt_c08_WAS_1, fval_c08_WAS_1, exitflag_c08_WAS_1] = fmincon(func_c08_WAS_1,sp_c08,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x_opt)nonlconstage1(x_opt,W1_ls_08_llm,t_c08_WAS));
xopt_c08_WAS_1_denorm = denormalizeuv(xopt_c08_WAS_1,lbx_c08,ubx_c08); 
xopt_c08_WAS_1_org = [0 xopt_c08_WAS_1_denorm(:,1:3) 0 xopt_c08_WAS_1_denorm(:,4) 0];
[psf_c08_WAS_1,pf_c08_WAS_1,re_c08_WAS_1,beta_c08_WAS_1] = mcspsfconstraint8c(xopt_c08_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c08_WAS_1_mm = srgtsWASEvaluate(t_c08_WAS,srgt_WAS_c08);

%Bounds in inactive subspace
W2_ls_08_llm = W_ls_08_llm(:,n_ls_08+1:m_c08);
t2_ls_08_WAS = xopt_c08_WAS_1*W2_ls_08_llm;

sp_c08_WAS = zeros(1,m_c08);
fmin_c08_WAS_1 = @(xmin_c08_1)(xmin_c08_1*W2_ls_08_llm(:,1));
[xmin_1_c08_WAS,fvalmin_1_c08_WAS,exitflagmin_1_c08_WAS]=fmincon(fmin_c08_WAS_1,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,2),W2_ls_08_llm(:,3),t2_ls_08_WAS(2),t2_ls_08_WAS(3)));

fmax_c08_WAS_1 = @(xmax_c08_1)(-(xmax_c08_1*W2_ls_08_llm(:,1)));
[xmax_1_c08_WAS,fvalmax_1_c08_WAS,exitflagmax_1_c08_WAS]=fmincon(fmax_c08_WAS_1,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,2),W2_ls_08_llm(:,3),t2_ls_08_WAS(2),t2_ls_08_WAS(3)));   

fmin_c08_WAS_2 = @(xmin_c08_2)(xmin_c08_2*W2_ls_08_llm(:,2));
[xmin_2_c08_WAS,fvalmin_2_c08_WAS,exitflagmin_2_c08_WAS]=fmincon(fmin_c08_WAS_2,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,3),t2_ls_08_WAS(1),t2_ls_08_WAS(3)));    

fmax_c08_WAS_2 = @(xmax_c08_2)(-(xmax_c08_2*W2_ls_08_llm(:,2)));
[xmax_2_c08_WAS,fvalmax_2_c08_WAS,exitflagmax_2_c08_WAS]=fmincon(fmax_c08_WAS_2,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,3),t2_ls_08_WAS(1),t2_ls_08_WAS(3)));    

fmin_c08_WAS_3 = @(xmin_c08_3)(xmin_c08_3*W2_ls_08_llm(:,3));
[xmin_3_c08_WAS,fvalmin_3_c08_WAS,exitflagmin_3_c08_WAS]=fmincon(fmin_c08_WAS_3,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,2),t2_ls_08_WAS(1),t2_ls_08_WAS(2)));       

fmax_c08_WAS_3 = @(xmax_c08_3)(-(xmax_c08_3*W2_ls_08_llm(:,3)));
[xmax_3_c08_WAS,fvalmax_3_c08_WAS,exitflagmax_3_c08_WAS]=fmincon(fmax_c08_WAS_3,sp_c08_WAS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive_c01_1(x,t_c08_WAS,W1_ls_08_llm,W2_ls_08_llm(:,1),W2_ls_08_llm(:,2),t2_ls_08_WAS(1),t2_ls_08_WAS(2)));    

x_bnd_c08_WAS_1 = [xmin_1_c08_WAS;xmax_1_c08_WAS;xmin_2_c08_WAS;xmax_2_c08_WAS;xmin_3_c08_WAS;xmax_3_c08_WAS];

ind_c08_WAS = 2*(m_c08-n_ls_08);

x_bnd_c08_WAS_org = denormalizeuv([zeros(ind_c08_WAS,1) x_bnd_c08_WAS_1(:,1:3) zeros(ind_c08_WAS,1) x_bnd_c08_WAS_1(:,4) zeros(ind_c08_WAS,1)],lbx,ubx);
[psf_c08_WAS_2,pf_c08_WAS_2,re_c08_WAS_2,beta_c08_WAS_2] = mcspsfconstraint8c(x_bnd_c08_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c08_WAS,std_c08_WAS] = momentsusingyounetal([psf_c08_WAS_2;psf_c08_WAS_1]);

% std_c08_WAS = std([psf_c08_WAS_2;psf_c08_WAS_1]);
y_c08_WAS = constraint8c([x_bnd_c08_WAS_org repmat([mup3 mup4],ind_c08_WAS,1)]);

%% Constraint 9
%Finding t at which PSF is one
sp_t_c09_WAS = zeros(1,n_ls_09);
func_t_c09_WAS =  @(t_opt) (t_opt*zeros(n_ls_09,1));
[t_c09_WAS, fval_t_c09_WAS, exitflag_t_c09_WAS] = fmincon(func_t_c09_WAS,sp_t_c09_WAS,[],[],[],[],t1_ls_09_min ,t1_ls_09_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c09,@srgtsWASEvaluate,1));

%Validation
func_c09_WAS_1 =  @(x_opt) (x_opt*zeros(m_c09,1));
[xopt_c09_WAS_1, fval_c09_WAS_1, exitflag_c09_WAS_1] = fmincon(func_c09_WAS_1,sp_c09,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x_opt)nonlconstage1(x_opt,W1_ls_09_llm,t_c09_WAS));
xopt_c09_WAS_1_denorm = denormalizeuv(xopt_c09_WAS_1,lbx_c09,ubx_c09); 
xopt_c09_WAS_1_org = [xopt_c09_WAS_1_denorm(:,1:4) 0 xopt_c09_WAS_1_denorm(:,5) 0];
[psf_c09_WAS_1,pf_c09_WAS_1,re_c09_WAS_1,beta_c09_WAS_1] = mcspsfconstraint9c(xopt_c09_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c09_WAS_1_mm = srgtsWASEvaluate(t_c09_WAS,srgt_WAS_c09);

%Bounds in inactive subspace
W2_ls_09_llm = W_ls_09_llm(:,n_ls_09+1:m_c09);
% t2_ls_09_min = -diag(sign(W2_ls_09_llm)'*W2_ls_09_llm);
% t2_ls_09_max = diag(sign(W2_ls_09_llm)'*W2_ls_09_llm);

sp_c09_WAS = zeros(1,m_c09);
fmin_c09_WAS_1 = @(xmin_c09_1)(xmin_c09_1*W2_ls_09_llm(:,1));
[xmin_1_c09_WAS,fvalmin_1_c09_WAS,exitflagmin_1_c09_WAS]=fmincon(fmin_c09_WAS_1,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(2),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmax_c09_WAS_1 = @(xmax_c09_1)(-(xmax_c09_1*W2_ls_09_llm(:,1)));
[xmax_1_c09_WAS,fvalmax_1_c09_WAS,exitflagmax_1_c09_WAS]=fmincon(fmax_c09_WAS_1,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(2),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmin_c09_WAS_2 = @(xmin_c09_2)(xmin_c09_2*W2_ls_09_llm(:,2));
[xmin_2_c09_WAS,fvalmin_2_c09_WAS,exitflagmin_2_c09_WAS]=fmincon(fmin_c09_WAS_2,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmax_c09_WAS_2 = @(xmax_c09_2)(-(xmax_c09_2*W2_ls_09_llm(:,2)));
[xmax_2_c09_WAS,fvalmax_2_c09_WAS,exitflagmax_2_c09_WAS]=fmincon(fmax_c09_WAS_2,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,3),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(3),t2_ls_09_llm(4)));

fmin_c09_WAS_3 = @(xmin_c09_3)(xmin_c09_3*W2_ls_09_llm(:,3));
[xmin_3_c09_WAS,fvalmin_3_c09_WAS,exitflagmin_3_c09_WAS]=fmincon(fmin_c09_WAS_3,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(4)));

fmax_c09_WAS_3 = @(xmax_c09_3)(-(xmax_c09_3*W2_ls_09_llm(:,3)));
[xmax_3_c09_WAS,fvalmax_3_c09_WAS,exitflagmax_3_c09_WAS]=fmincon(fmax_c09_WAS_3,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,4),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(4)));

fmin_c09_WAS_4 = @(xmin_c09_4)(xmin_c09_4*W2_ls_09_llm(:,4));
[xmin_4_c09_WAS,fvalmin_4_c09_WAS,exitflagmin_4_c09_WAS]=fmincon(fmin_c09_WAS_4,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(3)));

fmax_c09_WAS_4 = @(xmax_c09_4)(-(xmax_c09_4*W2_ls_09_llm(:,4)));
[xmax_4_c09_WAS,fvalmax_4_c09_WAS,exitflagmax_4_c09_WAS]=fmincon(fmax_c09_WAS_4,sp_c09_WAS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive_c03_4(x,t_c09_WAS,W1_ls_09_llm,W2_ls_09_llm(:,1),W2_ls_09_llm(:,2),W2_ls_09_llm(:,3),t2_ls_09_llm(1),t2_ls_09_llm(2),t2_ls_09_llm(3)));

x_bnd_c09_WAS_1 = [xmin_1_c09_WAS;xmax_1_c09_WAS;xmin_2_c09_WAS;xmax_2_c09_WAS;xmin_3_c09_WAS;xmax_3_c09_WAS;xmin_4_c09_WAS;xmax_4_c09_WAS];

ind_c09_WAS = 2*(m_c09-n_ls_09);
x_bnd_c09_WAS_org = denormalizeuv([x_bnd_c09_WAS_1(:,1:4) zeros(ind_c09_WAS,1) x_bnd_c09_WAS_1(:,5) zeros(ind_c09_WAS,1)],lbx,ubx);
[psf_c09_WAS_2,pf_c09_WAS_2,re_c09_WAS_2,beta_c09_WAS_2] = mcspsfconstraint9c(x_bnd_c09_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c09_WAS,std_c09_WAS] = momentsusingyounetal([psf_c09_WAS_2;psf_c09_WAS_1]);
% std_c09_WAS = std([psf_c09_WAS_2;psf_c09_WAS_1]);
y_c09_WAS = constraint9c([x_bnd_c09_WAS_org repmat([mup1 mup3],ind_c09_WAS,1)]);


%% Constraint 10
%Finding t at which PSF is one
sp_t_c10_WAS = zeros(1,n_ls_10);
func_t_c10_WAS =  @(t_opt) (t_opt*zeros(n_ls_10,1));
[t_c10_WAS, fval_t_c10_WAS, exitflag_t_c10_WAS] = fmincon(func_t_c10_WAS,sp_t_c10_WAS,[],[],[],[],t1_ls_10_min ,t1_ls_10_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_c10,@srgtsWASEvaluate,1));

%Validation
func_c10_WAS_1 =  @(x_opt) (x_opt*zeros(m_c10,1));
[xopt_c10_WAS_1, fval_c10_WAS_1, exitflag_c10_WAS_1] = fmincon(func_c10_WAS_1,sp_c10,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x_opt)nonlconstage1(x_opt,W1_ls_10_llm,t_c10_WAS));
xopt_c10_WAS_1_denorm = denormalizeuv(xopt_c10_WAS_1,lbx_c10,ubx_c10);
xopt_c10_WAS_1_org = [0 0 xopt_c10_WAS_1_denorm(:,1) 0 xopt_c10_WAS_1_denorm(:,2:4)];
[psf_c10_WAS_1,pf_c10_WAS_1,re_c10_WAS_1,beta_c10_WAS_1] = mcspsfconstraint10c(xopt_c10_WAS_1_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
psf_c10_WAS_1_mm = srgtsWASEvaluate(t_c10_WAS,srgt_WAS_c10);

%Bounds in inactive subspace
W2_ls_10_llm = W_ls_10_llm(:,n_ls_10+1:m_c10);
t2_ls_10_WAS = xopt_c10_WAS_1*W2_ls_10_llm;

sp_c10_WAS = zeros(1,m_c10);
fmin_c10_WAS_1 = @(xmin_c10_1)(xmin_c10_1*W2_ls_10_llm(:,1));
[xmin_1_c10_WAS,fvalmin_1_c10_WAS,exitflagmin_1_c10_WAS]=fmincon(fmin_c10_WAS_1,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,2),W2_ls_10_llm(:,3),t2_ls_10_WAS(2),t2_ls_10_WAS(3)));

fmax_c10_WAS_1 = @(xmax_c10_1)(-(xmax_c10_1*W2_ls_10_llm(:,1)));
[xmax_1_c10_WAS,fvalmax_1_c10_WAS,exitflagmax_1_c10_WAS]=fmincon(fmax_c10_WAS_1,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,2),W2_ls_10_llm(:,3),t2_ls_10_WAS(2),t2_ls_10_WAS(3)));   

fmin_c10_WAS_2 = @(xmin_c10_2)(xmin_c10_2*W2_ls_10_llm(:,2));
[xmin_2_c10_WAS,fvalmin_2_c10_WAS,exitflagmin_2_c10_WAS]=fmincon(fmin_c10_WAS_2,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,3),t2_ls_10_WAS(1),t2_ls_10_WAS(3)));    

fmax_c10_WAS_2 = @(xmax_c10_2)(-(xmax_c10_2*W2_ls_10_llm(:,2)));
[xmax_2_c10_WAS,fvalmax_2_c10_WAS,exitflagmax_2_c10_WAS]=fmincon(fmax_c10_WAS_2,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,3),t2_ls_10_WAS(1),t2_ls_10_WAS(3)));    

fmin_c10_WAS_3 = @(xmin_c10_3)(xmin_c10_3*W2_ls_10_llm(:,3));
[xmin_3_c10_WAS,fvalmin_3_c10_WAS,exitflagmin_3_c10_WAS]=fmincon(fmin_c10_WAS_3,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,2),t2_ls_10_WAS(1),t2_ls_10_WAS(2)));       

fmax_c10_WAS_3 = @(xmax_c10_3)(-(xmax_c10_3*W2_ls_10_llm(:,3)));
[xmax_3_c10_WAS,fvalmax_3_c10_WAS,exitflagmax_3_c10_WAS]=fmincon(fmax_c10_WAS_3,sp_c10_WAS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive_c01_1(x,t_c10_WAS,W1_ls_10_llm,W2_ls_10_llm(:,1),W2_ls_10_llm(:,2),t2_ls_10_WAS(1),t2_ls_10_WAS(2)));    

x_bnd_c10_WAS_1 = [xmin_1_c10_WAS;xmax_1_c10_WAS;xmin_2_c10_WAS;xmax_2_c10_WAS;xmin_3_c10_WAS;xmax_3_c10_WAS];

ind_c10_WAS = 2*(m_c10-n_ls_10);
x_bnd_c10_WAS_org = denormalizeuv([zeros(ind_c10_WAS,1) x_bnd_c10_WAS_1(:,1:3) zeros(ind_c10_WAS,1) x_bnd_c10_WAS_1(:,4) zeros(ind_c10_WAS,1)],lbx,ubx);
[psf_c10_WAS_2,pf_c10_WAS_2,re_c10_WAS_2,beta_c10_WAS_2] = mcspsfconstraint10c(x_bnd_c10_WAS_org,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[mean_c10_WAS,std_c10_WAS] = momentsusingyounetal([psf_c10_WAS_2;psf_c10_WAS_1]);

% mean_c10_WAS = mean([psf_c10_WAS_2;psf_c10_WAS_1]);
% std_c10_WAS = std([psf_c10_WAS_2;psf_c10_WAS_1]);
y_c10_WAS = constraint10c([x_bnd_c10_WAS_org(:,2:4) x_bnd_c10_WAS_org(:,6) repmat([mup2 mup3 mup4],ind_c10_WAS,1)]);

%% Weight

%Finding t at which PSF is one
sp_t_obj_WAS = zeros(1,n_obj);
func_t_obj_WAS =  @(t_opt) (t_opt*zeros(n_obj,1));
[t_obj_WAS, fval_t_obj_WAS, exitflag_t_obj_WAS] = fmincon(func_t_obj_WAS,sp_t_obj_WAS,[],[],[],[],t1_obj_min ,t1_obj_max ,@(t_opt)nonlcon_1(t_opt,srgt_WAS_obj_1,@srgtsWASEvaluate,mean(y_obj_new)));

%Validation
func_obj_WAS_1 =  @(x_opt) (x_opt*zeros(m_obj,1));
[xopt_obj_WAS_1, fval_obj_WAS_1, exitflag_obj_WAS_1] = fmincon(func_obj_WAS_1,sp_obj,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x_opt)nonlconstage1(x_opt,W1_obj_llm,t_obj_WAS));
xopt_obj_WAS_1_denorm = denormalizeuv(xopt_obj_WAS_1,lbx_obj,ubx_obj); 
xopt_obj_WAS_1_org = [xopt_obj_WAS_1_denorm(:,1:5) xopt_obj_WAS_1_denorm(:,6)];
obj_WAS_1 = Weightc(xopt_obj_WAS_1_org);
obj_WAS_1_mm = srgtsWASEvaluate(t_obj_WAS,srgt_WAS_obj_1);

%Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);
t2_obj_llm = xopt_obj_WAS_1*W2_obj_llm;

sp_obj_WAS = zeros(1,m_obj);
fmin_obj_WAS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,1));
[xmin_1_obj_WAS,fvalmin_1_obj_WAS,exitflagmin_1_obj_WAS]=fmincon(fmin_obj_WAS_1,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));     

fmax_obj_WAS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,1)));
[xmax_1_obj_WAS,fvalmax_1_obj_WAS,exitflagmax_1_obj_WAS]=fmincon(fmax_obj_WAS_1,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));

fmin_obj_WAS_2 = @(xmin_obj_2)(xmin_obj_2*W2_obj_llm(:,2));
[xmin_2_obj_WAS,fvalmin_2_obj_WAS,exitflagmin_2_obj_WAS]=fmincon(fmin_obj_WAS_2,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));    

fmax_obj_WAS_2 = @(xmax_obj_2)(-(xmax_obj_2*W2_obj_llm(:,2)));
[xmax_2_obj_WAS,fvalmax_2_obj_WAS,exitflagmax_2_obj_WAS]=fmincon(fmax_obj_WAS_2,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,3),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(3),t2_obj_llm(4),t2_obj_llm(5)));          

fmin_obj_WAS_3 = @(xmin_obj_3)(xmin_obj_3*W2_obj_llm(:,3));
[xmin_3_obj_WAS,fvalmin_3_obj_WAS,exitflagmin_3_obj_WAS]=fmincon(fmin_obj_WAS_3,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(4),t2_obj_llm(5))); 

fmax_obj_WAS_3 = @(xmax_obj_3)(-(xmax_obj_3*W2_obj_llm(:,3)));
[xmax_3_obj_WAS,fvalmax_3_obj_WAS,exitflagmax_3_obj_WAS]=fmincon(fmax_obj_WAS_3,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,4),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(4),t2_obj_llm(5)));

fmin_obj_WAS_4 = @(xmin_obj_4)(xmin_obj_4*W2_obj_llm(:,4));
[xmin_4_obj_WAS,fvalmin_4_obj_WAS,exitflagmin_4_obj_WAS]=fmincon(fmin_obj_WAS_4,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(5))); 

fmax_obj_WAS_4 = @(xmax_obj_4)(-(xmax_obj_4*W2_obj_llm(:,4)));
[xmax_4_obj_WAS,fvalmax_4_obj_WAS,exitflagmax_4_obj_WAS]=fmincon(fmax_obj_WAS_4,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,5),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(5)));    

fmin_obj_WAS_5 = @(xmin_obj_5)(xmin_obj_5*W2_obj_llm(:,5));
[xmin_5_obj_WAS,fvalmin_5_obj_WAS,exitflagmin_5_obj_WAS]=fmincon(fmin_obj_WAS_5,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4)));   

fmax_obj_WAS_5 = @(xmax_obj_5)(-(xmax_obj_5*W2_obj_llm(:,5)));
[xmax_5_obj_WAS,fvalmax_5_obj_WAS,exitflagmax_5_obj_WAS]=fmincon(fmax_obj_WAS_5,sp_obj_WAS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive_obj_1(x,t_obj_WAS,W1_obj_llm,W2_obj_llm(:,1),W2_obj_llm(:,2),W2_obj_llm(:,3),W2_obj_llm(:,4),t2_obj_llm(1),t2_obj_llm(2),t2_obj_llm(3),t2_obj_llm(4))); 


x_bnd_obj_WAS_1 = [xmin_1_obj_WAS;xmax_1_obj_WAS;xmin_2_obj_WAS;xmax_2_obj_WAS;xmin_3_obj_WAS;xmax_3_obj_WAS;xmin_4_obj_WAS;xmax_4_obj_WAS;xmin_5_obj_WAS;xmax_5_obj_WAS];

ind_obj_WAS = 2*(m_obj-n_obj);
x_bnd_obj_WAS_org = denormalizeuv([x_bnd_obj_WAS_1(:,1:5) x_bnd_obj_WAS_1(:,6)],lbx_obj,ubx_obj);
obj_WAS_2 = Weightc(x_bnd_obj_WAS_org);

std_obj_WAS = std([obj_WAS_2;obj_WAS_1]);
y_obj_WAS = Weightc(x_bnd_obj_WAS_org);