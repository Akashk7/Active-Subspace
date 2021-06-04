
%% Constraint 1

%Finding t at which PSF is one
sp_t_c01_PRS = zeros(1,n_ls_01);
func_t_c01_PRS =  @(t_opt) (t_opt*zeros(n_ls_01,1));
[t_c01_PRS, fval_t_c01_PRS, exitflag_t_c01_PRS] = fmincon(func_t_c01_PRS,sp_t_c01_PRS,[],[],[],[],t1_ls_01_min ,t1_ls_01_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c01,@srgtsPRSEvaluate,mm_psf_PRS_c01));

%Validation
func_c01_PRS_1 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xopt_c01_PRS_1, fval_c01_PRS_1, exitflag_c01_PRS_1] = fmincon(func_c01_PRS_1,sp_c01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t_c01_PRS));
xopt_c01_PRS_1_denorm = denormalizeuv(xopt_c01_PRS_1,lbx_c01,ubx_c01);
xopt_c01_PRS_1_org = xopt_c01_PRS_1_denorm;
[psf_c01_PRS_1,pf_c01_PRS_1,re_c01_PRS_1,beta_c01_PRS_1] = mcspsfconstraint1c(xopt_c01_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c01_PRS_1_mm = srgtsPRSEvaluate(t_c01_PRS,srgt_PRS_c01);

% Bounds in inactive subspace
W2_ls_01_llm = W_ls_01_llm(:,n_ls_01+1:m_c01);
t2_ls_01_PRS = xopt_c01_PRS_1*W2_ls_01_llm;

sp_c01_PRS = zeros(1,m_c01);
fmin_c01_PRS_1 = zeros(m_c01-n_ls_01,1);
fmax_c01_PRS_1 = zeros(m_c01-n_ls_01,1);
xmin_1_c01_PRS = zeros(m_c01-n_ls_01,m_c01);
xmax_1_c01_PRS = zeros(m_c01-n_ls_01,m_c01);
exitflagmin_1_c01_PRS = zeros(m_c01-n_ls_01,1);
exitflagmax_1_c01_PRS = zeros(m_c01-n_ls_01,1);

for i=1:m_c01-n_ls_01
    if i==1
        W2_nlcon_1_llm_c01 = W2_ls_01_llm(:,2:end);
        t2_nlcon_1_llm_c01 = t2_ls_01_PRS(:,2:end);

        fmin_c01_PRS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,i));
        [xmin_1_c01_PRS(i,:),fvalmin_1_c01_PRS,exitflagmin_1_c01_PRS(i,:)]=fmincon(fmin_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));

        fmax_c01_PRS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,i)));
        [xmax_1_c01_PRS(i,:),fvalmax_1_c01_PRS,exitflagmax_1_c01_PRS(i,:)]=fmincon(fmax_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));
                
    elseif i==m_c01-n_ls_01
        W2_nlcon_1_llm_c01 = W2_ls_01_llm(:,1:end-1);
        t2_nlcon_1_llm_c01 = t2_ls_01_PRS(:,1:end-1);

        fmin_c01_PRS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,i));
        [xmin_1_c01_PRS(i,:),fvalmin_1_c01_PRS,exitflagmin_1_c01_PRS(i,:)]=fmincon(fmin_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));

        fmax_c01_PRS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,i)));
        [xmax_1_c01_PRS(i,:),fvalmax_1_c01_PRS,exitflagmax_1_c01_PRS(i,:)]=fmincon(fmax_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));
        
    else 
        W2_nlcon_1_llm_c01 = [W2_ls_01_llm(:,1:i-1) W2_ls_01_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c01 = [t2_ls_01_PRS(:,1:i-1) t2_ls_01_PRS(:,i+1:end)];

        fmin_c01_PRS_1 = @(xmin_c01_1)(xmin_c01_1*W2_ls_01_llm(:,i));
        [xmin_1_c01_PRS(i,:),fvalmin_1_c01_PRS,exitflagmin_1_c01_PRS(i,:)]=fmincon(fmin_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));

        fmax_c01_PRS_1 = @(xmax_c01_1)(-(xmax_c01_1*W2_ls_01_llm(:,i)));
        [xmax_1_c01_PRS(i,:),fvalmax_1_c01_PRS,exitflagmax_1_c01_PRS(i,:)]=fmincon(fmax_c01_PRS_1,sp_c01_PRS,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x)nonlconinactive(x,t_c01_PRS,W1_ls_01_llm,W2_nlcon_1_llm_c01,t2_nlcon_1_llm_c01));
        
    end
end

x_bnd_c01_PRS_1 = [xmin_1_c01_PRS;xmax_1_c01_PRS];

ind_c01_PRS = 2*(m_c01-n_ls_01);
x_bnd_c01_PRS_org = denormalizeuv(x_bnd_c01_PRS_1,lbx,ubx);
[psf_c01_PRS_2,pf_c01_PRS_2,re_c01_PRS_2,beta_c01_PRS_2] = mcspsfconstraint1c(x_bnd_c01_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c01_PRS, std_c01_PRS] = momentsusingyounetal([psf_c01_PRS_2;psf_c01_PRS_1]);
% mean_c01_PRS = mean([psf_c01_PRS_2;psf_c01_PRS_1]);
% std_c01_PRS = std([psf_c01_PRS_2;psf_c01_PRS_1]);
y_c01_PRS = constraint1c(x_bnd_c01_PRS_org);

%% Constraint 2

%Finding t at which PSF is one
sp_t_c02_PRS = zeros(1,n_ls_02);
func_t_c02_PRS =  @(t_opt) (t_opt*zeros(n_ls_02,1));
[t_c02_PRS, fval_t_c02_PRS, exitflag_t_c02_PRS] = fmincon(func_t_c02_PRS,sp_t_c02_PRS,[],[],[],[],t1_ls_02_min ,t1_ls_02_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c02,@srgtsPRSEvaluate,mm_psf_PRS_c02));

%Validation
func_c02_PRS_1 =  @(x_opt) (x_opt*zeros(m_c02,1));
[xopt_c02_PRS_1, fval_c02_PRS_1, exitflag_c02_PRS_1] = fmincon(func_c02_PRS_1,sp_c02,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x_opt)nonlconstage1(x_opt,W1_ls_02_llm,t_c02_PRS));
xopt_c02_PRS_1_denorm = denormalizeuv(xopt_c02_PRS_1,lbx_c02,ubx_c02);
xopt_c02_PRS_1_org = xopt_c02_PRS_1_denorm;
[psf_c02_PRS_1,pf_c02_PRS_1,re_c02_PRS_1,beta_c02_PRS_1] = mcspsfconstraint2c(xopt_c02_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c02_PRS_1_mm = srgtsPRSEvaluate(t_c02_PRS,srgt_PRS_c02);

% Bounds in inactive subspace
W2_ls_02_llm = W_ls_02_llm(:,n_ls_02+1:m_c02);
t2_ls_02_PRS = xopt_c02_PRS_1*W2_ls_02_llm;

sp_c02_PRS = zeros(1,m_c02);
fmin_c02_PRS_1 = zeros(m_c02-n_ls_02,1);
fmax_c02_PRS_1 = zeros(m_c02-n_ls_02,1);
xmin_1_c02_PRS = zeros(m_c02-n_ls_02,m_c02);
xmax_1_c02_PRS = zeros(m_c02-n_ls_02,m_c02);
exitflagmin_1_c02_PRS = zeros(m_c02-n_ls_02,1);
exitflagmax_1_c02_PRS = zeros(m_c02-n_ls_02,1);

for i=1:m_c02-n_ls_02
    if i==1
        W2_nlcon_1_llm_c02 = W2_ls_02_llm(:,2:end);
        t2_nlcon_1_llm_c02 = t2_ls_02_PRS(:,2:end);

        fmin_c02_PRS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,i));
        [xmin_1_c02_PRS(i,:),fvalmin_1_c02_PRS,exitflagmin_1_c02_PRS(i,:)]=fmincon(fmin_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));

        fmax_c02_PRS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,i)));
        [xmax_1_c02_PRS(i,:),fvalmax_1_c02_PRS,exitflagmax_1_c02_PRS(i,:)]=fmincon(fmax_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));
                
    elseif i==m_c02-n_ls_02
        W2_nlcon_1_llm_c02 = W2_ls_02_llm(:,1:end-1);
        t2_nlcon_1_llm_c02 = t2_ls_02_PRS(:,1:end-1);

        fmin_c02_PRS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,i));
        [xmin_1_c02_PRS(i,:),fvalmin_1_c02_PRS,exitflagmin_1_c02_PRS(i,:)]=fmincon(fmin_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));

        fmax_c02_PRS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,i)));
        [xmax_1_c02_PRS(i,:),fvalmax_1_c02_PRS,exitflagmax_1_c02_PRS(i,:)]=fmincon(fmax_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));
    else 
        W2_nlcon_1_llm_c02 = [W2_ls_02_llm(:,1:i-1) W2_ls_02_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c02 = [t2_ls_02_PRS(:,1:i-1) t2_ls_02_PRS(:,i+1:end)];

        fmin_c02_PRS_1 = @(xmin_c02_1)(xmin_c02_1*W2_ls_02_llm(:,i));
        [xmin_1_c02_PRS(i,:),fvalmin_1_c02_PRS,exitflagmin_1_c02_PRS(i,:)]=fmincon(fmin_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));

        fmax_c02_PRS_1 = @(xmax_c02_1)(-(xmax_c02_1*W2_ls_02_llm(:,i)));
        [xmax_1_c02_PRS(i,:),fvalmax_1_c02_PRS,exitflagmax_1_c02_PRS(i,:)]=fmincon(fmax_c02_PRS_1,sp_c02_PRS,[],[],[],[],lbx_norm_c02,ubx_norm_c02,@(x)nonlconinactive(x,t_c02_PRS,W1_ls_02_llm,W2_nlcon_1_llm_c02,t2_nlcon_1_llm_c02));
    end
end

x_bnd_c02_PRS_1 = [xmin_1_c02_PRS;xmax_1_c02_PRS];

ind_c02_PRS = 2*(m_c02-n_ls_02);
x_bnd_c02_PRS_1_denorm = denormalizeuv(x_bnd_c02_PRS_1,lbx_c02,ubx_c02);
x_bnd_c02_PRS_org = x_bnd_c02_PRS_1_denorm;
[psf_c02_PRS_2,pf_c02_PRS_2,re_c02_PRS_2,beta_c02_PRS_2] = mcspsfconstraint2c(x_bnd_c02_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c02_PRS, std_c02_PRS] = momentsusingyounetal([psf_c02_PRS_2;psf_c02_PRS_1]);
y_c02_PRS = constraint2c(x_bnd_c02_PRS_org);

%% Constraint 3
%Finding t at which PSF is one
sp_t_c03_PRS = zeros(1,n_ls_03);
func_t_c03_PRS =  @(t_opt) (t_opt*zeros(n_ls_03,1));
[t_c03_PRS, fval_t_c03_PRS, exitflag_t_c03_PRS] = fmincon(func_t_c03_PRS,sp_t_c03_PRS,[],[],[],[],t1_ls_03_min ,t1_ls_03_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c03,@srgtsPRSEvaluate,mm_psf_PRS_c03));

%Validation
func_c03_PRS_1 =  @(x_opt) (x_opt*zeros(m_c03,1));
[xopt_c03_PRS_1, fval_c03_PRS_1, exitflag_c03_PRS_1] = fmincon(func_c03_PRS_1,sp_c03,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x_opt)nonlconstage1(x_opt,W1_ls_03_llm,t_c03_PRS));
xopt_c03_PRS_1_denorm = denormalizeuv(xopt_c03_PRS_1,lbx_c03,ubx_c03); 
xopt_c03_PRS_1_org = [xopt_c03_PRS_1_denorm(:,1:16) 0 xopt_c03_PRS_1_denorm(:,17)];
[psf_c03_PRS_1,pf_c03_PRS_1,re_c03_PRS_1,beta_c03_PRS_1] = mcspsfconstraint3c(xopt_c03_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c03_PRS_1_mm = srgtsPRSEvaluate(t_c03_PRS,srgt_PRS_c03);

% Bounds in inactive subspace
W2_ls_03_llm = W_ls_03_llm(:,n_ls_03+1:m_c03);
t2_ls_03_PRS = xopt_c03_PRS_1*W2_ls_03_llm;

sp_c03_PRS = zeros(1,m_c03);
fmin_c03_PRS_1 = zeros(m_c03-n_ls_03,1);
fmax_c03_PRS_1 = zeros(m_c03-n_ls_03,1);
xmin_1_c03_PRS = zeros(m_c03-n_ls_03,m_c03);
xmax_1_c03_PRS = zeros(m_c03-n_ls_03,m_c03);
exitflagmin_1_c03_PRS = zeros(m_c03-n_ls_03,1);
exitflagmax_1_c03_PRS = zeros(m_c03-n_ls_03,1);

for i=1:m_c03-n_ls_03
    if i==1
        W2_nlcon_1_llm_c03 = W2_ls_03_llm(:,2:end);
        t2_nlcon_1_llm_c03 = t2_ls_03_PRS(:,2:end);

        fmin_c03_PRS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,i));
        [xmin_1_c03_PRS(i,:),fvalmin_1_c03_PRS,exitflagmin_1_c03_PRS(i,:)]=fmincon(fmin_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));

        fmax_c03_PRS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,i)));
        [xmax_1_c03_PRS(i,:),fvalmax_1_c03_PRS,exitflagmax_1_c03_PRS(i,:)]=fmincon(fmax_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));
        
    elseif i==m_c03-n_ls_03
        W2_nlcon_1_llm_c03 = W2_ls_03_llm(:,1:end-1);
        t2_nlcon_1_llm_c03 = t2_ls_03_PRS(:,1:end-1);

        fmin_c03_PRS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,i));
        [xmin_1_c03_PRS(i,:),fvalmin_1_c03_PRS,exitflagmin_1_c03_PRS(i,:)]=fmincon(fmin_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));

        fmax_c03_PRS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,i)));
        [xmax_1_c03_PRS(i,:),fvalmax_1_c03_PRS,exitflagmax_1_c03_PRS(i,:)]=fmincon(fmax_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));

    else 
        W2_nlcon_1_llm_c03 = [W2_ls_03_llm(:,1:i-1) W2_ls_03_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c03 = [t2_ls_03_PRS(:,1:i-1) t2_ls_03_PRS(:,i+1:end)];

        fmin_c03_PRS_1 = @(xmin_c03_1)(xmin_c03_1*W2_ls_03_llm(:,i));
        [xmin_1_c03_PRS(i,:),fvalmin_1_c03_PRS,exitflagmin_1_c03_PRS(i,:)]=fmincon(fmin_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));

        fmax_c03_PRS_1 = @(xmax_c03_1)(-(xmax_c03_1*W2_ls_03_llm(:,i)));
        [xmax_1_c03_PRS(i,:),fvalmax_1_c03_PRS,exitflagmax_1_c03_PRS(i,:)]=fmincon(fmax_c03_PRS_1,sp_c03_PRS,[],[],[],[],lbx_norm_c03,ubx_norm_c03,@(x)nonlconinactive(x,t_c03_PRS,W1_ls_03_llm,W2_nlcon_1_llm_c03,t2_nlcon_1_llm_c03));
    end
end

x_bnd_c03_PRS_1 = [xmin_1_c03_PRS;xmax_1_c03_PRS];

ind_c03_PRS = 2*(m_c03-n_ls_03);
x_bnd_c03_PRS_1_denorm = denormalizeuv(x_bnd_c03_PRS_1,lbx_c03,ubx_c03);
x_bnd_c03_PRS_org = [x_bnd_c03_PRS_1_denorm(:,1:16) zeros(ind_c03_PRS,1) x_bnd_c03_PRS_1_denorm(:,17)];
[psf_c03_PRS_2,pf_c03_PRS_2,re_c03_PRS_2,beta_c03_PRS_2] = mcspsfconstraint3c(x_bnd_c03_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c03_PRS, std_c03_PRS] = momentsusingyounetal([psf_c03_PRS_2;psf_c03_PRS_1]);
% mean_c03_PRS = mean([psf_c03_PRS_2;psf_c03_PRS_1]);
% std_c03_PRS = std([psf_c03_PRS_2;psf_c03_PRS_1]);
y_c03_PRS = constraint3c(x_bnd_c03_PRS_org);

%% Constraint 4

%Finding t at which PSF is one
sp_t_c04_PRS = zeros(1,n_ls_04);
func_t_c04_PRS =  @(t_opt) (t_opt*zeros(n_ls_04,1));
[t_c04_PRS, fval_t_c04_PRS, exitflag_t_c04_PRS] = fmincon(func_t_c04_PRS,sp_t_c04_PRS,[],[],[],[],t1_ls_04_min ,t1_ls_04_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c04,@srgtsPRSEvaluate,mm_psf_PRS_c04));

%Validation
func_c04_PRS_1 =  @(x_opt) (x_opt*zeros(m_c04,1));
[xopt_c04_PRS_1, fval_c04_PRS_1, exitflag_c04_PRS_1] = fmincon(func_c04_PRS_1,sp_c04,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x_opt)nonlconstage1(x_opt,W1_ls_04_llm,t_c04_PRS));
xopt_c04_PRS_1_denorm = denormalizeuv(xopt_c04_PRS_1,lbx_c04,ubx_c04); 
xopt_c04_PRS_1_org = [xopt_c04_PRS_1_denorm(:,1:17) 0];
[psf_c04_PRS_1,pf_c04_PRS_1,re_c04_PRS_1,beta_c04_PRS_1] = mcspsfconstraint4c(xopt_c04_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c04_PRS_1_mm = srgtsPRSEvaluate(t_c04_PRS,srgt_PRS_c04);

% Bounds in inactive subspace
W2_ls_04_llm = W_ls_04_llm(:,n_ls_04+1:m_c04);
t2_ls_04_PRS = xopt_c04_PRS_1*W2_ls_04_llm;

sp_c04_PRS = zeros(1,m_c04);
fmin_c04_PRS_1 = zeros(m_c04-n_ls_04,1);
fmax_c04_PRS_1 = zeros(m_c04-n_ls_04,1);
xmin_1_c04_PRS = zeros(m_c04-n_ls_04,m_c04);
xmax_1_c04_PRS = zeros(m_c04-n_ls_04,m_c04);
exitflagmin_1_c04_PRS = zeros(m_c04-n_ls_04,1);
exitflagmax_1_c04_PRS = zeros(m_c04-n_ls_04,1);

for i=1:m_c04-n_ls_04
    if i==1
        W2_nlcon_1_llm_c04 = W2_ls_04_llm(:,2:end);
        t2_nlcon_1_llm_c04 = t2_ls_04_PRS(:,2:end);

        fmin_c04_PRS_1 = @(xmin_c04_1)(xmin_c04_1*W2_ls_04_llm(:,i));
        [xmin_1_c04_PRS(i,:),fvalmin_1_c04_PRS,exitflagmin_1_c04_PRS(i,:)]=fmincon(fmin_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));

        fmax_c04_PRS_1 = @(xmax_c04_1)(-(xmax_c04_1*W2_ls_04_llm(:,i)));
        [xmax_1_c04_PRS(i,:),fvalmax_1_c04_PRS,exitflagmax_1_c04_PRS(i,:)]=fmincon(fmax_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));   
    elseif i==m_c04-n_ls_04
        W2_nlcon_1_llm_c04 = W2_ls_04_llm(:,1:end-1);
        t2_nlcon_1_llm_c04 = t2_ls_04_PRS(:,1:end-1);

        fmin_c04_PRS_1 = @(xmin_c04_1)(xmin_c04_1*W2_ls_04_llm(:,i));
        [xmin_1_c04_PRS(i,:),fvalmin_1_c04_PRS,exitflagmin_1_c04_PRS(i,:)]=fmincon(fmin_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));

        fmax_c04_PRS_1 = @(xmax_c04_1)(-(xmax_c04_1*W2_ls_04_llm(:,i)));
        [xmax_1_c04_PRS(i,:),fvalmax_1_c04_PRS,exitflagmax_1_c04_PRS(i,:)]=fmincon(fmax_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));   
    else 
        W2_nlcon_1_llm_c04 = [W2_ls_04_llm(:,1:i-1) W2_ls_04_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c04 = [t2_ls_04_PRS(:,1:i-1) t2_ls_04_PRS(:,i+1:end)];

        fmin_c04_PRS_1 = @(xmin_c04_1)(xmin_c04_1*W2_ls_04_llm(:,i));
        [xmin_1_c04_PRS(i,:),fvalmin_1_c04_PRS,exitflagmin_1_c04_PRS(i,:)]=fmincon(fmin_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));

        fmax_c04_PRS_1 = @(xmax_c04_1)(-(xmax_c04_1*W2_ls_04_llm(:,i)));
        [xmax_1_c04_PRS(i,:),fvalmax_1_c04_PRS,exitflagmax_1_c04_PRS(i,:)]=fmincon(fmax_c04_PRS_1,sp_c04_PRS,[],[],[],[],lbx_norm_c04,ubx_norm_c04,@(x)nonlconinactive(x,t_c04_PRS,W1_ls_04_llm,W2_nlcon_1_llm_c04,t2_nlcon_1_llm_c04));   
    end
end

x_bnd_c04_PRS_1 = [xmin_1_c04_PRS;xmax_1_c04_PRS];

ind_c04_PRS = 2*(m_c04-n_ls_04);
x_bnd_c04_PRS_1_denorm = denormalizeuv(x_bnd_c04_PRS_1,lbx_c04,ubx_c04);
x_bnd_c04_PRS_org = [x_bnd_c04_PRS_1_denorm zeros(ind_c04_PRS,1)];
[psf_c04_PRS_2,pf_c04_PRS_2,re_c04_PRS_2,beta_c04_PRS_2] = mcspsfconstraint4c(x_bnd_c04_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c04_PRS, std_c04_PRS] = momentsusingyounetal([psf_c04_PRS_2;psf_c04_PRS_1]);
% mean_c04_PRS = mean([psf_c04_PRS_2;psf_c04_PRS_1]);
% std_c04_PRS = std([psf_c04_PRS_2;psf_c04_PRS_1]);
y_c04_PRS = constraint4c(x_bnd_c04_PRS_org);

%% Constraint 5
%Finding t at which PSF is one
sp_t_c05_PRS = zeros(1,n_ls_05);
func_t_c05_PRS =  @(t_opt) (t_opt*zeros(n_ls_05,1));
[t_c05_PRS, fval_t_c05_PRS, exitflag_t_c05_PRS] = fmincon(func_t_c05_PRS,sp_t_c05_PRS,[],[],[],[],t1_ls_05_min ,t1_ls_05_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c05,@srgtsPRSEvaluate,mm_psf_PRS_c05));

%Validation
func_c05_PRS_1 =  @(x_opt) (x_opt*zeros(m_c05,1));
[xopt_c05_PRS_1, fval_c05_PRS_1, exitflag_c05_PRS_1] = fmincon(func_c05_PRS_1,sp_c05,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x_opt)nonlconstage1(x_opt,W1_ls_05_llm,t_c05_PRS));
xopt_c05_PRS_1_denorm = denormalizeuv(xopt_c05_PRS_1,lbx_c05,ubx_c05); 
xopt_c05_PRS_1_org = xopt_c05_PRS_1_denorm;
[psf_c05_PRS_1,pf_c05_PRS_1,re_c05_PRS_1,beta_c05_PRS_1] = mcspsfconstraint5c(xopt_c05_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c05_PRS_1_mm = srgtsPRSEvaluate(t_c05_PRS,srgt_PRS_c05);

% Bounds in inactive subspace
W2_ls_05_llm = W_ls_05_llm(:,n_ls_05+1:m_c05);
t2_ls_05_PRS = xopt_c05_PRS_1*W2_ls_05_llm;

sp_c05_PRS = zeros(1,m_c05);
fmin_c05_PRS_1 = zeros(m_c05-n_ls_05,1);
fmax_c05_PRS_1 = zeros(m_c05-n_ls_05,1);
xmin_1_c05_PRS = zeros(m_c05-n_ls_05,m_c05);
xmax_1_c05_PRS = zeros(m_c05-n_ls_05,m_c05);
exitflagmin_1_c05_PRS = zeros(m_c05-n_ls_05,1);
exitflagmax_1_c05_PRS = zeros(m_c05-n_ls_05,1);

for i=1:m_c05-n_ls_05
    if i==1
        W2_nlcon_1_llm_c05 = W2_ls_05_llm(:,2:end);
        t2_nlcon_1_llm_c05 = t2_ls_05_PRS(:,2:end);

        fmin_c05_PRS_1 = @(xmin_c05_1)(xmin_c05_1*W2_ls_05_llm(:,i));
        [xmin_1_c05_PRS(i,:),fvalmin_1_c05_PRS,exitflagmin_1_c05_PRS(i,:)]=fmincon(fmin_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));

        fmax_c05_PRS_1 = @(xmax_c05_1)(-(xmax_c05_1*W2_ls_05_llm(:,i)));
        [xmax_1_c05_PRS(i,:),fvalmax_1_c05_PRS,exitflagmax_1_c05_PRS(i,:)]=fmincon(fmax_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));   
    elseif i==m_c05-n_ls_05
        W2_nlcon_1_llm_c05 = W2_ls_05_llm(:,1:end-1);
        t2_nlcon_1_llm_c05 = t2_ls_05_PRS(:,1:end-1);

        fmin_c05_PRS_1 = @(xmin_c05_1)(xmin_c05_1*W2_ls_05_llm(:,i));
        [xmin_1_c05_PRS(i,:),fvalmin_1_c05_PRS,exitflagmin_1_c05_PRS(i,:)]=fmincon(fmin_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));

        fmax_c05_PRS_1 = @(xmax_c05_1)(-(xmax_c05_1*W2_ls_05_llm(:,i)));
        [xmax_1_c05_PRS(i,:),fvalmax_1_c05_PRS,exitflagmax_1_c05_PRS(i,:)]=fmincon(fmax_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));   
    else 
        W2_nlcon_1_llm_c05 = [W2_ls_05_llm(:,1:i-1) W2_ls_05_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c05 = [t2_ls_05_PRS(:,1:i-1) t2_ls_05_PRS(:,i+1:end)];

        fmin_c05_PRS_1 = @(xmin_c05_1)(xmin_c05_1*W2_ls_05_llm(:,i));
        [xmin_1_c05_PRS(i,:),fvalmin_1_c05_PRS,exitflagmin_1_c05_PRS(i,:)]=fmincon(fmin_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));

        fmax_c05_PRS_1 = @(xmax_c05_1)(-(xmax_c05_1*W2_ls_05_llm(:,i)));
        [xmax_1_c05_PRS(i,:),fvalmax_1_c05_PRS,exitflagmax_1_c05_PRS(i,:)]=fmincon(fmax_c05_PRS_1,sp_c05_PRS,[],[],[],[],lbx_norm_c05,ubx_norm_c05,@(x)nonlconinactive(x,t_c05_PRS,W1_ls_05_llm,W2_nlcon_1_llm_c05,t2_nlcon_1_llm_c05));   
    end
end

x_bnd_c05_PRS_1 = [xmin_1_c05_PRS;xmax_1_c05_PRS];

ind_c05_PRS = 2*(m_c05-n_ls_05);
x_bnd_c05_PRS_1_denorm = denormalizeuv(x_bnd_c05_PRS_1,lbx_c05,ubx_c05);
x_bnd_c05_PRS_org = x_bnd_c05_PRS_1_denorm;
[psf_c05_PRS_2,pf_c05_PRS_2,re_c05_PRS_2,beta_c05_PRS_2] = mcspsfconstraint5c(x_bnd_c05_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c05_PRS, std_c05_PRS] = momentsusingyounetal([psf_c05_PRS_2;psf_c05_PRS_1]);
% mean_c05_PRS = mean([psf_c05_PRS_2;psf_c05_PRS_1]);
% std_c05_PRS = std([psf_c05_PRS_2;psf_c05_PRS_1]);
y_c05_PRS = constraint5c(x_bnd_c05_PRS_org);

%% Constraint 6
%Finding t at which PSF is one
sp_t_c06_PRS = zeros(1,n_ls_06);
func_t_c06_PRS =  @(t_opt) (t_opt*zeros(n_ls_06,1));
[t_c06_PRS, fval_t_c06_PRS, exitflag_t_c06_PRS] = fmincon(func_t_c06_PRS,sp_t_c06_PRS,[],[],[],[],t1_ls_06_min ,t1_ls_06_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c06,@srgtsPRSEvaluate,mm_psf_PRS_c06));

%Validation
func_c06_PRS_1 =  @(x_opt) (x_opt*zeros(m_c06,1));
[xopt_c06_PRS_1, fval_c06_PRS_1, exitflag_c06_PRS_1] = fmincon(func_c06_PRS_1,sp_c06,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x_opt)nonlconstage1(x_opt,W1_ls_06_llm,t_c06_PRS));
xopt_c06_PRS_1_denorm = denormalizeuv(xopt_c06_PRS_1,lbx_c06,ubx_c06); 
xopt_c06_PRS_1_org = xopt_c06_PRS_1_denorm;
[psf_c06_PRS_1,pf_c06_PRS_1,re_c06_PRS_1,beta_c06_PRS_1] = mcspsfconstraint6c(xopt_c06_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c06_PRS_1_mm = srgtsPRSEvaluate(t_c06_PRS,srgt_PRS_c06);

% Bounds in inactive subspace
W2_ls_06_llm = W_ls_06_llm(:,n_ls_06+1:m_c06);
t2_ls_06_PRS = xopt_c06_PRS_1*W2_ls_06_llm;

sp_c06_PRS = zeros(1,m_c06);
fmin_c06_PRS_1 = zeros(m_c06-n_ls_06,1);
fmax_c06_PRS_1 = zeros(m_c06-n_ls_06,1);
xmin_1_c06_PRS = zeros(m_c06-n_ls_06,m_c06);
xmax_1_c06_PRS = zeros(m_c06-n_ls_06,m_c06);
exitflagmin_1_c06_PRS = zeros(m_c06-n_ls_06,1);
exitflagmax_1_c06_PRS = zeros(m_c06-n_ls_06,1);

for i=1:m_c06-n_ls_06
    if i==1
        W2_nlcon_1_llm_c06 = W2_ls_06_llm(:,2:end);
        t2_nlcon_1_llm_c06 = t2_ls_06_PRS(:,2:end);

        fmin_c06_PRS_1 = @(xmin_c06_1)(xmin_c06_1*W2_ls_06_llm(:,i));
        [xmin_1_c06_PRS(i,:),fvalmin_1_c06_PRS,exitflagmin_1_c06_PRS(i,:)]=fmincon(fmin_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));

        fmax_c06_PRS_1 = @(xmax_c06_1)(-(xmax_c06_1*W2_ls_06_llm(:,i)));
        [xmax_1_c06_PRS(i,:),fvalmax_1_c06_PRS,exitflagmax_1_c06_PRS(i,:)]=fmincon(fmax_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));   
    elseif i==m_c06-n_ls_06
        W2_nlcon_1_llm_c06 = W2_ls_06_llm(:,1:end-1);
        t2_nlcon_1_llm_c06 = t2_ls_06_PRS(:,1:end-1);

        fmin_c06_PRS_1 = @(xmin_c06_1)(xmin_c06_1*W2_ls_06_llm(:,i));
        [xmin_1_c06_PRS(i,:),fvalmin_1_c06_PRS,exitflagmin_1_c06_PRS(i,:)]=fmincon(fmin_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));

        fmax_c06_PRS_1 = @(xmax_c06_1)(-(xmax_c06_1*W2_ls_06_llm(:,i)));
        [xmax_1_c06_PRS(i,:),fvalmax_1_c06_PRS,exitflagmax_1_c06_PRS(i,:)]=fmincon(fmax_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));   
    else 
        W2_nlcon_1_llm_c06 = [W2_ls_06_llm(:,1:i-1) W2_ls_06_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c06 = [t2_ls_06_PRS(:,1:i-1) t2_ls_06_PRS(:,i+1:end)];

        fmin_c06_PRS_1 = @(xmin_c06_1)(xmin_c06_1*W2_ls_06_llm(:,i));
        [xmin_1_c06_PRS(i,:),fvalmin_1_c06_PRS,exitflagmin_1_c06_PRS(i,:)]=fmincon(fmin_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));

        fmax_c06_PRS_1 = @(xmax_c06_1)(-(xmax_c06_1*W2_ls_06_llm(:,i)));
        [xmax_1_c06_PRS(i,:),fvalmax_1_c06_PRS,exitflagmax_1_c06_PRS(i,:)]=fmincon(fmax_c06_PRS_1,sp_c06_PRS,[],[],[],[],lbx_norm_c06,ubx_norm_c06,@(x)nonlconinactive(x,t_c06_PRS,W1_ls_06_llm,W2_nlcon_1_llm_c06,t2_nlcon_1_llm_c06));   
    end
end

x_bnd_c06_PRS_1 = [xmin_1_c06_PRS;xmax_1_c06_PRS];

ind_c06_PRS = 2*(m_c06-n_ls_06);
x_bnd_c06_PRS_1_denorm = denormalizeuv(x_bnd_c06_PRS_1,lbx_c06,ubx_c06);
x_bnd_c06_PRS_org = x_bnd_c06_PRS_1_denorm;
[psf_c06_PRS_2,pf_c06_PRS_2,re_c06_PRS_2,beta_c06_PRS_2] = mcspsfconstraint6c(x_bnd_c06_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c06_PRS, std_c06_PRS] = momentsusingyounetal([psf_c06_PRS_2;psf_c06_PRS_1]);

% mean_c06_PRS = mean([psf_c06_PRS_2;psf_c06_PRS_1]);
% std_c06_PRS = std([psf_c06_PRS_2;psf_c06_PRS_1]);
y_c06_PRS = constraint6c(x_bnd_c06_PRS_org);

%% Constraint 7
%Finding t at which PSF is one
sp_t_c07_PRS = zeros(1,n_ls_07);
func_t_c07_PRS =  @(t_opt) (t_opt*zeros(n_ls_07,1));
[t_c07_PRS, fval_t_c07_PRS, exitflag_t_c07_PRS] = fmincon(func_t_c07_PRS,sp_t_c07_PRS,[],[],[],[],t1_ls_07_min ,t1_ls_07_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c07,@srgtsPRSEvaluate,mm_psf_PRS_c07));

%Validation
func_c07_PRS_1 =  @(x_opt) (x_opt*zeros(m_c07,1));
[xopt_c07_PRS_1, fval_c07_PRS_1, exitflag_c07_PRS_1] = fmincon(func_c07_PRS_1,sp_c07,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x_opt)nonlconstage1(x_opt,W1_ls_07_llm,t_c07_PRS));
xopt_c07_PRS_1_denorm = denormalizeuv(xopt_c07_PRS_1,lbx_c07,ubx_c07); 
xopt_c07_PRS_1_org = xopt_c07_PRS_1_denorm;
[psf_c07_PRS_1,pf_c07_PRS_1,re_c07_PRS_1,beta_c07_PRS_1] = mcspsfconstraint7c(xopt_c07_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c07_PRS_1_mm = srgtsPRSEvaluate(t_c07_PRS,srgt_PRS_c07);

% Bounds in inactive subspace
W2_ls_07_llm = W_ls_07_llm(:,n_ls_07+1:m_c07);
t2_ls_07_PRS = xopt_c07_PRS_1*W2_ls_07_llm;

sp_c07_PRS = zeros(1,m_c07);
fmin_c07_PRS_1 = zeros(m_c07-n_ls_07,1);
fmax_c07_PRS_1 = zeros(m_c07-n_ls_07,1);
xmin_1_c07_PRS = zeros(m_c07-n_ls_07,m_c07);
xmax_1_c07_PRS = zeros(m_c07-n_ls_07,m_c07);
exitflagmin_1_c07_PRS = zeros(m_c07-n_ls_07,1);
exitflagmax_1_c07_PRS = zeros(m_c07-n_ls_07,1);

for i=1:m_c07-n_ls_07
    if i==1
        W2_nlcon_1_llm_c07 = W2_ls_07_llm(:,2:end);
        t2_nlcon_1_llm_c07 = t2_ls_07_PRS(:,2:end);

        fmin_c07_PRS_1 = @(xmin_c07_1)(xmin_c07_1*W2_ls_07_llm(:,i));
        [xmin_1_c07_PRS(i,:),fvalmin_1_c07_PRS,exitflagmin_1_c07_PRS(i,:)]=fmincon(fmin_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));

        fmax_c07_PRS_1 = @(xmax_c07_1)(-(xmax_c07_1*W2_ls_07_llm(:,i)));
        [xmax_1_c07_PRS(i,:),fvalmax_1_c07_PRS,exitflagmax_1_c07_PRS(i,:)]=fmincon(fmax_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));   
    elseif i==m_c07-n_ls_07
        W2_nlcon_1_llm_c07 = W2_ls_07_llm(:,1:end-1);
        t2_nlcon_1_llm_c07 = t2_ls_07_PRS(:,1:end-1);

        fmin_c07_PRS_1 = @(xmin_c07_1)(xmin_c07_1*W2_ls_07_llm(:,i));
        [xmin_1_c07_PRS(i,:),fvalmin_1_c07_PRS,exitflagmin_1_c07_PRS(i,:)]=fmincon(fmin_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));

        fmax_c07_PRS_1 = @(xmax_c07_1)(-(xmax_c07_1*W2_ls_07_llm(:,i)));
        [xmax_1_c07_PRS(i,:),fvalmax_1_c07_PRS,exitflagmax_1_c07_PRS(i,:)]=fmincon(fmax_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));   
    else 
        W2_nlcon_1_llm_c07 = [W2_ls_07_llm(:,1:i-1) W2_ls_07_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c07 = [t2_ls_07_PRS(:,1:i-1) t2_ls_07_PRS(:,i+1:end)];

        fmin_c07_PRS_1 = @(xmin_c07_1)(xmin_c07_1*W2_ls_07_llm(:,i));
        [xmin_1_c07_PRS(i,:),fvalmin_1_c07_PRS,exitflagmin_1_c07_PRS(i,:)]=fmincon(fmin_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));

        fmax_c07_PRS_1 = @(xmax_c07_1)(-(xmax_c07_1*W2_ls_07_llm(:,i)));
        [xmax_1_c07_PRS(i,:),fvalmax_1_c07_PRS,exitflagmax_1_c07_PRS(i,:)]=fmincon(fmax_c07_PRS_1,sp_c07_PRS,[],[],[],[],lbx_norm_c07,ubx_norm_c07,@(x)nonlconinactive(x,t_c07_PRS,W1_ls_07_llm,W2_nlcon_1_llm_c07,t2_nlcon_1_llm_c07));   
    end
end

x_bnd_c07_PRS_1 = [xmin_1_c07_PRS;xmax_1_c07_PRS];

ind_c07_PRS = 2*(m_c07-n_ls_07);
x_bnd_c07_PRS_1_denorm = denormalizeuv(x_bnd_c07_PRS_1,lbx_c07,ubx_c07);
x_bnd_c07_PRS_org = x_bnd_c07_PRS_1_denorm;
[psf_c07_PRS_2,pf_c07_PRS_2,re_c07_PRS_2,beta_c07_PRS_2] = mcspsfconstraint7c(x_bnd_c07_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c07_PRS, std_c07_PRS] = momentsusingyounetal([psf_c07_PRS_2;psf_c07_PRS_1]);
% mean_c07_PRS = mean([psf_c07_PRS_2;psf_c07_PRS_1]);
% std_c07_PRS = std([psf_c07_PRS_2;psf_c07_PRS_1]);
y_c07_PRS = constraint7c(x_bnd_c07_PRS_org);

%% Constraint 8
%Finding t at which PSF is one
sp_t_c08_PRS = zeros(1,n_ls_08);
func_t_c08_PRS =  @(t_opt) (t_opt*zeros(n_ls_08,1));
[t_c08_PRS, fval_t_c08_PRS, exitflag_t_c08_PRS] = fmincon(func_t_c08_PRS,sp_t_c08_PRS,[],[],[],[],t1_ls_08_min ,t1_ls_08_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c08,@srgtsPRSEvaluate,mm_psf_PRS_c08));

%Validation
func_c08_PRS_1 =  @(x_opt) (x_opt*zeros(m_c08,1));
[xopt_c08_PRS_1, fval_c08_PRS_1, exitflag_c08_PRS_1] = fmincon(func_c08_PRS_1,sp_c08,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x_opt)nonlconstage1(x_opt,W1_ls_08_llm,t_c08_PRS));
xopt_c08_PRS_1_denorm = denormalizeuv(xopt_c08_PRS_1,lbx_c08,ubx_c08); 
xopt_c08_PRS_1_org = xopt_c08_PRS_1_denorm;
[psf_c08_PRS_1,pf_c08_PRS_1,re_c08_PRS_1,beta_c08_PRS_1] = mcspsfconstraint8c(xopt_c08_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c08_PRS_1_mm = srgtsPRSEvaluate(t_c08_PRS,srgt_PRS_c08);

% Bounds in inactive subspace
W2_ls_08_llm = W_ls_08_llm(:,n_ls_08+1:m_c08);
t2_ls_08_PRS = xopt_c08_PRS_1*W2_ls_08_llm;

sp_c08_PRS = zeros(1,m_c08);
fmin_c08_PRS_1 = zeros(m_c08-n_ls_08,1);
fmax_c08_PRS_1 = zeros(m_c08-n_ls_08,1);
xmin_1_c08_PRS = zeros(m_c08-n_ls_08,m_c08);
xmax_1_c08_PRS = zeros(m_c08-n_ls_08,m_c08);
exitflagmin_1_c08_PRS = zeros(m_c08-n_ls_08,1);
exitflagmax_1_c08_PRS = zeros(m_c08-n_ls_08,1);

for i=1:m_c08-n_ls_08
    if i==1
        W2_nlcon_1_llm_c08 = W2_ls_08_llm(:,2:end);
        t2_nlcon_1_llm_c08 = t2_ls_08_PRS(:,2:end);

        fmin_c08_PRS_1 = @(xmin_c08_1)(xmin_c08_1*W2_ls_08_llm(:,i));
        [xmin_1_c08_PRS(i,:),fvalmin_1_c08_PRS,exitflagmin_1_c08_PRS(i,:)]=fmincon(fmin_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));

        fmax_c08_PRS_1 = @(xmax_c08_1)(-(xmax_c08_1*W2_ls_08_llm(:,i)));
        [xmax_1_c08_PRS(i,:),fvalmax_1_c08_PRS,exitflagmax_1_c08_PRS(i,:)]=fmincon(fmax_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));   
    elseif i==m_c08-n_ls_08
        W2_nlcon_1_llm_c08 = W2_ls_08_llm(:,1:end-1);
        t2_nlcon_1_llm_c08 = t2_ls_08_PRS(:,1:end-1);

        fmin_c08_PRS_1 = @(xmin_c08_1)(xmin_c08_1*W2_ls_08_llm(:,i));
        [xmin_1_c08_PRS(i,:),fvalmin_1_c08_PRS,exitflagmin_1_c08_PRS(i,:)]=fmincon(fmin_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));

        fmax_c08_PRS_1 = @(xmax_c08_1)(-(xmax_c08_1*W2_ls_08_llm(:,i)));
        [xmax_1_c08_PRS(i,:),fvalmax_1_c08_PRS,exitflagmax_1_c08_PRS(i,:)]=fmincon(fmax_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));   
    else 
        W2_nlcon_1_llm_c08 = [W2_ls_08_llm(:,1:i-1) W2_ls_08_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c08 = [t2_ls_08_PRS(:,1:i-1) t2_ls_08_PRS(:,i+1:end)];

        fmin_c08_PRS_1 = @(xmin_c08_1)(xmin_c08_1*W2_ls_08_llm(:,i));
        [xmin_1_c08_PRS(i,:),fvalmin_1_c08_PRS,exitflagmin_1_c08_PRS(i,:)]=fmincon(fmin_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));

        fmax_c08_PRS_1 = @(xmax_c08_1)(-(xmax_c08_1*W2_ls_08_llm(:,i)));
        [xmax_1_c08_PRS(i,:),fvalmax_1_c08_PRS,exitflagmax_1_c08_PRS(i,:)]=fmincon(fmax_c08_PRS_1,sp_c08_PRS,[],[],[],[],lbx_norm_c08,ubx_norm_c08,@(x)nonlconinactive(x,t_c08_PRS,W1_ls_08_llm,W2_nlcon_1_llm_c08,t2_nlcon_1_llm_c08));   
    end
end

x_bnd_c08_PRS_1 = [xmin_1_c08_PRS;xmax_1_c08_PRS];

ind_c08_PRS = 2*(m_c08-n_ls_08);
x_bnd_c08_PRS_1_denorm = denormalizeuv(x_bnd_c08_PRS_1,lbx_c08,ubx_c08);
x_bnd_c08_PRS_org = x_bnd_c08_PRS_1_denorm;
[psf_c08_PRS_2,pf_c08_PRS_2,re_c08_PRS_2,beta_c08_PRS_2] = mcspsfconstraint8c(x_bnd_c08_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c08_PRS, std_c08_PRS] = momentsusingyounetal([psf_c08_PRS_2;psf_c08_PRS_1]);
% mean_c08_PRS = mean([psf_c08_PRS_2;psf_c08_PRS_1]);
% std_c08_PRS = std([psf_c08_PRS_2;psf_c08_PRS_1]);
y_c08_PRS = constraint8c(x_bnd_c08_PRS_org);

%% Constraint 9

%Finding t at which PSF is one
sp_t_c09_PRS = zeros(1,n_ls_09);
func_t_c09_PRS =  @(t_opt) (t_opt*zeros(n_ls_09,1));
[t_c09_PRS, fval_t_c09_PRS, exitflag_t_c09_PRS] = fmincon(func_t_c09_PRS,sp_t_c09_PRS,[],[],[],[],t1_ls_09_min ,t1_ls_09_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c09,@srgtsPRSEvaluate,mm_psf_PRS_c09));

%Validation
func_c09_PRS_1 =  @(x_opt) (x_opt*zeros(m_c09,1));
[xopt_c09_PRS_1, fval_c09_PRS_1, exitflag_c09_PRS_1] = fmincon(func_c09_PRS_1,sp_c09,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x_opt)nonlconstage1(x_opt,W1_ls_09_llm,t_c09_PRS));
xopt_c09_PRS_1_denorm = denormalizeuv(xopt_c09_PRS_1,lbx_c09,ubx_c09); 
xopt_c09_PRS_1_org = xopt_c09_PRS_1_denorm;
[psf_c09_PRS_1,pf_c09_PRS_1,re_c09_PRS_1,beta_c09_PRS_1] = mcspsfconstraint9c(xopt_c09_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c09_PRS_1_mm = srgtsPRSEvaluate(t_c09_PRS,srgt_PRS_c09);

% Bounds in inactive subspace
W2_ls_09_llm = W_ls_09_llm(:,n_ls_09+1:m_c09);
t2_ls_09_PRS = xopt_c09_PRS_1*W2_ls_09_llm;

sp_c09_PRS = zeros(1,m_c09);
fmin_c09_PRS_1 = zeros(m_c09-n_ls_09,1);
fmax_c09_PRS_1 = zeros(m_c09-n_ls_09,1);
xmin_1_c09_PRS = zeros(m_c09-n_ls_09,m_c09);
xmax_1_c09_PRS = zeros(m_c09-n_ls_09,m_c09);
exitflagmin_1_c09_PRS = zeros(m_c09-n_ls_09,1);
exitflagmax_1_c09_PRS = zeros(m_c09-n_ls_09,1);

for i=1:m_c09-n_ls_09
    if i==1
        W2_nlcon_1_llm_c09 = W2_ls_09_llm(:,2:end);
        t2_nlcon_1_llm_c09 = t2_ls_09_PRS(:,2:end);

        fmin_c09_PRS_1 = @(xmin_c09_1)(xmin_c09_1*W2_ls_09_llm(:,i));
        [xmin_1_c09_PRS(i,:),fvalmin_1_c09_PRS,exitflagmin_1_c09_PRS(i,:)]=fmincon(fmin_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));

        fmax_c09_PRS_1 = @(xmax_c09_1)(-(xmax_c09_1*W2_ls_09_llm(:,i)));
        [xmax_1_c09_PRS(i,:),fvalmax_1_c09_PRS,exitflagmax_1_c09_PRS(i,:)]=fmincon(fmax_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));   
    elseif i==m_c09-n_ls_09
        W2_nlcon_1_llm_c09 = W2_ls_09_llm(:,1:end-1);
        t2_nlcon_1_llm_c09 = t2_ls_09_PRS(:,1:end-1);

        fmin_c09_PRS_1 = @(xmin_c09_1)(xmin_c09_1*W2_ls_09_llm(:,i));
        [xmin_1_c09_PRS(i,:),fvalmin_1_c09_PRS,exitflagmin_1_c09_PRS(i,:)]=fmincon(fmin_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));

        fmax_c09_PRS_1 = @(xmax_c09_1)(-(xmax_c09_1*W2_ls_09_llm(:,i)));
        [xmax_1_c09_PRS(i,:),fvalmax_1_c09_PRS,exitflagmax_1_c09_PRS(i,:)]=fmincon(fmax_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));   
    else 
        W2_nlcon_1_llm_c09 = [W2_ls_09_llm(:,1:i-1) W2_ls_09_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c09 = [t2_ls_09_PRS(:,1:i-1) t2_ls_09_PRS(:,i+1:end)];

        fmin_c09_PRS_1 = @(xmin_c09_1)(xmin_c09_1*W2_ls_09_llm(:,i));
        [xmin_1_c09_PRS(i,:),fvalmin_1_c09_PRS,exitflagmin_1_c09_PRS(i,:)]=fmincon(fmin_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));

        fmax_c09_PRS_1 = @(xmax_c09_1)(-(xmax_c09_1*W2_ls_09_llm(:,i)));
        [xmax_1_c09_PRS(i,:),fvalmax_1_c09_PRS,exitflagmax_1_c09_PRS(i,:)]=fmincon(fmax_c09_PRS_1,sp_c09_PRS,[],[],[],[],lbx_norm_c09,ubx_norm_c09,@(x)nonlconinactive(x,t_c09_PRS,W1_ls_09_llm,W2_nlcon_1_llm_c09,t2_nlcon_1_llm_c09));   
    end
end

x_bnd_c09_PRS_1 = [xmin_1_c09_PRS;xmax_1_c09_PRS];

ind_c09_PRS = 2*(m_c09-n_ls_09);
x_bnd_c09_PRS_1_denorm = denormalizeuv(x_bnd_c09_PRS_1,lbx_c09,ubx_c09);
x_bnd_c09_PRS_org = x_bnd_c09_PRS_1_denorm;
[psf_c09_PRS_2,pf_c09_PRS_2,re_c09_PRS_2,beta_c09_PRS_2] = mcspsfconstraint9c(x_bnd_c09_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c09_PRS, std_c09_PRS] = momentsusingyounetal([psf_c09_PRS_2;psf_c09_PRS_1]);
% mean_c09_PRS = mean([psf_c09_PRS_2;psf_c09_PRS_1]);
% std_c09_PRS = std([psf_c09_PRS_2;psf_c09_PRS_1]);
y_c09_PRS = constraint9c(x_bnd_c09_PRS_org);


%% Constraint 10
%Finding t at which PSF is one
sp_t_c10_PRS = zeros(1,n_ls_10);
func_t_c10_PRS =  @(t_opt) (t_opt*zeros(n_ls_10,1));
[t_c10_PRS, fval_t_c10_PRS, exitflag_t_c10_PRS] = fmincon(func_t_c10_PRS,sp_t_c10_PRS,[],[],[],[],t1_ls_10_min ,t1_ls_10_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c10,@srgtsPRSEvaluate,mm_psf_PRS_c10));

%Validation
func_c10_PRS_1 =  @(x_opt) (x_opt*zeros(m_c10,1));
[xopt_c10_PRS_1, fval_c10_PRS_1, exitflag_c10_PRS_1] = fmincon(func_c10_PRS_1,sp_c10,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x_opt)nonlconstage1(x_opt,W1_ls_10_llm,t_c10_PRS));
xopt_c10_PRS_1_denorm = denormalizeuv(xopt_c10_PRS_1,lbx_c10,ubx_c10);
xopt_c10_PRS_1_org = xopt_c10_PRS_1_denorm;
[psf_c10_PRS_1,pf_c10_PRS_1,re_c10_PRS_1,beta_c10_PRS_1] = mcspsfconstraint10c(xopt_c10_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c10_PRS_1_mm = srgtsPRSEvaluate(t_c10_PRS,srgt_PRS_c10);

% Bounds in inactive subspace
W2_ls_10_llm = W_ls_10_llm(:,n_ls_10+1:m_c10);
t2_ls_10_PRS = xopt_c10_PRS_1*W2_ls_10_llm;

sp_c10_PRS = zeros(1,m_c10);
fmin_c10_PRS_1 = zeros(m_c10-n_ls_10,1);
fmax_c10_PRS_1 = zeros(m_c10-n_ls_10,1);
xmin_1_c10_PRS = zeros(m_c10-n_ls_10,m_c10);
xmax_1_c10_PRS = zeros(m_c10-n_ls_10,m_c10);
exitflagmin_1_c10_PRS = zeros(m_c10-n_ls_10,1);
exitflagmax_1_c10_PRS = zeros(m_c10-n_ls_10,1);

for i=1:m_c10-n_ls_10
    if i==1
        W2_nlcon_1_llm_c10 = W2_ls_10_llm(:,2:end);
        t2_nlcon_1_llm_c10 = t2_ls_10_PRS(:,2:end);

        fmin_c10_PRS_1 = @(xmin_c10_1)(xmin_c10_1*W2_ls_10_llm(:,i));
        [xmin_1_c10_PRS(i,:),fvalmin_1_c10_PRS,exitflagmin_1_c10_PRS(i,:)]=fmincon(fmin_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));

        fmax_c10_PRS_1 = @(xmax_c10_1)(-(xmax_c10_1*W2_ls_10_llm(:,i)));
        [xmax_1_c10_PRS(i,:),fvalmax_1_c10_PRS,exitflagmax_1_c10_PRS(i,:)]=fmincon(fmax_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));   
    elseif i==m_c10-n_ls_10
        W2_nlcon_1_llm_c10 = W2_ls_10_llm(:,1:end-1);
        t2_nlcon_1_llm_c10 = t2_ls_10_PRS(:,1:end-1);

        fmin_c10_PRS_1 = @(xmin_c10_1)(xmin_c10_1*W2_ls_10_llm(:,i));
        [xmin_1_c10_PRS(i,:),fvalmin_1_c10_PRS,exitflagmin_1_c10_PRS(i,:)]=fmincon(fmin_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));

        fmax_c10_PRS_1 = @(xmax_c10_1)(-(xmax_c10_1*W2_ls_10_llm(:,i)));
        [xmax_1_c10_PRS(i,:),fvalmax_1_c10_PRS,exitflagmax_1_c10_PRS(i,:)]=fmincon(fmax_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));   
    else 
        W2_nlcon_1_llm_c10 = [W2_ls_10_llm(:,1:i-1) W2_ls_10_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c10 = [t2_ls_10_PRS(:,1:i-1) t2_ls_10_PRS(:,i+1:end)];

        fmin_c10_PRS_1 = @(xmin_c10_1)(xmin_c10_1*W2_ls_10_llm(:,i));
        [xmin_1_c10_PRS(i,:),fvalmin_1_c10_PRS,exitflagmin_1_c10_PRS(i,:)]=fmincon(fmin_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));

        fmax_c10_PRS_1 = @(xmax_c10_1)(-(xmax_c10_1*W2_ls_10_llm(:,i)));
        [xmax_1_c10_PRS(i,:),fvalmax_1_c10_PRS,exitflagmax_1_c10_PRS(i,:)]=fmincon(fmax_c10_PRS_1,sp_c10_PRS,[],[],[],[],lbx_norm_c10,ubx_norm_c10,@(x)nonlconinactive(x,t_c10_PRS,W1_ls_10_llm,W2_nlcon_1_llm_c10,t2_nlcon_1_llm_c10));   
    end
end

x_bnd_c10_PRS_1 = [xmin_1_c10_PRS;xmax_1_c10_PRS];

ind_c10_PRS = 2*(m_c10-n_ls_10);
x_bnd_c10_PRS_1_denorm = denormalizeuv(x_bnd_c10_PRS_1,lbx_c10,ubx_c10);
x_bnd_c10_PRS_org = x_bnd_c10_PRS_1_denorm;
[psf_c10_PRS_2,pf_c10_PRS_2,re_c10_PRS_2,beta_c10_PRS_2] = mcspsfconstraint10c(x_bnd_c10_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c10_PRS, std_c10_PRS] = momentsusingyounetal([psf_c10_PRS_2;psf_c10_PRS_1]);
% mean_c10_PRS = mean([psf_c10_PRS_2;psf_c10_PRS_1]);
% std_c10_PRS = std([psf_c10_PRS_2;psf_c10_PRS_1]);
y_c10_PRS = constraint10c(x_bnd_c10_PRS_org);

%% Constraint 11

%Finding t at which PSF is one
sp_t_c11_PRS = zeros(1,n_ls_11);
func_t_c11_PRS =  @(t_opt) (t_opt*zeros(n_ls_11,1));
[t_c11_PRS, fval_t_c11_PRS, exitflag_t_c11_PRS] = fmincon(func_t_c11_PRS,sp_t_c11_PRS,[],[],[],[],t1_ls_11_min ,t1_ls_11_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c11,@srgtsPRSEvaluate,mm_psf_PRS_c11));

%Validation
func_c11_PRS_1 =  @(x_opt) (x_opt*zeros(m_c11,1));
[xopt_c11_PRS_1, fval_c11_PRS_1, exitflag_c11_PRS_1] = fmincon(func_c11_PRS_1,sp_c11,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x_opt)nonlconstage1(x_opt,W1_ls_11_llm,t_c11_PRS));
xopt_c11_PRS_1_denorm = denormalizeuv(xopt_c11_PRS_1,lbx_c11,ubx_c11);
xopt_c11_PRS_1_org = [xopt_c11_PRS_1_denorm 0 0 0];
[psf_c11_PRS_1,pf_c11_PRS_1,re_c11_PRS_1,beta_c11_PRS_1] = mcspsfconstraint11c(xopt_c11_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c11_PRS_1_mm = srgtsPRSEvaluate(t_c11_PRS,srgt_PRS_c11);

% Bounds in inactive subspace
W2_ls_11_llm = W_ls_11_llm(:,n_ls_11+1:m_c11);
t2_ls_11_PRS = xopt_c11_PRS_1*W2_ls_11_llm;

sp_c11_PRS = zeros(1,m_c11);
fmin_c11_PRS_1 = zeros(m_c11-n_ls_11,1);
fmax_c11_PRS_1 = zeros(m_c11-n_ls_11,1);
xmin_1_c11_PRS = zeros(m_c11-n_ls_11,m_c11);
xmax_1_c11_PRS = zeros(m_c11-n_ls_11,m_c11);
exitflagmin_1_c11_PRS = zeros(m_c11-n_ls_11,1);
exitflagmax_1_c11_PRS = zeros(m_c11-n_ls_11,1);

for i=1:m_c11-n_ls_11
    if i==1
        W2_nlcon_1_llm_c11 = W2_ls_11_llm(:,2:end);
        t2_nlcon_1_llm_c11 = t2_ls_11_PRS(:,2:end);

        fmin_c11_PRS_1 = @(xmin_c11_1)(xmin_c11_1*W2_ls_11_llm(:,i));
        [xmin_1_c11_PRS(i,:),fvalmin_1_c11_PRS,exitflagmin_1_c11_PRS(i,:)]=fmincon(fmin_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));

        fmax_c11_PRS_1 = @(xmax_c11_1)(-(xmax_c11_1*W2_ls_11_llm(:,i)));
        [xmax_1_c11_PRS(i,:),fvalmax_1_c11_PRS,exitflagmax_1_c11_PRS(i,:)]=fmincon(fmax_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));   
    elseif i==m_c11-n_ls_11
        W2_nlcon_1_llm_c11 = W2_ls_11_llm(:,1:end-1);
        t2_nlcon_1_llm_c11 = t2_ls_11_PRS(:,1:end-1);

        fmin_c11_PRS_1 = @(xmin_c11_1)(xmin_c11_1*W2_ls_11_llm(:,i));
        [xmin_1_c11_PRS(i,:),fvalmin_1_c11_PRS,exitflagmin_1_c11_PRS(i,:)]=fmincon(fmin_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));

        fmax_c11_PRS_1 = @(xmax_c11_1)(-(xmax_c11_1*W2_ls_11_llm(:,i)));
        [xmax_1_c11_PRS(i,:),fvalmax_1_c11_PRS,exitflagmax_1_c11_PRS(i,:)]=fmincon(fmax_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));   
    else 
        W2_nlcon_1_llm_c11 = [W2_ls_11_llm(:,1:i-1) W2_ls_11_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c11 = [t2_ls_11_PRS(:,1:i-1) t2_ls_11_PRS(:,i+1:end)];

        fmin_c11_PRS_1 = @(xmin_c11_1)(xmin_c11_1*W2_ls_11_llm(:,i));
        [xmin_1_c11_PRS(i,:),fvalmin_1_c11_PRS,exitflagmin_1_c11_PRS(i,:)]=fmincon(fmin_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));

        fmax_c11_PRS_1 = @(xmax_c11_1)(-(xmax_c11_1*W2_ls_11_llm(:,i)));
        [xmax_1_c11_PRS(i,:),fvalmax_1_c11_PRS,exitflagmax_1_c11_PRS(i,:)]=fmincon(fmax_c11_PRS_1,sp_c11_PRS,[],[],[],[],lbx_norm_c11,ubx_norm_c11,@(x)nonlconinactive(x,t_c11_PRS,W1_ls_11_llm,W2_nlcon_1_llm_c11,t2_nlcon_1_llm_c11));   
    end
end

x_bnd_c11_PRS_1 = [xmin_1_c11_PRS;xmax_1_c11_PRS];

ind_c11_PRS = 2*(m_c11-n_ls_11);
x_bnd_c11_PRS_1_denorm = denormalizeuv(x_bnd_c11_PRS_1,lbx_c11,ubx_c11);
x_bnd_c11_PRS_org = [x_bnd_c11_PRS_1_denorm zeros(ind_c11_PRS,3)];
[psf_c11_PRS_2,pf_c11_PRS_2,re_c11_PRS_2,beta_c11_PRS_2] = mcspsfconstraint11c(x_bnd_c11_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c11_PRS, std_c11_PRS] = momentsusingyounetal([psf_c11_PRS_2;psf_c11_PRS_1]);

% mean_c11_PRS = mean([psf_c11_PRS_2;psf_c11_PRS_1]);
% std_c11_PRS = std([psf_c11_PRS_2;psf_c11_PRS_1]);
y_c11_PRS = constraint11c(x_bnd_c11_PRS_org);

%% Constraint 12

%Finding t at which PSF is one
sp_t_c12_PRS = zeros(1,n_ls_12);
func_t_c12_PRS =  @(t_opt) (t_opt*zeros(n_ls_12,1));
[t_c12_PRS, fval_t_c12_PRS, exitflag_t_c12_PRS] = fmincon(func_t_c12_PRS,sp_t_c12_PRS,[],[],[],[],t1_ls_12_min ,t1_ls_12_max ,@(t_opt)nonlcon_1(t_opt,srgt_PRS_c12,@srgtsPRSEvaluate,mm_psf_PRS_c12));

%Validation
func_c12_PRS_1 =  @(x_opt) (x_opt*zeros(m_c12,1));
[xopt_c12_PRS_1, fval_c12_PRS_1, exitflag_c12_PRS_1] = fmincon(func_c12_PRS_1,sp_c12,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x_opt)nonlconstage1(x_opt,W1_ls_12_llm,t_c12_PRS));
xopt_c12_PRS_1_denorm = denormalizeuv(xopt_c12_PRS_1,lbx_c12,ubx_c12);
xopt_c12_PRS_1_org = [xopt_c12_PRS_1_denorm 0 0 0];
[psf_c12_PRS_1,pf_c12_PRS_1,re_c12_PRS_1,beta_c12_PRS_1] = mcspsfconstraint12c(xopt_c12_PRS_1_org,lbx,ubx,sd,Nmcs,Pftarget);
psf_c12_PRS_1_mm = srgtsPRSEvaluate(t_c12_PRS,srgt_PRS_c12);

% Bounds in inactive subspace
W2_ls_12_llm = W_ls_12_llm(:,n_ls_12+1:m_c12);
t2_ls_12_PRS = xopt_c12_PRS_1*W2_ls_12_llm;

sp_c12_PRS = zeros(1,m_c12);
fmin_c12_PRS_1 = zeros(m_c12-n_ls_12,1);
fmax_c12_PRS_1 = zeros(m_c12-n_ls_12,1);
xmin_1_c12_PRS = zeros(m_c12-n_ls_12,m_c12);
xmax_1_c12_PRS = zeros(m_c12-n_ls_12,m_c12);
exitflagmin_1_c12_PRS = zeros(m_c12-n_ls_12,1);
exitflagmax_1_c12_PRS = zeros(m_c12-n_ls_12,1);

for i=1:m_c12-n_ls_12
    if i==1
        W2_nlcon_1_llm_c12 = W2_ls_12_llm(:,2:end);
        t2_nlcon_1_llm_c12 = t2_ls_12_PRS(:,2:end);

        fmin_c12_PRS_1 = @(xmin_c12_1)(xmin_c12_1*W2_ls_12_llm(:,i));
        [xmin_1_c12_PRS(i,:),fvalmin_1_c12_PRS,exitflagmin_1_c12_PRS(i,:)]=fmincon(fmin_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));

        fmax_c12_PRS_1 = @(xmax_c12_1)(-(xmax_c12_1*W2_ls_12_llm(:,i)));
        [xmax_1_c12_PRS(i,:),fvalmax_1_c12_PRS,exitflagmax_1_c12_PRS(i,:)]=fmincon(fmax_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));   
    elseif i==m_c12-n_ls_12
        W2_nlcon_1_llm_c12 = W2_ls_12_llm(:,1:end-1);
        t2_nlcon_1_llm_c12 = t2_ls_12_PRS(:,1:end-1);

        fmin_c12_PRS_1 = @(xmin_c12_1)(xmin_c12_1*W2_ls_12_llm(:,i));
        [xmin_1_c12_PRS(i,:),fvalmin_1_c12_PRS,exitflagmin_1_c12_PRS(i,:)]=fmincon(fmin_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));

        fmax_c12_PRS_1 = @(xmax_c12_1)(-(xmax_c12_1*W2_ls_12_llm(:,i)));
        [xmax_1_c12_PRS(i,:),fvalmax_1_c12_PRS,exitflagmax_1_c12_PRS(i,:)]=fmincon(fmax_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));   
    else 
        W2_nlcon_1_llm_c12 = [W2_ls_12_llm(:,1:i-1) W2_ls_12_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_c12 = [t2_ls_12_PRS(:,1:i-1) t2_ls_12_PRS(:,i+1:end)];

        fmin_c12_PRS_1 = @(xmin_c12_1)(xmin_c12_1*W2_ls_12_llm(:,i));
        [xmin_1_c12_PRS(i,:),fvalmin_1_c12_PRS,exitflagmin_1_c12_PRS(i,:)]=fmincon(fmin_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));

        fmax_c12_PRS_1 = @(xmax_c12_1)(-(xmax_c12_1*W2_ls_12_llm(:,i)));
        [xmax_1_c12_PRS(i,:),fvalmax_1_c12_PRS,exitflagmax_1_c12_PRS(i,:)]=fmincon(fmax_c12_PRS_1,sp_c12_PRS,[],[],[],[],lbx_norm_c12,ubx_norm_c12,@(x)nonlconinactive(x,t_c12_PRS,W1_ls_12_llm,W2_nlcon_1_llm_c12,t2_nlcon_1_llm_c12));   
    end
end

x_bnd_c12_PRS_1 = [xmin_1_c12_PRS;xmax_1_c12_PRS];

ind_c12_PRS = 2*(m_c12-n_ls_12);
x_bnd_c12_PRS_1_denorm = denormalizeuv(x_bnd_c12_PRS_1,lbx_c12,ubx_c12);
x_bnd_c12_PRS_org = [x_bnd_c12_PRS_1_denorm zeros(ind_c12_PRS,3)];
[psf_c12_PRS_2,pf_c12_PRS_2,re_c12_PRS_2,beta_c12_PRS_2] = mcspsfconstraint12c(x_bnd_c12_PRS_org,lbx,ubx,sd,Nmcs,Pftarget);
[mean_c12_PRS, std_c12_PRS] = momentsusingyounetal([psf_c12_PRS_2;psf_c12_PRS_1]);
y_c12_PRS = constraint12c(x_bnd_c12_PRS_org);

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

% Bounds in inactive subspace
W2_obj_llm = W_obj_llm(:,n_obj+1:m_obj);
t2_obj_PRS = xopt_obj_PRS_1*W2_obj_llm;

sp_obj_PRS = zeros(1,m_obj);
fmin_obj_PRS_1 = zeros(m_obj-n_obj,1);
fmax_obj_PRS_1 = zeros(m_obj-n_obj,1);
xmin_1_obj_PRS = zeros(m_obj-n_obj,m_obj);
xmax_1_obj_PRS = zeros(m_obj-n_obj,m_obj);
exitflagmin_1_obj_PRS = zeros(m_obj-n_obj,1);
exitflagmax_1_obj_PRS = zeros(m_obj-n_obj,1);

for i=1:m_obj-n_obj
    if i==1
        W2_nlcon_1_llm_obj = W2_obj_llm(:,2:end);
        t2_nlcon_1_llm_obj = t2_obj_PRS(:,2:end);

        fmin_obj_PRS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,i));
        [xmin_1_obj_PRS(i,:),fvalmin_1_obj_PRS,exitflagmin_1_obj_PRS(i,:)]=fmincon(fmin_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));

        fmax_obj_PRS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,i)));
        [xmax_1_obj_PRS(i,:),fvalmax_1_obj_PRS,exitflagmax_1_obj_PRS(i,:)]=fmincon(fmax_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));   
    elseif i==m_obj-n_obj
        W2_nlcon_1_llm_obj = W2_obj_llm(:,1:end-1);
        t2_nlcon_1_llm_obj = t2_obj_PRS(:,1:end-1);

        fmin_obj_PRS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,i));
        [xmin_1_obj_PRS(i,:),fvalmin_1_obj_PRS,exitflagmin_1_obj_PRS(i,:)]=fmincon(fmin_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));

        fmax_obj_PRS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,i)));
        [xmax_1_obj_PRS(i,:),fvalmax_1_obj_PRS,exitflagmax_1_obj_PRS(i,:)]=fmincon(fmax_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));   
    else 
        W2_nlcon_1_llm_obj = [W2_obj_llm(:,1:i-1) W2_obj_llm(:,i+1:end)] ;
        t2_nlcon_1_llm_obj = [t2_obj_PRS(:,1:i-1) t2_obj_PRS(:,i+1:end)];

        fmin_obj_PRS_1 = @(xmin_obj_1)(xmin_obj_1*W2_obj_llm(:,i));
        [xmin_1_obj_PRS(i,:),fvalmin_1_obj_PRS,exitflagmin_1_obj_PRS(i,:)]=fmincon(fmin_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));

        fmax_obj_PRS_1 = @(xmax_obj_1)(-(xmax_obj_1*W2_obj_llm(:,i)));
        [xmax_1_obj_PRS(i,:),fvalmax_1_obj_PRS,exitflagmax_1_obj_PRS(i,:)]=fmincon(fmax_obj_PRS_1,sp_obj_PRS,[],[],[],[],lbx_norm_obj,ubx_norm_obj,@(x)nonlconinactive(x,t_obj_PRS,W1_obj_llm,W2_nlcon_1_llm_obj,t2_nlcon_1_llm_obj));   
    end
end

x_bnd_obj_PRS_1 = [xmin_1_obj_PRS;xmax_1_obj_PRS];

ind_obj_PRS = 2*(m_obj-n_obj);
x_bnd_obj_PRS_1_denorm = denormalizeuv(x_bnd_obj_PRS_1,lbx_obj,ubx_obj);
x_bnd_obj_PRS_org = x_bnd_obj_PRS_1_denorm;
obj_PRS_2 = Weightc(x_bnd_obj_PRS_org);

mean_obj_PRS = mean(obj_PRS_2);
std_obj_PRS = std(obj_PRS_2);
y_obj_PRS = Weightc(x_bnd_obj_PRS_org);

