function [b_c01,C_c01,W_c01,Ev_c01,pv_c01, x_c01,...
          srgt_KRG_c01, PRESSRMS_KRG_c01, normalized_PRESSRMS_KRG_c01, eXV_KRG_c01, srgtOPT_KRG_c01, Y_hat_KRG_c01, predvar_KRG, R2_pred_KRG_c01,...
          srgt_PRS_c01, PRESSRMS_PRS_c01, normalized_PRESSRMS_PRS_c01, eXV_PRS_c01, srgtOPT_PRS_c01, Y_hat_PRS_c01, predvar_PRS, R2_pred_PRS_c01,...
          srgt_RBF_c01, PRESSRMS_RBF_c01, normalized_PRESSRMS_RBF_c01, eXV_RBF_c01, srgtOPT_RBF_c01, Y_hat_RBF_c01, R2_pred_RBF_c01,...
          srgt_WAS_c01, PRESSRMS_WAS_c01, normalized_PRESSRMS_WAS_c01, eXV_WAS_c01, srgtOPT_WAS_c01, Y_hat_WAS_c01, predvar_WAS, R2_pred_WAS_c01]=active_subspace_metamodel_c(n_design_var,scalenewDoe1,constraint,mcspsfconstraint,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget)

%Evaluate the limit state
x_cons = scalenewDoe1;
y_c01 = constraint(x_cons);

%Compute the active subspace
k_c=2;
alpha=10;
m=n_design_var;
npoints_active_susbpace = 10*m;
M_c01 = min(ceil(alpha*k_c*log(m)),npoints_active_susbpace-1);
p_c01 = 50;
[b_c01,C_c01,W_c01,Ev_c01]=llm(x_cons(:,1:n_design_var),y_c01,M_c01,p_c01);

%Define the dimension of active subspace
pv_c01=cumsum(diag(Ev_c01))/sum(diag(Ev_c01));
n_c01 = 2;

%DoE in the reduced space
n_doe_c01 = (10*n_c01);
W1_c01 = W_c01(:,1:n_c01);
t1_c01_min = -diag(sign(W1_c01)'*W1_c01);
t1_c01_max = diag(sign(W1_c01)'*W1_c01);
t1_c01_range = t1_c01_max'-t1_c01_min';
N_c01 = 140;
t1_c01_doe = t1_c01_min'+((t1_c01_range).*lhsdesign(N_c01,n_c01));

%Mapping to the original space
sp_c01=zeros(1,n_design_var);
x_c01_mean = mean(x_cons(:,1:n_design_var));
x_c01_std = std(x_cons(:,1:n_design_var));
lbx_c01_norm = zeros(1,n_design_var);
ubx_c01_norm = zeros(1,n_design_var);

for i=1:n_design_var
    lbx_c01_norm(1,i) = ((lbx(1,i)-x_c01_mean(1,i))./x_c01_std(1,i));
    ubx_c01_norm(1,i) = ((ubx(1,i)-x_c01_mean(1,i))./x_c01_std(1,i));
end

xmin_c01 = zeros(N_c01,n_design_var);
fval_c01 = zeros(N_c01,1);
exitflag_c01 = zeros(N_c01,1);

for i=1:N_c01
func =  @(x_opt) (x_opt*zeros(m,1));
[xmin_c01(i,:),fval_c01(i,1),exitflag_c01(i,1)] = fmincon(func,sp_c01,[],[],[],[],lbx_c01_norm,ubx_c01_norm,@(x_opt)nonlconstage1(x_opt,W1_c01,t1_c01_doe(i,:)));
end

% selecting any 10n points from new DoE
ind_c01=(exitflag_c01==1);
S_c01 = sum(ind_c01);
out_c01 = randperm(S_c01);
ind_c01_1=out_c01(1:n_doe_c01);
xmin_c01_1 = xmin_c01(ind_c01,:);
xmin_c01_2 = xmin_c01_1(ind_c01_1,:);

t1_c01_doe_1 = t1_c01_doe(ind_c01,:);
t1_c01_doe_2 = t1_c01_doe_1(ind_c01_1,:);

x_c01 = (xmin_c01_2.*x_c01_std)+x_c01_mean;

%Compute constraint 1
psf_c01 = mcspsfconstraint(x_c01,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget);
 
%Build the metamodel for constraint 1
t_c01 = t1_c01_doe_2;
y_c1 = psf_c01;

%Krgging
[srgt_KRG_c01, PRESSRMS_KRG_c01,  eXV_KRG_c01, srgtOPT_KRG_c01, Y_hat_KRG_c01, predvar_KRG, R2_pred_KRG_c01] = metamodel_KRG(t_c01,y_c1);
normalized_PRESSRMS_KRG_c01=PRESSRMS_KRG_c01/(max(y_c1)-min(y_c1));

%PRS
[srgt_PRS_c01, PRESSRMS_PRS_c01,  eXV_PRS_c01, srgtOPT_PRS_c01, Y_hat_PRS_c01, predvar_PRS, R2_pred_PRS_c01] = metamodel_PRS(t_c01,y_c1);
normalized_PRESSRMS_PRS_c01=PRESSRMS_PRS_c01/(max(y_c1)-min(y_c1));

%RBF
[srgt_RBF_c01, PRESSRMS_RBF_c01,  eXV_RBF_c01, srgtOPT_RBF_c01, Y_hat_RBF_c01, R2_pred_RBF_c01] = metamodel_RBF(t_c01,y_c1);
normalized_PRESSRMS_RBF_c01=PRESSRMS_RBF_c01/(max(y_c1)-min(y_c1));

%WAS
[srgt_WAS_c01, PRESSRMS_WAS_c01,  eXV_WAS_c01, srgtOPT_WAS_c01, Y_hat_WAS_c01, predvar_WAS, R2_pred_WAS_c01] = metamodel_WAS(t_c01,y_c1);
normalized_PRESSRMS_WAS_c01=PRESSRMS_WAS_c01/(max(y_c1)-min(y_c1));

end
