% clc;clear;close all;
%% Plot the contours of the objective function and the constraints
load('D:\MS Thesis\2 D problem\All Surrogate performance\Experiments_Algorithms\Algorithm_!\paper\Experiment 1.mat');
n=100;
x1 = lbx(1):(ubx(1)-lbx(1))/(n-1):ubx(1);
x2 = lbx(2):(ubx(2)-lbx(2))/(n-1):ubx(2);
[u1,u2]=meshgrid(x1,x2);
U = [u1(:) u2(:)];
U_norm = normalizeuv(U,lbx,ubx);
u1_norm = reshape(U_norm(:,1),n,n);
u2_norm = reshape(U_norm(:,2),n,n);
cont_mm_obj_KRG = srgtsKRGEvaluate(U, srgt_KRG_obj);
cont_mm_psf_KRG_c01 = srgtsKRGEvaluate(U, srgt_KRG_c01);
cont_mm_psf_KRG_c02 = srgtsKRGEvaluate(U, srgt_KRG_c02);
cont_mm_psf_KRG_c03 = srgtsKRGEvaluate(U, srgt_KRG_c03);

reshape_cont_mm_obj_KRG = reshape(cont_mm_obj_KRG,n,n);
reshape_cont_mm_psf_KRG_c01 = reshape(cont_mm_psf_KRG_c01,n,n);
reshape_cont_mm_psf_KRG_c02 = reshape(cont_mm_psf_KRG_c02,n,n);
reshape_cont_mm_psf_KRG_c03 = reshape(cont_mm_psf_KRG_c03,n,n);

%% Contours
h_obj = figure;
contour(u1,u2,reshape_cont_mm_obj_RBF,40,'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Objective Function');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
saveas(h_obj,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\obj.jpg');

h_c1 = figure;
contour(u1,u2,reshape_cont_mm_psf_RBF_c01,40,'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('PSF 1');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
saveas(h_c1,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1.jpg');

h_c2 = figure;
contour(u1,u2,reshape_cont_mm_psf_RBF_c02,40,'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('PSF 2');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
saveas(h_c2,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf2.jpg');

h_c3 = figure;
contour(u1,u2,reshape_cont_mm_psf_RBF_c03,40,'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('PSF 3');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
saveas(h_c3,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf3.jpg');

%% contour for PSF=1
h_psf_1 = figure;
contour(u1,u2,reshape_cont_mm_psf_RBF_c01,[1 1],'b','LineWidth',2);
hold on
contour(u1,u2,reshape_cont_mm_psf_RBF_c02,[1 1],'r','LineWidth',2);
contour(u1,u2,reshape_cont_mm_psf_RBF_c03,[1 1],'g','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Contours of PSF 1, PSF 2, PSF 3 at 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1', 'PSF 2', 'PSF 3');
saveas(h_psf_1,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf123.jpg');

%% contour using proposed algorithm
% clear;close all;
load('D:\MS Thesis\2 D problem\All Surrogate performance\Experiments_Algorithms\Algorithm_4\Uniform\1 dimensional active subspace\11 doe active subspace llm\Experiment 6.mat');
n=100;
x1 = lbx(1):(ubx(1)-lbx(1))/(n-1):ubx(1);
x2 = lbx(2):(ubx(2)-lbx(2))/(n-1):ubx(2);
[u1,u2]=meshgrid(x1,x2);
U = [u1(:) u2(:)];
U_norm = normalizeuv(U,lbx,ubx);
u1_norm = reshape(U_norm(:,1),n,n);
u2_norm = reshape(U_norm(:,2),n,n);
tu_norm_ls_01 = U_norm*W1_ls_01_llm;
tu_norm_ls_02 = U_norm*W1_ls_02_llm;
tu_norm_ls_03 = U_norm*W1_ls_03_llm;

cont_mm_psf_RBF_c01_llm = srgtsRBFEvaluate(tu_norm_ls_01, srgt_RBF_c01);
cont_mm_psf_RBF_c02_llm = srgtsRBFEvaluate(tu_norm_ls_02, srgt_RBF_c02);
cont_mm_psf_RBF_c03_llm = srgtsRBFEvaluate(tu_norm_ls_03, srgt_RBF_c03);

reshape_cont_mm_psf_RBF_c01_llm = reshape(cont_mm_psf_RBF_c01_llm,n,n);
reshape_cont_mm_psf_RBF_c02_llm = reshape(cont_mm_psf_RBF_c02_llm,n,n);
reshape_cont_mm_psf_RBF_c03_llm = reshape(cont_mm_psf_RBF_c03_llm,n,n);

h_psf_1_llm = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,[1 1],'b','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02_llm,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c03_llm,[1 1],'g','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Contours of PSF 1, PSF 2, PSF 3 at 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1', 'PSF 2', 'PSF 3');
saveas(h_psf_1_llm,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf123_llm.jpg');

%% Contours overlap

h_psf_1_overlap_llm = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,[1 1],'b--','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02_llm,[1 1],'r--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c03_llm,[1 1],'g--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01,[1 1],'b','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c03,[1 1],'g','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Contours of PSF 1, PSF 2, PSF 3 at 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1 Active Subspace', 'PSF 2 Active Subspace', 'PSF 3 Active Subspace','PSF 1','PSF 2','PSF 3');
saveas(h_psf_1_overlap_llm,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf123_llm_overlap.jpg');

xmin_act_norm = normalizeuv([3.44 3.28],lbx,ubx);
xmin_RBF_norm = normalizeuv(xmin_RBF,lbx,ubx);
h_psf_1_overlap_llm_journal = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,[1 1],'b--','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02_llm,[1 1],'r--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c03_llm,[1 1],'g--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01,[1 1],'b','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c03,[1 1],'g','LineWidth',2);
scatter(xmin_RBF_norm(1),xmin_RBF_norm(2),'y','filled','MarkerEdgeColor','k');
scatter(xmin_act_norm(1),xmin_act_norm(2),'g','filled','MarkerEdgeColor','k');
xlabel('Normalized x_1');
ylabel('Normalized x_2');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1 Active Subspace', 'PSF 2 Active Subspace', 'PSF 3 Active Subspace','PSF 1','PSF 2','PSF 3','Proposed Optima','Actual Optima');
saveas(h_psf_1_overlap_llm_journal,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf123_llm_overlap_journal.jpg');


%% Plotting 
w11 = W_ls_01_llm(1,1);
w12 = W_ls_01_llm(1,2);
w21 = W_ls_01_llm(2,1);
w22 = W_ls_01_llm(2,2);

h_psf1_w1 = figure;
plot([0 w11],[0 w21],'m','LineWidth',2);
hold on
plot([0 w12],[0 w22],'c','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
xlim([-1 1]);
ylim([-1 1]);
title('PSF 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('W1(active subspace)','W2(inactive subspace)');
saveas(h_psf1_w1,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_w1_w2.jpg');

h_psf1_doe = figure;
plot([0 w11],[0 w21],'m','LineWidth',2);
hold on
plot([0 w12],[0 w22],'c','LineWidth',2);
scatter(x_norm_c01_org_chose(1:6,1),x_norm_c01_org_chose(1:6,2),'y','filled','MarkerEdgeColor','k');
scatter(x_norm_c01_org_chose(8:11,1),x_norm_c01_org_chose(8:11,2),'y','filled','MarkerEdgeColor','k');
xlabel('x_1');
ylabel('x_2');
title('PSF 1');
set(gca,'FontWeight','bold','FontSize',12);
legend({'W1','W2','DoE in Active Subspace'},'Location','NorthWest');
saveas(h_psf1_doe,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_doe.jpg');

h_psf1_mcs=figure;
plot([0 w11],[0 w21],'m','LineWidth',2);
hold on
plot([0 w12],[0 w22],'c','LineWidth',2);
scatter(x_norm_c01_org_chose(1:6,1),x_norm_c01_org_chose(1:6,2),'y','filled','MarkerEdgeColor','k');
scatter(x_norm_c01_org_chose(8:11,1),x_norm_c01_org_chose(8:11,2),'y','filled','MarkerEdgeColor','k');
text([x_norm_c01_org_chose(1:6,1);x_norm_c01_org_chose(8:11,1)],[x_norm_c01_org_chose(1:6,2);x_norm_c01_org_chose(8:11,2)],num2str([y_c01(1:6);y_c01(8:11)]));
xlabel('x_1');
ylabel('x_2');
title('PSF Evaluation using MCS');
set(gca,'FontWeight','bold','FontSize',12);
legend({'W_1(Active Subspace)','W_2(Inactive Subspace)','DoE in Active Subspace'},'Location','NorthWest');
saveas(h_psf1_mcs,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_mcs.jpg');

w11_2 = W_ls_02_llm(1,1);
w12_2 = W_ls_02_llm(1,2);
w21_2 = W_ls_02_llm(2,1);
w22_2 = W_ls_02_llm(2,2);

xmin_norm_RBF = normalizeuv(xmin_RBF,lbx,ubx);
h_psf1_mcs_journal=figure;
plot([0 w11_2],[0 w21_2],'m','LineWidth',2);
hold on
plot([0 w12_2],[0 w22_2],'c','LineWidth',2);
scatter(x_norm_c02_org_chose(1:5,1),x_norm_c02_org_chose(1:5,2),'y','filled','MarkerEdgeColor','k');
scatter(xopt_c02_RBF_1(1),xopt_c02_RBF_1(2),'r','filled','MarkerEdgeColor','k');
scatter(x_bnd_c02_RBF_1(:,1),x_bnd_c02_RBF_1(:,2),'bs','filled','MarkerEdgeColor','k');
scatter(x_norm_c02_org_chose(7:11,1),x_norm_c02_org_chose(7:11,2),'y','filled','MarkerEdgeColor','k');
contour(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c02_llm,[1 1],'g--','LineWidth',2);
xlabel('Normalized x_1');
ylabel('Normalized x_2');
% title('Variation of PSF along inactive subspace');
set(gca,'FontWeight','bold','FontSize',12);
legend({'W_1(Active Subspace)','W_2(Inactive Subspace)','LHS DoE ','sample at PSF =1','bounds along W_2'},'Location','NorthWest');
saveas(h_psf1_mcs_journal,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_mcs_journal.jpg');


t1_ls_01_1 = t1_ls_01_min+0.2:(t1_ls_01_range-0.4)/(n-1):t1_ls_01_max-0.2;
t1_ls_01_2 = t1_ls_01_1';
sp_ls_01 = zeros(1,2);
xmin_ls_01 = zeros(n,2);
fval_ls_01 = zeros(n,1);
exitflag_ls_01 = zeros(n,1);

for i=1:n
func_ls_01 =  @(x_opt) (x_opt*zeros(m_c01,1));
[xmin_ls_01(i,:),fval_ls_01(i,1),exitflag_ls_01(i,1)] = fmincon(func_ls_01,sp_ls_01,[],[],[],[],lbx_norm_c01,ubx_norm_c01,@(x_opt)nonlconstage1(x_opt,W1_ls_01_llm,t1_ls_01_2(i,:)));
end

vis_mm_psf_RBF_c01_llm = srgtsRBFEvaluate(t1_ls_01_2, srgt_RBF_c01);

h_psf1_mm=figure;
plot3([0 w11],[0 w21],[max(y_c01) max(y_c01)],'m','LineWidth',2);
hold on
plot3([0 w12],[0 w22],[max(y_c01) max(y_c01)],'c','LineWidth',2);
scatter3([x_norm_c01_org_chose(1:6,1);x_norm_c01_org_chose(8:11,1)],[x_norm_c01_org_chose(1:6,2);x_norm_c01_org_chose(8:11,2)],[y_c01(1:6);y_c01(8:11)],'y','filled','MarkerEdgeColor','k');
scatter3(xmin_ls_01(:,1),xmin_ls_01(:,2),vis_mm_psf_RBF_c01_llm,20,vis_mm_psf_RBF_c01_llm,'filled');
xlabel('x_1');
ylabel('x_2');
zlabel('PSF 1');
title('PSF 1 metamodel');
set(gca,'FontWeight','bold','FontSize',12);
legend({'W1','W2','DoE in Active Subspace','Metamodel'},'Location','NorthWest');
colorbar;
grid on
saveas(h_psf1_mm,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_mm.jpg');

h_psf1_cont = figure;
plot3([0 w11],[0 w21],[max(y_c01) max(y_c01)],'m','LineWidth',2);
hold on
plot3([0 w12],[0 w22],[max(y_c01) max(y_c01)],'c','LineWidth',2);
scatter3(xmin_ls_01(:,1),xmin_ls_01(:,2),vis_mm_psf_RBF_c01_llm,20,vis_mm_psf_RBF_c01_llm,'filled');
contour3(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,[1 1],'b','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
zlabel('PSF 1');
title('Contour at PSF equals one');
set(gca,'FontWeight','bold','FontSize',12);
legend({'W1','W2','Metamodel','PSF Contour'},'Location','NorthWest');
colorbar;
grid on
saveas(h_psf1_cont,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_cont.jpg');

h_psf1_cont_3d = figure;
contour3(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,40,'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('PSF 1 Active Subspace Contour');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
legend({'Contour'},'Location','NorthWest');
grid on
saveas(h_psf1_cont_3d,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_cont_3d.jpg');

h_psf1_cont_3d_1 = figure;
plot3([0 w11],[0 w21],[max(y_c01) max(y_c01)],'m','LineWidth',2);
hold on
plot3([0 w12],[0 w22],[max(y_c01) max(y_c01)],'c','LineWidth',2);
scatter3(xmin_ls_01(:,1),xmin_ls_01(:,2),vis_mm_psf_RBF_c01_llm,20,vis_mm_psf_RBF_c01_llm,'filled');
contour3(u1_norm,u2_norm,reshape_cont_mm_psf_RBF_c01_llm,[1 1],'LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('PSF 1 Active Subspace Contour');
set(gca,'FontWeight','bold','FontSize',12);
colorbar;
legend({'W1','W2','Metamodel','PSF equals one contour'},'Location','NorthWest');
grid on
saveas(h_psf1_cont_3d_1,'D:\MS Thesis\2 D problem\All Surrogate performance\Individual Meeting-31 March,2021\Contours of original function\psf1_cont_3d_1.jpg');
