%% contour using proposed algorithm
% clear;close all;
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

load('D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Experiments_Algorithms\Algorithm_4\Uniform\1 dimensional active subspace\11 doe active subspace llm\Experiments\Experiment 5star.mat');
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

cont_mm_psf_KRG_c01_llm = srgtsKRGEvaluate(tu_norm_ls_01, srgt_KRG_c01);
cont_mm_psf_KRG_c02_llm = srgtsKRGEvaluate(tu_norm_ls_02, srgt_KRG_c02);
cont_mm_psf_KRG_c03_llm = srgtsKRGEvaluate(tu_norm_ls_03, srgt_KRG_c03);

reshape_cont_mm_psf_KRG_c01_llm = reshape(cont_mm_psf_KRG_c01_llm,n,n);
reshape_cont_mm_psf_KRG_c02_llm = reshape(cont_mm_psf_KRG_c02_llm,n,n);
reshape_cont_mm_psf_KRG_c03_llm = reshape(cont_mm_psf_KRG_c03_llm,n,n);

h_psf_1_llm = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c01_llm,1,'b','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c02_llm,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c03_llm,[1 1],'g','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Contours of PSF 1, PSF 2, PSF 3 at 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1', 'PSF 2', 'PSF 3');
saveas(h_psf_1_llm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Comparison\psf123_llm.jpg');

%% Contours overlap

h_psf_1_overlap_llm = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c01_llm,[1 1],'b--','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c02_llm,[1 1],'r--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c03_llm,[1 1],'g--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c01,[1 1],'b','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c02,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c03,[1 1],'g','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title('Contours of PSF 1, PSF 2, PSF 3 at 1');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1 Active Subspace', 'PSF 2 Active Subspace', 'PSF 3 Active Subspace','PSF 1','PSF 2','PSF 3');
saveas(h_psf_1_overlap_llm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Comparison\psf123_llm_overlap.jpg');

xmin_act_norm = normalizeuv([3.44 3.28],lbx,ubx);
xmin_KRG_norm = normalizeuv(xmin_KRG,lbx,ubx);
h_psf_1_overlap_llm_journal = figure;
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c01_llm,[1 1],'b--','LineWidth',2);
hold on
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c02_llm,[1 1],'r--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c03_llm,[1 1],'g--','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c01,[1 1],'b','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c02,[1 1],'r','LineWidth',2);
contour(u1_norm,u2_norm,reshape_cont_mm_psf_KRG_c03,[1 1],'g','LineWidth',2);
scatter(xmin_KRG_norm(1),xmin_KRG_norm(2),'y','filled','MarkerEdgeColor','k');
scatter(xmin_act_norm(1),xmin_act_norm(2),'g','filled','MarkerEdgeColor','k');
xlabel('Normalized x_1');
ylabel('Normalized x_2');
set(gca,'FontWeight','bold','FontSize',12);
legend('PSF 1 Active Subspace', 'PSF 2 Active Subspace', 'PSF 3 Active Subspace','PSF 1','PSF 2','PSF 3','Proposed Optima','Actual Optima');
saveas(h_psf_1_overlap_llm_journal,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Comparison\psf123_llm_overlap_journal.jpg');
