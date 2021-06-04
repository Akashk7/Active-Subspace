clc;close all;
%% Weight
t1_obj_mm_RBF = t1_obj_min:t1_obj_range/99:t1_obj_max;
obj_mm_RBF = srgtsRBFEvaluate(t1_obj_mm_RBF',srgt_RBF_obj_1);

h_obj_mm=figure;
scatter(t1_obj_1,y_obj_1,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_obj_mm_RBF,obj_mm_RBF,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('Weight');
%title('RBF Metamodel');
legend({'Objective Function','RBF Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_obj_mm,'D:\MS Thesis\2 D problem\All Surrogate Performance\Journal Figures\Metamodels\RBF\Weigth_RBF.jpg');
% close(h_obj_mm);

%% Constraint 1

t1_c01_mm_RBF = t1_ls_01_min:t1_ls_01_range/99:t1_ls_01_max;
c01_mm_RBF = srgtsRBFEvaluate(t1_c01_mm_RBF',srgt_RBF_c01);

h_c01_mm = figure;
scatter(t1_c01,y_c01,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c01_mm_RBF, c01_mm_RBF,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 1');
%title('PSF 1 RBF Metamodel');
legend({'PSF 1','RBF Metamodel'},'Location','NorthEast');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c01_mm,'D:\MS Thesis\2 D problem\All Surrogate Performance\Journal Figures\Metamodels\RBF\PSF 01_RBF.jpg');
%close(h_c01_mm)

%% Constraint 2

t1_c02_mm_RBF = t1_ls_02_min:t1_ls_02_range/99:t1_ls_02_max;
c02_mm_RBF = srgtsRBFEvaluate(t1_c02_mm_RBF',srgt_RBF_c02);

h_c02_mm = figure;
scatter(t1_c02,y_c02,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c02_mm_RBF, c02_mm_RBF,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 2');
%title('PSF 2 RBF Metamodel');
legend({'PSF 2','RBF Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c02_mm,'D:\MS Thesis\2 D problem\All Surrogate Performance\Journal Figures\Metamodels\RBF\PSF 02_RBF.jpg');
%close(h_c02_mm)

%% Constraint 3

t1_c03_mm_RBF = t1_ls_03_min:t1_ls_03_range/99:t1_ls_03_max;
c03_mm_RBF = srgtsRBFEvaluate(t1_c03_mm_RBF',srgt_RBF_c03);

h_c03_mm = figure;
scatter(t1_c03,y_c03,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c03_mm_RBF, c03_mm_RBF,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 3');
%title('PSF 3 RBF Metamodel');
legend({'PSF 3','RBF Metamodel'},'Location','NorthEast');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c03_mm,'D:\MS Thesis\2 D problem\All Surrogate Performance\Journal Figures\Metamodels\RBF\PSF 03_RBF.jpg');
%close(h_c03_mm)

