clc;close all;
%% Weight
t1_obj_mm_KRG = t1_obj_min:t1_obj_range/99:t1_obj_max;
obj_mm_KRG = srgtsKRGEvaluate(t1_obj_mm_KRG',srgt_KRG_obj_1);

h_obj_mm=figure;
scatter(t1_obj_1,y_obj_1,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_obj_mm_KRG,obj_mm_KRG,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('Weight');
%title('KRG Metamodel');
legend({'Objective Function','KRG Metamodel'},'Location','NorthEast');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_obj_mm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Metamodels\KRG\Weigth_KRG_2d.jpg');
% close(h_obj_mm);

%% Constraint 1

t1_c01_mm_KRG = t1_ls_01_min:t1_ls_01_range/99:t1_ls_01_max;
c01_mm_KRG = srgtsKRGEvaluate(t1_c01_mm_KRG',srgt_KRG_c01);

h_c01_mm = figure;
scatter(t1_c01,y_c01,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c01_mm_KRG, c01_mm_KRG,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 1');
%title('PSF 1 KRG Metamodel');
legend({'PSF 1','KRG Metamodel'},'Location','NorthEast');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c01_mm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Metamodels\KRG\PSF 01_KRG_2d.jpg');
%close(h_c01_mm)

%% Constraint 2

t1_c02_mm_KRG = t1_ls_02_min:t1_ls_02_range/99:t1_ls_02_max;
c02_mm_KRG = srgtsKRGEvaluate(t1_c02_mm_KRG',srgt_KRG_c02);

h_c02_mm = figure;
scatter(t1_c02,y_c02,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c02_mm_KRG, c02_mm_KRG,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 2');
%title('PSF 2 KRG Metamodel');
legend({'PSF 2','KRG Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c02_mm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Metamodels\KRG\PSF 02_KRG_2d.jpg');
%close(h_c02_mm)

%% Constraint 3

t1_c03_mm_KRG = t1_ls_03_min:t1_ls_03_range/99:t1_ls_03_max;
c03_mm_KRG = srgtsKRGEvaluate(t1_c03_mm_KRG',srgt_KRG_c03);

h_c03_mm = figure;
scatter(t1_c03,y_c03,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c03_mm_KRG, c03_mm_KRG,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 3');
%title('PSF 3 KRG Metamodel');
legend({'PSF 3','KRG Metamodel'},'Location','NorthEast');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c03_mm,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Metamodels\KRG\PSF 03_KRG_2d.jpg');
%close(h_c03_mm)

