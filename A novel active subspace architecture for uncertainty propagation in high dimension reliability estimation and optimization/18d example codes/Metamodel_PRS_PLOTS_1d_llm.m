clc;close all;
%% Weight
t1_obj_mm_PRS = t1_obj_min:t1_obj_range/99:t1_obj_max;
obj_mm_PRS = srgtsPRSEvaluate(t1_obj_mm_PRS',srgt_PRS_obj_1);

h_obj_mm=figure;
scatter(t1_obj_1,y_obj_1,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_obj_mm_PRS,obj_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('Objective Function');
%title('PRS Metamodel');
legend({'Objective Function','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_obj_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\Weigth_PRS_18d.jpg');
% close(h_obj_mm);

%% Constraint 1

t1_c01_mm_PRS = t1_ls_01_min:t1_ls_01_range/99:t1_ls_01_max;
c01_mm_PRS = srgtsPRSEvaluate(t1_c01_mm_PRS',srgt_PRS_c01);

h_c01_mm = figure;
scatter(t1_c01,y_c01,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c01_mm_PRS, c01_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 1');
%title('PSF 1 PRS Metamodel');
legend({'PSF 1','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c01_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 01_PRS_18d.jpg');
%close(h_c01_mm)

%% Constraint 2

t1_c02_mm_PRS = t1_ls_02_min:t1_ls_02_range/99:t1_ls_02_max;
c02_mm_PRS = srgtsPRSEvaluate(t1_c02_mm_PRS',srgt_PRS_c02);

h_c02_mm = figure;
scatter(t1_c02,y_c02,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c02_mm_PRS, c02_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 2');
%title('PSF 2 PRS Metamodel');
legend({'PSF 2','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c02_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 02_PRS_18d.jpg');
%close(h_c02_mm)

%% Constraint 3

t1_c03_mm_PRS = t1_ls_03_min:t1_ls_03_range/99:t1_ls_03_max;
c03_mm_PRS = srgtsPRSEvaluate(t1_c03_mm_PRS',srgt_PRS_c03);

h_c03_mm = figure;
scatter(t1_c03,y_c03,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c03_mm_PRS, c03_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 3');
%title('PSF 3 PRS Metamodel');
legend({'PSF 3','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c03_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 03_PRS_18d.jpg');
%close(h_c03_mm)

%% Constraint 4

t1_c04_mm_PRS = t1_ls_04_min:t1_ls_04_range/99:t1_ls_04_max;
c04_mm_PRS = srgtsPRSEvaluate(t1_c04_mm_PRS',srgt_PRS_c04);

h_c04_mm = figure;
scatter(t1_c04,y_c04,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c04_mm_PRS, c04_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 4');
%title('PSF 4 PRS Metamodel');
legend({'PSF 4','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c04_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 04_PRS_18d.jpg');
%close(h_c04_mm)

%% Constraint 5

t1_c05_mm_PRS = t1_ls_05_min:t1_ls_05_range/99:t1_ls_05_max;
c05_mm_PRS = srgtsPRSEvaluate(t1_c05_mm_PRS',srgt_PRS_c05);

h_c05_mm = figure;
scatter(t1_c05,y_c05,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c05_mm_PRS, c05_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 5');
%title('PSF 5 PRS Metamodel');
legend({'PSF 5','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c05_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 05_PRS_18d.jpg');
%close(h_c05_mm)

%% Constraint 6

t1_c06_mm_PRS = t1_ls_06_min:t1_ls_06_range/99:t1_ls_06_max;
c06_mm_PRS = srgtsPRSEvaluate(t1_c06_mm_PRS',srgt_PRS_c06);

h_c06_mm = figure;
scatter(t1_c06,y_c06,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c06_mm_PRS, c06_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 6');
%title('PSF 6 PRS Metamodel');
legend({'PSF 6','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c06_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 06_PRS_18d.jpg');
%close(h_c06_mm)

%% Constraint 7

t1_c07_mm_PRS = t1_ls_07_min:t1_ls_07_range/99:t1_ls_07_max;
c07_mm_PRS = srgtsPRSEvaluate(t1_c07_mm_PRS',srgt_PRS_c07);

h_c07_mm = figure;
scatter(t1_c07,y_c07,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c07_mm_PRS, c07_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 7');
%title('PSF 7 PRS Metamodel');
legend({'PSF 7','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c07_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 07_PRS_18d.jpg');
%close(h_c07_mm)

%% Constraint 8

t1_c08_mm_PRS = t1_ls_08_min:t1_ls_08_range/99:t1_ls_08_max;
c08_mm_PRS = srgtsPRSEvaluate(t1_c08_mm_PRS',srgt_PRS_c08);

h_c08_mm = figure;
scatter(t1_c08,y_c08,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c08_mm_PRS, c08_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 8');
%title('PSF 8 PRS Metamodel');
legend({'PSF 8','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c08_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 08_PRS_18d.jpg');
%close(h_c08_mm)

%% Constraint 9

t1_c09_mm_PRS = t1_ls_09_min:t1_ls_09_range/99:t1_ls_09_max;
c09_mm_PRS = srgtsPRSEvaluate(t1_c09_mm_PRS',srgt_PRS_c09);

h_c09_mm = figure;
scatter(t1_c09,y_c09,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c09_mm_PRS, c09_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 9');
%title('PSF 9 PRS Metamodel');
legend({'PSF 9','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c09_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 09_PRS_18d.jpg');
%close(h_c09_mm)

%% Constraint 10

t1_c10_mm_PRS = t1_ls_10_min:t1_ls_10_range/99:t1_ls_10_max;
c10_mm_PRS = srgtsPRSEvaluate(t1_c10_mm_PRS',srgt_PRS_c10);

h_c10_mm = figure;
scatter(t1_c10,y_c10,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c10_mm_PRS, c10_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 10');
%title('PSF 10 PRS Metamodel');
legend({'PSF 10','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c10_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 10_PRS_18d.jpg');
%close(h_c10_mm)
%% Constraint 11

t1_c11_mm_PRS = t1_ls_11_min:t1_ls_11_range/99:t1_ls_11_max;
c11_mm_PRS = srgtsPRSEvaluate(t1_c11_mm_PRS',srgt_PRS_c11);

h_c11_mm = figure;
scatter(t1_c11,y_c11,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c11_mm_PRS, c11_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 11');
%title('PSF 10 PRS Metamodel');
legend({'PSF 11','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c11_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 11_PRS_18d.jpg');
%close(h_c11_mm)

%% Constraint 12

t1_c12_mm_PRS = t1_ls_12_min:t1_ls_12_range/99:t1_ls_12_max;
c12_mm_PRS = srgtsPRSEvaluate(t1_c12_mm_PRS',srgt_PRS_c12);

h_c12_mm = figure;
scatter(t1_c12,y_c12,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c12_mm_PRS, c12_mm_PRS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 12');
%title('PSF 10 PRS Metamodel');
legend({'PSF 12','PRS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c12_mm,'D:\MS Thesis\18 d no dummy variables\All Surrogate Performance\Journal Figures\Metamodels\PRS\PSF 12_PRS_18d.jpg');
%close(h_c12_mm)
