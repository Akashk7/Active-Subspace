clc;close all;
%% Weight
t1_obj_mm_WAS = t1_obj_min:t1_obj_range/99:t1_obj_max;
obj_mm_WAS = srgtsWASEvaluate(t1_obj_mm_WAS',srgt_WAS_obj_1);

h_obj_mm=figure;
scatter(t1_obj_1,y_obj_1,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_obj_mm_WAS,obj_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('Weight');
%title('WAS Metamodel');
legend({'Weight','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_obj_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\Weigth_WAS.jpg');
% close(h_obj_mm);

%% Constraint 1

t1_c01_mm_WAS = t1_ls_01_min:t1_ls_01_range/99:t1_ls_01_max;
c01_mm_WAS = srgtsWASEvaluate(t1_c01_mm_WAS',srgt_WAS_c01);

h_c01_mm = figure;
scatter(t1_c01,y_c01,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c01_mm_WAS, c01_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 1');
%title('PSF 1 WAS Metamodel');
legend({'PSF 1','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c01_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 01_WAS.jpg');
%close(h_c01_mm)

%% Constraint 2

t1_c02_mm_WAS = t1_ls_02_min:t1_ls_02_range/99:t1_ls_02_max;
c02_mm_WAS = srgtsWASEvaluate(t1_c02_mm_WAS',srgt_WAS_c02);

h_c02_mm = figure;
scatter(t1_c02,y_c02,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c02_mm_WAS, c02_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 2');
%title('PSF 2 WAS Metamodel');
legend({'PSF 2','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c02_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 02_WAS.jpg');
%close(h_c02_mm)

%% Constraint 3

t1_c03_mm_WAS = t1_ls_03_min:t1_ls_03_range/99:t1_ls_03_max;
c03_mm_WAS = srgtsWASEvaluate(t1_c03_mm_WAS',srgt_WAS_c03);

h_c03_mm = figure;
scatter(t1_c03,y_c03,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c03_mm_WAS, c03_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 3');
%title('PSF 3 WAS Metamodel');
legend({'PSF 3','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c03_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 03_WAS.jpg');
%close(h_c03_mm)

%% Constraint 4

t1_c04_mm_WAS = t1_ls_04_min:t1_ls_04_range/99:t1_ls_04_max;
c04_mm_WAS = srgtsWASEvaluate(t1_c04_mm_WAS',srgt_WAS_c04);

h_c04_mm = figure;
scatter(t1_c04,y_c04,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c04_mm_WAS, c04_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 4');
%title('PSF 4 WAS Metamodel');
legend({'PSF 4','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c04_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 04_WAS.jpg');
%close(h_c04_mm)

%% Constraint 5

t1_c05_mm_WAS = t1_ls_05_min:t1_ls_05_range/99:t1_ls_05_max;
c05_mm_WAS = srgtsWASEvaluate(t1_c05_mm_WAS',srgt_WAS_c05);

h_c05_mm = figure;
scatter(t1_c05,y_c05,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c05_mm_WAS, c05_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 5');
%title('PSF 5 WAS Metamodel');
legend({'PSF 5','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c05_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 05_WAS.jpg');
%close(h_c05_mm)

%% Constraint 6

t1_c06_mm_WAS = t1_ls_06_min:t1_ls_06_range/99:t1_ls_06_max;
c06_mm_WAS = srgtsWASEvaluate(t1_c06_mm_WAS',srgt_WAS_c06);

h_c06_mm = figure;
scatter(t1_c06,y_c06,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c06_mm_WAS, c06_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 6');
%title('PSF 6 WAS Metamodel');
legend({'PSF 6','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c06_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 06_WAS.jpg');
%close(h_c06_mm)

%% Constraint 7

t1_c07_mm_WAS = t1_ls_07_min:t1_ls_07_range/99:t1_ls_07_max;
c07_mm_WAS = srgtsWASEvaluate(t1_c07_mm_WAS',srgt_WAS_c07);

h_c07_mm = figure;
scatter(t1_c07,y_c07,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c07_mm_WAS, c07_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 7');
%title('PSF 7 WAS Metamodel');
legend({'PSF 7','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c07_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 07_WAS.jpg');
%close(h_c07_mm)

%% Constraint 8

t1_c08_mm_WAS = t1_ls_08_min:t1_ls_08_range/99:t1_ls_08_max;
c08_mm_WAS = srgtsWASEvaluate(t1_c08_mm_WAS',srgt_WAS_c08);

h_c08_mm = figure;
scatter(t1_c08,y_c08,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c08_mm_WAS, c08_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 8');
%title('PSF 8 WAS Metamodel');
legend({'PSF 8','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c08_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 08_WAS.jpg');
%close(h_c08_mm)

%% Constraint 9

t1_c09_mm_WAS = t1_ls_09_min:t1_ls_09_range/99:t1_ls_09_max;
c09_mm_WAS = srgtsWASEvaluate(t1_c09_mm_WAS',srgt_WAS_c09);

h_c09_mm = figure;
scatter(t1_c09,y_c09,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c09_mm_WAS, c09_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 9');
%title('PSF 9 WAS Metamodel');
legend({'PSF 9','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c09_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 09_WAS.jpg');
%close(h_c09_mm)

%% Constraint 10

t1_c10_mm_WAS = t1_ls_10_min:t1_ls_10_range/99:t1_ls_10_max;
c10_mm_WAS = srgtsWASEvaluate(t1_c10_mm_WAS',srgt_WAS_c10);

h_c10_mm = figure;
scatter(t1_c10,y_c10,'r','filled','MarkerEdgeColor','k');
hold on
plot(t1_c10_mm_WAS, c10_mm_WAS,'b-','LineWidth',2);
xlabel('Active Variable');
ylabel('PSF 10');
%title('PSF 10 WAS Metamodel');
legend({'PSF 10','WAS Metamodel'},'Location','NorthWest');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_c10_mm,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Metamodels\WAS\PSF 10_WAS.jpg');
%close(h_c10_mm)
