clc;close all;
%% Weight
Mboot = 1e4;
[W1star_obj_llm,Lstar_obj_llm]=subspacevar(b_obj_llm,Mboot);
dim_obj_llm = [1:1:m_obj]'; 
min_Lstar_obj_llm = min(Lstar_obj_llm,[],2);
max_Lstar_obj_llm = max(Lstar_obj_llm,[],2);
x_polyshape_obj = [1:1:m_obj m_obj:-1:1];
y_polyshape_obj = ([min_Lstar_obj_llm' [max_Lstar_obj_llm(end:-1:1)]']);

h_obj = figure;
semilogy(dim_obj_llm,(diag(Ev_obj_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_obj,abs(y_polyshape_obj),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Weight Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_obj,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Weight_EV.jpg');
% close(h_obj);

%% Constraint 1
[W1star_ls_01_llm,Lstar_ls_01_llm]=subspacevar(b_ls_01_llm,Mboot);
dim_ls_01_llm = [1:1:m_c01]'; 
min_Lstar_ls_01_llm = min(Lstar_ls_01_llm,[],2);
max_Lstar_ls_01_llm = max(Lstar_ls_01_llm,[],2);
x_polyshape_ls_01 = [1:1:m_c01 m_c01:-1:1];
y_polyshape_ls_01 = ([min_Lstar_ls_01_llm' [max_Lstar_ls_01_llm(end:-1:1)]']);

h_ls_01 = figure;
semilogy(dim_ls_01_llm,(diag(Ev_ls_01_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_01,abs(y_polyshape_ls_01),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 1 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_01,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 1_EV.jpg');
% close(h_ls_01);

%% Constraint 2

[W1star_ls_02_llm,Lstar_ls_02_llm]=subspacevar(b_ls_02_llm,Mboot);
dim_ls_02_llm = [1:1:m_c02]'; 
min_Lstar_ls_02_llm = min(Lstar_ls_02_llm,[],2);
max_Lstar_ls_02_llm = max(Lstar_ls_02_llm,[],2);
x_polyshape_ls_02 = [1:1:m_c02 m_c02:-1:1];
y_polyshape_ls_02 = ([min_Lstar_ls_02_llm' [max_Lstar_ls_02_llm(end:-1:1)]']);

h_ls_02 = figure;
semilogy(dim_ls_02_llm,(diag(Ev_ls_02_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_02,abs(y_polyshape_ls_02),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 2 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_02,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 2_EV.jpg');
% close(h_ls_02);

%% Constraint 3

[W1star_ls_03_llm,Lstar_ls_03_llm]=subspacevar(b_ls_03_llm,Mboot);
dim_ls_03_llm = [1:1:m_c03]'; 
min_Lstar_ls_03_llm = min(Lstar_ls_03_llm,[],2);
max_Lstar_ls_03_llm = max(Lstar_ls_03_llm,[],2);
x_polyshape_ls_03 = [1:1:m_c03 m_c03:-1:1];
y_polyshape_ls_03 = ([min_Lstar_ls_03_llm' [max_Lstar_ls_03_llm(end:-1:1)]']);

h_ls_03 = figure;
semilogy(dim_ls_03_llm,(diag(Ev_ls_03_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_03,abs(y_polyshape_ls_03),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 3 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_03,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 03_EV.jpg');
% close(h_ls_03);

%% Constraint 4

[W1star_ls_04_llm,Lstar_ls_04_llm]=subspacevar(b_ls_04_llm,Mboot);
dim_ls_04_llm = [1:1:m_c04]'; 
min_Lstar_ls_04_llm = min(Lstar_ls_04_llm,[],2);
max_Lstar_ls_04_llm = max(Lstar_ls_04_llm,[],2);
x_polyshape_ls_04 = [1:1:m_c04 m_c04:-1:1];
y_polyshape_ls_04 = ([min_Lstar_ls_04_llm' [max_Lstar_ls_04_llm(end:-1:1)]']);

h_ls_04 = figure;
semilogy(dim_ls_04_llm,(diag(Ev_ls_04_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_04,abs(y_polyshape_ls_04),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 4 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_04,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 04_EV.jpg');
% close(h_ls_04);

%% Constraint 5

[W1star_ls_05_llm,Lstar_ls_05_llm]=subspacevar(b_ls_05_llm,Mboot);
dim_ls_05_llm = [1:1:m_c05]'; 
min_Lstar_ls_05_llm = min(Lstar_ls_05_llm,[],2);
max_Lstar_ls_05_llm = max(Lstar_ls_05_llm,[],2);
x_polyshape_ls_05 = [1:1:m_c05 m_c05:-1:1];
y_polyshape_ls_05 = ([min_Lstar_ls_05_llm' [max_Lstar_ls_05_llm(end:-1:1)]']);

h_ls_05 = figure;
semilogy(dim_ls_05_llm,(diag(Ev_ls_05_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_05,abs(y_polyshape_ls_05),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 5 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_05,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 05_EV.jpg');
% close(h_ls_05);

%% Constraint 06

[W1star_ls_06_llm,Lstar_ls_06_llm]=subspacevar(b_ls_06_llm,Mboot);
dim_ls_06_llm = [1:1:m_c06]'; 
min_Lstar_ls_06_llm = min(Lstar_ls_06_llm,[],2);
max_Lstar_ls_06_llm = max(Lstar_ls_06_llm,[],2);
x_polyshape_ls_06 = [1:1:m_c06 m_c06:-1:1];
y_polyshape_ls_06 = ([min_Lstar_ls_06_llm' [max_Lstar_ls_06_llm(end:-1:1)]']);

h_ls_06 = figure;
semilogy(dim_ls_06_llm,(diag(Ev_ls_06_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_06,abs(y_polyshape_ls_06),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 6 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_06,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 06_EV.jpg');
% close(h_ls_06);

%% Constraint 7

[W1star_ls_07_llm,Lstar_ls_07_llm]=subspacevar(b_ls_07_llm,Mboot);
dim_ls_07_llm = [1:1:m_c07]'; 
min_Lstar_ls_07_llm = min(Lstar_ls_07_llm,[],2);
max_Lstar_ls_07_llm = max(Lstar_ls_07_llm,[],2);
x_polyshape_ls_07 = [1:1:m_c07 m_c07:-1:1];
y_polyshape_ls_07 = ([min_Lstar_ls_07_llm' [max_Lstar_ls_07_llm(end:-1:1)]']);

h_ls_07 = figure;
semilogy(dim_ls_07_llm,(diag(Ev_ls_07_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_07,abs(y_polyshape_ls_07),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 7 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_07,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 07_EV.jpg');
% close(h_ls_07);

%% Constraint 8

[W1star_ls_08_llm,Lstar_ls_08_llm]=subspacevar(b_ls_08_llm,Mboot);
dim_ls_08_llm = [1:1:m_c08]'; 
min_Lstar_ls_08_llm = min(Lstar_ls_08_llm,[],2);
max_Lstar_ls_08_llm = max(Lstar_ls_08_llm,[],2);
x_polyshape_ls_08 = [1:1:m_c08 m_c08:-1:1];
y_polyshape_ls_08 = ([min_Lstar_ls_08_llm' [max_Lstar_ls_08_llm(end:-1:1)]']);

h_ls_08 = figure;
semilogy(dim_ls_08_llm,(diag(Ev_ls_08_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_08,abs(y_polyshape_ls_08),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 8 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_08,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 08_EV.jpg');
% close(h_ls_08);

%% Constraint 9

[W1star_ls_09_llm,Lstar_ls_09_llm]=subspacevar(b_ls_09_llm,Mboot);
dim_ls_09_llm = [1:1:m_c09]'; 
min_Lstar_ls_09_llm = min(Lstar_ls_09_llm,[],2);
max_Lstar_ls_09_llm = max(Lstar_ls_09_llm,[],2);
x_polyshape_ls_09 = [1:1:m_c09 m_c09:-1:1];
y_polyshape_ls_09 = ([min_Lstar_ls_09_llm' [max_Lstar_ls_09_llm(end:-1:1)]']);

h_ls_09 = figure;
semilogy(dim_ls_09_llm,(diag(Ev_ls_09_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_09,abs(y_polyshape_ls_09),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 9 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_09,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 09_EV.jpg');
% close(h_ls_09);

%% Constraint 10

[W1star_ls_10_llm,Lstar_ls_10_llm]=subspacevar(b_ls_10_llm,Mboot);
dim_ls_10_llm = [1:1:m_c10]'; 
min_Lstar_ls_10_llm = min(Lstar_ls_10_llm,[],2);
max_Lstar_ls_10_llm = max(Lstar_ls_10_llm,[],2);
x_polyshape_ls_10 = [1:1:m_c10 m_c10:-1:1];
y_polyshape_ls_10 = ([min_Lstar_ls_10_llm' [max_Lstar_ls_10_llm(end:-1:1)]']);

h_ls_10 = figure;
semilogy(dim_ls_10_llm,(diag(Ev_ls_10_llm)),'o-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
hold on
patch(x_polyshape_ls_10,abs(y_polyshape_ls_10),'y','EdgeColor','y');
alpha(0.3);
xlabel('Dimensions');
ylabel('Eigenvalues');
% title('Limit State 10 Eigen Value Distribution');
legend('EigenValue','Bootstrap Interval');
set(gca,'FontWeight','bold','FontSize',12);
grid on
saveas(h_ls_10,'D:\MS Thesis\Car crash worthiness - 7 dimension\All Surrogate Performance\Journal Figures\Eigen Values\Limit State 10_EV.jpg');
% close(h_ls_10);
