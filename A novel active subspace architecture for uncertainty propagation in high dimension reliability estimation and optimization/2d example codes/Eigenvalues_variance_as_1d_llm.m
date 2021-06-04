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
saveas(h_obj,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Eigen Values\Objective Function_EV_2d.jpg');
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
saveas(h_ls_01,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Eigen Values\Limit State 1_EV_2d.jpg');
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
saveas(h_ls_02,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Eigen Values\Limit State 2_EV_2d.jpg');
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
saveas(h_ls_03,'D:\MS_Thesis_2\Codes\2d example\All Surrogate performance\Journal Figures\Eigen Values\Limit State 03_EV_2d.jpg');
% close(h_ls_03);

