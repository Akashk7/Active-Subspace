function [c,ceq]= nonlcon_org_c(x,...
srgt_ls_1,srgt_ls_2,srgt_ls_3,srgt_ls_4,srgt_ls_5,srgt_ls_6,srgt_ls_7,srgt_ls_8,srgt_ls_9,srgt_ls_10,srgt_evaluate)

x_ls_01=[x(2:4) x(6)];                       
x_ls_02=x(1:3);
x_ls_03=[x(1:3) x(5) x(7)];
x_ls_04=[x(1:3) x(5:7)];
x_ls_05=[x(2:3) x(7)];
x_ls_06=[x(1:3) x(5:7)];
x_ls_07=[x(1:3) x(5:7)];
x_ls_08=[x(2:4) x(6)];
x_ls_09=[x(1:4) x(6)];
x_ls_10=[x(3) x(5:7)];

c(1)= 1-srgt_evaluate(x_ls_01, srgt_ls_1);
c(2)= 1-srgt_evaluate(x_ls_02, srgt_ls_2);
c(3)= 1-srgt_evaluate(x_ls_03, srgt_ls_3);
c(4)= 1-srgt_evaluate(x_ls_04, srgt_ls_4);
c(5)= 1-srgt_evaluate(x_ls_05, srgt_ls_5);
c(6)= 1-srgt_evaluate(x_ls_06, srgt_ls_6);
c(7)= 1-srgt_evaluate(x_ls_07, srgt_ls_7);
c(8)= 1-srgt_evaluate(x_ls_08, srgt_ls_8);
c(9)= 1-srgt_evaluate(x_ls_09, srgt_ls_9);
c(10)= 1-srgt_evaluate(x_ls_10, srgt_ls_10);
ceq=[];

end