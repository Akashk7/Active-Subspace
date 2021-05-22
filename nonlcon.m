function [c,ceq]= nonlcon(x,srgt_c1,srgt_c2,srgt_c3,srgt_evaluate)

x_c01=x;                       
x_c02=x;
x_c03=x;

c(1)= 1-srgt_evaluate(x_c01, srgt_c1);
c(2)= 1-srgt_evaluate(x_c02, srgt_c2);
c(3)= 1-srgt_evaluate(x_c03, srgt_c3);
ceq=[];

end