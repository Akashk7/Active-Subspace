function [c,ceq]= nonlcon(x,srgt_c1,srgt_c2,srgt_c3,srgt_c4,srgt_c5,srgt_c6,...
                           srgt_c7,srgt_c8,srgt_c9,srgt_c10,srgt_evaluate)
x_c01=[x(2:4) x(6)];                       
x_c02=x(1:3);
x_c03=[x(1:3) x(5) x(7)];
x_c04=[x(1:3) x(5:7)];
x_c05=[x(2:3) x(7)];
x_c06=[x(1:3) x(5:7)];
x_c07=[x(1:3) x(5:7)];
x_c08=[x(2:4) x(6)];
x_c09=[x(1:4) x(6)];
x_c10=[x(3) x(5:7)];

c(1)= 1-srgt_evaluate(x_c01, srgt_c1);
c(2)= 1-srgt_evaluate(x_c02, srgt_c2);
c(3)= 1-srgt_evaluate(x_c03, srgt_c3);
c(4)= 1-srgt_evaluate(x_c04, srgt_c4);
c(5)= 1-srgt_evaluate(x_c05, srgt_c5);
c(6)= 1-srgt_evaluate(x_c06, srgt_c6);
c(7)= 1-srgt_evaluate(x_c07, srgt_c7);
c(8)= 0.99-srgt_evaluate(x_c08, srgt_c8);
c(9)= 1-srgt_evaluate(x_c09, srgt_c9);
c(10)= 1-srgt_evaluate(x_c10, srgt_c10);
ceq=[];

end