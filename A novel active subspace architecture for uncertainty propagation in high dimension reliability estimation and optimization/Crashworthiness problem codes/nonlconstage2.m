function [c,ceq]= nonlconstage2(x,W1,t,W1_c01,W1_c02,W1_c03,W1_c04,W1_c05,W1_c06,W1_c07,W1_c08,W1_c09,W1_c10,...
 srgt_c01,srgt_c02,srgt_c03,srgt_c04,srgt_c05,srgt_c06,srgt_c07,srgt_c08,srgt_c09,srgt_c10,srgt_evaluate)

x_obj=[x(:,1:5) x(:,7)];
x_c01=([x(2:4) x(6)]);                       
x_c02=(x(1:3));
x_c03=([x(1:3) x(5) x(7)]);
x_c04=([x(1:3) x(5:7)]);
x_c05=([x(2:3) x(7)]);
x_c06=([x(1:3) x(5:7)]);
x_c07=([x(1:3) x(5:7)]);
x_c08=([x(2:4) x(6)]);
x_c09=([x(1:4) x(6)]);
x_c10=([x(3) x(5:7)]);

t_c01=x_c01*W1_c01;
t_c02=x_c02*W1_c02;
t_c03=x_c03*W1_c03;
t_c04=x_c04*W1_c04;
t_c05=x_c05*W1_c05;
t_c06=x_c06*W1_c06;
t_c07=x_c07*W1_c07;
t_c08=x_c08*W1_c08;
t_c09=x_c09*W1_c09;
t_c10=x_c10*W1_c10;

c(1)= 1-srgt_evaluate(t_c01, srgt_c01);
c(2)= 1-srgt_evaluate(t_c02, srgt_c02);
c(3)= 1-srgt_evaluate(t_c03, srgt_c03);
c(4)= 1-srgt_evaluate(t_c04, srgt_c04);
c(5)= 1-srgt_evaluate(t_c05, srgt_c05);
c(6)= 1-srgt_evaluate(t_c06, srgt_c06);
c(7)= 1-srgt_evaluate(t_c07, srgt_c07);
c(8)= 0.98-srgt_evaluate(t_c08, srgt_c08);
c(9)= 1-srgt_evaluate(t_c09, srgt_c09);
c(10)= 1-srgt_evaluate(t_c10, srgt_c10);

ceq = x_obj*W1-t;
end