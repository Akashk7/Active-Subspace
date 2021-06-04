function [c,ceq]= nonlcon_active_subspace(x,mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,...
               sd1,sd2,sd3,sd4,sd5,sd6,sd7,sd8,sd9,sd10,...
               W1_c01,W1_c02,W1_c03,W1_c04,W1_c05,W1_c06,W1_c07,W1_c08,W1_c09,W1_c10,...
srgt_c1,srgt_c2,srgt_c3,srgt_c4,srgt_c5,srgt_c6,srgt_c7,srgt_c8,srgt_c9,srgt_c10,srgt_evaluate)

% x_c01=[x(2) x(3) x(4) x(6)];
% x_c02=[x(1) x(2) x(3)];
% x_c03=[x(1) x(2) x(3) x(5) x(7)];
% x_c04=[x(1) x(2) x(3) x(5) x(6) x(7)];
% x_c05=[x(2) x(3) x(7)];
% x_c06=[x(1) x(2) x(3) x(5) x(6) x(7)];
% x_c07=[x(1) x(2) x(3) x(5) x(6) x(7)];
% x_c08=[x(2) x(3) x(4) x(6)];
% x_c09=[x(1) x(2) x(3) x(4) x(6)];
% x_c10=[x(3) x(5) x(6) x(7)];


x_c01=normalizerv([x(2:4) x(6)], mu1, sd1);                       
x_c02=normalizerv(x(1:3), mu2, sd2);
x_c03=normalizerv([x(1:3) x(5) x(7)], mu3, sd3);
x_c04=normalizerv([x(1:3) x(5:7)], mu4, sd4);
x_c05=normalizerv([x(2:3) x(7)], mu5, sd5);
x_c06=normalizerv([x(1:3) x(5:7)], mu6, sd6);
x_c07=normalizerv([x(1:3) x(5:7)], mu7, sd7);
x_c08=normalizerv([x(2:4) x(6)], mu8, sd8);
x_c09=normalizerv([x(1:4) x(6)], mu9, sd9);
x_c10=normalizerv([x(3) x(5:7)], mu10, sd10);

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

c(1)= 1-srgt_evaluate(t_c01, srgt_c1);
c(2)= 1-srgt_evaluate(t_c02, srgt_c2);
c(3)= 1-srgt_evaluate(t_c03, srgt_c3);
c(4)= 1-srgt_evaluate(t_c04, srgt_c4);
c(5)= 1-srgt_evaluate(t_c05, srgt_c5);
c(6)= 1-srgt_evaluate(t_c06, srgt_c6);
c(7)= 1-srgt_evaluate(t_c07, srgt_c7);
c(8)= 0.98-srgt_evaluate(t_c08, srgt_c8);
c(9)= 1-srgt_evaluate(t_c09, srgt_c9);
c(10)= 1-srgt_evaluate(t_c10, srgt_c10);
ceq=[];

end