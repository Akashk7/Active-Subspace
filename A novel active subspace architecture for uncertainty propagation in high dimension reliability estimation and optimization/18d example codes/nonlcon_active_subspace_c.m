function [c,ceq]= nonlcon_active_subspace_c(x,lbx,ubx,W1_c01,W1_c02,W1_c03,W1_c04,W1_c05,W1_c06,W1_c07,W1_c08,W1_c09,W1_c10,W1_c11,W1_c12,...
srgt_c1,srgt_c2,srgt_c3,srgt_c4,srgt_c5,srgt_c6,srgt_c7,srgt_c8,srgt_c9,srgt_c10,srgt_c11,srgt_c12,srgt_evaluate)

lbx_c01=lbx;
lbx_c02=lbx;
lbx_c03=[lbx(1:16) lbx(18)];
lbx_c04=lbx(1:17);
lbx_c05=lbx;
lbx_c06=lbx;
lbx_c07=lbx;
lbx_c08=lbx;
lbx_c09=lbx;
lbx_c10=lbx;
lbx_c11=lbx(1:15);
lbx_c12=lbx(1:15);

ubx_c01=ubx;
ubx_c02=ubx;
ubx_c03=[ubx(1:16) ubx(18)];
ubx_c04=ubx(1:17);
ubx_c05=ubx;
ubx_c06=ubx;
ubx_c07=ubx;
ubx_c08=ubx;
ubx_c09=ubx;
ubx_c10=ubx;
ubx_c11=ubx(1:15);
ubx_c12=ubx(1:15);

x_c01=normalizeuv(x, lbx_c01, ubx_c01);                       
x_c02=normalizeuv(x, lbx_c02, ubx_c02);
x_c03=normalizeuv([x(1:16) x(18)], lbx_c03, ubx_c03);
x_c04=normalizeuv(x(1:17), lbx_c04, ubx_c04);
x_c05=normalizeuv(x, lbx_c05, ubx_c05);
x_c06=normalizeuv(x, lbx_c06, ubx_c06);
x_c07=normalizeuv(x, lbx_c07, ubx_c07);
x_c08=normalizeuv(x, lbx_c08, ubx_c08);
x_c09=normalizeuv(x, lbx_c09, ubx_c09);
x_c10=normalizeuv(x, lbx_c10, ubx_c10);
x_c11=normalizeuv(x(1:15), lbx_c11, ubx_c11);
x_c12=normalizeuv(x(1:15), lbx_c12, ubx_c12);

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
t_c11=x_c11*W1_c11;
t_c12=x_c12*W1_c12;

c(1)= 1-srgt_evaluate(t_c01, srgt_c1);
c(2)= 1-srgt_evaluate(t_c02, srgt_c2);
c(3)= 1-srgt_evaluate(t_c03, srgt_c3);
c(4)= 1-srgt_evaluate(t_c04, srgt_c4);
c(5)= 1-srgt_evaluate(t_c05, srgt_c5);
c(6)= 1-srgt_evaluate(t_c06, srgt_c6);
c(7)= 1-srgt_evaluate(t_c07, srgt_c7);
c(8)= 1-srgt_evaluate(t_c08, srgt_c8);
c(9)= 1-srgt_evaluate(t_c09, srgt_c9);
c(10)= 1-srgt_evaluate(t_c10, srgt_c10);
c(11)= 1-srgt_evaluate(t_c11, srgt_c11);
c(12)= 1-srgt_evaluate(t_c12, srgt_c12);

ceq=[];

end