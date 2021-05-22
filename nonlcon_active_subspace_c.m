function [c,ceq]= nonlcon_active_subspace_c(x,lbx,ubx,W1_c01,W1_c02,W1_c03,...
srgt_c1,srgt_c2,srgt_c3,srgt_evaluate)

lbx_c01 = lbx;
lbx_c02 = lbx;
lbx_c03 = lbx;

ubx_c01 = ubx;
ubx_c02 = ubx;
ubx_c03 = ubx;

x_c01 = normalizeuv([x(1) x(2)], lbx_c01, ubx_c01);                       
x_c02 = normalizeuv([x(1) x(2)], lbx_c02, ubx_c02);
x_c03 = normalizeuv([x(1) x(2)], lbx_c03, ubx_c03);

t_c01 = x_c01*W1_c01;
t_c02 = x_c02*W1_c02;
t_c03 = x_c03*W1_c03;

c(1)= 1-srgt_evaluate(t_c01, srgt_c1);
c(2)= 1-srgt_evaluate(t_c02, srgt_c2);
c(3)= 1-srgt_evaluate(t_c03, srgt_c3);
ceq=[];

end