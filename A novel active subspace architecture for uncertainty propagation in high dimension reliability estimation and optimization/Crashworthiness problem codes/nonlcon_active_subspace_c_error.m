function [c,ceq]= nonlcon_active_subspace_c_error(x,lbx,ubx,W1_c01,W1_c02,W1_c03,W1_c04,W1_c05,W1_c06,W1_c07,W1_c08,W1_c09,W1_c10,...
srgt_c1,srgt_c2,srgt_c3,srgt_c4,srgt_c5,srgt_c6,srgt_c7,srgt_c8,srgt_c9,srgt_c10,srgt_evaluate,...
std_c1,std_c2,std_c3,std_c4,std_c5,std_c6,std_c7,std_c8,std_c9,std_c10)

lbx_c01=[lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c02=[lbx(1) lbx(2) lbx(3)];
lbx_c03=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(7)];
lbx_c04=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c05=[lbx(2) lbx(3) lbx(7)];
lbx_c06=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c07=[lbx(1) lbx(2) lbx(3) lbx(5) lbx(6) lbx(7)];
lbx_c08=[lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c09=[lbx(1) lbx(2) lbx(3) lbx(4) lbx(6)];
lbx_c10=[lbx(3) lbx(5) lbx(6) lbx(7)];

ubx_c01=[ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c02=[ubx(1) ubx(2) ubx(3)];
ubx_c03=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(7)];
ubx_c04=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c05=[ubx(2) ubx(3) ubx(7)];
ubx_c06=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c07=[ubx(1) ubx(2) ubx(3) ubx(5) ubx(6) ubx(7)];
ubx_c08=[ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c09=[ubx(1) ubx(2) ubx(3) ubx(4) ubx(6)];
ubx_c10=[ubx(3) ubx(5) ubx(6) ubx(7)];

x_c01=normalizeuv([x(2:4) x(6)], lbx_c01, ubx_c01);                       
x_c02=normalizeuv(x(1:3), lbx_c02, ubx_c02);
x_c03=normalizeuv([x(1:3) x(5) x(7)], lbx_c03, ubx_c03);
x_c04=normalizeuv([x(1:3) x(5:7)], lbx_c04, ubx_c04);
x_c05=normalizeuv([x(2:3) x(7)], lbx_c05, ubx_c05);
x_c06=normalizeuv([x(1:3) x(5:7)], lbx_c06, ubx_c06);
x_c07=normalizeuv([x(1:3) x(5:7)], lbx_c07, ubx_c07);
x_c08=normalizeuv([x(2:4) x(6)], lbx_c08, ubx_c08);
x_c09=normalizeuv([x(1:4) x(6)], lbx_c09, ubx_c09);
x_c10=normalizeuv([x(3) x(5:7)], lbx_c10, ubx_c10);

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

c(1)= (1-srgt_evaluate(t_c01, srgt_c1));
c(2)= (1-srgt_evaluate(t_c02, srgt_c2));
c(3)= (1-srgt_evaluate(t_c03, srgt_c3));
c(4)= (1-srgt_evaluate(t_c04, srgt_c4));
c(5)= (1-srgt_evaluate(t_c05, srgt_c5));
c(6)= (1-srgt_evaluate(t_c06, srgt_c6));
c(7)= (1-srgt_evaluate(t_c07, srgt_c7));
c(8)= (1-srgt_evaluate(t_c08, srgt_c8));
c(9)= (1-srgt_evaluate(t_c09, srgt_c9));
c(10)= (1-srgt_evaluate(t_c10, srgt_c10));
ceq=[];

end