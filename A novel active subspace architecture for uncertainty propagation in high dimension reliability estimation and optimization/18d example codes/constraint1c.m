function y =  constraint1c(x)
 
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
x5 = x(:,5);
x6 = x(:,6);
x7 = x(:,7);
x8 = x(:,8);
x9 = x(:,9);
x10 = x(:,10);
x11 = x(:,11);
x12 = x(:,12);
x13 = x(:,13);
x14 = x(:,14);
x15 = x(:,15);
x16 = x(:,16);
x17 = x(:,17);
x18 = x(:,18);

 y1 = 3*(x5.^2);
 y2 = x7.^2;
 y3 = 2*(x13.^2);
 y4 = (x1.*x10);
 y5 = 2*(x4.*x6);
 y6 = 3*(x2.*x9);
 y7 = 4*(x3.*x8);
 y8 = 5*(x11.*x12);
 y9 = 0.001*(x14+x15)-(x16.*x17)-x18;
 
 y = 210./(y1+y2+y3+y4+y5+y6+y7+y8+y9);
end