function y =  constraint2c(x)

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

 y1 = 2*(x2.^2);
 y2 = 4*(x3.^2);
 y3 = 5*(x1.*x10);
 y4 = 4*(x4.*x5);
 y5 = 3*(x6.*x7);
 y6 = 2*(x9.*x8);
 y7 = (x11.*x12);
 y8 = (x13.*x12);
 y9 = 0.001*(x14+x15)-x16-(2*(x17.^2))-(x18.^2);
 
 y = 200./(y1+y2+y3+y4+y5+y6+y7+y8+y9);
end