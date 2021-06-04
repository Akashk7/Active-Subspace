function y =  constraint4c(x)

x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
x5=x(:,4);
x6=x(:,5);
x7=x(:,6);
p1=x(:,7);
p2=x(:,8);
p3=x(:,9);
 y1 = 28.98;
 y2 = 3.818*(x3);
 y3 = -4.2*(x1.*x2);
 y4 = 0.0207*(x5.*p3); 
 y5 = 6.63*x6.*p2;
 y6 = -7.7*x7.*p1;
 y7 = 0.32*p2.*p3;
y = 32./(y1+y2+y3+y4+y5+y6+y7);
end