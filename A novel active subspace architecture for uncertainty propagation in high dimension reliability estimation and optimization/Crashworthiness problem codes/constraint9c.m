function y =  constraint9c(x)

x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
x4=x(:,4);
x6=x(:,5);
p1=x(:,6);
p3=x(:,7);

 y1 = 10.58;
 y2 = -0.674*(x1.*x2);
 y3 = -1.95*(x2.*p1);
 y4 = 0.02054*(x3.*p3);
 y5 = -0.0198*(x4.*p3);
 y6 = 0.028*(x6.*p3);
 y = (9.9)./(y1+y2+y3+y4+y5+y6);
end 