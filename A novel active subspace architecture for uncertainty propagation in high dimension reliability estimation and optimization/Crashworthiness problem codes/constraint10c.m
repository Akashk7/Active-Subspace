function y =  constraint10c(x)

x3=x(:,1);
x5=x(:,2);
x6=x(:,3);
x7=x(:,4);
p2=x(:,5);
p3=x(:,6);
p4=x(:,7);

 y1 = 16.45;
 y2 = -0.489*(x3.*x7);
 y3 = -0.843*(x5.*x6);
 y4 = 0.0432*(p2.*p3);
 y5 = -0.0556*(p2.*p4);
 y6 = -0.000786*(p4.^2);
y = (15.69)./(y1+y2+y3+y4+y5+y6);
end 