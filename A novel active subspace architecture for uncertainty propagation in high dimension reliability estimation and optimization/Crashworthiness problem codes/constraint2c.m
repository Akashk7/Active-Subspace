function y =  constraint2c(x)
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
p1=x(:,4);
p3=x(:,5);
y1 = 46.36;
y2 = -9.9*(x2);
y3 = -12.9*(x1.*p1);
y4 = 0.1107*(x3.*p3);
y = 32./(y1+y2+y3+y4);
end