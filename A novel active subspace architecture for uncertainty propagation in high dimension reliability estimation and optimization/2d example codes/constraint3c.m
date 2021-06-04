function y =  constraint3c(x)

x1 = x(:,1);
x2 = x(:,2);

y = 80./(((x1.^2)+(8*x2)+5));
end