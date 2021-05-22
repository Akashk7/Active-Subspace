function y =  constraint1c(x)
x1 = x(:,1);
x2 = x(:,2);

y = ((x1.^2).*x2)./20;
end