function y =  constraint2c(x)
x1 = x(:,1);
x2 = x(:,2);

y1 = ((x1 + x2 -5).^2)./30;
y2 = ((x1 - x2 - 12).^2)./120; 

y = (y1+y2);
end