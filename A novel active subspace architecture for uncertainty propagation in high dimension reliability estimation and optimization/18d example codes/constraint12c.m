function y =  constraint12c(x)

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
    
y1 = 0.001*(x1);
y2 = 0.001*(x2);
y3 = 0.001*(x3);
y4 = 0.001*(x4);
y5 = 0.001*(x5);
y6 = 0.001*(x6);
y7 = 0.001*(x7);
y8 = 0.001*(x8);
y9 = 0.001*(x9);
y10= 0.001*(x10);
y11= 0.001*(x11);
y12= 0.001*(x12);
y13= 0.001*(x13);
y14= -4*((x14+x15-5).^2);
y15= -((x14-x15-12).^2);

y = (y1+y2+y3+y4+y5+y6+y7+y8+y9+y10+y11+y12+y13+y14+y15)./(-120);

end 