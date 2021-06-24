function x = denormalizeuv(xnorm,lbx,ubx)
x = (xnorm+1).*((ubx-lbx)./2)+lbx;
end
