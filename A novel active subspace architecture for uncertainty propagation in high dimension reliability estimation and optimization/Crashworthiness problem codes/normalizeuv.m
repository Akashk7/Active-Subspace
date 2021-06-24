function Xnorm = normalizeuv(x,lbx,ubx)
Xnorm = ((x-lbx)./((ubx-lbx)./2))-1;
end