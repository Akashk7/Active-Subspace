function uc = plotasmetamodel(W1,tc12dmin,tc12drange,tc12dmax,srgtSRGTRBF,lbl)
Nsamples=1e4;
maxcount=1e5;
[t1,~]= zonotope_vertices(W1,Nsamples,maxcount);
tc1 = convhull(t1);
t = t1(tc1,:);
t1c1s = tc12dmin(1):tc12drange(1)/99:tc12dmax(1);
t2c1s = tc12dmin(2):tc12drange(2)/99:tc12dmax(2);
[u1,u2]=meshgrid(t1c1s,t2c1s);
U=[u1(:),u2(:)];
in = inhull(U,t1);
U1 = U(in,:);
psf1hat  = srgtsRBFEvaluate(U1, srgtSRGTRBF);
psf1corner  = srgtsRBFEvaluate(t, srgtSRGTRBF);
psf1hatnorm = norm01(psf1hat);
psf1cornernorm = norm01(psf1corner);
DT = delaunay(U1(:,1),U1(:,2));
fs =12; 
uc = figure;
trisurf(DT,U1(:,1),U1(:,2),psf1hatnorm,'EdgeColor','interp');
hold on
scatter3(t1(tc1,1),t1(tc1,2),psf1cornernorm,'r','filled','MarkerEdgeColor','k');
plot3(t1(tc1,1),t1(tc1,2),psf1cornernorm,'r-','LineWidth',2);
view([0 90]);
c = colorbar('FontWeight','bold','FontSize',fs);
c.Label.String=lbl;
% caxis([0 1]);
legend({'RBF Metamodel','Zonotope Vertices','ConvexHull'},'Location','SouthWest','FontWeight','bold','FontSize',fs);
grid on
xlabel('Active Variable 1','FontWeight','bold','FontSize',fs);
ylabel('Active Variable 2','FontWeight','bold','FontSize',fs);
set(gca,'FontWeight','bold','FontSize',fs);




