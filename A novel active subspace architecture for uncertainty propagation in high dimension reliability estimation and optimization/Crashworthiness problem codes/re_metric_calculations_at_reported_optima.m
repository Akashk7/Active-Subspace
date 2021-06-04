clc;clear;close all;
%% Defining the variables
t0=[0.8008490 1.35 0.7133922 1.5 0.875 1.2 0.4];
Nmcs=1e6;
Pftarget=0.00135;
lbx = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
ubx = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
sd=[0.03 0.03 0.03 0.03 0.05 0.03 0.03];
mup1 = 0.345;
mup2 = 0.192;
mup3 = 0;
mup4 = 0;
sdp1 = 0.006;
sdp2 = 0.006;
sdp3=7;
sdp4=7;
mup = [mup1 mup2 mup3 mup4];
sdp = [sdp1 sdp2 sdp3 sdp4];

%% PSF value at the reported optima by Method 2

[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);                                
[psf_c03,pf_c03,re_c03,beta_c03] = mcspsfconstraint3c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);                                
[psf_c04,pf_c04,re_c04,beta_c04] = mcspsfconstraint4c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c05,pf_c05,re_c05,beta_c05] = mcspsfconstraint5c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c06,pf_c06,re_c06,beta_c06] = mcspsfconstraint6c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c07,pf_c07,re_c07,beta_c07] = mcspsfconstraint7c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c08,pf_c08,re_c08,beta_c08] = mcspsfconstraint8c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c09,pf_c09,re_c09,beta_c09] = mcspsfconstraint9c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
[psf_c10,pf_c10,re_c10,beta_c10] = mcspsfconstraint10c(t0,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget);
