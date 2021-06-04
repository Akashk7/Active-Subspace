clc;clear;close all;
%% Defining the variables
t0=[2.12 3.39 2.00 3.75 2.40 2.00 2.46 4.00 3.66 4.00 4.00 2.00 3.52 3.47 3.29 3.27 2.00 2.00];
Nmcs=1e6;
Pftarget=0.0013;
lbx = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
ubx = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];
sd  = [0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];

%% PSF value at the reported optima by Method 2

[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2c(t0,lbx,ubx,sd,Nmcs,Pftarget);                                
[psf_c03,pf_c03,re_c03,beta_c03] = mcspsfconstraint3c(t0,lbx,ubx,sd,Nmcs,Pftarget);                                
[psf_c04,pf_c04,re_c04,beta_c04] = mcspsfconstraint4c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c05,pf_c05,re_c05,beta_c05] = mcspsfconstraint5c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c06,pf_c06,re_c06,beta_c06] = mcspsfconstraint6c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c07,pf_c07,re_c07,beta_c07] = mcspsfconstraint7c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c08,pf_c08,re_c08,beta_c08] = mcspsfconstraint8c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c09,pf_c09,re_c09,beta_c09] = mcspsfconstraint9c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c10,pf_c10,re_c10,beta_c10] = mcspsfconstraint10c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c11,pf_c11,re_c11,beta_c11] = mcspsfconstraint11c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c12,pf_c12,re_c12,beta_c12] = mcspsfconstraint12c(t0,lbx,ubx,sd,Nmcs,Pftarget);

f = Weightc(t0);