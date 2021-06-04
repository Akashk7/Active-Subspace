%clc;clear;close all;
%% Defining the variables
t0=[3.4405576 3.2799744];
Nmcs=1e6;
Pftarget=0.00135;
lbx = [0.1 0.1];
ubx = [10 10 ];
sd=[0.3 0.3];

%% PSF value at the reported optima by Method 2

[psf_c01,pf_c01,re_c01,beta_c01] = mcspsfconstraint1c(t0,lbx,ubx,sd,Nmcs,Pftarget);
[psf_c02,pf_c02,re_c02,beta_c02] = mcspsfconstraint2c(t0,lbx,ubx,sd,Nmcs,Pftarget);                                
[psf_c03,pf_c03,re_c03,beta_c03] = mcspsfconstraint3c(t0,lbx,ubx,sd,Nmcs,Pftarget);                                
