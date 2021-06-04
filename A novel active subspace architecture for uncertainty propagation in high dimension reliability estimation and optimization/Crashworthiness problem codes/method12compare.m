clc;clear all;close all;
%% Method1 and Method 2 comparison
t0=[0.511 1.417 0.500 1.352 0.658 1.473 0.500 0.345 0.192];
Nmcs=1e6;
Pftarget=0.0013;
lbx = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.192 0.192];
ubx = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 0.345 0.345];
for i=1:7
    sd(1,i)=(0.03);
end
sd(1,8) = 0.006;
sd(1,9) = 0.006;
mux10 = 0;
mux11 = 0;
sdx10=10;
sdx11=10;
mux = [mux10 mux11];
sdx = [sdx10 sdx11];
mu = [1 1 1 1 1 1 1 0.3 0.3];

R1 =ones(Nmcs,1)*1;
R2 =ones(Nmcs,1)*32;
R3 =ones(Nmcs,1)*32;
R4 =ones(Nmcs,1)*32;
R5 =ones(Nmcs,1)*0.32;
R6 =ones(Nmcs,1)*0.32;
R7 =ones(Nmcs,1)*0.32;
R8 =ones(Nmcs,1)*4;
R9 =ones(Nmcs,1)*9.9;
R10 =ones(Nmcs,1)*15.7;

%% Reliability values at the reported optima by Method 1

[Pf_c1_m1,Re_c1_m1]= mcsreconstraint1m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R1);
[Pf_c2_m1,Re_c2_m1]= mcsreconstraint2m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R2);
[Pf_c3_m1,Re_c3_m1]= mcsreconstraint3m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R3);
[Pf_c4_m1,Re_c4_m1]= mcsreconstraint4m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R4);
[Pf_c5_m1,Re_c5_m1]= mcsreconstraint5m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R5);
[Pf_c6_m1,Re_c6_m1]= mcsreconstraint6m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R6);
[Pf_c7_m1,Re_c7_m1]= mcsreconstraint7m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R7);
[Pf_c8_m1,Re_c8_m1]= mcsreconstraint8m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R8);
[Pf_c9_m1,Re_c9_m1]= mcsreconstraint9m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R9);
[Pf_c10_m1,Re_c10_m1]= mcsreconstraint10m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R10);

%% PSF value at the reported optima by Method 1

psf_c1_m1 = mcspsfconstraint1m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R1);
psf_c2_m1 = mcspsfconstraint2m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R2);                                
psf_c3_m1 = mcspsfconstraint3m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R3);                                
psf_c4_m1 = mcspsfconstraint4m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R4);
psf_c5_m1 = mcspsfconstraint5m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R5);
psf_c6_m1 = mcspsfconstraint6m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R6);
psf_c7_m1 = mcspsfconstraint7m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R7);
psf_c8_m1 = mcspsfconstraint8m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R8);
psf_c9_m1 = mcspsfconstraint9m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R9);
psf_c10_m1 = mcspsfconstraint10m(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R10);

%% Reliability values at the reported optima by Method 2

R=ones(Nmcs,1);
[Pf_c1_m2,Re_c1_m2]= mcsreconstraint1(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c2_m2,Re_c2_m2]= mcsreconstraint2(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c3_m2,Re_c3_m2]= mcsreconstraint3(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c4_m2,Re_c4_m2]= mcsreconstraint4(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c5_m2,Re_c5_m2]= mcsreconstraint5(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c6_m2,Re_c6_m2]= mcsreconstraint6(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c7_m2,Re_c7_m2]= mcsreconstraint7(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c8_m2,Re_c8_m2]= mcsreconstraint8(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c9_m2,Re_c9_m2]= mcsreconstraint9(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
[Pf_c10_m2,Re_c10_m2]= mcsreconstraint10(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);

%% PSF value at the reported optima by Method 2

psf_c1_m2 = mcspsfconstraint1(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c2_m2 = mcspsfconstraint2(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);                                
psf_c3_m2 = mcspsfconstraint3(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);                                
psf_c4_m2 = mcspsfconstraint4(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c5_m2 = mcspsfconstraint5(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c6_m2 = mcspsfconstraint6(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c7_m2 = mcspsfconstraint7(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c8_m2 = mcspsfconstraint8(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c9_m2 = mcspsfconstraint9(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);
psf_c10_m2 = mcspsfconstraint10(t0,lbx,ubx,sd,mux,sdx,Nmcs,Pftarget,R);

%% Comparison
Re=[Re_c1_m1,Re_c1_m2;Re_c2_m1,Re_c2_m2;Re_c3_m1,Re_c3_m2;Re_c4_m1,Re_c4_m2;Re_c5_m1,Re_c5_m2;Re_c6_m1,Re_c6_m2;Re_c7_m1,Re_c7_m2;Re_c8_m1,Re_c8_m2;Re_c9_m1,Re_c9_m2;Re_c10_m1,Re_c10_m2];
Pf=[Pf_c1_m1,Pf_c1_m2;Pf_c2_m1,Pf_c2_m2;Pf_c3_m1,Pf_c3_m2;Pf_c4_m1,Pf_c4_m2;Pf_c5_m1,Pf_c5_m2;Pf_c6_m1,Pf_c6_m2;Pf_c7_m1,Pf_c7_m2;Pf_c8_m1,Pf_c8_m2;Pf_c9_m1,Pf_c9_m2;Pf_c10_m1,Pf_c10_m2];
psf=[psf_c1_m1,psf_c1_m2;psf_c2_m1,psf_c2_m2;psf_c3_m1,psf_c3_m2;psf_c4_m1,psf_c4_m2;psf_c5_m1,psf_c5_m2;psf_c6_m1,psf_c6_m2;psf_c7_m1,psf_c7_m2;psf_c8_m1,psf_c8_m2;psf_c9_m1,psf_c9_m2;psf_c10_m1,psf_c10_m2];