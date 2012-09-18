% Wrapper script for generating synthetic data with space dependent
% diffusion constant. M.L. 2012-09-17

clear

%% simulation parameters
L=2000;
R=400;
stepT=1/300;
locAccuracy=20;
N=5e5 % number of time points

savefile='nucleoidmodel1_singletrace.mat';

% nucleoid model 
Rx=600;
Ryz=300;
% effective diffusion constants (including position noise)
Dcyt_eff=3e6;   % nm^2/s
Dnuc_eff=0.6e6; % nm^2/s
% intrinsic diffusion constants 
Dcyt=Dcyt_eff-locAccuracy^2/stepT
Dnuc=Dnuc_eff-locAccuracy^2/stepT
Dfun=@(x,y,z)(Dcyt+(Dnuc-Dcyt)*((x-L/2).^2/Rx^2+(y.^2+z.^2)/Ryz^2<1));

pause(1)

%% simulate and save results!
[timePoints,Traj,m] = simCell_nonconst_D(L,R,Dfun,N,stepT,locAccuracy);
save(savefile)
%% plot 
figure(1)
clf
subplot(2,1,1)
hold on
[X,Y]=meshgrid(-400:5:2400,-400:5:400);
%set(surf(X,Y,0*Y,Dfun(X,Y,0)*1e-6),'edgecolor','none')
set(pcolor(X,Y,Dfun(X,Y,0)*1e-6),'edgecolor','none')
axis([-R L+R -R R])
box on
subplot(2,1,2)
hold on

clf
hold on
[X,Y,Z]=ellipsoid(1000,0,0,Rx,Ryz,Ryz);
set(surf(X,Y,Z),'facealpha',0.3,'edgecolor','r','facecolor','r')
plot3(Traj(:,1),Traj(:,2),Traj(:,3),'-k')
axis([-R L+R -R R])
axis equal
box on
view([-40 20])
