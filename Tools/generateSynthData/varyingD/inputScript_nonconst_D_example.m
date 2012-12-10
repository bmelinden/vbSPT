% Wrapper script for generating synthetic data with space dependent
% diffusion constant. Will overwrite results by name 'savename'. 



%% simulation parameters
% Cell dimensions
L=2000;
R=400;
% Sampling interval
stepT=0.3e-3;
% trajectory lengths
trajLengths=round(5-10*log(1-rand(1,100)));
% Localization accuracy
locAccuracy=0;

savename = 'genData';

% an gaussian low-D region:
% Dfun=@(x,y,z)(1e6*(1-0.9*exp(-0.5*((x-1000).^2+y.^2+z.^2)/(300^2)))); % in nm^2/s

% two different domains:
% Dfun=@(x,y,z)(1e6*(1-0.9*(x>L/2)));

% simple nucleoid model
Dfun=@(x,y,z)(1e6*(1-0.9999*((x-L/2).^2/1e3^2+(y.^2+z.^2)/200^2<1)));

% linearly changing:
% Dfun=@(x,y,z)((1-0.99*(x/(L+R)))*1e6);


%% simulate!
for trajNr = 1:length(trajLengths)
    N=trajLengths(trajNr);
    [TimePoints, Traj, m] = simCell_nonconst_D(L, R, Dfun,  N, stepT,locAccuracy);
    data{trajNr} = Traj;
end
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
%plot(Traj(:,1),Traj(:,2),'-k')
for mm=1:length(trajLengths)
    %plot(data{mm}(:,1),data{mm}(:,2),'-k.')
    plot3(data{mm}(:,1),data{mm}(:,2),data{mm}(:,4),'-k.')
end
axis([-R L+R -R R])
box on
contour(X,Y,Dfun(X,Y,0)*1e-6)

% save generated data
save(savename);

% dx=diff(Traj(:,1));
% mean(dx(Traj(1:end-1,4)<1e6).^2/stepT/2)
% mean(dx(Traj(1:end-1,4)>1e6).^2/stepT/2)

