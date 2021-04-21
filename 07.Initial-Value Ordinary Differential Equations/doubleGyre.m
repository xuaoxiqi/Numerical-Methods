clc; close all; clear all
graphicsON = 1;   % flag for graphics
outputDatFile = 1;
tstart = tic;
A = 0.1;    % parameters from Shadden 2005 Physica D
eps = 0.1;
omega = 2.0*pi/10.0;  % frequency of gyre oscillations

%% Part 1 - Initialize grid of particles through vector field
dx = 0.05;
xvec = 0:dx:2;
yvec = 0:dx:1;
[x0,y0] = meshgrid(xvec,yvec);  % grid of particles
yIC(1,:,:) = x0';
yIC(2,:,:) = y0';

if(graphicsON)
    subplot(2,1,1)
    dy = doublegyreVEC(0,yIC,A,eps,omega);
    quiver(yIC(1,1:4:end,1:4:end),yIC(2,1:4:end,1:4:end),dy(1,1:4:end,1:4:end),dy(2,1:4:end,1:4:end));
    axis([0 2 0 1]), drawnow
    subplot(2,1,2)
    % plot initial conditions
    plot(yIC(1,:),yIC(2,:),'r.','LineWidth',2,'MarkerSize',4)
    axis([0 2 0 1]), drawnow
end
    

%% Part 2 - Compute trajectory 
dt =1e-1;     % timestep
T = 10.0;     % duration of integration

yin(1,:,:) = x0';  % 0.5;
yin(2,:,:) = y0';  % 0.25;
fileID = fopen('final.dat','wt+');
for i=1:T/dt
    time = i*dt;
    
    if(outputDatFile)
        fprintf(fileID, '%18.15f  %18.15f  %18.15f\n',yin, time); 
    end
    if(graphicsON)
        subplot(2,1,1)
        dy = doublegyreVEC(time,yIC,A,eps,omega);
        quiver(yIC(1,1:4:end,1:4:end),yIC(2,1:4:end,1:4:end),dy(1,1:4:end,1:4:end),dy(2,1:4:end,1:4:end));
        axis([0 2 0 1])
        drawnow
    end
    
    yout = rk4singlestep(@(t,y)doublegyreVEC(t,y,A,eps,omega),dt,time,yin);
    yin = yout;
    if(graphicsON)        
        subplot(2,1,2)
        plot(yout(1,:),yout(2,:),'r.','LineWidth',2,'MarkerSize',4)
        axis([0 2 0 1])
        hold on; %drawnow
    end    

end
fclose(fileID);

% %% Part 3 -  Compute the finite-time Lyapunov exponent (sigma)
% myBlack = [0 0 0];
% myGray = [192 192 192];
% myWhite = [255 255 255];
% myRed = [255 0 0];
% myGold = [255 215 0];
% myBlue = [0 0 200];
% myTangerine = [255 69 0];
% mycolorpoint=[myWhite
%     myWhite
%     myWhite
%     myBlue
%     myBlue];
% mycolorposition=[1 17 33 49 64];
% mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:64,'linear','extrap');
% mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:64,'linear','extrap');
% mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:64,'linear','extrap');
% mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/256.0;
 
% reshape 3-dim array into 2-dim array
xT = reshape(yout(1,:,:),length(xvec),length(yvec));
yT = reshape(yout(2,:,:),length(xvec),length(yvec));

% Finite difference to compute the gradient
[dxTdx0,dxTdy0] = gradient(xT,dx,dx);
[dyTdx0,dyTdy0] = gradient(yT,dx,dx);

% compute sigma: large sigma indicates large mixing!
for i=1:length(xvec)
    for j=1:length(yvec)
        D(1,1) = dxTdx0(i,j);
        D(1,2) = dxTdy0(i,j);
        D(2,1) = dyTdx0(i,j);
        D(2,2) = dyTdy0(i,j);
        sigma(i,j) = (0.5/time)*log(max(eig(D'*D)));
        % sigma(i,j) = (1/T)*sqrt(max(eig(D'*D)));
    end
end

if(graphicsON)
    %figure
    contourf(x0',y0',sigma,'LineColor','none')
    set(gcf,'Position',[700 300 600 300])
    colorbar
    % caxis([0, 500]);
    %colormap(mycolor);
    drawnow;
end

%%
telapsed = toc(tstart);
disp(['Successfully: simulation completed! ',num2str(telapsed),' seconds'])

%%
function yout = rk4singlestep(fun,dt,t0,y0)

f1 = fun(t0,y0);
f2 = fun(t0+dt/2,y0+(dt/2)*f1);
f3 = fun(t0+dt/2,y0+(dt/2)*f2);
f4 = fun(t0+dt,y0+dt*f3);
yout = y0+(dt/6)*(f1+2*f2+2*f3+f4);
end

%%
function [dy] = doublegyreVEC(t,yin,A,eps,om)
x = yin(1,:,:);
y = yin(2,:,:);
% A =  0.1; eps = 0.25; om = 2*pi/10;  % from Shadden 2005 Physica D

u = zeros(size(x)); v = u;

a = eps * sin(om * t);
b = 1.0 - 2.0 * a;

f = a * x.^2.0 + b * x;
df = 2.0 * a * x + b;

u = -pi * A * sin(pi * f) .* cos(pi * y);
v =  pi * A * cos(pi * f) .* sin(pi * y) .* df;

dy = [u;v];
end