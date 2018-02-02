clear all
close all
clc

if ispc == 1
    slash = '\';
else
    slash = '/';
end

pathDir = 'V:\Experiments10\Water Tunnel - JMB\Smooth_ZPG_adjusted_29Mp\Analysis\Results';
%pathDir = uigetdir();

fileName = 'EnsembleAverage_v2_3.dat';

[nc,I,J,dx,dy,X,Y,U,V,uu,uv,vv] = matrix([pathDir slash fileName]);

% Velocity conversion
U = -U;
uv = -uv;

% Apply the calibration
cal = 11.98e-3; % [mm/px]
dt = 200e-6; % [s]

X = X .* cal; % [mm]
Y = Y .* cal; % [mm]
U = U .* (cal/1000/dt); % [m/s]
u = sqrt(uu) .* (cal/1000/dt).^2; % [m/s]

% Choice Matrix
CHC = U == 0;
U(CHC) = NaN;
V(CHC) = NaN;
u(CHC) = NaN;

% Velocity profiles
y = Y(:,1);
Ul = nanmean(U,2);

figure(1)
plot(y,Ul,'ko',y,U(:,I/2),'rx')