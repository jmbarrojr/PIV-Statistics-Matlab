% Quadrant Analysis
% Function to compute the all Single Point Statistics
%
% Author: Julio Barros - UIUC 2014
% version: 1.2
% update for Matlab R2014a

clear
close all
clc

% Start Matlabpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool;
    cores = poolobj.NumWorkers;
else
    cores = poolobj.NumWorkers;
end

% Dealing with slash on all OS
if ispc == 1
    slash = '\';
else
    slash = '/';
end

%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%
type = input('Type if it is 2D, Stereo, Cross Plane, 3D: ','s');

pathDir = uigetdir('C:\','Select the dir which has the files to be processed');
%%%%%%%%
% DEBUG
% pathDir = 'I:\Experiments_Julio\WindTunnel_SmoothWall\30Hz_SW\Analysis\Fluctuation';
%%%%%%%%
ext = '*.dat';
files = dir(strcat(pathDir,slash,ext));
if isempty(files) == 1
    ext = '*.dat';
    files = dir(strcat(pathDir,slash,ext));
end
ResultsFol = strcat(pathDir,slash,'Results',slash);

disp('Do you want to apply any Spatial Calibration?')
calask = input('y or n?: ','s');
if strcmp(calask,'y') == 1
    dt = input('Input the dt in second [s]: ');
    mag = input('Input the magnification in [mm/px]: ');
    vel_factor = (mag/1000)/dt;
end

% Origin
X0 = input('X origin: ');
Y0 = input('Y origin: ');

%%%%%%%%%%%%%%%%%%%%%%
% files = files(1:5000); %DEBUG PURPOSE
%%%%%%%%%%%%%%%%%%%%%%
N = length(files);
%
% disp('Do you want to specify a Region of interest?');
% ROI = input('y , n?: ','s');
% if strcmp(ROI,'y')==1 || strcmp(ROI,'Y')==1
%     disp('Format: [left right bottom top]')
%     crop = input('Type the ROI: ');
% else
%     crop = [1 1 1 1];
% end

% tstart_p = tic; % Start computing Time

% Split the list of the files into the number of cores
[CompVecList] = DistVecContent(cores,files);

% Stack all the velocity into 3D matrix
spmd
    Np = length(CompVecList);
    for n = 1:Np
        
        vecfile = [pathDir slash CompVecList(n).name];
        if strcmp(type,'2D') == 1
            [~,I,J,~,~,X,Y,u,v,~,~,CHC] = matrix(vecfile);
            
        elseif strcmp(type,'Stereo') == 1
            [~,I,J,~,~,X,Y,Z,u,v,w,~,~,CHC] = matrix(vecfile);
            
        elseif strcmp(type,'Cross Plane') == 1
            [~,I,J,~,~,X,Y,Z,w,v,u,~,~,CHC] = matrix(vecfile);
        end
        
        if n == 1
            %u3d_p = zeros(J,I,Np);
            v3d_p = single(zeros(J,I,Np));
            u3d_p = single(zeros(J,I,Np));
        end
        
        u3d_p(:,:,n) = u;
        v3d_p(:,:,n) = v;
        
    end
end
% Gather data from all the cores
for c = 1:cores
    if c == 1;
        v3d = v3d_p{c};
    else
        v3d = cat(3,v3d,v3d_p{c});
    end
end
clear v3d_p
for c = 1:cores
    if c == 1;
        u3d = u3d_p{c};
    else
        u3d = cat(3,u3d,u3d_p{c});
    end
end
clear u3d_p

I = I{1};
J = J{1};
X = X{1};
Y = Y{1};
if strcmp(type,'Cross Plane') == 1 || strcmp(type,'Stereo') == 1
    Z = Z{1};
end
CHC = CHC{1};

u3d(CHC==0) = NaN;
v3d(CHC==0) = NaN;

%% Compute the Standard Deviation
sigma_u_2d = sqrt(nanmean(u3d.^2,3));
sigma_v_2d = sqrt(nanmean(v3d.^2,3));

sigma_u = sqrt(nanmean(nanmean(u3d.^2,3),2));
sigma_v = sqrt(nanmean(nanmean(v3d.^2,3),2));

%% Quadrant Analysis
H = input('What is the H you want?: ');
T_2d = H .* sigma_v_2d .* sigma_u_2d;

for i=1:I
    if i == 1
        T = zeros(J,I);
    end
    T(:,i) = H .* sigma_u .* sigma_v;
end

for n = 1:N
    if n == 1
        Q1 = zeros(J,I);
        Q2 = zeros(J,I);
        Q3 = zeros(J,I);
        Q4 = zeros(J,I);
        N1 = zeros(J,I);
        N2 = zeros(J,I);
        N3 = zeros(J,I);
        N4 = zeros(J,I);
    end
    
    %u = u3d(:,:,n);
    v = v3d(:,:,n);
    u = u3d(:,:,n);
    % Reynolds shear stress
    uv = u.*v;
    
    indexQ1 = u>0 & v>0 & abs(u.*v) >= T_2d;
    indexQ2 = u<0 & v>0 & abs(u.*v) >= T_2d;
    indexQ3 = u<0 & v<0 & abs(u.*v) >= T_2d;
    indexQ4 = u>0 & v<0 & abs(u.*v) >= T_2d;
    
    Q1(indexQ1) = Q1(indexQ1) + uv(indexQ1);
    Q2(indexQ2) = Q2(indexQ2) + uv(indexQ2);
    Q3(indexQ3) = Q3(indexQ3) + uv(indexQ3);
    Q4(indexQ4) = Q4(indexQ4) + uv(indexQ4);
    
    N1 = double(indexQ1) + N1;
    N2 = double(indexQ2) + N2;
    N3 = double(indexQ3) + N3;
    N4 = double(indexQ4) + N4;
    
end
% Ensemble Average
Q1 = Q1 ./ N;
Q2 = Q2 ./ N;
Q3 = Q3 ./ N;
Q4 = Q4 ./ N;

N1 = N1 ./ N;
N2 = N2 ./ N;
N3 = N3 ./ N;
N4 = N4 ./ N;

figure(1),
subplot(1,2,1),plot(Y(:,1),mean(Q1,2),'ro',Y(:,1),-mean(Q2,2),'bo')
subplot(1,2,2),plot(Y(:,1),mean(Q3,2),'ro',Y(:,1),-mean(Q4,2),'bo')

figure(2),
subplot(1,2,1),plot(Y(:,1),mean(N1,2),'ro',Y(:,1),mean(N2,2),'bo')
subplot(1,2,2),plot(Y(:,1),mean(N3,2),'ro',Y(:,1),mean(N4,2),'bo')

% Save the results
vel = mixing(I,J,X,Y,Q1,Q2,Q3,Q4,N1,N2,N3,N4,CHC);
vel = sortrows(vel,[2,1]);
TecplotHeader = ['VARIABLES="x", "y",'...
                 '"Q1", "Q2", "Q3", "Q4", '...
                 '"N1", "N2", "N3", "N4", '...
                 '"CHC", '...
                 'ZONE I=' num2str(I) ', J=' num2str(J) ', K=1, F=POINT'];

saver(ResultsFol,['QuadrantAnalysis_H=' num2str(H) '.dat'],TecplotHeader,vel);

% Save the Line Average
y = Y(:,1);
vel = [y nanmean(Q1,2) nanmean(Q2,2) nanmean(Q3,2) nanmean(Q4,2) nanmean(N1,2) nanmean(N2,2) nanmean(N3,2) nanmean(N4,2)];
TecplotHeader = ['VARIABLES="y", '...
                 '"Q1", "Q2", "Q3", "Q4", '...
                 '"N1", "N2", "N3", "N4", '...
                 'ZONE I=' num2str(1) ', J=' num2str(J) ', K=1, F=POINT'];
saver(ResultsFol,['LineAvg_QuadrantAnalysis_H=' num2str(H) '.dat'],TecplotHeader,vel);