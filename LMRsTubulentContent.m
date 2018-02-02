% Quadrant Analysis for Barnacle data
% Function to compute the all Single Point Statistics
%
% Author: Julio Barros - UIUC 2014
% version: 1.1
% update for Matlab R2014a

clear all %#ok<CLSCR>
close all
clc

% Dealing with slash on all OS
if ispc == 1
    slash = '\';
else
    slash = '/';
end

% Start Matlabpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool;
    cores = poolobj.NumWorkers;
else
    cores = poolobj.NumWorkers;
end

%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%
type = input('Type if it is 2D, Stereo, Cross Plane, 3D: ','s');

pathDir = uigetdir('C:\','Select the dir which has the files to be processed');
ext = '*.vec';
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
else
    mag = 1;
    vel_factor = 1;
end

% Origin
X0 = input('X origin: ');
Y0 = input('Y origin: ');
%%%%%%%%%%%%%%%%%%%%%%
% files = files(1:500); %DEBUG PURPOSE
%%%%%%%%%%%%%%%%%%%%%%
N = length(files);

%% Split the list of the files into the number of cores
[CompVecList] = DistVecContent(cores,files);

for n = 1:400
        
        vecfile = [pathDir slash files(n).name];
        if strcmp(type,'2D') == 1
            [~,I,J,~,~,X,Y,U,V,CHC] = matrix(vecfile);
            
        elseif strcmp(type,'Stereo') == 1
            [~,I,J,~,~,X,Y,Z,U,V,W,CHC] = matrix(vecfile);
            
        elseif strcmp(type,'Cross Plane') == 1
            [~,I,J,~,~,X,Y,Z,W,V,U,CHC] = matrix(vecfile);
        end
        
        U = U .* vel_factor;
        
        U(CHC<=0) = NaN;
        if n == 1
            Uvec = U(:);
        else
            Uvec = [Uvec;U(:)];
        end
end

% Normalize the U vel with the spatial mean
Uvec_mean = nanmean(Uvec);
% Calculate the p.d.f
[count,uvec_bin] = hist(Uvec,256);
figure(10),bar(uvec_bin./Uvec_mean,count./max(count))
drawnow

Thr = inputdlg('Enter the Threshold:','Input',1);
Thr = str2double(Thr);
% Un-normalize the Threshold
Thr = Thr .* Uvec_mean;

%% Stack all the velocity into 3D matrix
spmd
    Np = length(CompVecList);
    for n = 1:Np
        
        vecfile = [pathDir slash CompVecList(n).name];
        flucvecfile = [pathDir slash 'Fluctuation' slash 'Fluc_' CompVecList(n).name(1:end-4) '.dat'];
        if strcmp(type,'2D') == 1
            [~,~,~,~,~,~,~,U,V,CHC] = matrix(vecfile);
            [~,I,J,~,~,X,Y,u,v,~,~,CHCf] = matrix(flucvecfile);
            
        elseif strcmp(type,'Stereo') == 1
            [~,I,J,~,~,X,Y,Z,U,V,W,CHC] = matrix(vecfile);
            
        elseif strcmp(type,'Cross Plane') == 1
            [~,I,J,~,~,X,Y,Z,W,V,U,CHC] = matrix(vecfile);
        end
        
        U = U .* vel_factor;
        u = u .* vel_factor;
        V = V .* vel_factor;
        v = v .* vel_factor;
        
        U(CHC<=0) = NaN;
        u(CHCf<=0) = NaN;
        V(CHC<=0) = NaN;
        v(CHCf<=0) = NaN;
        
        % LMR's
        Ilmr = double(U < Thr );
        
        if n == 1
            Um_p = U;
            Vm_p = V;
            uu_p = u.^2;
            vv_p = v.^2;
            tke_p = 0.5 .* (u.^2 + v.^2);
            uv_p = u.*v;
            
            Um_lmr_p = U .* Ilmr;
            Vm_lmr_p = V .* Ilmr;
            uu_lmr_p = u.^2 .* Ilmr;
            vv_lmr_p = v.^2 .* Ilmr;
            tke_lmr_p = 0.5 .* (u.^2 + v.^2) .* Ilmr;
            uv_lmr_p = u.*v .* Ilmr;
            
            Nlmr_p = Ilmr;
            
        else
            Um_p = U + Um_p;
            Vm_p = V + Vm_p;
            uu_p = u.^2 + uu_p;
            vv_p = v.^2 + vv_p;
            tke_p = 0.5 .* (u.^2 + v.^2) + tke_p;
            uv_p = u .* v + uv_p;
            
            Um_lmr_p = U .* Ilmr + Um_lmr_p;
            Vm_lmr_p = V .* Ilmr + Vm_lmr_p;
            uu_lmr_p = (u.^2 .* Ilmr) + uu_lmr_p;
            vv_lmr_p = (v.^2 .* Ilmr) + vv_lmr_p;
            tke_lmr_p = (0.5 .* (u.^2 + v.^2) .*Ilmr )+ tke_lmr_p;
            uv_lmr_p = (u.*v .* Ilmr) + uv_lmr_p;
            
            Nlmr_p = Nlmr_p + Ilmr;
        end
            
    end
end
% Gather data from all the cores
for c = 1:cores
    if c == 1;
        Um = Um_p{c};
        Vm = Vm_p{c};
        uu = uu_p{c};
        vv = vv_p{c};
        tke = tke_p{c};
        uv = uv_p{c};
        
        Um_lmr = Um_lmr_p{c};
        Vm_lmr = Vm_lmr_p{c};
        uu_lmr = uu_lmr_p{c};
        vv_lmr = vv_lmr_p{c};
        tke_lmr = tke_lmr_p{c};
        uv_lmr = uv_lmr_p{c};
        
        Nlmr = Nlmr_p{c};
    else
        Um = Um + Um_p{c};
        Vm = Vm + Vm_p{c};
        uu = uu + uu_p{c};
        vv = vv + vv_p{c};
        tke = tke + tke_p{c};
        uv = uv + uv_p{c};
        
        Um_lmr = Um_lmr + Um_lmr_p{c};
        Vm_lmr = Vm_lmr + Vm_lmr_p{c};
        uu_lmr = uu_lmr + uu_lmr_p{c};
        vv_lmr = vv_lmr + vv_lmr_p{c};
        tke_lmr = tke_lmr + tke_lmr_p{c};
        uv_lmr = uv_lmr + uv_lmr_p{c};
        
        Nlmr = Nlmr + Nlmr_p{c};
    end
end
% Ensemble Average
Um = Um ./ N;
Vm = Vm ./ N;
uu = uu ./ N;
vv = vv ./ N;
tke = tke ./ N;
uv = uv ./ N;

Um_lmr = Um_lmr ./ N;
Vm_lmr = Vm_lmr ./ N;
uu_lmr = uu_lmr ./ N;
vv_lmr = vv_lmr  ./ N;
tke_lmr = tke_lmr  ./ N;
uv_lmr = uv_lmr  ./ N;

I = I{1};
J = J{1};
X = X{1};
Y = Y{1};
% Axis correction
Y = Y - Y0;
X = X - X0;
X = X .* mag;
Y = Y .* mag;
if strcmp(type,'Cross Plane') == 1 || strcmp(type,'Stereo') == 1
    Z = Z{1} .* mag;
end
CHC = CHC{1};
CHCf = CHCf{1};

Nlmr(CHC<=0) = NaN;

x = X(1,:);
y = Y(:,1);
figure(1), colormap jet
subplot(1,2,1),imagesc(x,y,Um),axis equal tight
set(gca,'Ydir','Normal')
subplot(1,2,2),imagesc(x,y,Um_lmr),axis equal tight
set(gca,'Ydir','Normal')

%% SAVE AVERAGE
Data = mixing(I,J,X,Y,Um_lmr./Um, Vm_lmr./Vm, uu_lmr./uu, uv_lmr./uv, vv_lmr./vv, Nlmr./N,CHCf);
Data = sortrows(Data,[2,1]);
Data = dealNaN(Data);

TecplotHeader = ['VARIABLES="x [mm]", "y [mm]",'...
                 '"<math>a</math>U<math>q</math><sub>lmr</sub>", '...
                 '"<math>a</math>V<math>q</math><sub>lmr</sub>", '...
                 '"<math>a</math>u''<sup>2</sup><math>q</math><sub>lmr</sub>", '...
                 '"<math>a</math>u''v''<math>q</math><sub>lmr</sub>", '...
                 '"<math>a</math>v''<sup>2</sup><math>q</math><sub>lmr</sub>", '...
                 '"Nlmr", "CHC", '...
                 'ZONE I=' num2str(I) ', J=' num2str(J) ', K=1, F=POINT'];

saver([pathDir slash 'Results' slash],'LMR_Content.dat',TecplotHeader,Data);