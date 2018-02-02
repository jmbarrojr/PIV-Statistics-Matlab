% PARALLEL
% Function to compute the all Single Point Statistics for velocity fields
% ca
%
% Author: Julio Barros - UIUC 2011
% version: 2.4
%
% v2.3: Added R2014a matlabpool commands
% v2.4: Added Insight binary vec files compatibility - 03/29/17

clear
close all
clc

%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%

%ext = input('Type, the extention of the files to be processed (eg: *.vec): ','s');
type = input('Type if it is 2D, Stereo, 3D: ','s');

if strcmp(type,'2D') == 1
    ext = '*.vec';
elseif strcmp(type,'Stereo') == 1
    ext = '*.v3d';
elseif strcmp(type,'3D') == 1
    ext = '*.v3v';
else
    error('You did not choose the correct option')
end

% Dealing with slash on all OS
if ispc == 1
    slash = '\';
else
    slash = '/';
end

pathDir = uigetdir('','Select the dir which has the files to be processed');
files = dir(strcat(pathDir,slash,ext));
if isempty(files) == 1
    ext = '*.dat';
    files = dir(strcat(pathDir,slash,ext));
    if isempty(files) == 1 && strcmp(type,'2D') == 1
        ext = '*.bvec';
        files = dir(strcat(pathDir,slash,ext));
    end 
end
ResultsFol = strcat(pathDir,slash,'Results',slash);
%%%%%%%%%%%%%%%%%%%%%%
%files = files(1:100); %DEBUG PURPOSE
%%%%%%%%%%%%%%%%%%%%%%

disp('Do you want to specify a Region of interest?');
ROI = input('y , n?: ','s');
if strcmp(ROI,'y')==1 || strcmp(ROI,'Y')==1
    %disp('Format: [left right bottom top]')
    %crop = input('Type the ROI: ');
    happyROI = 0;
    while happyROI == 0
        if strcmp(type,'2D') == 1
            [~,I,J,~,~,X,Y,U] = matrix([pathDir slash files(1).name]);
        else
            [~,I,J,~,~,X,Y,Z,U] = matrix([pathDir slash files(1).name]);
        end
        figure(1)
        subplot(1,2,1),contourf(U),axis equal,title('Original')
        rectROI = getrect;
        crop = round([rectROI(1) I-rectROI(3)-rectROI(1) rectROI(2) J-rectROI(4)-rectROI(2)]);
        [~,~,U] = WindowFile(crop,I,J,U);
        subplot(1,2,2),contourf(U),axis equal,title('Cropped')
        choice = questdlg('Happy with ROI?', ...
            'ROI Menu', ...
            'Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes'
                happyROI = 1;
            case 'No'
                happyROI = 0;
        end
    end
    close all
else
    crop = [1 1 1 1];
end

disp('Do you want to save the Fluctuation files?')
flucask = input('y or n?: ','s');

disp('Which direction do you want to caculate the line average?')
disp('1. Horizintal (across rows)')
disp('2. Vertical (across columns)')
direction = input(':? ');

disp('Do you want to apply any Spatial Calibration?')
calask = input('y or n?: ','s');
if strcmp(calask,'y') == 1
    dt = input('Input the dt in second [s]: ');
    mag = input('Input the magnification in [mm/px]: ');
    vel_factor = (mag/1000)/dt;
end

% Start Matlabpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool;
    cores = poolobj.NumWorkers;
else
    cores = poolobj.NumWorkers;
end

tstart_p = tic; % Start computing Time

% Split the list of the files into the number of cores
[CompVecList] = DistVecContent(cores,files);

disp('Calculating the Ensemble Average')
[AxisCell,Ucell,Lciavg,Omega_avg,CHCavg] = ensembleAverage(CompVecList,cores,crop,pathDir);

disp('')
disp('Calculating the Reynolds Stresses')
[uiujCell,lci,omega] = ReynoldsStress(CompVecList,cores,Ucell,Lciavg,Omega_avg,crop,pathDir,flucask);

% Determine the size of the final variables
[J,I] = size(Ucell{1});
N = length(Ucell);

% Applying the Calibration
if calask == 1
    % Axis
    for n=1:length(AxisCell)
        AxisCell{n} = AxisCell{n} .* mag;
    end
    % Velocities
    for n=1:length(Ucell)
        Ucell{n} = Ucell{n} .* vel_factor;
    end
    % Stresses
    for n=1:length(uiujCell)
        uiujCell{n} = uiujCell{n} .* (vel_factor.^2);
    end
end

% Putting all the variables in a Tecplot form
vel = mixingCell_Stat(I,J,AxisCell,Ucell,uiujCell,Lciavg,Omega_avg,lci,omega,CHCavg);
vel = dealNaN(vel);
vel = sortrows(vel, [2 , 1]);

if N == 2
    TecplotHeader = ['VARIABLES="x", "y",'...
        '"<math>a</math>U<math>q</math>", '...
        '"<math>a</math>V<math>q</math>", '...
        '"<math>a</math>u<sup>2</sup><math>q</math>", '...
        '"<math>a</math>uv<math>q</math>", '...
        '"<math>a</math>v<sup>2</sup><math>q</math>", '...
        '"<math>a</math>Lci<math>q</math>", '...
        '"<math>a</math>Omega<math>q</math>", '...
        '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
        '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
        '"CHC" ' ...
        'ZONE I=' num2str(I) ', J=' num2str(J) ', K=1, F=POINT'];
elseif N == 3
    TecplotHeader = ['VARIABLES="x", "y", "z",'...
        '"<math>a</math>U<math>q</math>", '...
        '"<math>a</math>V<math>q</math>", '...
        '"<math>a</math>W<math>q</math>", '...
        '"<math>a</math>u<sup>2</sup><math>q</math>", '...
        '"<math>a</math>uv<math>q</math>", '...
        '"<math>a</math>uw<math>q</math>", '...
        '"<math>a</math>v<sup>2</sup><math>q</math>", '...
        '"<math>a</math>vw<math>q</math>", '...
        '"<math>a</math>w<sup>2</sup><math>q</math>", '...
        '"<math>a</math>Lci<math>q</math>", '...
        '"<math>a</math>Omega<math>q</math>", '...
        '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
        '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
        '"CHC" ' ...
        'ZONE I=' num2str(I) ', J=' num2str(J) ', K=1, F=POINT'];
end

saver(ResultsFol,'EnsembleAverage_v2_3.dat',TecplotHeader,vel);

% Calculates the Line Average
if direction == 1
    if N == 2
        x = AxisCell{1}(1,:);
        U = nanmean(Ucell{1},direction);
        V = nanmean(Ucell{2},direction);
        uu = nanmean(uiujCell{1},direction);
        uv = nanmean(uiujCell{2},direction);
        vv = nanmean(uiujCell{3},direction);
        Lcil = nanmean(Lciavg,direction);
        Omegal = nanmean(Omega_avg,direction);
        lcil = nanmean(lci,direction);
        omegal = nanmean(omega,direction);
        tj = I;
        
        TecplotHeader2 = ['VARIABLES="x",'...
            '"<math>a</math>U<math>q</math>", '...
            '"<math>a</math>V<math>q</math>", '...
            '"<math>a</math>u<sup>2</sup><math>q</math>", '...
            '"<math>a</math>uv<math>q</math>", '...
            '"<math>a</math>v<sup>2</sup><math>q</math>", '...
            '"<math>a</math>Lci<math>q</math>", '...
            '"<math>a</math>Omega<math>q</math>", '...
            '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
            '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
            'ZONE I=1, J=' num2str(tj) ', K=1, F=POINT'];
        
        vel_l = [x' U' V' uu' uv' vv' Lcil' Omegal' lcil' omegal'];
        
    elseif N == 3
        
        x = AxisCell{1}(1,:);
        U = nanmean(Ucell{1},direction);
        V = nanmean(Ucell{2},direction);
        W = nanmean(Ucell{3},direction);
        uu = nanmean(uiujCell{1},direction);
        uv = nanmean(uiujCell{2},direction);
        uw = nanmean(uiujCell{3},direction);
        vv = nanmean(uiujCell{4},direction);
        vw = nanmean(uiujCell{5},direction);
        ww = nanmean(uiujCell{6},direction);
        Lcil = nanmean(Lciavg,direction);
        Omegal = nanmean(Omega_avg,direction);
        lcil = nanmean(lci,direction);
        omegal = nanmean(omega,direction);
        tj = I;
        
        TecplotHeader2 = ['VARIABLES="y", '...
            '"<math>a</math>U<math>q</math>", '...
            '"<math>a</math>V<math>q</math>", '...
            '"<math>a</math>W<math>q</math>", '...
            '"<math>a</math>u<sup>2</sup><math>q</math>", '...
            '"<math>a</math>uv<math>q</math>", '...
            '"<math>a</math>uw<math>q</math>", '...
            '"<math>a</math>v<sup>2</sup><math>q</math>", '...
            '"<math>a</math>vw<math>q</math>", '...
            '"<math>a</math>w<sup>2</sup><math>q</math>", '...
            '"<math>a</math>Lci<math>q</math>", '...
            '"<math>a</math>Omega<math>q</math>", '...
            '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
            '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
            'ZONE I=1, J=' num2str(J) ', K=1, F=POINT'];
        
        vel_l = [x' U' V' W' uu' uv' uw' vv' vw' ww' Lcil' Omegal' lcil' omegal'];
        
    end
    
    
elseif direction == 2
    
    if N == 2
        y = AxisCell{2}(:,1);
        U = nanmean(Ucell{1},direction);
        V = nanmean(Ucell{2},direction);
        uu = nanmean(uiujCell{1},direction);
        uv = nanmean(uiujCell{2},direction);
        vv = nanmean(uiujCell{3},direction);
        Lcil = nanmean(Lciavg,direction);
        Omegal = nanmean(Omega_avg,direction);
        lcil = nanmean(lci,direction);
        omegal = nanmean(omega,direction);
        tj = J;
        
        TecplotHeader2 = ['VARIABLES="y",'...
            '"<math>a</math>U<math>q</math>", '...
            '"<math>a</math>V<math>q</math>", '...
            '"<math>a</math>u<sup>2</sup><math>q</math>", '...
            '"<math>a</math>uv<math>q</math>", '...
            '"<math>a</math>v<sup>2</sup><math>q</math>", '...
            '"<math>a</math>Lci<math>q</math>", '...
            '"<math>a</math>Omega<math>q</math>", '...
            '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
            '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
            'ZONE I=1, J=' num2str(J) ', K=1, F=POINT'];
        
        vel_l = [y U V uu uv vv Lcil Omegal lcil omegal];
        
    elseif N == 3
        y = AxisCell{2}(:,1);
        U = nanmean(Ucell{1},direction);
        V = nanmean(Ucell{2},direction);
        W = nanmean(Ucell{3},direction);
        uu = nanmean(uiujCell{1},direction);
        uv = nanmean(uiujCell{2},direction);
        uw = nanmean(uiujCell{3},direction);
        vv = nanmean(uiujCell{4},direction);
        vw = nanmean(uiujCell{5},direction);
        ww = nanmean(uiujCell{6},direction);
        Lcil = nanmean(Lciavg,direction);
        Omegal = nanmean(Omega_avg,direction);
        lcil = nanmean(lci,direction);
        omegal = nanmean(omega,direction);
        tj = J;
        
        TecplotHeader2 = ['VARIABLES="y" '...
            '"<math>a</math>U<math>q</math>", '...
            '"<math>a</math>V<math>q</math>", '...
            '"<math>a</math>W<math>q</math>", '...
            '"<math>a</math>u<sup>2</sup><math>q</math>", '...
            '"<math>a</math>uv<math>q</math>", '...
            '"<math>a</math>uw<math>q</math>", '...
            '"<math>a</math>v<sup>2</sup><math>q</math>", '...
            '"<math>a</math>vw<math>q</math>", '...
            '"<math>a</math>w<sup>2</sup><math>q</math>", '...
            '"<math>a</math>Lci<math>q</math>", '...
            '"<math>a</math>Omega<math>q</math>", '...
            '"<math>a</math>lci<sub>rms</sub><math>q</math>", '...
            '"<math>a</math>omega<sub>rms</sub><math>q</math>", '...
            'ZONE I=1, J=' num2str(J) ', K=1, F=POINT'];
        
        vel_l = [y U V W uu uv uw vv vw ww Lcil Omegal lcil omegal];
    end
    
end
vel_l = dealNaN(vel_l);
vel_l = sortrows(vel_l,1);

% Save the file
saver(ResultsFol,'Line_Average_V2_3.dat',TecplotHeader2,vel_l);

tstop_p = toc(tstart_p);

disp('DONE')
disp([num2str(tstop_p/60) ' minutes taken'])

% Close Matlab Pool
%matlabpool close