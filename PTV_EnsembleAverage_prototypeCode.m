clc
close all
clear all

if ispc == 1
    slash = '\';
else
    slash = '/';
end

pathDir = 'C:\Experiments10\Water Tunnel - JMB\SmoothWall\Analysis';

files = dir([pathDir slash '*LA.par']);
N = length(files);

for n=1:N
    
    A = importdata([pathDir slash files(n).name]);
    A.data = sortrows(A.data,[2,1]);
    X = round(A.data(:,1));
    Y = round(A.data(:,2));
    U = A.data(:,3);
    V = A.data(:,4);
    CHC = A.data(:,5);
        
    CHC(CHC<0) = 0;
    
    % Global Filter
    ind = U >= 0;
    CHC(ind) = 0;
    ind = V < -5 | V > 5;
    CHC(ind) = 0;
    ind = V < -3*std(V) | V > 3*std(V);
    CHC(ind) = 0;
    
    U = U.*CHC;
    V = V.*CHC;
    
    % Velocity correction
    U = -U;
    
    i = 1;
    for j=1:2048
        
        if j == 1 && n == 1
            Up = zeros(2048,1);
            count = zeros(2048,1);
        end
        
        while Y(i) == j
            %disp([num2str(j) ' & ' num2str(i)])
            if CHC(i) > 0
                Up(j) = U(i) + Up(j);
                count(j) = count(j) + 1;
            end
            if i < length(Y)
                i = i + 1;
            else
                break
            end
        end
    end
end

%Ensemble Average
Up = Up ./ count;

semilogx(Up)
