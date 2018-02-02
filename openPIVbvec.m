function [NVar,I,J,K,VARIABLES,InsightHeader] = openPIVbvec(file)
%OPENPIVBVEC
% This function opens the new BINARY vector files on the new Insight V11.1
% 
% Output variables sequence is similar to MATRIX function for usability
% Author: Julio Barros
%         USNA - 2017
%
% VERSION 1.0
%
% See also: matrix

% Main function was prototyped for readbility
[NVar,I,J,K,VARIABLES,InsightHeader] = readTecplotBinary(file);

end

%%%%%%%%%%%%%%%
% Main Function
%%%%%%%%%%%%%%%
function [NVar,I,J,K,VARIABLES,InsightHeader] = readTecplotBinary(file)
fid = fopen(file,'r');
%[filename, permission, machineformat, encoding] = fopen(fid);

%%%%%%%%%%%%%%%%%%%
% I . Header Section
% Header section i. magic number 
% Header section starts with Magic number (i.e. Version number)
%
%  +-----------+
%  | "#!TDV112"|       8 Bytes, exact characters "#!TDV112".
%  +-----------+       Version number follows the "V" and
%                             consumes the next 3 characters (for
%                             example: "V75 ", "V101").
magic_number = fread(fid,8,'*char')';

% ii. Integer value of 1.
%
%  +-----------+
%  | INT32     |       This is used to determine the byte order
%  +-----------+       of the reader relative to the writer.
%
%
byteOrder = fread(fid,1,'*int32');

% iii. Title and variable names.
%
% 	      +-----------+ 
% 	      | INT32     |      FileType: 0 = FULL, 
% 	      +-----------+                1 = GRID, 
% 				                       2 = SOLUTION 
FileType = fread(fid,1,'*int32');

%         +-----------+
%         | INT32*N   |       The TITLE. (See note 1.)
%         +-----------+
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    Title(1,n) = char(ch);
    n = n + 1;
end

%         +-----------+
%         | INT32     |       Number of variables (NumVar) in the datafile.
%         +-----------+
NVar = fread(fid,1,'*int32');

%         +-----------+
%         | INT32*N   |       Variable names.  N =  L[1] + L[2] + .... L[NumVar]
%         +-----------+       where:
%                                    L[i] = length of the ith variable name + 1
%                                           (for the terminating 0 value).
%                                           (See note 1.)
for nv = 1:NVar
    n = 1;
    ch = 1;
    while ch~=0
        ch = fread(fid,1,'*int32');
        Var(1,n) = char(ch);
        n = n + 1;
    end
    
    VARIABLES{nv,1} = Var; %#ok<*AGROW>
end

% iv. Zones
%         +-----------+
%         | FLOAT32   |       Zone marker. Value = 299.0
%         +-----------+
ZoneMarker = fread(fid,1,'float32');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Various TSI Insight Header Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aplication
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{1,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32');
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{1,2}(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
fread(fid,1,'*int32'); % Null

% ImageWidth
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{2,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32'); % Null
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    ImageWidth(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
ImageWidth = str2double(ImageWidth);
InsightHeader{2,2} = ImageWidth;
fread(fid,1,'*int32'); % Null

% ImageHeight
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{3,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32'); % Null
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    ImageHeight(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
ImageHeight = str2double(ImageHeight);
InsightHeader{3,2} = ImageHeight;
fread(fid,1,'*int32'); % Null

% LengthUnit
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{4,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32');
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{4,2}(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
fread(fid,1,'*int32');

% OriginImageX
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{5,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32');
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    OriginImageX(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
OriginImageX = str2double(OriginImageX);
InsightHeader{5,2} = OriginImageX;
fread(fid,1,'*int32');

% OriginImageY
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{6,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32');
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    OriginImageY(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
OriginImageY = str2double(OriginImageY);
InsightHeader{6,2} = OriginImageY;
fread(fid,1,'*int32');

% MicrosecondsPerDeltaT
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{7,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32');
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    dt(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
dt = str2double(dt);
InsightHeader{7,2} = dt;
fread(fid,1,'*int32');

% TimeUnit
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{8,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32'); % NULL
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{8,2}(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
%disp(TimeUnit)
fread(fid,1,'*int32'); % NULL

% SecondaryPeakNumber
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{9,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32'); % NULL
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    SecondaryPeakNumber(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
SecondaryPeakNumber = str2double(SecondaryPeakNumber);
InsightHeader{9,2} = SecondaryPeakNumber;
fread(fid,1,'*int32'); % NULL

% DewarpImageSource
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    InsightHeader{10,1}(1,n) = char(ch);
    n = n + 1;
end
fread(fid,1,'*int32'); % NULL
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    Dewarp(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end
Dewarp = str2double(Dewarp);
InsightHeader{10,2} = Dewarp;
fread(fid,1,'*int32'); % NULL

%         +-----------+
%         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
%         +-----------+
% Zone Name
n = 1;
ch = 1;
while ch~=0
    ch = fread(fid,1,'*int32');
    ZoneName(1,n) = char(ch); %#ok<*AGROW>
    n = n + 1;
end

%         +-----------+
%         | INT32     |       ParentZone: Zero based zone number within this datafile
%         +-----------+                   to which this zone is a child.
ParentZone = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |       StrandID: -2 = pending strand ID for assignment by Tecplot,
%         +-----------+                 -1 = static strand ID
%                                        0 <= N < 32700 valid strand ID
StrandID = fread(fid,1,'*int32');

%         +-----------+
%         | FLOAT64   |       Solution time.
%         +-----------+
SolutionTime = fread(fid,1,'*float64');

%         +-----------+
%         | INT32     |       Zone Color (set to -1 if you want tecplot to
%         +-----------+       determine).
ZoneColor = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |       ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,
%         +-----------+                3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
%                                      6=POLYGON,7=POLYHEDRON
ZoneType = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |       DataPacking 0=Block, 1=Point
%         +-----------+
DataPacking = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |       Specify Var Location.  0 = Don't specify, all data is 
%         +-----------+       located at the nodes.  1 = Specify
%         if "specify var location" == 1
%           +-----------+
%           | INT32*NV  |     Variable Location (only specify if above is 1).  
%           +-----------+     0 = Node, 1 = Cell Centered (See note 5.)
VarLocation = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |     Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE)
%         +-----------+     These raw values are a compact form of the local 1-to-1 face
%                           neighbors are fully specified and therefore it will not
%                           perform auto face neighbor assignment, thereby improving
%                           Tecplot's time to first plot.
%                           See data section below for format details. ORDERED and
%                           FELINESEG zones must specify 0 for this value since raw
%                           face neighbors are not defined for these zone types.
RawLocal = fread(fid,1,'*int32');

%         +-----------+
%         | INT32     |     Number of miscellaneous user defined face neighbor connections 
%         +-----------+     (value >= 0) This value is in addition to the face neighbors
%                             supplied in the raw section.
%
%         if "number of miscellaneous user defined face neighbor connections" != 0
%           +-----------+
%           | INT32     |     User defined face neighbor mode
%           +-----------+     (0=Local 1-to-1, 1=Local 1-to-many, 2=Global 1-to-1, 
%                             3=Global 1-to-many)
%           if FE Zone:
%             +-----------+
%             | INT32     |     Indicates if the finite element face neighbors are
%             +-----------+     completely specified by the miscellaneous face neighbors
%                               given: (0=NO, 1=YES). If yes, then Tecplot will not perform
%                               auto assignment of face neighbors otherwise all faces not
%                               specified are considered boundaries. If no, then Tecplot will
%                               perform auto-assignment of the face neighbors unless the
%                               raw face neighbor array was supplied. This option is not
%                               valid for ORDERED zones.
%         if Ordered Zone:
%           +-----------+
%           | INT32*3   |     IMax,JMax,KMax
%           +-----------+
IJK = fread(fid,1*3,'*int32');
I = IJK(1);
J = IJK(2);
K = IJK(3);

%         For all zone types (repeat for each Auxiliary data name/value pair):
%         +-----------+
%         | INT32     |       1=Auxiliary name/value pair to follow
%         +-----------+       0=No more Auxiliar name/value pairs.
%         If the above is 1, then supply the following:
%           +-----------+
%           | INT32*N   |     name string (See note 1.)
%           +-----------+      
%           +-----------+
%           | INT32     |     Auxiliary Value Format (Currently only allow 0=AuxDataType_String)
%           +-----------+
%           +-----------+
%           | INT32*N   |     value string  (See note 1.)
%           +-----------+      
Auxiliary = fread(fid,1,'*int32');

EOHMARKER = fread(fid,1,'*float32');

% II.  DATA SECTION (don?t forget to separate the header from the data
%      with an EOHMARKER).  The data section contains all of the data
%      associated with the zone definitions in the header.
%   i. For both ordered and fe zones:
%           +-----------+
%           | FLOAT32   |   Zone marker  Value = 299.0
%           +-----------+
ZoneMarker = fread(fid,1,'*float32');

%           +-----------+
%           | INT32*N   | Variable data format, N=Total number of vars
%           +-----------+ 1=Float, 2=Double, 3=LongInt,
%                         4=ShortInt, 5=Byte, 6=Bit
VarDataFormat = fread(fid,NVar,'*int32');

%           +-----------+
%           | INT32     | Has passive variables: 0 = no, 1 = yes.
%           +-----------+
PassiveVar = fread(fid,1,'*int32');

%           +-----------+
%           | INT32     | Has variables sharing: 0 = no, 1 = yes.
%           +-----------+
SharedVar = fread(fid,1,'*int32');

%           +-----------+
%           | INT32     | Zero based zone number to share connectivity
%           +-----------+ list with (-1 = no sharing). FEPOLYGON and
%                         FEPOLYHEDRON zones use this zone number to
%                         share face map data.
ShareConnectivity = fread(fid,1,'*int32');

%Compressed list of min/max pairs for each non-shared and non-passive
%variable. For each non-shared and non-passive variable (as specified
%above):
%           +-----------+
%           | FLOAT64   | Min value
%           +-----------+
%           +-----------+
%           | FLOAT64   | Max value
%           +-----------+
for nv = 1 : NVar
    MinVal{nv} = fread(fid,1,'*float64'); %#ok<*AGROW>
    MaxVal{nv} = fread(fid,1,'*float64'); %#ok<*AGROW>
end

%           +-----------+
%           | xxxxxxxxxx| Zone Data. Each variable is in data format as
%           +-----------+ specified above.
Data = fread(fid,I*J*K*NVar,'*float');
for nv = 1 : NVar
    if nv == 1
        VARIABLES{nv,2} = reshape(Data(1:I*J*K),I,J,K)';
    else
        VARIABLES{nv,2} = reshape(Data(I*J*K*(nv-1)+1:I*J*K*nv),I,J,K)';
    end
end

fclose(fid);

end