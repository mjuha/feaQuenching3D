function outfile = readData(filename)

% global variables
global DBCSet HFCSet NBCSet TS MAT isTimeDBC nodeSet NLOPT TTableData
global HTableData isNBCTempDependent

% global rhoDataAust rhoDataPer rhoDataMar 
% global cpDataAust cpDataPer tauData kappaDataAust kappaDataPer kappaDataMar

% Open file
fileID = fopen(filename,'r');

% get first three lines and discard them
for i=1:3
    fgetl(fileID);
end
% get mesh input file
tline = fgetl(fileID);
tmp = strsplit(tline);
len = length(tmp);
if len > 4
    % concatenate name
    for i=4:len-1
        mshfile = tmp{i};
        s1 = tmp{i+1};
        mshfile = strcat(mshfile,{' '},s1);
    end
else
    mshfile = tmp{4};
end
% get output file location
tline = fgetl(fileID);
tmp = strsplit(tline);
len = length(tmp);
if len > 3
    % concatenate name
    for i=3:len-1
        outfile = tmp{i};
        s1 = tmp{i+1};
        outfile = strcat(outfile,{' '},s1);
    end
else
    outfile = tmp{3};
end
% get material set association
tline = fgetl(fileID);
tmp = strsplit(tline);
nms = str2double(tmp(4)); % number of material sets to read
if nms ~= 1
    error('Only one material is permitted!')
end
matSet = cell(nms,2);
if nms == 0
    fgetl(fileID); % dummy line
else
    for i=1:nms
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        phyN = str2double(tmp(1)); % physical entity number
        matSet(i,:) = {phyN, tmp{2}}; 
    end
end
% get side set association
tline = fgetl(fileID);
tmp = strsplit(tline);
nss = str2double(tmp(4)); % number of side sets to read
sideSet = cell(nss,2);
if nss == 0
    fgetl(fileID); % dummy line
else
    for i=1:nss
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        phyN = str2double(tmp(1)); % physical entity number
        sideSet(i,:) = {phyN, tmp{2}};
    end
end
% get node set association
tline = fgetl(fileID);
tmp = strsplit(tline);
nns = str2double(tmp(4)); % number of node sets to read
nodeSet = cell(nns,2);
if nns == 0
    fgetl(fileID); % dummy line
else
    for i=1:nns
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        phyN = str2double(tmp(1)); % physical entity number
        nodeSet(i,:) = {phyN, tmp{2}};
    end
end
% read Dirichlet (temperature) BCs
tline = fgetl(fileID);
tmp = strsplit(tline);
ndbc = str2double(tmp(3)); % number of DBC to read
DBCSet = cell(ndbc,3);
if ndbc == 0
    fgetl(fileID); % dummy line
else
    for i=1:ndbc
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        name = tmp{4};
        dof = sscanf(tmp{7},'%[TFunction]');
        if strcmp(dof,'T')
            value = str2double(tmp(8)); % value to assign
            DBCSet(i,:) = {name, 'T' ,value};
            isTimeDBC = false;
        elseif strcmp(dof,'TFunction')
            f = tmp{8};            
            if strcmp(f,'Linear')
                val = 1;
                value = [ val, str2double(tmp(9)), str2double(tmp(10)) ];
                DBCSet(i,:) = {name, 'TFunction' ,value};
            elseif strcmp(f,'Table')
                val = 2;
                tableName = tmp{9};
                TTableData = csvread(tableName);
                value = [ val, 0.0 ];
                DBCSet(i,:) = {name, 'TFunction' ,value};
            else
                error('Function must be Linear or Table');
            end
%             value = [ val, str2double(tmp(9)), str2double(tmp(10)) ];
%             DBCSet(i,:) = {name, 'TFunction' ,value};
            isTimeDBC = true;
        else
            error('DOF must be T or TFunction, please check')
        end
    end
end
% read Neumann (convection) BCs
tline = fgetl(fileID);
tmp = strsplit(tline);
nnbc = str2double(tmp(3)); % number of NBC to read
NBCSet = cell(nnbc,3);
if nnbc == 0
    fgetl(fileID); %dummy line
else
    for i=1:nnbc
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        name = tmp{4};
        dof = sscanf(tmp{6},'%[Fconvection]');
        if strcmp(dof,'convection')
            % value to assign to film coefficient
            coeff = str2double(tmp(7)); % value to assign to film coefficient
            val = str2double(tmp(8)); % ambient temperature
            value = [ coeff, val ]; % ambient temperature
            NBCSet(i,:) = {name, 'convection', value};
            isNBCTempDependent = false;
        elseif strcmp(dof,'Fconvection')
            f = tmp{7};            
            if strcmp(f,'Table')
                tableName = tmp{8};
                val = str2double(tmp(9)); % ambient temperature
                value = [ 1, val ];
                HTableData = csvread(tableName);
                NBCSet(i,:) = {name, 'Fconvection' ,value};
            else
                error('Function must be Table');
            end
            isNBCTempDependent = true;
        else
            error('NBC must be convection or Fconvection, please check')
        end
%         coeff = str2double(tmp(7)); % value to assign to film coefficient
%         value = str2double(tmp(8)); % ambient temperature
%         NBCSet(i,:) = {name, coeff, value};
    end
end
% read heat flux BCs
tline = fgetl(fileID);
tmp = strsplit(tline);
nhfc = str2double(tmp(4)); % number of heat flux bcs to read
HFCSet = cell(nhfc,2);
if nhfc == 0
    fgetl(fileID); %dummy line
else
    for i=1:nhfc
        tline = fgetl(fileID);
        tmp = strsplit(tline);
        name = tmp{4};
        dof = sscanf(tmp{6},'%[heatflux]');
        if ~strcmp(dof,'heatflux')
            error('HPC must be heatflux, please check')
        end
        value = str2double(tmp(7)); % value to assign
        HFCSet(i,:) = {name, value};
    end
end

% read initial condition
% get next two lines and discard them
for i=1:2
    fgetl(fileID);
end
IC = cell(1,2); % initialize
tline = fgetl(fileID);
tmp = strsplit(tline);
name = tmp{2};
if strcmp(name,'Constant')
    % next line
    tline = fgetl(fileID);
    tmp = strsplit(tline);
    value = str2double(tmp(2)); % value to assign
    IC(:) = {name, value};
else
    error('Function name must be Constant');
end

for i=1:2
    fgetl(fileID);
end

% read time stepping
TS = cell(1,4);
tline = fgetl(fileID);
tmp = strsplit(tline);
dt = str2double(tmp(3)); % value to assign
tline = fgetl(fileID);
tmp = strsplit(tline);
ft = str2double(tmp(3)); % value to assign
tline = fgetl(fileID);
tmp = strsplit(tline);
name = tmp{2}; % value to assign
if strcmp(name,'Trapezoidal')
    alpha = 0.5;
elseif strcmp(name,'Backward')
    alpha = 1.0;
elseif strcmp(name,'Forward')
    alpha = 0.0;
else
    warning('Using Trapezoidal rule!')
    alpha = 0.5;
end
tline = fgetl(fileID);
tmp = strsplit(tline);
nst = str2double(tmp(4)); % value to assign
TS(:) = {dt, ft, alpha, nst};

% read nonlinear options
% get next 2 lines and discard them
for i=1:2
    fgetl(fileID);
end
tline = fgetl(fileID);
tmp = strsplit(tline);
name = tmp{2}; % value to assign
if strcmp(name,'Full')
    opt = 1;
elseif strcmp(name,'Modified')
    opt = 2;
else
    warning('Using Full Newton!')
    opt = 1;
end
tline = fgetl(fileID);
tmp = strsplit(tline);
name = tmp{2}; % value to assign
if strcmp(name,'AbsoluteNorm')
    optNorm = 1;
elseif strcmp(name,'RelativeNorm')
    error('RelativeNorm not implemented yet')
else
    warning('Using absolute norm!')
    optNorm = 1;
end
tline = fgetl(fileID);
tmp = strsplit(tline);
value = str2double(tmp(2)); % value to assign
tline = fgetl(fileID); % number of iteration
tmp = strsplit(tline);
iterNum = str2double(tmp(4)); % value to assign
NLOPT = [opt, optNorm, value, iterNum];

% get next line and discard it
fgetl(fileID);
% Read material properties
tline = fgetl(fileID);
tmp = strsplit(tline);
nmat = str2double(tmp(3)); % number of materials
if nmat ~= 1
    error('Only one material is permitted!')
end
%propKey = cell(1,nmat);
%nameKey = zeros(1,nmat);
MAT = cell(nms,10);
for i=1:nmat
    tline = fgetl(fileID);
    tmp = strsplit(tline);
    name = tmp{2};
    if ~strcmp(name,matSet{i,2})
      error('Material does not match!')
    end
    % data
    tline = fgetl(fileID);
    tmp = strsplit(tline);
    if strcmp(tmp(2),'S1080')
        % read rho austenite
        s = strcat(pwd,'\data\S1080\densityAustenite.csv');
        rhoDataAust =  csvread(s);
        % read rho perlite
        s = strcat(pwd,'\data\S1080\densityPerlite.csv');
        rhoDataPer =  csvread(s);
        % read cp austenite
        s = strcat(pwd,'\data\S1080\specificHeatAustenite.csv');
        cpDataAust =  csvread(s);
        % read cp perlite
        s = strcat(pwd,'\data\S1080\specificHeatPerlite.csv');
        cpDataPer =  csvread(s);
        % read thermal conductivity austenite
        s = strcat(pwd,'\data\S1080\thermalConductivityAustenite.csv');
        kappaDataAust =  csvread(s);
        % read thermal conductivity perlite
        s = strcat(pwd,'\data\S1080\thermalConductivityPerlite.csv');
        kappaDataPer =  csvread(s);
        % read thermal conductivity martensite
        s = strcat(pwd,'\data\S1080\thermalConductivityAustenite.csv');
        kappaDataMar =  csvread(s);
        % read temperature transformation data
        s = strcat(pwd,'\data\S1080\tempData.csv');
        tempData =  csvread(s);
        % read IT data
        s = strcat(pwd,'\data\S1080\tauS.csv');
        tauSData = csvread(s);
        s = strcat(pwd,'\data\S1080\tauF.csv');
        tauFData = csvread(s);
    else
        error('Only data for AISI 1080 is available!')
    end
    MAT(i,:) = {rhoDataAust, rhoDataPer, cpDataAust, cpDataPer, ...
        kappaDataAust, kappaDataPer, kappaDataMar, tempData, ...
        tauSData, tauFData};
%    propKey(i) = {prop1};
    %MAT = containers.Map(nameKey,propKey);
end

fclose(fileID);

% read mesh
readMesh(mshfile)

% Initialize data
initialCondition(IC)

% parse data
parseData(nodeSet, sideSet);

% % ====================
% % Open file msh file
% % ====================
% 
% fileID = fopen(mshfile,'r');
% 
% % get first three lines and discard them
% for i=1:3
%     fgetl(fileID);
% end
% % Get physical names (mandatory)
% tline = fgetl(fileID);
% if ~strcmp(tline,'$PhysicalNames')
%     error('Input data MUST declare PhysicalNames. Please check.');
% end
% % get number of names
% nNames = str2double(fgetl(fileID));
% % get names
% phyNames = zeros(nNames,2);
% % each row contains: physical-dimension physical-number
% for i=1:nNames
%     tline = fgetl(fileID);
%     phyNames(i,:) = sscanf(tline,'%d %d %*s');
% end
% fgetl(fileID); % discard this line
% % Read nodes
% fgetl(fileID); % discard this line
% tline = fgetl(fileID);
% nn = str2double(tline);
% %
% coordinates = zeros(nn,3);
% for i=1:nn
%     tline = fgetl(fileID);
%     coordinates(i,:) = sscanf(tline,'%*d %f %f %f');
% end
% fgetl(fileID); % discard this line
% %
% % read elements
% fgetl(fileID); % discard this line
% tline = fgetl(fileID);
% nelT = str2double(tline);
% elementsT = cell(nelT,7);
% % count number of 1-node point
% pointCount = 0;
% % count number of 2-node line
% lineCount = 0;
% % count number of 3-node triangle
% triCount = 0;
% for i=1:nelT
%     tline = fgetl(fileID);
%     C = str2double(strsplit(tline));
%     switch C(2)
%         case 15
%             pointCount = pointCount + 1;
%         case 1
%             lineCount = lineCount + 1;
%         case 2
%             triCount = triCount + 1;
%         otherwise
%             error('Unknown element type. Please check.')
%     end
%     elementsT(i,:) = {C};
% end
% %close file
% fclose(fileID);
% % post-process data
% nel = triCount;
% elements = zeros(nel,4); % store number of physical entity, element tag
% pointNode = zeros(pointCount,2);
% lineNode = zeros(lineCount,3);
% %
% % count number of 1-node point
% pointCount = 0;
% % count number of 2-node line
% lineCount = 0;
% % count number of 3-node triangle
% triCount = 0;
% for i=1:nelT
%     % get array
%     v = elementsT{i};
%     switch v(2)
%         case 15
%             pointCount = pointCount + 1;
%             pointNode(pointCount,1) = v(4);
%             pointNode(pointCount,2) = v(6);
%         case 1
%             lineCount = lineCount + 1;
%             lineNode(lineCount,1) = v(4);
%             lineNode(lineCount,2:3) = v(6:7);
%         case 2
%             triCount = triCount + 1;
%             elements(triCount,1) = v(4);
%             elements(triCount,2:4) = v(6:8);
%         otherwise
%             error('Unknown element type. Please check.')
%     end
% end
% %
% clearvars elementsT
% 
% fprintf('************************\n')
% fprintf('Preparing data structures and printing initial mesh\n')
% fprintf('************************\n')
% fprintf('%s %d\n','Number of nodes     ......... ',nn)
% fprintf('%s %d\n','Number of elements  ......... ',nel)
% fprintf('%s %d\n\n','Number of equations ......... ',neq)
% parseData(DBCSet, NBCSet, PFCSet, nodeSet, sideSet)
% 
% WriteVTKFile(outfile,0)

end