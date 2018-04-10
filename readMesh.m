function readMesh(mshfile)

global coordinates elements nn nel triNode

% ====================
% Open file msh file
% ====================
fileID = fopen(mshfile,'r');

% get first three lines and discard them
for i=1:3
    fgetl(fileID);
end
% Get physical names (mandatory)
tline = fgetl(fileID);
if ~strcmp(tline,'$PhysicalNames')
    error('Input data MUST declare PhysicalNames. Please check.');
end
% get number of names
nNames = str2double(fgetl(fileID));
% get names
phyNames = zeros(nNames,2);
% each row contains: physical-dimension physical-number
for i=1:nNames
    tline = fgetl(fileID);
    phyNames(i,:) = sscanf(tline,'%d %d %*s');
end
fgetl(fileID); % discard this line
% Read nodes
fgetl(fileID); % discard this line
tline = fgetl(fileID);
nn = str2double(tline);
%
coordinates = zeros(nn,3);
for i=1:nn
    tline = fgetl(fileID);
    coordinates(i,:) = sscanf(tline,'%*d %f %f %f');
end
fgetl(fileID); % discard this line
%
% read elements
fgetl(fileID); % discard this line
tline = fgetl(fileID);
nelT = str2double(tline);
elementsT = cell(nelT,9);
% count number of 3-node triangle
triCount = 0;
% count number of 4-node tetrahedra
tetraCount = 0;
for i=1:nelT
    tline = fgetl(fileID);
    C = str2double(strsplit(tline));
    switch C(2)
        case 2
            triCount = triCount + 1;
        case 4
            tetraCount = tetraCount + 1;
        otherwise
            error('Unknown element type. Please check.')
    end
    elementsT(i,:) = {C};
end
%close file
fclose(fileID);
% post-process data
nel = tetraCount;
elements = zeros(nel,5); % store number of physical entity, element tag
triNode = zeros(triCount,4);
%
% count number of 3-node triangle
triCount = 0;
% count number of 4-node tetrahedra
tetraCount = 0;
for i=1:nelT
    % get array
    v = elementsT{i};
    switch v(2)
        case 4
            tetraCount = tetraCount + 1;
            elements(tetraCount,1) = v(4);
            elements(tetraCount,2:5) = v(6:9);
        case 2
            triCount = triCount + 1;
            triNode(triCount,1) = v(4);
            triNode(triCount,2:4) = v(6:8);
        otherwise
            error('Unknown element type. Please check.')
    end
end
%
clear elementsT

fprintf('************************\n')
fprintf('Preparing data structures and printing initial mesh\n')
fprintf('************************\n')
fprintf('%s %d\n','Number of nodes     ......... ',nn)
fprintf('%s %d\n','Number of elements  ......... ',nel)
fflush(stdout);
%parseData(DBCSet, NBCSet, PFCSet, nodeSet, sideSet)


% % ====================
% % Open msh file
% % ====================
% fileID = fopen(mshfile,'r');
% 
% % get first three lines and discard them
% for i=1:3
%     fgetl(fileID);
% end
% 
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
%
% % read elements
% fgetl(fileID); % discard this line
% tline = fgetl(fileID);
% nelT = str2double(tline);
% elementsT = zeros(nelT,7);
% for i=1:nelT
%     tline = fgetl(fileID);
%     C = str2double(strsplit(tline));
%     elementsT(i,:) = C;
% end
% %close file
% fclose(fileID);
% % post-process data
% % count the number of elements in the domain (using material set)
% matKeys = cell2mat(keys(MAT));
% nmat = length(matKeys);
% elementCount = 0;
% for i=1:nmat
%     for j=1:nelT
%         if elementsT(j,4) == matKeys(i)
%             elementCount = elementCount + 1;
%         end
%     end
% end
% 
% nel = elementCount;
% elements = zeros(nel,4); % store number of physical entity, material number, etc
% 
% % I do not know a better way to filter out the elements from gmsh.
% % I will loop several times over the element list. Not sure if the logic
% % below will work in general cases. Please check.
% 
% % fill-in elements material number and conectivity
% elementCount = 0;
% for i=1:nmat
%     mkey = matKeys(i);
%     for j=1:nelT
%         matnum = elementsT(j,4); 
%         if matnum == mkey
%             elementCount = elementCount + 1;
%             elements(elementCount,1) = matnum;
%             elements(elementCount,3:4) = elementsT(j,6:7);
%         end
%     end
% end
% 
% % fill-in direction
% dirKeys = cell2mat(keys(BDSet));
% ndir = length(dirKeys);
% for i=1:ndir
%     dkey = dirKeys(i);
%     for j=1:nelT
%         n1 = elementsT(j,6:7);
%         dirnum = elementsT(j,4);
%         if dirnum == dkey
%             index = find( (elements(:,3) == n1(1) ) & (elements(:,4) == n1(2) ) );
%             elements(index,2) = dirnum; %#ok<*FNDSB>
%         end
%     end
% end
% 
% % count the number of side loads in the domain (using NBCSet set)
% if ~isempty(NBCSet)
%     sideloadKeys = cell2mat(keys(NBCSet));
%     nsl = length(sideloadKeys);
%     count = 0;
%     for i=1:nsl
%         for j=1:nelT
%             if elementsT(j,4) == sideloadKeys(i)
%                 count = count + 1;
%             end
%         end
%     end
%     
%     % fill-in sideLoad
%     sideLoad = zeros(count,4);
%     elementCount = 0;
%     for i=1:nsl
%         slkey = sideloadKeys(i);
%         valueSet = cell2mat(NBCSet(slkey));
%         for j=1:nelT
%             slnum = elementsT(j,4);
%             if slnum == slkey
%                 elementCount = elementCount + 1;
%                 sideLoad(elementCount,3) = valueSet(1);
%                 sideLoad(elementCount,1:2) = elementsT(j,6:7);
%                 sideLoad(elementCount,4) = valueSet(2);
%             end
%         end
%     end
% end
% %
% clearvars elementsT
end