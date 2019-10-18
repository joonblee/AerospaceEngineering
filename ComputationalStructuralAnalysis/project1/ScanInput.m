function ScanInput(ifid)
%% Scan input file for Finite element analysis

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

%--- V A R I A B L E   D E F I N I T I O N---------------------------------
%   'title': title of project
%   'npoin': number of nodes
%   'nelem': number of elements
%   'nnode': number of nodes per element
%   'nmate': number of materials
%   'nsect': number of sections
%   'ndofn': number of DOFs per node
%   'nevab': size of elemental stiffness matrix
%   'lnods': [ node1 node2 node3 node4;
%            node1 node2 node3 node4;
%            ... ];  each row is for each element
%   Note size of 'lnods' is (number of nodes x number of nodes per element)
%--- V A R I A B L E   D E F I N I T I O N---------------------------------

global title
global npoin nelem nnode nmate nsect ndofn
global lnods
global COORD
global ELEM_DATA
global SECTION_DATA
global MATERIAL_DATA


% Scan control cards
tline = fgetl(ifid);

while ~feof(ifid) % if not the end of file, continue scanning
    
    % comment line
    if ~isempty(findstr(tline,'**')) ~= 0
        tline = fgetl(ifid);    
    elseif ~isempty(findstr(tline,'*tit'))~=0
        tline = fgetl(ifid);
        title = tline;
        tline = fgetl(ifid);
    elseif ~isempty(findstr(tline,'*contr'))~=0
        tline = fgetl(ifid);
        [token,rem] = strtok(tline, ', =');
        while ~isempty(token)   % to the end of this line
            if strcmp(token,'npoin')
                [token,rem] = strtok(rem,' ,=');
                npoin = str2double(token);
            elseif strcmp(token,'nelem')
                [token,rem] = strtok(rem,' ,=');
                nelem = str2double(token);
            elseif strcmp(token,'nsect')    % get the number of sections
                [token,rem] = strtok(rem,' ,=');
                nsect = str2double(token);
            elseif strcmp(token,'nmate')    % get the number of materials
                [token,rem] = strtok(rem,' ,=');
                nmate = str2double(token);
            elseif strcmp(token,'ndofn')    % get the number of DOFs per node
                [token,rem] = strtok(rem,' ,=');
                ndofn = str2double(token);
            end
            [token,rem] = strtok(rem,' ,=');
        end
        tline = fgetl(ifid);
    % read nodal coordinate values
    elseif ~isempty(findstr(tline,'*node'))~=0
        COORD = fscanf(ifid,'%d %f %f %f\n',[4,npoin]);
        COORD = COORD';
        tline = fgetl(ifid);
    % read elemental connectivity data
    % #, node1, node2, node3, node4, section #, material #
    elseif ~isempty(findstr(tline,'*elem'))~=0
        [token,rem] = strtok(tline, ', =');
        while ~strcmp(token,'nnode')
            [token,rem] = strtok(rem, ', =');
        end
        [token,rem] = strtok(rem, ', =');
        nnode = str2double(token);            
        ELEM_DATA = fscanf(ifid,'%d\n',[7,nelem]); 
        ELEM_DATA = ELEM_DATA';
        lnods = ELEM_DATA(:,2:(nnode+1));
        tline = fgetl(ifid);
    % read material data
    % #, Young's modulus, poisson ratio
    elseif ~isempty(findstr(tline,'*material')) ~=0
        MATERIAL_DATA =fscanf(ifid,'%i,%f,%f\n',[3,nmate]);
        MATERIAL_DATA = MATERIAL_DATA';
        tline = fgetl(ifid);
    % read sectional properties
    % #, thickness
    elseif ~isempty(findstr(tline,'*section')) ~=0    % read sectional data
        SECTION_DATA = fscanf(ifid,'%i,%f\n',[2,nsect]);
        SECTION_DATA = SECTION_DATA';
        tline = fgetl(ifid);           
    % Pass unidentified command
    else
        tline = fgetl(ifid);
    end    
end

fclose(ifid);
    
return