function [GL] = GetDLoad(ngaus,IDArray)

global nelem nnode ndofn nevab npoin neq
global lnods
global ELEM_DATA
global SECTION_DATA
global COORD

%%% conditions %%%
% ngaus == 2
% Only normal stress in x-direction is considered
% Q4 element
%%%%%%%%%%%%%%%%%%

% set uniform pressure
sigma0 = -10000;

% set variables
xipos     = 0.0;
etapos    = 0.0;
xiweight  = 0.0;
etaweight = 0.0;
detJ =0.0;

SFD = zeros(2,nnode);
coord2D = zeros(nnode,2);
Jmat = zeros(2,2);

% initialize global load vectors
GL = zeros(neq,1); % global load vector

% set node-coordinate matrix, NODE 
% node1: [x-coord  y-coord;
% node2:  x-coord  y-coord;
%         ...             ]
NODE=zeros(0,2);

for inode=1:npoin
	NODE(size(NODE,1)+1,1)=COORD(inode,2);
	NODE(size(NODE,1),2)=COORD(inode,3);
end

% find x-coord of left-end node, y-coord of bottom/top nodes
left_x = min(NODE(:,1)); % set left_x: x-coord of left end

%{
bot_y  = min(NODE(:,2)); % set bot_y: y-coord of bottom
top_y  = max(NODE(:,2));

% find left-end elements, left bottom node and left top node
left_nodes = zeros(0,1);
leftbot_node = 0;
lefttop_node = 0;
for inode=1:npoin
	if NODE(inode,1)==left_x
		left_nodes(size(left_nodes,1)+1,1) = inode;
		if NODE(inode,2)==bot_y
			leftbot_node = inode;
		elseif NODE(inode,2)==top_y
			lefttop_node = inode;
		end
	end
end
%}

left_elems = zeros(0,1); % storage of left-end elements

for ielem=1:nelem
	for inode=1:nnode
		if COORD(lnods(ielem,inode),2) == left_x
			left_elems(size(left_elems,1)+1,1) = ielem;
			break;
		end
	end
end

% loop
for i=1:size(left_elems,1)
	ielem = left_elems(i,1);
	no_sect = ELEM_DATA(ielem,6);
	thick = SECTION_DATA(no_sect,2);

	eload = zeros(ndofn*nnode,1); % initialize elemental load vector

	% find node order
	CASE = 0;
	nodes = zeros(4,3);
	for inode=1:nnode
		nod = lnods(ielem,inode);
		nodes(inode,1) = nod;
		nodes(inode,2) = COORD(nod,2);
		nodes(inode,3) = COORD(nod,3);
	end
	
	if nodes(1,2) == left_x
		if nodes(1,3) < nodes(4,3)
			CASE=1;
		else
			CASE=4;
		end
	else
		if nodes(2,3) < nodes(3,3)
			CASE=2;
		else
			CASE=3;
		end
	end

	for igaus=1:ngaus
		% find Jmat
		for jgaus=1:ngaus
			SF = zeros(1,nnode);
			SFD = zeros(2,nnode);
			coord2D = zeros(nnode,2);

			[xipos,etapos,xiweight,etaweight]= GetGPV(ngaus,igaus,jgaus,xipos,etapos,xiweight,etaweight);
			[SFD] = GetSFD(SFD,xipos,etapos);
			[coord2D] = GetCoord(ielem,coord2D);
			[detJ,Jmat] = GetDJabob(SFD,coord2D,detJ,Jmat);
		end

		% find Shape Functions
		xiweight = 1.0; etaweight = 1.0;
		if CASE == 1
			xipos = -1.0;
			etapos = (-1)^igaus/sqrt(3);
			s0 = sigma0;
		elseif CASE == 2
			xipos = (-1)^igaus/sqrt(3);
			etapos = 1.0;
			s0 = -sigma0;
		elseif CASE == 3
			xipos = 1.0;
			etapos = -(-1)^igaus/sqrt(3);
			s0 = -sigma0;
		else
			xipos = -(-1)^igaus/sqrt(3);
			etapos = -1.0;
			s0 = sigma0;
		end

		[SF] = GetSF(SF,xipos,etapos);

		fprintf('ielem: %d, igaus: %d, sigma0: %f, CASE: %d\n',ielem,igaus,sigma0,CASE);

		% find elemental load vector
		eload = eload + xiweight*etaweight*s0*thick*[SF(1)*(-Jmat(2,2));
		SF(1)*Jmat(2,1);
		SF(2)*(-Jmat(2,2));
		SF(2)*Jmat(2,1);
		SF(3)*(-Jmat(2,2));
		SF(3)*Jmat(2,1);
		SF(4)*(-Jmat(2,2));
		SF(4)*Jmat(2,1)]
	end

    % when using IDArray

	fprintf('\nelem: %d, CASE: %d\n',ielem,CASE);
	for inode=1:nnode
        nod = lnods(ielem,inode); % node number
		fprintf('%d(',nod);
        for jdofn=1:ndofn
            pos =(inode-1)*ndofn + jdofn; % numbering all dofn in this element
            KK(pos) = IDArray(jdofn, nod);
			fprintf('%d:%d ',pos,KK(pos));
            % 'IDArray': [ ndofn x npoin ] matrix that stores equation number
            % IDArray gives (1) 0 for Non-active DOF
            %               (2) numbering for all Active DOFs in the whole structure
            %                   from 0 to nevab
        end
		fprintf(')\n');
    end

	% Assemble element load vector into global load vector
	eload
	for ievab=1:nevab
		fprintf('ievab: %d, KK: %d, eload: %f\n',ievab,KK(ievab),eload(ievab));
		if( KK(ievab) <= 0 )
			continue;
		end
		I=KK(ievab);
		GL(I,1) = GL(I,1) + eload(ievab);
	end
	fprintf('\n -------------------------------------------------\n\n');
end

fprintf('\n +++ GL +++ \n');
GL'

quit();

