function Plot_result(visual_data,coord_r)

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

global COORD ELEM_DATA nelem nnode

dlmwrite('visual_data.dat',visual_data);

ifid=fopen('visual_data.dat','rt');
tline = fgetl(ifid);
i=1;
while ~feof(ifid)
    for j=1:nnode
        [token,rem]=strtok(tline,', =');
        [token,rem]=strtok(rem,', =');
        PostNode{i}(j)=str2double(token);
        [token,rem]=strtok(rem,', =');
        Stress_xx{i}(j)=str2double(token);
        [token,rem]=strtok(rem,', =');
        Stress_yy{i}(j)=str2double(token);
        [token,rem]=strtok(rem,', =');
        Stress_xy{i}(j)=str2double(token);
        [token,rem]=strtok(rem,', =');
        Disp_x{i}(j)=str2double(token);
        [token,rem]=strtok(rem,', =');
        Disp_y{i}(j)=str2double(token);
        % add 
        
        tline = fgetl(ifid);  
    end
    i=i+1;
end
fclose(ifid);

% dimens = (max(Node.cord) - min(Node.cord))/2;	
dimens = (max(COORD(:,2:end)) - min(COORD(:,2:end)))/2;	
%dimension of the mesh
origin = dimens+min(COORD(:,2:end));			
%origin(center) of the mesh
tmp = max(COORD(:,2:end));
origin(3) = tmp(3);


tmpCube.faces    =[1 2 3 4];

figure('position',[100 100 800 600])
title('xx stress')
axis off;
xlabel('x')
ylabel('y')
daspect([1 1 1]) % last term is scale factor for z direction
colormap('jet')
for ielem=1:nelem
    
    xx = [];
    
    for jnode=1:nnode
        ND_Ind = find(COORD(:,1)==ELEM_DATA(ielem,jnode+1));
        tmpCube.vertices(jnode,:) = coord_r(ND_Ind,2:3);
    end 
    
    xx = Stress_xx{ielem}(:);
    
    tmpCube.facecolor='interp';
    tmpCube.FaceVertexCData=xx;
    patch(tmpCube,'EdgeColor','b');
end
figure('position',[100 100 800 600])
title('yy stress')
axis off;
xlabel('x')
ylabel('y')
daspect([1 1 1])
colormap('jet')
for ielem=1:nelem
    
    yy = [];
    
    for jnode=1:nnode
        ND_Ind = find(COORD(:,1)==ELEM_DATA(ielem,jnode+1));
        tmpCube.vertices(jnode,:) = coord_r(ND_Ind,2:3);
    end 
    
    yy = Stress_yy{ielem}(:);
    
    tmpCube.facecolor='interp';
    tmpCube.FaceVertexCData=yy;
    patch(tmpCube,'EdgeColor','b');
end
figure('position',[100 100 800 600])
title('xy stress')
axis off;
xlabel('x')
ylabel('y')
daspect([1 1 1])
colormap('jet')
for ielem=1:nelem
    
    xy = [];
    
    for jnode=1:nnode
        ND_Ind = find(COORD(:,1)==ELEM_DATA(ielem,jnode+1));
        tmpCube.vertices(jnode,:) = coord_r(ND_Ind,2:3);
    end 
    
    xy = Stress_xy{ielem}(:);
    
    tmpCube.facecolor='interp';
    tmpCube.FaceVertexCData=xy;
    patch(tmpCube,'EdgeColor','b');
end
figure('position',[100 100 800 600])
title('x displacement')
axis off;
xlabel('x')
ylabel('y')
daspect([1 1 1])
colormap('jet')
for ielem=1:nelem
    
    x = [];
    
    for jnode=1:nnode
        ND_Ind = find(COORD(:,1)==ELEM_DATA(ielem,jnode+1));
        tmpCube.vertices(jnode,:) = coord_r(ND_Ind,2:3);
    end 
    
    x = Disp_x{ielem}(:);
    
    tmpCube.facecolor='interp';
    tmpCube.FaceVertexCData=x;
    patch(tmpCube,'EdgeColor','b');
end
figure('position',[100 100 800 600])
title('y displacement')
axis off;
daspect([1 1 1])
colormap('jet')
for ielem=1:nelem
    
    y = [];
    
    for jnode=1:nnode
        ND_Ind = find(COORD(:,1)==ELEM_DATA(ielem,jnode+1));
        tmpCube.vertices(jnode,:) = coord_r(ND_Ind,2:3);
    end 
    
    y = Disp_y{ielem}(:);
    
    tmpCube.facecolor='interp';
    tmpCube.FaceVertexCData=y;
    patch(tmpCube,'EdgeColor','b');
end
end