function main()
% written by Lehan Yao at 4/15/2020
% functions inpolyhedron.m and readObj.m are downloaded from MathWorks FileExchange
% Bernard Abayowa (2020). readObj (https://www.mathworks.com/matlabcentral/fileexchange/18957-readobj), MATLAB Central File Exchange. Retrieved April 15, 2020
% Sven (2020). inpolyhedron (https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume), MATLAB Central File Exchange. Retrieved April 15, 2020.
% Sometimes the inpolyhedron function can fail, resulting some atoms outside the particle. 
% When this happens, try to change the orientation of the object.
% inputs: 
    fName = 'examples/mid2.obj'; %name of the obj file to read (I only considered cases with only one object, please triangularize the meshes, non-triangular mesh is not supported)
    coarse_size = 0.667; %size of the "atoms" (length, unitless)
    box_size = 130;  %the volume to consider: box_size-by-box_size-by-box_size (length, unitless)
    ligand_density = 0.6; %ligands will be placed onto every meshes with this density (count*length-2, unitless)
    ligand_dist = 3; %ligand-to-surface distance (length, unitless)
% outputs:
    fileIDatom = fopen('perpendicular_64.51_10.68_atom_coordinate.txt','w');
    fileIDligand = fopen('0.6_perpendicular_64.51_10.68_ligand_coordinate.txt','w');
% atom_coordinate.txt, in which the coordinates of atoms arranged in a simple cube
% lattice inside the polyhedron volume confined by the obj file is recorded in a format of [x, y, z].
% ligand_coordinate.txt, positions of ligands, same as above
    disp(['points to consider: ', num2str(floor(box_size/coarse_size)^3)])
    shape = readObj(fName);
    object.vertices = shape.v;
    object.faces = shape.f.v;
    units= place_atoms(object,coarse_size,box_size);
    ligands = place_ligand(object,ligand_density,ligand_dist);
    figure(1);clf;
    hold on
    for i = 1 :size(object.faces,1)
        v_id = object.faces(i,:);
        mesh = [object.vertices(v_id(1),:);
                object.vertices(v_id(2),:);
                object.vertices(v_id(3),:)];
        fill3(mesh(:,1),mesh(:,2),mesh(:,3),1)
    end
    scatter3(units(:,1),units(:,2),units(:,3),'r')
    hold on
    axis equal
    scatter3(ligands(:,1),ligands(:,2),ligands(:,3),4,'b')
    xlabel("x")
    ylabel("y")
    zlabel("z")
    set(gcf,'Position', [100 100 512 512]);
    xlim([-box_size/2,box_size/2])
    ylim([-box_size/2,box_size/2])
    zlim([-box_size/2,box_size])
    fprintf(fileIDatom,'%f %f %f\n',units');
    fclose(fileIDatom);
    disp(['number of atoms: ', num2str(length(units))])

    fprintf(fileIDligand,'%f %f %f\n',ligands');
    fclose(fileIDligand);
    disp(['number of ligands: ', num2str(length(ligands))])
end

function units= place_atoms(object,coarse_size,box_size)
    [X, Y, Z] = meshgrid([-box_size/2 : coarse_size : box_size/2]);
    nGrid=length(X);
    QP = [X(:) Y(:) Z(:)];
    indexIntersect = inpolyhedron(object,QP);
    mask = logical(reshape(indexIntersect, [nGrid nGrid nGrid]));
    XL=reshape(X,[nGrid^3,1]);
    YL=reshape(Y,[nGrid^3,1]);
    ZL=reshape(Z,[nGrid^3,1]);
    maskL=reshape(mask,[nGrid^3,1]);
    units=[XL(maskL),YL(maskL),ZL(maskL)];
end

function ligands = place_ligand(object,density,dist)
ligands = [];
    for i = 1 :size(object.faces,1)
        v_id = object.faces(i,:);
        mesh = [object.vertices(v_id(1),:);
                object.vertices(v_id(2),:);
                object.vertices(v_id(3),:)];
        disp(mesh);
        l1 = mesh(1,:)-mesh(2,:);
        l2 = mesh(3,:)-mesh(2,:);
        disp(l1)
        n = cross(l1,l2);
        n = n/norm(n);
        x_ = n;
        y_ = l1/norm(l1);
        z_ = cross(x_,y_);
        z_ = z_/norm(z_);
        %quiver3(mesh(1,1),mesh(1,2),mesh(1,3),x_(1),x_(2),x_(3),'r')
        %quiver3(mesh(1,1),mesh(1,2),mesh(1,3),y_(1),y_(2),y_(3),'g')
        %quiver3(mesh(1,1),mesh(1,2),mesh(1,3),z_(1),z_(2),z_(3),'b')
        x = [1,0,0];
        y = [0,1,0];
        z = [0,0,1];
        M_trans = [x*x_', x*y_', x*z_';
                   y*x_', y*y_', y*z_';
                   z*x_', z*y_', z*z_'];
        M_trans_inv = inv(M_trans);
        mesh_ = M_trans_inv*mesh';
        mesh_ = mesh_';
        [grid_y_,grid_z_] = meshgrid(min(mesh_(:,2))+ 0.5/density: 1/density: max((mesh_(:,2))),min(mesh_(:,3))+ 0.5/density : 1/density: max((mesh_(:,3))));
        grid_y_ = grid_y_(:);
        grid_z_ = grid_z_(:);
        idx = inpolygon(grid_y_,grid_z_,mesh_(:,2),mesh_(:,3));
        if(sum(idx)>0)
            ligand_2D_ = [grid_y_(idx),grid_z_(idx)];
            ligand_ = [repmat(mesh_(1,1),length(grid_y_(idx)),1), ligand_2D_];
            ligand = M_trans*ligand_';
            ligand = ligand';
            indicators = inpolyhedron(object,ligand - repmat(n*dist,size(ligand,1),1));
            indicator = mean(indicators);
            if indicator<0.5
                ligand = ligand - repmat(n*dist,size(ligand,1),1); 
            else
                ligand = ligand + repmat(n*dist,size(ligand,1),1); 
            end
            ligands = [ligands;ligand];
        else
            disp('Warrning: No ligand placed onto this mesh')
        end
    end
end