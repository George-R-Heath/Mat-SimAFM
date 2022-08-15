clear 

PDBStruct = pdbread('Br.pdb');  %load pdb from local disk
%pdb_id = '2abm'; % Or load from PDB online 
% egs: GroEL = 1ss8; Actin = 6BNO; ATP Synthase = 1QO1; COVID spike = 6VXX;
%PDBStruct = pdbread(['http://www.rcsb.org/pdb/files/',pdb_id,'.pdb']);

%% Settings
r =1;             %tip radius(Å)
angle = 8;         %cone angle (o)
pix_per_ang =  1;  %sampling (pix/Å)
noise = 1;          %noise (rms Å)

%rotation about different axis:
theta_x = 0;    
theta_y = 0;
theta_z = 0;

z_thresh = 0.7;  % Set the fraction (0-1) of coordinates to exclude (eg membrane embedded fraction)

% Calculations
for i = 1:numel(PDBStruct.Model)
    if i ==1
    coords = [PDBStruct.Model(i).Atom.X ; PDBStruct.Model(i).Atom.Y; PDBStruct.Model(i).Atom.Z];
    else
    coords = [coords, [PDBStruct.Model(i).Atom.X ; PDBStruct.Model(i).Atom.Y; PDBStruct.Model(i).Atom.Z]];
    end
end
coords = coords';

if theta_x >0
coords_rx(:,2) = coords(:,2)*cos(theta_x*pi/180) - coords(:,3)*sin(theta_x*pi/180);
coords_rx(:,3) = coords(:,2)*sin(theta_x*pi/180) + coords(:,3)*cos(theta_x*pi/180);
coords(:,2:3) = coords_rx(:,2:3);
clear coords_rx
end

if theta_y >0
coords_ry(:,1) = coords(:,1)*cos(theta_y*pi/180) + coords(:,3)*sin(theta_y*pi/180);
coords_ry(:,3) = coords(:,3)*cos(theta_y*pi/180) - coords(:,1)*sin(theta_y*pi/180);
coords(:,1) = coords_ry(:,1);coords(:,3) = coords_ry(:,3);
clear coords_ry
end

if theta_z >0
v=[coords(:,1)';coords(:,2)'];
R = [cos(theta_z*pi/180) sin(theta_z*pi/180); -sin(theta_z*pi/180) cos(theta_z*pi/180)];
so = R*v;  
coords(:,1) = so(1,:);
coords(:,2) = so(2,:);
clear so R v
end


if z_thresh>0
 z_thresh = z_thresh.*(max(coords(:,3))-min(coords(:,3))) + min(coords(:,3));
 pos = coords(:,3)>z_thresh;
 coords = coords(pos,1:3);
end

coords(:,3) = coords(:,3)-min(coords(:,3));

rs = r*pix_per_ang;
fspace = (max(coords(:,3))-r)*tan(angle*pi/180)*pix_per_ang+1;
end_pos = [floor(min(coords(:,1).*pix_per_ang)-rs-fspace),ceil(max(coords(:,1).*pix_per_ang)+rs+fspace);floor(min(coords(:,2).*pix_per_ang)-rs-fspace),ceil(max(coords(:,2).*pix_per_ang)+rs+fspace)];

img_n=[];
    clear coords_s
coords_s(:,3) = coords(:,3);
coords_s(:,1) = coords(:,1).*pix_per_ang;
coords_s(:,2) = coords(:,2).*pix_per_ang;

img = zeros(end_pos(1,2)-end_pos(1,1)+1,end_pos(2,2)-end_pos(2,1)+1);
    
for i=1:numel(coords(:,3))
 offs_x =  round(coords_s(i,1))-coords_s(i,1);
 offs_y =  round(coords_s(i,2))-coords_s(i,2);
    [dx,dy]=ndgrid((-rs-fspace):(rs+fspace));
    dxl = dx(:)-offs_x;
    dyl = dy(:)-offs_y;
    
    for j = 1:numel(dxl)   
    dh(j) = real(sqrt(-((dxl(j)/pix_per_ang).^2)-((dyl(j)/pix_per_ang).^2)+r^2))-r;    %tip interaction
    if dh(j)>-r
    h(j) = coords_s(i,3) + dh(j);
    else
    di = real(sqrt(((dxl(j)/pix_per_ang).^2)+((dyl(j)/pix_per_ang).^2)))-r;
    h(j) = coords_s(i,3)-r - (di)/tan(angle*pi/180);
    end
    
    pos_x = round(coords_s(i,1) + dxl(j)-end_pos(1)+1);
    pos_y = round(coords_s(i,2) + dyl(j)-end_pos(2)+1);
        if img(pos_x, pos_y) < h(j)
        img(pos_x, pos_y) = h(j);
        end
    end
     
% end
end

img_n = img + noise*randn(size(img));

%plotting:
%plot3(coords(:,1),coords(:,2),coords(:,3),'.')
load('lutafm.mat')
imagesc(img_n)
colormap(AFM)
axis equal tight
