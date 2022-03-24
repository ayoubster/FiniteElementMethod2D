% Function
% Function project_FEM(file, BCvector)
%
% Takes in a .msh file and outputs the solution to the heat equation given
% the boundary conditions
%
% Required Input:
%   file        .msh file, general polygonal mesh 
%   BCvector    [] vector, specifying the boundary values at the 2 boundaries
%
% Produces the .dat file which contains the values of the heat flux at various 
% x,y coordinates.
%
function Temp = project_FEM(file, BCvector)
close all;

fname=strrep(file,'.msh','.dat');% Save the output in a .dat file format
grd = []; 
grd = read_gmsh_2d(file, grd);% Reading the GMSH File from a absolute path.
% [L,K] = boundary(grd, [1 2 3 4 5 6 7 8]); % Call the global stiffness matrix with boundary condition
[L,K] = boundary(grd, [1 2], BCvector);
Temp = K\L; % Reading Temperature

for i = 1:length(grd.x)
    grd.z(i) = Temp(i); % Drawing Temperature in Z axis
end

write_tecplot(fname, grd);
end

function [L,K] = boundary(grd, tags, BCvector)
K = assembly(grd);
L = zeros(length(grd.x), 1);

for i_bound = 1:size(grd.bn,1) % looping over boundary lines
    for i_pt = 1:2
        pt = grd.bn(i_bound,i_pt); % read the two boundary points
        K(pt,:) = 0.0; 
        K(pt,pt) = 1.0;
        
        if grd.tags(i_bound) == tags(end) % in boundary tag 2 has temperature of BC2
            bval = BCvector(2);
        else
            bval = BCvector(1); % in boundary tag 1 has temperature of BC2
        end
        
        L(pt,1) = bval;
    end
end       

end

function [psi,d_psi_dr,d_psi_ds] = lagrange_derev(grd,j)

if j == 1
    psi = 1 - grd.r - grd.s; %lagrange basis funstion
    d_psi_dr = -1; %derivative of lagrange basis funstion w.r to r
    d_psi_ds = -1;%derivative of lagrange basis funstion w.r to s

elseif j == 2
    psi = grd.r;
    d_psi_dr = 1;
    d_psi_ds = 0;
    
elseif j == 3
    psi = grd.s;
    d_psi_dr = 0;
    d_psi_ds = 1;
    
else
    error('Out of range')
end


end

function [J,J_inv,J_det] = jacobian(e,grd)

dx_dr=0;
dx_ds=0;
dy_dr=0;
dy_ds=0;


for i = 1:3  % looping over 1 to 3 points of each triangle. 
    
    pt = grd.icon(e,i); % read each of the three points per elements. 
    
    [~,d_psi_dr,d_psi_ds] = lagrange_derev(grd,i);
    dx_dr = dx_dr + (grd.x(pt) * d_psi_dr); % converting master coordinate r to physical coordinate x
    dx_ds = dx_ds + (grd.x(pt) * d_psi_ds);% converting master coordinate s to physical coordinate x
    dy_dr = dy_dr + (grd.y(pt) * d_psi_dr);% converting master coordinate r to physical coordinate y
    dy_ds = dy_ds + (grd.y(pt) * d_psi_ds);% converting master coordinate s to physical coordinate y
end
   
J = [dx_dr dx_ds;dy_dr dy_ds]; % Jacobian Matrix formation.
J_inv = inv(J);% Inverse of Jacobian Matrix.
J_det = det(J);% Determinant of Jacobian Matrix.

end

function [d_psi_dx,d_psi_dy] = mas_to_phy(e,j,grd)
[~,J_inv,~] = jacobian(e,grd);
[~,d_psi_dr,d_psi_ds] = lagrange_derev(grd,j);
[d_psi] =  J_inv * [d_psi_dr;d_psi_ds]; % Converting master coordinate to physical coordinate              
 d_psi_dx = d_psi(1);
 d_psi_dy = d_psi(2);
end

function Ke = stiffness(grd)

ne = size(grd.icon,1);
n = size(grd.icon,2);      
Ke=zeros(n,n,ne);

for e = 1:ne                 
    for d = 1:n
        [d_psi_dx_d,d_psi_dy_d] = mas_to_phy(e,d,grd);
        for m = 1:n        
            [d_psi_dx_m,d_psi_dy_m] = mas_to_phy(e,m,grd);
            ngauss =length(grd.r);    % length of gauss point r
            for k = 1:ngauss                                
                [~,~,J_det] = jacobian(e,grd);             
                C = 1;   %Conductivity of materials. For simplicity assumed that C(at gauss point)=1.
                Ke(d,m,e) = Ke(d,m,e) + C*((d_psi_dx_d * d_psi_dx_m)+(d_psi_dy_d * d_psi_dy_m))*abs(J_det)*grd.w(k) ;
            end
        end
    end
end
end

function Kglobe = assembly(grd)
% the function is for assembling stiffness elemental matrix
Ke = stiffness(grd);
n = length(grd.x);
Kglobe = zeros(n,n);

for e = 1:size(grd.icon,1)
    for d = 1:size(grd.icon,2)
         iglobe = grd.icon(e,d);
        for m = 1:size(grd.icon,2)
         jglobe=grd.icon(e,m);          
         Kglobe(iglobe,jglobe) = Kglobe(iglobe,jglobe) + Ke(d,m,e);
        end    
    end
end
end

function out=write_tecplot(fname, grd)

fid = fopen(fname,'W');
fprintf(fid, 'TITLE = " 2D Finite Element Analysis"\n');
fprintf(fid, 'VARIABLES = "X", "Y", "Z" \n');
fprintf(fid, 'ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n', length(grd.x), size(grd.icon,1) );

for i=1:length(grd.x)
    fprintf(fid, '%16.16f , %16.16f , %16.16f\n', grd.x(i), grd.y(i), grd.z(i));
end

for i=1:size(grd.icon, 1) %number of elements
    for j = 1:3
        fprintf(fid, '%d ', grd.icon(i, j));
    end
    fprintf(fid, '\n ');
end

fclose(fid);
end

function grd = read_gmsh_2d(file, grd)

% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range
% bounded by row offsets R1 and R2 and column offsets C1 and C2.


% no of nodes is mentioned in 5th row and first column

N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);

%------- 2D Geometry

two_d_nodes = nodes(:,1:2);
elem_type   = elements(:,2);

%--- find the starting indices of 2D elements
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end
%----------------------------------------------

two_d_elements(1:N_e-two_ind,1:3) = 0;
k = 1;
for i = two_ind:N_e
   two_d_elements(k,1:3) = elements(i,6:8);
   k = k+1;
end

% %---- visualize in matlab ---------------------
% 
% figure(1)
% triplot(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2))
% xlabel('X','fontsize',14)
% ylabel('Y','fontsize',14)
% title('Project 2D GMASH TO MATLAB','fontsize',12)
% fh = figure(1);
% set(fh, 'color', 'white'); 

% getting x and y first ..
grd.x = nodes(:,1);
grd.y = nodes(:,2);

% count the number of tris
ntri = 0;
for i = 1:size(elements,1)
    elem_type = elements(i, 2);
    if (elem_type == 2) % yes im getting  a triangle
       ntri = ntri + 1; 
    end
end

% allocate icon matrix
grd.icon = zeros(ntri, 3);

% fill icon matrix
jj=1;
for i = 1:size(elements,1)
    elem_type = elements(i, 2);
    if (elem_type == 2) % yes im getting  a triangle
       grd.icon(jj, :) = elements(i, 6:8);
       jj = jj + 1;
    end
end

% count the number of boundary lines (segments)
nbn = 0;
for i = 1:size(elements,1)
    elem_type = elements(i, 2);
    if (elem_type == 1) % yes im getting  a boundary
       nbn = nbn + 1; 
    end
end

% allocate bn matrix
grd.bn = zeros(nbn, 2);
grd.tags = zeros(nbn,1);
jj = 1;

% fill bn and matrix
for i = 1:size(elements,1)
    elem_type = elements(i, 2);
    if (elem_type == 1) % yes im getting  a boundary
       grd.bn(jj, :) = elements(i, 6:7);
       grd.tags(jj,1) = elements(i, 4);
       jj = jj + 1;
    end
end
% Read the gauss point for 2nd order lagrange.
grd.r =      [0.666666666666667   0.166666666666667   0.166666666666667];
grd.s =      [0.166666666666667   0.166666666666667   0.666666666666667];
grd.w = [0.333333333333333   0.333333333333333   0.333333333333333];
end

