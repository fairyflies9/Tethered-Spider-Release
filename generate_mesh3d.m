
%This file will generate a rectangular "leaf" that resists bending and stretching. 
%It will be held in place by a line of target points along its center

L = 3;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)

drag_length = 1.0;                  % drag line length (m)
npts_drag = ceil(2*(drag_length/L)*N);  % number of points along the span
npts = npts_drag;	    % total number points on each dragline 
ds = drag_length/(npts_drag-1);     % drag line mesh spacing (m)

mesh_name = 'spider3d_';      % structure name
offset = -0.25;	                    % plate offset from center of domain (m)

kappa_spring = 3e1;               % spring constant (Newton)
kappa_beam = 5.0e-6;                % beam stiffness constant (Newton m^2)
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information

vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');

fprintf(vertex_fid, '%d\n', npts);

for s = 0:npts_drag-1
   X(1) = offset;
   X(2) = 0.0;
   X(3) = 0.0 + 0.5*drag_length  - s*ds;
   fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));
end 

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the spring information
spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');

fprintf(spring_fid, '%d\n', npts-1);

for s = 0:npts-2
   fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, kappa_spring*ds/(ds^2), ds);
end


fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the beam information
beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');

fprintf(beam_fid, '%d\n', npts-2);

for s = 0:npts-2
   fprintf(beam_fid, '%d %d %d %1.16e\n', s, s+1, s+2, kappa_beam*ds/(ds^4));
end


fclose(beam_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information
    target_fid = fopen([mesh_name num2str(N) '.mass'], 'w');
    
    fprintf(target_fid, '%d\n', 1);  %note the number is only npts_chord since we only want a line of target points
    
    fprintf(target_fid, '%d %1.16e %1.16e\n', npts-1, 0.000002, kappa_target*ds/(ds^2));

    fclose(target_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
