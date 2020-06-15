function [] = plotting_plane(Nx, Ny, Nz)
    sensors = zeros(Nx, Ny, Nz);
    sensors(2:10:Nx-2, 2:10:Ny-2, gel_surface_z-1) = 1;

    % [x y] = meshgrid(2:1:Nx-2); % SAME AS SENSORS IN SIMULATION

    [x y] = meshgrid(1:Nx); % SAME AS SENSORS IN SIMULATION
    % x = reshape(x, [1 16]);
    % y = reshape(y, [1 16]);
    pos = vel_to_pos(sensor_data.uz(size(sensor_data.uz,1)/thickness_grid*2 + 1:end,:), dt);
    max_z = max(max(pos));
    min_z = min(min(pos));
    max_z_outer = 0;
    min_z_outer = 100;

    source_radius = source_radius +2;


    for i=1:10:4000 % length(pos)
        pos_z = pos(:,i);
        pos_plane = zeros(Nx, Ny);
        pos_plane(find(gel(:,:,6) == 1)) = pos_z;
        new_pos = pos_plane;
        new_pos(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius) = 0;
        pos_z = reshape(new_pos, [Nx Ny]);
        if(max_z_outer < max(max(new_pos)))
            max_z_outer = max(max(new_pos));
        end
        if(min_z_outer > min(min(new_pos)))
            min_z_outer = min(min(new_pos));
        end
        surf(x,y,new_pos);
        zlim([min_z max_z]);
        ts = ['timestep '  num2str(i)  ' of '  num2str(length(sensor_data.uz(1,:)))];
        title(ts);
    %     scatter3(x,y,u_z);
        pause(0.001);
    end
end