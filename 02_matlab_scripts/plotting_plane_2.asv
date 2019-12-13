uz = sensor_data.uz;
pos = vel_to_pos(sensor_data.uz, dt);
[x y] = meshgrid(1:Nx); % SAME AS SENSORS IN SIMULATION
% nbr_of_sensors = size(uz,1);
nbr_of_timesteps = size(uz,2);
ts = 1 * nbr_of_timesteps;
% pos_top = vel_to_pos(uz(end-Nx*Ny+1:end,:), dt);
% pos = vel_to_pos(uz, dt);

uz_3d = reshape(sensor_data.uz,Nx,Ny,Nz,kgrid.Nt);
uz_2d_top = zeros(Nx,Ny,kgrid.Nt);
uz_2d_top(:,:,:) = uz_3d(:,:,gel_surface_z,:);
uz_2d_top = uz_2d_top .* repelem(gel_cond,1,1, ts);
posz_2d_top = vel_to_pos_2d(uz_2d_top, dt);

uz_bottom = uz(1:Nx*Ny,:);
max_z = max(posz_2d_top,[], "all")
min_z = min(posz_2d_top,[], "all")
maxv_b = max(uz_bottom,[], 'all');
minv_b = min(uz_bottom,[], 'all');

for t=1:20:ts %nbr_of_timesteps
    pos_t = posz_2d_top(:,:,t);
    surf(x,y,pos_t);
    zlim([-0.001 0.001]); %     zlim([2*min_z 2*max_z]);%     zlim([2.5e-3 0.6*max_v]);
    tss = ['timestep '  num2str(t)  ' of '  num2str(length(uz(1,:))) ' (' num2str(t_end*1000) ' ms)'];
    title(tss);
    pause(0.1);
end

index = find(sensor_data.uz(:,2));
plot(kgrid.t_array, sensor_data.uz(index(3),:));