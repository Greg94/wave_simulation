clc;
%% parameter
% Actuator input:
sinusoidal = 0;
gaussian = 0;
impulse = 1;
% Boundary conditions:
reflecting = 0;
bc_vel_bottom = 0; % first two layer vel = 0
bc_free_surface = 0; % stress on the boundary = 0
c_compression = 25.0; %artificially low compression wave speed
c_shear = 4.0; % shear wave speed in gelatin 
% Other:
f_max = 300; %[Hz]
%% computational grid
% _________________________________________
% |     _____________________________ air | y1
% |  x1 |            gel            |     |
% """"""""""""""""""""""""""""""""""""""""" y2
x1_in_m = 0.010;            % [m]
y1_in_m = 0.006;            % [m]
y2_in_m = 0.010;            % [m]
diameter_gel_in_m = 0.15;   % [m]
thickness_gel_in_m = 0.014; % [m]
spacing = 5e-3;             % [m]
x_size = diameter_gel_in_m + 2* x1_in_m;            %[m] total size x-direction
y_size = x_size;                                    %[m] total size y-direction
z_size = thickness_gel_in_m + y1_in_m + y2_in_m;    %[m] total size z-direction
dx = spacing;                                       % grid point spacing in the x direction [m]
dy = spacing;                                       % grid point spacing in the y direction [m]
dz = spacing;                                       % grid point spacing in the z direction [m]
Nx = ceil(x_size/dx);                               % number of grid points in the x direction
Ny = ceil(y_size/dy);                               % number of grid points in the y direction
Nz = ceil(z_size/dz);                               % number of grid points in the z direction
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% simulation time and step size
dt = 0.05 * spacing / c_compression; % CFL*dx/c_max CFL 0.3 default
t_end = 0.010; %dt*500; % [s]
kgrid.t_array = 0:dt:t_end;
%% Geometry 
% 1. gel
thickness_gel = ceil(thickness_gel_in_m/dz); 
gel_radius = floor(ceil(diameter_gel_in_m/dx)/2); 
gel_center = floor(Nx/2); 
[x_gel,y_gel] = meshgrid(1:Nx, 1:Ny);
distance = (x_gel-gel_center).^2+(y_gel-gel_center).^2;
gel_cond = distance < gel_radius^2;
gel = zeros(Nx, Ny, Nz);
gel_height = ceil(y2_in_m/dz);
gel_surface_z = gel_height+thickness_gel-1; % z grid coordinates of top surface
gel(:,:,gel_height:gel_surface_z) = double(repelem(gel_cond,1,1,thickness_gel));
% 2. bottom
bottom = zeros(Nx, Ny, Nz);
bottom(:,:,1:gel_height) = 1;
% 3. air
air = ones(Nx, Ny, Nz) - gel - bottom;
%% properties of the propagation mediums
% 0. initialize 
medium.sound_speed_compression = zeros(Nx,Ny,Nz);                                       % [m/s]
medium.sound_speed_shear = zeros(Nx,Ny,Nz);                                             % [m/s]
medium.alpha_coeff_compression = zeros(Nx,Ny,Nz);                                       % [10–3 m–1] 
medium.alpha_coeff_shear = zeros(Nx,Ny,Nz);                                             % [10–3 m–1] 
medium.density = zeros(Nx,Ny,Nz);                                                       % [kg/m^3]
% 1. gel
medium.sound_speed_compression = medium.sound_speed_compression + gel * c_compression;  % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + gel * c_shear;                    % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + gel * 0.5;            % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + gel * 0.5;                        % [10–3 m–1] 
medium.density = medium.density + gel * 923.47;                                         % [kg/m^3]
% 2. air
medium.sound_speed_compression = medium.sound_speed_compression + air * c_compression;  % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + air * c_shear;                    % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + air * 0.5;            % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + air * 0.5;                        % [10–3 m–1] 
medium.density = medium.density + air * 1.225;                                          % [kg/m^3]
% 3. bottom
medium.sound_speed_compression = medium.sound_speed_compression + bottom * c_compression;        % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + bottom * c_shear;                 % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + bottom * 0.5;         % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + bottom * 0.5;                     % [10–3 m–1] 
medium.density = medium.density + bottom * 7700; % density steel                        % [kg/m^3]
%% boundary conditions
% defined by medium properties, PML at end of the script, for more bc look
% in earlier scripts
%% source
% define a square source element
source.u_mask = zeros(Nx,Ny,Nz);
source_radius = 0.5; %0.005/dx;  % [grid points]
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, gel_surface_z) = 1;
nbr_of_sources = sum(source.u_mask, "all");
if(sinusoidal)
    % define a time varying sinusoidal source
    source_freq = f_max; %2/(t_end);      % [Hz]
    source_mag = 0.001;           % [m]
    stimulation_z = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
    tp = find(kgrid.t_array < (0.5/f_max));
    sources_v = source_mag * cos(2 * pi * source_freq * kgrid.t_array) * 2 * pi * source_freq;
    smooth_beg = 0.5* sin(2*pi*source_freq*kgrid.t_array + pi/8*source_freq) + 0.5;
    smooth_beg(tp(end):end) = 1;
    v_smooth = sources_v .* smooth_beg;
    source.uz = repelem(v_smooth,nbr_of_sources, 1);
end
if(impulse)
    tp = find(kgrid.t_array < (0.5/f_max));
    u_imp = -0.0001*(0.5* sin(2*pi*f_max*kgrid.t_array + pi/8*f_max) + 0.5);
    u_imp(tp(end):end) = 0;
    v_imp = [0 diff(u_imp)/dt];
    v_imp(tp(end)) = 0;
    source.uz = repelem(v_imp,nbr_of_sources, 1);  
end
if(gaussian)
    pos = -0.01 * gaussmf(kgrid.t_array, [0.0005 0.5/f_max]);
    vel_gaussian = [0 diff(pos)./dt];
    source.uz = repelem(vel_gaussian,nbr_of_sources, 1);   
end
figure(11);
subplot(2,1,1);
plot(kgrid.t_array, source.uz(end,:));
title("Actuator input velocity");
subplot(2,1,2);
plot(kgrid.t_array, vel_to_pos(source.uz(end,:),dt));
title("Actuator input position");   
% source.u_mode = 'dirichlet';
%% sensor
% all points are being saved:
sensor.mask = ones(Nx,Ny,Nz);
% only surface is being saved:
% sensor.mask = zeros(Nx,Ny,Nz);
% sensor.mask(:,:,gel_surface_z);
%% plot gel, source and sensor
voxelPlot(gel);
% voxelPlot(air);
% voxelPlot(bottom);
start_ = [Nx/2 - source_radius, Ny/2 - source_radius, gel_surface_z];
size_ = [2*source_radius,2*source_radius,1];
voxel(start_, size_, 'b', 1);
%% run the simulation 
sensor.record = {'u'} % u particle velocity
% PML
% lambda=c0/f0;    % wavelength in lighter product (m)
% h=lambda/6;      % desired spatial step (m)
if(reflecting)
    PML_size = 0;             % size of the PML in grid points
    pml_alpha = 0;
    input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false, 'PMLAlpha', pml_alpha};
else
    PML_size = 5;             % size of the PML in grid points
    input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false};
end
fprintf('Number of grid points: %i \n', Nx*Ny*Nz);
fprintf('Number of time steps: %i \n', kgrid.Nt);
input("Hit [Enter] to start the Simulation...");
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});




