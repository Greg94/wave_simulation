%% sim_skin_water_int.m
% Simulation to try and identify couplings in the system with
% the hand immersed in water and stimulated with a sinusoidal input.
% Gelatin physical parameters used to model the finger and palm
% Written by Bharat Dandu, adapted from simulation7.m written by Nina Fiona

clc;
clearvars;
close all

%% parameter
% Actuator input:
sinusoidal = 1;

%gaussian = 0;
%impulse = 1;


% Boundary conditions:
reflecting = 0;
bc_vel_bottom = 0; % first two layer vel = 0
bc_free_surface = 0; % stress on the boundary = 0
c_compression = 25.0; % assuming same compression wave speed for water and gelatin
% artificially low compression wave speed can potentially be used to minimize computation and test, but unrealistic. 
c_shear = 5.0; % shear wave speed in gelatin assuming shear wave speed in water is 0 (non existent)

% Other:
f_max = 300; %[Maximum frequency in consideration, Hz]
% Currently, this is well under the limit of CFL/2dt

%% computational grid
% Side view (x-z plane) - 
% origin
% _________________________________________
% |     _____________________|S_|___water  | z1
% |     |_____________gel___________|      | z2
% |                                  water | z3
%
% Width of hand - 8 cm. Width of water tank - 16 cm
% Length of hand - 15 cm. Length of tank - 24 cm
% Height of hand - 2.4 cm.  Height of water above - 2.4 cm, below - 3.2 cm 
% Source is 1cm into the length of the hand, midpoint of width
%
% PML inside is set to false, specify accurate dimensions
% Otherwise, Note that 10 grid points in each direction required for PML, 
% which'll be in the water
% """"""""""""""""""""""""""""""""""""""""" 
x_gel = 0.15;            % Length of gelatin layer [m]
x_water = 0.045;           % Length of water layer (on each side)
y_gel = 0.08;            % [m]
y_water = 0.04;            % [m]
z_gel = 0.024;            % [m]
z_water_a = 0.024;            % [m]
z_water_b = 0.032;            % [m]

ds = 4e-3;             % 2 mm spacing of the grid, x10 = 4cm used for PML
%ds = 8e-3;             % 2 mm spacing of the grid, x10 = 2cm used for PML

x_size = x_gel+2*x_water;            %[m] total size x-direction
y_size = y_gel+2*y_water;                                    %[m] total size y-direction
z_size = z_gel+z_water_a+z_water_b;    %[m] total size z-direction

Nx = round(x_size/ds);                               % number of grid points in the x direction
Ny = round(y_size/ds);                               % number of grid points in the y direction
Nz = round(z_size/ds);                               % number of grid points in the z direction

kgrid = kWaveGrid(Nx, ds, Ny, ds, Nz, ds);

%% simulation time and step size
dt = 0.05 * ds / c_compression; % CFL*dx/c_max CFL 0.3 default
t_end = 0.0625; %dt*500; % [s] 5 cycles at 80 Hz
kgrid.t_array = 0:dt:t_end;

%% Geometry of system : Binary matrix indicating explicitly defined locations
% 1. gel

gel = zeros(Nx, Ny, Nz);
g_xR = round(x_water/ds):(round((x_water+x_gel)/ds)-1);
g_yR = round(y_water/ds):(round((y_water+y_gel)/ds)-1);
g_zR = round(z_water_a/ds):(round((z_water_a+z_gel)/ds)-1);

gel(g_xR,g_yR,g_zR) = 1;    

% 2. water
water = ones(Nx, Ny, Nz) - gel; % Everything else

%% properties of the propagation mediums
% 0. initialize 
medium.sound_speed_compression = zeros(Nx,Ny,Nz);                                       % [m/s]
medium.sound_speed_shear = zeros(Nx,Ny,Nz);                                             % [m/s]
medium.alpha_coeff_compression = zeros(Nx,Ny,Nz);                                       % [10–3 m–1] Loss parameters - tweak later
medium.alpha_coeff_shear = zeros(Nx,Ny,Nz);                                             % [10–3 m–1] 
medium.alpha_coeff_shear = zeros(Nx,Ny,Nz);                                             % [10–3 m–1] 
medium.density = zeros(Nx,Ny,Nz);                                                       % [kg/m^3]

% 1. gel
medium.sound_speed_compression = medium.sound_speed_compression + gel * c_compression;  % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + gel * c_shear;                    % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + gel * 0.8;            % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + gel * 0.8;                        % [10–3 m–1] 
medium.density = medium.density + gel * 923.47;                                         % [kg/m^3]

% 2. water
medium.sound_speed_compression = medium.sound_speed_compression + water * c_compression;  % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + water * 0;                    % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + water * 0.0022;            % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + water * 0;                        % [10–3 m–1] 
medium.density = medium.density + water * 1.225;                                          % [kg/m^3]

% % 3. bottom - consider if adding a plastic layer, I consider this
% unnecessary
% medium.sound_speed_compression = medium.sound_speed_compression + bottom * c_compression;        % [m/s]
% medium.sound_speed_shear = medium.sound_speed_shear + bottom * c_shear;                 % [m/s]
% medium.alpha_coeff_compression = medium.alpha_coeff_compression + bottom * 0.5;         % [10–3 m–1] 
% medium.alpha_coeff_shear = medium.alpha_coeff_shear + bottom * 0.5;                     % [10–3 m–1] 
% medium.density = medium.density + bottom * 7700; % density steel                        % [kg/m^3]

%% boundary conditions
% defined by medium properties, PML at end of the script, for more bc look
% in earlier scripts

%% source
% define a 8mm square source element 

source.u_mask = zeros(Nx,Ny,Nz);
source_side = 0.008; %0.005/dx;  % [grid points]
src_offset = .01; % 1cm into tip of finger
src_sd_ind = round(source_side/(2*ds));
src_x_center = round((x_water+src_offset)/ds);

source.u_mask(ceil(src_x_center - src_sd_ind):ceil(src_x_center + src_sd_ind - 1) , ceil(Ny/2 - src_sd_ind):ceil(Ny/2 + src_sd_ind - 1), g_zR(1)) = 1;
nbr_of_sources = sum(source.u_mask, "all");

if(sinusoidal)
    % define a time varying sinusoidal velocity source
    
    f_stim = 80; % Frequency of sinusoidal stimulation in Hz
    source_mag = 0.0045;           % [m/s] - calculated assuming .3g acceleration,
    source.uz = source_mag * sin(2 * pi * f_stim * kgrid.t_array) .* tukeywin(length(kgrid.t_array),.25)';
    
    %v_smooth = sources_v .* smooth_beg;
    %source.uz = repelem(v_smooth,nbr_of_sources, 1);
end

figure(1);
subplot(2,1,1);
plot(kgrid.t_array, source.uz(end,:));
title("Actuator input velocity");
subplot(2,1,2);
plot(kgrid.t_array, vel_to_pos(source.uz(end,:),dt));
title("Actuator input position");   
% source.u_mode = 'dirichlet';
%% sensor
% all points are being saved:
% sensor.mask = ones(Nx,Ny,Nz);
% sensor.record = {'u_split_field'};

% only surface, 2 layers above is being saved:
% sensor.mask = zeros(Nx,Ny,Nz);
% sensor.mask(:,:,(g_zR(1)-2):(g_zR(1)+1)) = 1;

% only surface
surface_mask = zeros(Nx,Ny,Nz);
surface_mask(:,:,g_zR(1)) = 1;

sensor.mask = surface_mask;
%sensor.mask = gel;

boundary_mask = zeros(Nx,Ny,Nz);
boundary_mask(g_xR,g_yR(1:2),g_zR(1:2)) = 1;
boundary_mask(g_xR,g_yR((end-1):end),g_zR(1:2)) = 1;
boundary_mask(g_xR,g_yR(1:2),g_zR((end-1):end)) = 1;
boundary_mask(g_xR,g_yR((end-1):end),g_zR((end-1):end)) = 1;
boundary_mask(g_xR(1:2),g_yR,g_zR(1:2)) = 1;
boundary_mask(g_xR((end-1):end),g_yR,g_zR(1:2)) = 1;
boundary_mask(g_xR(1:2),g_yR,g_zR((end-1):end)) = 1;
boundary_mask(g_xR((end-1):end),g_yR,g_zR((end-1):end)) = 1;
boundary_mask(g_xR(1:2),g_yR(1:2),g_zR) = 1;
boundary_mask(g_xR((end-1):end),g_yR(1:2),g_zR) = 1;
boundary_mask(g_xR(1:2),g_yR((end-1):end),g_zR) = 1;
boundary_mask(g_xR((end-1):end),g_yR((end-1):end),g_zR) = 1;

%% plot gel, source and sensor
% voxelPlot(gel);
% voxelPlot(water);

% start_ = [Nx/2 - source_radius, Ny/2 - source_radius, gel_surface_z];
% size_ = [2*source_radius,2*source_radius,1];
% voxel(start_, size_, 'b', 1);
%% run the simulation 
sensor.record = {'u','p'} % u particle velocity
% PML
% lambda=c0/f0;    % wavelength in lighter product (m)
% h=lambda/6;      % desired spatial step (m)
if(reflecting)
    PML_size = 0;             % size of the PML in grid points
    pml_alpha = 0;
    input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false, 'PMLAlpha', pml_alpha};
else
    PML_size = 10;             % size of the PML in grid points
    input_args = {'PMLSize', PML_size, 'PlotPML', false,'DisplayMask' , boundary_mask , 'PMLInside', false,'DataCast', 'single', 'PlotFreq', 20, 'RecordMovie', true,'MovieProfile','MPEG-4'};
end
fprintf('Number of grid points: %i \n', Nx*Ny*Nz);
fprintf('Number of time steps: %i \n', kgrid.Nt);
input("Hit [Enter] to start the Simulation...");
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});




