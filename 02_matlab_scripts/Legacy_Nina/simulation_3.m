%% simulation of experimental platform (gel, piezo actuators)
clc;

%% variables
f_max = 4; %[Hz]
%%  create the computational grid
% diameter: 15cm, thickness: 1.4cm
spacing = 5e-4;
x_size = 0.25;            %[m] total size x-direction
y_size = 0.25;            %[m] total size y-direction
z_size = 0.05;           %[m] total size z-direction
dx = spacing;             % grid point spacing in the x direction [m]
dy = spacing;             % grid point spacing in the y direction [m]
dz = spacing;             % grid point spacing in the z direction [m]
Nx = ceil(x_size/dx);     % number of grid points in the x direction
Ny = ceil(y_size/dy);     % number of grid points in the y direction
Nz = ceil(z_size/dz);     % number of grid points in the z direction
disp('Number of grid points: '); disp(Nx*Ny*Nz);
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% simulation time and step size
% dt = 0.05 * 1 / (2*f_max); % Nyquist
dt = 0.3 * spacing / 1440.3; % CFL*dx/c_max CFL 0.3 default 
t_end = dt*100; % [s]
kgrid.t_array = 0:dt:t_end;
%% Geometry
% binary gel mask 
diameter_gel = 0.15; % [m]
thickness_gel = 0.014; % [m]
radius = floor(ceil(diameter_gel/dx)/2); %ceil(0.15/2 / dx);
center = floor(Nx/2); 
[x_gel,y_gel] = meshgrid(1:Nx, 1:Ny);
distance = (x_gel-center).^2+(y_gel-center).^2;
gel_cond = distance < radius^2;
thickness_grid = ceil(thickness_gel/dz); 
gel = zeros(Nx, Ny, Nz);
gel(:,:,5:5+thickness_grid-1) = double(repelem(gel_cond,1,1,thickness_grid));
gel_surface_z = 5+thickness_grid-1; % z grid coordinates of top surface

%% define the properties of the propagation medium
% initialize 
medium.sound_speed_compression = zeros(Nx,Ny,Nz);                               % [m/s]
medium.sound_speed_shear = zeros(Nx,Ny,Nz);                                     % [m/s]
medium.alpha_coeff_compression = zeros(Nx,Ny,Nz);                               % [10–3 m–1] 
medium.alpha_coeff_shear = zeros(Nx,Ny,Nz);                                     % [10–3 m–1] 
medium.density = zeros(Nx,Ny,Nz);                                               % [kg/m^3]
% air
medium.sound_speed_compression = medium.sound_speed_compression + (1-gel) * 344;% [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + (1-gel) * 0.000001;              % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + (1-gel) * 0.5;% [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + (1-gel) * 0.5;            % [10–3 m–1] 
medium.density = medium.density + (1-gel) * 1.225;                              % [kg/m^3]
% gel
medium.sound_speed_compression = medium.sound_speed_compression + gel * 1440.3; % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + gel * 201.74;             % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + gel * 0.5;    % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + gel * 0.5;                % [10–3 m–1] 
medium.density = medium.density + gel * 923.47;                                 % [kg/m^3]
%% boundary conditions
% open boundaries 

% This is achieved by setting source.p_mode and source.u_mode to ‘dirichlet’. 
% In this case, at each time step, the input pressure and velocity values are 
% used to replace the existing values at the grid points specified by source.p_mask 
% and source.u_mask, rather than adding to them.
%% define the source
% define a square source element
source_radius = 0.005/dx;  % [grid points]
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, gel_surface_z) = 1;
% define a time varying sinusoidal source
source_freq = f_max; %2/(t_end);      % [Hz]
source_mag = 0.0000001;           % [m]
source.uz = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
% smooth the source
source.uz = filterTimeSeries(kgrid, medium, source.uz);
%% define the sensor
% define a series of Cartesian points to collect the data
% x = [-Nx/2+2:2:Nx/2-2] * dx;     % [m]
% y = x;                           % [m]
% z = ones(size(x)) * 0.5* (gel_surface_z * dz); % [m]
% sensor.mask = [x; y; z];
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(2:2:Nx-2, 2:2:Ny-2, gel_surface_z-1) = 1;
%% plot gel, source and sensor
voxelPlot(gel);
start_ = [Nx/2 - source_radius, Ny/2 - source_radius, gel_surface_z];
size_ = [2*source_radius,2*source_radius,1];
voxel(start_, size_, 'b', 1); 
voxelPlot(sensor.mask);

%% run the simulation 
sensor.record = {'u'} % u particle velocity
% PML
% lambda=c0/f0;    % wavelength in lighter product (m)
% h=lambda/6;                 % desired spatial step (m)
% PML_size=5;             % size of the PML in grid points
% input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false};
sensor_data = pstdElastic3D(kgrid, medium, source, sensor); %, input_args{:});