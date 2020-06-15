%% simulation of experimental platform (gel, piezo actuators)


% Quadratic gel simulation (no air), PML, no bc velocity, stress = 0 on the outside, 
% uz = gaussian
% dirichlet/ no dirichlet

clc;

%% variables
f_max = 300; %[Hz]
%%  create the computational grid
% diameter: 15cm, thickness: 1.4cm
spacing = 5e-3;
x_size = 0.15;            %[m] total size x-direction
y_size = 0.15;            %[m] total size y-direction
z_size = 0.014;            %[m] total size z-direction
dx = spacing;             % grid point spacing in the x direction [m]
dy = spacing;             % grid point spacing in the y direction [m]
dz = spacing;             % grid point spacing in the z direction [m]
Nx = ceil(x_size/dx)+2;   % number of grid points in the x direction
Ny = ceil(y_size/dy)+2;   % number of grid points in the y direction
Nz = ceil(z_size/dz)+2    % number of grid points in the z direction
disp('Number of grid points: '); disp(Nx*Ny*Nz);
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% simulation time and step size
% dt = 0.05 * 1 / (2*f_max); % Nyquist
dt = 0.1 * spacing / 1440.3; % CFL*dx/c_max CFL 0.3 default
t_end = 0.5/300; % dt*4000; % [s]
kgrid.t_array = 0:dt:t_end;
%% Geometry
% % binary gel mask 
% diameter_gel = 0.15; % [m]
thickness_gel = 0.014; % [m]
% radius = floor(ceil(diameter_gel/dx)/2); %ceil(0.15/2 / dx);
% center = floor(Nx/2); 
% [x_gel,y_gel] = meshgrid(1:Nx, 1:Ny);
% distance = (x_gel-center).^2+(y_gel-center).^2;
% gel_cond = distance < radius^2;
thickness_grid = ceil(thickness_gel/dz); 
% gel = zeros(Nx, Ny, Nz);
z_height = 2; % grid points
% gel(:,:,z_height:z_height+thickness_grid-1) = double(repelem(gel_cond,1,1,thickness_grid));
gel_surface_z = z_height+thickness_grid-1; % z grid coordinates of top surface
% % get outer layer of gel only:
% gel_outer = get_boundary(Nx,Ny,gel_cond); imagesc(gel_outer);
% gel_outer3d = zeros(Nx,Ny,Nz);
% gel_outer3d(:,:,[z_height z_height+thickness_grid-1]) = double(repelem(gel_cond,1,1,2));
% gel_outer3d(:,:,[z_height+1:z_height+thickness_grid-2]) = double(repelem(gel_outer,1,1,thickness_grid-2));
gel = ones(Nx, Ny, Nz);
%% define the properties of the propagation medium
% initialize 
medium.sound_speed_compression = zeros(Nx,Ny,Nz);                               % [m/s]
medium.sound_speed_shear = zeros(Nx,Ny,Nz);                                     % [m/s]
medium.alpha_coeff_compression = zeros(Nx,Ny,Nz);                               % [10–3 m–1] 
medium.alpha_coeff_shear = zeros(Nx,Ny,Nz);                                     % [10–3 m–1] 
medium.density = zeros(Nx,Ny,Nz);                                               % [kg/m^3]
% air
medium.sound_speed_compression = medium.sound_speed_compression + (1-gel) * 344;% [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + (1-gel) * 0;              % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + (1-gel) * 0.5;% [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + (1-gel) * 0.5;            % [10–3 m–1] 
medium.density = medium.density + (1-gel) * 1.225;                              % [kg/m^3]
% gel
medium.sound_speed_compression = medium.sound_speed_compression + gel * 1440.3; % [m/s]
medium.sound_speed_shear = medium.sound_speed_shear + gel * 201.74;             % [m/s]
medium.alpha_coeff_compression = medium.alpha_coeff_compression + gel * 1;    % [10–3 m–1] 
medium.alpha_coeff_shear = medium.alpha_coeff_shear + gel * 1;                % [10–3 m–1] 
medium.density = medium.density + gel * 923.47;                                 % [kg/m^3]
%% boundary conditions
gel_outer3d = ones(Nx,Ny,Nz);
gel_outer3d(2:Nx-1, 2:Ny-1, 2:Nz-1) = 0;
source.s_mask   = gel_outer3d;
source.sxx      = zeros(1, kgrid.Nt);
source.syy      = source.sxx;
source.sxy      = source.syy;
source.s_mode   = 'dirichlet';
% open boundaries, fixed at bottom
% source.u_mask = zeros(Nx, Ny, Nz);

% source.u_mode= 'dirichlet';

% source.u_mask(:,:,1:z_height) = 1;
% nbr_of_sources_1 = sum(sum(sum(source.u_mask)));
%source.uz(1:nbr_of_sources_1, 1:kgrid.Nt) = zeros(nbr_of_sources_1, kgrid.Nt);
% This is achieved by setting source.p_mode and source.u_mode to ‘dirichlet’. 
% In this case, at each time step, the input pressure and velocity values are 
% used to replace the existing values at the grid points specified by source.p_mask 
% and source.u_mask, rather than adding to them.
%% define the source
% define a square source element
source_radius = 0.005/dx;  % [grid points]
source.u_mask = zeros(Nx,Ny,Nz);
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, gel_surface_z) = 2;
mask = reshape(source.u_mask, 1, Nx*Ny*Nz);
nbr_of_sources_1 = 0;
nbr_of_sources_2 = 0.5 * (sum(sum(sum(source.u_mask))) - nbr_of_sources_1);
% define a time varying sinusoidal source
source_freq = f_max; %2/(t_end);      % [Hz]
source_mag = 0.001;           % [m]
stimulation_z = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
tp = find(kgrid.t_array == (0.5/source_freq));
sources_v = source_mag * cos(2 * pi * source_freq * kgrid.t_array) * 2 * pi * source_freq;
smooth_beg = 0.5* sin(2*pi*source_freq*kgrid.t_array + pi/8*source_freq) + 0.5;
smooth_beg(tp:end) = 1;
v_smooth = sources_v .* smooth_beg;
% source.uz(nbr_of_sources_1 + 1: nbr_of_sources_2 + nbr_of_sources_1, 1:kgrid.Nt) = repelem(v_smooth, nbr_of_sources_2, 1);

pos = gaussmf(kgrid.t_array, [t_end/20 t_end/6]);
vel_gaussian = [0 diff(pos)/dt];

mask = nonzeros(mask)';

source.uz(find(mask == 2),:) = repelem(vel_gaussian,nbr_of_sources_2, 1);
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, Nz) = 1;
figure(44)
% plot(kgrid.t_array, v_smooth);
% legend('Velocity Stimulation of Actuator');
% smooth the source
% source.uz = filterTimeSeries(kgrid, medium, source.uz);
%% define the sensor
% define a series of Cartesian points to collect the data
% x = [-Nx/2+2:2:Nx/2-2] * dx;     % [m]
% y = x;                           % [m]
% z = ones(size(x)) * 0.5* (gel_surface_z * dz); % [m]
% sensor.mask = [x; y; z];
sensor.mask = zeros(Nx, Ny, Nz);
% sensor.mask(2:1:Nx-2, 2:1:Ny-2, gel_surface_z-1) = 1;
sensor.mask = gel;
%% plot gel, source and sensor
voxelPlot(gel);
start_ = [Nx/2 - source_radius, Ny/2 - source_radius, gel_surface_z];
size_ = [2*source_radius,2*source_radius,1];
voxel(start_, size_, 'b', 1);
% voxelPlot(sensor.mask);

%% run the simulation 
sensor.record = {'u'} % u particle velocity
% PML
% lambda=c0/f0;    % wavelength in lighter product (m)
% h=lambda/6;                 % desired spatial step (m)
PML_size= 5;             % size of the PML in grid points
input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false};
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});