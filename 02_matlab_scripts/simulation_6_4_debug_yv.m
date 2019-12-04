
% Quadratic gel simulation (no air)
clc;
%% parameter
% Actuator input:
sinusoidal = 0;
gaussian = 1;
% Boundary conditions:
reflecting = 0;
bc_vel_bottom = 0; % first two layer vel = 0
bc_free_surface = 0; % stress on the boundary = 0
c_compression = 25.0; %artificially low compression wave speed
c_shear = 4.0; % shear wave speed in gelatin 
f_max = 300; %[Hz]
%%  create the computational grid
% diameter: 15cm, thickness: 1.4cm
spacing = 2.5e-3;
x_size = 0.15;            %[m] total size x-direction
y_size = 0.15;            %[m] total size y-direction
z_size = 0.014;           %[m] total size z-direction
dx = spacing;             % grid point spacing in the x direction [m]
dy = spacing;             % grid point spacing in the y direction [m]
dz = spacing;             % grid point spacing in the z direction [m]
Nx = ceil(x_size/dx)+2;   % number of grid points in the x direction
Ny = ceil(y_size/dy)+2;   % number of grid points in the y direction
Nz = ceil(z_size/dz)+2    % number of grid points in the z direction
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% simulation time and step size
dt = 0.1 * spacing / c_compression; %1440.3; % CFL*dx/c_max CFL 0.3 default
%dt = 0.1 * spacing / 1440.3; % CFL*dx/c_max CFL 0.3 default
t_end = dt*4000;% 0.25/300; % dt*4000; % [s]
kgrid.t_array = 0:dt:t_end;
%% Geometry see in other scripts to input gel geometry
%% define the properties of the propagation medium
medium.sound_speed_compression = ones(Nx,Ny,Nz) * c_compression; %1440.3;    % [m/s]
medium.sound_speed_shear = ones(Nx,Ny,Nz) * c_shear; %201.74;          % [m/s]
medium.alpha_coeff_compression = ones(Nx,Ny,Nz) * .5;%0.5       % [10–3 m–1] 
medium.alpha_coeff_shear = ones(Nx,Ny,Nz) * .5;%500;             % [10–3 m–1] 
medium.density = ones(Nx,Ny,Nz) * 923.47;                    % [kg/m^3]                               
%% boundary conditions
if (bc_vel_bottom)
    source.u_mask = zeros(Nx, Ny, Nz);
    source.u_mask(:,:,1) = 1; % HEIGHT = 2 of bottom to be zero
    nbr_of_sources_1 = sum(source.u_mask, "all");
    source.uz(1:nbr_of_sources_1, 1:kgrid.Nt) = zeros(nbr_of_sources_1, kgrid.Nt);
else
    nbr_of_sources_1 = 0;
    source.u_mask = zeros(Nx,Ny,Nz);
end
if (bc_free_surface)
    gel_outer3d = ones(Nx,Ny,Nz);
    gel_outer3d(2:Nx-1, 2:Ny-1, 2:Nz-1) = 0;
    source.s_mask   = gel_outer3d;
    source.sxx      = zeros(1, kgrid.Nt);
    source.syy      = source.sxx;
    source.sxy      = source.syy;
    source.s_mode   = 'dirichlet';
end
%% define the source
% define a square source element
source_radius = 1; %0.005/dx;  % [grid points]
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, end) = 2;
nbr_of_sources_2 = 0.5 * (sum(sum(sum(source.u_mask))) - nbr_of_sources_1);
mask = reshape(source.u_mask, 1, Nx*Ny*Nz);
mask = nonzeros(mask)';
if(sinusoidal)
    % define a time varying sinusoidal source
    source_freq = f_max; %2/(t_end);      % [Hz]
    source_mag = 0.001;           % [m]
    stimulation_z = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
    tp = find(kgrid.t_array == (0.5/source_freq));
    sources_v = source_mag * cos(2 * pi * source_freq * kgrid.t_array) * 2 * pi * source_freq;
    smooth_beg = 0.5* sin(2*pi*source_freq*kgrid.t_array + pi/8*source_freq) + 0.5;
    smooth_beg(tp:end) = 1;
    v_smooth = sources_v .* smooth_beg;
%     source.uz(nbr_of_sources_1 + 1: nbr_of_sources_2 + nbr_of_sources_1, 1:kgrid.Nt) = repelem(v_smooth, nbr_of_sources_2, 1);
    source.uz(find(mask == 2),:) = repelem(vel_gaussian,nbr_of_sources_2, 1);
end
if(gaussian)
    pos = 0.0003 * gaussmf(kgrid.t_array, [t_end/20 t_end/6]);
%    pos = 0.001 * gaussmf(kgrid.t_array, [t_end/20 t_end/6]);
    vel_gaussian = [0 diff(pos)./dt];
    source.uz(find(mask == 2),:) = repelem(vel_gaussian,nbr_of_sources_2, 1);
    figure(11);
    subplot(2,1,1);
    plot(kgrid.t_array, source.uz(end,:));
    title("Actuator input velocity");
    subplot(2,1,2);
    plot(kgrid.t_array, pos);
    title("Actuator input position");    
end
source.u_mask(find(source.u_mask == 2)) = 1;
source.u_mode= 'dirichlet';
%% define the sensor
sensor.mask = ones(Nx,Ny,Nz);
%% plot gel, source and sensor
voxelPlot(ones(Nx,Ny,Nz));
start_ = [Nx/2 - source_radius, Ny/2 - source_radius, Nz];
size_ = [2*source_radius,2*source_radius,1];
voxel(start_, size_, 'b', 1);
%% run the simulation 
sensor.record = {'u'} % u particle velocity
% PML
% lambda=c0/f0;    % wavelength in lighter product (m)
% h=lambda/6;                 % desired spatial step (m)
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