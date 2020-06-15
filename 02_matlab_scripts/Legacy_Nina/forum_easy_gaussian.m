%% simulation of experimental platform (gel, piezo actuators)


% Quadratic gel simulation (no air), PML, no bc velocity or stress, uz = gaussian
% dirichlet/ no dirichlet


clc;

%% variables
f_max = 300; %[Hz]
%%  create the computational grid
spacing = 5e-3;
x_size = 0.15;            %[m] total size x-direction
y_size = 0.15;            %[m] total size y-direction
z_size = 0.014;           %[m] total size z-direction
dx = spacing;             % grid point spacing in the x direction [m]
dy = spacing;             % grid point spacing in the y direction [m]
dz = spacing;             % grid point spacing in the z direction [m]
Nx = ceil(x_size/dx)+2;   % number of grid points in the x direction
Ny = ceil(y_size/dy)+2;   % number of grid points in the y direction
Nz = ceil(z_size/dz)+2    % number of grid points in the z direction
disp('Number of grid points: '); disp(Nx*Ny*Nz);
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% simulation time and step size
dt = 0.1 * spacing / 1440.3;    % CFL*dx/c_max CFL 0.3 default
t_end = 0.5/300;                % [s]
kgrid.t_array = 0:dt:t_end;
%% define the properties of the propagation medium
medium.sound_speed_compression = ones(Nx,Ny,Nz) * 1440.3;    % [m/s]
medium.sound_speed_shear = ones(Nx,Ny,Nz) * 201.74;          % [m/s]
medium.alpha_coeff_compression = ones(Nx,Ny,Nz) * 0.5;       % [10–3 m–1] 
medium.alpha_coeff_shear = ones(Nx,Ny,Nz) * 0.5;             % [10–3 m–1] 
medium.density = ones(Nx,Ny,Nz) * 923.47;                    % [kg/m^3]                      
%% define the source
% define a square source element
source_radius = 0.005/dx;  % [grid points]
source.u_mask = zeros(Nx,Ny,Nz);
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, Nz) = 1; 
nbr_of_sources = sum(source.u_mask,"all");
% define a gaussian impulse
pos = gaussmf(kgrid.t_array, [t_end/20 t_end/6]);
vel_gaussian = [0 diff(pos)/dt];
source.uz = repelem(vel_gaussian,nbr_of_sources, 1);
source.u_mode= 'dirichlet';
%% define sensor
sensor.mask = ones(Nx, Ny, Nz);
%% run the simulation 
sensor.record = {'u'} % u particle velocity
PML_size=5;             % size of the PML in grid points
input_args = {'PMLSize', PML_size, 'PlotPML', true, 'PMLInside', false};
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});