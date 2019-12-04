% create the computational grid 
Nx = 128;
Ny = 128;
Nz = 40;
dx = 0.1e-3; 
dy = 0.1e-3;
dz = 0.1e-3;
% [grid points] % [grid points] % [m] % [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);
dt = 1 * dx / 2000;
t_end = dt * 4000;
kgrid.t_array = 0:dt:t_end;
% define the compressional sound speed [m/s] 
medium.sound_speed_compression = 344*ones(Nx, Ny, Nz); 
medium.sound_speed_compression(Nx/2:end, :, :) = 1440.3;
% define the shear sound speed [m/s] 
medium.sound_speed_shear = zeros(Nx, Ny, Nz); 
medium.sound_speed_shear(Nx/2:end, :, :) = 202;
% define the mass density [kg/mˆ3] 
medium.density = 1.225*ones(Nx, Ny, Nz); 
medium.density(Nx/2:end, :, :) = 1000;
% define the absorption coefficients [dB/(MHzˆ2 cm)] 
medium.alpha_coeff_compression = 0.1; 
medium.alpha_coeff_shear = 0.5;
% define the initial pressure distribution 
disc_magnitude = 5; % [Pa] 
disc_x_pos = 40; 
disc_y_pos = 64; 
disc_radius = 5;
% [grid points] % [grid points] % [grid points]
%source.p0 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);
%% define the source
% define a square source element
source_radius = 0.005/dx;  % [grid points]
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(Nx/2 - source_radius:Nx/2 + source_radius, Ny/2 - source_radius:Ny/2 + source_radius, 5) = 1;
% define a time varying sinusoidal source
source_freq = 300; %2/(t_end);      % [Hz]
source_mag = 0.001;           % [m]
source.uz = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
% smooth the source
source.uz = filterTimeSeries(kgrid, medium, source.uz);
% define a circular binary sensor mask 
radius = 20;
% [grid points] 
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(2:10:Nx-2, 2:10:Ny-2, 5-1) = 1;
% run the simulation 
% sensor.record = {'u','p'}
sensor_data = pstdElastic3D(kgrid, medium, source, sensor);