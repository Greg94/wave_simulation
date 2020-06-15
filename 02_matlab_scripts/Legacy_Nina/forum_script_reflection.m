PML_size = 20; % [grid points]
Nx = 168 - 2*PML_size; % [grid points]
Ny = 168 - 2*PML_size; % [grid points]

dx = 0.1e-4; % grid point spacing in the x direction [m]
dy = 0.1e-4; % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed_compression = 10; % [m/s]
medium.sound_speed_shear = 3.5; % [m/s]
medium.density = 1000; % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.03; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear = 0.999; % [dB/(MHz^2 cm)]

% define a single source point
source.u_mask = zeros(Nx, Ny);
% source.u_mask(5, Ny/2) = 1;
source.u_mask(1:Nx/2, Ny/2-10:Ny/2+11) = 1;

% Defining Guassian source
sig= 5e-4;
xx= kgrid.y_vec(Ny/2-10:Ny/2+11)';
FGaus =zeros(Nx/2,Ny/2+11-(Ny/2-10)+1);
FGaus=repmat((-1/2/sig^2)*exp(-(xx./(2*sig)).^2),Nx/2,1);
% FGaus(:,Ny/2-10:Ny/2+10)=repmat((-1/2/sig^2)*exp(-(xx./(2*sig)).^2),Nx,1)

ux_F = FGaus.*kgrid.dx./medium.sound_speed_compression;
ux=reshape(ux_F,[],1);
source.ux=ux;

sensor.mask = zeros(Nx, Ny);
% sensor.mask(Nx/2-10:Nx/2+10, Ny/2) = 1;
sensor.mask(Nx/2, Ny/2-40:Ny/2+40) = 1;
% sensor.mask(pos2, Ny/2) = 1;
sensor.record = {'u','p',};

% set the CFL
cfl = 0.1;

% define the properties of the PML to allow plane wave propagation
pml_alpha = 2;
% pml_size = [10, 2];

% set the input arguments
% input_args = {'PlotScale', 'auto', 'PMLSize', pml_size,...
% 'PMLAlpha', pml_alpha, 'PlotPML', false,'PMLInside',false};
input_args = {'PlotScale', 'auto','PlotPML', false,'Smooth',true...
'PMLInside',false,'PMLAlpha', pml_alpha, 'RecordMovie',true,'DataCast', 'single'};

% set end time
t_end = 0.3e-3;

% create the time array
kgrid.makeTime(max(medium.sound_speed_compression, medium.sound_speed_shear), cfl, t_end);

% run the simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
sensor_data_reordered.p = reorderSensorData(kgrid, sensor, sensor_data.p);
sensor_data_reordered.ux = reorderSensorData(kgrid, sensor, sensor_data.ux);
sensor_data_reordered.uy = reorderSensorData(kgrid, sensor, sensor_data.uy);
%
figure;
imagesc(sensor_data_reordered.p,[-2000,2000]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;