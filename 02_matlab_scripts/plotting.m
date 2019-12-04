clf;
figure(1)
subplot(211)
title('Actuator');
% Stimulation
f = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
plot(kgrid.t_array,f);
% Response
subplot(212)
title('Response Acceleration');
u_z = sensor_data.uz(1,:);
t = linspace(0,t_end,t_end/dt+1);
% plot(t,u_z);
plot(kgrid.t_array, u_z);
hold off;