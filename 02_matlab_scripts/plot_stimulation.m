stimulation_z = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
tp = find(kgrid.t_array == (0.5/source_freq));
sources_v = source_mag * cos(2 * pi * source_freq * kgrid.t_array) * 2 * pi * source_freq;
smooth_beg = 0.5*max(sources_v) * sin(2*pi*source_freq*kgrid.t_array + pi/8*source_freq) + max(sources_v)*0.5;
smooth_beg(tp:end) = max(sources_v);
v_smooth = sources_v .* smooth_beg;
plot(kgrid.t_array, v_smooth);
legend('Acceleration Stimulation of Actuator');