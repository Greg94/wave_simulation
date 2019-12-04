function [pos] = vel_to_pos(vel, dt)
pos = zeros(size(vel));
    for datapoint = 1:length(pos(:,1))
        for time = 2:length(pos(1,:))
            pos(datapoint, time) = pos(datapoint, time - 1) + dt * vel(datapoint, time);
        end
    end
end