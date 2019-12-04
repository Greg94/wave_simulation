function [pos] = vel_to_pos_2d(vel, dt)
pos = zeros(size(vel));
    for x = 1:size(pos,1)
        for y = 1:size(pos,2)
            for time = 2:size(pos,3)
                pos(x, y , time) = pos(x, y, time - 1) + dt * vel(x, y, time);
            end
        end
    end
end