function [gel_outer] = get_boundary(Nx,Ny,gel_cond)
    gel_outer = zeros(Nx,Ny);
    for x = 1:Nx
        first = 0;
        for y = 1:Ny
            if (gel_cond(x,y) == 1 && first == 0)
                gel_outer(x,y:y+1) = 1;
                first = 1;
            elseif (gel_cond(x,y) == 0 && first > 0)
                gel_outer(x,y-2:y-1) = 1;
                break;
            end
        end
    end
    for x = 1:Nx
        first = 0;
        for y = 1:Ny
            if (gel_cond(y,x) == 1 && first == 0)
                gel_outer(y:y+1,x) = 1;
                first = 1;
            elseif (gel_cond(y,x) == 0 && first > 0)
                gel_outer(y-2:y-1,x) = 1;
                break;
            end
        end
    end
end