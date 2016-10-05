function [T,Y,DY] = do_nbody_leapfrog(t,x0,dx0,h,mus)
%Helper function for N-body leapfrog integration

n = length(mus); %this should be the number of bodies

%do the integration
[T,Y,DY] = leapfrog(@nbody_eq,t,x0,dx0,h);

    %second order nbody equation
    function ddx = nbody_eq(t,x)
        ddx = zeros(size(x));
        for j = 1:n
            indsj = (j-1)*3+1:(j-1)*3+3;
            for k = 1:n
                if k == j, continue; end
                indsk = (k-1)*3+1:(k-1)*3+3;
                % r_k/j = r_k/o - r_j/o
                rkj = x(indsk) - x(indsj);
                %d^2/dt^2(r_j/o) = sum_{k ~= j}(mu_k*r_k/j)/(r_k/j^3)
                ddx(indsj) = ddx(indsj)+mus(k)*rkj/norm(rkj)^3;
            end
        end
    end

end