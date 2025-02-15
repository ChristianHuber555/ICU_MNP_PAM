function Cav_Map = compute_PAM_RCB(RF, param, x, z, eps_RCB)
    RF = sign(RF).*sqrt(abs(RF));
    frames = size(RF,3);
    array_size = size(RF,2);
    
    N = size(RF, 2);
    I = ones(N);
    a = ones(N,1);

    Cav_Map = zeros(10,param.zdim,param.xdim);
    xdim = size(x,2);
    zdim = size(z,2);
    for i_frame = 1:frames
        parfor ix = 1:xdim
            xp = x(ix);
            for iz = 1:zdim
                zp = z(iz);
                dx = sqrt((param.xarray - xp).^2 + zp.^2);
                RF_b = zeros(size(RF,1),size(RF,2));
                for iarray = 1:array_size
                    t_new = param.t + dx(iarray)/param.c;
                    RF_b(:,iarray) = interp1(param.t,squeeze(RF(:,iarray,i_frame)),t_new,'linear',0);
                end
                R_cov = cov(RF_b);
                factor_R = 1e7 * trace(R_cov)/N;
                R_cov = R_cov .* factor_R;
                [U,V] = eig(R_cov);
    
                z_RCB = U'*a;

                lb = (norm(a) - sqrt(eps_RCB)) / (V(1,1)*sqrt(eps_RCB));
                valp = diag(V);
                lambda = 0.5*lb;

                % Find the optimal lambda
                while true
                    f = sum(abs(z_RCB).^2 ./ (1 + lambda * valp).^2)-eps_RCB;
                    fp = -2 * sum((valp .* abs(z_RCB).^2) ./ (1 + lambda * valp).^3);
                    if fp<=eps(1)
                        break
                    end
                    lambda = lambda - f / fp;
                end
                ahat = a - U * inv(I + lambda * V)*transpose(U)*a;
                ahat = ahat * N / norm(ahat);
                psi = (3*pi)/(1000 * param.c);
                psi_den = ahat' * U * V;
                psi_den_2 = inv(lambda^2*I + (2/lambda)*V +V.^2);
                psi_den_3 = U'*ahat;
                psi = psi/(psi_den * psi_den_2 * psi_den_3);
                Cav_Map(i_frame,iz,ix) = psi;
            end
            fprintf("Column %d from %d is finished", ix,param.xdim);
            fprintf('\n')
        end
        fprintf("RCB Image %d from %d is finished", i_frame,frames);
        fprintf('\n')
    end
end


