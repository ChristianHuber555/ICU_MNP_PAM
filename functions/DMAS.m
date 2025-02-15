function Cav_Map = compute_PAM_DMAS(RF, param, x, z)
    RF = sign(RF).*sqrt(abs(RF));
    frames = size(RF,3);
    array_size = size(RF,2);
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
                e1 = sum(RF_b,2);
                e2 = sum(RF_b.*RF_b,2);
                e1 = (e1.*e1-e2)/2;
                Cav_Map(i_frame,iz,ix) = sum(e1.*e1);
            end
            fprintf("Column %d from %d is finished", ix,param.xdim);
            fprintf('\n')
        end
        fprintf("DMAS Image %d from %d is finished", i_frame,frames);
        fprintf('\n')
    end
end