function back_proj = back_project(proj, u_cent, v_cent, du, dv, src_obj, det_src, x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz)

[Nv, Nu, Np] = size(proj); 

u_1 = ((0:Nu-1) - u_cent)*du;
v_1 = ((0:Nv-1) - v_cent)*dv;
[u, v] = meshgrid(u_1, v_1);

dLambda = (2*pi)/Np;
lambda=dLambda*(0:1:Np-1)';

eu = [-sin(lambda), cos(lambda), zeros(length(lambda),1)];
ev = [zeros(length(lambda),1), zeros(length(lambda),1), ones(length(lambda),1)];
ew = [cos(lambda), sin(lambda), zeros(length(lambda),1)];

xc = ((0:Nx-1)-x_cent) * dx;
yc = ((0:Ny-1)-y_cent) * dy;
zc = ((0:Nz-1)-z_cent) * dz;

[x, y] = meshgrid(xc, yc);

back_proj = zeros(Ny, Nx, Nz); 

for i = 1:Nz
    z = zc(i);
    
    for j = 1:length(lambda)
        
        lam_ = lambda(j);
        lam_x = src_obj .* cos(lam_) .* ones(size(x)); 
        lam_y = src_obj .* sin(lam_) .* ones(size(y)); 
        u_num = (x-lam_x).*eu(j,1) + (y-lam_y).*eu(j,2) + z.*eu(j,3);
        u_den = (x-lam_x).*ew(j,1) + (y-lam_y).*ew(j,2) + z.*ew(j,3);
        u_ = -det_src.*(u_num./u_den);
        v_num = (x-lam_x).*ev(j,1) + (y-lam_y).*ev(j,2) + z.*ev(j,3);
        v_den = u_den;
        v_ = -det_src.*(v_num./v_den);
        G_ = interp2(u, v, proj(:,:,j), u_, v_, 'linear');
        G_(isnan(G_)) = 0;
        d = ( (lam_x-x).*ew(j,1) + (lam_y-y).*ew(j,2) + (0-z).*ew(j,3) ).^2;
        b = (src_obj.*det_src) ./ d;
        back_proj(:,:,i) = back_proj(:,:,i) + b.*G_ .* dLambda;
        
    end
    
end
back_proj = back_proj ./ 2;
end