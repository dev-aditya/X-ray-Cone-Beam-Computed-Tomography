function proj = projections(u_cent, v_cent, du, dv, Nu, Nv, Np, src_obj, det_src, feat)

proj  = zeros(Nv, Nu, Np);

u_1=((0:Nu-1) - u_cent) * du;
v_1=((0:Nv-1) - v_cent) * dv;

[u , v]=meshgrid(u_1,v_1);

dlambda = 2*pi/(Np);
lambda=dlambda*(0:1:Np-1);

for i = 1:length(lambda)
    lam = lambda(i);
    
    for m = 1:size(feat,1)
        
        
        x  = feat(m,1);
        y  = feat(m,2);
        z  = feat(m,3);
        a   = feat(m,4);
        b   = feat(m,5);
        c   = feat(m,6);
        theta = feat(m,7)*pi/180;
        phi   = feat(m,8)*pi/180;
        psi = feat(m,9)*pi/180;
        mu  = feat(m,10);
        
        
        e1 = [cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta)];
        e2 = [-sin(phi), cos(phi), 0];
        
        alpha = zeros(Nv, Nu, 3);
        
        n1 = cos(psi)*e1+sin(psi)*e2;
        n2 = -sin(psi)*e1+cos(psi)*e2;
        n3 = [cos(phi)*sin(theta), sin(theta)*sin(phi), cos(theta)];
        
        a_sq = 1/a.^2;
        b_sq = 1/b.^2;
        c_sq = 1/c.^2;
      
        
        a = diag([a_sq, b_sq, c_sq]);
        
        q = [n1;n2;n3];
        
        gamma = q'*a*q;
        
        ew = [cos(lam) sin(lam) 0];
        eu = [-sin(lam) cos(lam) 0];
        ev = [0 0 1];
        for k = 1:3
            U=(u.*eu(k) + v*ev(k) - det_src*ew(k));
            alpha(:,:,k) = 1./sqrt(u.^2 + v.^2 + det_src.^2) .*U;
        end
        
        al1 = alpha(:,:,1);
        al2 = alpha(:,:,2);
        al3 = alpha(:,:,3);
        
        ahat = [src_obj*cos(lam), src_obj*sin(lam), 0];
        xhat = ahat - [x y z];
        
        el_0 = al1 .* ( al1.*gamma(1,1) + al2.*gamma(1,2) + al3.*gamma(1,3) ) + ...
               al2 .* ( al1.*gamma(2,1) + al2.*gamma(2,2) + al3.*gamma(2,3) ) + ...
               al3 .* ( al1.*gamma(3,1) + al2.*gamma(3,2) + al3.*gamma(3,3) );
        
        el_1 = al1 .* ( xhat(1).*gamma(1,1) + xhat(2).*gamma(1,2) + xhat(3).*gamma(1,3) ) + ...
               al2 .* ( xhat(1).*gamma(2,1) + xhat(2).*gamma(2,2) + xhat(3).*gamma(2,3) ) + ...
               al3 .* ( xhat(1).*gamma(3,1) + xhat(2).*gamma(3,2) + xhat(3).*gamma(3,3) );
        
        el2 = (xhat * gamma * xhat') - 1;
        el=1./el_0;
        
        t_u = (-el_1 - real(sqrt(el_1.^2-el_0.*el2))).*el;
        t_v = (-el_1 + real(sqrt(el_1.^2-el_0.*el2))).*el;
        
        proj(:,:,i) = proj(:,:,i) + mu.*(t_v-t_u);
        
        
    end
end