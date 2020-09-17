

function [P] = phantom(x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz, features)

P = zeros(Nx, Ny, Nz);  %% The size in X Y and Z direction

x = -x_cent*dx+(0:Nx-1)*dx;
y = -y_cent*dy+(0:Ny-1)*dy;
z = -z_cent*dz+(0:Nz-1)*dz;
[x,y,z] = meshgrid(x,y,z);

for k=1:length(features(:,1))
    
    x_  = features(k,1);
    y_  = features(k,2);
    z_  = features(k,3);
    a   = features(k,4);
    b   = features(k,5);
    c   = features(k,6);
    theta= features(k,7)*pi/180;
    phi= features(k,8)*pi/180;
    psi= features(k,9)*pi/180;
    mu  = features(k,10);
    
    e1 = [cos(phi)*cos(theta) sin(phi)*cos(theta) -sin(theta)];
    e2 = [-sin(phi) cos(phi) 0];
    
    n1 = cos(psi)*e1+sin(psi)*e2;
    n2 = -sin(psi)*e1+cos(psi)*e2;
    n3 = [cos(phi)*sin(theta) sin(theta)*sin(phi) cos(theta)];
    
    a_sq = 1/a^2;
    b_sq = 1/b^2;
    c_sq = 1/c^2;
    
    p1 = n1(1)^2*a_sq+n2(1)^2*b_sq+n3(1)^2*c_sq;
    p2 = n1(2)^2*a_sq+n2(2)^2*b_sq+n3(2)^2*c_sq;
    p3 = n1(3)^2*a_sq+n2(3)^2*b_sq+n3(3)^2*c_sq;
    p4 = n1(1)*n1(2)*a_sq+n2(1)*n2(2)*b_sq+n3(1)*n3(2)*c_sq;
    p5 = n1(1)*n1(3)*a_sq+n2(1)*n2(3)*b_sq+n3(1)*n3(3)*c_sq;
    p6 = n1(3)*n1(2)*a_sq+n2(3)*n2(2)*b_sq+n3(3)*n3(2)*c_sq;
    
    eq1 = p1*(x-x_).^2+p2*(y-y_).^2+p3*(z-z_).^2;
    eq2 = p4*(x-x_).*(y-y_)+p5*(x-x_).*(z-z_)+p6*(y-y_).*(z-z_);
    eq = eq1+2*eq2;
    
    ell    = zeros(size(eq));
    ell(eq<1.0) = 1;
 
    P      = P+mu*ell;
    
end
        
    
