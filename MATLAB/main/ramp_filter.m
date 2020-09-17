function filtered_back = ramp_filter(proj,u_cent,v_cent,du,dv,src_obj,src_det)

[Nv,Nu,Nw]=size(proj);
filtered_back = zeros(size(proj));
g = zeros(size(1,256));
for deg=1:Nw
    for i=1:Nu
        for j=1:Nv
            
            u = (i-1-u_cent)*du;
            v = (j-1-v_cent)*dv;
            
            
            for k=1:Nu
                
                u_1 = du*(k-1-u_cent);
                dist  = src_det/sqrt(u_1^2+v^2+src_det^2)*du;
                
                t_1 = pi*(u-u_1)/du;
                t_2= t_1+pi;
                t_3= t_1-pi;
                
                if abs(t_1) < 1e-4
                    r_1 = 0.5;
                else
                    r_1 =sin(t_1)/t_1+(cos(t_1)-1)/(t_1^2);
                end
                
                if abs(t_2) < 1e-4
                    r_2 = 0.5;
                else
                    
                  r_2 = sin(t_2)/t_2+(cos(t_2)-1)/(t_2^2);
                end
                
                if abs(t_3) < 1e-4
                    r_3 = 0.5;
                else
                    r_3 = sin(t_3)/t_3+(cos(t_3)-1)/(t_3^2);
                end
                
                ram = src_obj*r_1+0.5*(1-src_obj)*(r_2+r_3);
               
                ramp = 0.5/du^2*ram;
                
                g(k) = ramp*dist*proj(j,k,deg);
                
            end
            
            filtered_back(j,i,deg) = sum(g);
            
        end
    end
end
end



