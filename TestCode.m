%% Raw Code
clc;
clear;
nx = 20;
ny = 20;
mu = 1;
rho = 1;
dx = 0.01;
dy = 0.01;
anu = zeros(ny+2,nx+1);
asu = zeros(ny+2,nx+1);
awu = zeros(ny+2,nx+1);
aeu = zeros(ny+2,nx+1);
anv = zeros(ny+1,nx+2);
asv = zeros(ny+1,nx+2);
awv = zeros(ny+1,nx+2);
aev = zeros(ny+1,nx+2);
uold = zeros(ny+2,nx+1);
vold = zeros(ny+1,nx+2);
us = zeros(ny+2,nx+1);
v = zeros(ny+1,nx+2);
P = zeros(ny+2,nx+2);
anp = zeros(ny+2,nx+2);
asp = zeros(ny+2,nx+2);
awp = zeros(ny+2,nx+2);
aep = zeros(ny+2,nx+2);
alphaU = 0.9;
alphaP = 0.4;
b = zeros(ny+2,nx+2);
pprime = zeros(ny+2,nx+2);
v = zeros(ny+1,nx+2);
u = zeros(ny+2,nx+1);
uold = zeros(ny+2,nx+1);
for iter = 1:1000
for i = 2:ny+1
    for j = 2:nx
        % Defining flux of each cell for u
        awu(i,j) = max(0.5*rho*(uold(i,j)+uold(i,j-1)),0)+ mu*dy/dx;
        aeu(i,j) = max(-0.5*rho*(uold(i,j+1)+uold(i,j)),0)+ mu*dy/dx;
        anu(i,j) = max(-0.5*rho*(vold(i-1,j+1)+ vold(i-1,j)),0)+ mu*dx/dy;
        asu(i,j) = max(0.5*rho*(vold(i,j+1)+vold(i,j)),0)+ mu*dx/dy;
        
        % Boundary Values
        
        
       
               
        % Set Corner Values as zero, not supposed to enter the computation.
        awu(1,1) =0; awu(end,1) = 0; awu(1,end) =0; awu(end,end) = 0;
        aeu(1,1) =0; aeu(end,1) = 0; aeu(1,end) =0; aeu(end,end) = 0;
        anu(1,1) =0; anu(end,1) = 0; anu(1,end) =0; anu(end,end) = 0;
        asu(1,1) =0; asu(end,1) = 0; asu(1,end) =0; asu(end,end) = 0;
        
    end

end
for i = 2:ny
    for j = 2:nx+1
        % Define flux at each cell for v
        awv(i,j) = max(0.5*rho*(uold(i-1,j-1)+uold(i,j-1)),0) + mu*dy/dx;
        aev(i,j) = max(-0.5*rho*(uold(i-1,j)+uold(i,j)),0) + mu*dy/dx;
        asv(i,j) = max(0.5*rho*(vold(i+1,j)+vold(i,j)),0) + mu*dx/dy;
        anv(i,j) = max(-0.5*rho*(vold(i,j)+vold(i-1,j)),0) +mu*dx/dy;
        
        % Boundary Values
        
       
        awv(1,1) =0; awv(end,1) = 0; awv(1,end) =0; awv(end,end) = 0;
        aev(1,1) =0; aev(end,1) = 0; aev(1,end) =0; aev(end,end) = 0;
        anv(1,1) =0; anv(end,1) = 0; anv(1,end) =0; anv(end,end) = 0;
        asv(1,1) =0; asv(end,1) = 0; asv(1,end) =0; asv(end,end) = 0;
        
    end
end
au = awu+aeu+anu+asu;
%au(3,2) = 1e+100;
av = awv+aev+anv+asv;
for veliter = 1:20
for i = 2:ny+1
    for j = 2:nx
        u(i,j) = (alphaU/au(i,j))*(awu(i,j)*u(i,j-1)+aeu(i,j)*u(i,j+1)+asu(i,j)*u(i+1,j)+anu(i,j)*u(i-1,j)+...
            P(i,j-1)-P(i,j)+((1-alphaU)*au(i,j)/alphaU)*uold(i,j));
        u(2,:) = 0;u(end-1,:) = 0;
        u(1,j) = -u(2,j);
        u(i,1) = 1;
        u(end,j) = -u(end-1,j);
        u(i,end) = 1;
        
        
    end
   
end
for i = 2:ny
    for j = 2:nx+1
        v(i,j) = (alphaU/av(i,j))*(awv(i,j)*v(i,j-1)+aev(i,j)*v(i,j+1)+asv(i,j)*v(i+1,j)+anv(i,j)*v(i-1,j)+...
            P(i+1,j)-P(i,j)+((1-alphaU)*av(i,j)/alphaU)*vold(i,j));
        v(end,:) = 0;
        v(1,:) = 0;
        v(2:ny+1,end) = -v(2:ny+1,end-1);
        v(i,1) = -v(i,2);
        
    end
        
 
end
end
%% Defining Pressure correction loop
for i = 2:ny+1
    for j = 2:nx+1
        aep(i,j) = alphaP*rho*((dx*dy)^2)*(1/au(i,j));
        awp(i,j) = alphaP*rho*((dx*dy)^2)*(1/au(i,j-1));
        asp(i,j) = alphaP*rho*((dx*dy)^2)*(1/av(i,j));
        anp(i,j) = alphaP*rho*((dx*dy)^2)*(1/av(i-1,j));
        % Source Term
        b(i,j) = rho*(dx*dy)*(us(i,j-1)-us(i,j)+ v(i,j)- v(i-1,j));
        aep(:,nx+1) = 0;
        awp(:,2) = 0;
        anp(2,:)=0;
        asp(ny+1,:) = 0;
        
    end
    %Boundary Values

end
ap = awp+asp+aep+anp;
ap(:,nx+1) = 1e+100;
for piter = 1:1100
for i = 2:nx+1
    for j = 2:nx+1
           pprime(i,j) = pprime(i,j) + (1/ap(i,j))*...
                (aep(i,j)*pprime(i,j+1) + awp(i,j)*pprime(i,j-1)...
                + anp(i,j)*pprime(i-1,j) + asp(i,j)*pprime(i+1,j)...
                +b(i,j)-pprime(i,j)*ap(i,j));
            
            
    end
end
end
%% Updating velocities
for i = 2:ny+1
    for j = 2:nx
        u(i,j) = u(i,j)+((dx*dy)/au(i,j))*(pprime(i,j-1)-pprime(i,j)); 
        u(i,j) = alphaU*u(i,j)+(1-alphaU)*uold(i,j);
        u(2,:) = 0;u(end-1,:) = 0;
        u(end,j) = u(end-1,j);
        u(1,j) = -u(2,j);
        u(2:ny+1,end) = 1;
        u(i,1) = 1; 
    end
end
for i = 2:ny
    for j = 2:nx+1
        v(i,j) = v(i,j) + ((dx*dy)/av(i,j))*(pprime(i,j)-pprime(i-1,j));
        v(2:ny+1,end) = -v(2:ny+1,end-1);
        v(i,1) = -v(i,2);
        v(1,:) = 0;
        v(end,:) = 0;
        
    end
end
for i = 2:ny+1
    for j = 2:nx+1
        P(i,j) = P(i,j)+alphaP*(pprime(i,j));
        P(i,end) = 0-P(i,end-1);
        P(i,1) = P(i,2);
        P(1,j) = P(2,j);
        P(end,j) = P(end-1,j);
    end
end
uold = u;
vold = v;
end