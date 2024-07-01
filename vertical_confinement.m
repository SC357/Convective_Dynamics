

% ************************************************************
% Simulating the participate advection equation in 2D using
% pseudo-spectral code and ETD2 exponential time-stepping 
% Periodic boundary conditions are used.
% ***********************************************************

disp('*** 2D SHE SIMULATION ***');

u= VideoWriter('hor_density.avi');
open(u); 
figure('position', [200 200 950 800])
axis equal
colormap('jet')

% Set equation coefficients
eps=-0.2;
nu=0.1; 
q0=1;
Dr=0.4; 
Dc=10.; 
gam=0.02; 
beta=12; 
fac=-0.5; 
kap1=0.05; 
kap2=0.; 


% Set system parameters
L    = 60d0;        % Domain size (assume square container)
Tmax = 50d0;        % End time
N    = 512;       % Number of grid points
dT   = 0.05d0;      % Timestep
dps  = 10000;       % Number of stored times
      % Whether to continue simulation from final values of previous


betat=beta*dT; 
dx=L/N; 
adT=dT/(4*dx*dx); 
% Calculate some further parameters
nmax  = round(Tmax/dT);
XX    = (L/N)*(-N/2:N/2-1)'; 
[X,Y] = meshgrid(XX);
nplt  = floor(nmax/dps);

% Define initial conditions
A = double( 10^(-2)*randn(size(X)));
C = 1.+double( 10^(-2)*randn(size(X)));
%Rho = 1.+double( 10^(-2)*randn(size(X)));
Rho=exp(-0.001*(X.^2+Y.^2)); 
for i=1:N 
    for j=1:N 
        if i<N/4 & i>30 Rho(i,j)=Rho(N/4, j) ; end 
    end 
end 
% Set wavenumbers.
%Rho = sin(4*pi*X./L).*sin(4*pi*Y./L)+1+double( 10^(-2)*randn(size(X)));
k  				= double([0:N/2-1 0 -N/2+1:-1]'*(2*pi/L));
k2 				= double(k.*k);  
k2(N/2+1) = ((N/2)*(2*pi/L))^2;
[k2x,k2y] = meshgrid(k2); 
del2 			=double( k2x+k2y);

 
% Compute exponentials 
expA 	  	= double(exp(dT*(eps-nu.*del2)));
cor=fac*adT./del2; 
cor(1,1)=0; 

expR=exp(-dT*Dr.*del2);
expC=exp(-dT*Dc.*del2);
expB=betat*(1-exp(-kap1.*del2-kap2.*del2.*del2));


% Solve PDE
T=0;
for n = 1:nmax
	T = T + dT;

%calcualted the kernel 
w=A- real(fftn(expB.*ifftn(C.*Rho))); 


%updates A

	A= real(exp(-0.5*dT*A.^2).*fftn(expA.*ifftn(w)));
% potential function    
    phi= real(fftn(cor.*ifftn(A)));

%velocities
    for i=1:N
            ip=i+1; if i==N ip=1; end 
            im=i-1; if i==1 im=N; end 
        for j=1:N
            jp=j+1; if j==N jp=1; end 
            jm=j-1; if j==1 jm=N; end

       vx(i,j)=phi(ip,j)-phi(im,j);
       vy(i,j)=phi(i,jp)-phi(i,jm);
        end
    end

    %advection terms 

 for i=1:N
            ip=i+1; if i==N ip=1; end 
            im=i-1; if i==1 im=N; end 
        for j=1:N
            jp=j+1; if j==N jp=1; end 
            jm=j-1; if j==1 jm=N; end
cad(i,j)=vx(ip,j)*C(ip,j)-vx(im,j)*C(im,j)+vy(i,jp)*C(i,jp)...
-vy(i,jm)*C(i,jm); 
rad(i,j)=vx(ip,j)*Rho(ip,j)-vx(im,j)*Rho(im,j)+vy(i,jp)*Rho(i,jp)...
-vy(i,jm)*Rho(i,jm); 
        end
    end

    C= real(exp(-dT*gam*Rho.*C).*fftn(expC.*ifftn(C-cad)));
	Rho= real(fftn(expR.*ifftn(Rho-rad)));

    % Saving data
	
	
	% Commenting on time elapsed
	if mod(n,floor(nmax/200)) == 0
		amax=max(max(A)); 
        cmin=min(min(C)); 
    outp = strcat('  T=', num2str(T), ' amax=',num2str(amax), ' cmin=',num2str(cmin)); disp(outp);
    
	hold all

% Plot evolution
v = [0,0];
pcolor(X,Y, Rho')
cb = colorbar('eastoutside', FontSize=12, FontWeight='bold');
cb.Label.String = '\rho';
hold on
mp=16;

 
 vx1(1:N/mp,1:N/mp)=vx(1:mp:N,1:mp:N);
 vy1(1:N/mp,1:N/mp)=vy(1:mp:N,1:mp:N);
 [x1,y1]=meshgrid(0:mp:N-1, 0:mp:N-1);
%quiver(L*[-N/2:mp:N/2-0.5]/N, L*[-N/2:mp:N/2-0.5]/N,vx1',vy1',0.5,'w','linewidth',2)
% curvvec(L*x1/N-L/2,L*y1/N-L/2,vy1,vx1,'linewidth',1.7,'numpoints',30)
%hold on


%hold on 
%contour(X,Y, real(A),v,"r")
%hold on 
%contour(X,Y, imag(A),v,"b")
view(0,90), shading interp, axis tight
set(gca,'position', [0 0 1 1])
set(gca,...
	'xcolor',		[0.6 0.6 0.6],...
	'ycolor',		[0.6 0.6 0.6],...
	'fontsize',	6,...
	'fontname', 'courier')
xlabel('X',...
	'fontname', 'courier',...
	'fontsize', 6,...
	'color',		[0.6 0.6 0.6])
ylabel('Y',...
	'fontname', 'courier',...
	'fontsize', 6,...
	'color',		[0.6 0.6 0.6],...
	'rotation', 0)
title('final configuration of |A|',...
	'fontname', 'courier',...
	'fontsize', 6,...
	'color', 		[0.6 0.6 0.6])

drawnow

%hold on 
frame=getframe(gcf); 
writeVideo(u,frame); 
hold off 

end
end

close(u) 
