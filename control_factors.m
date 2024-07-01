

% ************************************************************
% Simulating coupled equation for convection and particulate transport 
% ***********************************************************


arrowThickness = 1; % Change this value to obtain different arrow thicknesses


disp('*** 2D chem convection equations ***');
velocityx=[];
velocityz=[];
com_x_total = [];
com_y_total = [];
for factor = 0.6:0.2:1.2
% figure('position', [200 200 1200 (h/4)*1200])
figure('position', [200 200 1350 300])
vidfile = VideoWriter([strcat('conv_video_chemfield',num2str(factor),'video.mp4')],'MPEG-4');
open(vidfile);



% Set equation coefficients
nu=0.4   % viscosity 
Dc=0.1  % fuel diffusion 
Drho=0.05 % particulate diffusion 
gamma=0.3   % fuel depletion rate 
epsilon=5  % effective buoyancy 
g= 0.5    % scaled gravity 
vs=0.02 % sedimentation rate 



% Set system parameters
% Lz    = (h/4)*4*40d0;     %domain height 
Lz    = 40d0;     %domain height 
Lx    = 4*40d0; % Domain size (assume rectangular  container)
Tmax = 25d0;    % End time
Nx    = 4*256;      % Number of grid points in x-direction 
Nz    = 256;      % Number of grid points in z-direction 
dT   = 0.0005d0;     % Timestep
dps  = 200;         % Number of stored times



% simulation constants 
dx=Lx/Nx; 
dz=Lz/Nz; 
dx2=2*dx; 
dz2=2*dz; 
dTdx=nu*dT/dx^2; 
dTdz=nu*dT/dz^2;
dTdxc=Dc*dT/dx^2; 
dTdzc=Dc*dT/dz^2;
dTdxr=Drho*dT/dx^2; 
dTdzr=Drho*dT/dz^2;
gammaT=dT*gamma; 
epsT=epsilon*dT; 
gT=g*dT;
vsT=vs*dT; 

T=0;     

% Calculate some further parameters
nmax  = round(Tmax/dT);
mp=32;   % contol number of arrows 
XX    = (Lx/Nx)*(-Nx/2:Nx/2-1); 
ZZ    = (Lz/Nz)*(-Nz/2:Nz/2-1);
[X,Z] = meshgrid(XX,ZZ);
R=4./cosh(0.2*sqrt((X').*(X')+(Z'+Lz/2-10).*(Z'+Lz/2-10))); % generates rho profile 
% R=R*factor;
nplt  = floor(nmax/dps);
cor=double(zeros(Nx,2*Nz));       % kernel of fft 
omega1=double(zeros(Nx,2*Nz));
psi1=double(zeros(Nx,2*Nz)); 

% Define initial conditions
psi = double(zeros(Nx,Nz) + 10^(-2)*randn(Nx,Nz)); % stream fucntion 
omega=double(zeros(Nx,Nz)+10^(-1)*randn(Nx,Nz))+ 0*sin(2*pi/Lx*X').*sin(pi/Lz*Z'); % vorticity 
omegan=omega; % updated omega 
c = double(ones(Nx,Nz));    % chem concentration 

c = c*factor;

% c = double(zeros(Nx,Nz));

% Gradient chemical field
% factor = -1;
% weight = 1;
% for i = 1:4*256
%     c(i,1:end) = c(i,1:end)+weight*factor;
%     factor = -1 + i*(2-2/(4*256))/(4*256);
% end
rho=double(zeros(Nx,Nz))+R;   % particulte density
% rho=double(zeros(Nx,Nz));
cn = c;    % chem new concentration
rhon=rho;  % density new

% Set wavenumbers. size doubled in z-direction to implement zero b.c in fft
k0x2= (2*pi/Lx)^2;
k0z2=(pi/Lz)^2;
for m=0:Nx-1
    kx2= k0x2*m^2;
            if m>Nx/2-1
                 kx2= k0x2*(Nx-m)^2; end;
    for n=0:2*Nz-1
    kz2= k0z2*n^2;
            if n>Nz-1                                                                                                                                                                                                                   
                  kz2= k0z2*(2*Nz-n)^2; end; 
              cor(m+1,n+1)= -1d0/(kx2+kz2); %/(2*Nx*Nz); 
    end 
end 
cor(1,1) = 0. ;



[nrows,ncols] = size(rho);
[Yw,Xw] = meshgrid(1:ncols,1:nrows);


 % Solve PDE
velox=[];
veloz=[];
com_x_separate=[];
com_y_separate=[];
density_projections=[];
distances=[];
old_distance=256;
for k = 1:nmax
T=T+dT;

	
  %symmetrization for the fft 
for n=1:Nz 
omega1(1:Nx,n)=omega(1:Nx,n); 
omega1(1:Nx,n+Nz)=-omega(1:Nx,Nz-n+1); 
end 
%calculating stream function 
psi1= -fftn(cor.*ifftn(omega1));
for n=1:Nz 
psi(1:Nx,n)=real(psi1(1:Nx,n)); 
end 
% calculating velocities 

for m=1:Nx 
    ip=m+1; if ip==Nx+1 ip=1; end; 
    im=m-1; if im==0 im=Nx; end; 
    for n= 2:Nz-1 
    jp=n+1; jm=n-1; 
    vz(m,n)= (psi(ip,n)-psi(im,n))/dx2; % vz = \partial psi/\partial  x 
    vx(m,n)=-(psi(m,jp)-psi(m,jm))/dz2;  % vx= - \partial psi/\partial z
    end
% non-slip conditions in z 
    n=1; jp=2; 
    vz(m,n)= (psi(ip,n)-psi(im,n))/dx2; 
    vx(m,n)=-(psi(m,jp)-psi(m,n))/dz2;  
    n=Nz;  jm=Nz-1; 
    vz(m,n)= (psi(ip,n)-psi(im,n))/dx2; 
    vx(m,n)=-(psi(m,n)-psi(m,jm))/dz2;
     end 
 % calculate gradients of omega, c. rho  



for m=1:Nx 
    ip=m+1; if ip==Nx+1 ip=1; end; 
    im=m-1; if im==0 im=Nx; end; 
    for n= 2:Nz-1 
    jp=n+1; jm=n-1; 
    omegax(m,n)= (omega(ip,n)-omega(im,n))/dx2; % omega_x 
    omegaz(m,n)= (omega(m,jp)-omega(m,jm))/dz2; % omega_z 
    rhox(m,n) = (rho(ip,n)-rho(im,n))/dx2 ; % gradient of rho 
    rhoz(m,n) = (rho(m,jp)-rho(m,jm))/dz2 ;
    rhoxc(m,n)=(rho(ip,n)*c(ip,n)-rho(im,n)*c(im,n))/dx2;
    rhoxv(m,n) = (vx(ip,n)*rho(ip,n)-vx(im,n)*rho(im,n))/dx2;
    rhozv(m,n) = (vz(m,jp)*rho(m,jp)-vz(m,jm)*rho(m,jm))/dz2;
    cvx(m,n) = (vx(ip,n)*c(ip,n)-vx(im,n)*c(im,n))/dx2 ; % gradient of c
    cvz(m,n) = (vz(m,jp)*c(m,jp)-vz(m,jm)*c(m,jm))/dz2 ;
    end 
    n=1; jp=2;
    omegax(m,n)=0; omegaz(m,n)=0; 
    rhox(m,n) = (rho(ip,n)-rho(im,n))/dx2 ;
    rhoxv(m,n) = (vx(ip,n)*rho(ip,n)-vx(im,n)*rho(im,n))/dx2; 
    rhozv(m,n) = (vz(m,jp)*rho(m,jp)+vz(m,n)*rho(m,n))/dz2;
    rhoxc(m,n)=(rho(ip,n)*c(ip,n)-rho(im,n)*c(im,n))/dx2;
    rhoz(m,n) = -g*rho(m,n)/Drho;
    cvx(m,n) = (vx(ip,n)*c(ip,n)-vx(im,n)*c(im,n))/dx2; 
    cvz(m,n) = (vz(m,jp)*c(m,jp)+vz(m,n)*c(m,n))/dz2 ;
    n=Nz; jm=Nz-1; 
    omegax(m,n)=0; omegaz(m,n)=0; 
    rhox(m,n) = (rho(ip,n)-rho(im,n))/dx2 ;
    rhoxc(m,n)=(rho(ip,n)*c(ip,n)-rho(im,n)*c(im,n))/dx2;
    rhoz(m,n) =  -g*rho(m,n)/Drho;
    rhoxv(m,n) = (vx(ip,n)*rho(ip,n)-vx(im,n)*rho(im,n))/dx2; 
    rhozv(m,n) = (-vz(m,n)*rho(m,n)-vz(m,jm)*rho(m,jm))/dz2;
    cvx(m,n) = (vx(ip,n)*c(ip,n)-vx(im,n)*c(im,n))/dx2 ; % gradient of c
    cvz(m,n) = (-vz(m,n)*c(m,n)-vz(m,jm)*c(m,jm))/dz2 ;
    end




%updates omega, c,  rho  

for m=1:Nx 
    ip=m+1; if m==Nx ip=1; end; 
    im=m-1; if m==1 im=Nx; end; 
    for n= 2:Nz-1 
    jp=n+1; jm=n-1; 
    omegan(m,n)=omega(m,n)+dTdx*(omega(ip,n)+omega(im,n)-2*omega(m,n))...
    +dTdz*(omega(m,jp)+omega(m,jm)-2*omega(m,n))-dT*vx(m,n)*omegax(m,n)-...
    dT*vz(m,n)*omegaz(m,n)-epsT*rhoxc(m,n)+gT*rhox(m,n);

    cn(m,n)=c(m,n)+dTdxc*(c(ip,n)+c(im,n)-2*c(m,n))...
    +dTdzc*(c(m,jp)+c(m,jm)-2*c(m,n))-dT*cvx(m,n)-...
    dT*cvz(m,n)-gammaT*rho(m,n)*c(m,n);

    rhon(m,n)=rho(m,n)+dTdxr*(rho(ip,n)+rho(im,n)-2*rho(m,n))...
    +dTdzr*(rho(m,jp)+rho(m,jm)-2*rho(m,n))-dT*rhoxv(m,n)-dT*rhozv(m,n)...
    +vsT*rhoz(m,n);



    end
    % vorticity, c, rho  b.c. 
    n=1;jp=2; 
    om0= -(8*psi(m,1)-psi(m,2))/dz2/2; 
    omegan(m,n)=omega(m,n)+dTdx*(omega(ip,n)+omega(im,n)-2*omega(m,n))+...
    dTdz*(omega(m,jp)-2*omega(m,n)+om0)-epsT*rhoxc(m,n)...
    +gT*rhox(m,n); 

    cn(m,n)=c(m,n)+dTdxc*(c(ip,n)+c(im,n)-2*c(m,n))...
    +dTdzc*(c(m,jp)+c(m,jp)-2*c(m,n))-gammaT*rho(m,n)*c(m,n)...
    -dT*cvx(m,n)-dT*cvz(m,n); 

     fac=(1d0+dz*vs/(2*Drho))/(1d0-dz*vs/(2*Drho));
    rm = fac*rho(m,n); 
    rhon(m,n)=rho(m,n)+dTdxr*(rho(ip,n)+rho(im,n)-2*rho(m,n))...
    +dTdzr*(rho(m,jp)+rm-2*rho(m,n))+vsT*(rho(m,jp)-rm)/dz2 ...
    -dT*rhoxv(m,n)-dT*rhozv(m,n); 





     n=Nz;jm=Nz-1;
    omN= -(8*psi(m,n)-psi(m,jm))/dz2/2; 
    omegan(m,n)=omega(m,n)+dTdx*(omega(ip,n)+omega(im,n)-2*omega(m,n))+...
    dTdz*(omega(m,jm)-2*omega(m,n)+omN)-epsT*rhoxc(m,n)...
    +gT*rhox(m,n); 

    cn(m,n)=c(m,n)+dTdxc*(c(ip,n)+c(im,n)-2*c(m,n))...
    +dTdzc*(c(m,jm)+c(m,jm)-2*c(m,n))-gammaT*rho(m,n)*c(m,n)...
    -dT*cvx(m,n)-dT*cvz(m,n); 

    rp = rho(m,n)/fac; 
    rhon(m,n)=rho(m,n)+dTdxr*(rho(ip,n)+rho(im,n)-2*rho(m,n))...
    +dTdzr*(rp+rho(m,jm)-2*rho(m,n))+vsT*(rp-rho(m,jm))/dz2...
    -dT*rhoxv(m,n)-dT*rhozv(m,n); 

end 

  omega=omegan; 
  c=cn; 
  rho=rhon; 

  total_weight = sum(rho(:));
  weighted_x_sum = sum(rho(:).*Xw(:));
  weighted_y_sum = sum(rho(:).*Yw(:));
  com_x = weighted_x_sum/total_weight;
  com_y = weighted_y_sum/total_weight;
	
	% Commenting on time elapsed
	if mod(k,floor(nmax/dps)) == 0
		outp = strcat('  T= ', num2str(T), ' completed'); disp(outp);
    omegamax=max(max(omega))
    cmin=min(min(c))
    totalrho=sum(sum(rho))

    com_x_separate = [com_x_separate,com_x];
    com_y_separate = [com_y_separate,com_y];
    
	hold all

% Plot evolution
 
%v = [0,0];
%pcolor(X,Z, omega')
 pcolor(X,Z, rho')
%  pcolor(X,Z,c')
 pbaspect([4 1 1])
 cb = colorbar('eastoutside', FontSize=12, FontWeight='bold')
 cb.Label.String = '\rho';

 % Binarize image to calculate velocity of the front
 binary = flip(rho>1.4,2);
 [op,I] = max(binary,[],2);
 distance = min(I(I>1));
 distances = [distances,distance];

%  velox = [velox,max(max(vx))];
 veloz = [veloz,old_distance-distance];
 old_distance = distance;
 density_projection = sum(rho,2); 
 density_projections = [density_projections, density_projection];

 vx1(1:Nx/mp,1:Nz/mp)=vx(1:mp:Nx,1:mp:Nz);
 vz1(1:Nx/mp,1:Nz/mp)=vz(1:mp:Nx,1:mp:Nz);

%  pcolor(X,Z, flip(double(binary), 2)')
%  pbaspect([4 1 1])

quiver(Lx*[-Nx/2:mp:Nx/2-1]/Nx, Lz*[-Nz/2+2.5:mp:Nz/2-1.5]/Nz,vx1',vz1','w', LineWidth=arrowThickness)
pbaspect([4 1 1])
hold on


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
title('Current  Configuration',...
	'fontname', 'courier',...
	'fontsize', 6,...
	'color', 		[0.6 0.6 0.6])

drawnow
hold on 
F = getframe(gcf); 
    writeVideo(vidfile,F);

end
end
close(vidfile)
velocityx=[velocityx,velox];
velocityz=[velocityz,veloz];
com_x_total = [com_x_total,com_x_separate];
com_y_total = [com_y_total,com_y_separate];

end


% figure(10)
t = Tmax/dps:Tmax/dps:Tmax;
% plot(t,velocityx(1,1:end))
% hold all
% plot(t,velocityx(2,1:end))
% hold all
% plot(t,velocityx(3,1:end))
% xlabel('time [s]')
% ylabel('Max velocity [m/s]')
% 
figure(11)
plot(t,velocityz(1,1:end))
hold all
plot(t,velocityz(2,1:end))
hold all
plot(t,velocityz(3,1:end))
xlabel('time [s]')
ylabel('Max velocity [m/s]')


figure(12)
plot(t,veloz(1,1:end))

xlabel('time [s]')
ylabel('Max velocity [m/s]')
