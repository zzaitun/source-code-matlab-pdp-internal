clc; clear;
clf

% =======================================
% Prepare the new file..
% file=input('Nama File= ','s');
% vidObj = VideoWriter(file);
% vidObj.Quality = 55;
% open(vidObj);
% =======================================

% Parameter
rho=0.9;
cp=1.9;
% K=0.4;
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
Lx=10;
Ly=10;

Nx=51;
Ny=51;
Nt=400;

dx=Lx/(Nx-1);
dy=Ly/(Ny-1);

c=1;
C=0.05;
dt=C*dx/c;

Tn=zeros(Ny,Nx);



% Panas
% K = ones(Ny,Nx);
% K(10:15, 20:45) = 0.00001;
% K(10:25, 40:45) = 0.00001;

% processor
K = ones(Ny,Nx);
K(19:33, 19:33) = 0.00001;
% K(10:25, 40:45) = 0.00001;


% Initial Condition
Tn(:,:)=0;
t=0;

name='Zaitun';

for n=1:Nt
    Tc = Tn;
    if (t*13/4 > 6)
        cp=2;
    end
    t = t + dt;
    for i=2:Nx-1
        for j=2:Ny-1
            Tn(j,i)= Tc(j,i) + ...
                dt * (K(j,i)/rho/cp) * ...
                ((Tc(j,i+1) + Tc(j+1,i)- 4*Tc(j,i) + ...
                Tc(j,i-1)+Tc(j-1,i))/(dx)^2);
        end
    end
%     Source
    % Sx = round(3*Nx/Lx);
    % Sy = round(3*Ny/Ly);
    fire="off";
    % Sources
    if (t*13/4 <= 6)
        fire="on";
        Tn(18:34,18) = Tn(18:34,18) + dt*500/rho/cp;
        Tn(34,34) = Tn(34,34) + dt*500/rho/cp;
        Tn(18:34,34) = Tn(18:34,34) + dt*500/rho/cp;
        Tn(34,18) = Tn(34,18) + dt*500/rho/cp;
        % Tn(Sy-7,Sx+32) = Tn(Sy,Sx) + dt*500/rho/cp;
        % Tn(Sy+15,Sx+27) = Tn(Sy+30,Sx+30) + dt*500/rho/cp;
        % Tn(Sy,Sx) = Tn(Sy,Sx) + dt*500/rho/cp;
        % Tn(Sy-7,Sx+32) = Tn(Sy,Sx) + dt*500/rho/cp;
        % Tn(Sy+15,Sx+27) = Tn(Sy+30,Sx+30) + dt*500/rho/cp;
    end
 
    % Boundary Conditions
    Tn(1,:) = Tn(2,:);
    Tn(end,:) = Tn(end-1,:);
    Tn(:,1) = Tn(:,2);
    Tn(:,end)=Tn(:,end-1);

    imagesc(x,y,Tn);
    fr=sprintf('Computer : %s || time = %.2f s',fire,t*13/4);
    title('Graph of Plate Temperature',fr,FontSize=15);
    xlabel('Length of Plate');
    ylabel('Width of Plate');
    axis([0 Lx 0 Ly 0 30]);
    set(gca,'ydir','normal');

    colorbar
    colorbar('Ticks',[0,1,2,3,4],...
        'TickLabels',{'Cool','Semi Cool','Neutral','Warm','Hot'});
    clim([0 4]);
    colormap jet
%     view([0 90]);

    drawnow
    % =======================================
    % Write each frame to the file..
        % currFrame = getframe(gcf);
        % writeVideo(vidObj,currFrame);
    % =======================================
end

% =======================================
% Close the file..
% close(vidObj);
% =======================================




