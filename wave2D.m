clc; clear;
clf

% =======================================
% Prepare the new file.
% file    = input('Nama File= ','s');
% vidObj  = VideoWriter(file);
% vidObj.Quality = 55;
% open(vidObj);
% =======================================

Lx = 10;
Ly = 10;
dx=0.1;
dy=dx;

nx=fix(Lx/dx);
ny=fix(Ly/dy);
x =linspace(0, Lx, nx);
y =linspace(0, Ly, ny);

T=6;

wn   = zeros(nx,ny);
wnm1 = wn;
wnp1 = wn;

CFL = 0.5;          %CFL = c.dt/dx
c   = 1;
dt  = CFL*dx/c;

% wn(49:51) = [0.1 0.2 0.1];
% wnp1(:) = wn(:);

t=0;

while(t<T*3/4)
    %     Reflecting Boundary Condition
        % wn(:,[1 end])=0;
        % wn([1 end],:)=0;
        % bound = 'Reflecting Boundary';

    %     Absorbsing boundary Conditions
    wnp1(1,:)   = wn(2,:) + ((CFL-1)/(CFL+1))*(wnp1(2,:) - wn(1,:));
    wnp1(end,:) = wn(end-1,:) + ((CFL-1)/(CFL+1)) * (wnp1(end-1,:)-wn(end,:));
    wnp1(:,1) = wn(:,2) + ((CFL-1)/(CFL+1))*(wnp1(:,2) - wn(:,1));
    % wnp1(:,1)   = 0;
    wnp1(:,end) = wn(:,end-1) + ((CFL-1)/(CFL+1)) * (wnp1(:,end-1)-wn(:,end));
    bound       = 'Combination Boundary';

    %       Iteration
    t  =t+dt;
    wnm1=wn; wn=wnp1;

    % source
    if t<=4/3 + 2/3
        wn(51,71) = -dt^2*6*sin((30*pi/20)*t);
    end

    for i=2:nx-1
        for j=2:ny-1
            wnp1(i,j) = 2*wn(i,j) - wnm1(i,j) ...
                        + CFL^2 * (wn(i+1,j) + wn(i,j+1) - 4*wn(i,j) ...
                        + wn(i-1,j) + wn(i,j-1));
        end
    end
    subplot(1,2,1)
    imagesc(x,y,wn'); colorbar; clim([-0.01 0.01]);
    set(gca,'ydir','normal');
    colormap hsv
    title(sprintf('With %s || Time = %.2f',bound,2*t/3),FontSize=15);


    subplot(1,2,2)
    surf(x, y, wn'); clim([-0.01 0.01]);
    axis([0 Lx 0 Ly -0.01 0.04]);
    colormap hsv
    title(sprintf('With %s || Time = %.2f',bound,2*t/3),FontSize=15);

    %     hold on
    %     plot(5,wn(50),'square',LineWidth=10,Color='red');
    %     ylabel('Amplitudo');
    %     xlabel('Distances');
    %     hold off

    grid on;
    drawnow

    % =======================================
    % Write each frame to the file.
    % currFrame = getframe(gcf);
    % writeVideo(vidObj,currFrame);
    % =======================================
end

% =======================================
% Close the file.
% close(vidObj);
% =======================================
