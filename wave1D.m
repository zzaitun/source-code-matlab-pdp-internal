clc; clear;
clf
% =======================================
% Prepare the new file.
file=input('Nama File= ','s');
vidObj = VideoWriter(file);
vidObj.Quality = 55;
open(vidObj);
% =======================================
Lx = 10;
dx=0.1;
nx=fix(Lx/dx);
x=linspace(0, Lx, nx);

T=15;

wn=zeros(nx,1);
wnm1 = wn;
wnp1 = wn;

CFL = 1;
c=1;
dt=CFL*dx/c;

% wn(49:51) = [0.1 0.2 0.1];
% wnp1(:) = wn(:);

t=0;

while(t<45)
%     Reflecting Boundary Condition
    % wn(1)=0;
    % wn(end)=0;
    % bound = 'Reflecting Boundary';

%     Absorbsing boundary Conditions
        wnp1(1) = wn(2) + ((CFL-1)/(CFL+1))*(wnp1(2)-wn(1));
        wnp1(end) = wn(end-1) + ((CFL-1)/(CFL+1)) * (wnp1(end-1)-wn(end));
        bound = 'Absorbsing Boundary';


    t  =t+dt;
    wnm1=wn; wn=wnp1;
    % source
    wn(50) =dt^2*20*sin(10*pi*t/T);

    for i=2:nx-1
        wnp1(i) = 2*wn(i) - wnm1(i) ...
            + CFL^2 * (wn(i+1) - 2*wn(i) + wn(i-1));
    end
    plot(x,wn,LineWidth=2);
    hold on
    plot(5,wn(50),'square',LineWidth=10,Color='red');
    plot(0,wn(1),'square',LineWidth=10,Color='black');
    plot(x(end),wn(end),'square',LineWidth=10,Color='black');
    ylabel('Amplitudo');
    xlabel('Distances');
    hold off
    grid on;
    xline(0);
    title(sprintf('With %s || Time = %.2f',bound,t/3),FontSize=15);
    axis([0 Lx -0.5 0.5]);
    drawnow
    % =======================================
    % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    % =======================================
end

% =======================================
% Close the file.
close(vidObj);
% =======================================