clc; clear;
clf
% =======================================
% Prepare the new file.
file=input('Nama File= ','s');
vidObj = VideoWriter(file);
vidObj.Quality = 55;
open(vidObj);
% =======================================

% Parameter
rho=1;
cp=0.6;
K=0.4;

A=K/rho/cp;

Lx=10;

Nx=31;
Nt=450;

dx=Lx/(Nx-1);
dy=dx;

c=1;
C=0.1;
dt=C*dx/c;

Tn=zeros(1,Nx);
x = linspace(0, Lx, Nx);

Tn(:)=0;
t=0;

xx=0:dx:Lx;             %Grid sumbu x
yy=0:dy:5;               %Grid sumbu y
CC=zeros(length(Tn),length(yy));
name='Zaitun';
for n=1:Nt
    Tc = Tn;

    t = t + dt;
    for i=2:Nx-1
        Tn(i)= Tc(i) + dt * A * ((Tc(i+1) - 2*Tc(i) + Tc(i-1))/(dx)^2);
    end
    Sx = round(7.5*Nx/Lx);
    fire="off";
    if (t<=4)
        fire="on";
        Tn(Sx) = Tn(Sx) + dt*10/rho/cp;
    end

    %     Tn(1) =0; Tn(end) = 0; %Dirichlet
    Tn(1) =Tn(2); Tn(end) = 0; %Neumann
    for i=1:1+5/dy
        CC(:,i)=Tn(:);
    end
    %     surf(x,y,transpose(C(:,:,am)));
    %                     view([0 90]);

    subplot(2,1,1);

    plot(x, Tn,'LineWidth',2);
    set(gca, 'ylim',[0 10]);
    grid on;
    xlabel('Length of Rod');
    ylabel('Temperature');
    title('Graph of Temperature Heat',FontSize=12);

    subplot(2,1,2)
    imagesc(xx,yy,CC');
    colorbar('Ticks',[0,1.75,3.5,5.25,7],...
        'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'});
    caxis([0 7]);
    colormap jet
    fr=sprintf('Fire : %s || time = %.2f s',fire,t);
    title(fr,FontSize=12);
    view([0 90]);

    namee=sprintf('Author Name : %s',name);
    sgtitle(namee,FontSize=20);
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




