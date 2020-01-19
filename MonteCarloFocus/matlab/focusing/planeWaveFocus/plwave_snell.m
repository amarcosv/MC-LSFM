% Simulation of a plane wave hitting a planar interface using snell's law

circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];

lambda = 500;
n0 = 1.45;
n1 = 1.0;
theta_inc = -42.6028; % in degrees 

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;
x = [-2000:2:2000];
y = [-2000:2:2000];
[X,Y] = meshgrid(x,y);

%Prepare the new file.
% folder = 'C:\Users\jripoll\OneDrive - Universidad Carlos III de Madrid\UC3M\Powerpoint\';
% vidObj = VideoWriter([folder 'plwave_inf.avi']);
% vidObj.FrameRate = 30;
% vidObj.Quality = 95;
% open(vidObj);

for itheta = 1:length(theta_inc)
    theta = theta_inc(itheta);
    
    k0_x = k0*sind(theta);
    k0_y = k0*cosd(theta);

    k1 = 2*pi*n1/lambda; % The wavelength in medium1 is lambda/n1
    theta1 = asind(k0/k1*sind(theta));
    k1_x = k1*sind(theta1);
    k1_y = k1*cosd(theta1);

    arci = circd(500,linspace(90-theta,90,10));  
    arcr = circd(500,linspace(-90-theta1,-90,10));  

    t = 2*pi/w*(0:0.05:2);
    for it=1:1%length(t)
        E = exp(1i*(k0_x*X + k0_y*Y).*(Y>0) + 1i*(k1_x*X + k1_y*Y).*(Y<=0)).*exp(1i*w*t(it));
        figure(2)
        hold off
        imagesc(x/1000,-y/1000,real(E));
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
        title(['\lambda_0 = ' num2str(lambda) 'nm. \theta_{inc}=' num2str(theta) '^o']); 
        axis ij
        hold on
        if n0~=n1
            plot(arci(1,:)/1000,-arci(2,:)/1000,'k-.','LineWidth',1)
            plot(arcr(1,:)/1000,-arcr(2,:)/1000,'k-.','LineWidth',1)
            h_hor = line([min(x) max(x)]/1000,[0 0],'color','w','LineWidth',1);
            h_ver = line([0 0],-[min(y) max(y)]/1000,'color','black','LineWidth',0.1,'LineStyle','-.');
            h_inc = line([0 max(x)*sind(theta)]/1000,-[0 max(y)*cosd(theta)]/1000,'color','k','LineWidth',1,'LineStyle',':');    
            h_ref = line([0 min(x)*sind(theta1)]/1000,-[0 min(y)*cosd(theta1)]/1000,'color','k','LineWidth',1,'LineStyle',':');
        end
        pause(0.1);
%         f = getframe(gcf);
%         writeVideo(vidObj,f);
    end
end
% Close the file.
% close(vidObj);
