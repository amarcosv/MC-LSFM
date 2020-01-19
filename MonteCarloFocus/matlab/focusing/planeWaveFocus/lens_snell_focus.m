% Simulation of a plane wave hitting a planar interface using snell's law

circd = @(radius,deg_ang,xo,yo)  [radius*cosd(deg_ang)+xo;  radius*sind(deg_ang)+yo];

lambda = 500;
n0 = 1.45;
n1 = 1.0;
theta_ref = 0;
lradius = 2000; % in nm

c = 3E+8; %c = 3E+8 m/s x 10+9 nm/m x 10-9 ns/s = 3E+8 nm/ns;
k0 = 2*pi*n0/lambda;
w = k0*c/n0;
wo = 2*pi/100*c;% the background freq. for time

x = [-lradius:20:lradius];
y = [-lradius*3:20:lradius];
[X,Y] = meshgrid(x,y);

k0_x = k0*sind(theta_ref);
k0_y = k0*cosd(theta_ref);
% The lens arc:
rcenter = [0 lradius];
XYlens = circd(lradius,linspace(180,360,181),rcenter(1),rcenter(2));

theta_inc = 90-180/pi*angle((XYlens(1,:)-rcenter(1)) + 1i*(XYlens(2,:)-rcenter(2)));


% Prepare the new file.
% folder = 'C:\Users\jripoll\OneDrive - Universidad Carlos III de Madrid\UC3M\Powerpoint\';
% vidObj = VideoWriter([folder 'plwave_lens_focus_500nm.avi']);
% vidObj.FrameRate = 30;
% vidObj.Quality = 95;
%
% open(vidObj);

t = 2*pi/wo*(0:0.05:4);
for it=1:1%length(t)
    
    Rott0 = zeros(size(X));
    E = zeros(size(X));
    
    for ii=1:1:length(theta_inc)
        
        Rott = (X-XYlens(1,ii))*sind(theta_inc(ii))+(Y-XYlens(2,ii))*cosd(theta_inc(ii));

        Rott0 = Rott0 + (Rott>0);
%         Rott0=(Rott>0);
        
        xn0 = min(x);
        xn1 = max(x);
        yn0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(-theta_inc(ii));
        yn1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(-theta_inc(ii));
        
        figure(1)
        %subplot(1,2,1)
        hold off
        imagesc(x/1000,-y/1000,Rott0);
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
        title(['Lens of radius ' num2str(lradius/1000) '\mum']);
        axis ij
        hold on
        plot(XYlens(1,:)/1000,-XYlens(2,:)/1000,'w.')
        axis equal tight
        
        
        theta = -(180-theta_inc(ii));
        k1 = 2*pi*n1/lambda; % The wavelength in medium1 is lambda/n1
        theta1 = asind(k0/k1*sind(theta));
        
        if isreal(theta1)
%             k1_x = k1*sind(theta1);
%             k1_y = k1*cosd(theta1);
            
%             arci = circd(500,linspace(90-theta,90,10),XYlens(1,ii),XYlens(2,ii));
%             arcr = circd(500,linspace(90-theta_inc(ii)+theta1,90-theta_inc(ii),10),XYlens(1,ii),XYlens(2,ii));
            
            new_angle = theta_inc(ii) - theta1;
            k1_xn = k1*sind(new_angle);
            k1_yn = k1*cosd(new_angle);
            
            deph1 = k1_xn*XYlens(1,ii) + k1_yn*XYlens(2,ii);
            deph0 = k0_x*XYlens(1,ii) + k0_y*XYlens(2,ii);
            

            
             E = E + (0*exp(1i*(k0_x*X + k0_y*Y ).*(Rott<=0)) + exp(- 1i*(k1_xn*X + k1_yn*Y - deph0 - deph1).*(Rott>0))).*exp(1i*w*t(it));
%             E =(0*exp(1i*(k0_x*X + k0_y*Y ).*(Rott<=0)) + exp(- 1i*(k1_xn*X + k1_yn*Y - deph0 - deph1).*(Rott>0))).*exp(1i*w*t(it));
            
%             xn0 = min(x);
%             xn1 = max(x);
%             yn0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(-theta_inc(ii));
%             yn1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(-theta_inc(ii));
%             yp0 = XYlens(2,ii) + (xn0-XYlens(1,ii))*tand(90-theta_inc(ii));
%             yp1 = XYlens(2,ii) + (xn1-XYlens(1,ii))*tand(90-theta_inc(ii));
%             if (theta_inc(ii)-180) > 0
%                 xr0 = XYlens(1,ii);
%                 xr1 = min(x);
%                 yr0 = XYlens(2,ii);
%                 yr1 = XYlens(2,ii) - (xr1-XYlens(1,ii))*tand(-90-theta1+theta_inc(ii));
%             else
%                 xr0 = XYlens(1,ii);
%                 xr1 = max(x);
%                 yr0 = XYlens(2,ii);
%                 yr1 = XYlens(2,ii) - (xr1-XYlens(1,ii))*tand(-90-theta1+theta_inc(ii));
%             end
%             
         end
        
        figure(2)
        subplot(1,2,1)
        hold off
        imagesc(x/1000,-y/1000,real(E.*(Rott0>0)));
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
        title({'${\it R_e}\{E\}$'},'Interpreter','latex');
        axis ij
        hold on
        plot(XYlens(1,:)/1000,-XYlens(2,:)/1000,'w-','LineWidth',5)
        axis equal tight
        subplot(1,2,2)
        hold off
        imagesc(x/1000,-y/1000,abs(E.*(Rott0>0)).^2);
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
        title({'$\mid E \mid^2$'},'Interpreter','latex');
        axis ij
        hold on
        plot(XYlens(1,:)/1000,-XYlens(2,:)/1000,'w-','LineWidth',5)
        axis equal tight
        
        
        %         figure(2)
        %         pause(0.1);
        %         f = getframe(gcf);
        %         writeVideo(vidObj,f);
        
        
    end
    figure(3)
    subplot(1,2,1)
    hold off
    imagesc(x/1000,-y/1000,real(E.*(Rott0>0)));
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    title({'${\it R_e}\{E\}$'},'Interpreter','latex');
    axis ij
    hold on
    plot(XYlens(1,:)/1000,-XYlens(2,:)/1000,'w-','LineWidth',5)
    axis equal tight
    subplot(1,2,2)
    hold off
    imagesc(x/1000,-y/1000,abs(E.*(Rott0>0)).^2);
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    title({'$\mid E \mid^2$'},'Interpreter','latex');
    axis ij
    hold on
    plot(XYlens(1,:)/1000,-XYlens(2,:)/1000,'w-','LineWidth',5)
    axis equal tight
end
% Close the file.
% close(vidObj);