clear
% Simulation of a plane wave hitting a planar interface using snell's law

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

% The lens arc:
rcenter = [0 lradius];
deg_ang=linspace(270,90,181); %angular sampling
XYlens = [lradius*sind(deg_ang)+rcenter(1);  lradius*cosd(deg_ang)+rcenter(2)];

%Create incoming angles:
theta_inc=deg_ang; %Contains angles between prop axis and lens points, center is located at lens radial center
%all rays come to the lens with direction parallell to the optical axis.


k0_x = k0*sind(theta_ref);
k0_y = k0*cosd(theta_ref);

%Start Simulation:
Rott0 = zeros(size(X));
E = zeros(size(X));

nAngles=1;
% nAngles=length(theta_inc);
% for ii=1:1:nAngles
     ii=70;
    % Rotate Y coordinates towards theta incoming angle
    Rott = (X-XYlens(1,ii))*sind(theta_inc(ii))+(Y-XYlens(2,ii))*cosd(theta_inc(ii));
    % Rott = (X-XYlens(1,ii))+(Y-XYlens(2,ii));
    % Rott=X+Y;
    %       Rott0 = Rott0 + (Rott>0);
    Rott0 = (Rott>0);
    
    %Plot wave space
    figure(1)
    %subplot(1,2,1)
    hold off
    imagesc(x/1000,y/1000,Rott0);
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    title(['Lens of radius ' num2str(lradius/1000) '\mum']);
    axis ij
    hold on
    plot(XYlens(1,:)/1000,XYlens(2,:)/1000,'w.')
    axis equal tight
    
    %calculate angle between surface normal and propagation axis. Now
    %center is at the lens curve
    theta = -(180-theta_inc(ii));
    k1 = 2*pi*n1/lambda; % The wavelength in medium1 is lambda/n1
    theta1 = asind(k0/k1*sind(theta)); %Out angle respect to the normal
    
    
      if isreal(theta1)
%             k1_x = k1*sind(theta1);
%             k1_y = k1*cosd(theta1);
            
%             arci = circd(500,linspace(90-theta,90,10),XYlens(1,ii),XYlens(2,ii));
%             arcr = circd(500,linspace(90-theta_inc(ii)+theta1,90-theta_inc(ii),10),XYlens(1,ii),XYlens(2,ii));
            
            new_angle = theta_inc(ii) - theta1;
            k1_xn = k1*sind(new_angle);
            k1_yn = k1*cosd(new_angle);
            
            %We calculate the dephasing. Then we just add them since the y
            %components have opposite signs.
            deph1 = k1_xn*XYlens(1,ii) + k1_yn*XYlens(2,ii);
            deph0 = k0_x*XYlens(1,ii) + k0_y*XYlens(2,ii);
            
%             E = E + (0*exp(1i*(k0_x*X + k0_y*Y ).*(Rott<=0)) + exp(- 1i*(k1_xn*X + k1_yn*Y - deph0 - deph1).*(Rott>0))).*exp(1i*w*t(it));
            E = E + (0*exp(1i*(k0_x*X + k0_y*Y ).*(Rott<=0)) + exp(- 1i*(k1_xn*X + k1_yn*Y - deph0 - deph1).*(Rott>0)));%.*exp(1i*w*t(it));

      end
% end