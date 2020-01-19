w0=0.1/2;

w=0.1:0.1:1.5;
lambda=500e-6;
f=pi*w0.*sqrt(w.^2-w0^2)./lambda;

figure, plot(w,f);
xlabel('Input beam waist w [mm]');
ylabel('Lens focal distance')
zR=pi*w0^2/lambda;

f=5;
w=w0*sqrt(1+(f/(pi*w0))^2);