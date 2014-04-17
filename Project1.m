%% Initialize
clear all
close all
%See wikipedia for soliton solutions 0.5*csech(sqrt(c)/2)(x-ct)-a)
%% Variable Declaration
tic
N = 512;
dt = 0.1/N^2;
% x=(2*pi/N)*(-N/2:N/2-1);%;
% k = [0:N/2-1 0 -N/2+1:-1];
x = linspace(-10,10,N);
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);
k = [0:delta_k:N/2*delta_k,-(N/2-1)*delta_k:delta_k:-delta_k];

a = 9; b = 16; c=2;
u = 3*a^2*sech(.5*(a*(x))).^2;%+3*b^2*sech(.5*(b*(x+2))).^2 +3*c*(sech(sqrt(c)/2*(x+8))).^2;%+

tmax = 0.6; nmax = round(tmax/dt);
U = fft(u);
udata=u'; tdata = 0;
nplt = floor((tmax/25)/dt);
ik3 = 1i*k.^3;
for n = 1:nmax
    t = n*dt; g = -.5i*dt*k;
    E = exp(dt*ik3/2); E2 = E.^2;
    a = g.*fft(real( ifft( U ) ).^2);
    b = g.*fft(real( ifft(E.*(U+a/2)) ).^2); % 4th-order
    c = g.*fft(real( ifft(E.*U + b/2) ).^2); % Runge-Kutta
    d = g.*fft(real( ifft(E2.*U+E.*c) ).^2);
    U = E2.*U + (E2.*a + 2*E.*(b+c) + d)/6;
    if mod(n,nplt) == 0
        u = real(ifft(U)); waitbar(n/nmax)
        udata = [udata u']; tdata = [tdata t];
    end
end
toc
% figure(6)
% set(gcf,'renderer','zbuffer')
% waterfall(x,tdata,udata')
%%
figure(1)

filename = 'test.gif';
for tt = 1:size(tdata,2)
    
    plot(x,udata(:,tt))
    axis([-10 10 0 800 ])
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if tt == 1;
        imwrite(A,map,filename,'gif', 'Loopcount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
