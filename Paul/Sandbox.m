% Sandbox for Paul to prove things to himself
clear


%% Match filtering on phase shifts
% xx = (-50:.03:50)';
% n = length(xx);
% lambda_c = pi/6;
% rr1 = sqrt(xx.^2      + (2e3*lambda_c).^2);
% rr2 = sqrt((xx+0).^2 + (2e3*lambda_c).^2);
% 
% % Generate phase delay at this distance along Az to filter with
% az = exp(1i.*(2.*pi*(2*rr1)./lambda_c));
% azPlus = exp(1i.*(2.*pi*(2*rr2)./lambda_c));
% % Shift to center frame
% w = exp(1i * 2 * pi * (n-1)/2 * [0:floor(n/2)-1 floor(-n/2):-1]'/ n);
% 
% azPlus = azPlus.*exp(1i*0);
% 
% mf = ifft( fft(azPlus)  .* conj( fft(az) ) .* w);
% 
% figure(5)
% clf
% subplot(311)
% plot(xx,real(az),'b-');
% hold on
% plot(xx,imag(az),'r-');
% plot(xx,real(azPlus),'b--');
% plot(xx,imag(azPlus),'r--');
% legend('real','imaginary')
% 
% subplot(312)
% plot(xx,real(mf),'b-');
% hold on
% plot(xx,imag(mf),'r-');
% plot(xx,abs(mf),'c-')
% legend('real','imaginary','abs')
% 
% subplot(313)
% plot(xx,angle(mf));

%% Range correct data
dx = .1;
dy = 1;

x = -50:dx:50;
r1 = sqrt(100.^2 + (x+15).^2);
r2 =  sqrt(125.^2 + (x).^2);
y = (1:dy:200);
y_sub = y(1:end-round(20/dy));
z = exp(-(y'-r1).^2) +  exp(-(y'-r2).^2);

tic
[yy,xx] = ndgrid(y,x);
lookup = griddedInterpolant(yy,xx,z);
z_shift = lookup(y'+x.^2./(2*y'),xx);
toc

figure(3)
subplot(211)
prettyPlot(z)
colorbar
subplot(212)
prettyPlot(z_shift)
colorbar

% %% Hyperbolas everywhere
% x = -5:.1:5;
% figure(2)
% clf 
% d = 100;
% for i = 1:20
%     y = sqrt(x.^2 + (d+i*1)^2)-i-d;
%     plot(x,y)
%     hold on
% end
% plot(x,(x.^2)/(2*d),'k-.')
%     

