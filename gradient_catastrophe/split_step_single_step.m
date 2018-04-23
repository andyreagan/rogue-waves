
% I'm going to write a split step code for the focusing NLS equation in
% hopes of reproducing the modulation instability


  eps=.02;

  L=4*pi; N=256; dt=0.01; tmax=2; nmax=round(tmax/dt);
  dx=L/N; x=[-L/2:dx:L/2-dx]'; k=[0:N/2-1 -N/2:-1]'*2*pi/L;

  u=exp(-(x.^2)).*exp(-i*log(cosh(x))/eps);
  
  udata=u; tdata=0;
  
  for nn=1:nmax
    u = ifft(exp(-dt*k.^2*i*eps).*fft(u.*exp(i*2*u.*conj(u)*dt/eps)));
    if mod(nn,round(nmax/25)) == 0
       udata=[udata u]; tdata=[tdata nn*dt];
    end
  end
  
  waterfall(x, tdata, abs(udata'));    % solution plotting
  colormap([0 0 0]); view(5, 60)
  text(-0.4,  -0.4, 'x', 'fontsize', 15)
  text(7, 1, 't', 'fontsize', 15)
  text(-3.7, -0.3, '-\pi', 'fontsize', 14)
  text(2.9, -0.3, '\pi', 'fontsize', 14)
  zlabel('|u|', 'fontsize', 15)
  axis([-L/2 L/2 0 tmax 0 2]); grid off
  set(gca, 'xtick', [-pi 0 pi], 'xticklabel',{'','0',''})
  set(gca, 'ytick', [0 1 2], 'yticklabel',{'0','','2'})
  set(gca, 'ztick', [0 1 2])