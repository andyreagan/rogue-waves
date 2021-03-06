function [] =P3(m,y)

% p3.m: fourth-order split-step method for solving the NLS equation
% i*eps*u_t+eps^2*u_{xx}+2|u|^2u=0.

  eps=.02; %epsilson
    
  L=4*pi; N=m; dt=.01/m; tmax=1; nmax=round(tmax/dt);
  dx=L/N; x=[-L/2:dx:L/2-dx]'; k=[0:N/2-1 -N/2:-1]'*2*pi/L;
  % u=sech(x); 
 
  % u=(1./cosh(x)).*exp(-2*i*log(cosh(x))/eps+i*pi);
  u=exp(-(x.^2)).*exp(-i*log(cosh(x))/eps);
  
  
    udata=u; tdata=0;
    
    
  c=1/(2-2^(1/3));                     % scheme coefficients
  a1=c/2; a2=(1-c)/2; a3=a2; a4=c/2;
  b1=c; b2=1-2*c; b3=c;
  E1=exp(-a1*dt*i*k.^2*eps);
  E2=exp(-a2*dt*i*k.^2*eps);
  E3=exp(-a3*dt*i*k.^2*eps);
  E4=exp(-a4*dt*i*k.^2*eps);
  for nn=1:nmax                        % integration begins
    v=ifft(fft(u).*E1);
    v=v.*exp(b1*dt*i*2*v.*conj(v)/eps);
    v=ifft(fft(v).*E2);
    v=v.*exp(b2*dt*i*2*v.*conj(v)/eps);
    v=ifft(fft(v).*E3);
    v=v.*exp(b3*dt*i*2*v.*conj(v)/eps);
    u=ifft(fft(v).*E4);
    if mod(nn,round(nmax/25)) == 0
       udata=[udata u]; tdata=[tdata nn*dt];
    end
  end                                  % integration ends
  
  y = N/m;
  
  x_small = x(1:y:end);
  u_small = udata(1:y:end,:);

  waterfall(x_small, tdata, abs(u_small)');    % solution plotting 
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
  
  fname = [num2str(N) 'pts_' num2str(m) 'plotted'];  %num2str(eps) 'eps_' 
  saveas(1,fname,'fig');       % save figure
  saveas(1,fname,'png');
