% p3.m: fourth-order split-step method for solving the NLS equation
% iu_t+u_{xx}+2|u|^2u=0.

  m=16384;
  eps=0.01;
  eps_str='point01eps';

  L=300; N=m; dt=(.1*L/N)^2; tmax=5; nmax=round(tmax/dt);
  dx=L/N; x=[-L/2:dx:L/2-dx]'; k=[0:N/2-1 -N/2:-1]'*2*pi/L;

  
  % Peregrine at t=0;
  t=-5;
  u=(1-(4*(1+2*i.*t))./(1+4*x.^2+4.*t.^2)).*exp(i.*t);
  
  % Adding noise

  u = u + eps*randn(m,1);
  
  RMSdata=0;
  udata=u; tdata=0;
  c=1/(2-2^(1/3));                     % scheme coefficients
  a1=c/2; a2=(1-c)/2; a3=a2; a4=c/2;
  b1=c; b2=1-2*c; b3=c;
  E1=exp(-a1*dt*i*k.^2);
  E2=exp(-a2*dt*i*k.^2);
  E3=exp(-a3*dt*i*k.^2);
  E4=exp(-a4*dt*i*k.^2);
  for nn=1:nmax                        % integration begins
    v=ifft(fft(u).*E1);
    v=v.*exp(b1*dt*i*2*v.*conj(v));
    v=ifft(fft(v).*E2);
    v=v.*exp(b2*dt*i*2*v.*conj(v));
    v=ifft(fft(v).*E3);
    v=v.*exp(b3*dt*i*2*v.*conj(v));
    u=ifft(fft(v).*E4);
    if mod(nn,round(nmax/25)) == 0
       anal=(1-(4*(1+2*i.*nn*dt))./(1+4*x.^2+4.*(nn*dt).^2)).*exp(i.*nn*dt);
       RMSE=sqrt(sum(abs(u(:))-abs(anal(:))).^2/numel(x));
       RMSdata= [RMSdata RMSE];
       udata=[udata u]; tdata=[tdata nn*dt];
    end
  end                                  % integration ends
  
  % display RMSE vector

  display(tdata);
  display(RMSdata);
  plot(tdata,RMSdata);

  fname=[num2str(m) 'pts_error_' eps_str];
  saveas(1,fname,'fig');       % save figure
  saveas(1,fname,'png');

  % plot with 512 grid points
  
  y=N/512;

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
  
  fname = ['peregrine_' num2str(N) 'pts_512plotted' eps_str '_normal'];
  saveas(1,fname,'fig');       % save figure
  saveas(1,fname,'png');
