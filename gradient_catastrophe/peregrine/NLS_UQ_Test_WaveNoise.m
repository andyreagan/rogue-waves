% p3.m: fourth-order split-step method for solving the NLS equation
% iu_t+u_{xx}+2|u|^2u=0.

  clear all;

  m=1500;
  eps=.0001;

  L=300; N=m; dt=(.1*L/N)^2; tmax=10; nmax=round(tmax/dt);
  dx=L/N; x=[-L/2:dx:L/2-dx]'; k=[0:N/2-1 -N/2:-1]'*2*pi/L;

  
  % Peregrine at t=0;
  t=-5;
  u=(1-(4*(1+2*i.*t))./(1+4*x.^2+4.*t.^2)).*exp(i.*t);
  
  % Adding noise
  RandWaves=zeros(length(x),1);
  NoisePower=1000; % How many waves to add up
  StorePsuedoRandomOutput=zeros(1,NoisePower);
  NoiseNicenessFactor=1; % 0 for very mean, 2 for very nice. really the wavelength minimum divisor
  for j=1:NoisePower % number of random waves to add
    psuedorandomoutput=rand(1)+NoiseNicenessFactor;
    StorePsuedoRandomOutput(j)=psuedorandomoutput;
    RandWaves=RandWaves+abs(randn(1))*exp(i*(x+rand(1)*300-150)/((abs(psuedorandomoutput))));
  end
  
  SignChanges=find(real(RandWaves(1:end-1)).*real(RandWaves(2:end))<0);
  yelp=floor(length(SignChanges)/2);
  Henry=zeros(yelp,1);
  for j=1:yelp-1
     Henry(j)=max(real(RandWaves(SignChanges(2*j):SignChanges(2*j+2))));
  end
  John=mean(Henry); % John is the mean wave height
  
  Ursula=eps/John; % Ursula will make John of RandWaves equal to eps
  figure;
  plot(x,Ursula*real(RandWaves));
  %figure;
  %plot(x,abs(u));
  u=u+Ursula*RandWaves;
  figure;
  plot(x,abs(u));
  
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
    if mod(nn,round(nmax/50)) == 0
       anal=(1-(4*(1+2*i.*nn*dt))./(1+4*x.^2+4.*(nn*dt).^2)).*exp(i.*nn*dt);
       RMSE=sqrt(sum((abs(u(:))-abs(anal(:))).^2)/numel(x));
       RMSdata= [RMSdata RMSE];
       udata=[udata u]; tdata=[tdata nn*dt];
    end
  end                                  % integration ends
  
  % save RMSE vector
  y=N/512;

  x_small = x(1:y:end);
  u_small = udata(1:y:end,:);
  figure;
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
