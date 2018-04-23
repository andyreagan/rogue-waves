% % p2.m: the pseudo-spectral method for solving the NLS equation
% % iu_t+u_{xx}+(2|u|^2 - 1)*u=0.
% 
%   L=80; N=256; %dt=0.02; tmax=4; nmax=round(tmax/dt);
%   dx=L/N; x=[-L/2:dx:L/2-dx]'; %k=[0:N/2-1 -N/2:-1]'*2*pi/L; k2=k.^2;
% %  u=1.2*sech(1.2*(x+20)).*exp(i*x)+0.8*sech(0.8*x);
% %  u = 1 - (4*(1-4*i*t))/(1+4*(x-0.5)^2 +16*t^2);
% 
%   u = 0.*x;
%   udata=0.*x; tdata=0;
%   for t=-1:0.5:1                               % integration begins
%     
%       for i=1:length(x)
%       
%           u = 1 - (4*(1-4*i*t))/(1+4*(x(i)-0.5)^2 +16*t^2);
%     
%       end
%       udata=[udata u]; 
%       tdata=[tdata t];
% %       if mod(nn,round(nmax/25)) == 0
% %        udata=[udata u]; tdata=[tdata t*dt];
% %     
% %       end
%     
%   end                                         % integration ends
%   waterfall(x, tdata, abs(udata));           % solution plotting
% %   colormap([0 0 0]); view(10, 60)
% %   text(-2,  -6, 'x', 'fontsize', 15)
% %   text(50, 5, 't', 'fontsize', 15)
% %   zlabel('|u|', 'fontsize', 15)
% %   axis([-L/2 L/2 0 tmax 0 2]); grid off
% %   set(gca, 'xtick', [-40 -20 0 20 40])
% %   set(gca, 'ytick', [0 10 20])
% %   set(gca, 'ztick', [0 1 2])
% 
% 
% 
% a3 = -1/12;
% [X,T] = meshgrid(-4:.2:4);
% X = X-0.5;
% 
% 
% %R = sqrt(X.^2 + Y.^2) + eps;
% %Z = abs(1 - (4*(1-4*i.*T))./(1+4.*(X-0.5)^2 +16.*T.^2));
% z = 1+(24.*((3.*X - 6.*X.^2 + 4.*X.^3 - 2.*X.^4 - 48.*T.^2 + 48.*X.*T.^2 - 48.*X.^2.*T.^2 - 160.*T.^4)+ i.*T.*(-12 + 12.*X - 19.*X.^3 + 8.*X.^4 + 32.*T.^2 - 64.*X.*T.^2 + 64.*X.^2.*T.^2 + 128.*T.^4)+ 6*a3.*(1 - 2.*X + X^2 - 4.*i.*T + 4.*i.*X.*T - 4.*T^.2)+ 6*conj(a3).*(-X.^2 + 4*i.*X.*T + 4.* T.^2 )))./((9 - 36.*X + 72.*X.^2 - 72.*X.^3 + 72.*X.^4 - 48.*X.^5 + 16.*X.^6)+96.*T.^2.*(3 + 3.*X - 4.*X.^3 +2.*X.^4) + 384.*T.^4.*(5 - 2.*X+2.*X.^2) +1024.* T.^6 + 24*(a3+conj(a3)).*(3.*X.^2 - 2.*X.^3 - 12.*T.^2 +24.*X.*T.^2)+48.i.*(a3 - conj(a3)).*(3.*T +6.*X.*T - 6.*X.^2.*T + 8.*T.^3)+144*a3*conj(a3));
% 
% 
% % 1ST ORDER SOLUTION
% 
% Z = zeros(100);
% for w = 1:100
%     for j = 1:100
%         X = 0; T = 0; z = 0;
%         X = 8*w/100-4;
%         T = 4*j/100 - 2;
%         z = 1 - (4*(1-4*sqrt(-1)*T))/(1+4*(X-0.5).^2 +16*T^2);
%         %z = 1+(24.*((3.*X - 6.*X.^2 + 4.*X.^3 - 2.*X.^4 - 48.*T.^2 + 48.*X.*T.^2 - 48.*X.^2.*T.^2 - 160.*T.^4)+ i.*T.*(-12 + 12.*X - 16.*X.^3 + 8.*X.^4 + 32.*T.^2 - 64.*X.*T.^2 + 64.*X.^2.*T.^2 + 128.*T.^4)+ 6*a3.*(1 - 2.*X + X^2 - 4*i.*T + 4.*i.*X.*T - 4.*T^.2)+ 6*conj(a3).*(-X.^2 + 4*i.*X.*T + 4.* T.^2 )))./((9 - 36.*X + 72.*X.^2 - 72.*X.^3 + 72.*X.^4 - 48.*X.^5 + 16.*X.^6)+96.*T.^2.*(3 + 3.*X - 4.*X.^3 +2.*X.^4) + 384.*T.^4.*(5 - 2.*X+2.*X.^2) +1024.* T.^6 + 24*(a3+conj(a3)).*(3.*X.^2 - 2.*X.^3 - 12.*T.^2 +24.*X.*T.^2)+48*i*(a3 - conj(a3)).*(3.*T +6.*X.*T - 6.*X.^2.*T + 8.*T.^3)+144*a3*conj(a3));    
%         Z(w,j) = z;
%     end
% end
% mesh(-4:8/100:4-8/100,-2:4/100:2-4/100,abs(Z))
% 
% 
% 
% % 2ND ORDER SOLUTION
% 
% a3 = 5/3;
% Z = zeros(100);
% for w = 1:100
%     for j = 1:100
%         X = 0; T = 0; z = 0;
%         X = 8*w/100-4;
%         T = 4*j/100 - 2;
%         z = 1+(24*((3*X - 6*X^2 + 4*X^3 - 2*X^4 - 48*T^2 + 48*X*T^2 - 48*X^2*T^2 - 160*T^4)+ i*T*(-12 + 12*X - 16*X^3 + 8*X^4 + 32*T^2 - 64*X*T^2 + 64*X^2*T^2 + 128*T^4)+ 6*a3*(1 - 2*X + X^2 - 4*i*T + 4*i*X*T - 4*T^2)+ 6*conj(a3)*(-X^2 + 4*i*X*T + 4*T^2 )))/((9 - 36*X + 72*X^2 - 72*X^3 + 72*X^4 - 48*X^5 + 16*X^6)+96*T^2*(3 + 3*X - 4*X^3 +2*X^4) + 384*T^4*(5 - 2*X+2*X^2)+1024*T^6 + 24*(a3+conj(a3))*(3*X^2 - 2*X^3 - 12*T^2 +24*X*T^2)+48*i*(a3 - conj(a3))*(3*T +6*X*T - 6*X^2*T + 8*T^3)+144*a3*conj(a3));    
%         %z = 1+(24.*((3.*X - 6.*X.^2 + 4.*X.^3 - 2.*X.^4 - 48.*T.^2 + 48.*X.*T.^2 - 48.*X.^2.*T.^2 - 160.*T.^4)+ i.*T.*(-12 + 12.*X - 16.*X.^3 + 8.*X.^4 + 32.*T.^2 - 64.*X.*T.^2 + 64.*X.^2.*T.^2 + 128.*T.^4)+ 6*a3.*(1 - 2.*X + X^2 - 4*i.*T + 4.*i.*X.*T - 4.*T^.2)+ 6*conj(a3).*(-X.^2 + 4*i.*X.*T + 4.* T.^2 )))./((9 - 36.*X + 72.*X.^2 - 72.*X.^3 + 72.*X.^4 - 48.*X.^5 + 16.*X.^6)+96.*T.^2.*(3 + 3.*X - 4.*X.^3 +2.*X.^4) + 384.*T.^4.*(5 - 2.*X+2.*X.^2) +1024.* T.^6 + 24*(a3+conj(a3)).*(3.*X.^2 - 2.*X.^3 - 12.*T.^2 +24.*X.*T.^2)+48*i*(a3 - conj(a3)).*(3.*T +6.*X.*T - 6.*X.^2.*T + 8.*T.^3)+144*a3*conj(a3));    
%         Z(w,j) = z;
%     end
% end
% mesh(-2:4/100:2-4/100,-4:8/100:4-8/100,abs(Z))
% 
% % solutions from txt files

% I = i;
% Z = zeros(50);
% z = importdata('order2_peaks_formatted.txt');
% s = z{1};
% for w = 1:50
%     for j = 1:50
%         x = 16*w/100-4;
%         t = 8*j/100 - 2;
%         Z(w,j) = eval(s);
%     end
%     w
% end
% mesh(-2:8/100:2-8/100,-4:16/100:4-16/100,abs(Z))

order=5;
I = i;
Z = zeros(50);
a = importdata('0__.txt');
b = importdata('1__.txt');
for w = 1:50
    for j = 1:50
        x = 16*w/100-4;
        t = 8*j/100 - 2;
        a_mat = zeros(order);
        b_mat = zeros(order);
        for c=1:order
            for d=1:order
                a_mat(c,d) = eval(a{(c-1)*order+d});
            end
        end
        for c=1:order
            for d=1:order
                b_mat(c,d) = eval(b{(c-1)*order+d});
            end
        end
        Z(w,j) = det(b_mat)/det(a_mat);
    end
    w
end
mesh(-2:8/100:2-8/100,-4:16/100:4-16/100,abs(Z))

