% 6/19/12 Andy Reagan
% Goal is to search the space of free parameters of general rogue wave
% solutions to the NLS. I'll be using analytical solutions from Ohta and Yang (2012)
% and searching the space using a "monte carlo" approach

% This is the 2nd order solution

search=200; %number of points to include

r=10; %radius of search

a=zeros(1,search); %store the values of the parameter choice

res=150; %number of points to plot surface with in each dimension

for q=1:search
    a3 = r*rand*exp(2*pi*i*rand); %random complex within desired radius
    a(q) = a3;
    Z = zeros(res);
    for w = 1:res
        for j = 1:res
            X = 8*w/res-4;
            T = 4*j/res - 2;
            Z(w,j) = 1+(24*((3*X - 6*X^2 + 4*X^3 - 2*X^4 - 48*T^2 + 48*X*T^2 - 48*X^2*T^2 - 160*T^4)+ i*T*(-12 + 12*X - 16*X^3 + 8*X^4 + 32*T^2 - 64*X*T^2 + 64*X^2*T^2 + 128*T^4)+ 6*a3*(1 - 2*X + X^2 - 4*i*T + 4*i*X*T - 4*T^2)+ 6*conj(a3)*(-X^2 + 4*i*X*T + 4*T^2 )))/((9 - 36*X + 72*X^2 - 72*X^3 + 72*X^4 - 48*X^5 + 16*X^6)+96*T^2*(3 + 3*X - 4*X^3 +2*X^4) + 384*T^4*(5 - 2*X+2*X^2)+1024*T^6 + 24*(a3+conj(a3))*(3*X^2 - 2*X^3 - 12*T^2 +24*X*T^2)+48*i*(a3 - conj(a3))*(3*T +6*X*T - 6*X^2*T + 8*T^3)+144*a3*conj(a3));    
        end
    end    
    % plot figure
    mesh(-2:4/res:2-4/res,-4:8/res:4-8/res,abs(Z));
    a_str=['a3 = ' num2str(a3)];
    title(a_str);
    fname = ['r' num2str(r) '_' num2str(q)];
    saveas(1,fname,'fig');
    q % output iteration
end

% plot points searched
figure;
plot(real(a),imag(a),'o');
axis([-2*r,2*r,-2*r,2*r]);
tit=['Radius ' num2str(r) ' sample points'];
title(tit);
fname = ['r' num2str(r) 'sampling'];
saveas(2,fname,'fig');