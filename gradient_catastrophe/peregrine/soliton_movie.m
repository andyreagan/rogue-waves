% 6/21/12 Eric Roma, Andy Reagan
% This code generates a movie of Peregrine Soliton solution in the 1+1D
% NLS equation

clear all;

a3 = 5*i/2;
a5 = 1/240;

res=100;

Z = zeros(res);

for w = 1:res
    for j = 1:res
        x = 8*w/res-4;
        t = 16*j/res-8;
        Z(j,w) = (1-(4*(1+2*i*x))/(1+4*t^2+4*x^2))*exp(i*x);
    end
end

for b=1:res
    c = Z(:,b);
    plot(-8:16/100:8-8/100,abs(c));
    axis([-8,8,-0.1,5]);
    F(b)=getframe;
end

movie(F,3,15);