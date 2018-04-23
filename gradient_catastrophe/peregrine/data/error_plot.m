% plot the error in the simulartions of the Peregrine solition under uncertain initial conditions

% pull in the 100 for 0 eps

big = zeros(100,50);

for i=1:100
    filename = ['16384_0_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps0(t) = sum(big(:,t))/100;
end




% tdata is known for each of them

tdata = linspace(0,10,50);

plot(tdata,eps0,'r');
axis([0 10 0 1]);






% pull in the 20 for 1 eps

big = zeros(20,50);

for i=1:20
    filename = ['16384_1_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps1(t) = sum(big(:,t))/20;
end

figure;
plot(tdata,eps1,'b');
axis([0 10 0 1]);





% pull in the 20 for 10 eps

big = zeros(20,50);

for i=1:20
    filename = ['16384_10_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps10(t) = sum(big(:,t))/20;
end

figure;
plot(tdata,eps10,'g');
axis([0 10 0 1]);






% pull in the 20 for 100 eps

big = zeros(20,50);

for i=1:20
    filename = ['16384_100_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps100(t) = sum(big(:,t))/20;
end

figure;
plot(tdata,eps100,'o');
axis([0 10 0 1]);




% pull in the 20 for 1000 eps

big = zeros(20,50);

for i=1:20
    filename = ['16384_1000_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps1000(t) = sum(big(:,t))/20;
end

figure;
plot(tdata,eps1000,'y');

axis([0 10 0 1]);


% pull in the 20 for 10000 eps

big = zeros(20,50);

for i=1:20
    filename = ['16384_10000_' num2str(i)];
    load(filename);
    display(RMSdata);
    for j=1:50
        big(i,j) = RMSdata(j);
    end
end

for t=1:50
    eps10000(t) = sum(big(:,t))/20;
end

figure;
plot(tdata,eps10000,'y');
axis([0 10 0 1]);





figure;

%plot(tdata,eps0,'r',tdata,eps1,'b',tdata,eps10,'g',tdata,eps100,'p',tdata,eps1000,'i',tdata,eps10000,'y');
plot(tdata,eps0,tdata,eps1,tdata,eps10,tdata,eps100,tdata,eps1000,tdata,eps10000);
%plot(tdata,eps0,'r',tdata,eps1,'b',tdata,eps10,'g',tdata,eps100,'p',tdata,eps1000,'i',tdata,eps10000,'y');
%plot(tdata,eps0,'r',tdata,eps1,'b',tdata,eps10,'g',tdata,eps100,'p',tdata,eps1000,'i',tdata,eps10000,'y');

%axis([0 10 0 1]);
