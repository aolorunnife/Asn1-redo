%part 1
clearvars; clearvars -GLOBAL
close all
%set(0,'DefaultFigureWindowStyle','docked')  % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)

%%scattering


%solve for Vth
mo= 9.1093837015E-31;
mn = 0.26*mo;
%may need to an array to do this
l= 200E-9;  %come back into the other side
h= 100E-9;  %bounce back
a= l*h;
%%% ask TA
T= 300;
k= 1.38064852E-23;  %for thermal velocity
Tmn= 0.2E-12;  %mean free path

Vth = sqrt((2*k*T)/(mn));
np = 1000; %number of particles
partWatch = randi(np,7,1);  %an array of indices to select a sample of particles
%mean free path
mfp = Vth*Tmn;

X = rand(np,1)*l;  %random number between 0 and1 multiplied by limits of silicon
Y = rand(np,1)*h;

X1 = X;
Y1 = Y;

Vx = Vth*(randn(np,1)-0.5);  %x component
Vy = Vth*(randn(np,1)-0.5);  %y component

dt = h/Vth/50;  %should be 1/100 of the region size


numit=1000;
Xp = X;  %when it hits the top, X remains the same
%when hitting the x max, come back from left so X-width
Yp = Y;   % when it hits the top, direction will just be -Vy

pathnum = 0;
distancesum = 0;
for i=1:numit  

    Xp = X;
    Yp = Y;

    X= X + dt*Vx;
    Y= Y + dt*Vy;


    ix = X < 0;
    X(ix) = X(ix)+l;
    Xp(ix) = Xp(ix) + l;
    Y = Y;

    ix = X > l;
    X(ix) = X(ix)-l;
    Xp(ix) = Xp(ix)-l;
    Y = Y;

    iy = Y<0 | Y > h;
    Vy(iy) = -Vy(iy);

    %scattering, ONLY PLOTS A FEW
    std=Vth/sqrt(2);
    Pscat = 1 - exp(-dt/Tmn);
    iscat = Pscat > rand(np,1);

    X = X + Vx *dt;
    Y = Y + Vy*dt;
    Vx(iscat) = std*randn(sum(iscat),1);
    Vy(iscat) = std*randn(sum(iscat),1);


    %mean free path
    pathnum = pathnum + sum(iscat);
    distance = sqrt((X1(iscat)-X(iscat)).^2 + ((Y1(iscat)-Y(iscat)).^2));
    distancesum = distancesum + sum(distance);
    X1(iscat) = X(iscat);

   
    %calculating temperarure using avg velocity

    avgV = mean(sqrt(Vx.^2 + Vy.^2));  %vx and vy are vectors
    semiT = (avgV).^2*mn/(2*k);

    figure(1);
    hold on
    axis([0 l 0 h])

    %titlep = sprintf('Temp = %.3d', semiT);

    title(sprintf('Temp = %.3d', semiT))
    for j=1:length(partWatch)

        plot([Xp(partWatch(j)),X(partWatch(j))]',[Yp(partWatch(j)),Y(partWatch(j))]','SeriesIndex',j)
        xlabel('x');
        ylabel('y');
    end
    figure (2);
    hold on
    plot(i, semiT, 'bo')
    title('Temperature Plot')
    xlabel('iteration')
    ylabel('temperature')
    pause(0.05)

end

hold off


figPlot = figure (3);
histvelmb = histogram(sqrt(Vx.^2 + Vy.^2),50);
title(' Maxwell Boltzman Distribution')
xlabel('Velocity in both directions (m/s)')
ylabel('# of Particles')

%summing for mean free path 
mfp = sum(distance());
avgmfp = distancesum/pathnum;
Tmn = avgmfp/Vth;

