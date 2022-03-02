%part 1
clearvars; clearvars -GLOBAL
close all
%set(0,'DefaultFigureWindowStyle','docked')  % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)
%%electron modeling

%solve for Vth
mo= 9.1093837015E-31;
mn = 0.26*mo;
%may need to an array to do this
l= 200E-9;  %come back into the other side
h= 100E-9;  %bounce back

%%% ask TA
T= 300;
k= 1.38064852E-23;  %for thermal velocity
Tmn= 0.2E-12;  %mean free path

Vth = sqrt((2*k*T)/(mn))
np = 1000; %number of particles
partWatch = randi(np,7,1);  %an array of indices to select a sample of particles
%mean free path
mfp = Vth*Tmn

X = rand(np,1)*l;  %random number between 0 and1 multiplied by limits of silicon
Y = rand(np,1)*h;

Vx = Vth*(randn(np,1));  %x component
Vy = Vth*(randn(np,1));  %y component

dt = h/Vth/100;  %should be 1/100 of the region size

numit=1000;
Xp = X;  %when it hits the top, X remains the same
%when hitting the x max, come back from left so X-width
Yp = Y;   % when it hits the top, direction will just be -Vy
for i=1:numit  %going from 1 iteraton to 5

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

    %calculating temperarure ising avg velocity

    avgV = mean(sqrt(Vx.^2 + Vy.^2));  %vx and vy are vectors
    semiT = (avgV).^2*mn/(2*k)*(1/sqrt(2));

    figure(1);
    hold on
    axis([0 l 0 h])

   %titlep = sprintf('Temp = %.3d', semiT);

    title(sprintf('Temp = %.3d Degrees Celcius', semiT))
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


