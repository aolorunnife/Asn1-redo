clearvars; clearvars -GLOBAL
close all
%set(0,'DefaultFigureWindowStyle','docked') % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)

%solve for Vth
mo= 9.1093837015E-31;
mn = 0.26*mo;
%may need to an array to do this
l= 200E-9; %come back into the other side
h= 100E-9; %bounce back
a= l*h;
%%% ask TA
T= 300;
k= 1.38064852E-23; %for thermal velocity
Tmn= 0.2E-12; %mean free path

Vth = sqrt((2*k*T)/(mn));
np = 1000; %number of particles
partWatch = randi(np,20,1); %an array of indices to select a sample of particles
%mean free path
mfp = Vth*Tmn;

X = zeros(np,1); %injecting from left = 0;
Y = rand(np,1)*h;

Vx = Vth*(abs(randn(np,1))) ; %x component
Vy = Vth*(randn(np,1)-0.5); %y component

dt = h/Vth/100; %should be 1/100 of the region size

numit=1000;
Xp = X; %when it hits the top, X remains the same
%when hitting the x max, come back from left so X-width
Yp = Y; % when it hits the top, direction will just be -Vy

figure(1);
hold on
axis([0 l 0 h])


box1 = rectangle('Position',[0.8E-7, 0, 0.4E-7, 0.4E-7]);
box2 = rectangle('Position',[0.8E-7, 0.6E-7 ,0.4E-7,0.4E-7]);

InBox = X > 0.8E-7 & X < 1.2E-7 & (Y > 0.6E-7 | Y < 0.4E-7);
while sum(InBox)> 0

    X(InBox) = rand(sum(InBox),1)*l;
    Y(InBox) = rand(sum(InBox),1)*h;
    InBox = X > 0.8E-7 & X < 1.2E-7 & (Y> 0.6E-7 | Y < 0.4E-7);
end

for i=1:numit %going from 1 iteraton to 5


    Xp = X;
    Yp = Y;

    X= X + dt*Vx;
    Y= Y + dt*Vy;

   
    %injection
    bcx = 1;

    if bcx == 1 % if conditions are on, run old code
        ix = X < 0;
        X(ix) = X(ix)+ l;
        Xp(ix) = Xp(ix) + l;
        Y = Y;

        ix2 = X > l;
        X(ix2) = X(ix2)-l;
        Xp(ix2) = Xp(ix2)-l;
        Y = Y;

    elseif bcx == 0 %if conditions are new, do new code, we dont want to loop left to right
        ix = X < 0;
        X(ix) = X(ix);
        Xp(ix)= Xp(ix);
        Y = Y;

        ix2 = X > l;
        X(ix2) = X(ix2);
        Xp(ix2) = Xp(ix2);
        Y = Y;
    end
    iy = Y<0 | Y > h;
    Vy(iy) = -Vy(iy);

    %scattering, ONLY PLOTS A FEW

    sct = 1; %off

    if sct == 1
        Pscat = 1 - exp((-dt/Tmn));
    elseif sct == 0
        Pscat = 0;
    end

    std=Vth/sqrt(2);
    iscat = Pscat > rand(np,1);
    X = X + Vx *dt;
    Y = Y + Vy*dt;
    Vx(iscat) = std*randn(sum(iscat),1);
    Vy(iscat) = std*randn(sum(iscat),1);

    %calculating temperarure ising avg velocity

    avgV = mean(sqrt(Vx.^2 + Vy.^2)); %vx and vy are vectors
    semiT = (avgV).^2*mn/(2*k);

    InBox1 = X > 0.8E-7 & X < 1.2E-7 & Y> 0.6E-7;
    outsidebox = Xp < 0.8E-7 | Xp > 1.2E-7;
    X(InBox1 & outsidebox) = Xp(InBox1 & outsidebox);
    Vx(InBox1 & outsidebox) = -Vx(InBox1 & outsidebox);
    Vy(InBox1 & ~outsidebox) = -Vy(InBox1 & ~outsidebox);
    Y(InBox1 & ~outsidebox) = Yp(InBox1 & ~outsidebox);

    InBox2 = X > 0.8E-7 & X < 1.2E-7 & Y < 0.4E-7;
    outsidebox = Xp < 0.8E-7 | Xp > 1.2E-7;
    X(InBox2 & outsidebox)= Xp(InBox2 & outsidebox);
    Vx(InBox2 & outsidebox) = -Vx(InBox2 & outsidebox);
    Vy(InBox2 & ~outsidebox) = -Vy(InBox2 & ~outsidebox);
    Y(InBox2 & ~outsidebox) = Yp(InBox2 & ~outsidebox);

    figure(1)
    hold on
    axis([0 l 0 h])

  
    %titlep = sprintf('Temp = %.3d', semiT);

    title(sprintf('Temp = %.3d Particles when Boundaries & Scattering are On', semiT))
    xlabel('x')
    ylabel('y')

    for j=1:length(partWatch)

        plot([Xp(partWatch(j)),X(partWatch(j))]',[Yp(partWatch(j)),Y(partWatch(j))]','SeriesIndex',j)
        
    end
    

end

%electron density

hold off
figPlot = figure (2);
histdens = hist3([X Y], [20 20]);
title('Density Plot')
xlabel(('m/s in x direction'))
zlabel('# of Particles')
ylabel('m/s in y direction')

%temperature map

binsx = discretize(X, 10); %20 bins;
binsy = discretize(Y, 10);

tempmap = zeros(10,10);
for x = 1:10

    for y= 1:10

        avgV = mean(sqrt(Vx(binsx==x& binsy==y).^2 + Vy(binsx==x & binsy==y).^2)); %vx and vy are vectors
        tempmap(x,y) = (avgV).^2*mn/(2*k);


    end

end

surf(tempmap)
title('Temperature of System (Boundaries and Scattering On) in Degrees Celsius')
xlabel('x')
ylabel('y')
zlabel('z')


avgV = mean(sqrt(Vx.^2 + Vy.^2)); %vx and vy are vectors
semiT = (avgV).^2*mn/(2*k);