clc
clear

fprintf(' Hello there friend! Help me do some physics real quick!\n')

v = input(' Please give me a velocity to start with.\n');

theta = input(' What angle are we shooting this thing off at?\n');

throwBallFunc(v,theta);

[X,~,~,~]=throwBallFunc(v,theta);

[~,Y,~,~]=throwBallFunc(v,theta);

[~,~,hitDistance,~]=throwBallFunc(v,theta);

[~,~,~,hitTime]=throwBallFunc(v,theta);

[~,~,~,~,xt,~]=throwBallFunc(v,theta);

[~,~,~,~,~,yt]=throwBallFunc(v,theta);

fprintf('The ball hits the ground at %f seconds.\n',hitTime);

fprintf('The ball hits the ground %f meters from the launch point.\n',X);


figure
hold on
title('Projectile Trajectory')
plot(xt,yt,'r')
xlabel('Distance(m)')
ylabel('Ball height(m)')
%plot the ground

xG=linspace(0,X,100);
ground=zeros(100);
plot(xG,ground,'c--');
xlim([0 3]);
ylim([-1 2]);
name='Ground';
text(xt(1)/2, -0.25,name);
