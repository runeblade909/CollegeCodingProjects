%% Subplots
titleSize = 20;
axisSize = 15;
sz2 = 100;

chosenWingload = 25.37;
[~,WLidx] = min(abs(chosenWingload - Wload));

subplot(2,3,1)
hold on
grid minor
yyaxis left
plot(Wload,E,'LineWidth',2)
yline(7,'g--','LineWidth',2)
yline(10,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Glide Endurance (s)','fontsize',axisSize)

yyaxis right
plot(Wload,CL_e,'LineWidth',2)
ylabel('Coefficient of Lift','fontsize',axisSize)

yyaxis left
scatter(chosenWingload,E(WLidx),sz2,'red')
t2 = text(chosenWingload,E(WLidx),sprintf('      Endurance = %f seconds' ,E(WLidx)));
t2.HorizontalAlignment = 'left';


legend('Emax Achieved','Min Emax','Max Emin','Chosen Wing Load','CL Req','Location','Northeast')
title('Max Endurance & CL vs Wing Load','fontsize',titleSize)
%xline(chosenWingload,'k:','linewidth',1.5)
xlim([10 45])



hold off

subplot(2,3,2)
hold on
grid minor
yyaxis left
plot(Wload,Grange,'LineWidth',2)
yline(70,'g--','LineWidth',2)
yline(100,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Glide Range (m)','fontsize',axisSize)

scatter(chosenWingload,Grange(WLidx),sz2,'red')
t3 = text(chosenWingload,Grange(WLidx),sprintf('     Range = %f meters' ,Grange(WLidx)));
t3.HorizontalAlignment = 'left';

yyaxis right
plot(Wload,CL_r,'LineWidth',2)
ylabel('Coefficient of Lift','fontsize',axisSize)



legend('Rmax Achieved','Min Rmax','Max Rmin','Chosen Wing Load','CL Req','Location','Northeast')
title('Max Glide Range & CL vs Wing Load','fontsize',titleSize)
xlim([10 45])
hold off

subplot(2,3,3)
hold on
grid minor
plot(Wload,Wto,'LineWidth',2)
scatter(chosenWingload,Wto(WLidx),sz2,'red')
t3 = text(chosenWingload,Wto(WLidx),sprintf('     Weight = %f Newtons' ,Wto(WLidx)));
t3.HorizontalAlignment = 'left';
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Aircraft Takeoff Weight (N)','fontsize',axisSize)
title('Varitation of Weight with Wing Load','fontsize',titleSize)
xlim([10 45])
hold off

subplot(2,3,4)
hold on
grid minor
plot(Wload,Vr,'LineWidth',2)
plot(Wload,Ve,'LineWidth',2)
plot(Wload,Vstall,'r--','LineWidth',2)
yline(7,'g--','LineWidth',2)
yline(12,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Velocity (m/s)','fontsize',axisSize)
scatter(chosenWingload,Vr(WLidx),sz2,'red')
t3 = text(chosenWingload,Vr(WLidx),sprintf('     Velocity = %f m/s' ,Vr(WLidx)));
t3.HorizontalAlignment = 'left';


legend('Max Range Velocity','Max Endurance Velocity','Stall Velocity',...
    'Minimum Velocity','Maximum Velocity','Chosen Wing Load','Location','Southeast')

title('Max E and Max R Velocity vs Wing Load','fontsize',titleSize)
xlim([10 45])
hold off

subplot(2,3,5)
hold on
grid minor
yyaxis left
plot(Wload,b,'LineWidth',2)
yline(1,'g--','LineWidth',2)
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Wing Span (m)','fontsize',axisSize)

scatter(chosenWingload,b(WLidx),sz2,'red')
t3 = text(chosenWingload,b(WLidx),sprintf('     Wingspan = %f meters' ,b(WLidx)));
t3.HorizontalAlignment = 'left';

yyaxis right
plot(Wload,Sref,'LineWidth',2)
ylabel('Wing Planar Area (m^2)','fontsize',axisSize)


legend('Span','Max Span','Chosen Wing Load','Planar Area','Location','East')
title('Variation of Wing Geometery','fontsize',titleSize)
xlim([10 45])
hold off

subplot(2,3,6)
hold on
grid minor
plot(Wload,Cost,'LineWidth',2)
scatter(chosenWingload,Cost(WLidx),sz2,'red')
t3 = text(chosenWingload,Cost(WLidx),sprintf('    Cost = %f $' ,Cost(WLidx)));
t3.HorizontalAlignment = 'left';
xlabel('Wing Load (N/m^2)','fontsize',axisSize)
ylabel('Aircraft Cost ($)','fontsize',axisSize)

title('Variation of Cost with Wing Load','fontsize',titleSize)
xlim([10 45])
hold off