function mountainSurfacePlot(Mountain_Range)

Mountain_Graph = fsurf(Mountain_Range,[0,5,0,5]);
Mountain_Graph.FaceAlpha = 0.5;
colorbar;
xlim([0,5]);
ylim([0,5]);
zlim([7000,12000]);
view(135,30);
title('Surface Plot');

end