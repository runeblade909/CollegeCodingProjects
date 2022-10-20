function mountainContourPlot(Mountain_Range,Mountain_Gradient)

syms x y;
[xs,ys] = meshgrid(0:.2:5, 0:.2:5);
Mountain_Contour = fcontour(Mountain_Range,[0,5,0,5]);
Mountain_Contour.LevelStep = 250;
hold on;
quiver(xs,ys,subs(Mountain_Gradient(1,1),{x,y},{xs,ys}),subs(Mountain_Gradient(2,1),{x,y},{xs,ys}));
title('Contour Plot');
hold off;

end