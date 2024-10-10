%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Example code 1: Comparison with matlab default dashed line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%   Data
x = [0, 600];
y = x;
y0 = 0 * x;

h = figure_setups;
axis([0, 600, 0, 15000])
hold on
p_0 = plot(x, y0 + 2, ':', 'DisplayName', 'Line 1 :');
p_1 = plot(y0 + 2, x, '-', 'DisplayName', 'Line 2 -');
p_2 = plot(x, x - 2, '-.', 'DisplayName', 'Line 3 -.');
p_3 = plot(x, x + 0, '--', 'DisplayName', 'Line 4 --');

xlabel('time, days')
ylabel('population')
ll = legend('Location','se');
getLegendPosition;
plotLegendBoundary(h, x0, x1, y0, y1, 'FaceColor', 'w');
fac = 0.02;
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*0, p_0, [1,0],  1);
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*1, p_1, [1,0],  1);
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*2, p_2, [1.0, 0.2, 0.5, 0.2],  9); % '-.'
plotLegendLine(h, x0+dx*fac, y1-0.04 - 0.08*3, p_3, [0.5, 0.2, 0.5,0.2], 8); %'--'
