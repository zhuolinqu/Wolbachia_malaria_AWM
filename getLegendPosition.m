pos_legend = ll.Position; % based on outposition
legend off
pos_canvas = get(gca,'Position');
dx = pos_legend(3)/pos_canvas(3);
dy = pos_legend(4)/pos_canvas(4);
x0 = (pos_legend(1)-pos_canvas(1))/pos_canvas(3);
x1 = x0+dx;
y0 = (pos_legend(2)-pos_canvas(2))/pos_canvas(4);
y1 = y0+dy;