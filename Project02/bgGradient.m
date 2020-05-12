function varargout = bgGradient(ax,leftcolor,rightcolor)

if isempty(ax)
    ax = gca;
end

lim = axis(ax);
xdata = [lim(1) lim(2) lim(2) lim(1)];
ydata = [lim(3) lim(3) lim(4) lim(4)];
cdata(1,2,:) = rightcolor;
cdata(1,4,:) = rightcolor;
cdata(1,1,:) = leftcolor;
cdata(1,3,:) = leftcolor;
p = patch(xdata,ydata,'k','Parent',ax);
set(p,'CData',cdata, ...
    'FaceColor','interp', ...
    'EdgeColor','none');
uistack(p,'bottom') % Put gradient underneath everything else
if nargout
    varargout{1} = p;
end
set(gca,'Layer','top')