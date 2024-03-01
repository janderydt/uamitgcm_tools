function PlotLengthScale(x0,y0,ticks,dy)

%% everything in km!!

dx = ticks(2:end) - ticks(1:end-1);

for ii=1:length(ticks)-1
    patch(gca,x0+ticks(ii)+[0 0 dx(ii) dx(ii)],y0+[0 dy dy 0],[0 0 0]+mod(ii,2),'edgecolor','k');   
end

for ii=1:length(ticks)
        text(x0+ticks(ii),y0+dy+10,num2str(ticks(ii)),'horizontalalignment','center','verticalalignment','bottom',...
            'fontsize',10,'color','k');
    if ii==length(ticks)
        text(x0+ticks(ii)+200,y0+dy+10,'km','horizontalalignment','left','verticalalignment','bottom',...
            'fontsize',10,'color','k');
    end
end


