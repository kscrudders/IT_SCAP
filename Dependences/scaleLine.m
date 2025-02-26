function [hl ht] = scaleLine(ha,len,pxlsize,loc,direction,color,width,fontsize,showtext)
% SCALELINE     make scalebar for axes
%
% [hl ht] = scaleLine(ha,len,pxlsize,loc,direction,color,width,fontsize)
% INPUTS
%   ha - axes handle. Default gca.
%   len - desired length of scale bar um. Default 5.
%   pxlsize - size of pixels in um. Default 0.267.
%   loc - desired location [x y]. x specifies the right edge of the
%         scalebar (Default sum(max(get(ha,'xlim'))*[1 -.1])). y specifies
%         the y position of the line sum(max(get(gca,'ylim'))*[1 -.1]).
%   direction - binary for directin of bar Default 1 - horizontal.
%   0=vertical
%   color - string or vector of length 3 indicating the desired color of
%           the scalebar. Default [1 1 1].
%   width - desired line width. Default 3.
%   fontsize - desired fontsize. Default 16.
%   showtext - binary show length in um. Default 1 - show text.
% OUPUTS
%   hl - line handle
%   ht - text handle
%% Default values
if nargin <1 || isempty(ha)
    ha = gca;
end
if nargin <2 || isempty(len)
    len = 5;
end
if nargin <3 || isempty(pxlsize)
    pxlsize = 0.267;
end
xlim = get(ha,'xlim');
ylim = get(ha,'ylim');
ydist = (ylim(2)-ylim(1))/200;
xdist = (xlim(2)-xlim(1))/100;
if nargin <4 || isempty(loc)
    loc = [max(xlim) max(ylim)];
end
if nargin <5 || isempty(direction)
    direction = 'hor';
end
if nargin <6 || isempty(color)
    color = [1 1 1];
end
if nargin <7 || isempty(width)
    width = 3;
end
if nargin <8 || isempty(fontsize)
    fontsize = 16;
end
if nargin <9 || isempty(showtext)
    showtext = 1;
end


%% Body
if direction ==1
    if loc(1) + len/pxlsize > max(xlim)
        loc(1) = max(xlim)-xdist;
    end
else
    if loc(2) + len/pxlsize > min(ylim)
        loc(2) = 5+ydist;
    end
end
if showtext == 1
    ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' µm'],'verticalalignment','top','horizontalalignment','center','fontweight','bold','fontsize',fontsize,'color',color,'parent',ha);
    a = get(ht,'extent');
    if a(1)+a(3) > max(xlim)
        loc(1) = max(xlim)-a(3)/2;
        delete(ht)
        ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' µm'],'verticalalignment','top','horizontalalignment','center','fontweight','bold','fontsize',fontsize,'color',color,'parent',ha);
    end
    if a(2)+a(4) > max(ylim)
        loc(2) =   max(ylim)-a(4)-ydist;
        delete(ht)
        ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' µm'],'verticalalignment','top','horizontalalignment','center','fontweight','bold','fontsize',fontsize,'color',color,'parent',ha);
    end
end
if direction == 1
    hl = line([loc(1)-len/pxlsize loc(1)],[loc(2) loc(2)],'color',color,'linewidth',width,'parent',ha);   
else
    hl = line([max(xlim)-8 max(xlim)-8],[loc(2) loc(2)+len/pxlsize],'color',color,'linewidth',width,'parent',ha);    
end
set(get(ha,'parent'),'inverthardcopy','off')
