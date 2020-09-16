function [h,dmax] = plotmts(xM,drawtype,standardize,maxcols,taus,nameM,figno)
% [h,dmax] = plotmts(xM,drawtype,standardize,maxcols,taus,nameM,figno)
% PLOTMTS make a figure of the columns of the input matrix 'xM', which is
% handy for making time history plors of multivariate time series. Two
% options can be given in the input: a) to standardize the time data in
% each column, so that the range is the same for all columns, b) to put a
% maximum number of columns (time series) at each figure, so that if the
% number of columns is larger more figures will be displayed. In that case
% a vector of figure handles is given to the output.
% INPUTS 
% - xM          : the data matrix of size n x m, where n denotes the time
%                 series length and m the number of time series.
% - drawtype    : the type of data drawing: 1->lines, 2->dots, 3->line and
%                 dots. Default is 1.
% - standardize : if set to 1, transform each time series (each column) to
%                 [0,1], if it is 0 leave them at the range they are given,
%                 and if it is a two component vector, use them as the range.
%                 Default is 0.
% - maxcols     : the maximum number of columns to be displayed in one
%                 figure. If m > maxcols more figures will be created.
%                 Default is 10.
% - taus        : the sampling time. If not know just set to 1. Default is
%                 1.
% - nameM       : a cell array or string matrix, its components will be
%                 used as labels for each displayed time series. 
% - figno       : handle of figure to use for the first plot (and if other
%                 plots increase figure handle number by one).
% OUTPUTS       
% - h           : the vector of handles for each of the displayed figures.
% - dmax        : the range of y-axis at each panel.

if nargin==6
    figno = 1;
elseif nargin==5
    figno = 1;
    nameM= [];
elseif nargin==4
    figno = 1;
    nameM= [];
    taus = 1;
elseif nargin==3
    figno = 1;
    nameM= [];
    taus = 1;
    maxcols = 10;
elseif nargin==2
    figno = 1;
    nameM= [];
    taus = 1;
    maxcols = 10;
    standardize = 0;
elseif nargin==1
    figno = 1;
    nameM= [];
    taus = 1;
    maxcols = 10;
    standardize = 0;
    drawtype = 1;
end
if isempty(taus), taus = 1; end
if isempty(maxcols), maxcols = 10; end
if isempty(standardize), standardize = 0; end
if isempty(drawtype), drawtype = 1; end

[nrow,ncol]=size(xM);
tV = [1:nrow]'*taus;

if ~isempty(nameM) && ~((size(nameM,1)==1 && size(nameM,2)==ncol) || (size(nameM,1)==ncol && size(nameM,2)==1))
    error('The list of labels for the time series does not match the number of time series.');
end
% if ~isempty(nameM) && size(nameM,1) ~= ncol
%     error('The list of labels for the time series does not match the number of time series.');
% end

xminV = NaN*ones(ncol,1);
xmaxV = NaN*ones(ncol,1);
for j=1:ncol
    xminV(j)=min(xM(:,j));
    xmaxV(j)=max(xM(:,j));
end
if standardize==1
    for j=1:ncol
        d = xmaxV(j) - xminV(j);
        xM(:,j)= (xM(:,j) - xminV(j))/d;
    end
    xmin = 0;
    xmax = 1;
elseif standardize==0
    xmin = min(xminV);
    xmax = max(xmaxV);
elseif length(standardize)==2
    xmin = standardize(1);
    xmax = standardize(2);
end
dmax = xmax - xmin;

if dmax==0
    dmax =1 ;
end
nfig = ceil(ncol/maxcols);
h = NaN*ones(nfig,1);
for ifig=1:nfig
    ncolend = min(ncol,ifig*maxcols);
    ncolstart = (ifig-1)*maxcols+1;
    ncolnow = ncolend-ncolstart+1;
    
    h(ifig) = figure(figno-1+ifig);
    clf;
    % paperfigure;
    set(gcf,'Position',[10 40 1000 600])
    fz = get(gca,'fontsize');
    hold on
    for j=1:ncolnow
        if drawtype==1
            plot(tV,xM(:,ncolend-j+1)+(j-1)*dmax-xmin,'-b','linewidth',1)
        elseif drawtype==2
            plot(tV,xM(:,ncolend-j+1)+(j-1)*dmax-xmin,'.b','Markersize',8)
        else
            plot(tV,xM(:,ncolend-j+1)+(j-1)*dmax-xmin,'.-b','Markersize',8,'linewidth',1)
        end
        plot([tV(1) tV(nrow)],j*dmax*[1 1],'k')
        if ~isempty(nameM)
            if iscell(nameM)
                nowname = nameM{ncolend-j+1};
            else
                nowname = deblank(nameM(ncolend-j+1,:));
            end
            text(tV(nrow),j*dmax,nowname,'HorizontalAlignment','right',...
                'VerticalAlignment','top','fontsize',fz-2)
        end
    end
    axis([tV(1) tV(nrow) 0 dmax*ncolnow])
    set(gca, 'YTickLabel', [])
    set(gca, 'YTick', [])
    set(gca, 'Box','On')
    xlabel(sprintf('time step t [sampling time=%f]',taus))
    if standardize
        ylabel('x(t) [standardized]')
    else
        ylabel('x(t)')
    end    
end
