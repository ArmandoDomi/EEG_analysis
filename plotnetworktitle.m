function h = plotnetworktitle(meaM,mealimV,nodenameM,tittxt,figno)
% h = plotnetworktitle(meaM,mealimV,nodenameM,figno)
% PLOTNETWORK draws a network of weighted or binary directed connections.
% If the adjacency matrix is symmetric the undirected connections are drawn 
% by straight lines, otherwise directed connections are drawn as arcs with
% arrows in the middle of the arcs. The weight of the connection is given 
% by line thickness. If the connections are binary, all arcs (or straight 
% lines if undirected) have the same line thickness.
% INPUTS:
% - meaM      : the adjacency (measure) matrix of size K x K, where K is 
%               the number of nodes.
% - mealimV   : the 2 x 1 vector of measure limits to be drawn. The two 
%               limits determines the smallest and largest thickness of the 
%               arcs or lines. The default is the minimum and maximum value
%               in the meaM.
% - nodenameM : a matrix of K rows (of varying length of strings of 
%               non-blank characters) of node names. Default is an ordered
%               index.
% - tittxt    : the given text to be displayed in the title, default is
%               empty.
% - figno     : the figure handle number. The default is the current
%               figure.

[nnode,tmp] = size(meaM);
if nnode~=tmp, error('The input adjacency matrix is not square.'); end
% Check the input arguments and set default values where not given.
if nargin==4
    figno = gcf;
elseif nargin==3
    figno = gcf;
    tittxt = [];
elseif nargin==2
    figno = gcf;
    tittxt = [];
    nodenameM = '1';
    for inode=2:nnode
        nodenameM = str2mat(nodenameM,int2str(inode));
    end
elseif nargin==1
    figno = gcf;
    tittxt = [];
    nodenameM = '1';
    for inode=2:nnode
        nodenameM = str2mat(nodenameM,int2str(inode));
    end
    mealimV = [];
end
if isempty(nodenameM)
    nodenameM = '1';
    for inode=2:nnode
        nodenameM = str2mat(nodenameM,int2str(inode));
    end
end
meaM = meaM - diag(diag(meaM)); % Force diagonal components to be zero

% Check whether all non-zero entries in the adjacency matrix have the same
% value and in this case assign the network as binary.
tmpV = meaM(:);
tmpV(tmpV==0)=[];
if min(tmpV)==max(tmpV)
    isbinary=1;
else
    isbinary=0;
end

% Internal fixed parameters for the sizes of lines, points and fonts.
if isbinary
    linewidthlimitV = [0 1]; % The limits from smallest to largest line width to be plotted
else
    linewidthlimitV = [0 4]; % The limits from smallest to largest line width to be plotted
end
dotnode = 50/sqrt(nnode); % The size of the dot indicating the node
fontsizenode = 50/sqrt(nnode); % The font size for the node name

% Check if the adjacency matrix is symmetric
if sum(sum(abs(meaM-meaM')))<100*eps
    symyesno = 'y'; %input('Is the input matrix symmetric? (y/n)>','s'); %Εβαλα y για συμμετρικό
    if symyesno=='y'
        issymmetric = 1;
    else
        issymmetric = 0;
        f=inline('2*1.03*R*sin(th/2)-R*th','th','R'); % For arc drawing
    end
else
    issymmetric = 0;
    f=inline('2*1.03*R*sin(th/2)-R*th','th','R'); % For arc drawing
end

meamin = min(min(meaM)); 
meamax = max(max(meaM));
% if no limits for the measure are given set them to min and max measure
% values
if isempty(mealimV) && meamin<meamax
    if meamax<=0
        mealimV = 0;
        maxminmealim = 0;
    elseif meamin>0
        mealimV = [meamin-0.05*(meamax-meamin) meamax];
        maxminmealim = mealimV(2) - mealimV(1);
    else
        mealimV = [0 meamax];
        maxminmealim = mealimV(2) - mealimV(1);
    end
elseif isempty(mealimV) && meamin>=meamax
    maxminmealim = 0;
    mealimV = 0;
else 
    maxminmealim = mealimV(2) - mealimV(1);
end    
% maxminmealim = mealimV(2) - mealimV(1);
maxminlinewidthlimit = linewidthlimitV(2) - linewidthlimitV(1);


% Put the nodes on a circle
thnodeV = linspace(0,2*pi,nnode+1);
locationsM = NaN*ones(nnode,2);
[locationsM(:,1),locationsM(:,2)] = pol2cart(thnodeV(1:nnode)',ones(nnode,1));


h = figure(figno);
clf
% First draw the nodes and the node names
plot(locationsM(:,1),locationsM(:,2),'k.','Markersize',dotnode);
hold on
text(1.05*locationsM(:,1),1.05*locationsM(:,2),deblank(nodenameM),...
    'FontSize',fontsizenode,'HorizontalAlignment','center');
if nansum(nansum(meaM))>0 && issymmetric
    for irow=1:nnode-1
        node1 = locationsM(irow,:);
        for icol=irow+1:nnode
            if meaM(irow,icol)-mealimV(1) > 0
                normmea=(meaM(irow,icol)-mealimV(1))/maxminmealim;
                normmea = min(max(0, normmea),1);
                linewidthnow = linewidthlimitV(1) + maxminlinewidthlimit * normmea;
                node2 = locationsM(icol,:);
                % Draw a straight line between the pair of nodes
                plot([node1(1) node2(1)],[node1(2) node2(2)],'b','linewidth',linewidthnow)
            end % if positive measure
        end % for icol
    end % for irow
elseif nansum(nansum(meaM))>0
    % First draw all the arcs and find the largest arc to use as reference
    % for the arrow size.
    xC = cell(nnode,nnode);
    yC = cell(nnode,nnode);
    dM = zeros(nnode,nnode);
    for irow=1:nnode
        node1 = locationsM(irow,:);
        for icol=1:nnode
            if icol~=irow && meaM(irow,icol)-mealimV(1) > 0
                normmea=(meaM(irow,icol)-mealimV(1))/maxminmealim;
                normmea = min(max(0, normmea),1);
                linewidthnow = linewidthlimitV(1) + maxminlinewidthlimit * normmea;
                node2 = locationsM(icol,:);
                % Drawing of arc
                difnode=node2-node1;
                [difth,difR]=cart2pol(difnode(1),difnode(2));
                options=optimset('Display','on');
                thnew=fzero(f,sqrt(24*(1-difR/(1.03*difR))),options,difR);
                Rnew=difR/(2*sin(thnew/2));
                [x,y]=pol2cart(difth+(pi-thnew)/2,Rnew);
                CC=node1+[x,y];
                newnode=node1-CC;
                thnode1=cart2pol(newnode(1),newnode(2));
                [x,y]=pol2cart(linspace(thnode1,thnode1+thnew),Rnew);
                x = x+CC(1);
                y = y+CC(2);
                plot(x,y,'b','linewidth',linewidthnow)
                xC{irow,icol}=x;
                yC{irow,icol}=y;
                dM(irow,icol)=norm([x(1)-x(100) y(1)-y(100)]);
            end % if positive measure
        end % for icol
    end % for irow
    % [a,b]=max(dM); % Alternative use of max instead of mean for arc size
    % [tmp,jbest]=max(a);
    % ibest=b(jbest);
    tmpV=dM(:);
    tmpV(tmpV==0)=[];
    dmean = mean(tmpV);
    tmpM = abs(dM-dmean*ones(size(dM)));
    [a,b]=min(tmpM);
    [tmp,jbest]=min(a);
    ibest=b(jbest);

    xbest = xC{ibest,jbest};
    ybest = yC{ibest,jbest};
    dbest = norm([xbest(45)-xbest(55) ybest(45)-ybest(55)]);
    % Having the reference length for the arrow, draw the arrow at the
    % middle of each arc
    for irow=1:nnode
        for icol=1:nnode
            if icol~=irow && meaM(irow,icol)-mealimV(1) > 0
                x = xC{irow,icol};
                y = yC{irow,icol};
                dtryV = NaN*ones(49,1);
                for ii=1:49
                    pstart = [x(50-ii) y(50-ii)];
                    pend = [x(50+ii) y(50+ii)];
                    dtryV(ii) = norm([pstart(1)-pend(1) pstart(2)-pend(2)]);
                end
                [tmp,iibest]=min(abs(dtryV-dbest));
                % Draw the arrow
                pstart = [x(50-iibest) y(50-iibest)];
                pend = [x(50+iibest) y(50+iibest)];
                ioffset = round(0.3*2*iibest);
                offset1(1) = x(50-iibest)-(y(50-iibest+ioffset)-y(50-iibest));
                offset1(2) = y(50-iibest)-(x(50-iibest)-x(50-iibest+ioffset));
                offset2(1) = x(50-iibest)+(y(50-iibest+ioffset)-y(50-iibest));
                offset2(2) = y(50-iibest)+(x(50-iibest)-x(50-iibest+ioffset));
                fill([offset1(1) pend(1) offset2(1) pstart(1) offset1(1)],...
                    [offset1(2) pend(2) offset2(2) pstart(2) offset1(2)],'k')
            end % if positive measure
        end % for icol
    end % for irow
end % if symmetric
% daspect([1,1,1])
set(gca,'Visible','Off')
if ~isempty(tittxt)
	text(0,1.10,tittxt,'HorizontalAlignment','Center')
end