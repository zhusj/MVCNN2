function fig = plotting_series(varargin)
%fig = plotting(A,L,AddMore,fig_in,estimated_L)
%	Function to do neat plotting of labelled data in Matlab figure
%	A: data to be visualized row-wise data organization
%	L: colum vector holds the set of labels
%	AddMore [optional]: boolean to plot more point to previous created figure
%	fig_in [optional except if AddMore = 1] : the previously created figure
%	estimated_L: compare the estimated Labels and the inserted labels and add red circle around the 
%		wrongly classified points

% Copy rights: Amr Bakry 2012

%% prepair parameters
A = varargin{1};
%L = varargin{2};
L = (1:size(A,1))';
if length(varargin)>1 && ~isempty(varargin{2})
    labels = varargin{2};
else
    labels = {'Start point'};%,'End point'};%char('A'+unique(L)-1);%{'Label 1','Label 2','Label 3'};
end
hold_previous = 1;
if length(varargin)>2
    AddMore = varargin{3};
    fig_in = varargin{4};
    if length(varargin)>4
        hold_previous = varargin{5};
    else
        hold_previous = 1;
    end
    if length(varargin)>5
        estimated_L = varargin{6};
        
        if ~exist('estimated_L','var')
            estimated_L = 1;
        end

        if ~exist('AddMore','var')
            AddMore = 0;
        end
    end
else
    AddMore = 0;
end

%% start plotting
colorOrder = zeros(19,3);
for i = 2:19%length(colorOrder)
    colorOrder(i,:) = mod(colorOrder(i-1,:) + [0.23,0.52,0.89],1.0);
end
n = size(A,1);
colorOrder = [linspace(1,0,n);linspace(0,1,n);linspace(0,0,n)]';

style = ['-','--',':','-.'];
colors = ['r','g','b','c','m','k'];

mc = length(colors);
markers = ['.','+','v','o','*','x','s','d','^','>','<','p','h'];
ml = length(markers);

if AddMore==0
    if exist('fig_in','var')
        fig = figure(fig_in);
    else
        fig = figure;
    end
    
    if ~hold_previous
        clf
    end
    hold on
    for i = 1:n-1%unique(L)'%length(A)
        %indx = strmatch(i, L, 'exact');%(L==i);
        plot3(A(i,1),A(i,2),A(i,3),'color' ,colorOrder(i,:),'marker','*','MarkerSize',15,'LineWidth',3);%'Color',colorOrder(i,:)
        plot3(A(i:i+1,1),A(i:i+1,2),A(i:i+1,3),'color' ,colorOrder(i,:),'MarkerSize',15,'LineWidth',3);
        %Arrow3d([ A(i,1),A(i,2),A(i,3)], [ A(i+1,1),A(i+1,2),A(i+1,3)],'color' ,colorOrder(i,:) );
    end
    plot3(A(end,1),A(end,2),A(end,3),'color' ,colorOrder(i,:),'marker','*','MarkerSize',15,'LineWidth',3);%,%,'MarkerSize',4);
    %l=char('A'+unique(L)-1);%{'Label 1','Label 2','Label 3'};
    %legend(labels);
    legend({'Start point'});
    hold off
else
    fig = figure(fig_in);
    hold on
    if length(L)>1
        for i = unique(L)'%length(A)
            indx = strmatch(i, L, 'exact');%(L==i);
            plot3(A(indx,1),A(indx,2),A(indx,3),[colors(mod(i,mc)+1),markers(mod(i,ml)+1)]);
        end
        indx = (estimated_L==L);%well classified samples
        plot3(A(indx,1),A(indx,2),A(indx,3),'co','MarkerSize',10);
		indx = (estimated_L~=L);%wrongly classified samples
        plot3(A(indx,1),A(indx,2),A(indx,3),'ro','MarkerSize',10);
    else
        plot3(A(1),A(2),A(3),[colors(mod(L,mc)+1),markers(mod(L,ml)+1)]);
		if estimated_L==L
			plot3(A(1),A(2),A(3),'co','MarkerSize',10);
		else
			plot3(A(1),A(2),A(3),[colors(mod(estimated_L,mc)+1),markers(mod(estimated_L,ml)+1)],'MarkerSize',10);
        end
    end
    %l=char('A'+unique(L)-1);%{'Label 1','Label 2','Label 3'};
    %legend(l);
    hold off
end

end