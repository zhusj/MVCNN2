function fig = plotting(varargin)
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
L = varargin{2};
L = L(:);
if length(varargin)>2 && ~isempty(varargin{3})
    labels = varargin{3};
else
    labels = char('A'+unique(L)-1);%{'Label 1','Label 2','Label 3'};
end

if length(varargin)>3
    AddMore = varargin{4};
    fig_in = varargin{5};
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
    hold on
    uL = unique(L);
    for i = 1:length(uL)%length(A)
        v = uL(i);
        indx = strmatch(v, L, 'exact');%(L==i);
        plot3(A(indx,1),A(indx,2),A(indx,3),[colors(mod(i,mc)+1),markers(mod(i,ml)+1)],'MarkerSize',10,'LineWidth',2);%,%,'MarkerSize',4);'Color',colorOrder(i,:)
    end
    %l=char('A'+unique(L)-1);%{'Label 1','Label 2','Label 3'};
    legend(labels);
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