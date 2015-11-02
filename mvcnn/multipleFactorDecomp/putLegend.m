hold on
hLegend = legend('a','b','c','d','e','f');          %# Create the legend
hKids = get(hLegend,'Children');    %# Get the legend children
hText = hKids(strcmp(get(hKids,'Type'),'text'));  %# Select the legend children
                                                  %#    of type 'text'
set(hText,{'Color'},{colorOrder(1,:); colorOrder(2,:); colorOrder(3,:); colorOrder(4,:); colorOrder(5,:);colorOrder(6,:)});    %# Set the colors