function result = growth_fn_2(mmol_m2_s)

if mmol_m2_s>800
    mmol_m2_s=800;
end

data=[0, 0
36, 0 %https://doi.org/10.1371/journal.pone.0278159 and DOI: 10.32615/ps.2020.013
100, 127.98
200, 145.65 
400, 158.45 
600, 162.89
800, 135.56];


pot_per_m2=20; %1 pot is 0.2*0.25 m^2 so 1m^2 holds 20 pots.
growth_perhour_per_second = (14*17);

data2=zeros(801,2);
data2(:,1)=0:800;
data2(37:end,2)=interp1(data(:,1),data(:,2),data2(37:end,1),'makima');


PPFD = data2(:,1);
M = data2(:,2)*pot_per_m2/growth_perhour_per_second;
% 
% set(0, 'DefaultLineLineWidth', 2); %set thickness of all the lines = 2
% figure('Renderer', 'painters', 'Position', [500 300 500 420]) % starting point and height - width of the frame
% set(gca, 'ColorOrder', [0 0 0;0 0 0; 1 0 0;0 0.5 0;0.75, 0.75, 0], 'NextPlot', 'replacechildren');% color of lines in the plot with the given order. remember it is periodic
% hAx=gca;
% yyaxis left
% plot(data2(:,1),data2(:,2))
% ylh=ylabel({'Fresh mass of a single';'above-ground lettuce [g]'});
% ylh.VerticalAlignment	= 'bottom'; %if it is not alligned well, try 'top' and 'bottom' too
% ylim([0 200] )
% yyaxis right
% plot(data(:,1),data(:,2)*pot_per_m2/growth_perhour_per_second,'o')
% hAx.XColor = [0 0 0];
% hAx.YColor = [0 0 0];
% hAx.LineWidth = 1.5;
% axis square
% ylh=ylabel({'Lettuce mass increase rate';'per area, M [g m^{-2} h^{-1}]'});
% xlabel({'Photosynthetic photon flux';'density, PPFD [\mumol m^{-2} s^{-1}]'}')
% ylh.VerticalAlignment	= 'top'; %if it is not alligned well, try 'top' and 'bottom' too
% xlim([0 800])
% ylim([0 200*pot_per_m2/growth_perhour_per_second] )
% set(gca,'FontSize',13)
% set(gca,'XMinorTick','on','YMinorTick','on')
% box on
% saveas(gcf,'growth.png')
% % saveas(gcf,'fig.emf')
% % saveas(gcf,'fig.fig')

result=interp1(PPFD,M,mmol_m2_s);% g lettuce