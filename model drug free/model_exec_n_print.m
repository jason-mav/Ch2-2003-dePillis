close all;
clear i j;

save = 1;

I0 = 0.15;
I0_ts = timeseries(I0);

simTime = tf;
sim('model_depillis_drugfree',simTime);

fig_drugfree = figure();
hold on;
plot(Drugfree_Cells_out.time,Drugfree_Cells_out.data(:,1), 'LineWidth',1)
plot(Drugfree_Cells_out.time,Drugfree_Cells_out.data(:,2), 'LineWidth',1)
plot(Drugfree_Cells_out.time,Drugfree_Cells_out.data(:,3), 'LineWidth',1)
axis([0 tf 0 1])
set(gca,'FontSize',11)
I_0 = int8(I0*100);
title(sprintf('Drug-free Cell Populations (I_0=0.%d)', I_0), 'fontsize',12)
xlabel('Days', 'fontsize',12)
ylabel('Cells (10^{11})', 'fontsize',12)
legend('N', 'T', 'I')

% Save the results
if (save == 1)    
    % Print - All populations
    saveas(fig_drugfree, sprintf('figures\\I_0=0%d_drugfree', I_0),'fig');
    print(fig_drugfree,'-dpng',sprintf('figures\\I_0=0%d_drugfree.png', I_0));    
end

