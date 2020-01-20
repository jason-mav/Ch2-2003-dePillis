close all;
clear i j period active_days;

simTime = tf;
sim('model_depillis',simTime);

I_0 = (Io_out.data(1)*100);

path_nodrug = sprintf('figures\\I(0)=0%d_nodrug', I_0);
path_pulsed_cells = sprintf('figures\\I(0)=0%d_pulsed', I_0);
path_pulsed_drug = sprintf('figures\\I(0)=0%d_pulsed-d', I_0);

%%% Drug-free -------------------
%%% Cells figure
fig_cells = figure(1);
plot(Cells_out.time,Cells_out.data,'LineWidth',1)

title('Drug-free populations', 'fontsize',12)
set(gca,'FontSize',11)
xlabel('Days','fontsize',12)
ylabel('Cells','fontsize',12)
legend('N','T','I')

saveas(fig_cells, path_nodrug, 'fig');
print(fig_cells,'-dpng',strcat(path_nodrug,'.png'));

%%% Traditionally pulsed -------------------
%%% Cells figure
fig_cells2 = figure(2);
plot(pulsed_Cells_out.time,pulsed_Cells_out.data,'LineWidth',1)

title('Traditional pulsed', 'fontsize',12)
set(gca,'FontSize',11)
xlabel('Days','fontsize',12)
ylabel('Cells','fontsize',12)
legend('N','T','I','M')

saveas(fig_cells2, path_pulsed_cells, 'fig');
print(fig_cells2,'-dpng',strcat(path_pulsed_cells,'.png'));

%%% Drug figure

fig_drug = figure(3);
stairs(pulsed_Drug_out.time,pulsed_Drug_out.data,'k','LineWidth',1)

title(sprintf('Drug limit: %g', max_Drug_out.data(1)),'fontsize',12)

set(gca,'FontSize',11)
xlabel('Days','fontsize',12)
ylabel('Drug','fontsize',12)
ylim([0 1.2])

saveas(fig_drug, path_pulsed_drug, 'fig');
print(fig_drug,'-dpng',strcat(path_pulsed_drug,'.png'));

