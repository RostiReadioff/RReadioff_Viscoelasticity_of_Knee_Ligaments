clear all
load('Recovery.mat')
plot(1:length(HysteresisCycles_Asc),HysteresisCycles_Asc,'-or',...
    1:length(HysteresisCycles_Desc),HysteresisCycles_Desc,'-ok', 'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',12)
ylabel('Dissipated Energy (MPa)','FontWeight','bold','FontName','Times New Roman','FontSize',14)
xlabel('Cycle Number','FontWeight','bold','FontName','Times New Roman','FontSize',14)
xticks(0:2:20);
grid on
legend('Ascending', 'Descending','FontName','Times New Roman','FontSize',14)