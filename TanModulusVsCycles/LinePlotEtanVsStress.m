clear all
load('TangentModulusVsCycle.mat')
%plotting stress_strain (using Sample 52 as an example)
plot(Stress_Tan_Asc(:,1),Stress_Tan_Asc(:,2),'-r',...
    Stress_Tan_Asc(:,3),Stress_Tan_Asc(:,4),'-.r',...
    Stress_Tan_Asc(:,5),Stress_Tan_Asc(:,6),'--r',...
    Stress_Tan_Desc(:,1),Stress_Tan_Desc(:,2),'-k',...
    Stress_Tan_Desc(:,3),Stress_Tan_Desc(:,4),'-.k',...
    Stress_Tan_Desc(:,5),Stress_Tan_Desc(:,6),'--k','LineWidth',1.5)

set(gca,'FontName','Times New Roman','FontSize',12)
ylabel('Tangent Modulus (MPa)','FontWeight','bold','FontName','Times New Roman','FontSize',14)
xlabel('Stress (MPa)','FontWeight','bold','FontName','Times New Roman','FontSize',14)
%xticks(0:2:20);
grid on
legend('0.1 %/min - Asc', '1 %/min - Asc','10 %/min - Asc',...
    '0.1 %/min - Desc', '1 %/min - Desc','10 %/min - Desc',...
    'FontName','Times New Roman','FontSize',14)