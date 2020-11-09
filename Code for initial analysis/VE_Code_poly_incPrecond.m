%% Code information
%This code is developed by Dr Rosti Readioff on 10/07/2019 to extract and
%analyse data collected from tensile testing machine. Please cite this code
%appropriately when used & contact r.readioff@keele.ac.uk for how to cite.
%To use this code properly some lines might need to be modified. The Excell
%sheet needs to be set up in the same manner as in the sample spreadsheet.
%All formatting in excel sheet needs to be cleared.
%
%Fitting stress-strain polynomial lines
%Fitting Etan-stress polynomial lines
clear all; close all; clc;
%% Import raw data
numData=xlsread('47L.xlsx'); %load Excel sheet - change name accordingly.
Time=numData(11:end,1); %There is an issue reading the 1st row hence 
%starting data at 11th cell
Extension=numData(11:end,2);
Load=numData(11:end,3);
CycleCount=numData(11:end,4);
Area=numData(5,2);
Length=numData(4,2);
RawData=[Time, Extension, Load, CycleCount]; %You can uncomment this when
%required
%% Balance parameters, calc stress, strain, Etan, Hysteresis
i=1;
for m=1:length(Extension) %Balance parameters, calc stress & strain
Extension_b(i)=Extension(m,end)- Extension(1,1); %balancing extension
Load_b(i)=Load(m,end)-Load(1,1); %balancing load
Strain(i)=Extension_b(end,m)/Length; %calc of strain
Stress(i)=Load_b(end,m)/Area; %calc of stress
Etan=diff(Stress)./diff(Strain); % calc of Etan. This calc causes delay
i=i+1;
end
Etan(isinf(Etan)|isnan(Etan)) = 0; %replacing Inf and NaN with zero

v=2;
for w=1:length(Stress)-1 %Calc energy (area under curve)
Energy(v)=abs(0.5*(Stress(w+1)+Stress(w))*(Strain(w+1)-Strain(w)));
v=v+1;
end

RawData_b=[Strain', Stress',[0;Etan'],Energy',CycleCount,Extension_b'];

x1=1;
for x=0:0.5:19 %Calc area under curve - load & unload
    AA{x1}=RawData_b(:,5)==x;
    AAA{x1}=RawData_b(AA{x1},4); %picking energy load&unload together
    BBB{x1}=RawData_b(AA{x1},6); %picking extension data
    x1=x1+1;
end 
x1=1;
for x=0:0.5:18.5
    if x==0
        AUC_k{x1}=sum(AAA{1,x1});
    Oddcells{x1}=AAA{1,x1+1}(1,1);
    AUC_L{x1}=AUC_k{1,x1}+Oddcells{1,x1};
    
    elseif floor(x)==x %AUC for loading
        AUC_k{x1}=sum(AAA{1,x1}(721:end,:));
    Oddcells{x1}=AAA{1,x1+1}(1,1);
    AUC_L{x1}=AUC_k{1,x1}+Oddcells{1,x1};
Extension_0{x1}=BBB{1,x1}(2,1);
Extension_360{x1}=BBB{1,x1}(720,1);
    else %AUC for unloading
        AUC_k{x1}=sum(AAA{1,x1}(2:end,:));
        Oddcells{x1}=AAA{1,x1+1}(1,1);
        AUC_UL{x1}=AUC_k{1,x1}+Oddcells{1,x1};
    end
    x1=x1+1;
end
L=39; % the last recovery cycle which is 39
AUC_L=transpose(cell2mat(AUC_L));
AUC_UL=transpose(cell2mat(AUC_UL));
Extension_0=cell2mat(Extension_0);
LastCycle_0=BBB{1,L}(2,1);
Extension_0=[Extension_0';LastCycle_0]; 
Extension_360=cell2mat(Extension_360);
LastCycle_360=BBB{1,L}(720,1);
Extension_360=[Extension_360';LastCycle_360]; 
Recover_aft_360=Extension_0-Extension_360;
for i = 1:1:19 %Calc Hysteresis
    Hysteresis(i) = AUC_L(i,1)-AUC_UL(i,1);
end
Hysteresis=transpose(Hysteresis);
%% Separate loading & unloading cycles, including precond
n=0; %starting load cycle no. 0
j=0.5; %starting unload cycle no. is 0.5
for k=1:1:20 %Chose 10 because we have 19 cycles in total
    %% Loading
    ind_load{k}=RawData_b(:,5)==n; %Find indices to elements in specified 
    %column of RawData that satisfy the equality
    LoadingCycle{k}=RawData_b(ind_load{k},:); %Use the logical indices to 
    %index into RawData to return required sub-matrices
    if k>1
    LoadingCycle_woRecovery{k}=LoadingCycle{1,k}(721:end,:);%Deleting the 
    %recovery time from the loading cycle (6minutes)
    else
        LoadingCycle_woRecovery{k}=LoadingCycle{1,k};%Not 
        %deleting any data as there is no recovery time.
    end

    %% Unloading
    ind_unload{k}=RawData_b(:,5)==j;%Find indices to elements in specified 
    %column of RawData that satisfy the equality
    UnloadingCycle{k}=RawData_b(ind_unload{k},:); %Use the logical indices 
    %to index into RawData to return required sub-matrices
    %Calc max stress and strain
    j=j+1;
    n=n+1;
end
%loop to correctly pick cyclic data (inc following cyclic value)
for kk=1:1:19 
    LoadingCycle_woRecovery{kk}=[LoadingCycle_woRecovery{1,kk}...
        ;UnloadingCycle{1,kk}(1,:)];
    UnloadingCycle{kk}=[UnloadingCycle{1,kk}(2:end,:)...
        ;LoadingCycle{1,kk+1}(1,:)];
    
    Stress_Lmax{kk}=max(LoadingCycle_woRecovery{1,kk}(:,2));
    Strain_Lmax{kk}=max(LoadingCycle_woRecovery{1,kk}(:,1));
    Stress_ULmax{kk}=max(UnloadingCycle{1,kk}(:,2));
    Strain_ULmax{kk}=max(UnloadingCycle{1,kk}(:,1));
end

%% Calculate polynomial coefficients of stress-strain & stress-Etan
    for C=1:1:19 %Loop to calc coeff
    %% loading
    coeff_Lss{C}=polyfitZero(LoadingCycle_woRecovery{1,C}(:,1),...
        LoadingCycle_woRecovery{1,C}(:,2),4); %polynomial coeff 
    %for stress-strain zero intercept - loading cycle

    coeff_Lets{C}=polyfit(LoadingCycle_woRecovery{1,C}(:,2),...
        LoadingCycle_woRecovery{1,C}(:,3),4); %polynomial coeff
    %for Etan-stress - loading cycle
    
    %% unloading
    coeff_ULss{C}=polyfit(UnloadingCycle{1,C}(:,1),...
        UnloadingCycle{1,C}(:,2),4'); %polynomial coeff 
    %for stress-strain - unloading cycle    
    
    coeff_ULets{C}=polyfit(UnloadingCycle{1,C}(:,2),...
        UnloadingCycle{1,C}(:,3),4); %polynomial coeff
    %for Etan-stress - uloading cycle
    end
%% Calc fitted Stress & Etan
for r=1:1:19
    Strain_presc{r}=(0.0001:0.0025:Strain_Lmax{1,r});
    Stress_presc{r}=(0:0.025:Stress_Lmax{1,r});
for t=1:length(Strain_presc{1,r}(1,:))
Stress_p{r}=((coeff_Lss{1,r}(1,1))*(Strain_presc{1,r}).^4)+...
    ((coeff_Lss{1,r}(1,2))*(Strain_presc{1,r}).^3)+...
    ((coeff_Lss{1,r}(1,3))*(Strain_presc{1,r}).^2)+...
    ((coeff_Lss{1,r}(1,4))*(Strain_presc{1,r}).^1)+...
    ((coeff_Lss{1,r}(1,5))*(Strain_presc{1,r}).^0);

Etan_p{r}=((coeff_Lets{1,r}(1,1))*(Stress_presc{1,r}).^4)+...
    ((coeff_Lets{1,r}(1,2))*(Stress_presc{1,r}).^3)+...
    ((coeff_Lets{1,r}(1,3))*(Stress_presc{1,r}).^2)+...
    ((coeff_Lets{1,r}(1,4))*(Stress_presc{1,r}).^1)+...
    ((coeff_Lets{1,r}(1,5))*(Stress_presc{1,r}).^0);
end
for t2=0.01:0.02:0.05 %condition to calc stress at 0.01,0.03,0.05 strains &
    %Etan at 0.1, 0.3, 0.5 stresses.
    if t2==0.01 
        Stress_01{r}=((coeff_Lss{1,r}(1,1))*(t2).^4)+...
    ((coeff_Lss{1,r}(1,2))*(t2).^3)+...
    ((coeff_Lss{1,r}(1,3))*(t2).^2)+...
    ((coeff_Lss{1,r}(1,4))*(t2).^1)+...
    ((coeff_Lss{1,r}(1,5))*(t2).^0);

        Etan_1{r}=((coeff_Lets{1,r}(1,1))*(0.1).^4)+...
    ((coeff_Lets{1,r}(1,2))*(0.1).^3)+...
    ((coeff_Lets{1,r}(1,3))*(0.1).^2)+...
    ((coeff_Lets{1,r}(1,4))*(0.1).^1)+...
    ((coeff_Lets{1,r}(1,5))*(0.1).^0);

    elseif t2==0.03
        Stress_03{r}=((coeff_Lss{1,r}(1,1))*(t2).^4)+...
    ((coeff_Lss{1,r}(1,2))*(t2).^3)+...
    ((coeff_Lss{1,r}(1,3))*(t2).^2)+...
    ((coeff_Lss{1,r}(1,4))*(t2).^1)+...
    ((coeff_Lss{1,r}(1,5))*(t2).^0);

        Etan_3{r}=((coeff_Lets{1,r}(1,1))*(0.3).^4)+...
    ((coeff_Lets{1,r}(1,2))*(0.3).^3)+...
    ((coeff_Lets{1,r}(1,3))*(0.3).^2)+...
    ((coeff_Lets{1,r}(1,4))*(0.3).^1)+...
    ((coeff_Lets{1,r}(1,5))*(0.3).^0);

    elseif t2==0.05
        Stress_05{r}=((coeff_Lss{1,r}(1,1))*(t2).^4)+...
    ((coeff_Lss{1,r}(1,2))*(t2).^3)+...
    ((coeff_Lss{1,r}(1,3))*(t2).^2)+...
    ((coeff_Lss{1,r}(1,4))*(t2).^1)+...
    ((coeff_Lss{1,r}(1,5))*(t2).^0);

        Etan_5{r}=((coeff_Lets{1,r}(1,1))*(0.5).^4)+...
    ((coeff_Lets{1,r}(1,2))*(0.5).^3)+...
    ((coeff_Lets{1,r}(1,3))*(0.5).^2)+...
    ((coeff_Lets{1,r}(1,4))*(0.5).^1)+...
    ((coeff_Lets{1,r}(1,5))*(0.5).^0);
    end
end
end
%% Plotting data
%plotting Stress vs Strain - loading
figure(1);
title('Stress vs Strain')
xlabel('Strain(mm/mm)', 'fontsize', 14)
%xlim([0 0.05])
ylabel('Stress (MPa)', 'fontsize', 14)
%ylim([0 0.5])
hold on
for u1=1:1:19 %19 no. of cycles
    Heading{u1}=char(strcat('Cycle',num2str(u1)));
    SS(u1)=plot(Strain_presc{1,u1}(1,:),Stress_p{1,u1}(1,:),'color', rand(1,3),'DisplayName',Heading{u1});
    legend('show')
end
savefig('StressVsStrain.fig')
hold off

%plotting Etan vs Stress - loading
figure(2);
title('Etan vs Stress')
xlabel('Stress (MPa)', 'fontsize', 14)
ylabel('Etan (MPa)', 'fontsize', 14)
hold on
for u1=1:1:19 %19 no. of cycles
    Heading{u1}=char(strcat('Cycle',num2str(u1)));
    Et_S(u1)=plot(Stress_presc{1,u1}(1,:),Etan_p{1,u1}(1,:),'color', rand(1,3),'DisplayName',Heading{u1});
    legend('show')
end
savefig('EtanVsStress.fig')
hold off

%plotting Stress vs Cycle No. - loading
figure(3);
title('Stress vs Cycle No.')
xlabel('Cycle No.', 'fontsize', 14)
ylabel('Stress (MPa)', 'fontsize', 14)
hold on
Stress_01=transpose(cell2mat(Stress_01));
Stress_03=transpose(cell2mat(Stress_03));
Stress_05=transpose(cell2mat(Stress_05));
u3=(1:1:19);
S_Cyc=plot(u3,Stress_01,'-*',u3,Stress_03,'-o',u3,Stress_05,'-s');
legend('Strain 0.01','Strain 0.03','Strain 0.05')
savefig('StressVsCycle.fig')
hold off

%plotting Etan vs Cycle No. - loading
figure(4);
title('Etan vs Cycle No.')
xlabel('Cycle No.', 'fontsize', 14)
ylabel('Etan (MPa)', 'fontsize', 14)
hold on
Etan_1=transpose(cell2mat(Etan_1));
Etan_3=transpose(cell2mat(Etan_3));
Etan_5=transpose(cell2mat(Etan_5));
u3=(1:1:19);
Et_Cyc=plot(u3,Etan_1,'-*',u3,Etan_3,'-o',u3,Etan_5,'-s');
legend('Stress 0.1','Stress 0.3','Stress 0.5')
savefig('EtanVsCycle.fig')
hold off

%plotting Energy vs Cycle No. - loading
figure(5);
title('Energy vs Cycle No.')
xlabel('Cycle No.', 'fontsize', 14)
ylabel('Energy (MPa)', 'fontsize', 14)
hold on
u3=(1:1:19);
Ener_Cyc=plot(u3,AUC_L,'-*',u3,AUC_UL,'-o',u3,Hysteresis,'-s');
legend('Area Under Loading','Area Under Unloading',...
    'Area Between Loading & Unloading')
savefig('EnergyVsCycle.fig')
hold off

%plotting Recovery vs Cycle No. - loading
figure(6);
title('Recovery vs Cycle No.')
xlabel('Cycle No.', 'fontsize', 14)
ylabel('Recovery (mm)', 'fontsize', 14)
hold on
u3=(1:1:19);
Recov_Cyc=plot(u3,Extension_0,'-*',u3,Extension_360,'-o',u3,...
    Recover_aft_360,'-s');
legend('Extension After Unloading','Extension After 360s Recovery',...
    'Total Recovery After 360s')
savefig('RecoveryVsCycle.fig')
hold off
%% Exporting Data
%tabulating some of the non-cell arrays
Cycles=u3'; 
Tab_Recovery=table(Cycles,Extension_0,Extension_360, Recover_aft_360);
Tab_AreaUnderCurve=table(Cycles,AUC_L,AUC_UL,Hysteresis);
Tab_Etan_Cycl=table(Cycles,Etan_1,Etan_3,Etan_5);
Tab_Stress_Cycl=table(Cycles,Stress_01,Stress_03,Stress_05);

%writing cells in an excel spreadsheet
writematrix(Cycles,'Data.xls','Sheet','Strain_presc', 'Range','A1')
writecell(Strain_presc','Data.xls','Sheet','Strain_presc', 'Range','B1');
writematrix(Cycles,'Data.xls','Sheet','Stress_calc', 'Range','A1')
writecell(Stress_p','Data.xls','Sheet','Stress_calc', 'Range','B1');
writematrix(Cycles,'Data.xls','Sheet','Stress_presc', 'Range','A1')
writecell(Stress_presc','Data.xls','Sheet','Stress_presc', 'Range','B1');
writematrix(Cycles,'Data.xls','Sheet','Etan_calc', 'Range','A1')
writecell(Etan_p','Data.xls','Sheet','Etan_calc', 'Range','B1');

%writing tabulated data in an excel spreadsheet
writetable(Tab_Recovery,'Data.xls','Sheet','Recovery');
writetable(Tab_AreaUnderCurve,'Data.xls','Sheet','Hysteresis');
writetable(Tab_Etan_Cycl,'Data.xls','Sheet','EtanVsCycl');
writetable(Tab_Stress_Cycl,'Data.xls','Sheet','StressVsCycl');

%% END %%
%
%