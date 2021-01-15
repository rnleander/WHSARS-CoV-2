%%%%%Write on Jul 6 2020
%%%%To solve the within host coronavirus model

clear all
format long
T1=15; %%T=maximum day
T2=30;


%plot(tdata/24, ydata)

%plot(SampleTData/24, SampleYdata, 'o')

%%%%set incubation period to be six days
incubation=5;

%%%%%%%Initial guesses for some parameters 


%%%%%Initial value for optimization
pp=.05;
A1_Initial=1.96*10^4;
A2_Initial=3.29*10^4;
M_Initial=5.99*10^3;
V_Initial=200*10^-6;
y0=[A1_Initial; pp*A2_Initial; (1-pp)*A2_Initial; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0; M_Initial; V_Initial];



rhoX=0.006;

hold on
%%%%use the estimated paprameters to graph the solutions 
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial)
[t,y] = ode15s(funfunfun, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, 6-2+log10(y(:, 14)),'LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model2(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y2] = ode15s(funfunfun2, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t/24-incubation, 6-2+log10(y2(:, 14)),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model3(t,y, rhoX, pp, A1_Initial, A2_Initial,incubation)
[t,y3] = ode15s(funfunfun2, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t/24-incubation, 6-2+log10(y3(:, 14)),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model4(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y4] = ode15s(funfunfun2, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p4=plot(t/24-incubation, 6-2+log10(y4(:, 14)),'--','LineWidth', 2)

%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model5(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y5] = ode15s(funfunfun2, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p5=plot(t/24-incubation, 6-2+log10(y5(:, 14)),'--','LineWidth', 2)

ylabel('Viral load (log_{10})')
xlabel('days')
legend([p1, p2, p3, p4, p5],{'No interferon treatment','\Delta F=0.1 K_F from day 1', '\Delta F=0.2 K_F from day 1', '\Delta F=0.1 K_F from day 3', '\Delta F=0.2 K_F from day 3'},'Location','southwest')
exportgraphics(gcf,'INFTreatment.eps')

hold off 

figure
hold on 
%%%%use the estimated paprameters to graph the solutions 
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial)
[t,y] = ode15s(funfunfun, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, y(:, 2)+y(:, 3)+y(:, 4)+y(:, 5),'LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model2(t,y, rhoX, pp, A1_Initial, A2_Initial,incubation)
[t,y2] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t/24-incubation, y2(:, 2)+y2(:, 3)+y2(:, 4)+y2(:, 5),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model3(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y3] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t/24-incubation, y3(:, 2)+y3(:, 3)+y3(:, 4)+y3(:, 5),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model4(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y4] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p4=plot(t/24-incubation, y4(:, 2)+y4(:, 3)+y4(:, 4)+y4(:, 5),'--','LineWidth', 2)

%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model5(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y5] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p5=plot(t/24-incubation, y5(:, 2)+y5(:, 3)+y5(:, 4)+y5(:, 5),'--','LineWidth', 2)

ylabel('A_2 cells (10^6 cells)')
xlabel('days')
legend([p1, p2, p3, p4, p5],{'No interferon treatment','\Delta F=0.1 K_F from day 1', '\Delta F=0.2 K_F from day 1', '\Delta F=0.1 K_F from day 3', '\Delta F=0.2 K_F from day 3'},'Location','northeast')
exportgraphics(gcf,'IFNTreatmentA1Cells.eps')

hold off 

figure
hold on 
%%%%use the estimated paprameters to graph the solutions 
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial)
[t,y] = ode15s(funfunfun, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, y(:, 6)+y(:, 7),'LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model2(t,y, rhoX, pp, A1_Initial, A2_Initial,incubation)
[t,y2] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t/24-incubation, y2(:, 6)+y2(:, 7),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model3(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y3] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t/24-incubation, y3(:, 6)+y3(:, 7),'--','LineWidth', 2)


%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model4(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y4] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p4=plot(t/24-incubation, y4(:, 6)+y4(:, 7),'--','LineWidth', 2)

%%%%use the estimated paprameters to graph the solutions 
funfunfun2=@(t, y) model5(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)
[t,y5] = ode15s(funfunfun2, [0 24*(T2+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p5=plot(t/24-incubation, y5(:, 6)+y5(:, 7),'--','LineWidth', 2)

ylabel('I cells (10^6 cells)')
xlabel('days')
legend([p1, p2, p3, p4, p5],{'No interferon treatment','\Delta F=0.1 K_F from day 1', '\Delta F=0.2 K_F from day 1', '\Delta F=0.1 K_F from day 3', '\Delta F=0.2 K_F from day 3'},'Location','northeast')
exportgraphics(gcf,'IFNTreatmentInfectedCells.eps')

hold off 




%plot(t/24-incubation, y5(:, 9),'--','LineWidth', 2)


%%%%%The model 
function dydt = model1(t,y, rhoX, pp, A1_Initial, A2_Initial)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KV=1000;
KX=500;
rhoT=0.12;
rhoV=3.1848;
gamma=7.7282;
sigmaA=0.00035;

AT_Initial=A1_Initial+A2_Initial;

delta=sigmaA*A1_Initial/A2_Initial;

KA=AT_Initial/(1-.01);
r2=(delta+sigmaA)/(1-AT_Initial/KA);

a2p=(1-pp)*gamma*(delta+sigmaA);
a2m=pp*gamma*(delta+sigmaA);

diff_rate=1/(7*24);
frac_prolif=delta*7*24;
KA1=A1_Initial/(1-frac_prolif);


mu=0.005;
sigmaM=0.0005;
sigmaMs=0.02;

%%%Immume cells
rM=3; 
rMs=350;

kM0=10^-4;
kM1=3*kM0;

%%%%Alveolar cells infection parameters
%beta=1/5;
%KV=10^5;
%KV=3*10^3;
sigmaI=1/72; 

%%%%%%%%%%%%%%%%%%%%%%
%%%%Fluid volume
C1=20;
C2=250*10^-3;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
KF=100;
%KF=0.008*100000;



%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines
%rhoT=0.08;
%rhoT=0.12;

r=0.1;

%KT=5*10^2;
KT=3*10^2;
%KX=50;
sigmaX=1;
%sigmaY=8.32; %%%%%%%%%%%Not sure 

rhoF1=0.01;
rhoF2=0.0;

qV=1;



%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

%dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2n+mu*A2sn;

dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*(T/(T+KT))*A1;

dydt(2)=r2*(1-AT/KA)*A2+a2m*A2n-(a+a2p+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*(T/(T+KT))*A2+mu*A2s;

dydt(3)=r2*(1-AT/KA)*A2n+a2p*A2-(a+a2m+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*(T/(T+KT))*A2n+mu*A2sn;

dydt(4)=alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-(a+mu+sigmaA)*A2s-r*(T/(T+KT))*A2s;

dydt(5)=alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-(a+mu+sigmaA)*A2sn-r*(T/(T+KT))*A2sn;

dydt(6)=beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-(sigmaI+sigmaA)*I-r*(T/(T+KT))*I;

dydt(7)=alpha*F*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-(sigmaI+sigmaA)*Is-r*(T/(T+KT))*Is;

dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*(T/(T+KT))*(I+Is)-(kM0*M+kM1*Ms)*D/C1;

dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end

 
function dydt = model2(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KV=1000;
KX=500;
rhoT=0.12;
rhoV=3.1848;
gamma=7.7282;
sigmaA=0.00035;

AT_Initial=A1_Initial+A2_Initial;

delta=sigmaA*A1_Initial/A2_Initial;

KA=AT_Initial/(1-.01);
r2=(delta+sigmaA)/(1-AT_Initial/KA);

a2p=(1-pp)*gamma*(delta+sigmaA);
a2m=pp*gamma*(delta+sigmaA);

diff_rate=1/(7*24);
frac_prolif=delta*7*24;
KA1=A1_Initial/(1-frac_prolif);


mu=0.005;
sigmaM=0.0005;
sigmaMs=0.02;

%%%Immume cells
rM=3; 
rMs=350;

kM0=10^-4;
kM1=3*kM0;

%%%%Alveolar cells infection parameters
%beta=1/5;
%KV=10^5;
%KV=3*10^3;
sigmaI=1/72; 

%%%%%%%%%%%%%%%%%%%%%%
%%%%Fluid volume
C1=20;
C2=250*10^-3;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
KF=100;
%KF=0.008*100000;



%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines
%rhoT=0.08;
%rhoT=0.12;

r=0.1;

%KT=5*10^2;
KT=3*10^2;
%KX=50;
sigmaX=1;
%sigmaY=8.32; %%%%%%%%%%%Not sure 

rhoF1=0.01;
rhoF2=0.0;

qV=1;



%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;

if t>(incubation+1)*24
    G=KF*0.1+F;
else
    G=F;
end

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

%dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2n+mu*A2sn;

dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*(T/(T+KT))*A1;

dydt(2)=r2*(1-AT/KA)*A2+a2m*A2n-(a+a2p+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2+mu*A2s;

dydt(3)=r2*(1-AT/KA)*A2n+a2p*A2-(a+a2m+sigmaA)*A2n-alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2n+mu*A2sn;

dydt(4)=alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2s-r*(T/(T+KT))*A2s;

dydt(5)=alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2sn-r*(T/(T+KT))*A2sn;

dydt(6)=beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*I-r*(T/(T+KT))*I;

dydt(7)=alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*Is-r*(T/(T+KT))*Is;

dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*(T/(T+KT))*(I+Is)-(kM0*M+kM1*Ms)*D/C1;

dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end



% function dydt = model2(t,y, rhoX)
% %%%%%%%%%This function is the model to be solved
% %%%%Cell parameters 
% a=0.0014;
% sigmaA=0.0023;
% mu=0.005;
% sigmaM=0.0005;
% sigmaMs=0.02;
% r2=0.074;
% KA=5.5*10^4;
% %%%Immume cells
% rM=3; 
% rMs=142;
% 
% kM0=10^-4;
% kM1=3*kM0;
% 
% %%%%Alveolar cells infection parameters
% beta=1/6;
% KV=10^3;
% sigmaI=1/72; 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% %%%%Fluid volume
% C1=36;
% C2=500;
% 
% 
% %%%%%%Interferon
% sigmaF=0.35;
% alpha=0.6;
% %KF=0.008*100000;
% 
% 
% 
% %%%%%%Virus
% rhoV=9.6;
% rhoVs=rhoV*10^-3;
% sigmaV=1/3;
% 
% %%%%Toxins and chemokines
% %rhoT=0.08;
% rhoT=0.5;
% sigmaT=0.29;
% r=0.04;
% 
% KT=10;
% KX=0.001;
% sigmaX=8.32;
% 
% rhoF1=0.05;
% rhoF2=0;
% 
% KF=4.9;
% 
% 
% %%%%%%%%%%%%
% A1=y(1); A2=y(2); A2n=y(3);
% A2s=y(4); A2sn=y(5);
% I=y(6); Is=y(7); D=y(8);
% F=y(9); X=y(10);  T=y(11); Ms=y(12); M=y(13); V=y(14); 
% AT=A1+A2+A2n+A2s+A2sn+I+Is;
% 
% if t>4*24
%     G=KF*0.05+F;
% else
%     G=F;
% end
% 
% %dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV)-alpha*G*A2/(KF+G)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*G*A2n/(KF+G)-r*T/(T+KT)*A2n+mu*A2sn;
% 
% dydt(4)=alpha*G*A2/(KF+G)-(mu+sigmaA)*A2s-r*T/(T+KT)*A2s;
% 
% dydt(5)=alpha*G*A2n/(KF+G)-(mu+sigmaA)*A2sn-r*T/(T+KT)*A2sn;
% 
% dydt(6)=beta*V*A2/(V+KV)-alpha*G*I/(KF+G)-(sigmaI+sigmaA)*I-r*T/(T+KT)*I;
% 
% dydt(7)=alpha*G*I/(KF+G)-(sigmaI+sigmaA)*Is-r*T/(T+KT)*Is;
% 
% dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*T/(T+KT)*AT-(kM0*M+kM1*Ms)*D/C1;
% 
% 
%  dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*G;
% 
% dydt(10)=rhoX*(I+Is+A2s+Ms)/C2-sigmaX*X;
% 
% 
% 
% dydt(11)=rhoT*Ms/C2-sigmaT*T;
% 
% dydt(12)=kM0*M*(V+D/C1)-r*T/(T+KT)*Ms-sigmaMs*Ms;
% 
% dydt(13)=rM+rMs*X/(X+KX)-r*T/(T+KT)*M-sigmaM*M-kM0*M*(V+D/C1);
% 
% 
% dydt(14)=(rhoVs*Is+rhoV*I+kM0*M*V+kM1*Ms*V)/C1-sigmaV*V;
% 
% dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
% end



function dydt = model3(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KV=1000;
KX=500;
rhoT=0.12;
rhoV=3.1848;
gamma=7.7282;
sigmaA=0.00035;

AT_Initial=A1_Initial+A2_Initial;

delta=sigmaA*A1_Initial/A2_Initial;

KA=AT_Initial/(1-.01);
r2=(delta+sigmaA)/(1-AT_Initial/KA);

a2p=(1-pp)*gamma*(delta+sigmaA);
a2m=pp*gamma*(delta+sigmaA);

diff_rate=1/(7*24);
frac_prolif=delta*7*24;
KA1=A1_Initial/(1-frac_prolif);


mu=0.005;
sigmaM=0.0005;
sigmaMs=0.02;

%%%Immume cells
rM=3; 
rMs=350;

kM0=10^-4;
kM1=3*kM0;

%%%%Alveolar cells infection parameters
%beta=1/5;
%KV=10^5;
%KV=3*10^3;
sigmaI=1/72; 

%%%%%%%%%%%%%%%%%%%%%%
%%%%Fluid volume
C1=20;
C2=250*10^-3;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
KF=100;
%KF=0.008*100000;



%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines
%rhoT=0.08;
%rhoT=0.12;

r=0.1;

%KT=5*10^2;
KT=3*10^2;
%KX=50;
sigmaX=1;
%sigmaY=8.32; %%%%%%%%%%%Not sure 

rhoF1=0.01;
rhoF2=0.0;

qV=1;



%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;

if t>(incubation+1)*24
    G=F+KF*0.2;
else
    G=F;
end

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

%dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2n+mu*A2sn;

dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*(T/(T+KT))*A1;

dydt(2)=r2*(1-AT/KA)*A2+a2m*A2n-(a+a2p+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2+mu*A2s;

dydt(3)=r2*(1-AT/KA)*A2n+a2p*A2-(a+a2m+sigmaA)*A2n-alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2n+mu*A2sn;

dydt(4)=alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2s-r*(T/(T+KT))*A2s;

dydt(5)=alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2sn-r*(T/(T+KT))*A2sn;

dydt(6)=beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*I-r*(T/(T+KT))*I;

dydt(7)=alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*Is-r*(T/(T+KT))*Is;

dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*(T/(T+KT))*(I+Is)-(kM0*M+kM1*Ms)*D/C1;

dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end




function dydt = model4(t,y, rhoX, pp, A1_Initial, A2_Initial,incubation)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KV=1000;
KX=500;
rhoT=0.12;
rhoV=3.1848;
gamma=7.7282;
sigmaA=0.00035;

AT_Initial=A1_Initial+A2_Initial;

delta=sigmaA*A1_Initial/A2_Initial;

KA=AT_Initial/(1-.01);
r2=(delta+sigmaA)/(1-AT_Initial/KA);

a2p=(1-pp)*gamma*(delta+sigmaA);
a2m=pp*gamma*(delta+sigmaA);

diff_rate=1/(7*24);
frac_prolif=delta*7*24;
KA1=A1_Initial/(1-frac_prolif);


mu=0.005;
sigmaM=0.0005;
sigmaMs=0.02;

%%%Immume cells
rM=3; 
rMs=350;

kM0=10^-4;
kM1=3*kM0;

%%%%Alveolar cells infection parameters
%beta=1/5;
%KV=10^5;
%KV=3*10^3;
sigmaI=1/72; 

%%%%%%%%%%%%%%%%%%%%%%
%%%%Fluid volume
C1=20;
C2=250*10^-3;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
KF=100;
%KF=0.008*100000;



%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines
%rhoT=0.08;
%rhoT=0.12;

r=0.1;

%KT=5*10^2;
KT=3*10^2;
%KX=50;
sigmaX=1;
%sigmaY=8.32; %%%%%%%%%%%Not sure 

rhoF1=0.01;
rhoF2=0.0;

qV=1;



%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;


if t>(incubation+3)*24
    G=F+KF*0.1;
else
    G=F;
end

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

%dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2n+mu*A2sn;

dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*(T/(T+KT))*A1;

dydt(2)=r2*(1-AT/KA)*A2+a2m*A2n-(a+a2p+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2+mu*A2s;

dydt(3)=r2*(1-AT/KA)*A2n+a2p*A2-(a+a2m+sigmaA)*A2n-alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2n+mu*A2sn;

dydt(4)=alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2s-r*(T/(T+KT))*A2s;

dydt(5)=alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2sn-r*(T/(T+KT))*A2sn;

dydt(6)=beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*I-r*(T/(T+KT))*I;

dydt(7)=alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*Is-r*(T/(T+KT))*Is;

dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*(T/(T+KT))*(I+Is)-(kM0*M+kM1*Ms)*D/C1;

dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end

function dydt = model5(t,y, rhoX, pp, A1_Initial, A2_Initial, incubation)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KV=1000;
KX=500;
rhoT=0.12;
rhoV=3.1848;
gamma=7.7282;
sigmaA=0.00035;

AT_Initial=A1_Initial+A2_Initial;

delta=sigmaA*A1_Initial/A2_Initial;

KA=AT_Initial/(1-.01);
r2=(delta+sigmaA)/(1-AT_Initial/KA);

a2p=(1-pp)*gamma*(delta+sigmaA);
a2m=pp*gamma*(delta+sigmaA);

diff_rate=1/(7*24);
frac_prolif=delta*7*24;
KA1=A1_Initial/(1-frac_prolif);


mu=0.005;
sigmaM=0.0005;
sigmaMs=0.02;

%%%Immume cells
rM=3; 
rMs=350;

kM0=10^-4;
kM1=3*kM0;

%%%%Alveolar cells infection parameters
%beta=1/5;
%KV=10^5;
%KV=3*10^3;
sigmaI=1/72; 

%%%%%%%%%%%%%%%%%%%%%%
%%%%Fluid volume
C1=20;
C2=250*10^-3;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
KF=100;
%KF=0.008*100000;



%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines
%rhoT=0.08;
%rhoT=0.12;

r=0.1;

%KT=5*10^2;
KT=3*10^2;
%KX=50;
sigmaX=1;
%sigmaY=8.32; %%%%%%%%%%%Not sure 

rhoF1=0.01;
rhoF2=0.0;

qV=1;



%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;


if t>(incubation+3)*24
    G=F+KF*0.2;
else
    G=F;
end

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

%dydt=zeros(1,12);
% dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*T/(T+KT)*A1;
% 
% dydt(2)=r2*(1-AT/KA)*A2-(a+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*F*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2+mu*A2s;
% 
% dydt(3)=r2*(1-AT/KA)*A2n-(a+sigmaA)*A2n-alpha*F*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+F)-r*T/(T+KT)*A2n+mu*A2sn;

dydt(1)=a*(A2+A2n+A2s+A2sn)-sigmaA*A1-r*(T/(T+KT))*A1;

dydt(2)=r2*(1-AT/KA)*A2+a2m*A2n-(a+a2p+sigmaA)*A2-beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2+mu*A2s;

dydt(3)=r2*(1-AT/KA)*A2n+a2p*A2-(a+a2m+sigmaA)*A2n-alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-r*(T/(T+KT))*A2n+mu*A2sn;

dydt(4)=alpha*G*A2/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2s-r*(T/(T+KT))*A2s;

dydt(5)=alpha*G*A2n/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(a+mu+sigmaA)*A2sn-r*(T/(T+KT))*A2sn;

dydt(6)=beta*V*A2/(V+KV+qV*A2/2/C1)-alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*I-r*(T/(T+KT))*I;

dydt(7)=alpha*G*I/(qF*(A2+A2n+A1+I)*10^-2/C1/6.02+KF+G)-(sigmaI+sigmaA)*Is-r*(T/(T+KT))*Is;

dydt(8)=(sigmaI+sigmaA)*(I+Is)+r*(T/(T+KT))*(I+Is)-(kM0*M+kM1*Ms)*D/C1;

dydt(9)=(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end




