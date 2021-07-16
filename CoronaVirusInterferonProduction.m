%%%%%Write on Jul 6 2020
%%%%To solve the within host coronavirus model

clear all
format long
%T=15;
T=30; %%T=maximum day


incubation=5;

%%%%%Paramter values from optimization
pp=.05;
A1_Initial=1.96*10^4;
A2_Initial=3.29*10^4;
M_Initial=5.99*10^3;
V_Initial=20*10^-6;
y0=[A1_Initial; pp*A2_Initial; (1-pp)*A2_Initial; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0; M_Initial; V_Initial];

rhoX=0.006;

%%%%use the estimated paprameters to graph the solutions 
rhoF21=0;
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial,rhoF21)
[t,y] = ode15s(funfunfun, [0 24*(T+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, 6-2+log10(y(:, 14)),'LineWidth', 2)
hold on 

%%%%use the estimated paprameters to graph the solutions 
rhoF22=.01;
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial,rhoF22)
y0=[A1_Initial; pp*A2_Initial; (1-pp)*A2_Initial; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0; M_Initial; V_Initial];
[t2,y2] = ode15s(funfunfun, [0 24*(T+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t2/24-incubation, 6-2+log10(y2(:, 14)),'--','LineWidth', 2)

%%%%use the estimated paprameters to graph the solutions 
rhoF23=.1;
funfunfun=@(t, y) model1(t,y, rhoX, pp, A1_Initial, A2_Initial,rhoF23)
y0=[A1_Initial; pp*A2_Initial; (1-pp)*A2_Initial; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0; M_Initial; V_Initial];
[t3,y3] = ode15s(funfunfun, [0 24*(T+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t3/24-incubation, 6-2+log10(y3(:, 14)),'--','LineWidth', 2)


ylabel('Viral load (log_{10})','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
L1=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF21);
L2=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF22);
L3=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF23);
legend([p1, p2, p3],{L1,L2,L3},'Location','southeast','FontSize', 16)
exportgraphics(gcf,'rhoF2.eps')

hold off 
 

figure
hold on 

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, y(:, 2)+y(:, 3)+y(:, 4)+y(:, 5),'LineWidth', 2)


%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t2/24-incubation, y2(:, 2)+y2(:, 3)+y2(:, 4)+y2(:, 5),'--','LineWidth', 2)


%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t3/24-incubation, y3(:, 2)+y3(:, 3)+y3(:, 4)+y3(:, 5),'--','LineWidth', 2)

ylabel('A_2 cells (10^6 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
L1=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF21);
L2=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF22);
L3=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF23);
legend([p1, p2, p3],{L1,L2,L3},'Location','southeast','FontSize', 16)
exportgraphics(gcf,'rhoF2_A2Cells.eps')

hold off 
 

figure
hold on 

%%%%We assume the viral load in the lung is two magnitue larger than serum
p1=plot(t/24-incubation, y(:,9),'LineWidth', 2)


%%%%We assume the viral load in the lung is two magnitue larger than serum
p2=plot(t2/24-incubation, y2(:, 9),'--','LineWidth', 2)


%%%%We assume the viral load in the lung is two magnitue larger than serum
p3=plot(t3/24-incubation, y3(:, 9),'--','LineWidth', 2)

ylabel('F (pM)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
L1=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF21);
L2=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF22);
L3=sprintf('rhoF2=%.2f pmol/10^6 cell hr', rhoF23);
legend([p1, p2, p3],{L1,L2,L3},'Location','northeast','FontSize', 16)
exportgraphics(gcf,'rhoF2_F.eps')


%%%%%The model 
function dydt = model1(t,y, rhoX, pp, A1_Initial, A2_Initial,rhoF2)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

beta=1/6;
KX=500;
KV=1000;
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

rhoF1=0.01;

%%%%%%Virus

rhoVs=rhoV*10^-3;
sigmaV=1/3;
qV=1;

%%Toxins
r=0.1;


KT=3*10^2;

%Chemokines
sigmaX=1;








%%%%%%%%%%%%
A1=y(1); A2=y(2); A2n=y(3);
A2s=y(4); A2sn=y(5);
I=y(6); Is=y(7); D=y(8);
F=y(9); X=y(10); T=y(11); Ms=y(12); M=y(13); V=y(14); 
AT=A1+A2+A2n+A2s+A2sn+I+Is;

%compute per capita rate of diff

a=diff_rate*(1-A1/KA1);

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