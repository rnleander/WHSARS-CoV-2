%%%%%Write on Dec 30 2020
%%%%To solve the within host coronavirus model

clear all
format long
T=10; %%T=maximum day
%to see what would happen in the long term without an adaptive immune
%response
T1=120;

%%%%read viral load data
viraldata=csvread('viraldata.csv', 0, 1);
SampleYdata=viraldata(1:(T+1), 1); %%%%%use the first 10 days data
%%%use interpolation to get the data for each hour
SampleTData=(0:1:T)'*24;
tdata=(0:1:T*24)';
%ydata=interp1(SampleTData, SampleYdata, tdata,'spline');
ydata=interp1(SampleTData, SampleYdata, tdata);
plot(tdata/24, ydata,'LineWidth', 2)
hold on
plot(SampleTData/24, SampleYdata, 'o')

%%%%set incubation period to be two days
incubation0=2;

%%%%%%%Initial guesses 
beta0=1/6;     %%%%%%%%
KV0=10^3;      %%%%%%%%
KX0=500;
%%%%%%%%%%
rhoX0=.006;  %%%%%%%%%%
rhoT0=0.04;   %%%%%%%%%
rhoV0=10;
gamma0=10;


%%%%%Initial value for optimization
%pp=fraction of A2 cells that are suceptible
pp=.05;
A1_Initial=1.96*10^4;
A2_Initial=3.29*10^4;
M_Initial=5.99*10^3;
V_Initial=200*10^-6;
y0=[A1_Initial; pp*A2_Initial; (1-pp)*A2_Initial; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0; M_Initial; V_Initial];

x0=[beta0, KV0, KX0, rhoX0, rhoT0, rhoV0, incubation0, gamma0];   %%%%%%%%%%
%LeastSquareError(x0)
%%%set upper and lower bound for the estimate
lb=[1/6, 10^3, 500, 0.006, 0.12, 1, 1, 1];     %%%%%%%%%
ub=[1/6, 10^3, 500, 0.006, 0.12, 100, 7, 100];    %%%%%%%%%%%%%

%%%%estimate the parameters 
LeastSquareError1=@(x) LeastSquareError(x, y0, T, pp, A1_Initial, A2_Initial);
x = fmincon(LeastSquareError1,x0, [], [], [], [],lb,ub)

beta=x(1);      %%%%%%%%%%%%
KV=x(2);      %%%%%%%%%
KX=x(3);    %%%%%%%%
rhoX=x(4); %%%%%%%%
rhoT=x(5);  %%%%%%%%%
rhoV=x(6);
incubation=x(7);
gamma=x(8);


fprintf('Estimated \\rho_X is %9.4f\n', rhoX)

%%%%use the estimated paprameters to graph the solutions 
funfunfun=@(t, y) model1(t,y, beta, KV, KX, rhoX, rhoT, rhoV, gamma, pp, A1_Initial, A2_Initial)    %%%%%%%%%%
[t,y] = ode15s(funfunfun, [0 24*(T+incubation)], y0);

[t1,y1] = ode15s(funfunfun, [0 24*(T1+incubation)], y0);

%%%%We assume the viral load in the lung is two magnitue larger than serum
plot(t/24-incubation, 6-2+log10(y(:, 14)),'LineWidth', 2)
ylabel('Viral load (log_{10})','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'viral.eps')



figure
plot(t1/24-incubation, y1(:, 1), 'LineWidth', 2)
ylabel('A_1 cells (10^6 cells)','FontSize', 20)
set(gca,'fontsize',16)
xlabel('Days','FontSize', 20)
exportgraphics(gcf,'A1.eps')


figure
plot(t1/24-incubation, y1(:, 2),'LineWidth', 2)
ylabel('A_2^+ cells (10^6 cells)','FontSize', 20, 'LineWidth', 2)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'A2+.eps')

figure
plot(t1/24-incubation, y1(:, 3),'LineWidth', 2)
ylabel('A_2^- cells (10^6 cells)','FontSize', 20, 'LineWidth', 2)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'A2-.eps')


figure
plot(t1/24-incubation, y1(:, 4),'LineWidth', 2)
ylabel('A_2^{*+} cells (10^6 cells)','FontSize', 20, 'LineWidth', 2)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'A2+s.eps')

figure
plot(t1/24-incubation, y1(:, 5),'LineWidth', 2)
ylabel('A_2^{*-} cells (10^6 cells)','FontSize', 20, 'LineWidth', 2)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'A2-s.eps')

figure
plot(t1/24-incubation, y1(:, 2)+y1(:,3)+y1(:,4)+y1(:,5), 'LineWidth', 2)
ylabel('Healthy A_2 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf, 'A2_total.eps')



figure
plot(t1/24-incubation, y1(:, 6),'LineWidth', 2)
ylabel('I cells','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'I.eps')

figure
plot(t1/24-incubation, y1(:, 7),'LineWidth', 2)
ylabel('I^* cells (10^6 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'Is.eps')

figure
plot(t1/24-incubation, y1(:, 8),'LineWidth', 2)
ylabel('D cells (10^6 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'D.eps')

figure
plot(t1/24-incubation, y1(:, 9),'LineWidth', 2)
ylabel('Interferons (pM)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'F.eps')

figure
plot(t1/24-incubation, y1(:, 10),'LineWidth', 2)
ylabel('Chemoxines (pM)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'X.eps')


figure
plot(t1/24-incubation, y1(:, 11),'LineWidth', 2)
ylabel('Toxins (10^6 NETs','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'T.eps')

figure
plot(t1/24-incubation, y1(:, 12),'LineWidth', 2)
ylabel('M^* cells (10^6 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'Ms.eps')

figure
plot(t1/24-incubation, y1(:, 13),'LineWidth', 2)
ylabel('M cells (10^6 cells)','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'M.eps')

figure
plot(t1/24-incubation, log10(y1(:, 14))+6-2,'LineWidth', 2)
ylabel('Log10 Virus','FontSize', 20)
xlabel('Days','FontSize', 20)
set(gca,'fontsize',16)
exportgraphics(gcf,'V_long_run.eps')

%compute R0
sigmaA=0.00035;

sigmaV=1/3;

mu=0.005;

kM0=10^-4;

sigmaI=1/72;

C1=20;

R0=(2*C1*beta*pp*A2_Initial*rhoV/((2*C1*KV+pp*A2_Initial)*(kM0*M_Initial+C1*sigmaV)*(sigmaA+sigmaI)))^.5;

R0_new_cells=2*C1^2*beta*pp*A2_Initial/((2*C1*KV+pp*A2_Initial)*(kM0*M_Initial+C1*sigmaV));

R0_new_virus=rhoV/(C1*(sigmaA+sigmaI));


function error=LeastSquareError(input, y0, T, pp, A1_Initial, A2_Initial)  
%%%for a particular parameter set 'input', this function computes the error
%%%between the simulation and real viral load data

%%%%%The following are parameters to be estimated
beta=input(1);      %%%%%%%%%%%%%%%
KV=input(2);        %%%%%%%%%%%%
KX=input(3);         %%%%%%%%%%
rhoX=input(4);      %%%%%%%%%%%
rhoT=input(5);    %%%%%%%%%%%%%%
rhoV=input(6);
incubation=input(7);
gamma=input(8);



%%read the data 
viraldata=csvread('viraldata.csv', 0, 1);
SampleYdata=viraldata(1:(T+1), 1); %%%%%use the first 10 days data

%%%use interpolation to get viral data for each hour
SampleTData=(0:1:T)'*24;
tdata=(0:1:T*24)';
ydata=interp1(SampleTData, SampleYdata, tdata);
%ydata=interp1(SampleTData, SampleYdata, tdata,'spline');

%%%%Solve the model and compute the error
funfunfun=@(t, y) model1(t,y, beta, KV, KX, rhoX, rhoT, rhoV, gamma, pp, A1_Initial, A2_Initial);    %%%%%%%%
[t,y] = ode45(funfunfun, [0 24*(T+incubation)], y0);
%yydata=interp1(t-incubation*24, 6+log10(y(:, 11))-2, tdata,'spline');
yydata=interp1(t-incubation*24, 6-2+log10(y(:, 14)), tdata);

%%% The following compute the error
error=sum((yydata-ydata).^2);
end



%%%%%The model 
function dydt = model1(t,y, beta, KV, KX, rhoX, rhoT, rhoV, gamma, pp, A1_Initial, A2_Initial)    %%%%%%%%%%%%%%
%%%%%%%%%This function is the model to be solved
%%%%Cell parameters 

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
C2=250;


%%%%%%Interferon
sigmaF=0.35;
alpha=0.6;
qF=20;
%KF=10 before correction to volume
KF=100;




%%%%%%Virus
%rhoV=10;
rhoVs=rhoV*10^-3;
sigmaV=1/3;


%%%%Toxins and chemokines


r=0.1;


KT=3*10^2;

sigmaX=1;


rhoF1=0.01;
%rhoF2=0.03;
%In case infected cells do not produce interferons
rhoF2=0.00;


qV=1;



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

dydt(9)=10^3*(rhoF1*Ms+rhoF2*I+rhoF2*Is)/C2-sigmaF*F;

dydt(10)=10^3*rhoX*(I+Is+A2s+Ms+D)/C2-sigmaX*X;

dydt(11)=rhoT*Ms-kM0*T*M-kM1*T*Ms;

dydt(12)=kM0*M*(V+D/C1)-rhoT*Ms-sigmaMs*Ms-kM1*T*Ms;

dydt(13)=rM+rMs*X/(X+KX)-sigmaM*M-kM0*M*(V+D/C1)-kM0*T*M;

dydt(14)=(rhoVs*Is+rhoV*I-(kM0*M*V-kM1*Ms*V))/C1-sigmaV*V;

dydt=[dydt(1); dydt(2); dydt(3); dydt(4); dydt(5); dydt(6); dydt(7); dydt(8); dydt(9); dydt(10); dydt(11); dydt(12);  dydt(13);  dydt(14)];
end










