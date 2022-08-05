
%This code runs a 2-dimensional of kinetic models in a continuous-flow simulation.
%each experiment with a unique combination of sorbent concentration (Cs) and initial sorbate concentration (C0)



%Adsorption kinetics are modelled using the modified pseudo-second order
%rate equation, being (a) second order with respect to the relative
%proportion of adsorption capacity remaining (with maximum adsorption
%capacity determined using the Freundlich adsorption isotherm), (b) first
%order with respect to sorbate concentration, and (c) zero-/first- order
%with respect to sorbent concentration when the rate is normalised to
%sorbent concentration and total volume, respectively.%

%filename,number_exps,j_flow,C_inf,Cs0
number_exps=1;
j_flow=0.1;
C_inf=2000;
Cs0=10;

number_of_experiments=number_exps;  %each experiment will have a different initial sorbate concentration (C0) 实验次数
 
number_of_variables = 8;            %counting how many rows in the array we need to store the input parameters for each kinetic plot 参数数量

experiments = zeros(number_of_experiments,number_of_variables);   %each experiment refers to a single kinetic plot  54*8的矩阵（行为第几次试验，列表示参数）

%creating the input variables for each kinetic plots

C_init = 1;     % ug L-1  - this is the initial concentration of aqueous sorbate in the suspension 这是悬浮液中山梨醇水溶液的初始浓度

q_init = 0;     % ug L-1  - this is the initial concentration of adsorbed sorbate in the suspension 吸收AS的浓度

Cs = Cs0;       %g L-1 - this is the concentration of sorbent  吸附剂浓度

Cinfluent = C_inf;  %ppb or ug L-1 - this is the concentration of sorbate in the influent (for continuous-flow modelling)

j = j_flow;          %this the turn-over frequency, i.e. bed volumes per minute

k = 0.1111;     %this is the value of normalised k' (L g-1 min-1)

KF = 5.10;      %this is the Freundlich constant (mg g-1) used to determine 'qe' at each time step

n = 2.63;       %this is the second parameter for the Freundlich adsorption isotherm (g L-1) used to determine 'qe' at each time step


%setting the time intervals upon which data is recorded. This will be
%overwritten later to avoid wasting computational time after breakthrough
%has already occurred.
t_end = 100000;   %1440 = 1 day
t_steps = 2000;
t_step = t_end/t_steps;
bv_end = t_end*j;

%setting up the arrays where calculated data is to be stored

global exp;

%赋值
for i = 1:number_of_experiments
    experiments(i,1) = C_init;
    experiments(i,2) = q_init;
    experiments(i,3) = j;
    experiments(i,4) = Cinfluent;
    experiments(i,5) = k;
    experiments(i,6) = Cs*10^(i-1); %exponentially increasing sorbent concentration 吸附剂浓度呈指数级增长的百分比
    experiments(i,7) = KF;
    experiments(i,8) = n;   
end


results_table=[];

for i = 1:number_of_experiments
    exp = zeros(1,number_of_variables);   
    exp(1,:) = experiments(i,:);
    exp_C_init = exp(1,1);
    exp_q_init = exp(1,2);
    %making sure that we model an appropriate length of time (duration) with
    %appropriate interval lengths for each simulation.
    %making sure that we model an appropriate length of time (duration) with
    %appropriate interval lengths for each simulation.
 
    t_end = 50000 * exp(1,6) / (exp(1,3) * exp(1,4));  % 50000*Cs/(j*Cinf)
    if t_end<1000
        t_end = t_end*20;
    end
    if t_end<500
        t_end=t_end*5;
    end
    if j<0.01
        t_end=t_end*10;
    end
    %t_step = t_end/t_steps;
    stepping = [1:1:t_steps+1];      %collect data at shorter time intervals in the initial stages of the simulation
    stepping(1)=0;
    gradient_1=1;
    gradient_2=60;
    for j = 2:(t_steps+1)
    %change the time intervals from being evenly spaced to having a smooth
    %transition from gradient_1 to gradient_2
       stepping(j)=stepping(j-1)+(gradient_2*(stepping(j)/t_steps))+(gradient_1*((t_steps-stepping(j))/t_steps));
    end
    for k = 2:(t_steps+1)
        %normalise the time intervals to 1 and then multiply out by the desired
        %final time
        stepping(k)=(stepping(k)/stepping(t_steps+1))*t_end;
        %"i is " + i +" and stepping(i) is " + stepping(i)
    end
    options = odeset('BDF',1e-4,'Stats','on','OutputFcn',@odeplot);   %need to increase the tolerance to avoid errors, was 1e-5 to begin with, tried changing to 1e-4
    [t,C]=ode15s(@DiffEq,stepping,[exp_C_init exp_q_init],options);   %call the ODE function
    [numRows,numCols]=size(t);
    results_t = zeros(number_of_experiments,numRows);
    results_Ct = zeros(number_of_experiments,numRows);
    results_qt = zeros(number_of_experiments,numRows);
    results_bv = zeros(number_of_experiments,numRows);
    results_t(i,:) = t;
    %results_t(i,:) = t/speed;
    results_bv(i,:) = results_t(i,:)*exp(1,3);
    results_Ct(i,:) = C(:,1);
    results_qt(i,:) = C(:,2);




end
plot( results_t ,results_Ct)

function dCdt = DiffEq(t,conditions)
global exp;
%global speed

time=t;

j = exp(1,3);
Cinfluent = exp(1,4);
k = exp(1,5);
Cs = exp(1,6);
KF = exp(1,7);
n = exp(1,8);

Ct = conditions(1);   %ppb
if Ct<0.001
   Ct=0.001;            %set a minimum concentration of sorbate in the reactor to prevent the Freundlich adsorption isotherm from running an error.
end
if Ct>Cinfluent
    Ct=Cinfluent;       %control incase ode15s time intervals are too large
end
qt = conditions(2);   %ppb
Ct_mgL = Ct/1000;         %mg L-1
qt_mgg = qt/(Cs*1000);    %mg g-1

%"Ct is " + Ct + " ug L-1 and qt is " + qt + "ug L-1";

%the rate equation we are using is
%dqt/dt = k' * Ct * (1-(qt/(Kf * Ce^(1/n))))^2
%please see our manuscript for derivation and further information

rate_ads = 1000*k*Ct_mgL*Cs*((1-(qt_mgg/(KF*(Ct_mgL^(1/n)))))^2); 
%calculate dq/dt in ppb L-1 min-1
rate_influx = j*Cinfluent; %calculate the rate of sorbate influx (continuous-flow systems only)
rate_outflux = j*Ct; %calculate the rate of sorbate outflux (continuous-flow systems only)
dCdt=zeros(2,1);
dCdt(1)=-rate_ads+rate_influx-rate_outflux;
dCdt(2) = rate_ads; 
%dCdt = [-rate_ads+rate_influx-rate_outflux;rate_ads]; %adjust the concentration of aqueous sorbate and adsorbed sorbate, respectively
%dCdt = [speed*(-rate_ads+rate_influx-rate_outflux);speed*(rate_ads)]; %adjust the concentration of aqueous sorbate and adsorbed sorbate, respectively
end






