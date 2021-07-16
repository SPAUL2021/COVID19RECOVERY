clearvars;close all;clc;

format long;

global Npop  

global EDT EDD EDR

global EDR1

global N N2

global nval

global nd

global maxnum

maxnum = 5000;

global tau
global value

global ns
global sigma
global del
global nt

global prelim
global npar

nd = 14 % total number of tau

% Time definition for model
dt = 1; % time step
time1 = datetime(2020,01,22,0,0,0):dt:datetime(2020,07,16,0,0,0);
N = numel(time1)

t = [0:N-1].*dt;

% Time definition for data
ddt = 1; % time step
time2 = datetime(2020,01,25,0,0,0):ddt:datetime(2020,07,16,0,0,0);
N2 = numel(time2)


%Generate the data
Npop= 37894799; % population 
%Q0 = 200; % Initial number of infectious that have bee quanrantined
%I0 = par0(1)*200; % Initial number of infectious cases non-quarantined
%R0 = 0; % Initial number of recovereds
%D0 = 0; % Initial number of deads

% Loading data
data = load('caronacanada.csv');
d11 = size(data) 

DT = zeros(N2,1);
DD = zeros(N2,1);
DR = zeros(N2,1);

DT = data(:,1);
DR = data(:,2);
DD = data(:,3);

EDT = zeros(N,1);
EDD = zeros(N,1);
EDR = zeros(N,1);

for k = N-N2+1:N
    EDT(k) = DT(k-N+N2);
    EDD(k) = DD(k-N+N2);
    EDR(k) = DR(k-N+N2);
end    
dq1 = EDT;

for i = 1:nd
    tau(i) = i; % (1/3)+(i-1);
end    

%fileID = fopen('canadaF1.txt','w');

% parameters for the first half
   
      value(1) = 0.020286916248476;   
      value(2) = 0.001781131418980;  
      value(3) = 0.657996331943508;  
      value(4) = 0.000111958007666;  
      value(5) = 0.015736796934949;  
      value(6) = 0.051856278402357;
      value(7) = 0.099711856222703;  
      value(8) = 0.063723914092713;  
      value(9) = 0.035145965096083;  
      value(10) = 0.075213864896830;  
      value(11) = 0.068595920318325;  
      value(12) = 0.083835584936427;
      value(13) = 0.111475107827384; 
      value(14) = 0.014211094023449;  
      value(15) = 0.021537268880149;  
      value(16) = 0.035561405417237;  
      value(17) = 0.014441955121292;  
      value(18) = 0.005835407229078;
      value(19) = 0.035639440313039;


% Choice of a particular form for lambda(t)


   par(1) = 0.001019813091030;
   par(2) = 0.000850036398043;
   par(3) = 0.733476973077193;
   par(4) = 0.021796876787298;
   par(5) = 0.024538732378698;
   par(6) = 0.041286612172221;
   par(7) = 0.059309854040723;
   par(8) = 0.113739382297089;
   par(9) = 0.146577417836896;
   par(10) = 0.139072739363188;
   par(11) = 0.107811439191143;
   par(12) = 0.089887511511006;
   par(13) = 0.067546137201272;
   par(14) = 0.051510938643997;
   par(15) = 0.038336611218976;
   par(16) = 0.032482166419202;
   par(17) = 0.020557109818042;
   par(18) = 99;
  
[TT,yy] = comp(par);

p1 = size(yy)

 y1 = yy(3,:);
  y2 = yy(4,:);
   y3 = yy(5,:);
    y4 = yy(6,:);
     y5 = yy(7,:);
      y6 = yy(8,:);
       y7 = yy(9,:);
        y8 = yy(10,:);
         y9 = yy(11,:);
          y10 = yy(12,:);
           y11 = yy(13,:);
            y12 = yy(14,:);
             y13 = yy(15,:);
              y14 = yy(16,:);
              
TTa = TT;
yya = y1;
 
 %
 multi = 1;
 
 TT1 = multi*TT;
 
 T1 = multi*y1;
 T2 = multi*y2;
 T3 = multi*y3;
 T4 = multi*y4;
 T5 = multi*y5;
 T6 = multi*y6;
 T7 = multi*y7;
 T8 = multi*y8;
 T9 = multi*y9;
 T10 = multi*y10;
 T11 = multi*y11;
 T12 = multi*y12;
 T13 = multi*y13;
 T14 = multi*y14;
 
 EDR1 = zeros(N,1);
 
 for i = 4:N
     we1(i) = T1(i)/TT1(i);
     EDR1(i)= we1(i)*EDR(i);
 end
 
 EDR11 = EDR1;
 
 
 sum = 0;
   
   for i = 1:N
       Q11 = (TT1(i) - EDT(i))^2;
       sum = sum + Q11;
   end

    f0 = sqrt(sum)/N; 

fval1 = f0


ns = 29;

for i = 1:ns
    sigma(i) = 14+(i-1); % + (1/3);
end  

del = zeros(nd+ns,1);

for i = 1:nd
    del(i)=tau(i);
end
for i = nd+1:nd+ns
    del(i)=sigma(i-nd);
end    

del1 = del;

prelim = load('canadaR1d.txt');
d11 = size(data); 

 
  pp = abs(prelim);   
    
rpt= zeros(ns,nd);    
summr = 0;
sumc = 0;
for ii = 1:nd   
    nt = ii;

for i = 1:ns
    spar0(i) =  pp((ii - 1)*ns + i);
end    

[R,uu] = scomp(spar0);

MR(:,ii) = R;
rr0(ii) = R(N);

uusize = size(uu);

sumc = sumc + uu(3,:);

for j= 1:ns
rpt(j,ii) = uu(3+j,N);
end

summr = summr+MR(:,ii);
end
rr0a = rr0;
TMR = summr;

 multi = EDR(N)/TMR(N)
 
 R1 = multi*TMR;
 
 scd = sumc;
 
 sum = 0;
   
   for i = 4:N
       Q11 = (R1(i) - EDR(i))^2;
       sum = sum + Q11;
   end

    f0 = sqrt(sum)/N 
    
    D99 = load('deaths.txt');

figure(1)
subplot(2,2,1)
 plot(time1,yy(2,:),'r-','LineWidth',3)
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(d)','FontSize',16);
xlabel('Day','FontSize',16);
 ylabel('Infected','FontSize',16)
 
 subplot(2,2,2)
 plot(time1,TT1,'b-','LineWidth',3)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 hold on
  plot(time1,EDT,'ro')
legend('Model','Data')
title('(e)','FontSize',16);
xlabel('Day','FontSize',16);
 ylabel('Total Cases','FontSize',16)
 
 subplot(2,2,3)
 plot(time1,R1,'b-','LineWidth',3)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 hold on
 plot(time1,EDR,'ro')
 legend('Model','Data')
 title('(f)','FontSize',16);
 xlabel('Day','FontSize',16)
 ylabel('Recovered','FontSize',16)
 
 subplot(2,2,4)
 plot(time1,D99,'b-','LineWidth',3)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 hold on
 plot(time1,EDD,'ro')
 legend('Model','Data')
 title('(g)','FontSize',16);
 xlabel('Day','FontSize',16)
 ylabel('Deaths','FontSize',16)
 
      rr = multi*rr0;
 
 rpt1 = size(rpt);
 rpt2 = rpt;
 
 for j1 = 1:ns
 sumr = 0;
 for i1 = 1:nd
     sumr = sumr + rpt(j1,i1);
 end
     rsig(j1) = multi*sumr;
 end    
 
 for i2 = 1:nd
 sumr = 0;
 for j2 = 1:ns
     sumr = sumr + rpt(j2,i2);
 end
     rtau(i2) = multi*sumr;
 end    
 rtau1 = rtau;
 
 rptm = zeros(ns,nd);
 
 for i = 1:nd
     for j = 1:ns
         rptm(j,i) = multi*rpt(j,i);
     end
 end   

 rptm1 = rptm;
 
 fileID88 = fopen('rts.txt','w')
 fprintf(fileID88,'%12.8f\n',rptm);
 

sum1 = 0;
for i = 1:nd
    sum1 = sum1+rtau(i);
end
rsum1 = sum1

sum2 = 0;
for i = 1:ns
    sum2 = sum2+rsig(i);
end
rsum2 = sum2

fileID = fopen('rsig.txt','w') 
fprintf(fileID,'%12.8f\n',rsig);

 
%--------------------------------------------------------------------------

%---------------------------------------------------------------------------
function [TT, yy] = comp(par2)

global N N2
global tau
global alpha2 nu2 beta2 delta2 gamma2
global alpha3 nu3 beta3 delta3 gamma3
global value
global nd
global I0

alpha2 = value(1);
nu2 = value(2);
beta2 = value(3);
  % I0 = 1.0;
I0 = value(4);

alpha3 = abs(par2(1));
nu3 = abs(par2(2));
beta3 = abs(par2(3));
 
for k = 1:nd 
    delta2(k) = beta2*value(4+k);
end    

for k = 1:nd 
    delta3(k) = beta3*abs(par2(3+k));
end   

gamma2 = abs(value(4+nd+1));
gamma3 = abs(par2(3+nd+1));

tspan = [0:0.1:N-1]';
opts = ddeset('RelTol',1e-2,'AbsTol',1e-2);
sol1 = dde23(@system,tau,@hist,tspan([1, end]),opts);
y_obs = deval(sol1, tspan);

T0 = sum(y_obs(3:nd+2,:));

for k = 1:N
    kk = 9*(k-1)+k;
    TT1(k)=T0(kk);
end    

for i = 1:17
for k = 1:N
    kk = 9*(k-1)+k;
    yy(i,k)=y_obs(i,kk);
end  
end

TT = TT1;

q = size(TT);
p = size(yy);

end
%--------------------------------------------------------------------------

function s = hist(t)
% Constant history function 
global Npop
global nd
global I0

s = zeros(nd+3,1);
s(2) = I0;
s(3:nd+2) = 0;
s(2+nd+1) = 0;
s(1) = Npop - sum(s(2:nd+3));
end
% --------------------------------------------------------------------------

function dydt = system(t,y,Z)
% Differential equations 
% S = y(1), I = y(2), T1 = y(3), T2 = y(4), T3 = y(5), L = y(6) 
global Npop
global alpha2 nu2 beta2 delta2 gamma2
global alpha3 nu3 beta3 delta3 gamma3
global nd
global tau

if (t < 178)
    alpha = alpha2;
       nu = nu2;
     beta = beta2;
    delta = delta2;
    gamma = gamma2;
else
    alpha = alpha3;
       nu = nu3;
     beta = beta3;
    delta = delta3;
    gamma = gamma3;
end    

for i = 1:nd
    
if (t < tau(i))
    d(i) = 0;
else 
    d(i) = 1;
end    

end

for i = 1:nd
    w(i) = d(i)*delta(i)*Z(1,i)*Z(2,i)/Npop;
end    

ylag1 = Z(:,1);
ylag2 = Z(:,2);
ylag3 = Z(:,3);
dydt = [ -beta*y(1)*y(2)/Npop - alpha*y(1) + nu*y(2+nd+1)
          beta*y(1)*y(2)/Npop - gamma*y(2) - sum(w(1:nd))
          w(1)
          w(2) % nd arries
          w(3)
          w(4)
          w(5)
          w(6)
          w(7)
          w(8)
          w(9)
          w(10)
          w(11)
          w(12)
          w(13)
          w(14)
          alpha*y(1) - nu*y(2+nd+1)                ];
end

%%
% --------------------------------------------------------------------------
function [R, uu] = scomp(par2)

global N N2
global del
global nd
global ns
global lambda

lambda = par2;

tspan = [0:0.1:N-1]';
opts = ddeset('RelTol',1e-2,'AbsTol',1e-2);
sol1 = dde23(@ssystem,del,@shist,tspan([1, end]),opts);
y_obs = deval(sol1, tspan);

R0 = sum(y_obs(4:ns+3,:));

for k = 1:N
    kk = 9*(k-1)+k;
    R1(k)=R0(kk);
end    

for i = 1:ns+4
for k = 1:N
    kk = 9*(k-1)+k;
    uu(i,k)=y_obs(i,kk);
end  
end

R = R1;

qs = size(R);
ps = size(uu);

end
%--------------------------------------------------------------------------

function s = shist(t)
% Constant history function 
global Npop
global ns
global I0

s = zeros(ns+4,1);
s(2) = I0;
s(3:ns+3) = 0;
s(3+ns+1) = 0;
s(1) = Npop - sum(s(2:ns+4));
end
% --------------------------------------------------------------------------

function dydt = ssystem(t,y,Z)
% Differential equations 
% S = y(1), I = y(2), T1 = y(3), T2 = y(4), T3 = y(5), L = y(6) 
global Npop
global alpha2 nu2 beta2 delta2 gamma2
global alpha3 nu3 beta3 delta3 gamma3
global ns
global nd
global del
global lambda
global nt

if (t < 178)
    alpha = alpha2;
       nu = nu2;
     beta = beta2;
    delta = delta2;
    gamma = gamma2;
else
    alpha = alpha3;
       nu = nu3;
     beta = beta3;
    delta = delta3;
    gamma = gamma3;
end    

for i = 1:nd
    
if (t < del(i))
    d(i) = 0;
else 
    d(i) = 1;
end    

end

for i = nd+1:ns+nd
    
if (t < del(i))
    dw(i-nd) = 0;
else 
    dw(i-nd) = 1;
end    

end

for i = 1:nd
    w(i) = d(i)*delta(i)*Z(1,i)*Z(2,i)/Npop;
end    

for i = 1:ns
    ww(i) = dw(i)*lambda(i)*delta(nt)*Z(1,i+nd)*Z(2,i+nd)/Npop;
end    

dydt = [ -beta*y(1)*y(2)/Npop - alpha*y(1) + nu*y(3+ns+1)
          beta*y(1)*y(2)/Npop - sum(w(1:nd)) - gamma*y(2) 
          w(nt)-sum(ww(1:ns))
          ww(1) % ns arries
          ww(2)
          ww(3)
          ww(4)
          ww(5)
          ww(6)
          ww(7)
          ww(8)
          ww(9)
          ww(10)
          ww(11)
          ww(12)
          ww(13)
          ww(14)
          ww(15)
          ww(16)
          ww(17)
          ww(18)
          ww(19)
          ww(20)
          ww(21)
          ww(22)
          ww(23)
          ww(24)
          ww(25)
          ww(26)
          ww(27)
          ww(28)
          ww(29)
          alpha*y(1) - nu*y(3+ns+1)                ];
end
