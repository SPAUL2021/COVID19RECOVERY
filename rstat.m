
nd = 14;

for i = 1:nd
    tau(i) = i; % (1/3)+(i-1);
end   

ns = 29;

for i = 1:ns
    sigma(i) = 14+(i-1); % + (1/3);
end  

data = load('rts.txt');

for i = 1:nd
    for j = 1:ns
        rptm(j,i) = data(ns*(i-1)+j);
    end
end    

dist = load('rsig.txt');

sumd = 0;
for j = 1:ns
    sumd = sumd+dist(j);
end   

dens = (1/sumd)*dist;

k = 18.62067;
theta = 1.18892;

x1 = 14:0.1:50;
y1 = distibu(x1,k,theta);

k2 = 34.55447;
theta2 = 0.60847;

x2 = 14:0.1:50;
y2 = distibu(x2,k2,theta2);

k3 = 226.40545;
theta3 = 0.17171;

x3 = 14:.1:50;
y3 = distibu(x3,k3,theta3);

lambda = 0.9365;

data = load('rdist.txt');

sd11 = size(data)

for i = 1:5587
    data2(i) = data(i);
end

for i = 1:368
    data3(i) = data(5575+i);
end    

params(1) = k;
params(2) = theta;
[~,ncov] = gamlike(params,data)

params2(1) = k2;
params2(2) = theta2;
[~,ncov2] = gamlike(params2,data2)

params3(1) = k3;
params3(2) = theta3;
[~,ncov3] = gamlike(params3,data3)

p = 0:0.0125:0.999;
[xp,xpl,xpu] = gaminv(p,k,theta,ncov);
[xp2,xp2l,xp2u] = gaminv(p,k2,theta2,ncov2);
[xp3,xp3l,xp3u] = gaminv(p,k3,theta3,ncov3);

xpbi = lambda*xp2+(1-lambda)*xp3;
xpbil = lambda*xp2l+(1-lambda)*xp3l;
xpbiu = lambda*xp2u+(1-lambda)*xp3u;

figure(2)
subplot(2,3,1)
imagesc(tau,sigma,rptm)
colormap(jet );
 colorbar
  set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(a)','FontSize',16);
xlabel('Incubation period (days)','FontSize',16);
ylabel('Recovery period (days)','FontSize',16);

subplot(2,3,[2,3])
w1 = 0.7;
bar(sigma,dens, w1, 'FaceColor',[0.2 0.2 0.5])
hold on
plot(x1,y1,'k','LineWidth',3)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 legend('Model','Gamma distribution')
 title('(b)','FontSize',16);
xlabel('Recovery period (days)','FontSize',16);
ylabel('Density','FontSize',16);

subplot(2,3,[4,5])
w1 = 0.7;
bar(sigma,dens, w1, 'FaceColor',[0.2 0.2 0.5])
hold on
plot(x3,0.9365*y2+0.0635*y3,'r-','LineWidth',3)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
  legend('Model','Bimodal gamma distribution')
 title('(c)','FontSize',16);
xlabel('Recovery period (days)','FontSize',16);
ylabel('Density','FontSize',16);

subplot(2,3,6)
plot(xp,100*p,'k-',xpbi,100*p,'r-','LineWidth',3)
legend('Unimodal','Bimodal')
hold on
plot(xpl,100*p,'k-',xpu,100*p,'k-','LineWidth',3)
hold on
plot(xpbil,100*p,'r-',xpbiu,100*p,'r-','LineWidth',3)
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(d)','FontSize',16);
xlabel('Recovery period (days)','FontSize',16);
ylabel('Percentile','FontSize',16);



function g = distibu(x,k,theta)

        g = gampdf(x,k,theta);
end        
       
       
       

