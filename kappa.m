

prelim = load('canadaD2.txt');
d11 = size(prelim) 

  pp = abs(prelim)   
    
  nd = 14;
  ns = 29;
  spar = zeros(nd,ns);
for ii = 1:nd   

for i = 1:ns
    spar(ii,i) =  pp((ii - 1)*ns + i);
end    

end

s0 = min(spar)
s00 = min(s0)
s9 = max(spar)
s99 = max(s9)

s = spar;

s1 = size(spar)


figure(1)
axes('position', [0.2  0.2 0.71 0.38]);
imagesc(spar)
colorMap = jet; %colorcube;
colormap(colorMap);
colorbar;
set(gca,'LineWidth',2,'FontSize',20,'Box','on');
xlabel('Column','FontSize',20)
ylabel('Row','FontSize',20)
