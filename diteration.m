
for i = 1:14
    x(i) = i;
end


y(1) = 51.8 ;
y(2) = 49.9 ;
y(3) = 46.2 ;
y(4) = 43.9 ;
y(5) = 42.9 ;
y(6) = 39.8 ;
y(7) = 37.14 ;
y(8) = 34.23 ;
y(9) = 29.5 ;
y(10) = 29.1 ;
y(11) = 28.33 ;
y(12) = 27.41 ;
y(13) = 26.9 ;
y(14) = 26.6 ;

plot(x,y,'b-','LineWidth',2)
axis([0 15 25 55]);
set(gca,'LineWidth',2,'FontSize',20,'Box','on');
hold on
plot(x,y,'r.','MarkerSize',25)
xlabel('Iteration','FontSize',20)
ylabel('E_D','FontSize',20)

