function examine_embedded_velocity(zValues)

figure(100)
plot(zValues)

zvel = sqrt(sum(abs(diff(zValues)).^2,2));

figure(97)
subplot(2,1,1)

plot(zvel)
subplot(2,1,2)

hist(zvel,logspace(-4,4,1000))
set(gca,'XScale','log')

badvals = find(zvel>30);
goodvals = find(zvel<30);

figure(101)
subplot(1,3,1)
plot(zValues(badvals(1:end),1),zValues(badvals(1:end),2),'bo','MarkerSize',1)
hold on
plot(zValues(goodvals(1:end),1),zValues(goodvals(1:end),2),'ro','MarkerSize',1)
hold off
subplot(1,3,2)
plot(zValues(badvals(1:end),1),zValues(badvals(1:end),2),'bo','MarkerSize',1)
subplot(1,3,3)
plot(zValues(goodvals(1:end),1),zValues(goodvals(1:end),2),'ro','MarkerSize',1)
