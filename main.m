% set different design points
function main()
thk_row = linspace(0.5,1.3,9);
% thk_row = linspace(0.5,0.6,2);
n = length(thk_row);
y = zeros(n,1);
adj = zeros(n,1);
devi = zeros(n,1);
for ithk = 1: n
    thk = thk_row(ithk);
    disp(thk);
    
    %% create new desvar.inp file
    fid = fopen('desvars.inp','wt');
    fprintf(fid, 'DESVAR  30      T.C       %4.3f', thk);
    fid = fopen('desvars1.inp','wt');
    fprintf(fid, 'DESVAR  30      T.C       %4.3f     YES', thk);
    %% call FLUT_thk_NILSS_FullAdjoint.m
    [y(ithk), adj(ithk), devi(ithk)] = FLUT_thk_NILSS_FullAdjoint(5001);
end
save('yAdj')

%% plot
figure;
p1 = plot(thk_row, y, 'sb','MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',5);
hold on;
scale = (max(y)-min(y))/(thk_row(end)-thk_row(1));
for i = 1: length(thk_row)
    width = (thk_row(2) - thk_row(1) )/ 2 / sqrt(1+(adj(i)/scale)^2);
    xx = linspace(thk_row(i)-width, thk_row(i)+width,3);
    yy = y(i)+adj(i)*(xx-thk_row(i));
    p2 = plot(xx,yy,'r');
end

% plot deviation bar
devi = devi * 3; % 95 confidence interval
for i = 1: length(thk_row)
    xx = ones(2,1) * thk_row(i);
    yy = [y(i) + devi(i), y(i)-devi(i)];
    plot(xx,yy,'k')
    
    width = (thk_row(2) - thk_row(1)) * 0.1;
    xx = [thk_row(i) - width , thk_row(i) + width];
    yy = [y(i) + devi(i), y(i)+devi(i)];
    plot(xx,yy,'k')
    
    xx = [thk_row(i) - width , thk_row(i) + width];    
    yy = [y(i) - devi(i), y(i)-devi(i)];
    plot(xx,yy,'k')
end
xlabel('thickness');
ylabel('averaged C_L^2');
ylim([3e-4,6e-4]);

% plot the linear regression
B = polyfit(thk_row', y, 1);
xx = linspace(thk_row(1), thk_row(end), 100);
yy = B(1) * xx + B(2);
p3 = plot(xx,yy,'--');
legend([p1,p2,p3],'average objective','sensitivity by NI-LSS','linear regression')

savefig('adj_thk.fig');
saveas(gcf,'adj_thk.png')
end
