% set different design points
function main()
%thk_row = linspace(0.8,1.0,21);
thk_row = linspace(0.5,1.5,11);
% thk_row = linspace(2.0,2.0,1);
n = length(thk_row);
y = zeros(n,1);
adj = zeros(n,1);
for ithk = 1: n
    thk = thk_row(ithk);
    disp(thk);
    
    %% create new desvar.inp file
    fid = fopen('desvars.inp','wt');
    fprintf(fid, 'DESVAR  30      T.C       %4.3f', thk);
    fid = fopen('desvars1.inp','wt');
    fprintf(fid, 'DESVAR  30      T.C       %4.3f     YES', thk);
    %% call FLUT_thk_NILSS_FullAdjoint.m
    [y(ithk), adj(ithk)] = FLUT_thk_NILSS_FullAdjoint();
end
save('yAdj')


% % % adjustment
% % adj = adj * 0.5*1.4*0.25^2;
%% plot
figure;
plot(thk_row, y);
hold on;
scale = (max(y)-min(y))/(thk_row(end)-thk_row(1));
for i = 1: length(thk_row)
    width = (thk_row(2) - thk_row(1) )/ 2 / sqrt(1+(adj(i)/scale)^2);
    xx = linspace(thk_row(i)-width, thk_row(i)+width,3);
    yy = y(i)+adj(i)*(xx-thk_row(i));
plot(xx,yy)
end

end
