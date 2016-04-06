% plot the diverging adjoint and the adjoint solution w.r.t. nus

close all
clear all
load('everything')

%% plot the diverging adjoint
% reshape VAL and VAU
v = reshape(vstar_obj,[],2);
VAL = reshape(v(:, 1), nsnc(3),nsnc(2), LSTEP);
VAU = reshape(v(:, 2), nsnc(3),nsnc(2), LSTEP);

% compute derivative
ns=1:nsnc(2);
np=1:nsnc(3);
adjgradi=zeros(LSTEP,1); 
for istep=1:LSTEP        
    adjgu=dnxu.*VAU(np,ns,istep); 
    adjgl=dnxl.*VAL(np,ns,istep);             
    adjgradi(istep)=(sum(sum(adjgu))+sum(sum(adjgl))); 
end

adjgrad = adjgradi / DV;
figure;
semilogy(abs(adjgrad));
xlim([0,5000]);
xlabel('time');
ylabel('dJ/ds');
savefig('divergingAdj.fig');
saveas(gcf,'divergingAdj.png')


%% plot adj w.r.t. nus
load('everything')
w_opt_all = w_opt;
w_obj_all = w_obj;
adjgrad_nus = zeros(1,5);
for nus = 1: 5
    
    w_opt = w_opt_all(:,1:nus);
    w_obj = w_obj_all(:,1:nus);
    M = zeros(nus, nus);
    for i = 1: nus
        for j = 1: nus
            M(i,j) = dot( w_opt(:,i) , w_opt(:,j));
        end
    end
    rhs = zeros(nus,1);
    for i = 1:nus
        rhs(i) = dot(vstar_opt, w_opt(:,i));
    end  
    v_opt = vstar_opt - ((M\rhs)' * w_opt')';
    v_obj = vstar_obj - ((M\rhs)' * w_obj')';
    v = v_obj;

    % reshape VAL and VAU
    v = reshape(v,[],2);
    VAL = reshape(v(:, 1), nsnc(3),nsnc(2), LSTEP);
    VAU = reshape(v(:, 2), nsnc(3),nsnc(2), LSTEP);


    % compute derivative
    ns=1:nsnc(2);
    np=1:nsnc(3);
    adjgradi=zeros(LSTEP,1); 
    for istep=1:LSTEP        
        adjgu=dnxu.*VAU(np,ns,istep);
        adjgl=dnxl.*VAL(np,ns,istep);             
        adjgradi(istep)=(sum(sum(adjgu))+sum(sum(adjgl))); 
    end

    % window function
    t = linspace(0,1,LSTEP);
    window = 2 * (sin(t*pi) .^2);
    adjgrad = sum( (adjgradi+adjgradi2).*window' )/DV;

    adjgrad_nus(nus) = adjgrad;
end

figure;
semilogy(linspace(1,5,5), abs(adjgrad_nus))
savefig('adj_nus')

