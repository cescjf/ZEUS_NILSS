%%%%%%%%%%%% Jack Optimization %%%%%%%%%%
% this function performs minimization on the adjoint of the whole field
function FLUT_thk_NILSS_FullAdjoint()
sdir='thksens';
resume = 0;
if resume == 0
    %% prepare the primal solution
%     scomd=sprintf('zeus dir=%s main_init.inp',sdir); 
%     dos(scomd);
%     pause(0.5);

    %% prepare the non-homogeneous solution
    % main0.inp file, homo =0, requires no input files
    dos('copy main0_nonhomo.inp main0.inp /Y');
    % compute v*
    [vstar_obj, vstar_opt] = run(2);

    %% prepare the homogeneous solutions
    nus = 2;
    w_obj = zeros(length(vstar_obj), nus);
    w_opt = zeros(length(vstar_opt), nus);
    for ius = 1: nus
        % main_0.inp file
        dos('copy main0_homo.inp main0.inp /Y');
        % copute homogeneous solutions: w
        [w_obj(:,ius), w_opt(:,ius)] = run(1);
    end
    % save data
    save('vstar_w');
elseif resume == 1
    load('vstar_w');
    scomd=sprintf('zeus dir=%s main0.inp',sdir); 
    dos(scomd);            
    pause(0.5);
end

%% construct left matrix and right hand side, compute v
disp('...construct left matrix and right hand side, compute v');
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
% save data
save('v_vstar_w');

%% read in coefficients
disp('read in coefficients')
%READ IN FUM0 AND FLM0
fid=fopen('ADJAFOIL2.DAT','r');   
linn=fgets(fid);
nht=sscanf(linn,'%d');
for n=1:nht
linn=fgets(fid);
nsnc=sscanf(linn,'%d');
FUM0=zeros(nsnc(3),nsnc(2));
FLM0=zeros(nsnc(3),nsnc(2));
FUY0=zeros(nsnc(3),nsnc(2));
FLY0=zeros(nsnc(3),nsnc(2));     
AREA0=zeros(nsnc(3),nsnc(2));  
for ns=1:nsnc(2)
    for np=1:nsnc(3)
        linn=fgets(fid);
        xy=sscanf(linn,'%f');
        FUM0(np,ns)=xy(1); 
        FLM0(np,ns)=xy(2);   
        FUY0(np,ns)=xy(3); 
        FLY0(np,ns)=xy(4);  
        AREA0(np,ns)=xy(5);            
    end
end
end
fclose(fid);

% read in steps
fid=fopen('ADJOUTP_UNST.DAT','r');
linn=fgets(fid);
xy=sscanf(linn,'%d %d %d %f');
NSTEP=xy(1);
NSTEP=NSTEP+1; %INCLUDING THE INITIAL STEADY ADJOINT
LSTEP=5001;
fclose(fid)

% read in FUM and FLM, first run main00.inp 
scomd=sprintf('zeus dir=%s main00.inp',sdir); 
dos(scomd);  
pause(0.5);        

fid=fopen('ADJAFOIL2.DAT','r');          
linn=fgets(fid);
nht=sscanf(linn,'%d');
for n=1:nht            %THIS LOOP IS NOT ACTIVATED YET
    linn=fgets(fid);
    nsnc=sscanf(linn,'%d');
    FUM=zeros(nsnc(3),nsnc(2));
    FLM=zeros(nsnc(3),nsnc(2));
    FUY=zeros(nsnc(3),nsnc(2));
    FLY=zeros(nsnc(3),nsnc(2)); 
    AREA=zeros(nsnc(3),nsnc(2));              
    for ns=1:nsnc(2)
        for np=1:nsnc(3)
            linn=fgets(fid);
            xy=sscanf(linn,'%f');
            FUM(np,ns)=xy(1); 
            FLM(np,ns)=xy(2);  
            FUY(np,ns)=xy(3); 
            FLY(np,ns)=xy(4);  
            AREA(np,ns)=xy(5);                         
        end
    end
end
fclose(fid);       

%% compute derivative
disp('compute derivative')
% reshape VAL and VAU
v = reshape(v,[],2);
VAL = reshape(v(:, 1), nsnc(3),nsnc(2), LSTEP);
VAU = reshape(v(:, 2), nsnc(3),nsnc(2), LSTEP);

% prepare coefficients
ns=1:nsnc(2);
np=1:nsnc(3);

adjgradi=zeros(LSTEP,1);       
adjgradi2=zeros(LSTEP,1);            
dnxu=FUM(np,ns)-FUM0(np,ns); 
dnxl=FLM(np,ns)-FLM0(np,ns);  
dnyu=FUY(np,ns)-FUY0(np,ns); 
dnyl=FLY(np,ns)-FLY0(np,ns); 
darea=AREA(np,ns)-AREA0(np,ns);          

% compute derivative
for istep=1:LSTEP        
    adjgu=dnxu.*VAU(np,ns,istep); %+ dnyu.*VBU(np,ns,istep);
    adjgl=dnxl.*VAL(np,ns,istep); %+ dnyl.*VBL(np,ns,istep);     
    %adjgu2=darea.*VCU(np,ns,istep);
    %adjgl2=darea.*VCL(np,ns,istep);             
    adjgradi(istep)=(sum(sum(adjgu))+sum(sum(adjgl))); 
    %adjgradi2(istep)=(sum(sum(adjgu2))+sum(sum(adjgl2))); 
end


fid=fopen('GRAD_GRID.DAT','r');          
linn=fgets(fid);
xy=sscanf(linn,'%f');
DV=xy(1);

% window function
t = linspace(0,1,LSTEP);
window = 2 * (sin(t*pi) .^2);
adjgrad = sum( (adjgradi+adjgradi2).*window' )/DV    

save('everything');
end

% function V = run()
%     V = [0;0];
% end

%% function run
function [V_obj, V_opt] = run(homo)
% homo = 0: non-homo with 0 terminal condition
% homo = 1: homo
V_opt = [];
% need an empty output file
dos('del ADJSOL_OUT.dat');            
dos('copy ADJSOL_EMPTY.dat ADJSOL_OUT.dat /Y');

% prepare random ADJSOL_IN.dat file
fid = fopen('ADJSOL_IN.dat','wt');
for i = 1: 2
    temp = rand(1, 6);
    if homo == 2
        temp = zeros(1, 6);
    end
    V_opt = [V_opt; temp'];
    fprintf(fid, '  %18.11e', temp);
    fprintf(fid, '\n');
end
for ii = 1: 15456
    temp = rand(1, 5);
    if homo == 2
        temp = zeros(1, 5);
    end
    V_opt = [V_opt; temp'];
    fprintf(fid,'  %18.11e', temp);
    fprintf(fid, '\n');            
end
fclose(fid);
% if homo is 0 then the ADJSOL_IN is meaningless
% if homo == 0
%     V_opt = zeros(size(V_opt));
% end

%GET TARGET AIRFOIL Cp
format long;
sdir='thksens';
scomd=sprintf('zeus dir=%s main0.inp',sdir); 
dos(scomd);            
pause(0.5);

% read in coefficients
fid=fopen('ADJAFOIL2.DAT','r');          
linn=fgets(fid);
nht=sscanf(linn,'%d');
for n=1:nht            %THIS LOOP IS NOT ACTIVATED YET
    linn=fgets(fid);
    nsnc=sscanf(linn,'%d');
end
fclose(fid);       

%READ IN values
fid=fopen('ADJOUTP_UNST.DAT','r');  
linn=fgets(fid);
xy=sscanf(linn,'%d %d %d %f');
NSTEP=xy(1);
NSTEP=NSTEP+1; %INCLUDING THE INITIAL STEADY ADJOINT
LSTEP = 5001;
%LSTEP=xy(2);
%LSTEP=LSTEP+1; %INCLUDING THE INITIAL STEADY ADJOINTLAEROL=xy(3);
%DT=xy(4);
%NBEG=xy(6)+1;
VAL=zeros(nsnc(3),nsnc(2),(LSTEP));
VAU=zeros(nsnc(3),nsnc(2),(LSTEP));
% VBL=zeros(nsnc(3),nsnc(2),(LSTEP));
% VBU=zeros(nsnc(3),nsnc(2),(LSTEP));
% VCL=zeros(nsnc(3),nsnc(2),(LSTEP));
% VCU=zeros(nsnc(3),nsnc(2),(LSTEP));
num  = 1;
for istep=LSTEP:-1:1
for ns=1:nsnc(2)
    for np=1:nsnc(3)
        num = num + 1;
        linn=fgets(fid);
        xy=sscanf(linn,'%f');
        VAL(np,ns,istep)=xy(1); 
        VAU(np,ns,istep)=xy(2);   
%         VBL(np,ns,istep)=xy(3); 
%         VBU(np,ns,istep)=xy(4);     
%         VCL(np,ns,istep)=xy(5); 
%         VCU(np,ns,istep)=xy(6);          
    end
end
end
fclose(fid);
V_obj = [reshape(VAL,[],1); reshape(VAU,[],1)];

%READ IN values for optimization: full adjoint solution at t = 0,end
fid=fopen('ADJSOL_OUT.dat','r');  
while ~feof(fid)
    linn=fgets(fid);
    xy=sscanf(linn,'%f');
    V_opt = [V_opt; xy];
end
fclose(fid);

end
