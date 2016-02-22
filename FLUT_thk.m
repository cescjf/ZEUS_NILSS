%%%%%%%%%%%% Jack Optimization %%%%%%%%%%
function FLUT_thk()
sdir='thksens';
resume = 0;
if resume == 0
    % prepare the primal solution
    scomd=sprintf('zeus dir=%s main_init.inp',sdir); 
    dos(scomd);
    pause(0.5);

    %% prepare the non-homogeneous solution
    % main0.inp file, homo =0, requires no input files
    dos('copy main0_nonhomo.inp main0.inp /Y');
    % compute v*
    vstar = run();

    % prepare the homogeneous solutions
    nus = 50;
    w = zeros(length(vstar), nus);
    for ius = 1: nus
        % main_0.inp file
        dos('copy main0_homo.inp main0.inp /Y');
        % ADJSOL_IN.dat file
        NMDS = 3;
        NXT = 139;
        NYT = 2;
        NZT = 57;
        fid = fopen('ADJSOL_IN.dat','wt');
        for i = 1: 2
            temp = rand(1, 2* NMDS);
            fprintf(fid, '  %13.7e', temp);
            fprintf(fid, '\n');
        end
        for ii = 1: 2
        for i = 2: NXT
            for j = 2: NYT
                for k = 2: NZT
                    temp = rand(1, 5);
                    fprintf(fid,'  %13.7e', temp);
                    fprintf(fid, '\n');            
                end
            end
        end
        end
        fclose(fid);
        % copute homogeneous solutions: w
        w(:,ius) = run();
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
M = zeros(nus, nus);
for i = 1: nus
    for j = 1: nus
    M(i,j) = dot( w(:,i) , w(:,j));
    end
end
rhs = zeros(nus,1);
for i = 1:nus
    rhs(i) = dot(vstar, w(:,i));
end 
v = vstar - ((M\rhs)' * w')';
% save data
save('v_vstar_w');


%% compute derivative
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
LSTEP=1501;
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


% fid=fopen('GRAD_GRID.DAT','r');          
% linn=fgets(fid);
% xy=sscanf(linn,'%f');
% DV=xy(1);
DV = 0.254000013E-03;

% window function
t = linspace(0,1,LSTEP);
window = 2 * (sin(t*pi) .^2);
adjgrad = sum( (adjgradi+adjgradi2).*window' )/DV    

save('everything');
end

% function V = run()
%     V = [0;0];
% end

function V = run()

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
LSTEP = 1501;
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

V = [reshape(VAL,[],1); reshape(VAU,[],1)];

end
