%%%%%%%%%%%% Jack Optimization %%%%%%%%%%
function thicksens
%GET TARGET AIRFOIL Cp
format long;
sdir='thksens';
scomd=sprintf('zeus dir=%s main_init.inp',sdir); 
dos(scomd);       
pause(0.5);
scomd=sprintf('zeus dir=%s main0.inp',sdir); 
dos(scomd);            
    
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

%READ IN values
fid=fopen('ADJOUTP_UNST.DAT','r');  
linn=fgets(fid);
xy=sscanf(linn,'%d %d %d %f');
NSTEP=xy(1);
NSTEP=NSTEP+1; %INCLUDING THE INITIAL STEADY ADJOINT
LSTEP=xy(2);
LSTEP=LSTEP+1; %INCLUDING THE INITIAL STEADY ADJOINTLAEROL=xy(3);
DT=xy(4);
NBEG=xy(6)+1;
VAL=zeros(nsnc(3),nsnc(2),NSTEP);
VAU=zeros(nsnc(3),nsnc(2),NSTEP);
VBL=zeros(nsnc(3),nsnc(2),NSTEP);
VBU=zeros(nsnc(3),nsnc(2),NSTEP);
VCL=zeros(nsnc(3),nsnc(2),NSTEP);
VCU=zeros(nsnc(3),nsnc(2),NSTEP);
for istep=NSTEP:-1:NBEG
for ns=1:nsnc(2)
    for np=1:nsnc(3)
        linn=fgets(fid);
        xy=sscanf(linn,'%f');
        VAL(np,ns,istep)=xy(1); 
        VAU(np,ns,istep)=xy(2);   
%        VBL(np,ns,istep)=xy(3); 
%        VBU(np,ns,istep)=xy(4);     
%        VCL(np,ns,istep)=xy(5); 
%        VCU(np,ns,istep)=xy(6);          
    end
end
end
fclose(fid);

scomd=sprintf('zeus dir=%s main00.inp',sdir); 
dos(scomd);  
pause(0.5);        

        fid=fopen('ADJAFOIL2.DAT','r');   
%        c4 = onCleanup(@()fclose(fid));        
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
        
        ns=1:nsnc(2);
        np=1:nsnc(3);
        
        adjgradi=zeros(LSTEP,1);       
        adjgradi2=zeros(LSTEP,1);            
        dnxu=FUM(np,ns)-FUM0(np,ns); 
        dnxl=FLM(np,ns)-FLM0(np,ns);  
        dnyu=FUY(np,ns)-FUY0(np,ns); 
        dnyl=FLY(np,ns)-FLY0(np,ns); 
        darea=AREA(np,ns)-AREA0(np,ns);          
        for istep=1:LSTEP        
        adjgu=dnxu.*VAU(np,ns,istep); %+dnyu.*VBU(np,ns,istep);
        adjgl=dnxl.*VAL(np,ns,istep); %+dnyl.*VBL(np,ns,istep);     
%        adjgu2=darea.*VCU(np,ns,istep);
%        adjgl2=darea.*VCL(np,ns,istep);             
        adjgradi(istep)=(sum(sum(adjgu))+sum(sum(adjgl))); 
%        adjgradi2(istep)=(sum(sum(adjgu2))+sum(sum(adjgl2))); 
        end

        
        fid=fopen('GRAD_GRID.DAT','r');          
        linn=fgets(fid);
        xy=sscanf(linn,'%f');
        DV=xy(1);
        fclose(fid);     
        adjgrad=(sum(adjgradi))/DV          
        
end
