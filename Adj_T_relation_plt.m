% plot adj w.r.t T
close all
clear all
dbstop if error

adjgrad_t = zeros(1,10);
for T = 1: 10
    %% generate new main0_homo file
    % Read txt into cell A
    fid = fopen('main0_homo.inp','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{36} = sprintf('DESADJ  100     max      CL              %0.1f           26.0     1', T + 15);
    % Write cell A into txt
    fid = fopen('main0_homo.inp', 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end

    %% generate new main0_nonhomo file
    % Read txt into cell A
    fid = fopen('main0_nonhomo.inp','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{36} = sprintf('DESADJ  100     max      CL              %0.1f            26.0    2', T + 15);
    % Write cell A into txt
    fid = fopen('main0_nonhomo.inp', 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end

    LSTEP = 500 * (11-T) + 1;
    [yy, adjgrad] = FLUT_thk_NILSS_FullAdjoint(LSTEP);
    adjgrad_t(T) = adjgrad;
end
save('adj_T_relation')
plot(linspace(10,1,10), adjgrad)
savefig('adj_T')



