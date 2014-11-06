log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun2'};

for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(8:11)
        fid = fopen([log_path 'HBA_' func_name '_' num2str(N) '.log'],'a+');
        for EPS = [2:2:10]
            runCbs(N, func_name, EPS, fid);
            fprintf('Func %s, N %4d, EPS %2d finished.\n',func_list{func_i},N,EPS);
        end
        fclose(fid);
    end
end