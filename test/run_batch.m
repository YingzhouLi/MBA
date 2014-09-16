log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'funF', 'fun0', 'fun1'};

for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(7:11)
        fid = fopen([log_path 'HBA_' func_name '_' num2str(N) '.log'],'a+');
        for EPS = [5 7 9 11]
            runCbs(N, func_name, EPS, fid);
        end
        fclose(fid);
    end
end