% -------------------------------------------------------------------------
% KF-GINS-Matlab: An EKF-based GNSS/INS Integrated Navigation System in Matlab
%
% Copyright (C) 2024, i2Nav Group, Wuhan University
%
%  Author : Liqiang Wang
% Contact : wlq@whu.edu.cn
%    Date : 2023.3.2
% -------------------------------------------------------------------------

clear;
clc;
% add function to workspace
addpath("function\");

%% define parameters and importdata process config
param = Param();
%cfg = ProcessConfig1();
% cfg = ProcessConfig2();
% cfg = ProcessConfig3();
cfg = ProcessConfig4();
%cfg = ProcessConfig6();
%% importdata data
% imudata
imudata = importdata(cfg.imufilepath);
imustarttime = imudata(1, 1);
imuendtime = imudata(end, 1);

% gnss data
gnssdata = importdata(cfg.gnssfilepath);
gnssdata(:, 2:3) = gnssdata(:, 2:3) * param.D2R;
if (size(gnssdata, 2) < 13)
    cfg.usegnssvel = false;
end
gnssstarttime = gnssdata(1, 1);
gnssendtime = gnssdata(end, 1);

% odo data
if (cfg.useodonhc)
    ododata = importdata(cfg.odofilepath);
end


%% save result
navpath = [cfg.outputfolder, '/NavResult_100HzTEST'];
if cfg.usegnssvel
    navpath = [navpath, '_GNSSVEL'];
    disp("use GNSS velocity!");
end
if cfg.useodonhc
    navpath = [navpath, '_ODONHC'];
    disp("use ODO velocity!");
end
navpath = [navpath, '.nav'];
navfp = fopen(navpath, 'wt');

imuerrpath = [cfg.outputfolder, '/ImuError.txt'];
imuerrfp = fopen(imuerrpath, 'wt');

stdpath = [cfg.outputfolder, '/NavSTD.txt'];
stdfp = fopen(stdpath, 'wt');


%% get process time
% start time and end time
if imustarttime > gnssstarttime
    starttime = imustarttime;
else
    starttime = gnssstarttime;
end
if imuendtime > gnssendtime
    endtime = gnssendtime;
else
    endtime = imuendtime;
end
if cfg.starttime < starttime
    cfg.starttime = starttime;
end
if cfg.endtime > endtime
    cfg.endtime = endtime;
end

if cfg.useodonhc
    % epoch to get odo vel
    EPOCH_TO_GETVEL = 20;
    ododatarate = 1.0 / mean(diff(ododata(:, 1)));
    if cfg.odoupdaterate > ododatarate / EPOCH_TO_GETVEL
        cfg.odoupdaterate = ododatarate / EPOCH_TO_GETVEL;
        disp("warning: set ODO udpaterate to " + num2str(cfg.odoupdaterate) + "Hz!");
    end

    % odo update time
    updateinterval = 1.0 / cfg.odoupdaterate;
    time_to_nextupdate = updateinterval - mod(cfg.starttime, updateinterval);
    odoupdatetime = cfg.starttime + time_to_nextupdate;
end

% data in process interval
imudata = imudata(imudata(:,1) >= cfg.starttime, :);
imudata = imudata(imudata(:,1) <= cfg.endtime, :);
gnssdata = gnssdata(gnssdata(:, 1) >= cfg.starttime, :);
gnssdata = gnssdata(gnssdata(:, 1) <= cfg.endtime, :);

%% ================= 核心修改：模拟 GNSS 失锁 (GNSS Outages) =================
    % 定义失锁时间段 [开始时间TOW, 结束时间TOW]
    % 注意：你的数据起步是 268858，所以我给你设定了下面两个测试区间
    outages = [
        %270370, 270430;
        269200, 269260;  % 第一次失锁：起步 5 分钟后，断网 60 秒 (模拟长隧道)
        %270000, 270030   % 第二次失锁：起步 19 分钟后，断网 30 秒 (模拟立交桥)
    ];

    % 遍历所有设定的失锁区间，把落在区间内的 GNSS 数据无情删掉！
    for i = 1:size(outages, 1)
        outage_start = outages(i, 1);
        outage_end = outages(i, 2);
        
        % 找出不在这个失锁区间内的数据（保留正常的，剔除失锁的）
        valid_idx = (gnssdata(:, 1) < outage_start) | (gnssdata(:, 1) > outage_end);
        gnssdata = gnssdata(valid_idx, :);
    end
    disp(['已人为注入 ', num2str(size(outages, 1)), ' 个 GNSS 失锁区间段！']);
    % =========================================================================


    %% ================= 核心修改：IMU 降采样实验 =================
    % 容错设计：如果配置文件里忘了写 downsample_factor，默认设为 1 (不降采样)
    if ~isfield(cfg, 'downsample_factor')
        cfg.downsample_factor = 1; 
    end

    if cfg.downsample_factor > 1
        factor = cfg.downsample_factor;
        old_len = size(imudata, 1);
        new_len = floor(old_len / factor);
        
        % 预分配内存，加速运行
        new_imudata = zeros(new_len, 7); 
        
        for k = 1:new_len
            idx_start = (k-1) * factor + 1;
            idx_end = k * factor;
            
            % 1. 时间戳：取当前压缩窗口的最后一个历元时间
            new_imudata(k, 1) = imudata(idx_end, 1);
            
            % 2. 惯性数据：因为 KF-GINS 输入的是角度/速度增量，必须用 sum() 严格累加能量！
            new_imudata(k, 2:4) = sum(imudata(idx_start:idx_end, 2:4), 1);
            new_imudata(k, 5:7) = sum(imudata(idx_start:idx_end, 5:7), 1);
        end
        
        % 偷梁换柱：用降采样后的数据覆盖原内存矩阵
        imudata = new_imudata;
        
        disp(['⚠️ 警告：正在执行 IMU 降采样实验！降采样倍率：', num2str(factor), ...
              '，等效频率约：', num2str(100/factor), 'Hz']);
    end
    % =================================================================
%% for debug
disp("Start GNSS/INS Processing!");
lastprecent = 0;


%% initialization 
[kf, navstate] = Initialize(cfg);
laststate = navstate;

% data index preprocess
lastimu = imudata(1, :)';
thisimu = imudata(1, :)';
imudt = thisimu(1, 1) - lastimu(1, 1);
gnssindex = 1;
while gnssdata(gnssindex, 1) < thisimu(1, 1)
    gnssindex = gnssindex + 1;
end

if cfg.useodonhc
    odoindex = 1;
    while ododata(odoindex, 1) < thisimu(1, 1) && odoindex < size(ododata, 1)
        odoindex = odoindex + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN PROCEDD PROCEDURE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================== 新增：开始纯算法计时 ==================
tic_algorithm = tic; 
% ==========================================================

for imuindex = 2:size(imudata, 1)-1


    %% set value of last state
    lastimu = thisimu;
    laststate = navstate;
    thisimu = imudata(imuindex, :)';
    imudt = thisimu(1, 1) - lastimu(1, 1);


    %% compensate IMU error
    thisimu(2:4, 1) = (thisimu(2:4, 1) - imudt * navstate.gyrbias)./(ones(3, 1) + navstate.gyrscale);
    thisimu(5:7, 1) = (thisimu(5:7, 1) - imudt * navstate.accbias)./(ones(3, 1) + navstate.accscale);

    
    %% adjust GNSS index
    while (gnssindex <= size(gnssdata, 1) && gnssdata(gnssindex, 1) < lastimu(1, 1))
        gnssindex = gnssindex + 1;
    end
    % check whether gnss data is valid
    if (gnssindex > size(gnssdata, 1))
        disp('GNSS file END!');
        break;
    end

    %% determine whether gnss update is required
    if lastimu(1, 1) == gnssdata(gnssindex, 1)
        % do gnss update for the current state
        thisgnss = gnssdata(gnssindex, :)';
        kf = GNSSUpdate(navstate, thisgnss, kf, cfg.antlever, cfg.usegnssvel, lastimu, imudt);
        [kf, navstate] = ErrorFeedback(kf, navstate);
        gnssindex = gnssindex + 1;
        laststate = navstate;
        
        % do propagation for current imu data
        imudt = thisimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, thisimu);
        kf = InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
    elseif (lastimu(1, 1) < gnssdata(gnssindex, 1) && thisimu(1, 1) > gnssdata(gnssindex, 1))
        % ineterpolate imu to gnss time
        [firstimu, secondimu] = interpolate(lastimu, thisimu, gnssdata(gnssindex, 1));
        % NOTE：内插之后采样间隔会变化，严格上不满足INSMech的假设，但影响较小暂时忽略
        
        % do propagation for first imu
        imudt = firstimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, firstimu);
        kf = InsPropagate(navstate, firstimu, imudt, kf, cfg.corrtime);

        % do gnss update
        thisgnss = gnssdata(gnssindex, :)';
        kf = GNSSUpdate(navstate, thisgnss, kf, cfg.antlever, cfg.usegnssvel, firstimu, imudt);
        [kf, navstate] = ErrorFeedback(kf, navstate);
        gnssindex = gnssindex + 1;
        laststate = navstate;
        lastimu = firstimu;

        % do propagation for second imu
        imudt = secondimu(1, 1) - lastimu(1, 1);
        navstate = InsMech(laststate, lastimu, secondimu);
        kf = InsPropagate(navstate, secondimu, imudt, kf, cfg.corrtime);
    else
        %% only do propagation
        % INS mechanization
        navstate = InsMech(laststate, lastimu, thisimu);
        % error propagation
        kf = InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
    end


    if cfg.useodonhc
        %% update odo index
        while ododata(odoindex, 1) < thisimu(1, 1) && odoindex < size(ododata, 1)
            odoindex = odoindex + 1;
        end

        %% odonhc udpate
        if (thisimu(1, 1) >= odoupdatetime)
            startindex = odoindex - round(EPOCH_TO_GETVEL / 2);
            endindex = odoindex + round(EPOCH_TO_GETVEL / 2);
            if (startindex < 1)
                startindex = 1;
            end
            if (endindex > size(ododata, 1))
                endindex = size(ododata, 1);
            end
           
            % get odovel and update
            [odovel, valid] = GetOdoVel(ododata(startindex:endindex, :), thisimu(1, 1));
            if valid
                odonhc_vel = [odovel; 0; 0];
                kf = ODONHCUpdate(navstate, odonhc_vel, kf, cfg, thisimu, imudt);
                [kf, navstate] = ErrorFeedback(kf, navstate);
            end
            odoupdatetime = odoupdatetime + 1 / cfg.odoupdaterate;
        end
    end


    %% save data
    % write navresult to file
    nav = zeros(11, 1);
    nav(2, 1) = navstate.time;
    nav(3:5, 1) = [navstate.pos(1) * param.R2D; navstate.pos(2) * param.R2D; navstate.pos(3)];
    nav(6:8, 1) = navstate.vel;
    nav(9:11, 1) = navstate.att * param.R2D;
    fprintf(navfp, '%2d %12.6f %12.8f %12.8f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', nav);

    % write imu error, convert to common unit
    imuerror = zeros(13, 1);
    imuerror(1, 1) = navstate.time;
    imuerror(2:4, 1) = navstate.gyrbias * param.R2D * 3600;
    imuerror(5:7, 1) = navstate.accbias * 1e5;
    imuerror(8:10, 1) = navstate.gyrscale * 1e6;
    imuerror(11:13, 1) = navstate.accscale * 1e6;
    fprintf(imuerrfp, '%12.6f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', imuerror);

    % write state std, convert to common unit
    std = zeros(1, 22);
    std(1) = navstate.time;
    for idx=1:21
        std(idx + 1) = sqrt(kf.P(idx, idx));
    end
    std(8:10) = std(8:10) * param.R2D;
    std(11:13) = std(11:13) * param.R2D *3600;
    std(14:16) = std(14:16) * 1e5;
    std(17:22) = std(17:22) * 1e6;
    fprintf(stdfp, '%12.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n', std);


    %% print processing information
    if (imuindex / size(imudata, 1) - lastprecent > 0.01) 
        disp("processing " + num2str(floor(imuindex * 100 / size(imudata, 1))) + " %!");
        lastprecent = imuindex / size(imudata, 1);
    end
end

% ================== 新增：结束计时并打印 ==================
algorithm_time = toc(tic_algorithm);  
disp('--------------------------------------------------');
disp(['✅ 主解算循环执行完毕！']);
if isfield(cfg, 'downsample_factor')
    disp(['⚙️ 当前 IMU 降采样倍率: ', num2str(cfg.downsample_factor), 'x']);
end
disp(['⏳ 核心算法解算总耗时: ', num2str(algorithm_time), ' 秒']);
disp('--------------------------------------------------');
% ==========================================================

% close file
fclose(imuerrfp);
fclose(navfp);
fclose(stdfp);

disp("GNSS/INS Integration Processing Finished!")