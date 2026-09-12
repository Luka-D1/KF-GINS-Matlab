% -------------------------------------------------------------------------
%  Author : Duan FangKai
%    Date : 2026.3.31
% -------------------------------------------------------------------------

function cfg = ProcessConfig4()

    param = Param();

    %% filepath
    cfg.imufilepath = 'dataset4\KF_IMU.txt';
    cfg.gnssfilepath = 'dataset4\KF_GNSS_PPP.txt';
    cfg.odofilepath = '';
    cfg.outputfolder = 'dataset4';

    %% configure
    cfg.usegnssvel = false;
    cfg.useodonhc = false;
    cfg.odoupdaterate = 100; % [Hz]

    % ======= 新增：IMU 降采样配置 =======
    cfg.downsample_factor = 10 % 1表示不降采样(100Hz)，2表示降到50Hz，5表示降到20Hz，10表示降到10Hz

    %% initial information
    cfg.starttime = 268858;
    cfg.endtime = inf;

    cfg.initpos = [40.070072333; 116.275650512; 36.7091]; % [deg, deg, m]
    cfg.initvel = [0.002; 0.004; -0.032]; % [m/s]   (NE地)
    cfg.initatt = [-0.429844; -0.937233 ; 250.298870]; % [deg]

    cfg.initposstd = [10.0 ; 10.0 ; 10.0]; %[m]
    cfg.initvelstd = [2.0 ; 2.0 ; 2.0]; %[m/s]
    cfg.initattstd = [5.0; 5.0; 0.8]; %[deg]

    cfg.initgyrbias = [0; 0; 0]; % [deg/h]
    cfg.initaccbias = [0; 0; 0]; % [mGal]
    cfg.initgyrscale = [0; 0; 0]; % [ppm]
    cfg.initaccscale = [0; 0; 0]; % [ppm]

    cfg.initgyrbiasstd = [20; 20; 20]; % [deg/h]
    cfg.initaccbiasstd = [50; 50; 50]; % [mGal]
    cfg.initgyrscalestd = [300; 300; 300]; % [ppm]
    cfg.initaccscalestd = [300; 300; 300]; % [ppm]

    cfg.gyrarw = 0.1; % [deg/sqrt(h)]
    cfg.accvrw = 0.1; % [m/s/sqrt(h)]
    cfg.gyrbiasstd = 25; % [deg/h]
    cfg.accbiasstd = 200; % [mGal]
    cfg.gyrscalestd = 3000; % [ppm]
    cfg.accscalestd = 3000; % [ppm]
    cfg.corrtime = 1; % [h]

    %% install parameters 安装参数
    cfg.antlever = [0.430; -0.340; -1.005]; % [m]
    cfg.odolever = [0; 0; 0]; %[m]
    cfg.installangle = [0; 0; 0]; %[deg]

    %% ODO/NHC measurement noise 观测噪声
    cfg.odonhc_measnoise = [0.1; 0.1; 0.1]; % [m/s]


    %% convert unit to standard unit (单位转换)
    cfg.initpos(1) = cfg.initpos(1) * param.D2R;
    cfg.initpos(2) = cfg.initpos(2) * param.D2R;
    cfg.initatt = cfg.initatt * param.D2R;

    cfg.initattstd = cfg.initattstd * param.D2R;

    cfg.initgyrbias = cfg.initgyrbias * param.D2R / 3600;
    cfg.initaccbias = cfg.initaccbias * 1e-5;
    cfg.initgyrscale = cfg.initgyrscale * 1e-6;
    cfg.initaccscale = cfg.initaccscale * 1e-6;
    cfg.initgyrbiasstd = cfg.initgyrbiasstd * param.D2R / 3600;
    cfg.initaccbiasstd = cfg.initaccbiasstd * 1e-5;
    cfg.initgyrscalestd = cfg.initgyrscalestd * 1e-6;
    cfg.initaccscalestd = cfg.initaccscalestd * 1e-6;

    cfg.gyrarw = cfg.gyrarw * param.D2R / 60;
    cfg.accvrw = cfg.accvrw / 60;
    cfg.gyrbiasstd = cfg.gyrbiasstd * param.D2R / 3600;
    cfg.accbiasstd = cfg.accbiasstd * 1e-5;
    cfg.gyrscalestd = cfg.gyrscalestd * 1e-6;
    cfg.accscalestd = cfg.accscalestd * 1e-6;
    cfg.corrtime = cfg.corrtime * 3600;

    cfg.installangle = cfg.installangle * param.D2R;
    cfg.cbv = euler2dcm(cfg.installangle);

end

