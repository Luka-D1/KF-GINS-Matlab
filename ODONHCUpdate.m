% -------------------------------------------------------------------------
% KF-GINS-Matlab: An EKF-based GNSS/INS Integrated Navigation System in Matlab
%
% Copyright (C) 2024, i2Nav Group, Wuhan University
%
%  Author : Liqiang Wang
% Contact : wlq@whu.edu.cn
%    Date : 2023.3.9
% -------------------------------------------------------------------------

function kf = ODONHCUpdate(navstate, odonhc_vel, kf, cfg, thisimu, dt)

    param = Param();

    %% measurement innovation
    wib_b = thisimu(2:4, 1) / dt;
    wie_n = [param.WGS84_WIE * cos(navstate.pos(1)); 0; -param.WGS84_WIE * sin(navstate.pos(1))];
    wen_n = [navstate.vel(2) / (navstate.Rn + navstate.pos(3)); 
            -navstate.vel(1) / (navstate.Rm + navstate.pos(3)); 
            -navstate.vel(2) * tan(navstate.pos(1)) / (navstate.Rn + navstate.pos(3))];
    win_n = wie_n + wen_n;
    wnb_b = wib_b - navstate.cbn' * win_n;

    vel_pre = cfg.cbv * (navstate.cbn' * navstate.vel + skew(wnb_b) * cfg.odolever);
    Z = vel_pre - odonhc_vel;

    %% measurement equation and noise
    H = zeros(3, kf.RANK);

    % Velocity error projected from the navigation frame to the vehicle
    % frame.
    H(:, 4:6) = cfg.cbv * navstate.cbn';

    % The attitude error changes both the projected IMU velocity and the
    % rotational velocity at the odometer lever arm. The small dependence
    % of win_n on position and velocity errors is neglected, consistently
    % with the observation model in the algorithm document.
    H(:, 7:9) = -cfg.cbv * navstate.cbn' * skew(navstate.vel) ...
                  -cfg.cbv * skew(cfg.odolever) * navstate.cbn' * skew(win_n);

    % Gyroscope bias and scale-factor errors affect the angular rate used
    % in the lever-arm compensation. Accelerometer errors have no direct
    % contribution to this velocity observation.
    H(:, 10:12) = -cfg.cbv * skew(cfg.odolever);
    H(:, 16:18) = -cfg.cbv * skew(cfg.odolever) * diag(wib_b);

    R = diag(cfg.odonhc_measnoise .^ 2);

    %% update
    K = kf.P * H' / (H * kf.P * H' + R);
    kf.x = kf.x + K*(Z - H*kf.x);
    kf.P=(eye(kf.RANK) - K*H) * kf.P * (eye(kf.RANK) - K*H)' + K * R * K';

end
