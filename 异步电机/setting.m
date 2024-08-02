close all;
clear;
clc;

% inverter parameters
% 直流电源电压
Vdc = 500;
% 开关频率
fsw = 10e3;
% 开关周期
Tsw = 1 / fsw;
pxyz = sqrt(3) *Tsw / Vdc;

% motor parameters
P_rated = 3700; % 额定功率
V_rated = 220; % 额定电压（线电压的有效值）
I_rated = 12.8; % 额定电流
f_rated = 50; % 额定频率
w_rated = 2 * pi * f_rated; % 额定角速度
rpm_rated = 1460; % 额定转速
np = 2; % 极对数

V_max = sqrt(2) / sqrt(3) * V_rated; % 相电压的最大值

rpm_syn = 60 * f_rated / np; % 同步转速
w_syn = rpm_syn * pi / 30; % 同步角速度
s_rated = (rpm_syn - rpm_rated) / rpm_syn; % 额定转差率


Rs = 0.295; % 定子电阻
Rr = 0.379; % 转子电阻
Lm = 59e-3; % 互感 Lm = 1.5 * Lms
Lls = 1.794e-3; % 定子漏感
Llr = 1.794e-3; % 转子漏感
Ls = Lm + Lls; % 定子电感
Lr = Lm + Llr; % 转子电感
sigma = 1 - Lm^2 / Ls / Lr; % 转子漏磁系数
R = Rs + Rr * (Lm / Lr)^2; % 等效电阻
L = sigma * Ls; % 等效电感
tau_r = Lr / Rr; % 转子时间常数

Te_rated = np * P_rated / (w_syn * (1 - s_rated)); % 额定转矩
Ids_rated =V_max / sqrt(Rs^2 + (w_rated * Ls)^2) ; % 额定定子电流的d轴分量
flux_rated = Lm * Ids_rated; % 额定转子磁链
Kt = 1.5 * np * Lm^2 * Ids_rated / Lr; % 转矩系数(用在转速控制器)
Iqs_rated = Te_rated / Kt; % 额定定子电流的q轴分量
wsl_rated = Iqs_rated / (tau_r * Ids_rated); % 转差频率的额定值(用在电流环的前馈解耦)

% Load parameters (负载参数)
J = 0.0252; % 负载转动惯量
B = 0; % 阻尼系数

% ACR parameters (电流控制器参数)
wcc = 0.0707 * fsw * 2 * pi; % 电流控制器的带宽
Kpc = L * wcc; % 电流控制器的比例增益
Kic = R * wcc; % 电流控制器的积分增益
Kac = 1 / Kpc; % 电流控制器的前馈增益

% ASR parameters (速度控制器参数)
wcs = wcc / 5; % 速度控制器的带宽
Kps =  J * wcs / Kt; % 速度控制器的比例增益
Kis = Kps * wcs / 5; % 速度控制器的积分增益
Kas = 1 / Kps; % 速度控制器的前馈增益
disp('complete the parameters setting');