close all;
clear;
clc;

% 逆变器参数
Vdc = 28;
fsw = 421;
Tsw = 1 / fsw;
pxyz = sqrt(3) * Tsw / Vdc;

% 电机参数
P_rated = 64; % 额定功率
V_rated = 24; % 额定电压（线电压有效值）
I_rated = 3.13; % 额定电流
f_rated = 200; % 额定频率 根据3000rpm和4p计算
w_rated = 2 * pi * f_rated; % 额定角速度
np = 4; % 极对数

V_max = sqrt(2) / sqrt(3) * V_rated; % 相电压最大值

rpm_syn = 60 * f_rated / np; % 同步转速
w_syn = rpm_syn * pi / 30; % 同步角速度
rpm_rated = rpm_syn; % 对于PMSM，转子速度等于同步速度

% PMSM参数
Rs = 0.89; % 定子电阻
Ld = 0.00062; % d轴电感
Lq = 0.00062; % q轴电感
Psi_f = 0.005926786; % 永磁磁链

% 转矩计算
Kt = (3/2) * np * Psi_f; % 转矩常数
Te_rated = P_rated / w_syn; % 额定转矩
Ids_rated = 0; % d轴电流设为零
Iqs_rated = Te_rated / Kt; % 计算q轴电流

% 负载参数
J = 0.000028; % 转动惯量
B = 0.00008; % 阻尼系数

% 电流控制器参数
wcc = 0.0707 * fsw * 2 * pi; % 电流控制器带宽
Kpc = 2769; % 比例增益
Kic = 2394; % 积分增益
Kac = 1 / Kpc; % 前馈增益

% 速度控制器参数
wcs = wcc / 5; % 速度控制器带宽
Kps = 10; % 比例增益
Kis = 700; % 积分增益
Kas = 0; % 前馈增益

disp('完成PMSM参数设置');
