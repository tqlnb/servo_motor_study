% Inveter parameters (逆变器参数)
Vdc = 200;
fsw = 20e3; % 开关频率
Tsw = 1 / fsw;


% DC motor parameters
power_rated = 3336; % 额定功率
Va_rated = 140; % 额定电压
Ia_rated = 25; % 额定电流
rpm_rated = 3000; % 额定转速rpm
wm_rated = rpm_rated * 2 * pi / 60; % 机械转速
Te_rated = power_rated / wm_rated;    % 额定转矩
Kt = Te_rated / Ia_rated; % 由电流得到转矩的系数


La = 1.7e-3; % 电枢电感？？
Ra = 0.26; % 电枢电阻


% Load Parameters （负载参数）
TL = 1; % 负载转矩
J = 0.00252; % 转动惯量
B = 0; % 阻尼系数（为0 不考虑）


% ACR parameters （电流控制器参数）
wcc = 0.0707 * fsw * 2 *pi;  % 截止频率
Kpc = La * wcc;
Kic = Ra * wcc;
Kac = 1 / Kpc;


% ASR parameters （速度控制器参数）
wcs = wcc / 5; % 截止频率
Kps = J * wcs / Kt;
Kis = Kps * wcs / 5;
Kas = 1 / Kps;