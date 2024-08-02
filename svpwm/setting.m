% inverter parameters
% 直流电源电压
Vdc = 500;
% 开关频率
fsw = 10e3;
% 开关周期
Tsw = 1 / fsw;
pxyz = sqrt(3) *Tsw / Vdc;
