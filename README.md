PMSM的dq轴电压方程

![image](https://github.com/user-attachments/assets/5d4791bb-362d-44a4-8fc2-64acc3b5e093)

![image](https://github.com/user-attachments/assets/4487af95-a6c0-4e75-abe5-679e2fd3e7b8)

最终传递函数：

$$
G_(s)=\frac{k_ps+k_i}{L_ds^2+(k_p+R)s+k_i}
$$

![image](https://github.com/user-attachments/assets/4d8e4301-8358-4402-a557-e1f99b62e3eb)

**带宽**：

带宽就是一个用来评价系统动态性能的参数。

控制系统的性能一般从两个方面分析，也就是稳态有无静差以及动态响应的快慢。评价一个系统动
态响应的快慢最简单的方式是给一个阶跃输入，看系统的调节时间。但是这种方式虽然直观但不方
便使用。于是，大佬们就从时域转向频域来考虑这个问题。比如都给100Hz变化的正弦信号，看两
个系统的失真情况。这个失真一般就是从幅值衰减和相位滞后两个角度考虑，因此，一般就规定
3dB衰减或者-45度相移作为带宽的评价标准。给一个变化速度超过带宽的输入信号，
系统就无法完
美的跟随，而存在一定的失真;但变化速度小于带宽的话，则系统可以基本跟随上。

至于电流环带宽100Hz，就可以理解为，如果输入电流指令信号是90Hz的正弦信号，电流环可以完
全跟随上。

但是在实际电机的电流控制中，由于采用了坐标变换，dq轴电流在稳态情况下都是直流给定，而只
在加减速或是突加减负载的时候，速度环会等效输出一个阶跃指令，**那么这里电流环的带宽就是电
流环能够以多快的速度去响应这个指令**。

![image](https://github.com/user-attachments/assets/8ac80eab-6208-4102-bf52-1ed6174afcdc)

然后

$$
k_p = L_d\omega_c 
$$
$$
k_i = R\omega_c
$$

ST指定PI参数：

![image](https://github.com/user-attachments/assets/d7f9c276-2ec3-4f4d-a80b-0812be50d5a3)

![image](https://github.com/user-attachments/assets/2e41a19e-2e8e-4cde-a493-0ca5125fa2c7)

电机运动方程：

![image](https://github.com/user-attachments/assets/5630f20c-5c03-4f0a-a8a1-79c811bf6a9e)


