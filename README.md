README

THIS IS A READ ME OF AtmosRefraction

##### **大气折射作业要求-2018**

一、功能：大气折射引起的太阳系外天体本征方向到观测方向的计算

二、要求：

1. ##### 这里只需考虑大气折射效应。

2. 这里计算蒙气差和表示天体方向的公式参照<天体测量学导论>教科书第55页和59页的内容。
   这里假设我们进行光学观测的波段靠近红端（中心波长约0.7微米），而不是通常55页公式所对应的0.5微米。

3. 迭代是否可退出的标准为：高度或天顶距的变化小于3.6毫角秒。

4. 输出结果的单位为度，有效位数保留到毫角秒（即精确到1毫角秒）。

5. 编程语言：C或Fortran或Python

三、已知参数：

1. 天体在台站处地平坐标系中的本征方向的地平坐标：方位A=30度 高度H=50度 （方位：南=0度，东=90度 右手系）
2. 台站处的气温和气压：T=20摄氏度， P=980百帕 
3. 观测波段的中心波长约0.7微米

四、输出量：此天体观测方向的方位和高度角（单位：角度，精度到1毫角秒）