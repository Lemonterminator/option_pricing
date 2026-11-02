# 06 Vega、Theta、Rho：对哪个输入求偏导？

目标：复用同一个密度恒等式推导其他一阶 Greeks，再把它们接到数值求导与 BS PDE。用时约 90 分钟。前置：02；数值部分前置 03-05。

## 1. 先做一张偏导表

记 $m=\log(S/K)+(r-q)T$、$\tau=\sigma\sqrt T$、$d_1=m/\tau+\tau/2$。一般变量 $x$ 下，商法则给出

$$
(d_1)_x=\frac{m_x}{\tau}+\left(\frac12-\frac{m}{\tau^2}\right)\tau_x,\quad(d_2)_x=(d_1)_x-\tau_x.
$$

| 变量 x | m 的导数 | tau 的导数 | d1 的导数 |
| --- | --- | --- | --- |
| S | $1/S$ | $0$ | $1/(S\tau)$ |
| sigma | $0$ | $\sqrt T$ | $-d_2/\sigma$ |
| r | $T$ | $0$ | $\sqrt T/\sigma$ |
| q | $-T$ | $0$ | $-\sqrt T/\sigma$ |
| T | $r-q$ | $\tau/(2T)$ | $(r-q)/\tau-d_2/(2T)$ |

所有参数独立；求 Rho 时固定 $q$，不能同时令 $q=r$ 跟着一起变。若系统中输入的是远期 $F$ 而不是现货 $S$，所谓“Rho”还必须说明固定哪一个。

**纸笔停顿 A：**自己推导 $d_{1,\sigma}$ 与 $d_{1,T}$，不要直接背表。

### 停顿 A 答案

$d_{1,\sigma}=\sqrt T(1/2-m/\tau^2)=\sqrt T-d_1/\sigma=-d_2/\sigma$。期限导数为 $(r-q)/\tau+\tau/(4T)-m/(2T\tau)=(r-q)/\tau-d_2/(2T)$。

## 2. Vega：CDF 导数项这次没有全消掉

对 $\sigma$ 求导时 $A,B$ 不变。由密度恒等式令 $H=A\phi(d_1)=B\phi(d_2)$：

$$
\nu=C_\sigma=H((d_1)_\sigma-(d_2)_\sigma)=H\tau_\sigma=A\phi(d_1)\sqrt T.
$$

因为 $S,\sigma,T>0$ 且 $\phi>0$，Vega 为正。可从 Gamma 公式消去同一个密度因子：

$$
\nu=\sigma T S^2\Gamma.
$$

这把“对扩散强度求导”和“对状态求两次导数”连接起来，适用于此处的常参数欧式 BS Call/Put，不是任意定价模型的普遍恒等式。

**纸笔停顿 B：**场景 1 的 Vega 为多少？波动率从 0.20 到 0.21 时，线性价格变化是多少？

### 停顿 B 答案

$\nu\approx39.69525475$，一个百分点是 $\delta\sigma=0.01$，价格线性变化约 0.39695255。若误乘 1，就把变化夸大 100 倍。场景 2 的 Vega 约为 2.08815926；Gamma 较大并不表示 Vega 也较大，因为它们相差 $\sigma T S^2$。

## 3. Rho 与股息敏感度

对 $r$ 求导，$A_r=0$，$B_r=-TB$，且 $d_{1,r}=d_{2,r}$。CDF 导数项消去，留下

$$
\rho_C=C_r=TB\Phi(d_2).
$$

对 $q$ 求导，$A_q=-TA$，$B_q=0$，同理

$$
\psi_C=C_q=-TA\Phi(d_1)=-TS\Delta_C.
$$

这两个导数取决于“其余变量固定”的约定。现实中期限结构、股息和远期可能一起变动；真实组合冲击需要把输入变化链式传播。

## 4. Theta：先求剩余期限，再变号

$A_T=-qA$、$B_T=-rB$、$d_{1,T}-d_{2,T}=\tau/(2T)$。所以

$$
C_T=-qA\Phi(d_1)+rB\Phi(d_2)+H\frac{\tau}{2T}.
$$

日历时间向前，剩余期限缩短，故 $\Theta=-C_T$：

$$
\Theta_C=-\frac{A\phi(d_1)\sigma}{2\sqrt T}-rB\Phi(d_2)+qA\Phi(d_1).
$$

第一项是扩散/曲率贡献，第二项与现金贴现有关，第三项与股息有关。Theta 并非所有参数下都为负；例如股息项可能改变符号，不能把“通常时间损耗”当作符号定理。

将它和 Delta/Gamma 代回 PDE，应满足

$$
\Theta+\frac12\sigma^2S^2\Gamma+(r-q)S\Delta-rV=0.
$$

**纸笔停顿 C：**场景 1 中 $r=q=0$，用 PDE 而不是 Theta 闭式式子求 Theta。

### 停顿 C 答案

$\Theta=-\sigma^2S^2\Gamma/2\approx-3.969525475$（每年）。每日线性近似为 $\Theta/365\approx-0.01087541$，不是把原公式直接改成每天到期；有限一天的价格差还包括高阶期限效应。

## 5. 数值方法如何统一处理这些 Greeks？

将价格视为六个输入的函数 $V(S,K,r,q,\sigma,T)$。对选定输入 $x$，定义中心差分与 CS：

$$
V_x\approx\frac{V(x+h_x)-V(x-h_x)}{2h_x},\quad V_x\approx\frac{\operatorname{Im}V(x+ih_x)}{h_x}.
$$

Vega bump $\sigma$；Rho bump $r$；股息敏感度 bump $q$；Theta 则在 bump $T$ 得到的导数前加负号。每种输入的量纲不同，不能统一用同一个绝对步长。

对于 $\sigma,T$ 的中心差分，必须让 $x-h_x>0$。在边界附近可选单边公式，但要重新分析阶数。CS 在数学上可避开实轴跨越，却仍依赖内部复数 sqrt/log 与分派，而不是给 double 参数强塞复数就一定能工作。

普通一阶偏导用完整解析价格时是 $O(h_x^2)$。有限截断 CDF 经新参数的复合链后，应逐方向验证；对 Spot Delta 验证过，不表示 Vega/Theta 自动也验证过。

## 6. 建立三类一致性检查

解析式对照数值偏导，是第一类。PDE 残差 $R_{PDE}$，是第二类。Vega-Gamma 关系残差，是第三类：

$$
R_\nu=\nu-\sigma TS^2\Gamma.
$$

这些残差能定位单位、符号和公式错误，但共享代码错误可能同时让多个残差变小。至少搭配独立高精度价格求导和非零 $r,q$ 场景；零利率零股息会隐藏很多错误。

![Greeks 随期限变化](../../artifacts/bs-math/greeks_term_structure.png)

图取零 carry、ATM、$\sigma=0.2$，展示 Gamma、Vega 和年度 Theta 随期限改变。不要比较纵轴数值大小来判断哪个风险“更大”，因为单位不同。

## 7. 公式状态与下一步

本篇推导的 Vega、Theta、Rho、股息敏感度用于教学参考脚本；当前生产批量接口只返回 Call price/Delta/Gamma。扩展生产时先明确单位、边界和输出字段，再逐层暴露接口。

年度 Theta 的软件单位换算可对照 [MathWorks 文档](https://www.mathworks.com/help/finance/blstheta.html)。本文始终按 $-V_T$ 定义；输出“每百分点”的系统必须显式作 0.01 缩放。下一篇利用收益平价把整套推导传给 Put。
