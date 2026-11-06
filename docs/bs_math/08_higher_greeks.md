# 08 高阶 Greeks：Hessian、时间方向与 IV 条件数

目标：用微积分重建高阶 Greeks，把多个风险输入的变化统一成 Taylor；扩展到 IV 的数值问题。用时约 90-120 分钟。前置：02、06、07。

## 1. 多变量价格不是只有 Spot 曲率

先只考虑 $S,\sigma$ 两个输入。记 $s=\delta S$、$v=\delta\sigma$，价格变化到二阶为

$$
\delta V=\Delta s+\nu v+\frac12\Gamma s^2+V_{S\sigma}sv+\frac12 V_{\sigma\sigma}v^2+O((|s|+|v|)^3).
$$

这里各输入的单位已固定；严格比较“大小”时应先按各自尺度无量纲化。混合项前没有 $1/2$：Hessian 二次型里对称的两个交叉项各有 $1/2$，加起来正好是 1。

Vanna 定义为 $V_{S\sigma}=\Delta_\sigma=\nu_S$；Vomma/Volga 定义为 $V_{\sigma\sigma}=\nu_\sigma$。名称可能有不同软件约定，本教材始终以偏导式定义为准。

## 2. Vanna：先从 Delta 求 sigma 导数

第 06 篇已得 $(d_1)_\sigma=-d_2/\sigma$。因此

$$
\operatorname{Vanna}=e^{-qT}\phi(d_1)(d_1)_\sigma=-\frac{e^{-qT}\phi(d_1)d_2}{\sigma}.
$$

也可从 $\nu=A\phi(d_1)\sqrt T$ 对 $S$ 求导，用 $\phi'=-d_1\phi$ 得到

$$
\operatorname{Vanna}=\frac{\nu}{S}\left(1-\frac{d_1}{\tau}\right).
$$

两式相等因为 $\tau-d_1=-d_2$。在光滑区域混合偏导对称，提供第二条检查路线。

## 3. Vomma：对 Vega 再求一次 sigma 导数

Vega 中 $A\sqrt T$ 不随 $\sigma$ 变，只需求密度导数：

$$
\operatorname{Vomma}=\nu\frac{d_1d_2}{\sigma}.
$$

Vega 始终为正，Vomma 的符号却由 $d_1d_2$ 决定，所以价格对波动率不必全局凸。以第 07 篇的参数卡为例，$d_2=0$，Vanna 与 Vomma 都为零；这是特殊点，不代表一般值为零。

**纸笔停顿 A：**在零 carry ATM，$d_1=\tau/2,d_2=-\tau/2$，Vomma 是正还是负？

### 停顿 A 答案

$d_1d_2=-\tau^2/4$，故 $\operatorname{Vomma}=-\nu\sigma T/4<0$。波动率增加会提高价格，但增加带来的边际价格收益可以下降。

## 4. Speed 与 Zomma：Gamma 对输入怎样变化？

对 Gamma 先取对数：$\log\Gamma=-qT+\log\phi(d_1)-\log S-\log\sigma-\log T/2$。对 $S$ 求导，得到 Speed：

$$
\operatorname{Speed}=V_{SSS}=-\frac{\Gamma}{S}\left(1+\frac{d_1}{\tau}\right).
$$

对 $\sigma$ 求导，则得到 Zomma：

$$
\operatorname{Zomma}=V_{SS\sigma}=\frac{\Gamma}{\sigma}(d_1d_2-1).
$$

Speed 是第 03 篇前向 Gamma 领先误差 $h f^{(3)}$ 的系数。因此一个“高阶 Greek”同时也是数值分析里决定截断误差的量。Delta 的 CS 领先误差则是 $-h^2\operatorname{Speed}/6$。

## 5. 时间相关 Greeks：统一 calendar-time 符号

记

$$
g=(d_1)_T=\frac{r-q}{\tau}-\frac{d_2}{2T}.
$$

Charm 定义为 Delta 对日历时间的偏导 $\partial_t\Delta=-\Delta_T$。对 Call/Put 分别有

$$
\operatorname{Charm}_C=q\Delta_C-e^{-qT}\phi(d_1)g.
$$

$$
\operatorname{Charm}_P=q\Delta_P-e^{-qT}\phi(d_1)g.
$$

两者之差为 $q e^{-qT}$。这个差也可直接对 Delta 平价求日历时间导数得到。

Color 定义为 $\partial_t\Gamma=-\Gamma_T$。对数求导给出

$$
\operatorname{Color}=\Gamma\left(q+\frac1{2T}+d_1g\right).
$$

Veta 定义为 $\partial_t\nu=-\nu_T$，同理

$$
\operatorname{Veta}=\nu\left(q-\frac1{2T}+d_1g\right).
$$

这些时间导数固定当前 Spot、其他模型输入不变，不等于沿一条随机市场路径观察到的 Delta/Gamma/Vega 总变化。有些资料采用剩余期限导数定义 Charm/Color/Veta，符号会与这里相反。

**纸笔停顿 B：**从 $\log\Gamma$ 推 Color，再从 $\log\nu$ 推 Veta。为什么一个出现 $+1/(2T)$，另一个是负号？

### 停顿 B 答案

$\partial_T\log\Gamma=-q-d_1g-1/(2T)$，$\partial_T\log\nu=-q-d_1g+1/(2T)$。日历时间导数取负号；根源是 Gamma 含 $T^{-1/2}$，Vega 含 $T^{1/2}$。

## 6. Ultima、Dual Greeks 与弹性

设 $p=d_1d_2$。由 $d_{1,\sigma}=-d_2/\sigma$、$d_{2,\sigma}=-d_1/\sigma$，得到 $p_\sigma=-(d_1^2+d_2^2)/\sigma$。再对 Vomma 求导：

$$
\operatorname{Ultima}=V_{\sigma\sigma\sigma}=\frac{\nu}{\sigma^2}\left(p(p-1)-d_1^2-d_2^2\right).
$$

对执行价 $K$ 求导则是 Dual Delta，不应与 Spot Delta 混同：

$$
C_K=-e^{-rT}\Phi(d_2),\quad P_K=e^{-rT}\Phi(-d_2).
$$

$$
C_{KK}=P_{KK}=\frac{e^{-rT}\phi(d_2)}{K\sigma\sqrt T}.
$$

最后，价格弹性 Lambda 为 $\lambda=S\Delta/V$，表示局部百分比价格变化率。若 $V$ 很小，这个比率会很大且病态；$V=0$ 时不应直接相除。

所有含 $\sigma$ 或至少二次 $S$ 导数的这些高阶 Greeks，在同参数 Call/Put 下相同，因为平价差 $A-B$ 与波动率无关且对 Spot 线性。Charm 是一阶 Spot 加时间导数，因而保留 Call/Put 差异。

## 7. 怎样数值核验混合导数？

一个不用新增生产接口的中心混合差分是

$$
V_{S\sigma}\approx\frac{V(S+h,\sigma+k)-V(S+h,\sigma-k)-V(S-h,\sigma+k)+V(S-h,\sigma-k)}{4hk}.
$$

局部截断误差为 $O(h^2+k^2)$，舍入误差可被 $1/(hk)$ 放大。也可对解析 Vega 的 $S$ 方向差分、对解析 Delta 的 $\sigma$ 方向差分来交叉检查；两者更省价格调用，但若解析式共同错误仍可能相互支持。

不要简单嵌套同一个普通虚数 $i$ 的 CS 来求混合导数：内层的“取虚部”通常已经破坏外层所需的解析结构。高阶 CS、multicomplex、自动微分可作为后续路线，需要明确表示法和理论。

![Delta-Gamma-Vega 近似与混合项](../../artifacts/bs-math/taylor_risk_approximation.png)

图在一个非特殊参数点同时改变 $S,\sigma$，比较价格的精确重新计算、线性近似、二阶近似。它是模型内部的局部误差，不是市场预测表现。

## 8. 把 Greeks 接到模型风险与曲面输入

若外部系统输入的波动率本身是 $\sigma(S)$，价格的总 Spot 导数就不是固定波动率 Delta：

$$
\frac{dV}{dS}=\Delta+\nu\sigma_S.
$$

$$
\frac{d^2V}{dS^2}=\Gamma+2\operatorname{Vanna}\sigma_S+\operatorname{Vomma}\sigma_S^2+\nu\sigma_{SS}.
$$

这解释了为什么系统应同时说明“模型输入固定的 Greeks”与“按某个波动率曲面规则移动的风险”。本课程前面的公式均为固定其余 BS 输入的偏导。

## 9. IV：Vega 是反问题的局部条件数

给定合规市场价格 $V_m$，求 $V(\sigma)=V_m$。在 $S,K,T>0$、有限正 $\sigma$ 下，Vega 大于 0，所以内部可达价格对应唯一正 IV。边界价格可能只对应零波动率或无穷波动率的极限；越界报价没有此模型下的解。

微分隐函数得

$$
\delta\sigma\approx\frac{\delta V_m}{\nu},\quad \frac{d\sigma}{dV_m}=\frac{1}{\nu}.
$$

深 OTM 或极短期限时 Vega 可很小，微小报价误差就对应很大 IV 不确定性。Newton 更新为

$$
\sigma_{n+1}=\sigma_n-\frac{V(\sigma_n)-V_m}{\nu(\sigma_n)}.
$$

大步可能产生负波动率，且 Vomma 会改变价格对波动率的曲率。可靠实现应先检查无套利边界，建立 bracket，牛顿越界或 Vega 太小时回退二分。价格残差小不自动保证 IV 准：还应考虑残差除以 Vega、bracket 宽度，以及 bid/ask 对应的 IV 区间。

**纸笔停顿 C：**若价格误差为 0.01，而 Vega 分别为 40 与 0.1，估计 IV 误差。

### 停顿 C 答案

分别为 0.00025（0.025 个百分点）与 0.1（10 个百分点）。绝对价格精度相同，反求波动率的精度却不同。本仓库当前未发布 IV 求解器。

## 10. 最后一张验收卡

遮住正文，从 $\phi'=-x\phi$、$d_{1,\sigma}=-d_2/\sigma$ 和 Gamma/Vega 的对数表达式，重建 Vanna、Vomma、Speed、Zomma、Charm、Color、Veta、Ultima。再说明它们怎样进入价格 Taylor、差分截断误差、IV Newton 分母，以及系统的风险输入链。

本套课程到这里完成原作业的全部数学任务，并扩展了 Put 与常用高阶 Greeks；生产接口是否实现仍以课程首页的状态表为准。
