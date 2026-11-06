# 07 Put：从收益平价到价格、Greeks 与边界

目标：把 Call 的推导扩展到 Put，而不是重新背另一套公式。用时约 75 分钟。前置：01、02、06。

## 1. 到期收益恒等式

对任意终值 $s$，分 $s\ge K$ 与 $s<K$ 两种情况，可检查

$$
(s-K)^+-(K-s)^+=s-K.
$$

在同一标的、执行价、期限、连续股息与欧式行权约定下，取风险中性期望并折现：

$$
C-P=e^{-rT}(F-K)=A-B.
$$

这是 put-call parity。它来自终值收益完全一致，以及相同收益应有相同现值的无套利原则。Put 是

$$
P=C-A+B=B\Phi(-d_2)-A\Phi(-d_1).
$$

后一步使用 $1-\Phi(x)=\Phi(-x)$。也能直接积分 $Z<-d_2$：现金部分为 $B\Phi(-d_2)$，资产加权部分通过配平方变为 $A\Phi(-d_1)$。

**纸笔停顿 A：**分两种终值情况验证收益平价，再从 Call 闭式式子手动推 Put。

### 停顿 A 答案

$s\ge K$ 时是 $(s-K)-0=s-K$；$s<K$ 时是 $0-(K-s)=s-K$。由 $C-A+B$ 展开并将 $1-\Phi$ 改为负参数 CDF，即得 Put。

## 2. 对平价求导：Delta/Gamma/Vega

对 $S$ 求导：$\Delta_C-\Delta_P=e^{-qT}$，所以

$$
\Delta_P=e^{-qT}(\Phi(d_1)-1)=-e^{-qT}\Phi(-d_1).
$$

平价右边是 $S$ 的线性函数，第二次求导为 0：

$$
\Gamma_P=\Gamma_C=\frac{e^{-qT}\phi(d_1)}{S\sigma\sqrt T}.
$$

右边不依赖 $\sigma$，所以 $\nu_P=\nu_C$。Call 与 Put 斜率不同，但在同一参数下曲率和波动率敏感度相同。

## 3. Rho、股息敏感度、Theta

固定 $S,q,\sigma,T$ 求 $r$ 导数：

$$
\rho_C-\rho_P=TB,\quad\rho_P=-TB\Phi(-d_2).
$$

固定 $S,r,\sigma,T$ 求 $q$ 导数：

$$
\psi_C-\psi_P=-TA,\quad\psi_P=TA\Phi(-d_1).
$$

先对 $T$ 求导再变号：$-\partial_T(A-B)=qA-rB$，所以

$$
\Theta_C-\Theta_P=qA-rB.
$$

$$
\Theta_P=-\frac{A\phi(d_1)\sigma}{2\sqrt T}+rB\Phi(-d_2)-qA\Phi(-d_1).
$$

**纸笔停顿 B：**在非零 $r,q$ 下从平价推 Theta 差。常见错误是把 $\Theta=-V_T$ 的负号漏掉。

### 停顿 B 答案

$\partial_T(A-B)=-qA+rB$，变号后是 $qA-rB$。因此 $\Theta_P=\Theta_C-qA+rB$，代入 Call Theta 得上述式子。

## 4. 一张共同的数值卡

取 $S=K=100,r=0.03,q=0.01,\sigma=0.2,T=1$。此时 $F=100e^{0.02}$、$d_1=0.2,d_2=0$，$A\approx99.00498337$、$B\approx97.04455335$。

| 项目 | Call | Put |
| --- | --- | --- |
| Price | 约 8.82732123 | 约 6.86689121 |
| Delta | 约 0.57349598 | 约 -0.41655385 |
| Gamma | 约 0.01935759 | 相同 |
| Vega，按 1.00 | 约 38.71517542 | 相同 |
| Rho，按 1.00 | 约 48.52227668 | 约 -48.52227668 |
| 股息敏感度，按 1.00 | 约 -57.34959790 | 约 41.65538547 |
| Theta，每年 | 约 -4.75368986 | 约 -2.83240310 |

价格差约 1.96043002，Delta 差约 0.99004983。表的值四舍五入显示；可复现的高精度值见 `artifacts/bs-math/greeks_reference.csv`。不要用显示后被截断的数值要求机器级残差。

## 5. 两条 Put 计算路线与相消

“先算 Call 再做 $P=C-A+B$”适合构造数学一致性对照，但若 Put 非常小，可能减去两个大的量而丢失相对精度。直接式使用 $\Phi(-d_1),\Phi(-d_2)$，能避免某些平价变换相消；不过两项本身也可能接近，极端参数仍需专门稳定算法。

类似地，不要把极端尾部的 $\Phi(-d)$ 总是写成 $1-\Phi(d)$。应使用直接负参数 CDF/erfc，或相应尾部原语。

数值 Greeks 既可对 Put 价格独立 bump，也可从 Call 与平价推得。两条路线交叉比较时，不应让它们完全共享同一段求导代码，否则不能充分发现错误。

## 6. 无套利边界与退化极限

对非负终值，$0\le(s-K)^+\le s$、$0\le(K-s)^+\le K$。再用凸性得到下界：

$$
\max(A-B,0)\le C\le A,\quad\max(B-A,0)\le P\le B.
$$

$\sigma\to0$ 时，终值变成确定的远期 $F$；价格极限为折现内在价值 $C=(A-B)^+$、$P=(B-A)^+$。$T=0$ 时价格等于即期到期收益。

在折点 $A=B$，Delta 的经典导数不存在。零波动率下的 Gamma 在折点不能被一个普通有限数表示。工程接口可以约定 Delta 取 midpoint、Gamma 取 0，但那是返回规则。甚至 Vega 的边界也需独立分析：零 carry ATM、固定正期限时，$\sigma\downarrow0$ 的右侧 Vega 极限为 $S\sqrt T/\sqrt{2\pi}$，不为 0。

**纸笔停顿 C：**固定零 carry ATM，证明 Call 价格在小 $\sigma$ 下约为 $S\sigma\sqrt T/\sqrt{2\pi}$；解释为什么 blanket “sigma=0 所有 Greeks 都是 0”不成立。

### 停顿 C 答案

$C=S(2\Phi(\sigma\sqrt T/2)-1)$。用 $\Phi(x)=1/2+\phi(0)x+O(x^3)$，得到 $C=S\phi(0)\sigma\sqrt T+O(\sigma^3)$，右侧 Vega 非零。边界导数须分变量、方向与是否位于折点讨论。

## 7. 把平价变成检查体系

价格平价、Delta 差、Gamma 相等、Vega 相等、Rho 差和 Theta 差可形成一组残差。再检查两个价格的 PDE 残差和无套利区间，以及独立高精度价格的数值偏导。

Put 价格仍由本篇的教学脚本提供。Python 和 C++ 的 [Greeks 接口](../bs_greeks.md) 可计算欧式 Put 的 Vega、Theta、Rho。欧式、连续收益率软件约定可对照 [MathWorks blsprice](https://www.mathworks.com/help/finance/blsprice.html)。美式行权与离散现金股息需要各自的定价模型。

## 8. Put 的四张独立冲击卡

回到第 4 节的非零利率/股息参数。保持其他输入固定，分别研究四个独立场景；不要把它们未经混合导数修正就合并成一张联合冲击。

**纸笔停顿 D：**预测 Put 的价格变化与符号：标的下跌 1（保留二阶 Spot 项）；波动率上升一个百分点（只留 Vega）；利率上升 10 个基点（只留 Rho）；股息率上升 10 个基点（只留股息敏感度）。先写出每个输入的实际小数变化，再计算。

### 停顿 D 答案

Spot 变化为 $s=-1$，$\delta P\approx\Delta_Ps+\Gamma s^2/2\approx0.42623265$。波动率变化为 0.01，价格线性变化约 0.38715175。10 个基点为 0.001，因此利率冲击约 -0.04852228，股息冲击约 +0.04165539。

计算完再用教学 Put 价格逐项重新定价，观察近似与精确价格变化的差。第一张卡已经保留 Gamma；其他卡只保留一阶敏感度，因此不要求有限冲击下完全相等。能先判断四个符号，再说明被省略的项，是本篇的纸笔验收。
