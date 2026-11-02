# Black-Scholes 与 Greeks：纸笔推导课程

这套资料把原作业的全部五页改写为八篇数学课，再扩展到 Put、Vega、Theta、Rho 和高阶 Greeks。你已经略懂 BS 与扩散模型，因此从“扩散怎样变成定价积分”开始，不重复完整的随机微积分课程。写法沿用 [复数 CDF 研究报告](../complex_cdf_series.md)：先定义、逐步推导、纸笔停顿、答案、数值对照。

## 怎样学

准备一个计算器和八组纸。通常每篇约 60-90 分钟，第 01 篇增补后约 90-120 分钟。先遮住“停顿答案”，再推导。第一轮重建公式；第二轮根据误差曲线解释算法。所有价格与 Greeks 默认按一单位标的计，合约乘数由调用方处理。

| 篇 | 学习资料 | 你要在纸上完成的结果 |
| --- | --- | --- |
| 01 | [P/Q 测度、扩散与 BS 价格](01_model_and_price.md) | 区分预测概率与定价权重，推出 Call 积分与 PDE |
| 02 | [Delta、Gamma 与曲率](02_delta_gamma.md) | 证明密度恒等式，求一阶和二阶导 |
| 03 | [有限差分与浮点误差](03_finite_differences.md) | 展开四种差分，平衡截断和舍入误差 |
| 04 | [Complex-step 与 45° Gamma](04_complex_step.md) | 修正原作业，推出二阶与四阶收敛 |
| 05 | [复数 CDF：Hermite 到 P3](05_cdf_and_order.md) | 解释 P1、P2、P3 为何给出不同 Gamma |
| 06 | [Vega、Theta、Rho 与股息敏感度](06_first_order_greeks.md) | 推出其公式、单位与 PDE 关系 |
| 07 | [Put 定价、平价与全部一阶 Greeks](07_put_and_parity.md) | 从到期收益推出 Put 与 Call 的关系 |
| 08 | [高阶 Greeks、多变量 Taylor 与 IV](08_higher_greeks.md) | 推出 Vanna、Vomma、Speed、Charm、Color 等 |

可打印版本位于 `output/pdf/black_scholes_math/`：每篇独立 PDF，另有合订本 `black_scholes_paper_workbook.pdf`。Markdown 是可维护的正文；PDF 由正文生成，附纸笔留白。重建：`python -X utf8 scripts/build_bs_math_workbook.py`。生成器需要 Matplotlib、ReportLab、Pillow、pypdf 和中文字体；当前 Windows 运行时会查找 Codex 的捆绑依赖。教学数值脚本另需 mpmath，以及已有的 NumPy/SciPy。

第 01 篇的换测度权重可在第二轮阅读，正文说明为什么不能把 $\Phi(d_2)$ 或 Delta 直接当成 P 下的事件发生概率。只刷新第 01 篇与合订本：`python -X utf8 scripts/build_bs_math_workbook.py --lessons 01`；其余分篇 PDF 须已存在。

## 统一符号与边界

到期日固定为 $t_e$，当前日历时间为 $t$，剩余期限 $T=t_e-t$，按年计。$S,K>0$、$\sigma,T>0$；$r,q$ 为固定、连续复利年率。$q$ 是连续股息率。计算某个偏导时，除该变量外其余输入固定。

$$
F=Se^{(r-q)T},\quad D=e^{-rT},\quad A=Se^{-qT},\quad B=Ke^{-rT},\quad \tau=\sigma\sqrt{T}.
$$

$$
m=\log(F/K),\quad d_1=\frac{m}{\tau}+\frac{\tau}{2},\quad d_2=d_1-\tau.
$$

$\phi$ 是标准正态密度，$\Phi$ 是标准正态 CDF。$\Theta=\partial V/\partial t=-V_T$。Vega 记为 $\nu=V_\sigma$（Vega 本身不是一个希腊字母名）；Rho 为 $\rho=V_r$，股息敏感度记为 $\psi=V_q$。这里 Vega、Rho 都按输入改变 1.00 计；“每增加一个百分点”须乘 0.01。

价格公式与 $d_1,d_2$ 默认采用货币市场账户为计价基准的 Q 测度。P 是现实预测的目标测度；第 01 篇另定义 $\mu$ 为含股息预期总收益率，P 下除息漂移为 $\mu-q$。Q 下到期价内概率与 P 下到期价内概率必须分开标记。

这些是欧式、连续股息 BS 的数学资料。$T=0$ 或 $\sigma=0$ 的折点不能直接套光滑导数公式。仓库对退化输入的 midpoint Delta、zero Gamma 是接口约定，见第 07 篇；它们不是折点存在经典导数的证明。

## 对原作业的修正

| 原作业位置 | 原文问题 | 本课程采用的修正 |
| --- | --- | --- |
| 第 1 页，密度计算 | 对数计算被描述为可以避免 underflow | 先在对数域合并比例因子可减少中间下溢；最终结果小于可表示范围仍会下溢 |
| 第 2 页，CS 实部二阶导 | 从展开移项到真导数的四阶导项符号有误 | 估计量等于真二阶导减去 $h^2 f^{(4)}/12$；反解真导数时应加回来 |
| 第 2 页，45° 对称 Taylor | 写成 $\omega^{2j}+\omega^{-2j}$ | $x-h\omega$ 的偶数幂仍为 $\omega^{2j}$；正确系数是 $2\omega^{2j}$ |
| 第 2 页，CDF 分派 | 一阶截断被描述为使价格解析延拓 | 有限阶实虚部分派一般不解析；P1 对 Delta 的局部二阶精度不能推广为所有 Gamma 的精度 |
| 第 2 页，45° 余项 | 原文给出较宽松的 $O(h^6)$ | 对完整解析且实轴取实值的函数，筛掉八阶项后是 $O(h^8)$；P3 截断价格不自动继承这一余项 |

第 04-05 篇用展开和反例逐条解释，不要求你接受一个没有证明的“修正结论”。当前生产 CDF 默认 P3，可选 `order=2..6`；本课程不改变生产代码，也不把尚未实现的 Put/其他 Greeks 写成已有 API。

## 原作业逐页覆盖

| 页 | 全部任务内容 | 教材位置 |
| --- | --- | --- |
| 1 | 稳定实数价格、解析 Delta/Gamma、前向 FD、Taylor 起点 | 01-04；密度数值问题见 02-03 |
| 2 | CS 一阶、实部二阶、45° 二阶、实/复 CDF 分派 | 04-05，含上表修正 |
| 3 | 类型分派价格、三种 CS Greeks、两个场景、相对步长网格 | 04-05；03 的实验流程 |
| 4 | CSV 精确列名、误差图/统计、稳定性解释、方法建议、代码清单 | 03 的实验交付清单与 05 的方法选择 |
| 5 | 两个 CSV 与报告提交清单 | 03 的实验交付清单 |

## 仓库入口与数值证据

| 现有相对路径 | 当前职责 |
| --- | --- |
| `Python_prototype/black_scholes.py` | 实数/复数 Call、Delta/Gamma、FD、CS、P2-P6 |
| `Python_prototype/bs_batch.py` | 校验输入后的 Call 批量接口 |
| `BS_cpp/black_scholes.h`、`BS_cpp/complex_step.h`、`BS_cpp/bs_batch.h` | 对应 C++ 接口；复数接口目前只把 S 暴露为复数 |
| `experiments/assignment/black_scholes_sweep.py` | 原作业的步长扫描与 CSV |
| `experiments/complex_cdf/study_series.py` | 独立 80 位 CDF 截断研究 |
| `experiments/bs_math/numerical_examples.py` | 本课程的教学参考计算、图和数列，含 Put/扩展 Greeks |
| `artifacts/bs-math/` | 教学脚本生成的数值表和曲线 |
| `scripts/build_bs_math_workbook.py` | Markdown 到可打印 PDF |

运行 `python -X utf8 experiments/bs_math/numerical_examples.py` 可刷新教学图与数值表。它是独立的教学脚本，并非新增生产 API 或全套生产测试。课程中的“精确/解析”指模型的闭式表达，不代表计算机输出无舍入误差。

本次 80 位公式对照、扫描统计与输出边界见 [教学数值对照记录](numerical_notes.md)。

## 出处与阅读范围

原作业为仓库根目录的 `Quant_Black_Scholes_Greeks_Assignment.pdf`，共五页。本套推导为重新编写；没有照抄整篇外部教材。网页访问日期为 2026-10-04。

- [NIST DLMF 7.2](https://dlmf.nist.gov/7.2)：erf/erfc 的定义与复变量性质；[7.10.1](https://dlmf.nist.gov/7.10#E1)：Hermite 与高阶导数。用于第 05 篇。
- [Martins、Sturdza、Alonso，2003，作者机构出版记录](https://mdolab.engin.umich.edu/bibliography/Martins2003a)：complex-step 论文出处，TOMS 29(3), 245-262。用于第 04 篇的方法背景；本教材的公式推导逐行给出。
- [Oxford Financial Derivatives 课程纲要](https://courses.maths.ox.ac.uk/course/info.php?id=6476)：扩散、风险中性、PDE、Greeks、股息的学习体系；只作体系定位，不替代推导。
- [Utah 金融数学第 6 章讲义](https://www.math.utah.edu/~zhu/5760.13f/chapter06note.pdf)：第 2.8、3.5 节的等价测度与 Girsanov。第 01 篇按本课程的连续股息约定推导 P/Q 差别。
- [Engstrom，2014，FEDS Note](https://www.federalreserve.gov/econresdata/notes/feds-notes/2014/forecasting-stock-market-crashes-is-hard-especially-future-ones-20140507.html)：期权隐含概率、风险厌恶与现实预测的区别。第 01 篇用作实证背景，不将特定样本的偏差方向推广为普遍定理。
- [MathWorks blsprice](https://www.mathworks.com/help/finance/blsprice.html)：欧式 Call/Put 与连续收益率输入约定；[blstheta](https://www.mathworks.com/help/finance/blstheta.html)：年度 Theta 的单位换算。用于 06-07 的接口单位对照。本文以偏导定义为准，不混用不同软件对“每百分点”的输出约定。

## 建议节奏

第一周 01-02：先分清 P/Q，再从风险中性终值积分走到 Delta/Gamma。第二周 03-05：把 Taylor、浮点误差、复数 CDF 连接起来。第三周 06-07：手算全部一阶 Greeks 与 Put。第四周 08：混合导数、多变量风险近似和 IV 条件数。每周挑一个数值例子，先预测符号和阶数，再看脚本结果。
