# 教学数值对照记录

日期：2026-10-04。此记录说明数学资料中数字与图片的来源；不是生产接口全参数域的准确性证书。

## 1. 独立价格求导

脚本：`experiments/bs_math/numerical_examples.py`。使用 mpmath 80 位十进制算术，完整 `erfc` CDF，独立 Call 与 Put 闭式价格。对输入作高精度数值求导，与教材解析 Greeks 对照。

四组参数为原作业场景 1、场景 2、带非零利率/股息的 ATM，以及带非零利率/股息的 OTM。共两种期权，每组包含 16 个导数及一个 PDE 残差，共 136 项。导数包括 Delta、Gamma、Vega、Theta、Rho、股息敏感度、Vanna、Vomma、Speed、Zomma、Charm、Color、Veta、Ultima、Dual Delta、Dual Gamma。

本次最大绝对差约为 `1.3494e-78`。这表示这几个样例中的公式与高精度价格偏导一致，不能把有限样例的一致性误写成整个定义域的证明。数学证明在每篇正文中给出。

可审查结果：`artifacts/bs-math/greeks_reference.csv`、`derivative_comparisons.csv`。前者给出 24 位显示值；后者保留高精度比较差。比较时未先将参考值四舍五入成显示值。

## 2. 原作业扫描

调用现有 Python 的解析、FD、CS 内核，P3 CDF，每个场景 24 个 `h/S`。两份 CSV 严格使用原作业表头。数据文件为 `scenario1.csv`、`scenario2.csv`；全网格统计为 `error_summary.csv`。本次 240 个数值估计均为有限数，但含错误为零的 Delta 对照点与不能分辨的 Gamma 点。

此处参考是生产解析 float64 值。扫描检查的是方法差异，不把解析 float64 当成无限精度真值。零误差点在 log 图上省略；统计保留其真实零值，不填造数值。最佳值只表示网格中的最佳值，不能声称为连续最优步长。

## 3. P1 的错误极限

独立高精度价格使用原作业 P1 CDF，CS-real Gamma 与真 Gamma 的比值随 h 缩小趋于 2；P2/P3 比值趋于 1。完整数列为 `p1_gamma_limit.csv`。P1 仅存在于教学脚本，生产接口仍只接受 2-6。

## 4. 图片、打印版和实现边界

六张图分别展示价格/Delta/Gamma、短期限 Gamma 宽度、原作业全部误差曲线、Greeks 期限形状、多变量 Taylor 近似、P/Q 对数收益密度。`taylor_shocks.csv` 保留 Taylor 图的原始数值。PDF 使用 Markdown 正文、中文嵌入字体、400 dpi 数学公式图片，以及纸笔停顿留白。

课程补充的 Put、其他 Greeks 与 IV 讨论是教学内容；生产库本次没有增加这些 API。没有运行完整生产测试或修改原有测试结果。生成的 PDF 沿用仓库现有忽略规则；Markdown、生成器、教学 CSV 和图可维护在仓库中。

原版八份分篇 PDF 和 37 页合订本已逐页检查。增补后第 01 篇为 8 页，合订本为 41 页，保留九个章节书签，有 27 个外部链接标注。已查看增补分篇全部八页与更新目录四页的渲染；本次其余分篇沿用已检查版本。合订本全部页面的文字和公式图片边界检查未发现越界，并对照分篇确认合并内容一致。PDF 中的数学式以高分辨率图片呈现，维护与复制公式应使用 Markdown 源文件。

## 5. P/Q 同一事件对照

第 01 篇事件固定为 `S_expiry > K`，参数为 `S=K=100,r=0.03,q=0.01,vol=0.20,T=1`。`mu` 表示 P 下含股息预期总收益率，分别假设为 8%、3%、1%；不是行情估计。用 mpmath 80 位算术计算 Q/P 的终值阈值、概率、Delta 与价格，保留 30 位显示值。

输出为 `measure_probabilities.csv` 与 `measure_distributions.png`。Q 概率保持 50%；三种假设的 P 概率依次约 59.8706%、50%、46.0172%，BS Call 价格均为 8.827321。图片只画 8% 这一组。概率与导数都使用完整实数 CDF，不依赖复数 CDF 截断阶数。

单独刷新此例：`python -X utf8 -c "from experiments.bs_math.numerical_examples import measure_example; measure_example()"`。本次增补没有重新运行第 1-3 节的计算，旧数值对照的边界与记录继续保留。
