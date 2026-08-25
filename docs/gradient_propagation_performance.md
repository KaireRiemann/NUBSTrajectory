# NUBS 梯度传导性能问题与优化总结

## 问题定义

NUBS 与 MINCO 使用相同的外部决策变量：中间航点 `P` 和分段时长
`T`。NUBS 的控制点由

\[
A(T)C=b(P)
\]

恢复，因此能量梯度需要同时包含直接项和隐式控制点变化：

\[
\frac{dE}{dT_k}=
\left.\frac{\partial E}{\partial T_k}\right|_C-
\lambda^\mathsf{T}\frac{\partial A}{\partial T_k}C,
\qquad A^\mathsf{T}\lambda=\frac{\partial E}{\partial C}.
\]

早期实现把每个物理 B-spline span 的所有局部时长方向都携带在
`LocalJet` 中。虽然一个 span 最多只依赖 `2p-1` 个时长，B-spline
基函数、导数递推、结点差分和能量累加的每一次基础运算都要更新整条
jet。构造项中的 `M-1` 个 waypoint 行也重复执行同一模式。故而，
**控制点系统维度更小并不自动意味着梯度更快**：原始瓶颈是每个局部
方向上的微分工作，而不是线性系统的行列数。

## 当前实现

固定阶 `NUBSTrajectoryT<Dim, s>` 的生产路径已采用如下方案：

1. 对控制点重复应用导数控制点递推

   \[
   D_i^{(r)}=
   \frac{p-r+1}{u_{i+p+1}-u_{i+r}}
   \left(D_{i+1}^{(r-1)}-D_i^{(r-1)}\right),
   \]

   以获得次数为 `s-1` 的曲线 `p^(s)`；其平方范数由 `s` 点 Gauss
   积分精确计算。
2. 显式反传上述差分递推及低阶基函数的 Algorithm A2.2 递推，直接得到
   局部结点/时长导数；不再让前向 jet 穿过每个基础运算。
3. 先通过一次 `A^T` 求解获得构造伴随变量，再对每个内部 waypoint 行
   作标量基函数反传。端点导数行数量为 `O(s)`，仍保留小型 `LocalJet`
   路径。

该实现在 [NUBSTrajectory.hpp](../include/NUBSTrajectory.hpp) 中由
`getEnergyPartialGradDerivativeControl` 和
`propagateEnergyGradDerivativeControl` 提供。

## 测量结果

以下是 Release、严格数值验证、每个配置 500 次重复调用的最终结果。
`direct` 是内部控制变量和时长的直接梯度，`full` 是传播到共同的航点/
时长决策变量后的完整梯度。单位为微秒。

| 阶数 s | 段数 M | NUBS direct | MINCO direct | NUBS full | MINCO full | NUBS/MINCO full |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 64 | 4.230 | 0.564 | 7.513 | 2.671 | 2.81x |
| 3 | 64 | 10.712 | 1.070 | 19.951 | 4.618 | 4.32x |
| 4 | 64 | 22.795 | 1.573 | 41.802 | 7.324 | 5.71x |

对于最具代表性的 `s=4, M=64`，旧标量 Local-AD 基线的直接/完整梯度
分别为 1026.9 / 1125.7 微秒；当前实现为 22.8 / 41.8 微秒，分别降低约
45 倍和 27 倍。与 MINCO 的完整梯度差仍为约 5.7 倍，因此论文中应表述为
“显著消除了 NUBS 自身的 AD 开销”，而不能表述为“已与 MINCO 同速”。

最终基准中的最大 NUBS--MINCO 误差（`s=4, M=64`）为：航点梯度
`1.714e-08`，时长梯度 `1.686e-07`。完整测试套件 26/26 通过，包含
局部解析梯度、稠密解析梯度和有限差分交叉验证。

`build-refactor/validation/gradient_kernel_comparison.svg` 展示了旧标量
Local-AD、融合 LocalJet、精确 LocalJet 和最终标量反传内核的演进。运行时结果以
`build-refactor/validation/gradient_propagation_benchmark.csv` 为准。

## 对论文表述的影响

- 较小的 NUBS 线性系统仍是正向构造上的优势，但不是完整梯度速度优势的
  充分条件。
- B-spline 的实际工程价值应表述为局部支撑、非均匀时长的自然表达、较小
  的前向带状系统及可局部化的微分；本实现没有把凸包性质作为碰撞安全
  证书，不能以此宣称直接安全保证。
- 若要回应端到端规划收益，必须在相同目标函数、相同决策变量和相同约束
  下报告轨迹构造、目标/导数/各惩罚项、梯度回传和线搜索的分项时间；本
  文档的微基准不能替代该规划级实验。

## 复现

```bash
bash scripts/run_validation.sh build-refactor 500
```

该命令会以严格数值模式构建项目、运行全部 26 项测试，并导出梯度、强非
均匀时长和构造基准结果及 SVG 图。
