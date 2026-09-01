# CRESTU-1F Phase 3.2B 预声明响应收敛门槛

```ini
GATE_FREEZE_TIME = BEFORE_ANY_PHASE3_2B_RESPONSE_SOLVE
TOPOLOGY_MODE = QUALITY_CONTROLLED_V2_EXPLICIT
QUALITY_CONTROLLED_V2_DEFAULT = OFF
RAW_RCOND_MINIMUM = 1.0e-12
SCALED_RCOND_MINIMUM = 1.0e-12
RELATIVE_RESIDUAL_MAXIMUM = 1.0e-10
MAXIMUM_ABS_OUTER_NTHETA = 1.0e-10
MAXIMUM_OUTER_COLUMN_RESIDUAL = 1.0e-14
MAX_CONCURRENT_L3_BEM_SOLVES = 1
```

本文是 Phase 3.2B 的冻结判据。判据依据为：Phase 3.2A 已接受的几何质量与 1.5 rad/s bridge、Phase 3.1 已冻结的代数门槛、三等级工程网格收敛的常用相邻级比较，以及共振峰离散跟踪所需的频率分辨率。本文生成后不得依据 Phase 3.2B 新响应结果修改；如判据本身存在缺陷，只能将阶段标为 `UNCERTAIN` 或另开有版本、有理由、重新运行全部受影响 case 的后续阶段，不能覆盖本文件。

## 1. 前置、几何与代数门槛

每个 case 必须先在唯一目录内重新生成网格，并在组装矩阵前与 `Phase3_2B_V5_Frozen_Mesh_Manifest.csv` 核对：频率、mesh level、显式 topology mode、四组件 panel count、总未知数、四组件 hash、combined geometry hash、top/bottom outer extent、水深、production-source hash、configured/effective control hash，以及可在组装前得到的 operator-state hash。任何一项不匹配即 `V5_MESH_REPRODUCIBILITY=FAIL`，该 case 不得进入矩阵组装。

每个 clean solve 必须同时满足：

- raw rcond 与 scaled rcond 均不小于 `1.0e-12`；
- radiation/diffraction 的 linear residual 与 algebraic Frobenius relative residual 均不大于 `1.0e-10`；
- 继承 bridge 的更严格 operator audit：maximum absolute outer `n_theta` 不大于 `1.0e-10`，outer source-column reconstruction relative residual 不大于 `1.0e-14`；
- A33、B33、F3、RAO3、D3 和门槛字段全部有限；
- cache 状态为 `MISS`，矩阵维数等于 frozen manifest 的总 panel count；
- solver warning count 为零，或警告被逐条保存且经审核判定不影响代数解。未被审核的 warning 使 case 为 `NUMERICAL_FAILURE`。

代数门槛任一失败即 `SENTINEL_ALGEBRAIC_GATE=FAIL`。真实内存/运行资源不足单独记为 `BLOCKED_RESOURCE`，不得改写成物理或收敛 FAIL。

## 2. 响应量、单位、floor 与相邻敏感度

主量和单位为：A33 `[kg]`、B33 `[kg/s]`、F3 `[N/m wave amplitude]`、RAO3 `[m/m]`、D3 `[N/m]`。采用 `exp(+iωt)` 约定，单体全局 DOF 3 为 heave，力和位移正方向均为全局 `+z`，D3 为 RAO 动态矩阵的 `(3,3)` 元素。

相邻级 `coarse→fine` 的主敏感度定义为：

```text
S(Q) = |Q_fine - Q_coarse| / max(|Q_fine|, Q_floor)
```

固定 floors：A33=`1 kg`，B33=`1 kg/s`，F3=`1 N/m`，RAO3=`1e-6 m/m`，D3=`1 N/m`。这些值仅防止低幅值除零，不是精度放宽。每一行必须保存实际使用的 floor，并以 `0.1×floor` 和 `10×floor` 重算；若 PASS/FAIL 或分类改变，标记 `FLOOR_SENSITIVE`，不得单独以 phase 判定。

F3 与 RAO3 必须保存 real、imaginary、magnitude、wrapped phase、沿频率排序的 unwrapped phase、absolute complex difference 和 relative complex difference。总体收敛以 complex norm 为主。只有当 coarse/fine 两端幅值均大于 `100×Q_floor` 时，phase 才可作为辅助连续性诊断。

L2→L3 的预声明主门槛：

| 量 | 最大相对敏感度 |
|---|---:|
| A33 | 5% |
| B33 | 5% |
| complex F3 | 5% |
| complex RAO3 | 10% |

L1→L2 不直接决定最终 PASS，但必须保存并用于判断趋势。`convergence ratio = absolute difference(L2,L3) / absolute difference(L1,L2)`；分母低于对应 floor 的 `1e-6` 时记 `NOT_DEMONSTRABLE`。

## 3. 分类规则

- `CONVERGED`：L2→L3 满足量的门槛、floor 扫描不改变结论，且适用时可证明渐近区。
- `CONVERGED_WITHOUT_ASYMPTOTIC_PROOF`：L2→L3 满足门槛，但三点 GCI/observed order 因非等比、跨零、非单调或阶次不可解而不适用。
- `OSCILLATORY_CONVERGENCE`：相邻实增量反号；对复数以连续复平面轨迹和主 complex sensitivity 报告，不将复数直接做 GCI。
- `PRE_ASYMPTOTIC`：L2→L3 小于 L1→L2 但仍超过预声明门槛，或 observed order/GCI 尚未稳定。
- `RESONANCE_SHIFT`：1.30–1.50 rad/s 内的固定频点误差可由峰位置变化解释，且代数门槛通过；必须转入峰跟踪，不能把该单点直接标为普通收敛 PASS。
- `NONCONVERGENT`：L2→L3 超过门槛且不满足预声明 resonance-shift 延后条件，或差异随细化不减且无低幅值解释。
- `NUMERICAL_FAILURE`：几何前置失败之外的非有限量、conditioning/residual/cache/solver failure。
- `BLOCKED_RESOURCE`：真实资源不足阻止合法 solve。

单个 sentinel 量在 1.4 或 1.5 rad/s 可暂记 `RESONANCE_SHIFT` 而不使总体 sentinel 立即 FAIL，仅当：全部代数门槛 PASS、L2→L3 complex RAO sensitivity 不超过 25%、L2→L3 absolute complex difference 小于 L1→L2，且无复平面不连续。该状态强制运行 resonance tracking；峰跟踪失败后总体必须 FAIL。

`SENTINEL_RESPONSE_GATE=PASS` 要求所有强制量满足主门槛或被分类为 `RESONANCE_SHIFT`；存在 `PRE_ASYMPTOTIC`、`NONCONVERGENT`、`NUMERICAL_FAILURE` 时为 FAIL；存在资源阻断时为 BLOCKED。

## 4. Observed order、Richardson 与 GCI

real scalar 的三点计算只在以下条件全部满足时适用：三值有限、没有跨零、两个相邻增量同号、两个 realized refinement ratio 均大于 1、generalized nonuniform-ratio order 方程在 `0.1≤p≤10` 有唯一有限解。使用每一级冻结的 `eta_control`：`r21=h1/h2`、`r32=h2/h3`，不假设两比值相等。

满足条件时，fine-grid extrapolation 为 `Q_ext=Q3+(Q3-Q2)/(r32^p-1)`，fine-grid GCI 使用安全因子 `Fs=1.25`。A33、B33 可直接计算；F3、RAO3 只能分别对 Re 与 Im 在满足条件时计算，且必须独立列出。不得对 complex value 本身直接应用实数 GCI。任一适用条件不满足时写：

```ini
ASYMPTOTIC_RANGE = NOT_DEMONSTRATED
```

GCI 不可证明不自动等于 mesh convergence FAIL；但不能替代 L2→L3 主门槛，也不得删除任一级来制造阶次。

## 5. Resonance tracking 判据

冻结 coarse grid 为 `1.30:0.02:1.50 rad/s`，三等级分别 clean solve。比较 RAO3 magnitude peak、complex trajectory、D3 minimum、A33/B33、damping contribution 与 unwrapped phase continuity。L2→L3 门槛为：

- peak frequency 差不超过 `0.02 rad/s`；
- peak RAO3 magnitude 相对差不超过 `10%`；
- D3 minimum frequency 差不超过 `0.02 rad/s`；
- minimum `abs(D3)` 相对差不超过 `15%`；
- peak 不得落在 bracket 端点，所有点通过代数门槛且 complex trajectory 无非数值跳变。

允许局部 `0.01 rad/s` 加密，但只能在下列客观条件之一出现时触发，并且三等级使用同一加密点集：L2/L3 coarse peak 相差恰好一个 `0.02` 网格步；peak 相邻点任一达到 peak 的 98%；相邻 RAO unwrapped phase 变化超过 45°；D3 minimum 相邻点在最小值的 5% 以内。加密范围为触发候选频率并集的 `±0.02 rad/s` 与原 bracket 的交集。不得因看到峰形后任意移动 bracket 或选择性加点。

必须区分 `peak-frequency shift`、`peak-height nonconvergence`、`damping nonconvergence`、`fixed-frequency samples on opposite sides of peak` 与 `true numerical discontinuity`。

## 6. Resource、formal continuation 与降级

资源表必须先记录每个 sentinel 的实际 N、dense complex 单矩阵字节数 `16N²`、至少 A/L/U 三个 N×N 数组的下界，并引用 bridge 的 wall-clock。由于 bridge 未保存可靠 peak memory，正式批量前只允许一个最大 N (`ω=0.5,L3`) pilot。L3 最大并发固定为 1；其他并发也不得与 L3 重叠，除非新增可靠内存证据并另行审核，本阶段默认不启用。

formal sweep 仅在以下条件全部满足后继续：

```yaml
V5_MESH_REPRODUCIBILITY: PASS
SENTINEL_ALGEBRAIC_GATE: PASS
SENTINEL_RESPONSE_GATE: not FAIL and not BLOCKED
RESONANCE_CONVERGENCE: not FAIL and not BLOCKED
REVIEWER_PRECHECK_BLOCKERS: 0
```

正式频率固定为 `0.5:0.1:2.0 rad/s`。默认运行 `3 levels × 16 frequencies`。若资源 pilot 或 sentinel 的实测证据表明 L3 formal 总成本超过 8 小时，或单 case peak working set 超过当时可用物理内存的 70%，则使用预声明降级：三等级六个 sentinel 加两个 finest levels 的 16-frequency formal。降级不放宽门槛、不生成 L3 缺失值、不声称完整三点 Richardson/GCI。

## 7. L4 触发与 production-mesh candidate

只在以下任一项出现时触发 L4 设计请求；本阶段不得自动改变 topology 或凭空生成 L4：非共振强制量为 `PRE_ASYMPTOTIC`；L2→L3 差异不减；resonance peak/D3 minimum 超门槛；GCI 指示 L3 误差超过对应主门槛且结果 otherwise valid；或 L3 是唯一可能候选但没有渐近证据。L4 需要单独资源评审和明确 mesh controls，未授权时状态为 `L4_REQUIRED_NOT_RUN`。

production-mesh candidate 必须同时通过：可复现、代数、sentinel、resonance、formal、运行时间/内存证据。优先推荐满足全部门槛的最粗级。L2 只有在全部 mandatory L2→L3 比较通过时才可推荐；若 L2 不通过，L3 只有在有效 GCI 或后续 L4 证明其误差满足门槛时才可推荐。本阶段绝不自动切换 production default。

## 8. 禁止的事后操作

禁止在看到结果后改变 floor、门槛、频率 bracket、peak 定义、加密触发、资源降级条件或分类逻辑；禁止删除某一级/某频点、调 physics/topology/参数、使用旧 cache 或旧结果、把资源失败伪装为收敛失败、把 GCI 缺失等同于 FAIL、对 complex value 直接做实数 GCI、复制 WAMIT raw、用 WAMIT 反调门槛、把 mesh convergence 宣称为完整 physical validation 或 production ready。
