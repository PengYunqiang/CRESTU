```markdown
===============================================================================
=                                                                             =
=            CCCC   RRRR    EEEEE   SSSSS  TTTTT  U   U                       =
=           C       R   R   E       S        T    U   U                       =
=           C       RRRR    EEEE    SSSSS    T    U   U                       =
=           C       R  R    E           S    T    U   U                       =
=            CCCC   R   R   EEEEE   SSSSS    T     UUU                        =
=                                                                             =
=                            C R E S T U - 1 F                                =
=     (Computational Rankine-source Engine for Sea Technology, SJTU)          =
=                 3D Rankine Source BEM Hydrodynamic Solver                   =
=                                                                             =
=            Project Name : CRESTU-1F                                         =
=            Lead Authors : Yunqiang Peng, Zhentao Jiang                      =
=            Affiliation  : Shanghai Jiao Tong University                     =
=            Status       : Research Software — Active Development            =
=            Start Date   : Aug. 24, 2026                                     =
=                                                                             =
===============================================================================
```
# CRESTU-1F v1.0 User Manual

> **[Status: Active Development / Early Research Alpha]** 
> 本仓库核心一阶框架持续更新中。部分高阶模块和测试算例正在逐步完善并随开发进度推送。
> *Core first-order modules are under active development. Additional verification cases and modules will be pushed incrementally.*

> 中文说明在前，English documentation follows below.

---

## 1. 项目简介

**CRESTU-1F**（**C**omputational **R**ankine-source **E**ngine for **S**ea **T**echnology, **SJTU**）是由上海交通大学流固耦合团队博士生开发的三维 Rankine 源边界元（BEM）水动力学计算程序。当前版本基于 MATLAB 平台架构开发，主要功能包括：

- 定常均匀流 Hess–Smith 面元法计算；
- 单体及多体频域一阶辐射、绕射求解；
- 附加质量、辐射阻尼、波浪激励力、静水恢复力及 6N 自由度运动响应（RAO）；
- `ISX/ISY` 几何对称加速；
- 浮体物面压力导出；
- 二阶平均漂移力计算（Pinkster 近场法与远场法）；
- 标准浮体（单球/多球）基准算例与收敛性验证。

> **提示**：本程序当前处于科研开发与验证阶段（Active Research Development），适用于学术研究、算法验证与教学。

---

## 2. 核心符号约定

- **时谐因子**： $e^{\mathrm{i}\omega t}$
- **线性动水压力**： $p=-\mathrm{i}\omega\rho\Phi$
- **面元法向**： 物面面元法向向量一律指向流体域
- **自由度定义**： 1–6 依次为 Surge, Sway, Heave, Roll, Pitch, Yaw
- **多体自由度排列**： 第 $k$ 个浮体占据全局自由度序号 $6(k-1)+(1:6)$
- **单位系统**： 国际标准单位制（SI units）
- **浪向角**： 波浪传播方向，以全局 $+x$ 轴为基准逆时针计量
- **无量纲漂移力系数**： $C_d = \frac{F^{(2)}}{0.5\rho g A^2 L}$

---

## 3. 目录架构

```text
CRESTU/
├─ 0.Tools Code/       网格与几何可视化工具
├─ 1.Input/            算例配置文件解析
├─ 2.Mesh/             BMF 网格读取、生成、变换与水线处理
├─ 3.HessSmith/        定常源面元速度核函数
├─ 4.Potential/        Rankine 边界积分方程（BIE）组装与势求解
├─ 5.Force/            一阶水动力系数、激励力、静水力与 RAO 求解
├─ 6.MeanDriftLoads/   二阶平均漂移力计算模块
├─ Case_UniformFlow/   定常均匀流球体测试算例
├─ Case_Wave/          频域波浪主求解器与基准测试算例
├─ CRESTU-1F.cfg       顶层参数配置模板
└─ README.md           项目说明文档

```

---

## 4. BMF 网格文件规范

`.bmf`（Boundary Mesh File）采用标准四节点面元格式，第二行包含几何特征尺度与网格类型标识：

```text
Line 1 : Header text
Line 2 : ULEN   PanelType
Line 3 : ISX    ISY
Line 4 : NPAN
Next   : 每个面元连续 4 行，每行给出 X Y Z 节点坐标

```

**面元类型标识（PanelType）：**

* `1`：物面（Body surface, `*_body.bmf`）
* `2`：自由液面（Free surface, `*_fs.bmf`）
* `3`：自由面阻尼层（Damping lid, `*_damping.bmf`）
* `4`：海底边界（Seabed, `*_seabed.bmf`）
* `5`：外域圆柱截断控制面（Far-field boundary, `*_farfield.bmf`）

---

## 5. 控制参数说明（`.cfg`）

| 参数段 | 功能描述 |
| --- | --- |
| `PARA1` | 算例名称及输出文件前缀 |
| `PARA2` | 执行控制开关：`IPOTEN`（势流求解）、`IFORCE`（载荷计算）、`IRAD`（辐射）、`IDIFF`（绕射）、`IDRIFT`（漂移力） |
| `PARA3` | 入射波频率范围、步长及浪向角数组 |
| `PARA4` | 水深 $h$、重力加速度 $g$、流体密度 $\rho$ |
| `PARA5` | 浮体数量及各浮体空间位置、初始偏航角 |
| `PARA6` | 浮体质量、全局重心位置与惯性张量 |
| `PARA7` | 激活的辐射模态（1–6） |
| `PARA8` | 几何对称标志 `ISX ISY` |
| `PARA9` | 水线容差、近场与吸收层网格离散参数 |
| `PARA10` | 自由面人工阻尼系数 $\mu_0$ |

---

## 6. 快速开始

### 6.1 环境要求

* MATLAB R2022b 及以上版本（已在 Windows 平台验证）；
* 建议配备 Parallel Computing Toolbox 以支持矩阵装配并行计算。

### 6.2 运行算例

在 MATLAB 命令行中进入工程目录并执行：

```matlab
% 1. 进入主求解器目录
cd('Case_Wave');

% 2. 运行默认配置算例
results = HydroMain();

% 3. 或指定配置文件路径运行
cfg_path = fullfile(pwd, 'Case1_SingleSphere_Convergence', 'Mesh_Fine', 'Case1_Fine.cfg');
results = HydroMain(cfg_path);

```

### 6.3 运行基准验证

```matlab
% 运行全部验证集并导出对比图表
validation = Run_All_Validation(false);

```

---

## 7. 输出数据格式

| 文件名 | 内容描述 |
| --- | --- |
| `*_PotCache.mat` | 辐射势、入射势与绕射势函数解缓存 |
| `*_Results.mat` | 附加质量、阻尼、激励力、静水恢复刚度、RAO 与压力场 |
| `*_MeanDrift.mat/.csv` | 二阶平均漂移力及无量纲系数 |
| `*.png` | 水动力系数与运动响应验证曲线图 |

---

## 8. 开源许可与免责声明

* **开源协议**：本项目基于 **GNU General Public License v3.0 (GPLv3)** 开源。
* **免责声明**：本程序为科研开发版本，仅供学术交流、理论研究与教学使用。用于实际工程设计时应结合物理模型试验或权威基准进行独立校核。

---

---

# CRESTU-1F v1.0 Package Guide — English

## 1. Overview

**CRESTU-1F** (**C**omputational **R**ankine-source **E**ngine for **S**ea **T**echnology, **SJTU**) is a 3D Rankine Source Boundary Element Method (BEM) hydrodynamic solver developed at Shanghai Jiao Tong University. Implemented in MATLAB, the current framework provides:

* Steady uniform-flow Hess–Smith source-panel calculations;
* Single- and multi-body first-order radiation and diffraction analysis;
* Added mass, radiation damping, wave excitation forces, hydrostatics, and 6N-DOF RAOs;
* `ISX/ISY` geometric symmetry acceleration;
* Complex wetted-surface dynamic pressure evaluation;
* Near-field (Pinkster) and far-field second-order mean drift load calculations;
* Verification benchmark suite (floating single-/two-sphere cases).

> **Note**: This software is an active research release intended for academic research, numerical testing, and education.

---

## 2. Mathematical & Coordinate Conventions

* **Time-harmonic factor**: $\mathrm{e}^{\mathrm{i}\omega t}$
* **Linear dynamic pressure**: $p = -\mathrm{i}\omega\rho\Phi$
* **Panel normal vector**: Positive outwards from body into fluid domain
* **DOF definitions**: 1–6 denote Surge, Sway, Heave, Roll, Pitch, Yaw
* **Multi-body DOF sequence**: Body $k$ occupies global degrees of freedom $6(k-1)+(1:6)$
* **Units**: International System of Units (SI)
* **Wave heading**: Counter-clockwise angle measured from the global $+x$ axis
* **Nondimensional drift coefficient**: $C_d = \frac{F^{(2)}}{0.5 \rho g A^2 L}$

---

## 3. Directory Layout

```text
CRESTU/
├─ 0.Tools Code/       Geometry and mesh inspection utilities
├─ 1.Input/            Configuration parser and verification
├─ 2.Mesh/             BMF I/O, generation, spatial transforms, and waterline tools
├─ 3.HessSmith/        Steady source-panel velocity kernels
├─ 4.Potential/        Rankine BIE assembly and linear potential solvers
├─ 5.Force/            First-order hydrodynamic coefficients, excitation, and RAOs
├─ 6.MeanDriftLoads/   Second-order mean drift load computation
├─ Case_UniformFlow/   Steady uniform flow test cases
├─ Case_Wave/          Wave solver main scripts and validation benchmarks
├─ CRESTU-1F.cfg       Top-level configuration template
└─ README.md           Project documentation

```

---

## 4. Boundary Mesh File (.bmf)

The `.bmf` format follows a standard 4-node panel convention with panel-type identification on line 2:

```text
Line 1 : Header text
Line 2 : ULEN   PanelType
Line 3 : ISX    ISY
Line 4 : NPAN
Next   : 4 rows per panel specifying nodal coordinates (X Y Z)

```

**PanelType Identifiers:**

* `1`: Body surface (`*_body.bmf`)
* `2`: Free surface (`*_fs.bmf`)
* `3`: Damping surface (`*_damping.bmf`)
* `4`: Seabed boundary (`*_seabed.bmf`)
* `5`: Far-field cylindrical control boundary (`*_farfield.bmf`)

---

## 5. Configuration Parameters (`.cfg`)

| Section | Description |
| --- | --- |
| `PARA1` | Case name and prefix for generated files |
| `PARA2` | Execution switches: `IPOTEN`, `IFORCE`, `IRAD`, `IDIFF`, `IDRIFT` |
| `PARA3` | Frequency range, spacing, and wave heading angles |
| `PARA4` | Water depth $h$, gravity $g$, and fluid density $\rho$ |
| `PARA5` | Number of bodies, positions, and yaw angles |
| `PARA6` | Body mass, center of gravity, and inertia tensor |
| `PARA7` | Active radiation modes (1–6) |
| `PARA8` | Symmetry flags (`ISX`, `ISY`) |
| `PARA9` | Waterline tolerance and mesh domain sizing |
| `PARA10` | Free-surface Rayleigh damping coefficient $\mu_0$ |

---

## 6. Quick Start

### 6.1 Prerequisites

* MATLAB R2022b or later (tested on Windows);
* Parallel Computing Toolbox recommended for parallel matrix assembly.

### 6.2 Running a Simulation

```matlab
% 1. Navigate to solver directory
cd('Case_Wave');

% 2. Run default configuration
results = HydroMain();

% 3. Run a specific case
cfg_path = fullfile(pwd, 'Case1_SingleSphere_Convergence', 'Mesh_Fine', 'Case1_Fine.cfg');
results = HydroMain(cfg_path);

```

### 6.3 Automated Benchmarking

```matlab
% Execute complete verification suite and generate validation figures
validation = Run_All_Validation(false);

```

---

## 7. Output Files

| File | Contents |
| --- | --- |
| `*_PotCache.mat` | Cached radiation, incident, and diffraction velocity potentials |
| `*_Results.mat` | Hydrodynamic matrices ($A, B$), excitation forces, hydrostatics, and RAOs |
| `*_MeanDrift.mat/.csv` | Second-order mean drift forces and coefficients |
| `*.png` | High-resolution hydrodynamic response curves |

---

## 8. License & Disclaimer

* **License**: This project is distributed under the **GNU General Public License v3.0 (GPLv3)**.
* **Disclaimer**: This package is research software intended for academic and educational purposes. Engineering applications require independent verification against physical experiments or established commercial benchmarks.

```

```
