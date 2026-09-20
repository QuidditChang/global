# Dannberg (2024) Appendix C：cmbhf_EBA 方法核对、审计与改造方案

日期：2026-09-19。本文保存实施前的方法、审计与计划；最新实现和实测结果见 [VALIDATION_REPORT.md](VALIDATION_REPORT.md)，续接状态见 [IMPLEMENTATION_STATUS.md](IMPLEMENTATION_STATUS.md)。

## 1. 范围与证据版本

- 指定论文：Dannberg et al., *Changes in core–mantle boundary heat flux patterns throughout the supercontinent cycle*, GJI 237, 1251–1274, DOI 10.1093/gji/ggae075。
- 已读取本地 PDF，目视核对 PDF 第 21–22 页（印刷页 1271–1272）C1–C5；读取第 23–24 页的基准表 C1/C2。
- solver 工作树：`worktrees/src-global-cmbhf_EBA`，分支 `cmbhf_EBA`，基线 `b23c2b118f5717d7642b164f347a9bf815791dd0`。
- runs 工作树：`worktrees/runs-cmbhf_EBA`，分支 `cmbhf_EBA`，基线 `650f9f0240f6aed2eaed731eb3f2a2687a950f74`。
- 两工作树初始均干净。主 `src/global` 为 `cmbhf_ALA_strict`，不切换、不修改。
- 根目录历史 `CURRENT_CBF_ARCHITECTURE.md`、`STRICT_CBF_*` 针对其他历史版本，不能作为 EBA 当前实现或论文的权威。

## 2. 先确定论文方法

### 2.1 真正的 CBF

CBF 在给定已经计算出的温度、时间导数、速度和热源后，把 Dirichlet 边界节点原先被消去的能量方程重新用于恢复边界热流。它不是边界温度梯度的简单差分。

从能量方程

`rho Cp (T_t + u.grad T) - div(k grad T) = H`

与物理外向导热通量 `q_out = -k grad(T).n` 出发，采用一致的物理能量弱式：

`B q_out = b`

`B_ij = integral_GammaD N_i N_j dA`

`b_i = integral_Omega [N_i H - N_i rho Cp(T_t + u.grad T) - k grad N_i.grad T] dV - integral_GammaN N_i qN_out dA`。

EBA 相变项 `L` 若位于能量方程左边，须在 RHS 再减 `integral N_i L`。输入 Neumann 值若按向内正定义，已知通量项相应为加号。必须转换定义后再比较，不能仅照抄符号。

论文 C3 使用除以 rho Cp 的简写；在可变 rho/Cp 时，不能把强式简单除法后忽略系数梯度，并称为等价的守恒能量弱式。移植应以 CitcomS 实际物理能量方程为权威，明确归一化，输出 W/m²，不把温度通量和能量通量混用。

### 2.2 必须纠正的旧方案结论

论文第 1272 页明确说明：边界采用 Gauss–Lobatto–Legendre (GLL) 积分，积分点与有限元支撑节点重合，因此边界质量矩阵为对角矩阵。

因此“严格 CBF 必须保留非对角项并用 PCG 求解”不符合本次指定论文的具体实现。完整一致质量矩阵也是一种 CBF 离散，但不是本文采用的边界积分方案。不能仅因现有代码做逐节点除法就判定它不是 CBF；须核对其分子、分母、几何和积分规则。

论文温度元为 Q2，边界用 3×3 GLL，单方向点为 -1,0,1，权重为 1/3,4/3,1/3。CitcomS 当前温度元为 Q1，严格移植“与温度元同阶、支撑点共置的 GLL”原则应使用 2×2 端点规则：-1,+1，权重均为 1。不能声称 Q1 实现复现了论文 Q2 收敛阶。

Q1 每个面节点 i 的局部对角质量为 `D_i = J_surface(xi_i,eta_i)`。仿射平面矩形上 `D_i=A/4`，与行和集中质量相同；一般非仿射/曲面几何上不保证与高斯积分 `integral N_i dA` 相同。

### 2.3 印刷公式疑点与处理

目视可见 C3 的物质导数是 `T_t + u.grad T`，C5 外部负号内却印为 `T_t - u.grad T`，不符合从 C3 移项的代数。C3/C5 的 Neumann 积分还未显式保留测试函数。按能量方程和 C3 的物质导数推导执行：源项正、瞬态项负、对流项负、导热刚度项负。报告中保留此说明；不得悄悄按疑似排印错误实现反号对流。

### 2.4 Galerkin 与 SUPG 不是同一件事

论文 C3 写普通温度基函数 N_i。CitcomS 时间推进使用 Petrov/SUPG 测试函数。应共享热物性、相变项、几何、速度处理和所有能量项的生产内核，但明确选择：

- 论文 C3 产品：普通 Galerkin N_i 的物理弱式 + GLL 边界质量；这是本任务主目标。
- SUPG 边界反力：用 PG_i 的离散残差，是另一个诊断，可用于衡量稳定化影响；不能在没有说明时当作论文原式。

不应改变温度推进算法以迁就后处理。共享 API 必须防止诊断调用改写 heating_phase 或其他求解状态。

## 3. 当前 EBA 代码审计

以下定位针对上述基线，后续编辑可能移动行号。

| 项目 | 证据 | 结论与影响 |
|---|---|---|
| 已有 CBF 入口 | `CitcomS/Controller.py:151` → `Solver/Solver.py:248` → `module/outputs.c` → `lib/Output_gzdir.c:1262` | 确实已实现旧诊断；不是从零开始 |
| cfg 实际关闭 | `runs/cmbhf_EBA.cfg:19–22` monitoringFrequency=50、CBF=0 | 每 50 步普通监测并不会自动输出 CBF |
| 球坐标梯度错误 | `Process_buoyancy.c:361–399`；对照 `Size_does_matter.c:50`、`Advection_diffusion.c:982` | sphere=1 返回球坐标导数；旧注释误称 Cartesian，遗漏 1/r、1/(r sin theta)；横向导热及对流错误 |
| 遗漏相变能量 | `Advection_diffusion.c:1024–1098` 对照旧 CBF | 当前 EBA 使用显式 phase_energy，旧 CBF 仍镜像旧 heating_latent 形式；非零熵跳时方程不一致 |
| 可关闭对流 | `Instructions.c:473`、`Process_buoyancy.c:405` | 完整 CBF 不能任意省略实际激活的对流项 |
| 速度处理不同 | `Advection_diffusion.c:764–775` | 求解器限制过大的局部速度，旧 CBF 不限制；残差可能来自不同有效速度 |
| 边界规则不符 | `Process_buoyancy.c:427–436` | 2×2 Gauss 积分 N_i 的行和；不是支撑点 GLL 的 D_i |
| 面几何近似 | `Size_does_matter.c:378` | 旋转后只用二维投影行列式；新方案需与 Q1 三维面映射一致的叉积面积元，并验几何误差 |
| 径向 MPI 风险 | `Process_buoyancy.c:287,443–444`，`Full_parallel_related.c:1046` 后的 vertical Sendrecv | 非边界 rank 先返回，边界 rank 却调用包含径向通信的完整交换；nprocz>1 可发生通信不匹配/挂起。当前 cfg nprocz=2 |
| RHS 支撑 | 仅边界层单元组装所有 8 个体节点 | 边界通量仅需对应边界行；不要向内部径向层传播无用残差 |
| 单精度累计 | `float *b_CBF,*A_CBF`、exchange_node_f | 小差值相消、MPI 划分敏感；应改 double |
| 非法面积静默处理 | A<=0 时写 0 | 掩盖网格/装配错误；应所有 rank 一致失败，不发布伪零通量 |
| 已知通量边界 | 旧 CBF 无 Neumann 项、无 Dirichlet 检查 | 当前上下定温可支持；混合/定通量情形需显式支持或拒绝，不能静默当定温 |
| 热源时刻 | 时间推进最终刷新 adi/visc，CBF 独立重建 | 必须明确初始 step 0、温度更新与 Stokes 更新之后的取值时刻，不能沿用未初始化热源 |
| 同化/过滤 | `PG_timestep_solve:375–436` | T 被后处理而 Tdot 未同幅修正；不能无条件宣称闭域能量守恒。需保存/标记热推进态并独立报告同化增量 |
| 监测均值 | q 四节点平均乘 eco.area | 不是所选 GLL 的直接积分；应以局部 D_i q_i 求和，每个面只算一次 |
| 输出不足 | `Output_gzdir.c:1218–1275` | 只有 q_surf/q_botm_CBF.<rank>.<step>.gz；无 lon/lat、无全局 GRD |
| 老脚本证据失效 | 历史文档所称 scripts/make_q_CBF_grd.py 当前不存在 | 不把历史脚本结论冒充当前审计；无现成可信 GRD 桥接 |
| 基准覆盖不足 | 当前 tests/phase_energy 有生产相变块测试 | 未见 CBF 边界 GLL、球壳热流和径向 MPI 专项验收 |

当前结果不能作为经验证的 Appendix C 生产热流产品；尤其不能仅打开频率并给 gzip 改扩展名。

## 4. 输出契约

默认产品（同一运行目录的相对路径）：

```
PostProc/HF_CBF/q_CBF_<step>.grd
PostProc/HF_CBF/eshf_CBF_<step>.grd
```

- 保留用户原命名；cmbhf 为底部、eshf 为顶部。
- 两者单位 W/m²。顶部正=地幔向外，底部正=核向地幔。内部统一 outward，再对底部翻号一次。
- 与 `monitoringFrequency` 同步；独立 CBF 周期仅为显式覆盖。建议旧键 -1=继承、0=禁用、正整数=独立周期；默认继承，EBA 主 cfg 设置继承。
- 重启按实际 global step 命名；拒绝时间标签与状态不符。step 0 若无可靠 Tdot/热源不能冒充瞬态 CBF；须初始化一致导数或明确标记初始诊断。
- GMT 可读 NetCDF classic，x=longitude、y=latitude、z(y,x)，推荐 0.5°、0..360°、-90..90°，721×361，gridline registration。
- 0°/360°同值；每个极点各 longitude 同值；不得用任意最近邻填补未覆盖点而隐藏拓扑缺口。
- 元数据：method=CBF_GLL_Q1、weak_form、step、time_nd、time_seconds、boundary、positive_direction、units、radius_m、k0、DeltaT、R0、native_integrated_heat_W、mapping_method、initial/transient 状态、同化/过滤状态。
- GRD 是经纬采样场。原生 FEM/GLL 积分是热流总量权威；普通点插值不保证经纬栅格保守，必须标注并测量映射积分误差。不要人为整体乘系数迫使守恒而改变局地物理值。

## 5. 详细改造顺序

### A. 共享物理 RHS 内核

1. 在 EBA `Advection_diffusion.c` 内部提取可调用的单元热量残差 API，保持现有时间推进调用语义。
2. 提供测试权重参数/模式 N 或 PG；CBF 主产品调用 N；有效速度准备、rho/Cp、k(T,d,C)、源项、phase_energy 完全共享。
3. `element_residual` 当前会写 heating_phase：增加显式诊断只读模式，或把该赋值移回推进包装层。
4. CBF 输出时刷新或复用有明确状态标签的热源；所有全局 collective 必须由所有 rank 进入。
5. 同化和过滤设计先落实状态契约。至少同步记录热推进后的未修改态与最终态的差异；若仅输出最终态瞬时弱式，明确它不是严格时间离散反力。不能伪造 Tdot 或默认零瞬态项。
6. 删除/拒绝 CBF_use_advection=off；当前 cfg on 无迁移损失。

### B. GLL 面质量和 MPI

1. Q1 面四角使用与体单元相同的三维等参映射；两切向导数叉积模长得到面积元。
2. 在四个 GLL 端点计算 D_i，不用现有 Gauss 行和替代。
3. 上下边界独立累计 RHS 和 D，全程 double。仅装配边界支撑行。
4. 简洁可靠的首版可让所有 rank 以零初始化数组参与已有完整 exchange_node_d，避免径向提前返回；只在真实物理边界求商/输出。后续如需优化再专用 horizontal exchange，不改全局交换实现。
5. 所有本地异常先 Allreduce 错误标志，再统一失败；禁止一个 rank 退出其余 rank 等待。
6. `q_out=b/D`；D 必须正且有限，q 必须有限。
7. 热流积分采用每个拥有面上的局部 D_i 乘已组装 q_i；面遍历天然唯一，不能把已全局相加的 nodal D 在重复节点再求和。

### C. 经纬映射与 GRD

1. 映射消费原生 FEM 通量，绝不从输出的温度/速度重新计算另一套 CBF。
2. 采用面内 Q1 插值：经纬网格方向与三维边界面求交，逆等参求局部坐标，再按 N_i 插值。避免 lon/lat 平面三角剖分在日期线和极点出错。
3. 并行方案：每 rank 处理本地边界面覆盖的目标点，归约贡献、覆盖计数及接缝差；重复命中只允许在共享边/顶点并要求值一致。
4. 只向写者归约固定大小 GRD 数组，避免收集整个高分辨率三维模型。边界拓扑或空间索引可缓存。
5. 首选已有可用 NetCDF 库；若增加构建依赖需同步 configure/Makefile 与 HPC 配方。不能为写两张固定表引入整个 Python 后处理环境。
6. 写到临时文件，关闭并验证成功后 rename；双文件失败时不宣称本步完整。输出目录创建由唯一 rank 执行，结果广播。
7. 保留旧绑定名称可减少调用层迁移，但改为新全球产品，旧 rank gzip 若保留只能作为显式调试开关。

### D. 控制层和配置

- `CitcomS/Controller.py`：解析继承监测频率，在所有 rank 同一步调用。
- `CitcomS/Solver/Solver.py`：同步 inventory 与 C 参数，保证关闭模式不分配/通信。
- `lib/Instructions.c`, `lib/global_defs.h`：路径、grid spacing、错误验证、double 存储与释放。
- `lib/Output_gzdir.c`, `module/outputs.c`：同步生成双边界 GRD，保持正常温度/速度输出无变化。
- `lib/Makefile.am`/构建生成文件：登记新文件和必要依赖。
- `runs-cmbhf_EBA/cmbhf_EBA.cfg`：只修改 CBF 相关键；其他 EBA 变体 cfg 单独列出是否启用，避免默认意外覆盖历史实验。

### E. 验证（必须记录实际执行结果）

1. Q1 单面：矩形 D=A/4；一般三维扭曲面用独立叉积/高阶参考验证；GLL B 非对角元为零；方向翻转不改正面积。
2. RHS：直接调用生产内核，逐项测试瞬态、横向对流、径向/横向导热、可变 k、内部热源、adi/visc、非零相变熵。对照普通 N 与 PG 差异，不把两者混称。
3. 符号与尺寸：稳态球壳 T=A+B/r，q_top=kB/r_top²、q_bottom=kB/r_bottom²，两输出正；积分相等。用至少三套网格报告误差阶，不套论文 Q2 阶。
4. 变 k 球壳：满足 r² k T'=常数；检验热流尺度 k0 DeltaT/R0 只乘一次。
5. 制造瞬态/带源问题；明确解析采样和数值 PDE 求解是不同测试，不拿一个代替另一个。
6. MPI：同一全局网格改变水平及径向分解，至少 nprocz=1/2；设置超时捕捉死锁；比较原生节点、积分、GRD 接缝。全局 12-cap 的最小 rank 数按实现约束，不盲目要求不可用的 1-rank 配置。
7. I/O：独立 NetCDF/GMT 读取尺寸、坐标、单位和符号；验证 poles、seam、无 NaN、实际 cfg 步 0/50/100 与 restart 非零步命名。模拟写失败，检查无已发布截断文件。
8. 不侵入性：同一小算例开启/关闭 CBF，T/u/Tdot/热源/步长一致；检测诊断函数副作用。
9. 生产前：小型 EBA 跑若干接受时间步并输出真实 GRD，再 HPC 同配置验证。用含储热、全部热源、同化/过滤和边界平流的收支解释总量；不能规定未定义的闭域 1% 阈值掩盖外部强迫。

## 6. 完成门槛和非目标

完成必须同时包括：论文方法有明确移植说明、旧几何/相变/MPI 错误已修、监测步实际产生两张可读 GRD、数学单测与 MPI/小模型证据、没有修改 ALA 分支或其他求解模块。

不把全套 Q2 求解器重写、PCG 边界系统、确定性逐位 MPI 求和或新地核耦合协议列为本任务必需项。它们不是用户本次要求，也不是 Appendix C GLL 方法成立的前提。

## 7. 续接记录

已建立本任务每小时自动续接。额度恢复时先查本文件与 git diff，再继续未完成阶段；不重复创建自动任务，不兑换额度。完成后停用 automation `cbf`。

截至本报告写入：尚未修改生产 C/Python/cfg，尚未运行 CBF 数值验收，尚未生成用户要求的真实模型 GRD。下一步是落实 RHS/输出状态契约并实施 A/B；不得将本方案标记为实现完成。
