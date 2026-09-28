# P0：EBA-PICES 路线 B 数值契约 v1

基线：`cmbhf_EBA` / `3f7f4de2ab22b0cf45d18ffc829259c1fda4211b`。
实现分支：`cmbhf_EBA_PICES`。runs 分支：`cmbhf_EBA`。
日期：2026-09-28。

**本提交只冻结设计并建立 PG 验收基线，不实现 PICES 时间推进。**
`lib/`、`bin/`、`module/`、`CitcomS/` 及构建入口保持基线内容。
下面的 API、字段、检查点扩展、拒绝条件均为后续实现契约，当前不是可用功能。
P0 必须经用户 HPC 结果审计通过后才进入 P1；本地通过不自动开启下一阶段。

## 1. 范围及阶段限制

采用 Dong et al. (2024), DOI 10.26464/epp2024021 的粒子平流、网格扩散/源项及子网格增量扣除流程，复用 EBA 物理核，不合并旧 `cmbhf_PICES`。

P1 只开放全局 12-cap 网格、每 rank 一个 cap 或其子域、正向推进、常正热导率、常 rho/Cp、无相变/体热源/TA、固定径向 Dirichlet。粒子迁移基础复用，但跨缝和重启在 P2 独立验收。P1 不声称具备 P3/P4 的完整 EBA 物理。纯平流测试允许显式测试配置令 k=0。

后续参数 `energy_solver=pg|pices` 默认 pg，非法值报错。P0 cfg 不填写此尚未注册的参数。启用 PICES 时对阶段尚不支持的 restart、区域、多 cap/rank、backward、filter、rheo.dat 温度修改、边界或物理项应明确拒绝，不能静默走一半流程。

## 2. 状态、索引与变量语义

- 粒子编号始终 `1..ntracers`；槽位下标可从 0 开始，二者不可混淆。分配长度必须严格大于 ntracers。第 0 粒子不读、不写、不参与归约。
- 持久新属性先仅 `Tp`，在 tracer 分配前注册，放在已有 flavor 等属性之后。迁移、删除交换、扩容、读取检查点都由同一属性表负责。不可用固定槽号假定 flavor 存在。
- `Tp` 是无量纲完整温度，与 `E->T` 同一归一化，不是温度异常，也不是 Kelvin。采用 double 粒子存储；现有网格 float 精度不在 P0 改动。
- `Tg` 为该步网格热求解/同化的发布结果；`Tp` 保留子网格结构。允许 `P(Tp) != Tg`，但差额必须诊断和入账，不能额外投影后无声覆盖 Tg。
- `deltaT_heat`、`deltaT_sub_p`、`deltaT_assim`、重映射差额独立存储/统计。deltaT_sub_p 是本阶段临时量，在其计算和使用之间禁止移动或重排粒子。
- `Tdot_heat=deltaT_heat/dt` 为热阶段变化率，不能直接冒充 PG 的 Eulerian `Tdot`。P3/P4 的 CBF 与相变诊断使用明确的材料导数/阶段残差接口。

## 3. P/Q 映射及子网格公式

首版仅采用已有全局 tracer 的 gnomonic wedge × 径向线性局部形函数。插值助手从 Full_tracer_advection.c 小范围导出或抽取供共享使用，不拷贝一套独立几何公式。保留现有节点映射规则，验证每个 wedge、径向边界、cap 缝与常数保持。

P：累加 `sum(Tp*w)` 和 `sum(w)`，两者均使用 finest-level MPI 节点交换后再相除。Q：同一组非负、单位和权重插值节点值。仅在舍入级容差 `1e-12` 内夹断负权重再归一化；超界先重定位一次，仍失败则报错。不能用任意凸包夹断掩盖错误 host element。

初版若活动内部节点全局 `sum(w)<=0` 则终止并报告粒子/单元/节点/进程信息；不采用旧温度回填。固定温度边界可直接用已知边界值，并记录零覆盖数量。空单元与零节点权重分别统计。补粒子/合并不是 P1 的隐式副作用，需 P2 单独设计与能量账本。

一次热子步中粒子位置固定，冻结本子步的系数，首版用显式 Euler：

```text
Tg0      = P(Tp)，再施加该子步 Dirichlet
Tg1      = Tg0 + dt_s * M_L^-1 * [ -K(k_gp) Tg0 + F_sources + F_boundary ]
deltaTg  = Tg1 - Tg0
Tenv_p   = Q(Tg1)                  # 明确采用已求解环境温度
deltaTp_s = (Tenv_p-Tp) * [-expm1(-dt_s/tau_p)]
deltaTg_s = P(deltaTp_s)
Tp       = Tp + deltaTp_s + Q(deltaTg-deltaTg_s)
```

原始 deltaTp_s 不得被 `Q(P(deltaTp_s))` 覆盖，否则子网格项会抵消。
多子步时从上一子步 Tg1 继续网格求解，不在每个热子步重新 P(Tp) 覆盖网格解；仅 deltaTp_s 每子步重新映射。由此产生的 grid/particle mismatch 独立计量。

P0 冻结 tau 为 `Le²/kappa_eff`，不使用最近邻、应变修正或 TSC。
物理 Le 选为三组相对方向的平均边长中的最小值：每组四条平行拓扑边先计算 Cartesian 弦长平方的平均并开方，再取三组最小值。使用 E->x/Jacobian 物理几何，禁止直接对 `(theta,phi,r)` 差做欧氏范数。单位与半径归一化一致。
这是论文“单元长度”在本项目各向异性网格上的明确选择，不是声称论文规定了该定义。长度灵敏度比较在 P5 进行，变更定义必须改变记录的方法版本。

P1 单元 kappa 常数；P3 采用同一物理 kernel 的局部扩散率，不能引入覆盖 EBA 物性的自由 `pices_kappa`。几何可缓存，kappa/tau 每次热子步求值；首次热更新及 restart 后无条件构造。kappa=0 时松弛因子严格为零；Le<=0、负/非有限 kappa 或非法容量均报错。

## 4. 步序与稳定性：冻结为一次移动、热子步、不跨 MPI 回滚

首版采用一阶 Lie 分裂，已知速度 u^n 驱动轨迹 predictor/corrector：

1. 在已接受状态 n 选择 dt，受粒子 CFL 限制；P1 采用最小物理边长的 0.25 倍作为单步位移上限估计。
2. 粒子移动一次、MPI 迁移、重新定位，组分重建，Tp 暂不变。
3. P(Tp) 得到已平流场；覆盖检查、边界、系数与源更新。
4. 热区间 dt 内按稳定限制分热子步，完成网格热求解与粒子增量修正。
5. P4 才在整个 dt 完成后施加一次 TA、记录并传递实际增量；物理边界时间使用 t^n+dt，不提前提交全局时钟。
6. 发布 Tg/Tp，时钟与计数提交一次；使用同一时层的温度/组分求 Stokes，再输出和检查点。

为避免在已迁移粒子后套用只恢复 E->T 的 PG retry，首版不支持完整物理步自动回滚。固定 dt 若违反粒子 CFL，在移动前拒绝；热步使用子步。子步遇非法值、覆盖丢失、容量失败或超过子步数上限时终止，不写“已接受”检查点。不能继续用已失败的部分状态。以后添加回滚必须包含粒子 owner、坐标、属性、组分、全部阶段增量和时钟。

固定系数热算子 `A=M_L^-1 K`，先在组装且交换后的自由节点计算 `R=max_i(sum_j |Aij|)`，用 `dt_s <= 0.8/R` 作为首版保守谱稳定限制；R=0 无传导限制。跨 MPI 的物理条目不能重复计数。正质量、对称半正定 K 是这一界的前提，必须验证；它不保证离散最大值原理。还要取源项/相变非线性变化限制，P3 启用前验证。不能沿用“fixed_timestep 直接返回”而跳过 PICES 安全检查。

当子步剩余时间接近浮点尾差时按统一容差收尾，提交的总 dt 必须等于轨迹所用 dt。Pyre/C 两入口使用同一 C 编排，不能再次调用 advectTracers。粒子 predictor/corrector 不提升 Lie 分裂整体一阶时间精度。

## 5. 路线 B 共享物理 API 的职责

以下为接口设计，P0 不新增空实现或链接到生产二进制：

| 拟定接口 | 输入/输出契约 | 实施阶段 |
|---|---|---|
| `thermal_transport_at_gp` | E,cap,element,完整 Tgp -> rho_gp,Cp_gp,k_gp,kappa_gp；只读，无热状态副作用 | P1 常物性适配，P3 完整导出 |
| `thermal_phase_coefficients` | Tgp,物理 u_r,reference/phase -> 附加容量与压力源，另提供诊断值 | P3 |
| `thermal_assemble_diffusion_source` | 显式 field 与 source-stage -> K作用量、F、容量与边界贡献；不包含网格平流/SUPG | P1/P3 |
| `tracer_temperature_weights` | cap,host,粒子 Cartesian/球面坐标 -> 节点 IDs/weights；沿用 tracer 几何 | P1 |
| `pices_advance` | 已接受状态和 dt -> 完整接受状态或明确终止；唯一移动/时钟入口 | P1 |
| `pices_apply_assimilation_increment` | 完整步实际网格 deltaT -> 粒子增量及映射残差；不重新计算 TA | P4 |
| `pices_restore` | 已恢复 legacy tracer 属性 + 版本元数据 -> 分配/重建缓存；不初始化/覆盖 Tp | P2 |

PG 必须保留原来的算式、求和顺序和默认路径。共享抽取分小步进行，并以 k 非恒定及非零相变的生产核测试保护。不能通过清零 VV 关闭网格平流，那会误删绝热及径向穿相界源。

现有 EBA 热方程按其当前符号整理为：

```text
Ceff = rho * (Cp + Tabs * sum(ds * dX/dT))
Sphase_pressure = -rho * Tabs * sum(ds * dX/dr_pressure) * ur
Ceff * DT/Dt = div(k grad T) + rho Q - Hadi + Hvis + Sphase_pressure
```

P3 的 PICES 质量矩阵用 GP Ceff 的 lumped 装配，逐 GP 验证正且有限，冻结系数于热子步起点并在下子步更新。不再乘旧 heating_latent。PG 原 TMass 的节点 rho*Cp 加权近似保持原样；两种离散质量不同，不能声称非恒定参数下逐位一致，应对照同一制造解做收敛验证。常 rho/Cp、无相变、u=0 时作为直接残差对照。

## 6. 同化、CBF 与能量权重

P4 保留原 TA 的目标/深度支持/指数松弛，完整接受 dt 只执行一次。保存实际网格 deltaT_assim 后 Q 到粒子；另记录 `P(Q(deltaT_assim))-deltaT_assim`。该量不收敛时必须升级投影，不重复做 TA。

网格传热账本采用同一 GP/质量离散；物理体积分只数一次共享节点/边界面。P1 无代表体积/质量时，只可报告粒子温度和采样误差，不能把粒子温度总和当作物理能量。P2/P3 若引入代表体积，需版本化并验证迁移及重采样守恒。

P3/P4 CBF 用 PICES 网格热阶段相同残差。若需 Eulerian 形式，显式由同状态材料导数减 u.gradT；不得把 Tdot_heat 直接送入旧残差后又加 u.gradT。当前 PG 日志保持 schema=2；PICES 增加方法/导数语义标记后再开放此诊断。

账本分开列：网格求解储热差、物理源/边界、TA、粒子重映射差额、补粒子/合并差额。归一化权重与子网格扣除不构成能量守恒证明。禁止通过全域温度平移或 Lenardic filter 让预算看似通过。

## 7. 检查点格式与转换（P2 实现）

保持 PG legacy general header 不变，不在其后直接插入无版本字段。PICES 的 Tp 作为已注册 extraq 由 legacy tracer 段写出；每 rank 同步写强制伴随 metadata 文件 `*.pices.json`。

metadata v1 至少包含 magic=`CITCOMS_EBA_PICES`、schema=1、solver commit、phase=`accepted`、step/time/dt、MPI 分解与 cap IDs、粒子数量、extraq 名称/槽位表、Tp 归一化、Le/插值方法版本、基准参数摘要、主 chkpt 文件 SHA256。JSON 采用明确字段名，浮点按足够精度序列化；读取时检查主文件时钟与 metadata 精度一致。所有 rank 文件先以临时名写，所有 rank 成功后再发布完成 manifest；残缺集合拒绝恢复。

恢复先确认 cfg/metadata、注册属性、分配，再读 legacy 文件；随后重建几何/host、tau 和临时缓存。不能走 fresh N→P 初始化。旧 EBA -> PICES 用显式转换模式/独立工具从网格初始化 Tp，记录 conversion；这是新数值实验，不是等价续算。无伴随文件的 PICES restart、属性数/顺序不符、checksum 不符或跨分区恢复均拒绝。Pg 读写仍完全沿用原格式。

## 8. 审计风险到契约的对应

| 审计项 | 本契约处理 |
|---|---|
| F01 索引 | 第 2 节：1-based 粒子及独立属性槽注册 |
| F02 几何 | 第 3 节：Cartesian 物理长度，固定定义 |
| F03 启动 tau | 第 3 节：首次/每子步重算，k=0 明确 |
| F04 时间层滞后 | 第 4 节：先移动/投影再热求解及 Stokes |
| F05 EBA 物理丢失 | 第 5 节：共享物性/源核及 Ceff；阶段未支持则拒绝 |
| F06 TA 不一致 | 第 6 节：接受 dt 一次同化及独立增量 |
| F07 重启 | 第 7 节：主文件不插字段，版本化伴随元数据 |
| F08 MPI/空覆盖 | 第 3 节：交换后归一、无覆盖显式失败 |
| F09 非守恒 | 第 6 节：物理权重与重映射差额分开，不作保证 |
| F10 子网格抵消 | 第 3 节：保留原始粒子子网格增量 |

## 9. P0 验收与停止点

运行说明见 runs 仓库 `PICES_P0_HPC.md`。P0 使用现有 C 主程序的 PG 路径：12 ranks、每 cap 5³ 节点、两步，Q0=0.3，Di=0.1，恒定 k，phase/TA/tracer 关闭，输出温度/速度、CBF、热源日志和检查点。无需 production age/velocity/refstate 数据。

这个小案例验证部署、PG、热源/CBF及输出完整性，不是 PICES、Pyre 或生产 EBA 全功能验收。特别不能从其无 tracer 推论迁移已经验证。P0 生产源码与基线相同的 Git 比较是 PG 未变的主要证据。

验收要求：记录 solver/runs commit、build log、binary/input hashes；12 rank 的 step 0/1/2 节点文件齐全且有限、径向边界温度精确；各 rank 两次 thermal_exit，无非法温度；时钟为 float 存储的 1e-6/2e-6；step 0 和 2 检查点齐全；native CBF 积分/极值与日志相符，Qinternal SI 标度和 Qtotal 符号和相符。

旧终止函数正常运行也返回 8；不能仅凭该返回码通过或失败。审计脚本接受 0/8，但其余数据全部必须通过。源码未加入“P0通过”的自动升级机制。HPC 回传后人工审计结果，明确通过后才开始 P1。
