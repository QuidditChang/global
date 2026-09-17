# 相变全链路复核：密度、潜热、几何和诊断

日期：2026-09-17。对象：`src/global` 的 `cmbhf_ALA_strict`，HEAD `324ecda15d06bc2f3bd4af50f12782cf97a2e780`，加本次尚未提交的 GP 深度修复；当前 runs 配置和参考文件。

## 结论与证据级别

**确认的是能量调用的几何参数错误，不是共享相变函数整体错误。** 含相变跳跃的参考密度按用户选定模型保留，未将其当作本次缺陷修改。密度、beta、Tref、Gamma、DeltaS、宽度、黏度及求解参数均未改变。

| 检查对象 | 结果 | 证据范围 |
|---|---|---|
| `rtf[3]` 实际语义和能量调用 | 原调用错误；已修复 | 生产者/消费者追踪、实际 C 球面单元、旧调用负对照 |
| `q/X` 公式、温度和径向导数 | 通过 | 独立中心差分与实际 C 能量执行 |
| 相变密度异常和浮力符号 | 通过 | 实际 C 节点相变、`apply_one_phase`、当前参数和参考表 |
| `X-Xref` 参考扣除 | 通过 | T=Tref 时实际节点浮力严格为零；热/冷扰动响应正确 |
| 单熵潜热项及三个相变叠加 | 通过 | 完整热残差，分别开启三项与同时开启相符 |
| 归一化、参数传递、参考文件 | 通过当前有效输入检查 | 最大归一化误差约 1.00e-15；两种输入路径逐项核对；参考文件哈希未变 |
| 旧一维潜热标定目标 | 可复现 | +60/+43/-34 K；是标定重现，不是三维演化独立验证 |
| `Qphase` 空间平均 | 存在诊断积分偏差 | 算术 GP 平均不等于 Jacobian 加权平均 |
| `Qphase` 与最终温度的时间对应 | 不严格同步 | 值来自最后一次 corrector 之前；其他热源在步末重算 |
| 旧 `.heating` 第三列 | 不能当潜热 | 仍输出重置为 1 的 legacy `heating_latent`；当前 cfg 未启用此可选文件 |
| 全局热演化、时间/网格收敛 | 尚未验证 | 本轮无完整 Python2/Pyre/MPI 模型运行 |

因此不能把结果写成“全部相变物理错误”，也不能写成“全部三维能量闭合已证明”。

## 1. 几何错误的直接证据

真实调用链为：

```text
pg_solver
  -> get_global_shape_fn(..., pressure=0, sphere=1, rtf, ...)
       -> form_rtf_bc(k, Cartesian_GP, rtf, bc)
  -> pg_shape_fn(..., rtf, ...)  [只读 rtf]
  -> element_residual(..., rtf, ...)
       -> phase_change_state(..., depth, ...)
```

`Size_does_matter.c:210` 明确写入 `rtf[3]=1/sqrt(x*x+y*y+z*z)`。这符合角向梯度所需的 1/r。`pg_shape_fn` 不会把它改回 r。原热残差用 `ro-rtf[3]`，实际是 `ro-1/r`，现改为 `ro-1/rtf[3]`。

节点相变的 `sx[...,3]` 存储 r，节点调用 `ro-sx[...,3]` 是正确的。不能改动 `rtf` 生产者，也不能把节点调用再取倒数。

在单位外半径、R0=6371 km 下，410/520/660 km 原先对应 -438.200/-566.214/-736.274 km。原调用把相变函数推向饱和区，严重抑制相变导数；此前节点浮力仍能正常显示相变。

局部 Git 历史中，错误调用出现在 `5bb15b4`（2026-08-17，Finalize strict-ALA thermodynamic closure）。此前运行没有报错不能排除它：tanh 仍返回有限的 [0,1] 数值，不会产生崩溃或非法值。

## 2. 共享相变函数与完整导数

约定：z=ro-r 向下为正，ur 向外为正，X 是高压相分数，DeltaS 是高压减低压熵差。

```text
q = (z-z0) rho_ref*g_ref - Gamma (T-T0)
X = (1+tanh(q/w))/2
L = 2 X(1-X)/w
dX/dT = -L Gamma
dX/dr |_T = L [-rho_ref*g_ref + (z-z0) d(rho_ref*g_ref)/dr]
DX/Dt = X_T (partial_t T + u.grad T) + X_r|_T ur
```

对当前固定径向参考系数，这是所选 reduced-pressure proxy 的链式求导。`phase_change_state` 正确实现这些导数；热核的 `material_dT` 已含时间导数和平流，不能再次加一遍温度平流。GP 使用当前传入的 thermal-corrector `field` 插值得到的 T，而不是旧 `phase_B` 缓存。

GP 插值采用 `rho_g=sum(N*rho_node*g_node)`，同时用 GNx 的径向导数求 `d(rho_g)/dr`。`phase_change_state` 的参数传法 `rho=1, gravity=rho_g` 是传入同一个乘积，不是漏掉密度。热容量前面的 `rho_gp` 则单独由节点密度插值得到。

这仍是局部 `(z-z0)rho*g` 压力代理，不是完整静水压积分。三维 Q1 单元内参考量的几何插值误差与真实锐界面问题，均不能由此局部函数测试消除；本轮不更改这些模型约定。

## 3. 密度变化、参考扣除、重力及机械功

当前动态相变密度异常为

```text
delta_rho_phase = sum Delta_rho_i (X_i-Xref_i)
Xref_i = X_i(z,Tref(z))
```

`apply_one_phase` 对 buoy 加 `-phase_Ra*(X-Xref)`，`get_buoyancy` 随后统一乘 `g_ref/g0` 并去水平平均。`phase_Ra=Delta_rho*g0*R0^3/(eta0*kappa0)`，无需再乘一次 rho_ref；密度对比度本身已经是 kg/m3。

因此 X 增加意味着更密、更向下的浮力；正 Gamma 的相变在同深度变热时 X 减少，负 Gamma 则相反。实际 C 节点测试确认这两种响应。T=Tref 时 X=Xref，动态相变异常为零，没有将整个静态 Delta_rho 再加一遍。

`Profile_output.c` 和 `Drive_solvers.c` 的相变浮力/机械功分解使用同一 `X-Xref`，各分量去水平平均后乘径向 g；在 g 仅随半径变化的当前参考态下，顺序与总浮力一致。`Solver/Solver.py::advance` 先推进热/示踪，再解速度；组装力时刷新节点 X。没有发现这条正常路径把旧节点相态用于下一次动量方程的错误。

连续性和压力浮力继续消费已有 beta；本次没有加入相变导数、改变 H 或修改任何参考密度跳跃。rheol7 也没有新增使用 `phase_B` 的路径，其已有深度/Tref 黏度分支保持原样。

## 4. 潜热、热容量、符号与重复计算

现有 reduced entropy 模型的相变项在热残差左侧：

```text
rho Cp DT/Dt - alpha T_abs DP/Dt
    + rho T_abs sum(DeltaS_i DX_i/Dt) = conduction + other heating
```

故右侧净热源是 **-Qphase**。下沉形成高压相时：410/520 的 DeltaS<0，释放热；660 的 DeltaS>0，吸收热。反向相变反转符号。输出中 Qphase>0 表示左侧吸热贡献，不能直接读成“正加热”。

`phase_delta_s=DeltaS/Cp0`；温度使用 `T+surface_temp=T_abs/DeltaT`。相变项在归一化后无需额外乘 Di、phase_Ra 或局部 Cp。热容量修正可写作 `Cp + T_abs*sum(DeltaS*X_T)`；当前 Gamma 与 DeltaS 的符号使额外容量非负（正绝对温度、正宽度），没有从该项产生负热容量的问题。

参考 Cp 由源表分支拟合并归一化，生成器没有向其加入 DeltaS*dX/dT 潜热峰。旧 `D_phi/H_phi`、逆 D 缩放扩散和基于旧 phase_B 的潜热路径已退出当前热残差。三个相变共享一次循环中的单一 DeltaS 项；完整 C 单元逐项开启测试确认没有重复叠加。含相变偏移的 Tref 是背景/初值，不是每个时间步另加的热源。

原 Stage 7 已明确：Gamma=[4,6,-2] MPa/K 与 DeltaS 标定是约定的 reduced model，并非要求由 Katsura Delta_rho 重构完整 Gibbs EOS。本轮保持这个选择，不因两组参数独立就擅自重新标定。

## 5. 归一化、输入和参考标定复核

| 量 | 生产归一化 |
|---|---|
| z0 | z0_dim/R0 |
| Gamma | Gamma_dim*DeltaT/(rho0*g0*R0) |
| w | w_dim/(rho0*g0*R0)，不是一个固定几何厚度 |
| T0 | (Tref(z0)-Ttop)/DeltaT |
| DeltaS | DeltaS_dim/Cp0 |
| phase_Ra | Delta_rho*g0*R0^3/(eta0*kappa0) |

读取当前文件后独立重算，五个 cfg 向量的最大相对差约 1.00e-15。C 独立解析和 Pyre `Phase_set_properties` 均向相同的 phase 结构传递 depth/density/entropy/width/Gamma/T0，width 在 C 中取倒数，phase_Ra 在启动验证时由密度差重新计算。当前输入有限、width>0、熵差非零，未发现当前有效参数被漏传。

runtime 参数使用 float，因此中心/导数在真实 C 中存在约单精度参数误差；不能用 Python 双精度中心精确命中作为 C 零误差声明。

参考文件哈希仍为：

```text
refstate_ALA_strict.txt  36a43a62688644b23040287269d0c80d7cbc3d2db5e3c5d0a35637304d992850
interval_ALA_strict.txt   c7bbe4347d235c36ce5f1fbdf55647e1de2c2d6878a741810e71f90567e0c255
```

原一维参考轨迹积分仍给出 60、43、-34 K。这证明标定未被此次修改破坏，但由于 DeltaS 原本就是按这个积分标定的，它不是独立的三维时步验收。

## 6. 新确认的输出限制

### 空间平均

`element_residual` 当前保存 `heating_phase=sum(Q_gp)/8`，而它组装到节点残差时使用 `Q_gp*Jacobian_gp`。球面六面体的 Jacobian 通常随 GP 变化。因此用该算术均值再乘单元体积，不严格等于热残差中的相变积分。

在本轮采用当前径向单元、角向半宽 0.012/0.014 rad、T 偏移 0.01、Tdot=0.03、ur=0.2 的三个隔离测试单元中，算术均值相对体积加权均值的差为 **0.0346%、0.0691%、0.2069%**。这是特定测试场的偏差，不是全局误差界。真正的节点残差和独立相变积分相对差约 1.3–1.5e-9（含有限差分误差）。

### 时间对应

每个 thermal corrector 的顺序是：计算热源和残差/缓存 Qphase → 更新 T/Tdot → 温度边界处理。循环结束后还可能过滤温度、同化，再重算普通绝热/黏性热，最后打印预算；`process_heating` 不重新计算 Qphase。

所以打印 Qphase 是最后一个残差评估状态的值，不能声称它来自最终 T/Tdot，也不能把这个混合时刻的 Qtotal 直接用来证明一步积分能量守恒。是否重算步末瞬时值，还是保存与时间离散一致的积分，应先确定输出含义；本轮不悄悄改变该语义。

### legacy 输出

`Output.c::output_heating` 的第三列仍是 `heating_latent`，该字段被重置为 1，已不进入热残差。正确相变量来自 profile `qphase` / THERMAL_BUDGET `Qphase`；`qphase_adi` 是兼容别名。当前 cfg 的 `output_optional=comp_nd,surf,botm,k` 未开启旧 `.heating`，但分析历史文件时必须避免误读。

## 7. 为什么“以前都通过”

旧测试的强项是方程符号、独立 Python 参考轨迹、参数标定与源码结构。`test_active_operator_uses_current_state_and_complete_derivative` 只检查调用存在和几个字符串；一维积分直接使用正确的 depth，没有运行 `get_global_shape_fn` 的 rtf。C 语法检查也不可能判断一个 double 是否代表 r 或 1/r。

本轮新增的完整单元测试执行真实 `construct_shape_functions`、`get_global_shape_fn`、`form_rtf_bc`、`pg_shape_fn`、`element_thermal_transport`、`element_residual` 和节点相变/浮力函数；仅将无关导热系数固定，关闭边界通量，不使用完整 MPI mesh/solver。

- 修复版本：完整几何/热残差、参考扣除、三个相变叠加测试通过。
- 在内存中恢复旧调用：12 个完整单元能量案例全部失败，而节点密度与参考扣除测试仍通过。
- 前一轮 36 个局部材料导数案例同样能够捕获旧错误。
- 完整 strict 套件 342 项：341 通过；仅已有的 rheol7/cold_scale 配置对比测试失败，与本轮相变改动无关。随后增加真实 PG 形函数调用后，相关单元测试再次通过。
- 两个相关 C 翻译单元语法检查通过，保留既有旧式声明等警告。

新增测试不等于完整时间积分：默认 thermal corrector 为固定两次，TMass 是基础 rho*Cp lumped mass，不包含状态相关潜热容量；后者通过残差迭代进入。该安排并不自动错误，但修复恢复潜热后，需要实际时步/校正次数收敛和能量审查。

## 8. 交付与后续边界

量化结果见同目录 `phase_full_audit_2026-09-17.json`。可复跑：

```bash
python3 -m unittest discover -s tests/ala_strict -p 'test_phase_*.py' -v
```

本轮没有进一步改变生产物理或诊断实现，没有改参考态或参数，没有提交 HPC 作业、commit 或 push。唯一生产补丁仍是前一轮的 GP 半径倒数修复。

下一次真实热演化验收应比较同一初始/重启状态下的 dt、dt/2 和校正次数，分别报告相变层温度、X、相变积分及敏感性；不能把修复前后产生不同温度轨迹当作补丁失败，也不能用本轮局部通过代替全局物理验证。

外部方法交叉核对：ASPECT 的 [LatentHeat 材料文档](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1MaterialModel_1_1LatentHeat.html) 提供相变函数/熵导数方法背景；本报告的具体通过与缺陷结论均以本地 CitcomS 代码和执行证据为准。
