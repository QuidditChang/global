# P4：接受步 TA、材料热残差 CBF 和重启

2026-10-01。P3 HPC 作业 12249416 审计通过后进入本阶段。
`pices_p4=on` 是显式开关，依赖 `energy_solver=pices`、`pices_eba=on`。
本地验证记录见 `P4_LOCAL_RESULT.json`。2026-10-01，HPC 作业 12252234
通过阶段验收；完整报告位于 runs 的 `PICES_P4_HPC_AUDIT_12252234.md`。

## TA 编排

移动粒子一次、投影一次、完成全部热子步后，才调用原有
`assimilate_lith_conform_bcs` 的 relaxed TA 分支一次。强制
`lith_age_asml=1`、`lith_age_time=1`、`temperature_bound_adj=0`。
目标仍是 `Tref - Tref_surface*erfc(depth/(2*sqrt(age_nd)))`，板龄上限、
有符号海沟距离决定的支持深度和指数空间权重均沿用 `Lith_age.c`。
使用完整外步 dt 的原有指数松弛，不按热子步次数乘 TA。

读取板龄及计算目标时临时采用接受端点 `t+dt`；随后恢复旧时钟，最终
发布 Tg/Tp 时统一提交一次。实际增量是边界处理完成后的
`deltaTA = Tg_after - Tg_before`；仅以 `Q(deltaTA)` 更新持久粒子温度。
不会重新由目标计算第二次粒子松弛，也不把映射残差反馈到 Tg。

`PICES_TA` 逐步报告调用次数、接受时间、dt、最大实际温度增量、
`max|P(Q(deltaTA))-deltaTA|`（含边界节点）及独立 TA 储热项。
储热采用 TA 前冻结的有效容量积分，是离散诊断，不是完整非线性焓差。
`PICES_EBA` 保持纯热子步账本；`PICES_STEP remap_energy` 单独报告重映射。
未引入物理粒子质量，因此这些量不能证明粒子物理能量守恒。

standalone 参数入口补齐已有 Pyre 的 `flag_depth_file`、
`flag_depth_new_file`、`tf_file` 字符串，使用现有强迫文件读取器。

## CBF 的精确定义

热子步 s 的单元节点残差采用该子步实际使用的源、扩散和 lumped 容量：

```
Rbar_e,a = sum_s [dt_s * (F_e,a^s - K_e,ab^s * T_b^s)
                 - M_e,a^s * (T_a^(s+1) - T_a^s)] / dt_outer
```

在热阶段内缓存并累积，不在输出时用最终 Tg 或新 Stokes 速度重算。
该量经过原有 Q1/GLL 边界面组装和面积归一后得到通量；顶部正方向为
地幔向地表，底部正方向为地核向地幔。独立面积分应满足

```
Q_bottom_W - Q_surface_W
  = PICES_EBA.boundary_reaction / dt * k0 * deltaT_K * radius_m
```

这是 **TA 前、完整热区间平均的材料热残差通量**。粒子已承担平流，
不再添加 Eulerian `u.gradT`，不包含 TA 或粒子重映射的能量项。
文件明确写入 `state=PICES_heat_stage_average_before_TA`、
`derivative=material_heat_only advection=particles TA=excluded`。
`CBF_use_advection` 是 PG 残差选项，在 P4 缓存分支不使用；P4 cfg 设为 off。
`write_q_files` 的旧输出仍拒绝。PG CBF 算式和默认行为不变。

新算例第 0 步没有先前热区间，记录 `unavailable_initial`，不制造 CBF
文件或 `CBF_NATIVE_COMPLETE`。恢复第 2 步时缓存有效，可立即复现通量。

## 检查点 schema 2

P4 的 legacy 主二进制布局不变。JSON `schema=2`，其 `.pices.state`
按 native ABI 顺序为：三个 nno 长度 float 速度数组、一个 int
`cbf_valid`、nel*8 个 double 单元节点 `Rbar`。原有
`accepted_velocity_sha256` 字段在 schema 2 中覆盖整个伴随文件；长度、
有限性、有效标志及集合 manifest 都要通过。manifest 自身仍为 schema 1。
P2/P3 继续使用原 schema 1。

物理指纹加入 TA 参数、Tref、当前时刻实际读取的板龄与几何数据；只把
真实径向边界 TB 纳入边界指纹，排除动态内部 TA 目标。恢复时刷新同一
接受时刻的强迫数据，重建热系数，保留 Tp、CBF 缓存，不施加 TA。
未来未读取时刻的文件不属于检查点指纹；HPC 归档另保存完整输入 SHA256，
审计还比较各段实际输入。只支持相同可执行提交、ABI、分区的续算。
必须新建 P4 算例，不能把 P3 检查点当作 P4 检查点。

## 本地验证与边界

- 12-rank 真实 Stokes：连续 4 步、独立 2 步、恢复 2→4；检查
  216 份解码场输出、120 份 CBF 文件及有效二进制状态逐位一致。
- 强源制造解触发 3 个热子步，TA 恰好一次；独立指数公式核对网格增量，
  用输出粒子权重核对 Q，按共享 Cartesian 节点汇总独立核对 P(Q) 诊断。
  网格 TA 与粒子 Q 误差分别为 5.53e-17、6.25e-17；独立映射诊断误差为 0。
  TA 开/关两算例的 CBF 文件逐位一致，边界收支用面权重独立积分。
- 合成 HSC 层、连续初始 T 的 5/9/17/25 节点加密测试。
  两个最粗网格尚未进入误差下降区间；不能据此宣称所有网格单调收敛。
  最大映射误差依次为 0.0158121 / 0.0195817 / 0.0136692 / 0.00991013；
  17→25 节点继续下降约 27.5%。
  测试使用合成 10000 Ma 板龄来展宽可解析的 HSC 层，只检验数值一致性，
  不是地球物理参数选择或收敛阶证明。实际板龄的边界层仍须用足够细网格。
- 5 个重启拒绝测试：TA 时间尺度、指数、Tref、板龄变化和缓存损坏。
- P3 完整回归通过；180 份输出和检查点有效状态与改造前逐位相同。
  11 项 PICES 内核测试和 6 项 PG CBF 内核测试通过。

阶段范围仍为均匀 Newtonian 黏度、无组分反馈、固定分区、无补粒子。
Pyre 参数已接线，但本阶段运行验证使用 standalone CitcomSFull；没有宣称
在本地验证 Python 2.6/Pyre 安装或真实板块重建数据。HPC 作业 12252234 已完成三段算例并通过审计；本次 HPC 每步只有一个热子步，
多子步证据仍来自本地制造解。P5 可行性实验见 `P5_IMPLEMENTATION.md`。

## 重复验证

```bash
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py --build-dir /tmp/p4-build
python3 tests/pices/build_manufactured.py /tmp/p4-build
python3 tests/pices/run_p2_smoke.py /tmp/p4-build RUNS_PATH --stage P4 --output /tmp/p4-smoke
python3 tests/pices/run_p4_local.py /tmp/p4-build RUNS_PATH/cmbhf_EBA_PICES_P4.cfg --output /tmp/p4-local
python3 tests/pices/run_p4_guards.py /tmp/p4-build /tmp/p4-smoke --output /tmp/p4-guards
```
