# CBF GLL Q1 实现与验证报告

日期：2026-09-19。分支：solver/runs 均为 `cmbhf_EBA`。

## 当前交付状态

CBF 数值计算、双边界运行时 GRD、频率继承及回归测试已实现。
用户授权后，已修复 `Parsing.c` 两处未初始化指针写入，并正确跳过分隔逗号解析上下限。最终构建直接使用仓库源码，不再使用临时解析器补丁；12/24 ranks 均完成一步推进、输出 step 0/1 双边界 GRD，8 个文件均通过尺寸、有限值、接缝和极点检查。完整生产 HPC 作业仍未执行。

方法以 `APPENDIX_C_AUDIT_PLAN.md` 为准：物理能量 Galerkin RHS、与 Q1 温度元共置的 GLL 面质量、外向导热通量、CMB 写出时翻号。没有采用旧设计中的 PCG/完整边界矩阵要求。

## 实际改动

- `Advection_diffusion.c` 新增只读诊断包装，调用原 `element_residual`，选择普通 N 测试函数；共享球坐标梯度、rho/Cp、可变 k、全部显式相变能量项和源项。
- `Process_buoyancy.c` 替换旧镜像残差，按实际 Q1 面四角 GLL 面积元组装；所有 rank 参与现有双精度节点交换，取消危险的径向提前返回。
- 输出时重新计算当前 T/u 对应的 adi/visc 到临时数组，退出恢复原指针、phase 诊断值和应变几何工作缓存。
- `cbf_geometry.h` 提供三维面 Jacobian 与三维射线/Q1 逆映射，避免经纬平面插值的接缝/极点问题。
- `cbf_output.h` 原生 C NetCDF 输出 `PostProc/HF_CBF/cmbhf_CBF_<step>.grd` 和 `eshf_CBF_<step>.grd`。无需 Python 后处理。
- `Controller.py` 的 `monitoringFrequency_cmbhf_CBF=-1` 继承普通监测频率；0 禁用；正数独立周期。主 EBA cfg 已设为 -1，因此每 50 步输出。
- standalone `bin/Citcom.c` 支持 `cmbhf_CBF_freq` 同样语义，便于本地独立验证；其默认 0 保持兼容。
- configure 检查 NetCDF C 库。库不可用而启用 CBF 时明确失败，不能悄悄改为 gzip。

## 已完成测试

| 验证 | 实际结果 |
|---|---|
| 参数控制字符串回归 | 1/1 通过，覆盖 8 种默认值/上下限组合 |
| 生产核/几何单测 | 5/5 通过 |
| 原有生产相变几何回归 | 3/3 通过 |
| 合成闭合六面网格 NetCDF/MPI 回归 | 1/1 通过，内部运行 1 和 2 ranks |
| 全部 standalone 源码编译及链接 | 通过；使用 legacy C 警告降级编译选项 |
| 最小真实 EBA 时间推进 | 12/24 ranks 均实际写出 step 0 和 1 的双边界 GRD |
| 真实 CitcomS 全球网格，解析径向导热 | 12 ranks 与 24 ranks (nprocz=2) 的两个 GRD **逐位相同** |
| CBF 开/关数值侵入性 | 12 ranks 相同 cfg，96 个常规 gzip 数据文件解压内容完全相同 |
| 三档真实网格解析场细化 | 顶/底局地 RMS 与积分误差均随细化降低，见下表 |
| git diff --check | 通过 |

### GLL 及 RHS 专项

单测检查仿射面 A/4、旋转/反转不改变质量、扭曲三维面端点 Jacobian 与普通 Gauss 面积不同；使用非单位半径的球坐标导数，独立验证瞬态、横向对流、导热、源项符号；显式非零相变熵项进入 RHS；诊断不改变 heating_phase。

### 独立 GRD 读取

用 SciPy NetCDF 读取 C 写出的文件，检查 361×721、单位/方法元数据、所有点有限、0°/360°完全相等、极点各经度完全相等。
合成网格的解析点值最大误差为顶部 `1.1921e-7`、底部 `2.3841e-7`，属于 float32 存储误差；原生面积分分别为 `72e6` 和 `42e6 W`，精确一致。
这些合成文件不是地幔模拟成果。

### 真实网格解析导热细化

把解析 `T=B(1/r-1/r_top)` 采样到实际 CitcomS 网格，置 u=0、Tdot=0、源项=0、k=4 W/(m K)。顶部半径 6371 km、底部 3504.05 km（测试 cfg 默认 ri=0.55），DeltaT=3400 K。
解析热流总量 `1.330780882082e12 W`；底部/顶部密度分别 `0.00862494035`、`0.00260904445 W/m²`，均为正。

| 每 cap 节点数/轴 | 底部经纬采样相对 RMS | 顶部经纬采样相对 RMS | 底部原生积分相对误差 | 顶部原生积分相对误差 |
|---|---:|---:|---:|---:|
| 5 | 3.9510e-2 | 3.1735e-2 | 1.2340e-2 | 5.5101e-3 |
| 9 | 1.2073e-2 | 1.0142e-2 | 3.3339e-3 | 1.2881e-3 |
| 17 | 3.8099e-3 | 3.3890e-3 | 8.7071e-4 | 3.1275e-4 |

RMS 是规则经纬点等权统计，不是假称面积加权 FE 范数。此测试是**解析场采样后的提取/几何验证**，不是求解制造 PDE，也不能视为复现论文 Q2 的基准收敛阶。

## 必须保留的解释边界

1. 实际生产重建的板块同化与温度过滤可在时间推进后修改 T，而 Tdot 没有同步包含这些外部增量。文件记录 `output T + solver Tdot` 及相关开关；不得称之为封闭系统的精确离散能量反力。若需求升级为含同化的完整时间离散能量收支，需要独立定义状态/增量归属。
2. step 0 使用初始化 Tdot（通常零），在元数据中明确标为初始诊断；不是已经接受的瞬态时间步。
3. GRD 是 Q1 原生场的点采样，不是保守重映射。总功率应读取原生 GLL 积分元数据，不用经纬网格算术均值代替。
4. 当前保证全球球壳、上下定温边界；不支持悄悄把混合或定通量边界当定温恢复。缺失覆盖/接缝不一致/非法面积均失败。
5. 演化小模型在 12/24 ranks 上得到的流场和 GRD 并非逐位一致，不能把不同解的差别称为 CBF 装配误差；冻结相同解析场的 CBF 提取已验证逐位一致。
6. 小模型是功能 smoke test，不是生产物理校准。未执行用户完整重建数据集的 HPC 作业或 Blankenbach/King 全套基准。
7. 原有 `parallel_process_termination()` 即使正常结束也返回 8。最小 evolving 测试须检查已完成步、GRD 和日志，而不能把退出 8 单独当作 CBF 失败；manufactured 驱动显式返回 0。

## 复现入口

- `tests/cbf/test_cbf_kernel.py`
- `tests/cbf/test_grid_output.py`
- `tests/cbf/build_validation.py --build-dir /tmp/cbf-validation`
- `tests/cbf/smoke.cfg`、`refstate.txt`
- `tests/cbf/manufactured_state.h`
- `tests/cbf/verify_real_outputs.py`

测试产物分别位于 `/tmp/cbf-smoke-{12,24}`、`/tmp/cbf-smoke-off`、`/tmp/cbf-manufactured-{12,24}`、`/tmp/cbf-refine-{9,17}`。这些目录可能被系统清理；关键误差和结论保存在本报告中。

## 共享解析器修复后的最终复验

`lib/Parsing.c` 将两处写入未初始化指针的表达式改为 `strchr` 指针赋值，检查空指针后越过逗号读取数值。不改变调用接口。回归入口为 `tests/cbf/test_parser_control.py`。

最终源码构建位于 `/tmp/cbf-final-build`，记录 `parser_workaround=False`。最终 12/24 ranks 日志与文件位于 `/tmp/cbf-final-smoke-{12,24}`，两者日志均为 `cycles=1`。正常结束沿用原有退出码 8。此前细化、冻结场分区一致性和开关对照使用的是等价临时解析器修复构建；这次最终复验验证仓库正式修复的启动及实际时间推进，没有重复宣称完整 HPC 物理验证。
