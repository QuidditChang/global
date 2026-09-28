# P0 本地交付与验证记录

2026-09-28。源码分支已创建为 `cmbhf_EBA_PICES`，基线 `3f7f4de`。
**本地准备完成，等待用户 HPC 测试及回传审计；未进入 P1。**

## 交付

- `P0_CONTRACT.md`：冻结索引、映射、Le/tau、子网格公式、先平流后热子步、共享 API、相变/TA/CBF 边界、检查点方案及阶段拒绝条件。
- `tests/pices/verify_p0.py`：HPC 下载结果的只读验收入口；输出 PASS/FAIL JSON。
- `tests/pices/test_verify_p0.py`：正常、NaN、非有限速度、缺行、错时钟、错上下边界共 7 项。
- runs 的 `cmbhf_EBA` 分支新增 standalone cfg、lsf、5 行 refstate、HPC/commit/push 说明。

## 已执行

1. 用现有 `tests/cbf/build_validation.py` 编译全部生产 C 库与原 standalone 驱动，成功。该构建没有替换 parser；运行的是 CitcomSFull，不是制造解驱动。不是 HPC Intel/Pyre 安装验证。
2. `/usr/local/bin/mpiexec --oversubscribe -n 12` 运行 P0 cfg，完成 step 0/1/2；正常 legacy exit=8。12 rank 全部温度/速度文件有限，上下边界为 300/3700 K；所有输出温度范围 300–3700 K；两步时间约 1e-6/2e-6。
3. 全部 3 个时刻、上下两边界的 native CBF 面积分/极值与日志相符；Qinternal=0.3 的 SI 标度及 Qtotal 符号求和通过；每 rank step 0/2 检查点 header 通过。
4. 相变几何 3 项、TA 3 项、CBF kernel 6 项，以及新增输出验收 7 项：**19/19 通过**。
5. 在临时副本中删去一个 rank、截断 step-2 gzip、把 Qtotal 行改成 NaN，三种损坏均被拒绝。未修改真实验收数据。
6. `bash -n` 检查 LSF，通过；Git 比较确认生产源码/构建入口与 EBA 基线无差异。尚未在实际 LSF 调度器执行 launcher。
7. 用明确标记的临时合成 metadata 验证 provenance 分支：字段完整时接受，输入 SHA256 被篡改时拒绝。此项只测试检查器结构，不是 HPC 构建来源验证，合成目录已自动清理。

机器可读摘要：`P0_LOCAL_RESULT.json`。本地 MPI 原始目录在 `/tmp/pices-p0-local-final`，仅为当前会话临时数据；以后以保存的摘要及用户 HPC 归档为审计材料。

## 过程中修正

初始 cfg 将 maxtotstep 设为 2，只跑了一个演化步。确认 legacy total_timesteps 初值为 1 后，将 cfg 改为 `maxstep=2, maxtotstep=3`，在全新目录重新跑通。未改求解器的计数行为。
本机沙盒最初阻止 MPI 打开本地 socket；在获得该本地测试的执行权限后成功运行，未连接 HPC。

## 结果解释与下一关口

本案例为部署与日志/热流离散一致性的 smoke test。小网格 CBF 热流可能有负值和较强时变；当前检查的是其残差、native 积分、单位与输出一致性，不能解释为生产热流已收敛或热预算闭合。没有宣称 Stokes 生产精度、Pyre、PICES 平流、MPI 粒子迁移、完整 EBA 相变/TA 或长期守恒通过。

用户运行 runs 的 `cmbhf_EBA_PICES_P0.lsf`，下载完整 tar.gz 与 LSF out/err 后，使用不带 `--local` 的 verifier 并人工审阅 build/provenance。HPC P0 审计通过后才启动 P1；不自动 commit/push，不自动提交 HPC，也不提前实现 P1。
