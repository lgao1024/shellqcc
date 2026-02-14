# shellqcc

[English](README.md) | [Chinese](README.zh-CN.md)

用于 Quantum ESPRESSO 与 ORCA 工作流的 Shell 工具集。

## 项目简介

本仓库提供了一组 Bash 脚本，用于简化常见计算化学流程：

- 生成与转换 Quantum ESPRESSO 输入文件（`qei.sh`）
- 执行 QE 工作流并自动生成后处理输入（`qec.sh`）
- 批量准备与运行 ORCA 任务（`orcai.sh`）
- 从 ORCA 频率计算输出提取 Gibbs 自由能（`orcac.sh`）

此外，仓库还包含 `tools/` 下的辅助脚本与 `examples/` 示例文件。

## 目录结构

- `qei.sh`: QE 输入生成 / 转换 / 输出分析工具
- `qec.sh`: QE 任务流程脚本（`pw`、`dos`、`bands`、`w90`、`phonon` 等）
- `orcai.sh`: ORCA 输入准备与运行助手
- `orcac.sh`: ORCA Gibbs 自由能提取脚本
- `examples/`: 分子与晶体结构、ORCA 模板示例
- `tools/`: 环境设置、输入生成、CIF 处理等辅助工具

## 依赖要求

### 系统依赖

- Bash
- `awk`、`sed`、`grep`、`bc`
- `mpirun`（用于 QE 相关任务）

### `qei.sh` 依赖

- `cif2cell`
- ASE 命令行工具（`ase`）
- 可访问赝势库（脚本内路径可配置）

### `qec.sh` 依赖

- Quantum ESPRESSO 可执行文件（按任务需要：`pw.x`、`dos.x`、`projwfc.x`、`pp.x`、`bands.x`、`ph.x`、`q2r.x`、`matdyn.x` 等）
- 可选：`gnuplot`、Wannier90（`wannier90.x`、`pw2wannier90.x`）、`open_grid.x`、`sumpdos.x`、`thermo_pw.x`

### `orcai.sh` 依赖

- ORCA（脚本内配置二进制路径）
- ASE 命令行工具（在 xtb 模式下用于 `.mol` 到 `.xyz` 转换）

## 快速开始

```bash
# 1) 从 CIF 生成 QE 输入
./qei.sh structure.cif qe scf

# 2) 运行 QE 默认任务（pw）
./qec.sh structure.scf.in

# 3) 串行执行多个 QE 任务
./qec.sh structure.relax.in pw dos bands

# 4) 基于模板准备 ORCA 任务
./orcai.sh xtb.tmp

# 5) 提取 ORCA 频率输出中的 Gibbs 自由能
./orcac.sh
```

## 脚本用法

### `qei.sh`

```bash
# CIF -> 输入（默认: qe scf）
./qei.sh file.cif [qe|vasp|cp2k] [calc_type]

# QE 输入 -> cif/conv
./qei.sh file.in [cif|conv]

# QE 输出 -> 分析
./qei.sh file.out [scf|relax|vc-relax]
```

说明：
- 当目标 `.cif` 文件不存在时，`qei.sh` 可进入交互式 CIF 生成向导。
- 脚本会交互询问赝势选择及相关参数。

### `qec.sh`

```bash
./qec.sh input.in [tasks...]
```

可用任务：
- `pw`（默认）、`dos`、`co`、`esp`、`hp`、`bands`、`w90`、`phonon`、`thermo`、`dry`、`all`

示例：

```bash
./qec.sh prefix.relax.in pw dos co esp bands
```

### `orcai.sh`

```bash
./orcai.sh <template_file>
```

行为：
- 使用 `xtb.tmp` 时：如有需要先将 `.mol` 转换为 `.xyz`，然后生成 ORCA 输入、创建任务目录、运行 ORCA，并将生成的 `.xtb.xyz` 拷回上级目录。
- 使用其他模板时：基于当前目录 `.xyz` 文件生成 ORCA 输入文件。

### `orcac.sh`

```bash
./orcac.sh
```

会扫描 `*/*.out` 中的 `Final Gibbs free energy`，并输出为 `Gibbs.txt`。

## 示例数据

- 晶体 CIF：`examples/crystal/`
- 分子结构：`examples/molecular/`
- ORCA 模板：`examples/orcai_tmp/`

## 许可证

本项目使用 MIT 许可证，详见 `LICENSE`。
