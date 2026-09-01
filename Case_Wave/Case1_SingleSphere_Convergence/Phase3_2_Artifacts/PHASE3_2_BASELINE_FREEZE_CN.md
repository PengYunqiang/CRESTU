# CRESTU-1F Phase 3.2 基线冻结与来源审计

生成日期：2026-09-01（Asia/Shanghai）

```ini
GIT_REPO_ROOT = E:\GitWorktrees\CRESTU-Phase3\_2
BRANCH = phase3-2-mesh-convergence
INITIAL_HEAD = 7f2d391d1b42366fdb04e1b493439dab438584eb
INITIAL_WORKTREE_STATUS = CLEAN
AUDIT_RUNTIME_HEAD_RELATION = CHECKPOINT_OR_DESCENDANT
MATLAB_PWD = E:\GitWorktrees\CRESTU-Phase3\_2
PHASE3_2_RUNTIME_PROVENANCE = PASS
PRODUCTION_PHYSICS_CHANGED = NO
WAMIT_REFERENCE_PATH_RECOVERED = YES
WAMIT_RAW_REDISTRIBUTION = UNCONFIRMED_DO_NOT_COMMIT
```

## 冻结边界

冻结 checkpoint 为 `7f2d391d1b42366fdb04e1b493439dab438584eb`。本审计没有修改 production physics、mesh topology 或 WAMIT 数值。

文本输入同时记录 raw-byte SHA-256、canonical-LF SHA-256 和 Git blob ID；没有统一换行、renormalize 或修改 `.gitattributes`。清单共 54 个已跟踪文件。

## MATLAB runtime path

要求的 9 个 production function 均精确解析到当前 worktree；没有路径解析到原仓库。

## WAMIT provenance

目标 worktree 推导路径 `E:\GitWorktrees\WAMIT` 不存在。只读调用链恢复出 Phase 3.1 历史外部目录 `E:\CRESTU-1F_v1.0_20260824\WAMIT`；6 个关键文件已记录 hash。它们不属于 Git，公开再分发权限未确认，禁止复制或提交。

## 可信度边界

本文件只证明 checkpoint/source/runtime provenance。它不证明 Phase 3.2 mesh convergence，不把 WAMIT 作为调参目标，也不提升 `PHYSICAL_VALIDATION=PARTIAL` 或 `READY_FOR_PRODUCTION=NO`。
