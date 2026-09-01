# Phase 3.2 v1 协同网格 family 失败记录

## 结论

```ini
MESH_FAMILY_VERSION = v1
MESH_GEOMETRY_GATE = FAIL
MESH_CONTROL_ISOLATION_GEOMETRY = FAIL
ALGEBRAIC_GATE = NOT_RUN
PRODUCTION_PHYSICS_CHANGED = NO
FAIL_CLASS = MESH_FAMILY_INCONSISTENCY
```

v1 在任何 BEM 组装前完成 18 个 `6 sentinel frequencies × 3 levels` 的 production mesh generation。所有 case 的 runtime/config/BMF 都位于当前 worktree 的唯一临时目录，结束后已清理；没有读取 cache，也没有访问原仓库 production 路径。

## 通过的项目

- body、FS、bottom、FARFIELD panel count 和 component/combined mesh hash 均实际变化；
- 同频率三等级的 1.5R top/bottom outer extent 完全相同；
- physical/source hash 固定；
- duplicate panel count 为 0；
- closed-domain edge multiplicity 全部为 2，nonmanifold edge count 为 0，Euler characteristic 为 2；
- normal orientation、body/FS waterline、FS/wall、bottom/wall continuity 和 analytic component-separation audit 均通过；
- 临时目录残留为 0。

## 原始 FAIL

v1 explicit outer theta controls 56/80/108（bottom 44/64/84）与 body divisions 3/5/7 产生的 legacy ring counts 24/40/56 不匹配。outer annulus count transition 因此产生 sliver triangles；观测到 maximum aspect ratio 122.20、minimum corner angle 4.78 deg。L1 body maximum warp 8.99 deg，poor-panel fraction 5.56%，也超过预声明 8 deg/5% gate。0.5/0.6 rad/s 的部分 median component h 非单调，故不能使用 v1 eta 进行 Richardson/GCI。

完整失败行未删除、未平滑、未隐藏，保存在 `Phase3_2_V1_*.csv`。

## v2 处置

v2 不放宽 threshold，不修改 production source。它把 body divisions 改为 4/6/8，并把每级 FS/bottom/wall theta count 固定为对应 legacy ring 的 32/48/64；同时按 6/9/12 协同细化 FS/bottom radial total 和 3/4/6 wall vertical control。全部 18 个受影响 mesh-only case 必须重跑，v2 通过 geometry gate 后才允许 sentinel BEM。
