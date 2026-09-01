# Phase 3.2A Topology Design Decision v3（冻结）

冻结时间：2026-09-02；在实现 bottom-core schedule 和查看其 endpoint/v5 结果之前。

原始 design、connector-gap v2 design、initial/gap-v2 FAIL evidence 全部保留。gap-v2 的 omega=1.0 唯一失败仍是 bottom median h：`5.439254 → 5.750831 → 4.075027 m`。quality、mean h、connector quality、counts、transition h、topology 和 extent 全 PASS。

Population audit：

- L2 structured 5-block core median h=`0.481 m`，不是 quality blocker；
- L2 overall median 落在 connector population minimum h=`5.750831 m`；
- 继续缩小 connector gap 会消耗 35 deg angle margin，因此不再改变 connector geometry。

## 显式 bottom core radial schedule

新增仅 QUALITY_CONTROLLED_V2 使用的 named control：

```ini
BOTTOM_CORE_RADIAL_LAYERS = 4  # L1
BOTTOM_CORE_RADIAL_LAYERS = 7  # L2
BOTTOM_CORE_RADIAL_LAYERS = 12 # L3
```

`Nx=Ny=N_waterline/4` 保持不变，因此 waterline/core/annulus boundary nodes 与 segment counts 完全不变。只增加四个 structured side blocks 的径向层数 `Nr`；central patch 仍为 `Nx×Ny`，没有 single-center fan。

| level | N waterline | Nx=Ny | Nr requested | side panels `4*Nx*Nr` | central panels `Nx*Ny` | total core panels |
|---|---:|---:|---:|---:|---:|---:|
| L1 | 32 | 8 | 4 | 128 | 64 | 192 |
| L2 | 40 | 10 | 7 | 280 | 100 | 380 |
| L3 | 56 | 14 | 12 | 672 | 196 | 868 |

该 schedule 明确、单调、满足边界 integer contract，不改变 outer annulus、connector、wall、body 或 physics。缺省/LEGACY mode 不读取该 control，继续使用原 `Nr=max(3,round((Nx+Ny)*0.25))`。

必须重跑 legacy 5/5 exact、六个 endpoint 和完整18-case v5。所有原 hard gates 不变。

```ini
THRESHOLD_RELAXATION = NO
WATERLINE_NODE_CHANGE = NO
CENTRAL_PATCH_TOPOLOGY_CHANGE = NO
SINGLE_CENTER_FAN = NO
PHYSICS_CHANGE = NO
```
