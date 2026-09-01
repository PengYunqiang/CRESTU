# CRESTU-1F Phase 3.2B v5 基线冻结

```ini
V5_MESH_REPRODUCIBILITY_BASELINE = PASS
FROZEN_MESH_CASES = 18
FROZEN_COMPONENT_QUALITY = 72/72 PASS
FROZEN_FREQUENCY_ISOLATION = 6/6 PASS
QUALITY_CONTROLLED_V2_DEFAULT = OFF
LEGACY_Q_ONLY_GEOMETRY = 5/5 EXACT_HASH_PASS
BASELINE_COMMIT = 4a1bb81f6db6e577a1b45263d9e224a9260a5153
BASELINE_TAG_OBJECT = e761035e156e54814b751b6b80f40e25803c6667
PRODUCTION_SOURCE_SHA256 = 7dee2e9eed9c3b80755559f913b90a93e3ba59eca9bafb621c86fe58c82f3783
```

本阶段冻结 Phase 3.2A 已接受的 18 行 v5 manifest、72 行组件质量、6 行频率隔离、三份 topology design、realized spacing、outer extent、组件/组合网格哈希和 source/config/control 哈希。没有重新设计 topology，也没有把 V2 改为默认模式。

## 冻结网格族

| Level | eta | body division | FS radial | bottom radial | wall nz | theta | bottom core |
|---|---:|---:|---:|---:|---:|---:|---:|
| MESH_L1 | 0.0707009403590185 | 4 | 4+2 | 4+2 | 3 | 64 | 4 |
| MESH_L2 | 0.0550177552570247 | 5 | 5+3 | 5+3 | 4 | 80 | 7 |
| MESH_L3 | 0.0412542094946932 | 7 | 6+4 | 6+4 | 5 | 112 | 12 |

每个 Phase 3.2B case 都必须在矩阵组装前与冻结行逐字段核验；不一致即停止，变化后的网格不得称为 v5。
