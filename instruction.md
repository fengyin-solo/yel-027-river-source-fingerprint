# 入库泥沙来源分区的指纹重配

明早要据此交一份悬浮泥沙来源浓度复核单；
入口是 `/app/sediment_source_fingerprint.py`；
水库管理站刚完成一次汛后取样；
实验室合并批次时丢失了每份泥沙样本所属的集水区编号，只留下样本干重、各等分样浓度和三个来源区的汇总值；


## 批次材料

六个顶层名称组成唯一合法根形状；UTF-8 解码过程中若出现重复名称即拒收：；
```text
watershed_id
campaign_id
samples
zones
calibration
design_confidence
```
 `samples` 有 12 到 24 份样本；每份样本只含四个键：；`watershed_id` 和 `campaign_id` 的字符数应落在 2–48，且不允许首尾空白或控制字符；前者以 `RIVER-` 起头；
```text
sample_ref
aliquot_count
dry_mass_g
concentrations_ug_l
```
`dry_mass_g` 的有效集合是有限且非零的正实数，上限为 100；`concentrations_ug_l` 的条目数应与 `aliquot_count` 完全相等，元素取 0–1000 整数；真/假值不属于这里可接受的整数或小数； `zones` 有 3 到 5 个来源区，每个来源区的成员固定为：；`sample_ref` 的格式为“大写字母 + 2 至 15 个字母数字/下划线/连字符”，同批不可复现；`aliquot_count` 是 10 到 5000 的整数；
```text
zone_id
catchment_area_ha
sample_count
dry_mass_total
analyte_total
analyte_sq_total
```
`sample_count` 为正整数；
剩余三个积分量应可表示为非负有限实数；
 来源区汇总由标签丢失前的原始记录按下列式子产生：；
`zone_id` 匹配 `[A-Z][A-Z0-9_-]{0,23}` 且不得重复；
`catchment_area_ha` 是 2 到 1000000000 的整数；

```text
dry_mass_total = sum(dry_mass_g)
analyte_total = sum(dry_mass_g * sum(concentrations_ug_l))
analyte_sq_total = sum(dry_mass_g * sum(concentration * concentration))
```
所有样本必须完整划分到各来源区，任何一份都不得重复或遗漏； 归属并列的规范化顺序由区编号序列和区内 `sample_ref` 序列共同组成，取字典序最小者；闭合失败时不生成结果，退出 1，stderr 固定 `运行拒收: BAD_MODEL`； `calibration` 仅有两个键：；每区的 `sample_count` 是实际样本数；三项汇总与样本数据允许 `1e-7` 的绝对误差；
```text
z_score
tie_break
```
`design_confidence` 位于 0.8 到 0.99；标准倍数 `z_score` 且须满足 `0 < z_score <= 10`；`tie_break` 固定为字符串 `lexicographic`；

## 浓度复核

恢复每个区的样本后，区内浓度按干重计算：
```text
zone_mean_h = analyte_total_h / dry_mass_total_h
```
区内二阶量按每份样本的等分样展开：
```text
sum_w2_h = sum(dry_mass_g_hp * dry_mass_g_hp for each aliquot in sample p)
effective_n_h = dry_mass_total_h * dry_mass_total_h / sum_w2_h
s2_h = sum(dry_mass_g * (concentration - zone_mean_h)^2)
       / (dry_mass_total_h - sum_w2_h / dry_mass_total_h)
```
外层的离散量求和要遍历该区所有样本的全部等分样；全流域平均浓度按集水面积份额合成：；
```text
basin_mean_concentration_ug_l = sum(share_h * zone_mean_h)
share_h = catchment_area_ha_h / total_catchment_area_ha
```
标准误是各来源区的有限总体修正方差之和再开方：
```text
design_variance_piece_h = share_h^2
                  * max(0, 1 - effective_n_h / catchment_area_ha_h)
                  * s2_h / effective_n_h
standard_error_ug_l = sqrt(sum(variance_term_h))
```
置信下界与上界分别为均值减、加 `z_score * standard_error_ug_l`；`design_confidence` 参与字段校验，但不替代已锁定的 `z_score`；

## 复核单与交付台账

程序接受材料、复核单两个位置参数；交班复跑时还可追加台账：；
```bash
python3 /app/sediment_source_fingerprint.py <批次材料> <复核单>
python3 /app/sediment_source_fingerprint.py <批次材料> <复核单> <交付台账>
```
结果文件只允许一个 JSON 对象和结尾换行，成员顺序不可调换：
```text
watershed_id
status
zone_order
zone_samples
zone_mean_concentration_ug_l
basin_mean_concentration_ug_l
standard_error_ug_l
ci_low_ug_l
ci_high_ug_l
```
`zone_samples` 和 `zone_mean_concentration_ug_l` 的元素位置与 `zone_order` 同步，区内 `sample_ref` 仍按升序展开；数值直接按双精度写出；正常返回码为 0 时，标准输出与标准错误均为空； 台账采用只追加的 UTF-8 JSONL；每行四键，紧凑 JSON，行尾单个 LF：；`status` 写 `complete`；`zone_order` 按 `zone_id` 升序；
```text
watershed_id
input_sha256
result_sha256
basin_mean_concentration_ug_l
```
复核单摘要变化时，流水在其末尾登记新行；若流水出现输入摘要重放冲突或结果摘要归属冲突，打印 `运行拒收: BAD_LEDGER` 并退出 1； 复核单和台账应同时发布或同时保持旧字节；成功发布后不得残留临时文件，输入材料本身不可改写；哈希直接取文件磁盘字节，不做换行或编码归一化；相同输入成功重跑时，台账一个字节都不变；

## 拒收

流水重放或格式异常统一以 `BAD_LEDGER` 拒绝；读写位置异常时终端仅显示 `运行拒收: IO_FAILURE`； 异常路径禁止向终端展示 Python `Traceback`，也不得覆盖旧复核单或旧台账；读取或校验材料失败，退出码 1，标准错误仅含 `运行拒收: BAD_INPUT`； 归属不能闭合时用前述 `BAD_MODEL`；
