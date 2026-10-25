# 排名对照输出接口

来源：已冻结的排名对照结果。

每例 `cases/{code}.npz`：`leads[4]`，`ids[n]`，`error_xy[4,n,2]`，`target[4,2]`，`translated_xy[4,n,2]`；`early_order[n]`，`late_orders[4,n]`，`random_orders[128,n]`。

对 early/late/random family：`{family}_points[4,repeats,2,2]`、`{family}_weights[4,repeats,2,n]`、`{family}_distance[4,repeats,2]`、`{family}_M_direction[4,repeats]`。两方向索引 D=0、I=1；repeats=1/1/128。权重均按原成员索引。M 逐例逐时效逐排列选 D/I 最小距离，精确相等取 D。`H_points[4,2]`、`H_weights[4,n]`、`H_distance[4]`。

case_results.csv.gz：3224 行，中心/台风/起报/code/lead/n；D/I/M × early/late/random 的距离（random为128次各自最优后均值）、H；D/M 的 early_minus_late、late_minus_H、early_minus_random、early_minus_H、random_minus_H。所有对照另存大于+1km、小于−1km指示，避免把数值近零当方向结论。I_better为反向胜正向的比例。

summary.csv：center/lead/metric/mean/ci_low/ci_high/cases/storms。按起报台风内均值、台风等权，5000次原配对整台风抽样；lead包含48、72、96、120和四时效均值。未做多重性修正。

random_mc.csv：D/M 的总体 Monte Carlo SE，以及前64与后64排列的加权均值和差。仅诊断，不据结果增加次数。

