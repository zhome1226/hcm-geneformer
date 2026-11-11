# CausalBench2023：原始计数与扰动索引

来源为Replogle Perturb-seq的K562/RPE1细胞系，官方CausalBench固定代码指向同两个Figshare文件；不是HCM队列。X在官方归一化步骤前读取，本次完整核查数值域并保留全部发布基因。

- K562：310385细胞、8563基因；可用310385，隔离0，non-targeting对照10691，48捕获。
- RPE1：247914细胞、8749基因；可用247914，隔离0，non-targeting对照11485，56捕获。

计数原件只读，以真实行列索引交付。扰动基因Ensembl ID、gene/transcript、sgID_AB和gem_group完整保留；双sgRNA不作双细胞，NANOG不作对照。Ensembl列唯一，重复symbol不合并；缺少目标基因表达列的细胞单独标记，但不伪造列或删除其有效表达。

不套HCM阈值，不做归一化/插补/整合或网络推断；原mitopercent字段单位未独立核实而不擅自换算。当前对象无固定split，原官方参数化分层划分及可选strong-perturbation过滤配方存档，本次不新建训练划分或按扰动强度筛选。所有逐细胞/逐基因文件只留服务器。

准备后须独立验收并核查退出码、进程和源码一致性。
