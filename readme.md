以 $\Theta$ 为中心，向左右建立nodes， 间隔BatchSize

每经过一个node， 计算一次 BL , BR


创建变长数组 subnodes 记录 $\tau, \Theta$ 间的 nodes

BLM,BRM 的长度为 length(subnodes)


