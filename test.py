import matplotlib.pyplot as plt
import numpy as np

def read_xvg(filename):
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):
                values = line.split()
                if len(values) >= 2:
                    x.append(float(values[0]))
                    y.append(float(values[1]))
    return np.array(x), np.array(y)

residue, rmsf = read_xvg('rmsf.xvg')
rmsf_angstrom = rmsf * 10

# 检查数据
print("总残基数:", len(residue))
print("前20个残基编号:", residue[:20])
print("残基编号的差值（前50个）:", np.diff(residue[:50]))

# 找出所有不连续的地方
gaps = np.where(np.diff(residue) != 1)[0]
print(f"\n找到 {len(gaps)} 个不连续点")
print("不连续位置的前10个:", gaps[:10])
print("对应的残基号:", residue[gaps[:10]], "→", residue[gaps[:10] + 1])


