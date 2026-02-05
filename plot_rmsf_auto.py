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

# 检测链断点
breaks = np.where(np.diff(residue) < 0)[0] + 1
chain_indices = np.split(np.arange(len(residue)), breaks)

print("="*70)
print("所有链的残基范围信息")
print("="*70)

chain_names = ['Chain A (Connexin)', 'Chain B (Connexin)', 'Chain C (Connexin)', 
               'Chain D (Connexin)', 'Chain E (Connexin)', 'Chain F (Connexin)',
               'Chain G (Light Chain)', 'Chain H (Heavy Chain)']

for i, (idx_list, name) in enumerate(zip(chain_indices, chain_names)):
    res_start = int(residue[idx_list[0]])
    res_end = int(residue[idx_list[-1]])
    length = len(idx_list)
    
    print(f"\n{name}:")
    print(f"  残基范围: {res_start} - {res_end}")
    print(f"  残基数量: {length}")
    print(f"  数据索引: {idx_list[0]} - {idx_list[-1]}")

print("\n" + "="*70)
