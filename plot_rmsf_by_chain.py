
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

breaks = np.where(np.diff(residue) < 0)[0] + 1
print(f"检测到断点位置: {breaks}")
print(f"总共 {len(breaks) + 1} 条链")

chain_indices = np.split(np.arange(len(residue)), breaks)

for i, chain in enumerate(chain_indices):
    print(f"链 {i+1}: 索引 {chain[0]:4d}-{chain[-1]:4d}, "
          f"残基号 {int(residue[chain[0]]):3d}-{int(residue[chain[-1]]):3d}, "
          f"长度 {len(chain)}")

connexin_chains = chain_indices[:6]
antibody_chains = chain_indices[6:8]

print(f"\nConnexin: {len(connexin_chains)} 条链")
print(f"Antibody: {len(antibody_chains)} 条链")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

chain_labels = ['A', 'B', 'C', 'D', 'E', 'F']
colors_connexin = ['#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7']

offset = 0
for i, idx_list in enumerate(connexin_chains):
    x_plot = np.arange(len(idx_list)) + offset
    ax1.plot(x_plot, rmsf_angstrom[idx_list], 
             color=colors_connexin[i], linewidth=1.5, 
             label=f'Chain {chain_labels[i]}', alpha=0.9)
    offset += len(idx_list)

ax1.set_xlabel('Residue Index (Sequential across all chains)', fontsize=13, fontweight='bold')
ax1.set_ylabel('RMSF (Å)', fontsize=13, fontweight='bold')
ax1.set_title('Connexin43 (Chains A-F) - Residue Fluctuation', fontsize=15, fontweight='bold')
ax1.legend(loc='upper right', fontsize=10, ncol=2)
ax1.grid(True, alpha=0.3, axis='y')

offset = 0
for i in range(len(connexin_chains) - 1):
    offset += len(connexin_chains[i])
    ax1.axvline(x=offset, color='gray', linestyle='--', alpha=0.3, linewidth=0.8)

all_connexin = np.concatenate([rmsf_angstrom[idx] for idx in connexin_chains])
textstr = f'All Connexin:\nAvg: {np.mean(all_connexin):.2f} Å\nMax: {np.max(all_connexin):.2f} Å\nStd: {np.std(all_connexin):.2f} Å'
ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

antibody_labels = ['Light Chain (G)', 'Heavy Chain (H)']
antibody_colors = ['#d7301f', '#2ca02c']

offset = 0
for i, idx_list in enumerate(antibody_chains):
    x_plot = np.arange(len(idx_list)) + offset
    ax2.plot(x_plot, rmsf_angstrom[idx_list], 
             color=antibody_colors[i], linewidth=2, 
             label=antibody_labels[i], alpha=0.9)
    offset += len(idx_list)

if len(antibody_chains) > 1:
    ax2.axvline(x=len(antibody_chains[0]), color='gray', linestyle='--', alpha=0.5, linewidth=1)

ax2.set_xlabel('Residue Index (Sequential)', fontsize=13, fontweight='bold')
ax2.set_ylabel('RMSF (Å)', fontsize=13, fontweight='bold')
ax2.set_title('M1 Antibody (Chains G-H) - Residue Fluctuation', fontsize=15, fontweight='bold')
ax2.legend(loc='upper right', fontsize=11)
ax2.grid(True, alpha=0.3, axis='y')

all_antibody = np.concatenate([rmsf_angstrom[idx] for idx in antibody_chains])
textstr = f'All Antibody:\nAvg: {np.mean(all_antibody):.2f} Å\nMax: {np.max(all_antibody):.2f} Å\nStd: {np.std(all_antibody):.2f} Å'
ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

plt.tight_layout()
plt.savefig('rmsf_by_chain_final.png', dpi=300, bbox_inches='tight')

print("\n" + "="*70)
print("RMSF ANALYSIS BY CHAIN")
print("="*70)

print("\n📊 Connexin43 Chains (A-F):")
for i, (idx_list, label) in enumerate(zip(connexin_chains, chain_labels)):
    chain_rmsf = rmsf_angstrom[idx_list]
    print(f"   Chain {label}: Avg={np.mean(chain_rmsf):5.2f} Å, "
          f"Max={np.max(chain_rmsf):5.2f} Å, "
          f"Std={np.std(chain_rmsf):5.2f} Å")

print(f"\n   Overall Connexin: Avg={np.mean(all_connexin):5.2f} Å, "
      f"Max={np.max(all_connexin):5.2f} Å")

print("\n📊 M1 Antibody Chains (G-H):")
for i, (idx_list, label) in enumerate(zip(antibody_chains, antibody_labels)):
    chain_rmsf = rmsf_angstrom[idx_list]
    print(f"   {label}: Avg={np.mean(chain_rmsf):5.2f} Å, "
          f"Max={np.max(chain_rmsf):5.2f} Å, "
          f"Std={np.std(chain_rmsf):5.2f} Å")

print(f"\n   Overall Antibody: Avg={np.mean(all_antibody):5.2f} Å, "
      f"Max={np.max(all_antibody):5.2f} Å")

max_idx = np.argmax(rmsf_angstrom)
max_chain = None
for i, chain in enumerate(chain_indices):
    if max_idx in chain:
        max_chain = i
        break

chain_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G (LC)', 'H (HC)']
print(f"\n🔥 Most flexible residue: "
      f"Index {max_idx}, "
      f"Chain {chain_names[max_chain]}, "
      f"RMSF = {rmsf_angstrom[max_idx]:.2f} Å")

print("\n" + "="*70)
print("✅ RMSF图已保存: rmsf_by_chain_final.png")
print("="*70)

plt.show()
