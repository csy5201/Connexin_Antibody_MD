import matplotlib.pyplot as plt
import numpy as np

def read_xvg(filename):
    """读取GROMACS生成的xvg文件"""
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):
                values = line.split()
                if len(values) >= 2:
                    x.append(float(values[0]))
                    y.append(float(values[1]))
    return np.array(x), np.array(y)

# 1. 读取数据
residue, rmsf = read_xvg('rmsf.xvg')
rmsf_angstrom = rmsf * 10  # 将单位从 nm 转换为 Å

# 2. 检测链断点 (基于残基号回落到1的逻辑)
breaks = np.where(np.diff(residue) < 0)[0] + 1
chain_indices = np.split(np.arange(len(residue)), breaks)

# 3. 数据分组
connexin_chains = chain_indices[:6]  # Chains A-F: 对应残基 1-257
light_chain = chain_indices[6]        # Chain G: 对应残基 1-113
heavy_chain = chain_indices[7]        # Chain H: 对应残基 1-119

# ===== 创建绘图画布 (1x3 布局) =====
fig = plt.figure(figsize=(18, 5.5))

# ========== 图1: Connexin (Chains A-F) - 6条曲线叠加 ==========
ax1 = plt.subplot(131)

chain_labels = ['Chain A', 'Chain B', 'Chain C', 'Chain D', 'Chain E', 'Chain F']
# 使用渐变蓝色系
colors_connexin = ['#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7']

for i, idx_list in enumerate(connexin_chains):
    ax1.plot(residue[idx_list], rmsf_angstrom[idx_list], 
             color=colors_connexin[i], linewidth=1.8, 
             label=chain_labels[i], alpha=0.85)

# 添加平均黑虚线
all_connexin_rmsf = np.array([rmsf_angstrom[idx] for idx in connexin_chains])
avg_connexin_rmsf = np.mean(all_connexin_rmsf, axis=0)
ax1.plot(residue[connexin_chains[0]], avg_connexin_rmsf, 
         'k--', linewidth=2.5, label='Average', alpha=0.8, zorder=10)

ax1.set_xlabel('Residue Number', fontsize=13, fontweight='bold')
ax1.set_ylabel('RMSF (Å)', fontsize=13, fontweight='bold')
ax1.set_title('Connexin43 (Chains A-F)\nResidues 1-257', 
              fontsize=14, fontweight='bold')
ax1.legend(loc='upper right', fontsize=9, framealpha=0.95, ncol=2)
ax1.grid(True, alpha=0.3, axis='y')
ax1.set_xlim(0, 260)

# 统计信息框 (左上角)
cx_avg = np.mean(avg_connexin_rmsf)
cx_max = np.max(avg_connexin_rmsf)
cx_max_pos = residue[connexin_chains[0]][np.argmax(avg_connexin_rmsf)]
textstr = f'Average RMSF:\nMean: {cx_avg:.2f} Å\nMax: {cx_max:.2f} Å\n(at residue {int(cx_max_pos)})'
ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
         verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9))

# ========== 图2: Light Chain (Chain G) ==========
ax2 = plt.subplot(132)

ax2.plot(residue[light_chain], rmsf_angstrom[light_chain], 
         color='#238b45', linewidth=2.5, label='Light Chain (G)', alpha=0.9)

ax2.set_xlabel('Residue Number', fontsize=13, fontweight='bold')
ax2.set_ylabel('RMSF (Å)', fontsize=13, fontweight='bold')
ax2.set_title('Antibody Light Chain (G)\nResidues 1-113', 
              fontsize=14, fontweight='bold')
ax2.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_xlim(0, 115)

# 统计信息框 (修正位置：从右上移到左上，对齐方式改为left)
lc_avg = np.mean(rmsf_angstrom[light_chain])
lc_max = np.max(rmsf_angstrom[light_chain])
lc_max_pos = residue[light_chain][np.argmax(rmsf_angstrom[light_chain])]
textstr = f'Mean: {lc_avg:.2f} Å\nMax: {lc_max:.2f} Å\n(at residue {int(lc_max_pos)})'
ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, fontsize=10,
         verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# ========== 图3: Heavy Chain (Chain H) ==========
ax3 = plt.subplot(133)

ax3.plot(residue[heavy_chain], rmsf_angstrom[heavy_chain], 
         color='#cb181d', linewidth=2.5, label='Heavy Chain (H)', alpha=0.9)

ax3.set_xlabel('Residue Number', fontsize=13, fontweight='bold')
ax3.set_ylabel('RMSF (Å)', fontsize=13, fontweight='bold')
ax3.set_title('Antibody Heavy Chain (H)\nResidues 1-119', 
              fontsize=14, fontweight='bold')
ax3.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax3.grid(True, alpha=0.3, axis='y')
ax3.set_xlim(0, 122)

# 统计信息框 (修正位置：从右上移到左上，对齐方式改为left)
hc_avg = np.mean(rmsf_angstrom[heavy_chain])
hc_max = np.max(rmsf_angstrom[heavy_chain])
hc_max_pos = residue[heavy_chain][np.argmax(rmsf_angstrom[heavy_chain])]
textstr = f'Mean: {hc_avg:.2f} Å\nMax: {hc_max:.2f} Å\n(at residue {int(hc_max_pos)})'
ax3.text(0.02, 0.98, textstr, transform=ax3.transAxes, fontsize=10,
         verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.9))

# 调整整体布局，防止标签重叠
plt.tight_layout()

# 保存高质量图片
plt.savefig('rmsf_final_superimposed_clean.png', dpi=300, bbox_inches='tight')

# ===== 打印终端统计报告 =====
print("="*80)
print("RMSF ANALYSIS - FINAL CORRECTED VERSION")
print("="*80)

print(f"\n📊 Panel 1: Connexin43 (Chains A-F)")
print(f"   Avg RMSF: {cx_avg:.2f} Å | Max: {cx_max:.2f} Å at res {int(cx_max_pos)}")

print(f"\n📊 Panel 2: Light Chain (G)")
print(f"   Avg RMSF: {lc_avg:.2f} Å | Max: {lc_max:.2f} Å at res {int(lc_max_pos)}")

print(f"\n📊 Panel 3: Heavy Chain (H)")
print(f"   Avg RMSF: {hc_avg:.2f} Å | Max: {hc_max:.2f} Å at res {int(hc_max_pos)}")

print("\n" + "="*80)
print("✅ 图片已完美保存: rmsf_final_superimposed_clean.png")
print("="*80)

plt.show()
