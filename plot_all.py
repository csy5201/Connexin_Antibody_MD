import matplotlib.pyplot as plt
import numpy as np

def read_xvg(filename):
    """读取xvg文件"""
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            # 跳过注释行
            if not line.startswith('#') and not line.startswith('@'):
                values = line.split()
                if len(values) >= 2:
                    x.append(float(values[0]))
                    y.append(float(values[1]))
    return np.array(x), np.array(y)

# 创建图形
fig = plt.figure(figsize=(15, 5))

# ===== 1. Energy Minimization =====
ax1 = plt.subplot(131)
steps, energy = read_xvg('energy_minimization.xvg')
ax1.plot(steps, energy, 'b-', linewidth=2)
ax1.set_xlabel('Step', fontsize=12, fontweight='bold')
ax1.set_ylabel('Potential Energy (kJ/mol)', fontsize=12, fontweight='bold')
ax1.set_title('Step 6: Energy Minimization', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
# 添加文本框显示关键信息
textstr = f'Initial: {energy[0]:.0f}\nFinal: {energy[-1]:.0f}\nΔE: {energy[0]-energy[-1]:.0f} kJ/mol'
ax1.text(0.65, 0.95, textstr, transform=ax1.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# ===== 2. RMSD =====
ax2 = plt.subplot(132)
time, rmsd = read_xvg('rmsd.xvg')
time_ns = time / 1000  # ps转ns
rmsd_angstrom = rmsd * 10  # nm转Å
ax2.plot(time_ns, rmsd_angstrom, 'g-', linewidth=2)
ax2.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
ax2.set_ylabel('RMSD (Å)', fontsize=12, fontweight='bold')
ax2.set_title('Step 8: Backbone RMSD', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
textstr = f'Average: {np.mean(rmsd_angstrom):.2f} Å\nMax: {np.max(rmsd_angstrom):.2f} Å'
ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

# ===== 3. RMSF =====
ax3 = plt.subplot(133)
residue, rmsf = read_xvg('rmsf.xvg')
rmsf_angstrom = rmsf * 10  # nm转Å
ax3.plot(residue, rmsf_angstrom, 'r-', linewidth=1.5)
ax3.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
ax3.set_ylabel('RMSF (Å)', fontsize=12, fontweight='bold')
ax3.set_title('Step 8: Residue Fluctuation', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')
textstr = f'Average: {np.mean(rmsf_angstrom):.2f} Å\nMax: {np.max(rmsf_angstrom):.2f} Å'
ax3.text(0.65, 0.95, textstr, transform=ax3.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

plt.tight_layout()
plt.savefig('all_results.png', dpi=300, bbox_inches='tight')
print("✅ 图片已保存: all_results.png")

# 打印统计信息
print("\n" + "="*70)
print("MEMBRANE EQUILIBRATION ANALYSIS RESULTS")
print("="*70)
print(f"\n📊 Step 6 - Energy Minimization:")
print(f"   Initial Energy: {energy[0]:>15,.2f} kJ/mol")
print(f"   Final Energy:   {energy[-1]:>15,.2f} kJ/mol")
print(f"   Energy Drop:    {energy[0]-energy[-1]:>15,.2f} kJ/mol")
print(f"   Total Steps:    {int(steps[-1]):>15,}")

print(f"\n📊 Step 8 - RMSD Analysis:")
print(f"   Simulation Time: {time_ns[-1]:.2f} ns")
print(f"   Average RMSD:    {np.mean(rmsd_angstrom):>10.2f} Å")
print(f"   Min RMSD:        {np.min(rmsd_angstrom):>10.2f} Å")
print(f"   Max RMSD:        {np.max(rmsd_angstrom):>10.2f} Å")
print(f"   Final RMSD:      {rmsd_angstrom[-1]:>10.2f} Å")
print(f"   Std Dev:         {np.std(rmsd_angstrom):>10.2f} Å")

print(f"\n📊 Step 8 - RMSF Analysis:")
print(f"   Total Residues:  {len(residue):>10}")
print(f"   Average RMSF:    {np.mean(rmsf_angstrom):>10.2f} Å")
print(f"   Max RMSF:        {np.max(rmsf_angstrom):>10.2f} Å (Residue {int(residue[np.argmax(rmsf_angstrom)])})")
print(f"   Std Dev:         {np.std(rmsf_angstrom):>10.2f} Å")

# 找出最灵活的残基
threshold = np.mean(rmsf_angstrom) + np.std(rmsf_angstrom)
flexible = [(int(residue[i]), rmsf_angstrom[i]) for i in range(len(residue)) if rmsf_angstrom[i] > threshold]
if flexible:
    print(f"\n   Most Flexible Residues (RMSF > {threshold:.2f} Å):")
    for res_id, res_rmsf in sorted(flexible, key=lambda x: x[1], reverse=True)[:5]:
        print(f"      Residue {res_id:>4}: {res_rmsf:>6.2f} Å")

print("\n" + "="*70)
print("✅ Analysis Complete! Check 'all_results.png' for visualization.")
print("="*70)

plt.show()