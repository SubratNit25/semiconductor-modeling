import pandas as pd
import matplotlib.pyplot as plt

try:
    df = pd.read_csv('rd_results.csv')
except FileNotFoundError:
    print("Error: 'rd_results.csv' not found. Run the C program first.")
    exit()

stress_df = df[df['Phase'] == 'Stress']
recovery_df = df[df['Phase'] == 'Recovery']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot Phase 1: Stress
ax1.plot(stress_df['Time'], stress_df['NIT1'], label=r'$N_{IT,1}$ (Channel/IL)', linewidth=2.5, color='#1f77b4')
ax1.plot(stress_df['Time'], stress_df['NIT2'], label=r'$N_{IT,2}$ (IL/HK)', linewidth=2.5, linestyle='--', color='#ff7f0e')
ax1.plot(stress_df['Time'], stress_df['Total_NIT'], label='Total $N_{IT}$', color='black', linewidth=3)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('Phase 1: Stress Kinetics ($V_G$ = -1.5V)', fontsize=14)
ax1.set_xlabel('Stress Time (s)', fontsize=12)
ax1.set_ylabel(r'Trap Density ($cm^{-2}$)', fontsize=12)
ax1.grid(True, which="both", ls="--", alpha=0.5)
ax1.legend(fontsize=11)

# Plot Phase 2: Recovery
ax2.plot(recovery_df['Time'], recovery_df['NIT1'], label=r'$N_{IT,1}$', linewidth=2.5, color='#1f77b4')
ax2.plot(recovery_df['Time'], recovery_df['NIT2'], label=r'$N_{IT,2}$', linewidth=2.5, linestyle='--', color='#ff7f0e')
ax2.plot(recovery_df['Time'], recovery_df['Total_NIT'], label='Total $N_{IT}$', color='black', linewidth=3)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title('Phase 2: Recovery Kinetics ($V_G$ = 0V)', fontsize=14)
ax2.set_xlabel('Recovery Time (s)', fontsize=12)
ax2.set_ylabel(r'Trap Density ($cm^{-2}$)', fontsize=12)
ax2.grid(True, which="both", ls="--", alpha=0.5)
ax2.legend(fontsize=11)

plt.tight_layout()
plt.show()
