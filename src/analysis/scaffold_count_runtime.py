
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

in_dir = Path("in/")
out_dir = Path("out/")

fig_path = out_dir / 'figs'
fig_path.mkdir(exist_ok= True)

# ––––– Model time ~ num ref chromosomes

# 14min w/ [0]; 43 min with 2; 121min with 3
X = np.array([[1,2, 2,3],[14,39,43,130]]).T
x = X[:,0]
y = np.log(X[:,1])
n = len(X)
p = X.shape[1]


MM = np.column_stack([np.ones(n), x])
COV = MM.T @ MM
pars = np.linalg.solve(COV , MM.T @ y)
B0,B1 = pars

resid = y - MM @ pars
rmse = np.sqrt(np.mean(resid**2))

# unbiased residual variance: RSS / df
df = (n-p)
sigma2 = np.sum(resid**2) / df


cov_pars = sigma2 * np.linalg.inv(COV) # (2,2) covariance of [B0,B1]
SE = np.sqrt(np.diag(cov_pars))              # [SE_B0, SE_B1]


t_stats = pars / SE                          # t = estimate / SE
p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stats), df=df))

print(f"RMSE = {rmse}")
print(f"** LM(time ~ num chromosomes)**",
      f"\tp-value = {p_vals}",
      f"\tB = {pars}",
      sep = "\n")

# ––––– predict time with 28 chromosomes
n_chr = 28
E_t = np.exp(np.array([1,n_chr]) @ pars)
print(f"Expected runtime for {n_chr} chromosomes = {np.floor(E_t/60):,.0f} hr {(E_t/60)%1 * 60:.0f} min")

N_new = 500
_max = 32
new_dat = np.column_stack([np.ones(N_new), np.linspace(0,_max, num=N_new)])
pred_t_at_n_chr = np.exp(new_dat @ pars)[np.argmin(abs(new_dat[:,1] - n_chr))]

fig, ax = plt.subplots(figsize = (6,6))
ax.plot(new_dat[:,1], np.exp(new_dat @ pars), linewidth = .7, c = "steelblue")
ax.scatter(x,np.exp(y))
ax.fill_between(new_dat[:,1], np.exp(new_dat @ (pars + SE)), np.exp(new_dat @ (pars - SE)), alpha = .4, color = "steelblue")
ax.axvline(x = n_chr, c = "#5a0d87", linestyle = "--")
ax.axhline(y = pred_t_at_n_chr, c = "#5a0d87", linestyle = "--")
ax.scatter(n_chr, pred_t_at_n_chr , c= "#5a0d87", s = 100, marker = 'x')
ax.text(0, pred_t_at_n_chr + pred_t_at_n_chr/2 , f"{int(pred_t_at_n_chr/60):.2e} hours", c = "#5a0d87")
ax.set_yscale("log")
ax.set_xticks(range(0,_max,2))
ax.set_ylabel("Time (mins)")
ax.set_xlabel("Number of scaffolds")
ax.set_title("`filter_reads()` runtime with increasing reference scaffold count")
plt.tight_layout()
plt.savefig(fig_path / "runtime_versus_index_scaffold_number.png", dpi = 150)

