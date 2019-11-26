"""
Annotated heatmaps for spectral clustering result
==================

"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import matplotlib as mpl


# Load the example flights dataset and conver to long-form
flights_long = sns.load_dataset("clustering_sc_vec25")
flights = flights_long.pivot("Valley", "Cluster","Overlap")

# Draw a heatmap with the numeric values in each cell
font = {'family' : 'Arial', 'size' : 8}
sns.set(font_scale = 0.4)
#cmapNew = cm.get_cmap('Default')
#cmapNew_r = reverse_colourmap(cmapNew)


#f, ax = plt.subplots(figsize=(4.8, 2.36),dpi = 300)
grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
f,(ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws, figsize=(2.8, 3),dpi = 600)
#cmap_hm = sns.cubehelix_palette(as_cmap=True);
cmap_hm = sns.cubehelix_palette(light=1, as_cmap=True);
ax = sns.heatmap(flights, annot=False, fmt="1.1f", linecolor='black', linewidths=.1, ax=ax, cmap=cmap_hm, vmin=0, vmax=100, square=True, cbar=1, cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal"})
ax.invert_yaxis()
#f.subplots_adjust(top = 1.1)
#f.subplots_adjust(bottom=0.15)
plt.savefig('/Users/Ahmed/Downloads/figure_sc.pdf', bbox_inches='tight')
plt.show()  
plt.close()