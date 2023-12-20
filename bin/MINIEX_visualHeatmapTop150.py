import pandas,seaborn,math,sys,natsort
from matplotlib.patches import Patch
from matplotlib.pyplot import gcf

DF_IN=sys.argv[1]

OUT_SPEC=sys.argv[2]
OUT_DE=sys.argv[3]
TOPS=int(sys.argv[4])
df=pandas.read_excel(DF_IN)

df_top=df.sort_values('borda_rank',ascending=True).head(TOPS)

filterToUse=list(set(df_top['alias']))

mino=list(set(df['qval_cluster']))
if 0.0 in mino:
    mino.remove(0.0)
mino=min(mino)

col_row=[]

# If the tissue name and cluster name are the same, simplify the redundant name
def simplify_cluster_name(cluster_name):
    tissue, cluster = cluster_name.split("_Cluster_")
    if tissue == cluster:
        return f"Cluster_{cluster}"
    else:
        return cluster_name.replace("_Cluster_", "-")
df["cluster"]=df["cluster"].astype(str).apply(lambda x: simplify_cluster_name(x))

plot_qval=pandas.DataFrame(0.0,index=list(set(df_top['alias'])),columns=natsort.natsorted(list(set(df['cluster'])), alg=natsort.IGNORECASE))
plot_de=pandas.DataFrame(float('nan'),index=list(set(df_top['alias'])),columns=natsort.natsorted(list(set(df['cluster'])), alg=natsort.IGNORECASE))

for index,row in df.iterrows():
    if row['alias'] in filterToUse:
        if row['qval_cluster'] != 0.0:
            plot_qval[row['cluster']][row['alias']]=round(-math.log10(row['qval_cluster']))
        else:
            plot_qval[row['cluster']][row['alias']]=round(-math.log10(mino))
        
        plot_de[row['cluster']][row['alias']]=row['isTF_DE']

        if row['hasTFrelevantGOterm'] == 'relevant_known_TF':
            col_row.append([row['alias'],'green'])
        elif row['hasTFrelevantGOterm'] == 'known_TF':
            col_row.append([row['alias'],'orange'])     
        else:
            col_row.append([row['alias'],'gray']) 

plot_qval=plot_qval.fillna(0.0).sort_index()
col_row=pandas.DataFrame(col_row)
col_row=col_row.drop_duplicates().set_index(0).sort_index()

tf2alias=pandas.Series(df.alias.values,index=df.TF).to_dict()

def get_x_font_size_for_heatmap(dataframe):
    font_size = None
    if len(dataframe.columns) > 20:
        font_size = 5
    else:
        font_size = 7
    return font_size

def get_y_font_size_for_heatmap(dataframe):
    font_size = None
    if len(dataframe) > 300:
        font_size = 1
    elif len(dataframe) > 200:
        font_size = 2
    elif len(dataframe) > 100:
        font_size = 4
    elif len(dataframe) > 20:
        font_size = 8
    else:
        font_size = 12
    return font_size


# cluster specificity heatmap
plot_qval=plot_qval.loc[:, (plot_qval != 0).any(axis=0)] ###remove columns with only zeros
plot_qval[plot_qval > 20] = 20

ax=seaborn.clustermap(plot_qval,cmap='Blues',xticklabels=True,yticklabels=True,row_colors=col_row,col_cluster=False,row_cluster=True,mask=(plot_qval==0),linewidths=0.01,linecolor='dimgrey')
ax1 = ax.ax_heatmap
ax1.set_facecolor("#b3b3b3ff")
ax.ax_cbar.set_title('-log10(qval)')
heatmap_pos = ax.ax_heatmap.get_position()
ax.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width*0.5, heatmap_pos.height])
fontsize = min([get_x_font_size_for_heatmap(plot_qval),get_y_font_size_for_heatmap(plot_qval)]) # enforce the same x and y font size
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=fontsize)
ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=fontsize)
ax.ax_row_colors.tick_params(bottom=False)
ax.ax_row_colors.set_xticklabels('')

# Add legend to plots
lut = {"Associated with GO term of interest": "#008000ff", "Associated with any other GO term": "#ffa500ff", "No GO info known": "#808080ff"}
handles = [Patch(facecolor=lut[name], linewidth=0.5, edgecolor="#303030ff") for name in lut]
ax.ax_col_dendrogram.legend(handles, lut, title='Legend',
           bbox_to_anchor=(1, 1), bbox_transform=gcf().transFigure, loc='upper right')


###get the order from the dendogram
indiNew=ax.dendrogram_row.reordered_ind
indiNewNames=[]
for i in range(len(indiNew)):
    indiNewNames.append(plot_qval.index.tolist()[indiNew[i]])

# DE calls
plot_de=plot_de.reindex(indiNewNames)
plot_de = plot_de[plot_de.columns.intersection(plot_qval.columns)]

bx=seaborn.clustermap(plot_de,cmap='Blues',xticklabels=True,yticklabels=True,row_colors=col_row,col_cluster=False,row_cluster=False,vmax=1,vmin=0,linewidths=0.01,linecolor='dimgrey')
bx.cax.set_visible(False) #remove cbar
bx1 = bx.ax_heatmap
heatmap_pos = bx.ax_heatmap.get_position()
bx.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width*0.5, heatmap_pos.height])
fontsize = min([get_x_font_size_for_heatmap(plot_de),get_y_font_size_for_heatmap(plot_de)]) # enforce the same x and y font size
bx1.set_xticklabels(bx1.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=fontsize)
bx1.set_yticklabels(ax1.get_yticklabels(), fontsize=fontsize)
bx1.set_facecolor("#b3b3b3ff")
bx.ax_row_colors.tick_params(bottom=False) 
bx.ax_row_colors.set_xticklabels('')

# Add legend to plot
lut = {"Associated with GO term of interest": "#008000ff", "Associated with any other GO term": "#ffa500ff", "No GO info known": "#808080ff",  \
        "Differentially expressed (DE) TF": "#08306bff", "Expressed (in 10% of cells) TF, not DE": "#f7fbffff", "Not DE, not expressed": "#b3b3b3ff"}
handles = [Patch(facecolor=lut[name], linewidth=0.5, edgecolor="#303030ff") for name in lut]
bx.ax_col_dendrogram.legend(handles, lut, title='Legend',
           bbox_to_anchor=(1, 1), bbox_transform=gcf().transFigure, loc='upper right')


ax.savefig(f"{OUT_SPEC}.svg")
ax.savefig(f"{OUT_SPEC}.png", dpi=600)
bx.savefig(f"{OUT_DE}.svg")
bx.savefig(f"{OUT_DE}.png", dpi=600)



