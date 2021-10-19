import pandas,seaborn,math,sys

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

plot_qval=pandas.DataFrame(0.0,index=list(set(df_top['alias'])),columns=sorted(list(set(df['cluster']))))
plot_de=pandas.DataFrame(float('nan'),index=list(set(df_top['alias'])),columns=sorted(list(set(df['cluster']))))

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

plot_qval=plot_qval.loc[:, (plot_qval != 0).any(axis=0)] ###remove columns with only zeros
plot_qval[plot_qval > 20] = 20                

ax=seaborn.clustermap(plot_qval,cmap='Blues',xticklabels=True,yticklabels=True,row_colors=col_row,col_cluster=False,row_cluster=True,mask=(plot_qval==0),linewidths = 0.01, linecolor='black')
ax1 = ax.ax_heatmap
ax1.set_facecolor("#b3b3b3ff")
ax.ax_cbar.set_title('-log10(qval)')
heatmap_pos = ax.ax_heatmap.get_position()
ax.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width*0.25, heatmap_pos.height])
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment='right')
ax.ax_row_colors.tick_params(bottom=False)
ax.ax_row_colors.set_xticklabels('')


###get the order from the dendogram
indiNew=ax.dendrogram_row.reordered_ind
indiNewNames=[]
for i in range(len(indiNew)):
    indiNewNames.append(plot_qval.index.tolist()[indiNew[i]])

plot_de=plot_de.reindex(indiNewNames)


plot_de = plot_de[plot_de.columns.intersection(plot_qval.columns)]

bx=seaborn.clustermap(plot_de,cmap='Blues',xticklabels=True,yticklabels=True,row_colors=col_row,col_cluster=False,row_cluster=False,vmax=1,vmin=0,linewidths = 0.01, linecolor='black')
bx.cax.set_visible(False) #remove cbar
bx1 = bx.ax_heatmap
heatmap_pos = bx.ax_heatmap.get_position()
bx.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width*0.25, heatmap_pos.height])
bx1.set_xticklabels(bx1.get_xticklabels(), rotation=90, horizontalalignment='right')
bx1.set_facecolor("#b3b3b3ff")
bx.ax_row_colors.tick_params(bottom=False) 
bx.ax_row_colors.set_xticklabels('')


ax.savefig(OUT_SPEC)
bx.savefig(OUT_DE)
