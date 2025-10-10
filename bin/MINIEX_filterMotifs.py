import sys,pandas

TF_INFO_FILE         = sys.argv[1]
ENRICHED_MOTIFS_FILE = sys.argv[2]
OUTPUT_FILE          = sys.argv[3]
TF_FAMILY_EXTENSION  = sys.argv[4]

# read in the tf info file as a dataframe
header = ["transcription_factor", "family", "has_motif", "motifs"]
tf_info_df = pandas.read_csv(TF_INFO_FILE, sep='\t', header=None, names=header).fillna("")
# extract motif info
motif_df = tf_info_df[tf_info_df["has_motif"] == "Y"].copy()
motif_df["motifs"] = motif_df["motifs"].str.split(",")
# expand each motif in the motif list to a separate row and filter out empty motifs
motif_df = motif_df.explode("motifs")
motif_df = motif_df[motif_df["motifs"].str.strip() != ""]
# convert dataframes into dictionaries
tf_to_family = tf_info_df.set_index("transcription_factor")["family"].to_dict()
motif_to_family = motif_df.groupby("motifs")["family"].unique().apply(list).to_dict()
motif_to_tf = motif_df.groupby("motifs")["transcription_factor"].unique().apply(list).to_dict()

# read and filter the enricher output
with open(OUTPUT_FILE, 'w') as output_file:
    with open(ENRICHED_MOTIFS_FILE, 'r') as enriched_motif_file:
        for line in enriched_motif_file:
            if line.startswith("#") or line.startswith("set_id"):
                pass # skip the comments and the header
            else:
                transcription_factor, motif, _, _, _, _, _, _, _, target_genes = line.strip().split("\t")
                family = tf_to_family[transcription_factor]
                # without TF family extension (never extend if family is "unknown")
                if TF_FAMILY_EXTENSION == "TF_motifs" or family == "Unknown":
                    keep = transcription_factor in motif_to_tf[motif]
                # with TF family extension
                elif TF_FAMILY_EXTENSION == "TF-F_motifs":
                    keep = transcription_factor in tf_to_family and family in motif_to_family[motif]
                # if passed filter, write to output file
                if keep:
                    for target_gene in target_genes.split(","):
                        output_file.write(f"{transcription_factor}\t{target_gene}\n")
