import pandas
from functools import reduce

gi=pandas.read_csv("DESeq2_DEG_GI_T-GI_N_all_genes.txt",sep="\t",header=0,usecols=["gene","log2fc","fdr"])
gi.columns=["gene","GI.log2fc","GI.fdr"]
skin=pandas.read_csv("DESeq2_DEG_Skin_T-Skin_N_all_genes.txt",sep="\t",header=0,usecols=["gene","log2fc"
,"fdr"])
skin.columns=["gene","Skin.log2fc","Skin.fdr"]

mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="gene"),[gi,skin])
mergeddf.fillna(' ',inplace=True)
mergeddf.drop_duplicates(inplace=True)

mergeddf.to_csv("iDEP_input.txt",sep="\t",index=False)
