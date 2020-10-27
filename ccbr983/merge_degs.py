import pandas
from functools import reduce

def readdf(filename,prefix):
	data=pandas.read_csv(filename,sep="\t",header=0,usecols=["ensid_gene","fc","log2fc","pvalue","fdr"])
	data=data[["ensid_gene","fc","log2fc","pvalue","fdr"]]
	data.columns=["ensid_gene",prefix+"_fc",prefix+"_log2fc",prefix+"_pvalue",prefix+"_fdr"]
	return data

gi=readdf('./GIT_vs_GIN/DESeq2_DEG_GI_T-GI_N_patient_covariate_all_genes.v2.txt',"GI")
skin=readdf('./Skin_T_vs_Skin_N/DESeq2_DEG_Skin_T-Skin_N_patient_covariate_all_genes.v2.txt',"Skin")
skingi=readdf('./T_vs_N/DESeq2_DEG_T-N_patient_covariate_all_genes.v2.txt',"SkinGI")

mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="ensid_gene"),[gi,skin,skingi])
mergeddf.fillna(' ',inplace=True)
mergeddf.drop_duplicates(inplace=True)

mergeddf.to_csv("merged_DEGs.v2.txt",sep="\t",index=False)
