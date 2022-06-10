#!/usr/bin/env python

import pandas as pd

KMM1 = pd.read_csv("KMM1-A_SNP-InDel.csv")
#KMM1[0:5]

KMM3 = pd.read_csv("KMM3-A_Snp-InDel.csv")
KMM4 = pd.read_csv("KMM4_Snp-InDel.csv")

# Remove KMM1 background from KMM3

merged_KMM1KMM3 = KMM1.merge(KMM3, on =['Chromosome','Region','Allele'],indicator=True,how='outer')

merged_KMM1KMM3_II = merged_KMM1KMM3[merged_KMM1KMM3['_merge'] == 'right_only']

merged_KMM1KMM3_III = merged_KMM1KMM3_II[['Chromosome', 'Region','Type_y','Reference_y', 'Allele', 'Reference allele_y', 'Length_y', 'Linkage_y','Zygosity_y', 'Count_y', 'Coverage_y', 'Frequency_y','Probability_y', 'Forward read count_y', 'Reverse read count_y','Forward/reverse balance_y', 'Average quality_y', 'Read count_y','Read coverage_y', '# unique start positions_y','# unique end positions_y', 'Read position test probability_y','Read direction test probability_y', 'Hyper-allelic_y','BaseQRankSum_y', 'Homopolymer_y', 'Homopolymer length_y']]

merged_KMM1KMM3_III.to_csv('KMM3_A_Snp_InDel_minusKMM1.csv',index=False)

merged_KMM1KMM3_IV = merged_KMM1KMM3_III[merged_KMM1KMM3_III['Reference_y'] != merged_KMM1KMM3_III['Allele']]

merged_KMM1KMM3_IV.to_csv('KMM3_A_Snp_InDel_minusKMM1_noMNVwithSame.csv',index=False)

# Remove KMM1 background from KMM4

merged_KMM1KMM4 = KMM1.merge(KMM4, on =['Chromosome','Region','Allele'],indicator=True,how='outer')
#merged_KMM1KMM4[0:6]

merged_KMM1KMM4_II = merged_KMM1KMM4[merged_KMM1KMM4['_merge'] == 'right_only']

merged_KMM1KMM4_III = merged_KMM1KMM4_II[['Chromosome', 'Region','Type_y','Reference_y', 'Allele', 'Reference allele_y', 'Length_y', 'Linkage_y','Zygosity_y', 'Count_y', 'Coverage_y', 'Frequency_y','Probability_y', 'Forward read count_y', 'Reverse read count_y','Forward/reverse balance_y', 'Average quality_y', 'Read count_y','Read coverage_y', '# unique start positions_y','# unique end positions_y', 'Read position test probability_y','Read direction test probability_y', 'Hyper-allelic_y','BaseQRankSum_y', 'Homopolymer_y', 'Homopolymer length_y']]

merged_KMM1KMM4_III.to_csv('KMM4_A_Snp_InDel_minusKMM1.csv',index=False)

merged_KMM1KMM4_IV = merged_KMM1KMM4_III[merged_KMM1KMM4_III['Reference_y'] != merged_KMM1KMM4_III['Allele']]

merged_KMM1KMM4_IV.to_csv('KMM4_A_Snp_InDel_minusKMM1_noMNVwithSame.csv',index=False)


# Common SNPs between KMM3clean to KMM4clean
merged_KMM3KMM4 = merged_KMM1KMM3_III.merge(merged_KMM1KMM4_III, on =['Chromosome','Region','Allele'],indicator=True,how='outer')

#merged_KMM3KMM4[0:6]

merged_KMM3KMM4_II = merged_KMM3KMM4[merged_KMM3KMM4['_merge'] == 'both']

merged_KMM3KMM4_III = merged_KMM3KMM4_II[['Chromosome', 'Region','Type_y_y','Reference_y_y', 'Allele', 'Reference allele_y_y', 'Length_y_y', 'Linkage_y_y','Zygosity_y_y', 'Count_y_y', 'Coverage_y_y', 'Frequency_y_y','Probability_y_y', 'Forward read count_y_y', 'Reverse read count_y_y','Forward/reverse balance_y_y', 'Average quality_y_y', 'Read count_y_y','Read coverage_y_y', '# unique start positions_y_y','# unique end positions_y_y', 'Read position test probability_y_y','Read direction test probability_y_y', 'Hyper-allelic_y_y','BaseQRankSum_y_y', 'Homopolymer_y_y', 'Homopolymer length_y_y']]
merged_KMM3KMM4_III.to_csv('KMM3KMM4_A_Snp_InDel_shared.csv',index=False)


#KMM3 with NCBI Gene Key
KMM3_cleaned = pd.read_csv("KMM3_A_Snp_InDel_minusKMM1_noMNVwithSame_a.csv")
NCBI_info = pd.read_csv('Scenedesmus_NCBI_Info_a.txt', index_col=False, sep='\t')
NCBI_info_cut.columns = ['seqname','Chromosome']
KMM3_cleaned_wNamesKey = KMM3_cleaned.merge(NCBI_info_cut, on='Chromosome')

KMM3_cleaned_wNamesKey_II = KMM3_cleaned_wNamesKey

KMM3_cleaned_wNamesKey_II['Region'] = KMM3_cleaned_wNamesKey_II.Region.map(lambda x: x[0: x.find('..')] if '..' in x else x)

KMM3_cleaned_wNamesKey_II['Region'] = KMM3_cleaned_wNamesKey_II.Region.map(lambda x: x[0: x.find('^')] if '^' in x else x)

#Import gff

Sobliq_annot = pd.read_csv('6thrun_+AugPreds_Sobliq_maker_1-21-16.gff',index_col=False, skiprows=1,sep='\t',names=['seqname','source','feature','start','end','score','strand','frame','attribute'])

Sobliq_annot_CDS = Sobliq_annot[Sobliq_annot['feature'] == 'CDS']

KMM3_cleaned_wNamesKey_II_annot_CDS = KMM3_cleaned_wNamesKey_II.merge(Sobliq_annot_CDS, on='seqname')
KMM3_cleaned_wNamesKey_II_annot_CDS_II[['Region']] = KMM3_cleaned_wNamesKey_II_annot_CDS_II[['Region']].astype(float)

KMM3_cleaned_wNamesKey_II_annot_CDS_III = KMM3_cleaned_wNamesKey_II_annot_CDS_II[(KMM3_cleaned_wNamesKey_II_annot_CDS_II['Region'] >= KMM3_cleaned_wNamesKey_II_annot_CDS_II['start']) & (KMM3_cleaned_wNamesKey_II_annot_CDS_II['Region'] <= KMM3_cleaned_wNamesKey_II_annot_CDS_II['end'])]

KMM4_cleaned_wNamesKey_II_annot_CDS_IV = KMM4_cleaned_wNamesKey_II_annot_CDS_III.drop_duplicates(subset='attribute',keep='first')
