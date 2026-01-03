import pandas as pd
import numpy as np
import wquantiles as wq # weighted-median


class ArmWeightedMedian:
    """Responsible for handling chromosome arm weighted median"""
    def __init__(self, cnr, centm_sorted):
        self.cnr = cnr
        self.centm_sorted = centm_sorted
    
    def chr_arm_weighted_median(self):
        """Calculate weighted median for each arm """
        # Sometimes, the cnr file has contigs, but centromere-bed only has standard chr info, so to only keep standard chr info in cnr
        # Defined standard human chromosomes 1-22 and X (as 23)
        chrs_standard = [str(i) for i in range(1, 23)] + ['23']
        self.cnr['chromosome'] = self.cnr['chromosome'].astype(str)
        chrs = self.cnr['chromosome'].unique()
        chrs_list = [c for c in chrs if c in chrs_standard]

        # Initialize empty dataframe for later use
        df_wmedian = pd.DataFrame(columns=self.cnr.columns).astype(self.cnr.dtypes) # to store weighted-median
        df_arms = pd.DataFrame(columns=["chromosome", "p_Genes", "q_Genes", "pArmEnds", "CM_row_pos"])
        qArm=[]
        for ch in chrs_list:
            cnr1 = self.cnr[self.cnr['chromosome']==ch].copy()
            cnr1['log2'] = cnr1['log2'].astype(float)
            ctm1 = self.centm_sorted[self.centm_sorted['chromosome']==ch]
            df_arms, qArm = self._find_arms(ch, cnr1, ctm1, df_arms, qArm)
            df_wmedian=self._cal_arm_weighted_median(ch, cnr1, qArm,df_arms, df_wmedian)
        
        return df_arms, df_wmedian

    def _find_arms(self, ch, cnr1, ctm1, df_arms, qArm):
        """find the row position where the p-arm ends/q-arm starts"""
        for i, pos in enumerate(cnr1['start']):
            if pos > ctm1['start'].iloc[0]:
                pos_idx = i
                if pos_idx == 0:
                    qArm.append(ch)

                p_genes_count = pos_idx
                p_arm_ends = cnr1.iloc[pos_idx-1]['end']
                q_genes_count = len(cnr1)-pos_idx
                cm_row_pos_idx = pos_idx
                new_row_arms = pd.DataFrame([{
                    "chromosome":ch,
                    "p_Genes":p_genes_count,
                    "q_Genes":q_genes_count,
                    "pArmEnds":p_arm_ends,
                    "CM_row_pos":cm_row_pos_idx
                }])
                df_arms = pd.concat([df_arms, new_row_arms], ignore_index=True)
                break
        return df_arms, qArm

    def _cal_arm_weighted_median(self, ch, cnr1, qArm, df_arms, df_wmedian):
        if ch in qArm:
            dat = np.array(cnr1['log2'])
            w = np.array(cnr1['weight'])
            m = wq.median(dat,w)
            cnr1['log2'] = m
            df_wmedian = pd.concat([df_wmedian, cnr1],ignore_index=True)
        else: # chromosome with both arms
            p_genes = df_arms.loc[df_arms['chromosome']==ch,'p_Genes'].item()
            # Process p-arms
            dat1 = np.array(cnr1.iloc[0:p_genes]['log2'])
            w1 = np.array(cnr1.iloc[0:p_genes]['weight'])
            m1 = wq.median(dat1, w1)
            cnr1.iloc[0:p_genes, cnr1.columns.get_loc('log2')] = m1
            # Process q-arms
            dat2 = np.array(cnr1.iloc[p_genes:]['log2'])
            w2 = np.array(cnr1.iloc[p_genes:]['weight'])
            m2 = wq.median(dat2, w2)
            cnr1.iloc[p_genes:, cnr1.columns.get_loc('log2')] = m2
            df_wmedian = pd.concat([df_wmedian, cnr1],ignore_index=True)
        return df_wmedian


