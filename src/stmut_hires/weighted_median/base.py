from abc import ABC, abstractmethod
import pandas as pd

class BaseWeightedMedian(ABC):

    def __init__(self, cnr, centm_sorted):
        self.cnr = cnr.copy()
        self.centm_sorted = centm_sorted

        self.cnr["chromosome"] = (self.cnr["chromosome"].astype(str))

    @abstractmethod
    def _calculate(self, ch, cnr1, oneArm, df_arms, df_wmedian):
        pass

    # cal weighted median 
    def run_smoother(self):
        """ 
        Shared orchestration logic -- Calculate weighted median.
        """
        # Sometimes, the cnr file has contigs, but centromere-bed only has standard chr info, so to only keep standard chr info in cnr
        # Defined standard human chromosomes 1-22 and X (as 23)
    
        chrs_standard = [str(i) for i in range(1, 23)] + ['23']
        self.cnr['chromosome'] = self.cnr['chromosome'].astype(str)
        chrs = self.cnr['chromosome'].unique()
        chrs_list = [c for c in chrs if c in chrs_standard]

        # Initialize empty dataframe to store smoothed wtmedian and gene counts in each arm per chr.
        df_wmedian = pd.DataFrame(columns=self.cnr.columns).astype(self.cnr.dtypes) # to store weighted-median
        df_arms = pd.DataFrame(columns=["chromosome", "p_Genes", "q_Genes", "pArmEnds", "CM_row_pos"])


        # Store single-arm only data
        oneArm=[] # store pArm-only, qArm-only chromosome 

        for ch in chrs_list:
            cnr1 = self.cnr[self.cnr['chromosome']==ch].copy()
            cnr1['log2'] = cnr1['log2'].astype(float)
            ctm1 = self.centm_sorted[self.centm_sorted['chromosome']==ch]
            
            if len(ctm1) == 0:
                print(f"[WARN] No centromere info for chromosome {ch}, skipping", flush=True)
                continue

            df_arms, oneArm = self._find_arms(ch, cnr1, ctm1, df_arms, oneArm)
            df_wmedian=self._calculate(ch, cnr1, oneArm, df_arms, df_wmedian)
        
        return df_arms, df_wmedian

    def _find_arms(self, ch, cnr1, ctm1, df_arms, oneArm):
        """Count gene numbers per arm """
        if cnr1["end"].values[-1] < ctm1["start"].values[0]:
            oneArm.append(ch)
            new_row_arms = pd.DataFrame([{
            "chromosome": ch,
            "p_Genes": len(cnr1),
            "q_Genes": 0,
            "pArmEnds": cnr1.iloc[-1]['end'],
            "CM_row_pos": len(cnr1)
            }])
            df_arms = pd.concat([df_arms, new_row_arms], ignore_index=True)

        # q-arm only
        elif cnr1["start"].values[0] > ctm1["end"].values[0]:
            oneArm.append(ch)
            new_row_arms = pd.DataFrame([{
            "chromosome": ch,
            "p_Genes": 0,
            "q_Genes": len(cnr1),
            "pArmEnds": ctm1["start"].values[0],
            "CM_row_pos": 0
            }])
            df_arms = pd.concat([df_arms, new_row_arms], ignore_index=True)

        else: # both arms
            centromere_start = ctm1['start'].iloc[0]
            for i, pos in enumerate(cnr1['start']):
                if pos > centromere_start:
                    pos_idx = i
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
        return df_arms, oneArm














