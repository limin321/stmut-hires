import pandas as pd
import numpy as np
import wquantiles as wq # weighted-median
from .base import BaseWeightedMedian


class ArmWeightedMedian(BaseWeightedMedian): 
    """Responsible for handling chromosome arm weighted median"""
    def __init__(self, cnr, centm_sorted, **kwargs):
        super().__init__(cnr, centm_sorted)
    
    def _calculate(self, ch, cnr1, oneArm, df_arms, df_wmedian):
        # Case1: q-arm only or p-arm only
        if ch in oneArm:
            dat = cnr1['log2'].to_numpy()
            w = cnr1['weight'].to_numpy()

            m = wq.median(dat,w)
            cnr1['log2'] = m
            df_wmedian = pd.concat([df_wmedian, cnr1],ignore_index=True)
                
        else: 
            # Case 2 Both arms present
            matches = df_arms.loc[df_arms['chromosome'] == ch, 'p_Genes']
            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly 1 arm row for chromosome {ch!r}, got {len(matches)}. "
                    f"df_arms chromosomes: {df_arms['chromosome'].tolist()}"
                )
            p_genes = int(matches.iloc[0])

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
