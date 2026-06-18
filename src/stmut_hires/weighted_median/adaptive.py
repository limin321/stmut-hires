import pandas as pd
import numpy as np
import wquantiles as wq # weighted-median
from .base import BaseWeightedMedian


class AdaptiveWeightedMedian(BaseWeightedMedian): 
    """Calculate regional weighted median"""
    
    def __init__(
        self, 
        cnr, 
        centm_sorted,
        target_weight = 25
    ):

        # calling the parent constructor
        super().__init__(cnr, centm_sorted)

        self.target_weight = target_weight
    


    """  
    grows a genomically local adaptive window around each bin until sufficient cumulative confidence/coverage is reached, then uses a weighted median to robustly smooth the signal.
    """
    # Calculate weighted median per arm
    @staticmethod
    def _local_wtmedian(armdata, target_weight):
        """Calculate local weighted-median given provided target_weights """
        starts = armdata["start"].to_numpy()
        ends = armdata["end"].to_numpy()
        weights = armdata["weight"].to_numpy()
        values = armdata["log2"].to_numpy()
        n = len(armdata)
        smoothed = np.zeros(len(armdata))

        # Edge case, total arm weight below target - use arm-level weighted median
        if weights.sum() < target_weight:
            smoothed[:] = wq.median(values, weights)
            return smoothed

        for i in range(n):
            left, right = i, i
            cumulative_weight = weights[i]

            # expand outward until target weight reached
            while cumulative_weight < target_weight:
                can_expand_left = left > 0
                can_expand_right = right < (n - 1)

                # Only stop when neither side can expand anymore. If one-side can still expand, continue...
                if not can_expand_left and not can_expand_right:
                    break

                # Prefer genomically closer side
                left_distance = starts[i] - ends[left - 1] if can_expand_left else np.inf
                right_distance = starts[right + 1] - ends[i] if can_expand_right else np.inf

                if left_distance <= right_distance:
                    left -= 1
                    cumulative_weight += weights[left]

                else:
                    right += 1
                    cumulative_weight += weights[right]


            # local window
            local_values = values[left:right + 1]
            local_weights = weights[left:right + 1]
            smoothed[i] = wq.median(local_values, local_weights)
        
        return smoothed


    # Calculate adaptive weighted median 
    def _calculate(self, ch, cnr1, oneArm, df_arms, df_wmedian):
        """
        Adaptive local weighted-median smoothing.

        For each gene:
            - expand left/right symmetrically
            - accumulate neighboring weights
            - stop when cumulative weight >= target_weight
            - assign weighted median to focal gene
        """
        # Case1: q-arm only or p-arm only
        if ch in oneArm:
            print("working on single-arm")
            m=self._local_wtmedian(cnr1,self.target_weight)
            cnr1['log2'] = m
            df_wmedian = pd.concat([df_wmedian, cnr1],ignore_index=True)
                
        else: 
            print("working on both arm")
            # Case 2 Both arms present
            matches = df_arms.loc[df_arms['chromosome'] == ch, 'p_Genes']
            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly 1 arm row for chromosome {ch!r}, got {len(matches)}. "
                    f"df_arms chromosomes: {df_arms['chromosome'].tolist()}"
                )
            p_genes = int(matches.iloc[0])

            # Process p-arms
            print(f"p-genes is {p_genes}")

            parmdata=cnr1.iloc[0:p_genes]
            print(f"p-arm has shape: {parmdata.shape}")
            m1=self._local_wtmedian(parmdata,self.target_weight)
            cnr1.iloc[0:p_genes, cnr1.columns.get_loc('log2')] = m1
            
            # Process q-arms
            qarmdata=cnr1.iloc[p_genes:]
            print(f"q-arm has shape {qarmdata.shape}")
            m2=self._local_wtmedian(qarmdata,self.target_weight)
            cnr1.iloc[p_genes:, cnr1.columns.get_loc('log2')] = m2
            df_wmedian = pd.concat([df_wmedian, cnr1],ignore_index=True)
        return df_wmedian














