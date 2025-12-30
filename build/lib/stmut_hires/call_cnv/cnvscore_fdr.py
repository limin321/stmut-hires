import pandas as pd
import numpy as np

class CNVScoreFDR:
    """
    Calculates observed & Permuted ranks and fractions for CNV scores in preparation for FDR analysis.
    Manages the calculation of statistical significance and False Discovery Rates (FDR) for CNV scores.
    """


    # Prepare CaseCNVscore for FDR
    @staticmethod
    def calculate_observed_fractions(caseCNVscore):
        """Calculates the proportions and ranks of observed data points."""
        caseCNVscore0 = caseCNVscore.copy()
        total_n = caseCNVscore0.shape[0]
        caseCNVscore0['Fraction_at_this_scor_or_higher_obs'] = caseCNVscore0['Rank']/total_n
        caseCNVscore0.sort_values(by='CNVScore', ascending=True, inplace=True)
        caseCNVscore0['Rank_reverse'] = caseCNVscore0['CNVScore'].rank(method='dense', ascending=True).astype(int)
        caseCNVscore0['Fraction_at_that_scor_or_lower_obs'] = caseCNVscore0['Rank_reverse']/total_n
        caseCNVscore0.sort_values(by='CNVScore', ascending=False, inplace=True)
        caseCNVscore0['mean_diff'] = abs(caseCNVscore0['CNVScore'] - caseCNVscore0['CNVScore'].mean())
        caseCNVscore0.sort_values(by='mean_diff', ascending=False, inplace=True)
        caseCNVscore0['mean_diff_rank'] = caseCNVscore0['mean_diff'].rank(method='dense', ascending=False).astype(int)
        caseCNVscore0['Fraction_at_this_diff_or_higher_obs'] = caseCNVscore0['mean_diff_rank']/total_n
        caseCNVscore0.sort_values(by='Rank', ascending=True, inplace=True)
        return caseCNVscore0


    @staticmethod
    def calculate_permuted_fractions(permuts_df):
        """Calculates the proportions and ranks of permuted scores."""
        flattened_arr = permuts_df.values.flatten()
        permuts_long = pd.DataFrame(flattened_arr, columns=['permutScore'])
        permuts_long.sort_values(by='permutScore', ascending=False, inplace=True)
        permuts_long['Rank'] = permuts_long['permutScore'].rank(method='dense', ascending=False).astype(int)

        nr = permuts_long.shape[0]
        permuts_long['Fraction_at_this_scor_or_higher_permut'] = permuts_long['Rank']/nr
        permuts_long.sort_values(by='permutScore', ascending=True, inplace=True)
        permuts_long['Rank_reverse'] = permuts_long['permutScore'].rank(method='dense', ascending=True).astype(int)
        permuts_long['Fraction_at_that_scor_or_lower_permut'] = permuts_long['Rank_reverse']/nr
        permuts_long.sort_values(by='permutScore', ascending=False, inplace=True)
        permuts_long['mean_diff'] = abs(permuts_long['permutScore'] - permuts_long['permutScore'].mean())
        permuts_long.sort_values(by='mean_diff', ascending=False, inplace=True)
        permuts_long['mean_diff_rank'] = permuts_long['mean_diff'].rank(method='dense', ascending=False).astype(int)
        permuts_long['Fraction_at_this_diff_or_higher_permut'] = permuts_long['mean_diff_rank']/nr
        permuts_long.sort_values(by='Rank', ascending=True, inplace=True)
        return permuts_long

    # Calculate FDR
    def calculate_observed_FDR(permuts_long, caseCNVscore):
        """Calcuate the FDR of case CNVScores"""
        # 1. Prepare data (use negative values to treat large scores as small values during sort)
        # and sort ascendingly based on these negative values
        permuts_long_sorted = permuts_long.copy()
        permuts_long_sorted['neg_score'] = -permuts_long_sorted['permutScore']
        permuts_long_sorted = permuts_long_sorted.sort_values(by='neg_score', ascending=True).reset_index(drop=True)
        
        case_scores_neg = -caseCNVscore['CNVScore'].values
        permut_scores_neg = permuts_long_sorted['neg_score'].values
        
        # 2. Perform binary search using the negative scores
        # We want to find the first index where the permuted score is less than the case score.
        # In negative space, this corresponds to finding the first index where 
        # the negative permuted score is GREATER than the negative case score ('left').
        indices = np.searchsorted(permut_scores_neg, case_scores_neg, side='left')
        
        # 3. Assign values using the found indices
        max_idx = len(permuts_long_sorted) - 1
        # Ensure indices don't go out of bounds (should rarely happen here if ranges overlap)
        safe_indices = np.clip(indices, 0, max_idx) 

        caseCNVscore['Fraction_at_this_scor_or_higher_permut'] = permuts_long_sorted.iloc[safe_indices]['Fraction_at_this_scor_or_higher_permut'].values
        caseCNVscore['Fraction_at_that_scor_or_lower_permut'] = permuts_long_sorted.iloc[safe_indices]['Fraction_at_that_scor_or_lower_permut'].values
        caseCNVscore['Fraction_at_this_diff_or_higher_permut'] = permuts_long_sorted.iloc[safe_indices]['Fraction_at_this_diff_or_higher_permut'].values
        caseCNVscore['FDR'] = caseCNVscore['Fraction_at_this_scor_or_higher_permut'] / caseCNVscore['Fraction_at_this_scor_or_higher_obs']

        return caseCNVscore
