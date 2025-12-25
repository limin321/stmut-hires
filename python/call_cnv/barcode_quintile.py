import os
import pandas as pd

class AssignBarcodeQuintile:
    """Assign quintile values (Q1...Q5, out) to individual barcode based on ranked CNVScore"""
    @staticmethod
    def barcode_quintils(merged_df,caseCNVscore_FDR):
        """Assign quintile value to each barcode."""
        merged_df_sub = merged_df[['main_barcode', 'merged_barcodes']].rename(
            columns={
                'main_barcode':'barcode',
                'merged_barcodes': 'merged_barcodes'
            }
        )

        barcode_in_scores = merged_df_sub['barcode'].isin(caseCNVscore_FDR['barcode'])
        outs=merged_df_sub[~barcode_in_scores].copy()
        outs['quintile'] = 'out'
        outs

        keep = merged_df_sub[barcode_in_scores]
        keep_quintile = pd.merge(keep, caseCNVscore_FDR[['barcode', 'CNVScore','Rank']].copy(), left_on='barcode', right_on='barcode')
        keep_quintile = keep_quintile.sort_values(by='Rank')
        keep_quintile

        # calculate quintiles labels
        a = round(len(keep_quintile)/5)
        if not keep_quintile.empty:
            bins = [0, a, 2*a, 3*a, 4*a, len(keep_quintile)]
            labels = ["Q1","Q2","Q3","Q4","Q5"]
            keep_quintile['quintile'] = pd.cut(keep_quintile.reset_index(drop=True).index, bins=bins, labels=labels, right=False, include_lowest=True)
        else:
            keep_quintile['quintile'] = pd.NA


        keep_quintile = keep_quintile[['barcode', 'quintile']]
        keep_quintile

        outs = outs[['barcode', 'quintile']]
        all_combined = pd.concat([keep_quintile,outs])
        all_barcodes_quintile = AssignBarcodeQuintile._expand_barcode(merged_df_sub,all_combined)
        return all_barcodes_quintile

    @staticmethod
    def _expand_barcode(merged_df_sub,all_combined):
        """Expand merged_barcodes and assign quintile"""
        merged_expand = (
            merged_df_sub
            .rename(columns={'barcode':'parent_barcode'})
            .assign(barcode = lambda df: df['merged_barcodes'].str.split(","))
            .explode("barcode")
        )
        merged_expand['barcode'] = merged_expand['barcode'].str.strip()

        final_quintile = (
            merged_expand
            .merge(
                all_combined,
                left_on = 'parent_barcode',
                right_on = 'barcode',
                how = 'left'
            )
            [['barcode_x', 'quintile']]
            .rename(columns={'barcode_x':'barcode'})
        )
        return final_quintile








