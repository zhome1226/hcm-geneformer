"""Summarize overlapping QC review reasons without changing any membership."""
from pathlib import Path
import pandas as pd
root=Path('/home/zhome/hcm-geneformer-data/processed')
for ds in ['GSE255296','GSE292754']:
    out=root/(ds+'_qc_20261114')
    x=pd.read_csv(out/'cell_identity_audit.tsv.gz',sep='\t',index_col=0,low_memory=False).fillna({'reference_pruned':''})
    rows=[]
    for sample,group in x.groupby('sample_id'):
        pending=group.loc[group.review_action.eq('pending_CM')]
        reasons={
            'pending_total':pd.Series(True,index=pending.index),
            'raw_SingleR_not_CM':pending.reference_label.ne('Cardiomyocyte'),
            'pruned_SingleR_not_CM':pending.reference_pruned.ne('Cardiomyocyte'),
            'fewer_than_2_CM_markers':pending.Cardiomyocyte_markers.lt(2),
            'CM_program_not_dominant':~pending.cm_program_dominates,
            'reference_tail_and_context_fail':~(pending.reference_tail_pass | (pending.Cardiomyocyte_markers.ge(3)&pending.cm_anchor_neighbor_fraction.ge(.5))),
            'strong_mixed_lineage_flag':pending.coherent_mixed_panel_flag,
            'algorithmic_doublet_call':pending.scDblFinder_class.eq('doublet'),
        }
        for reason,mask in reasons.items():rows.append({'dataset':ds,'sample_id':sample,'reason_nonexclusive':reason,'count':int(mask.sum())})
    pd.DataFrame(rows).to_csv(out/'PENDING_REASON_SUMMARY.tsv',sep='\t',index=False)
    print(ds,pd.DataFrame(rows).groupby('reason_nonexclusive')['count'].sum().to_dict(),flush=True)
