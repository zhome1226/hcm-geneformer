"""Preserve strict geometry quarantine and validate a separate voxel-grid task view."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import nibabel as nib
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
OUT=ROOT/'processed/ACDC2017_qc_20261111'


def digest(path,kind='sha256'):
    with path.open('rb') as f: return hashlib.file_digest(f,kind).hexdigest()


def main():
    assert not (OUT/'GEOMETRY_COMPLETION_GATE.json').exists(),'Successful geometry view exists'
    gate=json.loads((OUT/'FINAL_COMPLETION_GATE.json').read_text())
    verification=json.loads((OUT/'INDEPENDENT_VERIFICATION.json').read_text())
    assert gate['status']==verification['status']=='passed'
    assert gate['verification_SHA256']==digest(OUT/'INDEPENDENT_VERIFICATION.json')
    for row in pd.read_csv(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
        p=Path(row.server_path)
        assert p.stat().st_size==row.bytes and digest(p,'md5')==row.MD5 and digest(p)==row.SHA256
    pairs=pd.read_csv(OUT/'IMAGE_MASK_PAIR_QC.tsv',sep='\t',keep_default_na=False)
    file_qc=pd.read_csv(OUT/'NIFTI_FILE_QC.tsv',sep='\t',keep_default_na=False).set_index('source_path')
    checks,rows,conflicts=[],[],[]
    for row in pairs.itertuples(index=False):
        image,mask=nib.load(row.image_path),nib.load(row.mask_path)
        a,b=np.asarray(image.dataobj),np.asarray(mask.dataobj)
        tests={'image_pixels_unchanged':hashlib.sha256(np.ascontiguousarray(a)).hexdigest()==row.image_pixel_SHA256,
            'mask_pixels_unchanged':hashlib.sha256(np.ascontiguousarray(b)).hexdigest()==row.mask_pixel_SHA256,
            '3d_image_mask_shape_matches':a.ndim==b.ndim==3 and a.shape==b.shape,
            'image_mask_native_affine_matches':np.allclose(image.affine,mask.affine,rtol=1e-6,atol=1e-5),
            'all_phase_image_and_mask_pixels_finite':np.isfinite(a).all() and np.isfinite(b).all(),
            'segmentation_labels_valid':np.equal(b,np.rint(b)).all() and set(np.unique(b))<={0,1,2,3} and np.any(b>0),
            'phase_link_to_4d_exact_pixels_already_independently_verified':bool(row.phase_exactly_matches_4d_frame)}
        for name,result in tests.items():
            checks.append({'name':row.patient_id+'_'+row.phase+'_'+name,'passed':bool(result)})
        # Source geometry disagreement remains a restriction for spatial/physical use.
        voxel_valid=all(tests.values())
        output=row._asdict()
        output['usable_voxel_grid_segmentation_pair']=voxel_valid
        output['cross_file_physical_geometry_consistent']=bool(row.usable_pair)
        output['physical_geometry_review_required']=not row.image_cine_affine_matches
        output['voxel_grid_scope']='index_aligned_3D_image_and_mask_only_no_cross_file_physical_coordinate_claim'
        rows.append(output)
        if not row.image_cine_affine_matches:
            conflict={'patient_id':row.patient_id,'phase':row.phase,'official_split':row.official_split,'diagnostic_group':row.diagnostic_group,
                'image_path':row.image_path,'mask_path':row.mask_path,'cine4d_path':row.cine4d_path,
                'image_mask_affine_matches':bool(row.image_mask_affine_matches),'phase_pixels_equal_4d_frame':bool(row.phase_exactly_matches_4d_frame),
                'phase_image_affine':file_qc.loc[row.image_path,'affine'],'mask_affine':file_qc.loc[row.mask_path,'affine'],'cine4d_affine':file_qc.loc[row.cine4d_path,'affine'],
                'action':'retain_native_headers; usable_for_voxel_grid_segmentation; physical_geometry_requires_source_review'}
            conflicts.append(conflict)
    data=pd.DataFrame(rows)
    data.to_csv(OUT/'VOXEL_GRID_PAIR_QC.tsv',sep='\t',index=False)
    data.loc[data.usable_voxel_grid_segmentation_pair].to_csv(OUT/'USABLE_VOXEL_GRID_IMAGE_MASK_PAIRS.tsv',sep='\t',index=False)
    data.loc[~data.usable_voxel_grid_segmentation_pair].to_csv(OUT/'QUARANTINED_VOXEL_GRID_PAIRS.tsv',sep='\t',index=False)
    pd.DataFrame(conflicts).to_csv(OUT/'GEOMETRY_CONFLICT_AUDIT.tsv',sep='\t',index=False)
    assert all(c['passed'] for c in checks),'Inspect failed source image/mask checks before certifying voxel-grid view'
    report=f"# ACDC2017 空间头冲突与用途索引\n\n原严格几何验收{verification['checks_passed']}项通过；本次再从原路径读取全部300对3D图像/标注，{len(checks)}项逐对检查通过。\n\n体素网格分割可用{int(data.usable_voxel_grid_segmentation_pair.sum())}对；3D/4D空间头一致{int(data.cross_file_physical_geometry_consistent.sum())}对；空间头冲突{len(conflicts)}对（{len(set(r['patient_id'] for r in conflicts))}位患者）。冲突病例的图像/标注彼此对齐且原像素精确等于4D对应帧，但原3D头空间与4D头不同。\n\nUSABLE_VOXEL_GRID_IMAGE_MASK_PAIRS用于原3D图像/标注体素网格上的分割任务；原USABLE_IMAGE_MASK_PAIRS仅含空间头一致对象。72条原严格隔离记录仍保留，不能把体素可用性解释为物理几何已经修复。GEOMETRY_CONFLICT_AUDIT给出每对原始仿射；不覆盖原头、不填造几何、不计算临床体积、不重采样。\n"
    (OUT/'GEOMETRY_POLICY.md').write_text(report)
    now=datetime.now(timezone.utc).isoformat()
    evidence={'status':'passed','checked_utc':now,'checks_passed':len(checks),'checks':checks,'strict_verification_SHA256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),
        'voxel_grid_usable_pairs':int(data.usable_voxel_grid_segmentation_pair.sum()),'strict_geometry_usable_pairs':int(data.cross_file_physical_geometry_consistent.sum()),'source_geometry_conflict_pairs':len(conflicts),'original_headers_modified':False}
    (OUT/'GEOMETRY_VIEWS_VERIFICATION.json').write_text(json.dumps(evidence,indent=2)+'\n')
    names=['VOXEL_GRID_PAIR_QC.tsv','USABLE_VOXEL_GRID_IMAGE_MASK_PAIRS.tsv','QUARANTINED_VOXEL_GRID_PAIRS.tsv','GEOMETRY_CONFLICT_AUDIT.tsv','GEOMETRY_POLICY.md','GEOMETRY_VIEWS_VERIFICATION.json']
    manifest=[]
    for name in names:
        p=OUT/name
        manifest.append({'file':name,'server_path':str(p),'bytes':p.stat().st_size,'MD5':digest(p,'md5'),'SHA256':digest(p)})
    pd.DataFrame(manifest).to_csv(OUT/'GEOMETRY_VIEWS_MANIFEST.tsv',sep='\t',index=False)
    completion={'status':'passed','checked_utc':now,'verification_SHA256':digest(OUT/'GEOMETRY_VIEWS_VERIFICATION.json'),'manifest_SHA256':digest(OUT/'GEOMETRY_VIEWS_MANIFEST.tsv'),'original_strict_gate_SHA256':digest(OUT/'FINAL_COMPLETION_GATE.json')}
    (OUT/'GEOMETRY_COMPLETION_GATE.json').write_text(json.dumps(completion,indent=2)+'\n')
    print('Geometry views:',evidence['voxel_grid_usable_pairs'],'voxel-grid pairs;',evidence['strict_geometry_usable_pairs'],'strict geometry;',len(checks),'checks passed.',flush=True)


if __name__=='__main__':
    main()
