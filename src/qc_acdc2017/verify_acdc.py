"""Independent file-based NIfTI reader and exhaustive image/pair/index audit."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import traceback
import nibabel as nib
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
SOURCE=ROOT/'benchmarks/competitions/acdc2017'
WORK=ROOT/'_processing/all_data_cleaning_20261111/ACDC2017'
OUT=ROOT/'processed/ACDC2017_qc_20261111'
CHECKS=[]


def check(name,result):
    CHECKS.append({'name':name,'passed':bool(result)})
    if not result: raise AssertionError(name)


def digest(path,kind='sha256'):
    h=hashlib.new(kind)
    with path.open('rb') as f:
        while chunk:=f.read(4*1024*1024): h.update(chunk)
    return h.hexdigest()


def table(name):
    return pd.read_csv(OUT/name,sep='\t',keep_default_na=False,float_precision='round_trip')


def equal_rows(a,b):
    return a.columns.tolist()==b.columns.tolist() and sorted(map(tuple,a.astype(str).to_numpy()))==sorted(map(tuple,b.astype(str).to_numpy()))


def main():
    assert not (OUT/'FINAL_COMPLETION_GATE.json').exists(),'Existing successful verification must not be overwritten'
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    original=table('ORIGINAL_FILES.tsv').set_index('source_path')
    check('all_904_downloaded_files_accounted',len(baseline)==len(original)==904 and original.index.is_unique)
    check('no_missing_or_extra_source_files',set(original.index)=={str(p) for p in SOURCE.rglob('*') if p.is_file()})
    for item in baseline:
        path=Path(item['server_path'])
        row=original.loc[str(path)]
        name=str(path.relative_to(SOURCE))
        check(name+'_source_bytes',path.stat().st_size==int(item['bytes'])==row.bytes)
        check(name+'_source_MD5',digest(path,'md5')==item['md5']==row.MD5)
        check(name+'_source_SHA256',digest(path)==item['sha256']==row.SHA256)
    print('All 904 originals independently rehashed.',flush=True)
    patients=table('PATIENT_QC.tsv').set_index('patient_id')
    pairs=table('IMAGE_MASK_PAIR_QC.tsv')
    file_qc=table('NIFTI_FILE_QC.tsv').set_index('source_path')
    check('exact_150_patients_300_pairs_750_volumes',len(patients)==150 and patients.index.is_unique and len(pairs)==300 and not pairs[['patient_id','phase']].duplicated().any() and len(file_qc)==750 and file_qc.index.is_unique)
    cfg_files=sorted(SOURCE.glob('*/patient*/Info.cfg'))
    check('patient_index_covers_all_source_cases',set(patients.index)=={p.parent.name for p in cfg_files})
    check('file_index_covers_all_NIfTI_images',set(file_qc.index)=={str(p) for p in SOURCE.rglob('*.nii.gz')})
    expected_patient_good={}
    seen_pairs=[]
    for config_path in cfg_files:
        patient=config_path.parent.name
        split=config_path.parent.parent.name
        config=dict(line.strip().split(': ',1) for line in config_path.read_text().splitlines() if line.strip())
        pr=patients.loc[patient]
        check(patient+'_original_diagnosis_and_split',pr.official_split==split and pr.diagnostic_group==config['Group'])
        check(patient+'_original_ED_ES_NbFrame',pr.ED_one_based==int(config['ED']) and pr.ES_one_based==int(config['ES']) and pr.NbFrame==int(config['NbFrame']))
        check(patient+'_original_height_weight',float(pr.Height_source_value)==float(config['Height']) and float(pr.Weight_source_value)==float(config['Weight']))
        check(patient+'_original_config_path',pr.info_path==str(config_path))
        loaded={}
        for image_path in sorted(config_path.parent.glob('*.nii.gz')):
            # Different path from preparation: NiBabel opens the original gzip file.
            image=nib.load(image_path,mmap=False,keep_file_open=False)
            values=np.asarray(image.dataobj)
            recorded=file_qc.loc[str(image_path)]
            kind='cine4d' if image_path.name.endswith('_4d.nii.gz') else 'segmentation' if image_path.name.endswith('_gt.nii.gz') else 'phase_image'
            check(image_path.name+'_native_shape_dtype_and_kind',json.loads(recorded['shape'])==list(values.shape) and recorded['dtype']==str(values.dtype) and recorded.kind==kind)
            check(image_path.name+'_native_affine_zooms',np.array_equal(np.asarray(json.loads(recorded.affine)),image.affine) and tuple(json.loads(recorded.zooms))==tuple(image.header.get_zooms()))
            check(image_path.name+'_native_coordinate_codes',recorded.qform_code==int(image.header['qform_code']) and recorded.sform_code==int(image.header['sform_code']) and recorded.axis_codes==''.join(nib.aff2axcodes(image.affine)))
            check(image_path.name+'_every_pixel_exact',recorded.pixel_SHA256==hashlib.sha256(np.ascontiguousarray(values)).hexdigest())
            finite=np.isfinite(values)
            finite_values=values[finite]
            check(image_path.name+'_nonfinite_pixel_count',recorded.nonfinite_pixels==int((~finite).sum()))
            check(image_path.name+'_minimum_maximum',recorded.minimum==float(finite_values.min()) and recorded.maximum==float(finite_values.max()))
            check(image_path.name+'_patient_split_group',recorded.patient_id==patient and recorded.official_split==split and recorded.diagnostic_group==config['Group'])
            reasons=[]
            if values.ndim!=(4 if kind=='cine4d' else 3): reasons.append('wrong_dimensions')
            if not finite.all(): reasons.append('nonfinite_pixels')
            if not np.isfinite(image.affine).all() or abs(np.linalg.det(image.affine[:3,:3]))<1e-9: reasons.append('invalid_affine')
            if not all(np.isfinite(z) and z>0 for z in image.header.get_zooms()): reasons.append('invalid_zooms')
            if not np.any(values!=0): reasons.append('all_zero_image_or_mask')
            labels=np.unique(values).tolist() if kind=='segmentation' else []
            if kind=='segmentation' and (not np.equal(values,np.floor(values)).all() or not set(labels)<={0,1,2,3}): reasons.append('invalid_segmentation_label')
            check(image_path.name+'_segmentation_labels',json.loads(recorded.label_values)==labels)
            check(image_path.name+'_file_eligibility_and_reasons',recorded.usable_file==(not reasons) and recorded.file_reasons=='|'.join(reasons))
            check(image_path.name+'_full_gzip_CRC_was_verified',recorded.gzip_CRC_verified)
            loaded[image_path.name]=(image,values,not reasons)
        cine_path=config_path.parent/(patient+'_4d.nii.gz')
        cine,cine_values,cine_valid=loaded[cine_path.name]
        n_frames=int(config['NbFrame'])
        meta_ok=cine_values.ndim==4 and cine_values.shape[3]==n_frames and 1<=int(config['ED'])<=n_frames and 1<=int(config['ES'])<=n_frames
        patient_good=[]
        for phase in ['ED','ES']:
            frame=int(config[phase])
            image_path=config_path.parent/f'{patient}_frame{frame:02d}.nii.gz'
            mask_path=config_path.parent/f'{patient}_frame{frame:02d}_gt.nii.gz'
            image,values,image_valid=loaded[image_path.name]
            mask,labels,mask_valid=loaded[mask_path.name]
            recorded=pairs.loc[pairs.patient_id.eq(patient)&pairs.phase.eq(phase)]
            check(patient+'_'+phase+'_one_pair_row',len(recorded)==1)
            row=recorded.iloc[0]
            check(patient+'_'+phase+'_exact_image_mask_cine_paths',row.image_path==str(image_path) and row.mask_path==str(mask_path) and row.cine4d_path==str(cine_path))
            check(patient+'_'+phase+'_source_phase_numbering',row.frame_number_one_based==frame and row.frame_index_zero_based==frame-1)
            shape_ok=values.shape==labels.shape==cine_values.shape[:3]
            image_mask=np.allclose(image.affine,mask.affine,rtol=1e-6,atol=1e-5)
            image_cine=np.allclose(image.affine,cine.affine,rtol=1e-6,atol=1e-5)
            exact_frame=bool(meta_ok and shape_ok and np.array_equal(values,cine_values[:,:,:,frame-1]))
            check(patient+'_'+phase+'_shape_and_spatial_matches',row.shape_matches==shape_ok and row.image_mask_affine_matches==image_mask and row.image_cine_affine_matches==image_cine)
            check(patient+'_'+phase+'_all_pixels_compared_to_4d_phase',row.phase_exactly_matches_4d_frame==exact_frame)
            check(patient+'_'+phase+'_image_mask_pixel_fingerprints',row.image_pixel_SHA256==hashlib.sha256(np.ascontiguousarray(values)).hexdigest() and row.mask_pixel_SHA256==hashlib.sha256(np.ascontiguousarray(labels)).hexdigest())
            check(patient+'_'+phase+'_each_segmentation_label_voxel_count',all(row[key]==int(np.count_nonzero(labels==label)) for key,label in [('background_voxels',0),('RV_voxels',1),('MYO_voxels',2),('LV_voxels',3)]))
            reasons=[]
            if not image_valid: reasons.append('invalid_phase_image')
            if not mask_valid: reasons.append('invalid_segmentation')
            if not cine_valid: reasons.append('invalid_4d_sequence')
            if not meta_ok: reasons.append('frame_metadata_conflict')
            if not shape_ok: reasons.append('shape_mismatch')
            if not image_mask or not image_cine: reasons.append('affine_mismatch')
            if not exact_frame: reasons.append('phase_values_differ_from_4d_frame')
            good=not reasons
            check(patient+'_'+phase+'_eligibility_and_quarantine_reasons',row.usable_pair==good and row.quarantine_reason=='|'.join(reasons))
            check(patient+'_'+phase+'_diagnosis_and_split_preserved',row.official_split==split and row.diagnostic_group==config['Group'])
            patient_good.append(good)
            seen_pairs.append((patient,phase))
        expected_patient_good[patient]=all(patient_good)
        check(patient+'_both_phase_patient_eligibility',pr.usable_both_phases==all(patient_good) and pr.usable_phase_pairs==sum(patient_good))
        print('Independently verified case:',patient,split,'usable phases',sum(patient_good),flush=True)
    check('all_300_phase_pairs_covered_once',len(seen_pairs)==len(set(seen_pairs))==300)
    check('training_and_testing_patients_unchanged',set(patients.loc[patients.official_split.eq('training')].index)=={f'patient{i:03d}' for i in range(1,101)} and set(patients.loc[patients.official_split.eq('testing')].index)=={f'patient{i:03d}' for i in range(101,151)})
    check('exact_usable_pair_index',equal_rows(table('USABLE_IMAGE_MASK_PAIRS.tsv'),pairs.loc[pairs.usable_pair]))
    check('exact_quarantine_pair_index',equal_rows(table('QUARANTINED_IMAGE_MASK_PAIRS.tsv'),pairs.loc[~pairs.usable_pair]))
    expected=patients.reset_index().groupby(['official_split','diagnostic_group']).agg(patients=('patient_id','size'),usable_both_phases=('usable_both_phases','sum'),usable_phase_pairs=('usable_phase_pairs','sum')).reset_index()
    check('exact_group_and_split_summary',equal_rows(table('GROUP_SPLIT_SUMMARY.tsv'),expected))
    policy=json.loads((OUT/'qc_policy.json').read_text())
    check('no_pixel_or_training_operations',all(policy[key] is False for key in ['normalization','resampling','new_segmentation','intensity_changes','model_training']) and policy['read_only_originals'])
    check('original_license_and_mandatory_citation',all((OUT/name).read_bytes()==(SOURCE/name).read_bytes() for name in ['LICENSE_TERMS.md','MANDATORY_CITATION.md']))
    for row in table('OUTPUT_MANIFEST.tsv').itertuples(index=False):
        path=Path(row.server_path)
        check('output_bytes_'+row.file,path.stat().st_size==row.bytes)
        check('output_MD5_'+row.file,digest(path,'md5')==row.MD5)
        check('output_SHA256_'+row.file,digest(path)==row.SHA256)
    now=datetime.now(timezone.utc).isoformat()
    verification={'status':'passed','checked_utc':now,'checks_passed':len(CHECKS),'source_files':904,'patients':150,'NIfTI_files_fully_read':750,'phase_pairs':300,'usable_pairs':int(pairs.usable_pair.sum()),'quarantined_pairs':int((~pairs.usable_pair).sum()),'checks':CHECKS}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
    gate={'status':'passed','checked_utc':now,'independent_checks':len(CHECKS),'verification_SHA256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),'prepared_manifest_SHA256':digest(OUT/'OUTPUT_MANIFEST.tsv'),'source_immutable':True,'scope':'all 904 source files, all 750 NIfTI volumes and 300 image-mask phase pairs, exact official patient splits'}
    (OUT/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    records=[]
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name!='VERIFIED_OUTPUT_MANIFEST.tsv': records.append({'file':path.name,'server_path':str(path),'bytes':path.stat().st_size,'MD5':digest(path,'md5'),'SHA256':digest(path)})
    pd.DataFrame(records).to_csv(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    print('All',len(CHECKS),'independent checks passed.',flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as error:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/'VERIFICATION_FAILURE.json').write_text(json.dumps({'error':repr(error),'traceback':traceback.format_exc(),'checks':CHECKS},indent=2)+'\n')
        raise
