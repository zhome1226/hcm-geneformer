"""Validate immutable ACDC NIfTI images, annotations, cardiac phases and splits."""
from datetime import datetime,timezone
from pathlib import Path
import gzip
import hashlib
import importlib.metadata
import json
import shutil
import nibabel as nib
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
SOURCE=ROOT/'benchmarks/competitions/acdc2017'
WORK=ROOT/'_processing/all_data_cleaning_20261111/ACDC2017'
OUT=ROOT/'processed/ACDC2017_qc_20261111'


def digest(path,kind='sha256'):
    with path.open('rb') as f:
        return hashlib.file_digest(f,kind).hexdigest()


def pixels_digest(a):
    return hashlib.sha256(np.ascontiguousarray(a)).hexdigest()


def read_image(path,kind):
    # gzip.decompress validates the entire stream and trailer CRC before parsing.
    payload=gzip.decompress(path.read_bytes())
    img=nib.Nifti1Image.from_bytes(payload)
    data=np.asanyarray(img.dataobj)
    affine=img.affine
    zooms=img.header.get_zooms()
    finite=np.isfinite(data)
    values=data[finite]
    reasons=[]
    if data.ndim!=(4 if kind=='cine4d' else 3): reasons.append('wrong_dimensions')
    if not finite.all(): reasons.append('nonfinite_pixels')
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3,:3]))<1e-9: reasons.append('invalid_affine')
    if not all(np.isfinite(z) and z>0 for z in zooms): reasons.append('invalid_zooms')
    if not np.any(data!=0): reasons.append('all_zero_image_or_mask')
    labels=[]
    if kind=='segmentation':
        labels=np.unique(data).tolist()
        if not np.equal(data,np.rint(data)).all() or not set(labels)<={0,1,2,3}: reasons.append('invalid_segmentation_label')
    record={'source_path':str(path),'kind':kind,'shape':json.dumps(list(data.shape)),'dtype':str(data.dtype),
        'zooms':json.dumps(list(map(float,zooms))),'affine':json.dumps(affine.tolist()),
        'qform_code':int(img.header['qform_code']),'sform_code':int(img.header['sform_code']),
        'axis_codes':''.join(nib.aff2axcodes(affine)),'nonfinite_pixels':int((~finite).sum()),
        'minimum':float(values.min()) if values.size else None,'maximum':float(values.max()) if values.size else None,
        'pixel_SHA256':pixels_digest(data),'gzip_CRC_verified':True,'label_values':json.dumps(labels),
        'usable_file':not reasons,'file_reasons':'|'.join(reasons)}
    return img,data,record


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'EXPORT_SUCCESS.json').exists(),'Successful preparation exists; do not repeat'
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline)==904
    listed={str(Path(item['server_path']).relative_to(SOURCE)) for item in baseline}
    actual={str(p.relative_to(SOURCE)) for p in SOURCE.rglob('*') if p.is_file()}
    assert listed==actual,'Unaccounted for or missing files'
    originals=[]
    for item in baseline:
        path=Path(item['server_path'])
        md5,sha=digest(path,'md5'),digest(path)
        assert path.stat().st_size==int(item['bytes']) and md5==item['md5'] and sha==item['sha256']
        originals.append({'relative_path':str(path.relative_to(SOURCE)),'source_path':str(path),'bytes':path.stat().st_size,'MD5':md5,'SHA256':sha,'source_url':item['source_url']})
    pd.DataFrame(originals).to_csv(OUT/'ORIGINAL_FILES.tsv',sep='\t',index=False)
    print('All 904 source files matched download byte sizes/MD5/SHA256.',flush=True)
    patients,pairs,files=[],[],[]
    cfgs=sorted(SOURCE.glob('*/patient*/Info.cfg'))
    assert len(cfgs)==150
    for info_path in cfgs:
        patient=info_path.parent.name
        split=info_path.parent.parent.name
        cfg={key.strip():value.strip() for key,value in (line.split(':',1) for line in info_path.read_text().splitlines() if line.strip())}
        assert set(cfg)=={'ED','ES','Group','Height','NbFrame','Weight'}
        ed,es,frames=int(cfg['ED']),int(cfg['ES']),int(cfg['NbFrame'])
        assert cfg['Group'] in {'DCM','HCM','MINF','NOR','RV'}
        assert split in {'training','testing'} and ed!=es
        height,weight=float(cfg['Height']),float(cfg['Weight'])
        assert np.isfinite(height) and np.isfinite(weight) and height>0 and weight>0
        images={}
        paths=sorted(info_path.parent.glob('*.nii.gz'))
        assert len(paths)==5
        for path in paths:
            kind='cine4d' if path.name.endswith('_4d.nii.gz') else 'segmentation' if path.name.endswith('_gt.nii.gz') else 'phase_image'
            img,data,record=read_image(path,kind)
            record.update(patient_id=patient,official_split=split,diagnostic_group=cfg['Group'])
            files.append(record)
            images[path.name]=(img,data,record)
        cine_path=info_path.parent/f'{patient}_4d.nii.gz'
        cine,cine_data,cine_qc=images[cine_path.name]
        frame_meta_ok=cine_data.ndim==4 and cine_data.shape[3]==frames and 1<=ed<=frames and 1<=es<=frames
        patient_good=[]
        for phase,frame in [('ED',ed),('ES',es)]:
            image_path=info_path.parent/f'{patient}_frame{frame:02d}.nii.gz'
            mask_path=info_path.parent/f'{patient}_frame{frame:02d}_gt.nii.gz'
            image,values,image_qc=images[image_path.name]
            mask,labels,mask_qc=images[mask_path.name]
            shape_ok=values.shape==labels.shape==cine_data.shape[:3]
            affine_image_mask=np.allclose(image.affine,mask.affine,rtol=1e-6,atol=1e-5)
            affine_image_cine=np.allclose(image.affine,cine.affine,rtol=1e-6,atol=1e-5)
            same_frame=bool(frame_meta_ok and shape_ok and np.array_equal(values,cine_data[:,:,:,frame-1]))
            reasons=[]
            if not image_qc['usable_file']: reasons.append('invalid_phase_image')
            if not mask_qc['usable_file']: reasons.append('invalid_segmentation')
            if not cine_qc['usable_file']: reasons.append('invalid_4d_sequence')
            if not frame_meta_ok: reasons.append('frame_metadata_conflict')
            if not shape_ok: reasons.append('shape_mismatch')
            if not affine_image_mask or not affine_image_cine: reasons.append('affine_mismatch')
            if not same_frame: reasons.append('phase_values_differ_from_4d_frame')
            row={'patient_id':patient,'official_split':split,'diagnostic_group':cfg['Group'],'phase':phase,'frame_number_one_based':frame,
                'frame_index_zero_based':frame-1,'image_path':str(image_path),'mask_path':str(mask_path),'cine4d_path':str(cine_path),
                'image_pixel_SHA256':image_qc['pixel_SHA256'],'mask_pixel_SHA256':mask_qc['pixel_SHA256'],
                'shape_matches':shape_ok,'image_mask_affine_matches':bool(affine_image_mask),'image_cine_affine_matches':bool(affine_image_cine),
                'phase_exactly_matches_4d_frame':same_frame,'background_voxels':int((labels==0).sum()),'RV_voxels':int((labels==1).sum()),
                'MYO_voxels':int((labels==2).sum()),'LV_voxels':int((labels==3).sum()),'usable_pair':not reasons,'quarantine_reason':'|'.join(reasons)}
            pairs.append(row)
            patient_good.append(not reasons)
        patients.append({'patient_id':patient,'official_split':split,'diagnostic_group':cfg['Group'],'ED_one_based':ed,'ES_one_based':es,
            'NbFrame':frames,'Height_source_value':cfg['Height'],'Weight_source_value':cfg['Weight'],'info_path':str(info_path),
            'cine4d_path':str(cine_path),'usable_both_phases':all(patient_good),'usable_phase_pairs':sum(patient_good)})
        print('Patient checked:',split,patient,cfg['Group'],'usable phases',sum(patient_good),flush=True)
    ptab,itab,ftab=pd.DataFrame(patients),pd.DataFrame(pairs),pd.DataFrame(files)
    assert ptab.patient_id.is_unique and len(ptab)==150 and len(itab)==300 and len(ftab)==750
    assert ptab.groupby('official_split').size().to_dict()=={'testing':50,'training':100}
    assert set(ptab.loc[ptab.official_split.eq('training'),'patient_id'])=={f'patient{i:03d}' for i in range(1,101)}
    assert set(ptab.loc[ptab.official_split.eq('testing'),'patient_id'])=={f'patient{i:03d}' for i in range(101,151)}
    ptab.to_csv(OUT/'PATIENT_QC.tsv',sep='\t',index=False)
    itab.to_csv(OUT/'IMAGE_MASK_PAIR_QC.tsv',sep='\t',index=False)
    ftab.to_csv(OUT/'NIFTI_FILE_QC.tsv',sep='\t',index=False)
    itab.loc[itab.usable_pair].to_csv(OUT/'USABLE_IMAGE_MASK_PAIRS.tsv',sep='\t',index=False)
    itab.loc[~itab.usable_pair].to_csv(OUT/'QUARANTINED_IMAGE_MASK_PAIRS.tsv',sep='\t',index=False)
    ptab.groupby(['official_split','diagnostic_group']).agg(patients=('patient_id','size'),usable_both_phases=('usable_both_phases','sum'),usable_phase_pairs=('usable_phase_pairs','sum')).reset_index().to_csv(OUT/'GROUP_SPLIT_SUMMARY.tsv',sep='\t',index=False)
    policy={'dataset':'ACDC2017','source_root':str(SOURCE),'source_files':904,'patients':150,'NIfTI_files':750,'phase_pairs':300,
        'read_only_originals':True,'primary_delivery':'exact original image/mask paths plus usable/quarantine indices',
        'phase_numbering':'Info.cfg and filenames one-based; 4D array index=frame_number-1',
        'segmentation_labels':{'0':'background','1':'right_ventricular_cavity','2':'myocardium','3':'left_ventricular_cavity'},
        'affine_match_tolerance':{'rtol':1e-6,'atol':1e-5},'phase_pixel_comparison':'exact equality to corresponding 4D array frame',
        'official_split':'training patients001-100; testing patients101-150 preserved, no cross-split patient',
        'diagnostic_groups':'DCM/HCM/MINF/NOR/RV retained as published; no invented binary single-cell labels',
        'normalization':False,'resampling':False,'new_segmentation':False,'intensity_changes':False,'model_training':False,
        'license':'CC BY-NC-SA 4.0 and source mandatory citation retained','nibabel_version':importlib.metadata.version('nibabel')}
    (OUT/'qc_policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    for name in ['LICENSE_TERMS.md','MANDATORY_CITATION.md']:
        shutil.copyfile(SOURCE/name,OUT/name)
    report=f"# ACDC2017：图像、标注与病例索引\n\n904个原件与下载大小/MD5/SHA256一致；完整解压读取750个NIfTI，150位患者（训练100/测试50）、300个ED/ES图像-标注配对。可用配对{int(itab.usable_pair.sum())}，隔离{int((~itab.usable_pair).sum())}；双相位均可用患者{int(ptab.usable_both_phases.sum())}。\n\n核查每个像素的有限性、空间几何、尺寸、分割标签及ED/ES对应4D原帧逐值相同；gzip全流CRC通过。记录0背景/1右室腔/2心肌/3左室腔；来源诊断分组和官方划分原样保留。HCM只是其中的MRI诊断分组，不把影像当作单细胞计数或补造RNA疾病标签。\n\n原像素、强度和标注保持不变，不做重采样、归一化、自动分割或训练，以真实可用索引交付。异常配对有独立隔离清单及原因。原许可CC BY-NC-SA 4.0和强制引用已归档。准备后须独立验收，再检查退出码、进程结束和源码一致。\n"
    (OUT/'REPORT.md').write_text(report)
    records=[]
    for path in sorted(OUT.iterdir()):
        if path.is_file(): records.append({'file':path.name,'server_path':str(path),'bytes':path.stat().st_size,'MD5':digest(path,'md5'),'SHA256':digest(path)})
    pd.DataFrame(records).to_csv(OUT/'OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    (OUT/'EXPORT_SUCCESS.json').write_text(json.dumps({'status':'prepared_not_independently_verified','completed_utc':datetime.now(timezone.utc).isoformat(),'usable_pairs':int(itab.usable_pair.sum())},indent=2)+'\n')
    print('MRI preparation complete; independent verification required.',flush=True)


if __name__=='__main__':
    main()
