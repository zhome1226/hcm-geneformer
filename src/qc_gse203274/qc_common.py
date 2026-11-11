"""Frozen marker panels; all outputs retain the complete published feature axis."""
PANELS={
 'Cardiomyocyte':['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],
 'Endothelial':['PECAM1','VWF','KDR','EMCN','CLDN5'],
 'Immune':['PTPRC','LST1','TYROBP','CD3D','CD79A'],
 'Mural':['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN'],
 'Neural':['NRXN1','SNAP25','SYT1','SLC1A3','PLP1'],
 'Adipocyte':['ADIPOQ','FABP4','PLIN1','LPL','CIDEC'],
 'Mesothelial':['WT1','MSLN','UPK3B','KRT8','KRT18']}
PAIRS=[('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),
 ('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),('Endothelial','Immune'),('Immune','Mural')]
