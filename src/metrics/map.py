import copairs.map as copairs
import pandas as pd

from preprocessing.io import split_parquet


def _index(meta, plate_types, ignore_dmso=False):
    '''Select samples to be used in mAP computation'''
    index = meta['Metadata_PlateType'].isin(plate_types)
    index &= (meta['Metadata_PertType'] != 'poscon')
    valid_cmpd = meta.loc[index, 'Metadata_JCP2022'].value_counts()
    valid_cmpd = valid_cmpd[valid_cmpd > 1].index
    index &= meta['Metadata_JCP2022'].isin(valid_cmpd)
    # TODO: This compound has many more replicates than any other. ignoring it
    # for now. This filter should be done early on.
    index &= (meta['Metadata_JCP2022'] != 'JCP2022_033954')
    if ignore_dmso:
        index &= (meta['Metadata_JCP2022'] != 'DMSO')
    return index.values


def _group_negcons(meta: pd.DataFrame):
    '''
    Hack to avoid mAP computation for negcons. Assign a unique id for every
    negcon so that no pairs are found for such samples.
    '''
    negcon_ix = (meta['Metadata_JCP2022'] == 'DMSO')
    n_negcon = negcon_ix.sum()
    negcon_ids = [f'DMSO_{i}' for i in range(n_negcon)]
    pert_id = meta['Metadata_JCP2022'].astype("category").cat.add_categories(negcon_ids)
    pert_id[negcon_ix] = negcon_ids
    meta['Metadata_JCP2022'] = pert_id


def average_precision_negcon(parquet_path, ap_path, plate_types):
    meta, vals, _ = split_parquet(parquet_path)
    ix = _index(meta, plate_types)
    meta = meta[ix].copy()
    vals = vals[ix]
    _group_negcons(meta)
    result = copairs.average_precision(
        meta,
        vals,
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=['Metadata_Well'],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_PertType', 'Metadata_JCP2022'],
        batch_size=20000)
    result = result.query('Metadata_PertType!="negcon"')
    result.reset_index(drop=True).to_parquet(ap_path)


def average_precision_nonrep(parquet_path, ap_path, plate_types):
    meta, vals, _ = split_parquet(parquet_path)
    ix = _index(meta, plate_types, ignore_dmso=True)
    meta = meta[ix].copy()
    vals = vals[ix]
    result = copairs.average_precision(
        meta,
        vals,
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=[],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_JCP2022'],
        batch_size=20000,
    )
    result.reset_index(drop=True).to_parquet(ap_path)


def mean_average_precision(ap_path, map_path, threshold=0.05):
    ap_scores = pd.read_parquet(ap_path)

    map_scores = copairs.mean_average_precision(ap_scores,
                                                'Metadata_JCP2022',
                                                threshold=threshold,
                                                null_size=10000,
                                                seed=0)
    map_scores.to_parquet(map_path)
