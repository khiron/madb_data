import pathlib
from time import sleep
import cogent3
import madb 
from typing import Dict, List
import dataclasses
import shutil
from thesis_rec import create_thesis_rec, thesis_rec, thesis_rec_from_alignment
from tqdm import tqdm 


print(f"generating thesis data")
source_alignment_path = pathlib.Path('~/source/ensembl/primates100')
thesis_data = pathlib.Path('~/source/ensembl/thesis_data')

primates = {'gorilla': 'gorilla_gorilla', 'chimp': 'pan_troglodytes', 'macaque': 'macaca_mulatta'}
max_workers = 12

@cogent3.app.composable.define_app
def rename(align: cogent3.app.typing.AlignedSeqsType)->cogent3.app.typing.AlignedSeqsType:
    sleep(1)
    return align.rename_seqs(lambda x: x.split(':')[0])

with tqdm(total=len(primates), desc="Generating primate pair alignments", unit="step") as pbar:
    for primate in primates:
        pair_name = 'human_'+primate
        pair_path_root = thesis_data / pair_name
        pbar.set_description(f"Generating {pair_name} pair alignments")

        in_dstore = cogent3.open_data_store(source_alignment_path, suffix='fa') 
        out_dstore = cogent3.open_data_store(pair_path_root/'ensembl_alignments', suffix='fa', mode='w') 

        loader = cogent3.get_app('load_aligned', moltype='dna')
        select_4_primates = cogent3.get_app('take_named_seqs','homo_sapiens',*primates.values())
        select_pair = cogent3.get_app('take_named_seqs','homo_sapiens',primates[primate])
        omit_gaps = cogent3.get_app('omit_gap_pos', moltype="dna")
        writer = cogent3.get_app('write_seqs', data_store = out_dstore)
        app = loader + rename() + select_4_primates + select_pair + omit_gaps + writer 
        app.apply_to(in_dstore, show_progress=True, parallel=True, par_kw=dict(max_workers=max_workers))
        out_dstore.describe
        pbar.update(1)
        sleep(10)

@cogent3.app.composable.define_app
def rec_from_alignment(aln: cogent3.app.typing.AlignedSeqsType)->cogent3.app.typing.SerialisableType:
    sleep(1)
    return thesis_rec_from_alignment(aln)

@cogent3.app.composable.define_app
def cogent3(rec: cogent3.app.typing.SerialisableType)->cogent3.app.typing.SerialisableType:
    rec = thesis_rec.from_rich_dict(rec)
    sleep(1)
    return rec.align_cogent3()

with tqdm(total=len(primates), desc="Generating cogent3 alignments", unit="step") as pbar:
    for primate in primates:
        pair_name = 'human_'+primate
        pair_path_root = thesis_data / pair_name
        pbar.set_description(f"Generating {pair_name} cogent3 alignments")

        in_dstore = cogent3.open_data_store(pair_path_root/'ensembl_alignments', suffix='fa') 
        out_dstore = cogent3.open_data_store(pair_path_root/'cogent3', suffix='json', mode='w') 

        loader = cogent3.get_app('load_aligned', moltype='dna')
        writer = cogent3.get_app('write_json', data_store=out_dstore)
        app = loader + rec_from_alignment() + cogent3() + writer 
        app.apply_to(in_dstore, show_progress=True, parallel=True, par_kw=dict(max_workers=max_workers))
        print(out_dstore.describe)
        pbar.update(1)
        sleep(10)

@cogent3.app.composable.define_app
def sw(rec: cogent3.app.typing.SerialisableType)->cogent3.app.typing.SerialisableType:
    rec = thesis_rec.from_rich_dict(rec)
    sleep(1)
    return rec.calc_sw()

with tqdm(total=len(primates), desc="Generating ungapped smith-waterman", unit="step") as pbar:
    for primate in primates:
        pair_name = 'human_'+primate
        pair_path_root = thesis_data / pair_name
        pbar.set_description(f"Generating {pair_name} cogent3 alignments")

        in_dstore = cogent3.open_data_store(pair_path_root/'cogent3', suffix='json') 
        out_dstore = cogent3.open_data_store(pair_path_root/'smith_waterman', suffix='json', mode='w') 

        loader = cogent3.get_app('load_json', moltype='dna')
        writer = cogent3.get_app('write_json', data_store = out_dstore)
        app = loader + sw() + writer 
        app.apply_to(in_dstore, show_progress=True, parallel=True, par_kw=dict(max_workers=max_workers))
        print(out_dstore.describe)
        pbar.update(1)
        sleep(10)