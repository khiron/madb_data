import dataclasses
import json
import time
from typing import Any, Dict, Optional
from tqdm import tqdm 
import multiprocessing as mp
import warnings
import time

import cogent3
import madb

@dataclasses.dataclass
class thesis_rec:
    unique_id : str
    unaligned_seqs : Dict[str,str]
    ensembl_alignment : Dict[str,str]
    ensemble_score : float
    ensembl_pd : float
    cogent3_alignment : Dict[str,str]
    cogent3_score : float
    cogent3_time : float
    cogent3_pd : float
    jaccard_distance : Dict[int, float]
    madb_alignment : Dict[int,Dict[str,str]]
    madb_time : Dict[int,float]
    madb_score : Dict[int, float]    
    madb_distance : Dict[int, float]
    madb_bubbles : Dict[int,int]
    madb_braids : Dict[int, int]
    madb_cycles : Dict[int, bool]
    madb_longest_braid_length : Dict[int, int]
    madb_ungapped_smith_waterman_length : int
    madb_ungapped_smith_waterman : Optional[str] = None 
    ungapped_smith_waterman_time : Optional[float] = None
    def to_rich_dict(self) -> Dict[str, Any]:
        result = {}

        for field in dataclasses.fields(self):
            name = field.name
            value = getattr(self, name)

            # If the value has a .to_dict() method, call it
            if hasattr(value, "to_dict"):
                try:
                    result[name] = value.to_dict()
                except Exception as e:
                    result[name] = f"<<Failed to convert {name} via to_dict(): {e}>>"
            # If the value is a dict of objects with .to_dict()
            elif isinstance(value, dict):
                try:
                    result[name] = {
                        k: v.to_dict() if hasattr(v, "to_dict") else v
                        for k, v in value.items()
                    }
                except Exception as e:
                    result[name] = f"<<Failed to convert dict field {name}: {e}>>"
            # Default: copy value directly
            else:
                result[name] = value

        try:
            json.dumps(result)
        except TypeError as e:
            raise TypeError(f"{self.__class__.__name__}.to_rich_dict() returned non-serializable value: {e}")

        return result

    @classmethod
    def from_rich_dict(cls, d: Dict[str, Any]) -> "thesis_rec":
        if isinstance(d, thesis_rec):
            return d 
        # Get all field names and defaults
        field_names = {f.name for f in dataclasses.fields(cls)}
        init_kwargs = {}

        for field in field_names:
            if field in d:
                init_kwargs[field] = d[field]
            else:
                # Use the default or default_factory if provided
                dataclass_field = next(f for f in dataclasses.fields(cls) if f.name == field)
                if dataclass_field.default is not dataclasses.MISSING:
                    init_kwargs[field] = dataclass_field.default
                elif dataclass_field.default_factory is not dataclasses.MISSING:  # type: ignore
                    init_kwargs[field] = dataclass_field.default_factory()  # type: ignore
                else:
                    init_kwargs[field] = None  # fallback default

        return cls(**init_kwargs)

    def align_cogent3(self) -> "thesis_rec":
        if '' in self.unaligned_seqs:
            raise ValueError("Unaligned sequences contain empty strings.")
        # Perform pairwise alignment using Cogent3
        unaligned = cogent3.make_unaligned_seqs(data = self.unaligned_seqs, moltype="dna")
        # Perform pairwise alignment
        kmer_sizes = range(10, 100, 5)
    
        for kmer_size in kmer_sizes:
            jaccard = cogent3.get_app('jaccard_dist', k=kmer_size)
            jcdists = jaccard(unaligned)
            self.jaccard_distance[kmer_size] = jcdists[unaligned.names[0], unaligned.names[1]]

        start_time = time.time()
        # Perform global pairwise alignment
        score_matrix = cogent3.align.align.make_dna_scoring_dict(match=5, transition=-2, transversion=-4)
        gap_penalty = 4
        gap_extend = 1 
        alignment = cogent3.align.global_pairwise(unaligned.seqs[0], unaligned.seqs[1], score_matrix, gap_penalty, gap_extend)
        self.cogent3_alignment = alignment.to_dict()
        self.cogent3_time = time.time() - start_time
        self.cogent3_pd = alignment.distance_matrix('pdist')[alignment.names[0], alignment.names[1]]
        return self
    
    def align_madb(self, timeout_sec: int = 60) -> "thesis_rec":
        kmer_sizes = range(10, 100, 5)

        def madb_worker(seqs, k, queue):
            try:
                import cogent3
                import madb
                import time
                start_time = time.time()
                graph = madb.make_graph(seqs, kmer_size=k, moltype=cogent3.DNA)
                aligned = graph.align()
                elapsed = time.time() - start_time
                braid_diff, braid_nucleotides, total_braids, total_bubbles, has_cycles = graph.distance(braid_bubble_counts=True)
                longest_braid = graph.longest_braid().length
                queue.put({
                    "aligned": aligned.to_dict(),
                    "time": elapsed,
                    "graph": graph.to_dict(),
                    "distance": braid_diff / braid_nucleotides,
                    "bubbles": total_bubbles,
                    "braids": total_braids,
                    "cycles": has_cycles,
                    "longest_braid": longest_braid
                })
            except Exception as e:
                queue.put(e)

        with tqdm(total=len(kmer_sizes), desc="Aligning madb", unit="step") as pbar:
            for k in kmer_sizes:
                pbar.set_description(f"Aligning madb (k-mer size: {k})")

                queue = mp.Queue()
                p = mp.Process(target=madb_worker, args=(self.unaligned_seqs, k, queue))
                p.start()
                p.join(timeout_sec)

                if p.is_alive():
                    p.terminate()
                    p.join()
                    warnings.warn(f"MADB alignment for {self.unique_id}, k={k} timed out after {timeout_sec}s")
                    self.madb_alignment[k] = None
                    self.madb_time[k] = None
                    self.madb_distance[k] = None
                    self.madb_bubbles[k] = None
                    self.madb_braids[k] = None
                    self.madb_cycles[k] = None
                    self.madb_longest_braid_length[k] = None
                else:
                    result = queue.get()
                    if isinstance(result, Exception):
                        warnings.warn(f"MADB alignment for {self.unique_id}, k={k} failed: {result}")
                        self.madb_alignment[k] = None
                        self.madb_time[k] = None
                        self.madb_distance[k] = None
                        self.madb_bubbles[k] = None
                        self.madb_braids[k] = None
                        self.madb_cycles[k] = None
                        self.madb_longest_braid_length[k] = None
                    else:
                        self.madb_alignment[k] = result["graph"]
                        self.madb_time[k] = result["time"]
                        self.madb_score = {}  # optional: cleared, computed later
                        self.madb_distance[k] = result["distance"]
                        self.madb_bubbles[k] = result["bubbles"]
                        self.madb_braids[k] = result["braids"]
                        self.madb_cycles[k] = result["cycles"]
                        self.madb_longest_braid_length[k] = result["longest_braid"]

                pbar.update(1)

        return self    

    def calc_sw(self, timeout_sec: int = 60) -> "thesis_rec":

        def sw_worker(seqs, queue):
            try:
                from cogent3 import make_unaligned_seqs, get_app
                start_time = time.time()
                unaligned = make_unaligned_seqs(data=seqs, moltype="dna")
                sw = get_app('smith_waterman', moltype="dna", insertion_penalty=10_000)
                local_alignment = sw(unaligned)
                elapsed_time = time.time() - start_time
                length = len(local_alignment.seqs[0])
                queue.put((length, elapsed_time,local_alignment.to_dict()))
            except Exception as e:
                queue.put(e)

        queue = mp.Queue()
        p = mp.Process(target=sw_worker, args=(self.unaligned_seqs, queue))
        p.start()
        p.join(timeout_sec)

        if p.is_alive():
            p.terminate()
            p.join()
            warnings.warn(f"SW alignment for {self.unique_id} timed out after {timeout_sec}s")
            self.madb_ungapped_smith_waterman_length = None
            self.madb_ungapped_smith_waterman = None
            self.ungapped_smith_waterman_time = None
        else:
            result = queue.get()
            if isinstance(result, Exception):
                warnings.warn(f"SW alignment failed for {self.unique_id}: {result}")
                self.madb_ungapped_smith_waterman_length = 0
                self.madb_ungapped_smith_waterman = None
                self.ungapped_smith_waterman_time = None
            else:
                length, elapsed, local_alignment = result
                self.madb_ungapped_smith_waterman_length = length
                self.madb_ungapped_smith_waterman = local_alignment
                self.ungapped_smith_waterman_time = elapsed

        return self


    def score_alignments(self) -> "thesis_rec":

        alignment = cogent3.make_aligned_seqs(self.ensembl_alignment, moltype="dna")
        self.ensemble_score = alignment.alignment_quality('sp_score')

        alignment = cogent3.make_aligned_seqs(self.cogent3_alignment, moltype="dna")
        self.cogent3_score = alignment.alignment_quality('sp_score')
    
        kmer_sizes = range(10, 100, 5)
        with tqdm(total=len(kmer_sizes), desc="Aligning madb", unit="step") as pbar:
            for kmer_size in kmer_sizes:
                # Perform pairwise alignment using MADB
                alignment = cogent3.make_aligned_seqs(self.madb_alignment[kmer_size], moltype="dna")
                self.madb_score[kmer_size] = alignment.alignment_quality('sp_score')
                pbar.update(1)
        return self

    def describe(self) -> str:
        """
        Describe the object.
        """
        return f"thesis_rec({self.unique_id})"
    
def thesis_rec_from_alignment(align: cogent3.app.typing.AlignedSeqsType)-> thesis_rec:
    ensembl_pd = align.take_seqs(align.names)
    ensembl_alignment = align.to_dict()
    unaligned = align.degap()
    if len(unaligned) != 2:
        raise ValueError("Alignment must contain exactly two sequences.")
    if len(unaligned.seqs[0]) == 0 or len(unaligned.seqs[1]) == 0:
        raise ValueError("Unaligned sequences must not be empty.") 
    distance_mat = align.distance_matrix('pdist')
    ensembl_pdist = distance_mat[align.names[0],align.names[1]]
    new_rec = thesis_rec(
        unique_id = cogent3.app.data_store.get_data_source(align),
        unaligned_seqs = unaligned.to_dict(),
        ensembl_alignment = ensembl_alignment,
        ensemble_score = 0.0,
        ensembl_pd = ensembl_pdist,
        cogent3_alignment = {},
        cogent3_score = 0.0,
        cogent3_time = 0.0,
        cogent3_pd = 0.0,
        jaccard_distance = {},
        madb_alignment = {},
        madb_time = {},
        madb_score = {},
        madb_distance = {},
        madb_bubbles = {},
        madb_braids = {},
        madb_cycles = {},
        madb_longest_braid_length= {},
        madb_ungapped_smith_waterman_length = {},
        madb_ungapped_smith_waterman = None,
        ungapped_smith_waterman_time = None
    )
    return new_rec

@cogent3.app.composable.define_app
def create_thesis_rec(align: cogent3.app.typing.AlignedSeqsType)-> cogent3.app.typing.SerialisableType:
    """
    Create a thesis_rec object from an alignment.
    """
    return thesis_rec_from_alignment(align=align)
