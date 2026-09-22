#!/usr/bin/env python3
"""
generate_embeddings.py
----------------------

This script computes protein sequence embeddings for a small FASTA file
using three pretrained protein language models used by the FANTASIA pipeline:

  * **ProtT5**: `Rostlab/prot_t5_xl_uniref50`
  * **Ankh3‑Large**: `ElnaggarLab/ankh3-large`

For each sequence, the script loads each model, tokenizes the sequence,
runs a forward pass to obtain the per-residue hidden states, and averages
the last hidden layer across residues to produce a fixed‑length embedding.
The results are saved in an HDF5 file with one group per model.  Each
group contains two datasets:

  * ``ids``: a list of sequence identifiers from the FASTA file
  * ``embeddings``: a 2‑D array of shape ``(n_sequences, embedding_dim)``

**Note**: These models are large (hundreds of millions to billions of
parameters) and may require a GPU and significant memory.  The script
processes sequences individually to avoid padding artifacts, as done in
FANTASIA, but this means runtime increases linearly with the number of
sequences.  For just a handful of proteins it should be manageable on a
modern workstation.

Dependencies:

  * `transformers` ≥ 4.29
  * `torch`
  * `h5py`
  * `biopython`

Install these with pip if necessary:

```
pip install transformers torch h5py biopython
```

Example usage:

```
python generate_embeddings.py \
  --fasta my_sequences.fasta \
  --output embeddings.h5 \
  --models esm2 prot_t5 ankh3
```

This will write the embeddings to ``embeddings.h5``.  You can then
provide this file to the database‑free FANTASIA lookup script as
``embeddings_path`` in your YAML configuration.
"""

import argparse
import csv
import gzip
import json
import os
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    import h5py  # type: ignore
except ImportError:
    h5py = None  # type: ignore
from Bio import SeqIO  # type: ignore
import numpy as np
import torch  # type: ignore
from transformers import AutoModel, AutoTokenizer  # type: ignore
try:
    from tqdm import tqdm  # type: ignore
except ImportError:  # pragma: no cover
    class _TqdmShim:
        def __init__(self, iterable, **_kwargs):
            self._iterable = iterable

        def __iter__(self):
            return iter(self._iterable)

        def close(self):
            pass

    def tqdm(iterable, **kwargs):  # type: ignore
        return _TqdmShim(iterable, **kwargs)

import re

# Mapping from user‑friendly model names to the underlying HuggingFace model
# identifiers and lookup keys.  FANTASIA’s lookup table uses the
# HuggingFace model ID as the key for ESM‑2, ProtT5, and Ankh.  Each
# entry in this registry has two fields:
#
#   ``hf_id`` – the HuggingFace identifier used to download the model and
#               tokenizer (passed to ``AutoModel.from_pretrained`` and
#               ``AutoTokenizer.from_pretrained``).
#   ``key`` – the key used in the lookup table and when writing the
#             embeddings.  This should match the names provided in the
#             distributed lookup bundle so that the lookup script can
#             associate query embeddings with the corresponding reference
#             model.
MODEL_REGISTRY: Dict[str, Dict[str, object]] = {
    # ProtT5 model.  The lookup table refers to this model as ``Prot-T5``.
    "prot_t5": {
        "hf_id": "Rostlab/prot_t5_xl_uniref50",
        "key": "Prot-T5",
        "tokenizer_use_fast": False,
    },
    # Ankh3-Large model.  The lookup table key is ``Ankh3-Large``.
    "ankh3": {
        "hf_id": "ElnaggarLab/ankh3-large",
        "key": "Ankh3-Large",
        "tokenizer_use_fast": False,
    },
}


def parse_fasta(fasta_path: Path) -> Tuple[List[str], List[str]]:
    """Parse a FASTA file and return sequence IDs and sequences.

    Parameters
    ----------
    fasta_path : Path
        Path to the input FASTA file.

    Returns
    -------
    ids : List[str]
        A list of sequence identifiers.
    seqs : List[str]
        A list of amino acid sequences (upper‑case, no whitespace).
    """
    import re
    
    AA_RE = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYBXZJOUacdefghiklmnpqrstvwybxzjou]+$')
    
    records = []
    cur_header = None
    cur_seq_parts = []
    seq_started = False

    if fasta_path.suffix in {".gz", ".gzip"}:
        handle = gzip.open(fasta_path, "rt")
    else:
        handle = fasta_path.open()
    
    with handle:
        for raw in handle:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                # flush previous
                if cur_header is not None:
                    records.append((cur_header, ''.join(cur_seq_parts)))
                cur_header = line[1:].strip()
                cur_seq_parts = []
                seq_started = False
                continue

            s = line.strip()
            if not seq_started:
                # If line looks like a sequence (only AA letters), start sequence.
                if AA_RE.match(s):
                    seq_started = True
                    cur_seq_parts.append(s.upper())
                else:
                    # header continuation (contains spaces, '=', digits, parentheses, etc.)
                    if cur_header is None:
                        cur_header = ''
                    cur_header += ' ' + s
            else:
                # sequence continuation lines
                cur_seq_parts.append(s.replace(' ', '').upper())

        # flush last
        if cur_header is not None:
            records.append((cur_header, ''.join(cur_seq_parts)))

    # Convert to the format your main() expects
    ids = [record[0] for record in records]
    seqs = [record[1] for record in records]
    return ids, seqs


RESIDUE_SPACE_MODELS = {"prot_t5"}
SPECIAL_TOKEN_MODELS = {"prot_t5", "ankh3"}
ANKH_MODELS = {"ankh3"}
T5_AMINO_ACID_REMAP = str.maketrans({"U": "X", "Z": "X", "O": "X", "B": "X"})
DEFAULT_EMBED_BATCH_SIZE = 8
DEFAULT_MAX_BATCH_RESIDUES = 12000
DEFAULT_QUEUE_BATCH_SIZE = 100


@dataclass
class EmbeddingFailure:
    sequence_id: str
    model_name: str
    error_type: str
    error_category: str
    error_message: str
    sequence_length: int
    fasta_path: str


@dataclass
class EmbeddingResult:
    embeddings: Dict[str, np.ndarray]
    failures: List[EmbeddingFailure]
    failed_ids: Set[str]
    embedding_dim: int
    layer_embeddings: Dict[int, Dict[str, np.ndarray]] = field(default_factory=dict)


def normalize_sequence_for_embedding(
    sequence: str,
    *,
    max_sequence_length: int,
) -> str:
    """Apply upstream-compatible pre-embedding length filtering."""
    normalized = sequence.replace(" ", "").replace("\n", "").upper()
    if max_sequence_length > 0:
        return normalized[:max_sequence_length]
    return normalized


def preprocess_sequence(seq: str, model_name: str) -> str:
    """Format the raw FASTA sequence for the requested model."""
    stripped = seq.replace(" ", "").replace("\n", "")
    
    if model_name in RESIDUE_SPACE_MODELS:
        if stripped.isupper():
            # Uppercase sequences (amino acids)
            # Replace ambiguous amino acids with X
            remapped = re.sub(r"[UZOB]", "X", stripped)
            spaced = " ".join(list(remapped))
            return spaced
        else:
            # Lowercase sequences (3Di/fold tokens)
            spaced = " ".join(list(stripped))
            return f"<fold2AA> {spaced}"

    if model_name in ANKH_MODELS:
        return "[NLU]" + stripped

    return stripped.upper()


def _pool_hidden_states(
    hidden: torch.Tensor,
    attention_mask: torch.Tensor,
    *,
    model_name: str,
) -> torch.Tensor:
    """Pool token-level hidden states into one embedding per sequence."""
    if model_name in ANKH_MODELS:
        pooled: List[torch.Tensor] = []
        valid_lengths = attention_mask.sum(dim=1).tolist()
        for row_idx, valid_length in enumerate(valid_lengths):
            length = max(int(valid_length) - 2, 1)
            pooled.append(hidden[row_idx, 1 : 1 + length].mean(dim=0))
        return torch.stack(pooled, dim=0)

    weights = attention_mask.unsqueeze(-1).to(hidden.dtype)
    lengths = weights.sum(dim=1).clamp_min(1.0)
    return (hidden * weights).sum(dim=1) / lengths


def _collect_hidden_state_layers(
    outputs: Any,
    model: AutoModel,
) -> Tuple[torch.Tensor, Dict[int, torch.Tensor]]:
    """Return the last hidden state plus a normalized layer-index mapping."""
    if hasattr(outputs, "last_hidden_state"):
        hidden = outputs.last_hidden_state  # type: ignore[attr-defined]
    elif isinstance(outputs, tuple):
        hidden = outputs[0]
    else:
        raise RuntimeError("Unknown output structure from model")

    hidden_states = getattr(outputs, "hidden_states", None)
    if hidden_states is None:
        return hidden, {}

    hidden_state_list = list(hidden_states)

    layer_map: Dict[int, torch.Tensor] = {}
    for layer_idx, layer_hidden in enumerate(hidden_state_list):
        layer_map[layer_idx] = layer_hidden
    return hidden, layer_map


def resolve_layer_indices(
    requested_layers: Optional[Sequence[int]],
    *,
    export_all_layers: bool,
    available_layers: Dict[int, torch.Tensor],
) -> List[int]:
    """Resolve the final set of layer indices to export."""
    if export_all_layers:
        return sorted(available_layers.keys())
    if not requested_layers:
        return []

    missing = [layer_idx for layer_idx in requested_layers if layer_idx not in available_layers]
    if missing:
        preview = ", ".join(str(layer) for layer in missing[:10])
        raise ValueError(f"Requested layer(s) not available from model output: {preview}")
    return sorted(dict.fromkeys(int(layer_idx) for layer_idx in requested_layers))


def compute_embeddings_batch(
    sequences: Sequence[str],
    model: AutoModel,
    tokenizer: AutoTokenizer,
    device: torch.device,
    *,
    model_name: str,
    layer_indices: Optional[Sequence[int]] = None,
    export_all_layers: bool = False,
) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    """Compute mean pooled embeddings for a batch of sequences.

    This function tokenizes the input amino acid sequence, runs a forward
    pass through the transformer model, and returns one pooled vector per
    sequence. Sequences may be padded within the batch; pooling respects
    the attention mask so that padding does not affect the output.

    Parameters
    ----------
    sequences : Sequence[str]
        Prepared model inputs.
    model : AutoModel
        A pretrained transformer model.
    tokenizer : AutoTokenizer
        The corresponding tokenizer for the model.
    device : torch.device
        Device to run inference on (``torch.device('cuda')`` or ``'cpu'``).

    Returns
    -------
    Tuple[np.ndarray, Dict[int, np.ndarray]]
        The first item is the default pooled last-layer embedding matrix.
        The second item optionally contains pooled embeddings for selected
        hidden layers keyed by layer index.
    """
    add_special_tokens = model_name in SPECIAL_TOKEN_MODELS
    inputs = tokenizer(
        list(sequences),
        return_tensors="pt",
        add_special_tokens=add_special_tokens,
        truncation=False,
        padding=True,
    )
    inputs = {key: val.to(device) for key, val in inputs.items()}
    with torch.inference_mode():
        if getattr(model.config, "is_encoder_decoder", False) and hasattr(model, "encoder"):
            encoder = model.encoder  # type: ignore[attr-defined]
            outputs = encoder(**inputs, output_hidden_states=True)  # type: ignore[arg-type]
        else:
            outputs = model(**inputs, output_hidden_states=True)
        hidden, available_layers = _collect_hidden_state_layers(outputs, model)
    attention_mask = inputs.get("attention_mask")
    if attention_mask is None:
        attention_mask = torch.ones(hidden.shape[:2], device=hidden.device, dtype=torch.long)
    embeddings = _pool_hidden_states(hidden, attention_mask, model_name=model_name)

    pooled_layers: Dict[int, np.ndarray] = {}
    selected_layers = resolve_layer_indices(
        layer_indices,
        export_all_layers=export_all_layers,
        available_layers=available_layers,
    )
    for layer_idx in selected_layers:
        pooled_layers[layer_idx] = (
            _pool_hidden_states(
                available_layers[layer_idx],
                attention_mask,
                model_name=model_name,
            )
            .detach()
            .cpu()
            .numpy()
        )
    return embeddings.detach().cpu().numpy(), pooled_layers


def compute_embedding(
    sequence: str,
    model: AutoModel,
    tokenizer: AutoTokenizer,
    device: torch.device,
    *,
    model_name: str,
    layer_indices: Optional[Sequence[int]] = None,
    export_all_layers: bool = False,
) -> np.ndarray:
    """Backward-compatible single-sequence wrapper around batch inference."""
    embedding, _layer_embeddings = compute_embeddings_batch(
        [sequence],
        model,
        tokenizer,
        device,
        model_name=model_name,
        layer_indices=layer_indices,
        export_all_layers=export_all_layers,
    )
    return embedding[0]


def classify_embedding_error(exc: Exception) -> Optional[str]:
    """Return a short category string if the error is recoverable."""
    message = str(exc)
    lowered = message.lower()

    if isinstance(exc, torch.cuda.OutOfMemoryError):
        return "cuda_out_of_memory"
    if "out of memory" in lowered or "cuda error: out of memory" in lowered:
        return "cuda_out_of_memory"
    if "token indices sequence length" in lowered and "longer than the specified maximum" in lowered:
        return "sequence_too_long"
    if "is longer than the maximum length" in lowered or "sequence length" in lowered and "maximum" in lowered:
        return "sequence_too_long"
    if "input_ids" in lowered and "sequence length" in lowered and "exceeds" in lowered:
        return "sequence_too_long"
    if "too many positional tokens" in lowered or "index out of range in self" in lowered:
        return "sequence_too_long"

    return None


def infer_embedding_dim(model: AutoModel) -> int:
    """Best-effort extraction of the model's embedding dimensionality."""
    for attr in ("hidden_size", "d_model", "embed_dim", "projection_dim"):
        value = getattr(model.config, attr, None)
        if isinstance(value, int) and value > 0:
            return value
    # Fall back to parameter inspection
    for param in model.parameters():
        if param.ndim >= 2:
            return int(param.shape[-1])
    raise RuntimeError("Could not determine embedding dimensionality for model.")


def iter_embedding_batches(
    records: Sequence[Tuple[str, str, str]],
    *,
    max_batch_size: int,
    max_batch_residues: int,
) -> Iterable[List[Tuple[str, str, str]]]:
    """Yield batches of similar-length prepared sequences."""
    if max_batch_size <= 0:
        raise ValueError("max_batch_size must be positive.")
    if max_batch_residues <= 0:
        raise ValueError("max_batch_residues must be positive.")

    sorted_records = sorted(records, key=lambda item: len(item[2]))
    batch: List[Tuple[str, str, str]] = []
    batch_residue_budget = 0

    for record in sorted_records:
        prepared_len = max(len(record[2].split()), 1)
        projected_budget = max(batch_residue_budget, prepared_len) * (len(batch) + 1)
        if batch and (
            len(batch) >= max_batch_size or projected_budget > max_batch_residues
        ):
            yield batch
            batch = []
            batch_residue_budget = 0

        batch.append(record)
        batch_residue_budget = max(batch_residue_budget, prepared_len)

    if batch:
        yield batch


def iter_queue_packages(
    records: Sequence[Tuple[str, str]],
    *,
    queue_batch_size: int,
) -> Iterable[List[Tuple[str, str]]]:
    """Yield upstream-style queue packages before model forward batching."""
    if queue_batch_size <= 0:
        raise ValueError("queue_batch_size must be positive.")
    for idx in range(0, len(records), queue_batch_size):
        yield list(records[idx : idx + queue_batch_size])


def load_model_components(
    model_name: str, device: torch.device
) -> tuple[AutoTokenizer, AutoModel, int]:
    """Load tokenizer and model for ``model_name`` and move to ``device``."""
    model_info = MODEL_REGISTRY[model_name]
    hf_model_id = model_info["hf_id"]
    # BRAKER4: if BRAKER4_HF_MODEL_PATH is set to a local snapshot directory,
    # load directly from it instead of via HF Hub.  This bypasses all HF cache
    # machinery and works without internet or HF_HOME propagation into the container.
    local_path = os.environ.get("BRAKER4_HF_MODEL_PATH", "")
    if local_path:
        hf_model_id = local_path
    print(f"Loading model {model_name} ({hf_model_id})…")
    tokenizer_kwargs = dict(model_info.get("tokenizer_kwargs", {}))
    model_kwargs = dict(model_info.get("model_kwargs", {}))
    tokenizer_use_fast = bool(model_info.get("tokenizer_use_fast", False))
    tokenizer = AutoTokenizer.from_pretrained(
        hf_model_id,
        use_fast=tokenizer_use_fast,
        **tokenizer_kwargs,
    )
    model = AutoModel.from_pretrained(hf_model_id, **model_kwargs)
    model.to(device)
    model.eval()
    embedding_dim = infer_embedding_dim(model)
    return tokenizer, model, embedding_dim


def embed_sequences_for_model(
    records: Sequence[Tuple[str, str]],
    model_name: str,
    device: torch.device,
    *,
    skip_ids: Set[str],
    fasta_path: Path,
    embed_batch_size: int,
    max_batch_residues: int,
    layer_indices: Optional[Sequence[int]],
    export_all_layers: bool,
    model_components: Optional[Tuple[AutoTokenizer, AutoModel, int]] = None,
) -> EmbeddingResult:
    """Compute embeddings for all sequences using a given model.

    Parameters
    ----------
    records : Sequence[Tuple[str, str]]
        Sequence identifiers paired with their amino acid strings.
    model_name : str
        Key of the model in ``MODEL_REGISTRY``.
    device : torch.device
        Device to run the model on.
    skip_ids : Set[str]
        Sequence identifiers that should be skipped (failed previously).
    fasta_path : Path
        The FASTA file currently being processed (for reporting).
    model_components : Optional[Tuple[AutoTokenizer, AutoModel, int]]
        Preloaded tokenizer/model tuple.  When provided the caller retains ownership
        of the model and it will be reused across chunks.

    Returns
    -------
    EmbeddingResult
        Contains successful embeddings, any failures, and identifiers that failed for this model.
    """
    owns_model = False
    if model_components is None:
        tokenizer, model, embedding_dim = load_model_components(model_name, device)
        owns_model = True
    else:
        tokenizer, model, embedding_dim = model_components

    embeddings: Dict[str, np.ndarray] = {}
    layer_embeddings: Dict[int, Dict[str, np.ndarray]] = {}
    failures: List[EmbeddingFailure] = []
    newly_failed: Set[str] = set()

    pending_records = [
        (seq_id, raw_sequence, preprocess_sequence(raw_sequence, model_name))
        for seq_id, raw_sequence in records
        if seq_id not in skip_ids
    ]

    progress = tqdm(
        total=len(pending_records),
        desc=f"{model_name} embeddings",
        unit="seq",
        leave=False,
    )

    def handle_single_failure(seq_id: str, raw_sequence: str, exc: Exception) -> None:
        category = classify_embedding_error(exc)
        if category is None:
            raise exc
        error_type = exc.__class__.__name__
        message = str(exc)
        print(
            f"  Skipping sequence {seq_id} for model {model_name}: {category} ({message})"
        )
        failures.append(
            EmbeddingFailure(
                sequence_id=seq_id,
                model_name=model_name,
                error_type=error_type,
                error_category=category,
                error_message=message,
                sequence_length=len(raw_sequence),
                fasta_path=str(fasta_path),
            )
        )
        newly_failed.add(seq_id)
        if device.type == "cuda":
            torch.cuda.empty_cache()

    for batch in iter_embedding_batches(
        pending_records,
        max_batch_size=embed_batch_size,
        max_batch_residues=max_batch_residues,
    ):
        batch_ids = [seq_id for seq_id, _raw_sequence, _prepared in batch]
        batch_raw_sequences = [raw_sequence for _seq_id, raw_sequence, _prepared in batch]
        batch_prepared = [prepared for _seq_id, _raw_sequence, prepared in batch]
        try:
            batch_embeddings, batch_layer_embeddings = compute_embeddings_batch(
                batch_prepared,
                model,
                tokenizer,
                device,
                model_name=model_name,
                layer_indices=layer_indices,
                export_all_layers=export_all_layers,
            )
            for seq_id, embedding in zip(batch_ids, batch_embeddings, strict=False):
                embeddings[seq_id] = embedding
            for layer_idx, layer_matrix in batch_layer_embeddings.items():
                layer_store = layer_embeddings.setdefault(layer_idx, {})
                for seq_id, embedding in zip(batch_ids, layer_matrix, strict=False):
                    layer_store[seq_id] = embedding
            progress.update(len(batch))
            continue
        except Exception:  # pylint: disable=broad-except
            if device.type == "cuda":
                torch.cuda.empty_cache()
            for seq_id, raw_sequence, prepared in batch:
                try:
                    single_embedding, single_layer_embeddings = compute_embeddings_batch(
                        [prepared],
                        model,
                        tokenizer,
                        device,
                        model_name=model_name,
                        layer_indices=layer_indices,
                        export_all_layers=export_all_layers,
                    )
                    embeddings[seq_id] = single_embedding[0]
                    for layer_idx, layer_matrix in single_layer_embeddings.items():
                        layer_store = layer_embeddings.setdefault(layer_idx, {})
                        layer_store[seq_id] = layer_matrix[0]
                except Exception as exc:  # pylint: disable=broad-except
                    handle_single_failure(seq_id, raw_sequence, exc)
                finally:
                    progress.update(1)

    progress.close()

    result = EmbeddingResult(
        embeddings=embeddings,
        layer_embeddings=layer_embeddings,
        failures=failures,
        failed_ids=newly_failed,
        embedding_dim=embedding_dim,
    )
    if owns_model and device.type == "cuda":
        torch.cuda.empty_cache()
    return result


def write_failure_report(path: Path, failures: Sequence[EmbeddingFailure]) -> None:
    """Write a CSV report describing sequences that were skipped."""
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sequence_id",
        "model_name",
        "error_type",
        "error_category",
        "error_message",
        "sequence_length",
        "fasta_path",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for failure in failures:
            writer.writerow(asdict(failure))


def stack_embeddings_for_model(
    model_name: str,
    success_ids: Sequence[str],
    per_model_embeddings: Dict[str, Dict[str, np.ndarray]],
    per_model_dims: Dict[str, int],
) -> np.ndarray:
    embeddings_for_model = per_model_embeddings[model_name]
    if not success_ids:
        return np.empty((0, per_model_dims[model_name]), dtype=np.float32)
    missing = [seq_id for seq_id in success_ids if seq_id not in embeddings_for_model]
    if missing:
        preview = ", ".join(missing[:5])
        raise RuntimeError(
            f"Missing embeddings for model '{model_name}' and sequence(s): {preview}"
        )
    stacked = np.stack([embeddings_for_model[seq_id] for seq_id in success_ids]).astype(
        np.float32
    )
    return stacked


def stack_layer_embeddings_for_model(
    model_name: str,
    layer_idx: int,
    success_ids: Sequence[str],
    per_model_layer_embeddings: Dict[str, Dict[int, Dict[str, np.ndarray]]],
    per_model_dims: Dict[str, int],
) -> np.ndarray:
    model_layers = per_model_layer_embeddings.get(model_name, {})
    embeddings_for_layer = model_layers.get(layer_idx, {})
    if not success_ids:
        return np.empty((0, per_model_dims[model_name]), dtype=np.float32)
    missing = [seq_id for seq_id in success_ids if seq_id not in embeddings_for_layer]
    if missing:
        preview = ", ".join(missing[:5])
        raise RuntimeError(
            f"Missing layer embeddings for model '{model_name}', layer {layer_idx}, "
            f"and sequence(s): {preview}"
        )
    return np.stack([embeddings_for_layer[seq_id] for seq_id in success_ids]).astype(np.float32)


def save_embeddings_archive(
    output_path: Path,
    model_names: Sequence[str],
    per_model_embeddings: Dict[str, Dict[str, np.ndarray]],
    per_model_layer_embeddings: Dict[str, Dict[int, Dict[str, np.ndarray]]],
    per_model_dims: Dict[str, int],
    success_ids: Sequence[str],
    success_sequences: Sequence[str],
) -> None:
    suffix = output_path.suffix.lower()
    if suffix in {".h5", ".hdf5"}:
        if h5py is None:
            raise RuntimeError(
                "h5py is not installed.  Please install h5py or use a .npz output."
            )
        with h5py.File(output_path, "w") as h5f:
            dt = h5py.string_dtype(encoding="utf-8")
            for model_name in model_names:
                model_group = h5f.create_group(model_name)
                stacked = stack_embeddings_for_model(
                    model_name, success_ids, per_model_embeddings, per_model_dims
                )
                model_group.create_dataset("ids", data=np.array(success_ids, dtype=dt))
                model_group.create_dataset(
                    "embeddings", data=stacked, compression="gzip"
                )
                print(
                    f"Saved embeddings for {stacked.shape[0]} sequences with {model_name} to group '{model_name}'"
                )
                layer_map = per_model_layer_embeddings.get(model_name, {})
                if layer_map:
                    layers_group = model_group.create_group("layers")
                    for layer_idx in sorted(layer_map):
                        layer_stacked = stack_layer_embeddings_for_model(
                            model_name,
                            layer_idx,
                            success_ids,
                            per_model_layer_embeddings,
                            per_model_dims,
                        )
                        layer_group = layers_group.create_group(f"layer_{layer_idx}")
                        layer_group.create_dataset(
                            "embeddings", data=layer_stacked, compression="gzip"
                        )
                        layer_group.attrs["layer_index"] = int(layer_idx)
        print(f"All embeddings saved to {output_path}")
        return
    if suffix != ".npz":
        raise ValueError(f"Unsupported output file extension '{suffix}'. Use .h5 or .npz.")

    embeddings_dict: Dict[str, np.ndarray] = {}
    npz_dict: Dict[str, Any] = {}
    layer_embeddings_dict: Dict[str, Dict[str, np.ndarray]] = {}
    for model_name in model_names:
        model_info = MODEL_REGISTRY[model_name]
        model_key = model_info.get("key", model_name)
        stacked = stack_embeddings_for_model(
            model_name, success_ids, per_model_embeddings, per_model_dims
        )
        embeddings_dict[model_key] = stacked
        npz_dict[f"{model_key}_ids"] = np.array(success_ids, dtype=object)
        npz_dict[f"{model_key}_embeddings"] = stacked
        print(f"Computed embeddings for {stacked.shape[0]} sequences with {model_name}")
        layer_map = per_model_layer_embeddings.get(model_name, {})
        if layer_map:
            export_layers: Dict[str, np.ndarray] = {}
            for layer_idx in sorted(layer_map):
                layer_stacked = stack_layer_embeddings_for_model(
                    model_name,
                    layer_idx,
                    success_ids,
                    per_model_layer_embeddings,
                    per_model_dims,
                )
                export_layers[f"layer_{layer_idx}"] = layer_stacked
                npz_dict[f"{model_key}_layer_{layer_idx}_embeddings"] = layer_stacked
            layer_embeddings_dict[str(model_key)] = export_layers
    accessions_arr = np.array(success_ids, dtype=object)
    sequences_arr = np.array(success_sequences, dtype=object)
    npz_dict["accessions"] = accessions_arr
    npz_dict["sequences"] = sequences_arr
    npz_dict["embeddings"] = np.array(embeddings_dict, dtype=object)
    if layer_embeddings_dict:
        npz_dict["layer_embeddings"] = np.array(layer_embeddings_dict, dtype=object)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, **npz_dict)
    print(f"All embeddings saved to {output_path} (NumPy .npz format)")


def process_chunk(
    chunk_path: Path,
    output_path: Path,
    model_names: Sequence[str],
    device: torch.device,
    *,
    failure_report: Optional[Path],
    model_cache: Dict[str, Tuple[AutoTokenizer, AutoModel, int]],
    queue_batch_size: int,
    max_sequence_length: int,
    embed_batch_size: int,
    max_batch_residues: int,
    model_batch_sizes: Dict[str, int],
    layer_indices: Optional[Sequence[int]],
    export_all_layers: bool,
) -> Tuple[int, int]:
    ids, sequences = parse_fasta(chunk_path)
    if not ids:
        raise ValueError(f"No sequences found in {chunk_path}")
    records: List[Tuple[str, str]] = [
        (
            seq_id,
            normalize_sequence_for_embedding(
                seq,
                max_sequence_length=max_sequence_length,
            ),
        )
        for seq_id, seq in zip(ids, sequences)
    ]
    failed_ids: Set[str] = set()
    all_failures: List[EmbeddingFailure] = []
    per_model_embeddings: Dict[str, Dict[str, np.ndarray]] = {}
    per_model_layer_embeddings: Dict[str, Dict[int, Dict[str, np.ndarray]]] = {}
    per_model_dims: Dict[str, int] = {}

    for model_name in model_names:
        tokenizer, model, embedding_dim = model_cache[model_name]
        model_embeddings: Dict[str, np.ndarray] = {}
        model_layer_embeddings: Dict[int, Dict[str, np.ndarray]] = {}
        model_failures: List[EmbeddingFailure] = []
        model_failed_ids: Set[str] = set()
        effective_batch_size = model_batch_sizes.get(model_name, embed_batch_size)

        for queue_records in iter_queue_packages(
            records,
            queue_batch_size=queue_batch_size,
        ):
            result = embed_sequences_for_model(
                queue_records,
                model_name,
                device,
                skip_ids=failed_ids | model_failed_ids,
                fasta_path=chunk_path,
                embed_batch_size=effective_batch_size,
                max_batch_residues=max_batch_residues,
                layer_indices=layer_indices,
                export_all_layers=export_all_layers,
                model_components=(tokenizer, model, embedding_dim),
            )
            model_embeddings.update(result.embeddings)
            for layer_idx, layer_values in result.layer_embeddings.items():
                model_layer_embeddings.setdefault(layer_idx, {}).update(layer_values)
            if result.failures:
                model_failures.extend(result.failures)
            if result.failed_ids:
                model_failed_ids.update(result.failed_ids)

        per_model_embeddings[model_name] = model_embeddings
        per_model_layer_embeddings[model_name] = model_layer_embeddings
        per_model_dims[model_name] = embedding_dim
        if model_failures:
            all_failures.extend(model_failures)
        if model_failed_ids:
            failed_ids.update(model_failed_ids)

    success_ids: List[str] = [seq_id for seq_id in ids if seq_id not in failed_ids]
    success_sequences: List[str] = [
        seq for seq_id, seq in records if seq_id not in failed_ids
    ]

    save_embeddings_archive(
        output_path,
        model_names,
        per_model_embeddings,
        per_model_layer_embeddings,
        per_model_dims,
        success_ids,
        success_sequences,
    )

    if failure_report is not None:
        write_failure_report(failure_report, all_failures)
        if not all_failures:
            print(f"No sequences skipped. Wrote empty report to {failure_report}.")
        else:
            print(
                f"Wrote failure report for {len(all_failures)} sequence(s) to {failure_report}"
            )

    return len(records), len(failed_ids)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate embeddings for protein sequences using multiple models.")
    parser.add_argument(
        "--fasta",
        type=Path,
        help="Input FASTA file containing protein sequences. Required unless --chunks-file is provided.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file to store embeddings. Required unless --chunks-file is provided. Extension determines format (.h5 or .npz)",
    )
    parser.add_argument(
        "--chunks-file",
        type=Path,
        default=None,
        help=(
            "JSON file describing multiple chunks to process in a single run. Each entry "
            "must contain 'fasta' and 'output' fields and may include 'failure_report'."
        ),
    )
    parser.add_argument(
        "--models",
        nargs="+",
        choices=list(MODEL_REGISTRY.keys()),
        default=list(MODEL_REGISTRY.keys()),
        help="Subset of models to run (default: all)",
    )
    parser.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="PyTorch device to run inference on (e.g., 'cuda' or 'cpu')",
    )
    parser.add_argument(
        "--queue-batch-size",
        "--sequence-queue-package",
        type=int,
        default=DEFAULT_QUEUE_BATCH_SIZE,
        help=(
            "Outer sequence package size before per-model forward batching "
            "(upstream embedding.queue_batch_size; default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--max-sequence-length",
        "--length-filter",
        type=int,
        default=0,
        help=(
            "Truncate sequences before embedding. "
            "0 disables truncation (upstream embedding.max_sequence_length)."
        ),
    )
    parser.add_argument(
        "--embed-batch-size",
        type=int,
        default=DEFAULT_EMBED_BATCH_SIZE,
        help="Maximum number of sequences to embed at once (default: %(default)s).",
    )
    parser.add_argument(
        "--max-batch-residues",
        type=int,
        default=DEFAULT_MAX_BATCH_RESIDUES,
        help=(
            "Upper bound on padded residues per embedding batch. "
            "Larger values can be faster but use more memory (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--model-batch-sizes",
        nargs="*",
        default=(),
        help=(
            "Optional per-model forward batch overrides in the form "
            "'prot_t5=8 ankh3=4'."
        ),
    )
    parser.add_argument(
        "--layer-indices",
        nargs="*",
        type=int,
        default=None,
        help=(
            "Optional transformer layer indices to export in addition to the default "
            "last-layer embeddings. Intended for embedding-only workflows."
        ),
    )
    parser.add_argument(
        "--all-layers",
        action="store_true",
        help=(
            "Export pooled embeddings for all available transformer layers in addition "
            "to the default last-layer embeddings."
        ),
    )
    parser.add_argument(
        "--failure-report",
        type=Path,
        default=None,
        help=(
            "Optional CSV report listing sequences skipped due to excessive length or memory usage "
            "(not supported when --chunks-file is provided; specify per-chunk paths in the JSON instead)."
        ),
    )
    args = parser.parse_args()

    if args.chunks_file is not None:
        if args.fasta is not None or args.output is not None or args.failure_report is not None:
            parser.error("--chunks-file cannot be combined with --fasta/--output/--failure-report.")
    else:
        if args.fasta is None or args.output is None:
            parser.error("--fasta and --output are required when --chunks-file is not provided.")

    if args.chunks_file is not None:
        with args.chunks_file.open("r", encoding="utf-8") as handle:
            raw_specs = json.load(handle)
        if not isinstance(raw_specs, list) or not raw_specs:
            raise ValueError(
                f"Chunks file {args.chunks_file} must contain a non-empty list of chunk specifications."
            )
        chunk_specs: List[Dict[str, Optional[Path]]] = []
        for idx, entry in enumerate(raw_specs):
            if not isinstance(entry, dict):
                raise ValueError("Each chunk specification must be a JSON object.")
            try:
                fasta_path = Path(entry["fasta"])
                output_path = Path(entry["output"])
            except KeyError as exc:  # pragma: no cover - defensive
                raise ValueError(
                    f"Chunk specification at index {idx} is missing required fields 'fasta'/'output'."
                ) from exc
            failure_path = entry.get("failure_report")
            chunk_specs.append(
                {
                    "fasta": fasta_path,
                    "output": output_path,
                    "failure_report": Path(failure_path) if failure_path else None,
                }
            )
    else:
        chunk_specs = [
            {
                "fasta": args.fasta,
                "output": args.output,
                "failure_report": args.failure_report,
            }
        ]

    device = torch.device(args.device)
    print(f"Using device: {device}")

    model_batch_sizes: Dict[str, int] = {}
    for item in args.model_batch_sizes:
        if "=" not in item:
            raise ValueError(
                f"Invalid --model-batch-sizes entry {item!r}; expected model=value."
            )
        model_name, raw_value = item.split("=", 1)
        model_name = model_name.strip()
        if model_name not in MODEL_REGISTRY:
            raise ValueError(
                f"Unknown model in --model-batch-sizes: {model_name!r}."
            )
        batch_value = int(raw_value)
        if batch_value <= 0:
            raise ValueError("Per-model batch sizes must be positive integers.")
        model_batch_sizes[model_name] = batch_value
        if batch_value >= args.queue_batch_size:
            raise ValueError(
                f"Per-model batch size for {model_name!r} ({batch_value}) must be "
                f"strictly smaller than queue_batch_size ({args.queue_batch_size})."
            )

    if args.embed_batch_size >= args.queue_batch_size:
        raise ValueError(
            "embed_batch_size must be strictly smaller than queue_batch_size "
            "to preserve upstream batching semantics."
        )

    model_cache: Dict[str, Tuple[AutoTokenizer, AutoModel, int]] = {}
    for model_name in args.models:
        model_cache[model_name] = load_model_components(model_name, device)

    total_sequences = 0
    total_failures = 0
    for spec in chunk_specs:
        chunk_total, chunk_failed = process_chunk(
            spec["fasta"],
            spec["output"],
            args.models,
            device,
            failure_report=spec["failure_report"],
            model_cache=model_cache,
            queue_batch_size=args.queue_batch_size,
            max_sequence_length=args.max_sequence_length,
            embed_batch_size=args.embed_batch_size,
            max_batch_residues=args.max_batch_residues,
            model_batch_sizes=model_batch_sizes,
            layer_indices=args.layer_indices,
            export_all_layers=args.all_layers,
        )
        total_sequences += chunk_total
        total_failures += chunk_failed

    if total_failures:
        print(
            f"Skipped {total_failures} sequence(s) that exceeded model limits or memory constraints."
        )
    else:
        print(
            f"Successfully embedded {total_sequences} sequence(s) across {len(chunk_specs)} chunk(s)."
        )

    if device.type == "cuda":
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
