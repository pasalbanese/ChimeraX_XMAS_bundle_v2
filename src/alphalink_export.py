"""AlphaLink2 Colab export without ChimeraX/Qt dependencies.

Positions are indices in chain.characters, never PDB residue numbers.
"""
from dataclasses import dataclass
import math
import string

CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')


@dataclass(frozen=True)
class AlphaLinkExport:
    crosslinks: str
    fasta: str
    mapping: str
    count: int


def linked_chains(pseudobonds):
    """Return distinct chain objects in stable model/chain order."""
    chains = {}
    for pb in pseudobonds:
        for atom in pb.atoms:
            chain = atom.residue.chain
            if chain is None:
                raise ValueError('A crosslinked residue has no polymer chain.')
            chains[id(chain)] = chain
    return sorted(chains.values(), key=lambda c: (
        c.structure.id_string, c.chain_id))


def build_alphalink_export(pseudobonds, chains, fdr):
    """Build a matched FASTA, whitespace CSV, and residue mapping table.

    chains is the explicit FASTA order. fdr is a user-supplied probability
    shared by every restraint; search-engine scores are not probabilities.
    Reversed/replicate/atom-level duplicates collapse to one residue pair.
    """
    try:
        fdr = float(fdr)
    except (TypeError, ValueError):
        raise ValueError('Enter FDR as a probability between 0 and 1.') from None
    if not math.isfinite(fdr) or not 0 <= fdr <= 1:
        raise ValueError('Enter FDR as a probability between 0 and 1.')
    pseudobonds, chains = list(pseudobonds), list(chains)
    if not pseudobonds:
        raise ValueError('No crosslinks to export.')
    if not chains or len(chains) > len(CHAIN_IDS):
        raise ValueError('AlphaLink2 supports 1 to 62 chains (A-Z, a-z, 0-9).')
    if len({id(c) for c in chains}) != len(chains):
        raise ValueError('Each chain must appear exactly once in FASTA order.')

    positions = {}
    fasta = []
    mapping = ['alphalink_chain\tsequence_position\tmodel_id\tchain_id\tresidue_number\tinsertion_code']
    for index, chain in enumerate(chains):
        label = CHAIN_IDS[index]
        sequence = chain.characters.upper()
        residues = list(chain.residues)
        if not sequence or set(sequence) - AMINO_ACIDS:
            raise ValueError('Chain #%s/%s must contain only the 20 standard amino acids.'
                             % (chain.structure.id_string, chain.chain_id))
        if len(residues) != len(sequence):
            raise ValueError('Chain sequence and residue list lengths differ; cannot map safely.')
        fasta.extend(['>%s model_%s chain_%s' % (
            label, chain.structure.id_string, chain.chain_id), sequence])
        chain_positions = {}
        for position, residue in enumerate(residues, 1):
            if residue is None:
                number, insertion = '', ''
            else:
                chain_positions[id(residue)] = position
                number = str(residue.number)
                insertion = str(getattr(residue, 'insertion_code', '') or '').strip()
            mapping.append('\t'.join([label, str(position), chain.structure.id_string,
                                       chain.chain_id, number, insertion]))
        positions[id(chain)] = (label, chain_positions)

    pairs = set()
    for pb in pseudobonds:
        endpoints = []
        for atom in pb.atoms:
            residue = atom.residue
            try:
                label, chain_positions = positions[id(residue.chain)]
                position = chain_positions[id(residue)]
            except KeyError:
                raise ValueError('A crosslinked residue is absent from the exported chains.') from None
            endpoints.append((label, position))
        if len(endpoints) != 2:
            raise ValueError('Each crosslink must have exactly two endpoints.')
        if endpoints[0] == endpoints[1]:
            raise ValueError('A crosslink connects a residue to itself; resolve the mapping first.')
        pairs.add(tuple(sorted(endpoints)))
    lines = ['%s %s %s %s %s' % (p1, c1, p2, c2, format(fdr, '.12g'))
             for (c1, p1), (c2, p2) in sorted(pairs)]
    return AlphaLinkExport('\n'.join(lines) + '\n', '\n'.join(fasta) + '\n',
                           '\n'.join(mapping) + '\n', len(lines))
