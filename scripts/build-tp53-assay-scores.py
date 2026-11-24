#!/usr/bin/env python3
"""Build the TP53 score index from the published, CC0-cross-checked data."""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from itertools import product
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree

from openpyxl import load_workbook


REFSEQ = ('NM_000546.5', 'NM_000546.6')
GIACOMELLI = tuple(f'urn:mavedb:00000068-{screen}-1' for screen in 'abc')
KOTLER = 'urn:mavedb:00000059-a-1'
KOTLER_OFFSET = 101  # UniProt P04637 residue 102 is fragment residue 1.
MAVEDB_API = 'https://api.mavedb.org/api/v1/score-sets'
GENETIC_CODE = dict(zip(map(''.join, product('TCAG', repeat=3)),
                        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
AMINO_ACIDS = dict(zip('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter'.split(),
                       'ARNDCQEGHILKMFPSTWYV*'))
PUBLISHED_ALTERNATES = {*AMINO_ACIDS.values(), 'Z'}
SOURCE_URLS = {
    'tp53-refseq.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_000546.5,NM_000546.6&rettype=gb&retmode=xml',
    'giacomelli-table3.xlsx': 'https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-018-0204-y/MediaObjects/41588_2018_204_MOESM6_ESM.xlsx',
    'kotler-table-s3.xlsx': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115072/suppl/GSE115072_H1299_RFS_sequence_variants.xlsx',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source(cache, name, url):
    path = cache / name
    if not path.exists():
        with urlopen(Request(url, headers={'User-Agent': 'VariantSearch TP53 assay index'}), timeout=60) as response:
            contents = response.read()
        path.write_bytes(contents)
    return path


def translate(cds):
    return ''.join(GENETIC_CODE[cds[i:i + 3]] for i in range(0, len(cds), 3))


def load_cds(path):
    records = {}
    for record in ElementTree.parse(path).findall('./GBSeq'):
        accession = record.findtext('GBSeq_accession-version')
        feature = next(feature for feature in record.findall('GBSeq_feature-table/GBFeature')
                       if feature.findtext('GBFeature_key') == 'CDS')
        start, end = map(int, feature.findtext('GBFeature_location').split('..'))
        records[accession] = record.findtext('GBSeq_sequence')[start - 1:end].upper()
    require(set(records) == set(REFSEQ), 'Expected the two versioned TP53 RefSeq records')
    cds = records[REFSEQ[1]]
    require(len(cds) == 1182 and cds == records[REFSEQ[0]], 'TP53 RefSeq coding sequences differ')
    return cds


def load_mave(cache, urn, pmid, offset=0):
    basename = urn.rsplit(':', 1)[1]
    data = json.loads(source(cache, f'{basename}.json', f'{MAVEDB_API}/{urn}').read_text())
    require(data['urn'] == urn and data['license']['shortName'] == 'CC0', f'Unexpected score set or license for {urn}')
    require(any(pub['identifier'] == pmid for pub in data['primaryPublicationIdentifiers']), f'Unexpected paper for {urn}')
    target, = data['targetGenes']
    sequence = target['targetSequence']
    require(target['mappedHgncName'] == 'TP53' and sequence['sequenceType'] == 'dna', f'Unexpected target for {urn}')
    if offset:
        require(any(item['offset'] == offset and item['identifier']['identifier'] == 'P04637'
                    for item in target['externalIdentifiers']), f'Unexpected fragment offset for {urn}')
    protein = translate(sequence['sequence'].upper())
    scores = defaultdict(list)
    csv_path = source(cache, f'{basename}.csv', f'{MAVEDB_API}/{urn}/scores')
    with csv_path.open(newline='') as handle:
        for row in csv.DictReader(handle):
            match = re.fullmatch(r'p\.([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2}|\*)', row['hgvs_pro'])
            if not match or row['score'] in ('', 'NA'):
                continue
            old, residue, new = match.groups()
            residue = int(residue)
            old = AMINO_ACIDS[old]
            new = AMINO_ACIDS.get(new, new)
            require(residue <= len(protein) and protein[residue - 1] == old, f'Protein reference mismatch in {row["accession"]}')
            score = float(row['score'])
            require(math.isfinite(score), f'Non-finite score in {row["accession"]}')
            scores[f'{old}{residue + offset}{new}'].append(score)
    return protein, scores


def workbook_rows(path, skip=0):
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        rows = workbook.active.values
        for _ in range(skip):
            next(rows)
        header = next(rows)
        yield from (dict(zip(header, row, strict=True)) for row in rows)
    finally:
        workbook.close()


def build_giacomelli(cache, reference):
    protein = translate(reference)
    score_sets = []
    for urn in GIACOMELLI:
        target, scores = load_mave(cache, urn, '30224644')
        # The assay backbone has Arg72. RefSeq has Pro72.
        require(protein[71] == 'P' and target == protein[:71] + 'R' + protein[72:], f'Unexpected sequence for {urn}')
        score_sets.append(scores)
    table = source(cache, 'giacomelli-table3.xlsx', SOURCE_URLS['giacomelli-table3.xlsx'])
    columns = ('A549_p53WT_Nutlin-3_Z-score', 'A549_p53NULL_Nutlin-3_Z-score', 'A549_p53NULL_Etoposide_Z-score')
    index = {}
    for row in workbook_rows(table, skip=1):
        residue, old, new = row['Position'], row['AA_wt'], row['AA_variant']
        if not isinstance(residue, int) or not 1 <= residue <= 393 or residue == 72:
            continue
        if new not in PUBLISHED_ALTERNATES or new == old:
            continue
        require(old == protein[residue - 1] and row['Allele'] == f'{old}{residue}{new}', f'Unexpected Giacomelli allele {row["Allele"]}')
        key = f'{old}{residue}{new.replace("Z", "*")}'
        values = [row[column] for column in columns]
        for value, scores in zip(values, score_sets, strict=True):
            matched = scores.get(key, [])
            require(isinstance(value, (float, int)) and len(matched) == 1 and math.isclose(value, matched[0], abs_tol=1e-12, rel_tol=0),
                    f'Giacomelli score for {key} differs from the MaveDB record')
        require(key not in index, f'Duplicate Giacomelli allele {key}')
        index[key] = values
    require(len(index) == 7840, 'Unexpected Giacomelli coverage')
    return index


def build_kotler(cache, reference):
    target, scores = load_mave(cache, KOTLER, '29979965', KOTLER_OFFSET)
    require(target == translate(reference)[KOTLER_OFFSET:292], 'Kotler fragment does not match human TP53 residues 102–292')
    table = source(cache, 'kotler-table-s3.xlsx', SOURCE_URLS['kotler-table-s3.xlsx'])
    rows = list(workbook_rows(table))
    libraries = {}
    for row in rows:
        if row['Mut_type'] == f'wtp53_{row["SubLib"]}_seq':
            sequence = row['Var_seq']
            offset = reference.find(sequence)
            require(offset >= 0 and offset == reference.rfind(sequence), f'Cannot place {row["SubLib"]} on NM_000546.6')
            libraries[row['SubLib']] = (offset, sequence)
    require(set(libraries) == {'DBDA', 'DBDB', 'DBDC', 'DBDD'}, 'Missing Kotler sublibrary reference row')

    index = {}
    secondary_fields = ('Sec_mut_position', 'Sec_seq_change', 'Sec_codon_num', 'Sec_AA_change')
    for row in rows:
        if row['Backbone'] != 'wt' or row['Mut_type'] != 'Sub' or any(row[field] is not None for field in secondary_fields):
            continue
        construct = row['Identifier']
        change = re.fullmatch('([ACGT])>([ACGT])', row['Seq_change'])
        require(change is not None, f'Unexpected Kotler substitution in construct {construct}')
        ref, alt = change.groups()
        position = row['Position'] - 1
        start = position // 3 * 3
        within_codon = position - start
        codon = reference[start:start + 3]
        old = GENETIC_CODE[codon]
        new = GENETIC_CODE[codon[:within_codon] + alt + codon[within_codon + 1:]]
        if old == new:
            continue
        residue = start // 3 + 1
        offset, sequence = libraries[row['SubLib']]
        within_library = position - offset
        require(ref == reference[position] and 0 <= within_library < len(sequence)
                and row['Var_seq'] == sequence[:within_library] + alt + sequence[within_library + 1:]
                and row['Codon_num'] == residue and row['AA_change'].replace('Z', '*') == f'{old}>{new}',
                f'Kotler construct {construct} does not match the advertised single-base change')
        protein = f'{old}{residue}{new}'
        score = row['RFS_H1299']
        require(isinstance(score, (int, float)) and sum(math.isclose(score, value, abs_tol=1e-9, rel_tol=0)
                                                       for value in scores[protein]) == 1,
                f'Kotler score for construct {construct} differs from the MaveDB record')
        hgvs = f'C.{position + 1}{ref}>{alt}'
        require(hgvs not in index, f'Multiple Kotler constructs for {hgvs}')
        index[hgvs] = {'protein': protein, 'id': construct, 'score': score}
    require(len(index) == 1312, 'Unexpected Kotler coverage')
    return index


def js_export(name, index):
    rows = ',\n'.join(f'    {json.dumps(key)}: {json.dumps(value, allow_nan=False)}' for key, value in sorted(index.items()))
    return f'export const {name} = {{\n{rows}\n}};\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path, help='directory for the published input files')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'api/_tp53-published-scores.js')
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    reference = load_cds(source(args.cache, 'tp53-refseq.xml', SOURCE_URLS['tp53-refseq.xml']))
    giacomelli = build_giacomelli(args.cache, reference)
    kotler = build_kotler(args.cache, reference)
    contents = ('// Generated by scripts/build-tp53-assay-scores.py. See README.md for sources.\n'
                '// Giacomelli score order: WT nutlin-3, null nutlin-3, null etoposide.\n'
                + js_export('GIACOMELLI_SCORES', giacomelli) + '\n' + js_export('KOTLER_CONSTRUCTS', kotler))
    pending = args.output.with_suffix('.tmp')
    pending.write_text(contents, encoding='utf-8')
    pending.replace(args.output)
    print(f'Indexed {len(giacomelli):,} Giacomelli protein changes and {len(kotler):,} Kotler constructs.')


if __name__ == '__main__':
    main()
