# ClinVar and TP53 evidence fixtures

Retrieved November 24, 2026. The fixtures keep the fields exercised by the tests.

| Fixture | Source and scope |
| --- | --- |
| `clinvar-2774649.xml` | [NCBI VCV EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=2774649). VCV002774649.2, including the complete SCV004359989.2 submission. |
| `clinvar-142714.xml` | [NCBI VCV EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=142714). VCV000142714.41 with five of the 17 SCVs, including an expert panel, a dissenting germline assertion, a somatic clinical-impact assertion, and oncogenicity. |
| `clinvar-375941.xml` | [NCBI VCV EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=375941). VCV000375941.8 with all three SCVs and two somatic indication records. |
| `pubmed-tp53-functional.xml` | [PubMed EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12826609,30224644,29979965&retmode=xml). Bibliographic metadata and short abstract excerpts for all three records. |
| `tp53-mutationview-r21.csv` | [NCI TP53 MutationView r21](https://tp53.cancer.gov/static/data/MutationView_r21.csv). Four source rows at codon 237, including the two M237L alleles, and two R175H rows with Kato, Kakudo and NCI cellular functional results. |
| `myvariant-m237l.json` | [MyVariant.info](https://myvariant.info/v1/variant/chr17:g.7577572T%3EA). Selected annotation fields for the browser replay. |

The Chromium smoke test also replays the [NCBI Variation Services contextual](https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/NM_000546.6%3Ac.709A%3ET/contextuals) and [equivalent](https://api.ncbi.nlm.nih.gov/variation/v0/spdi/NM_000546.6%3A850%3AA%3AT/all_equivalent_contextual) SPDI responses in the two human reference assemblies. The test adds a literal markup string to the replayed rationale to exercise safe rendering. The live NCBI ClinVar VCV for Variation ID [127821](https://www.ncbi.nlm.nih.gov/clinvar/variation/127821/) was also checked with the parser.
