# Known limitations

## Annotation

- HMDB XML is parsed in memory with `xml2::read_xml`. Memory requirements increase with export size.
- Adduct expansion generates hypotheses without formula, isotope, or coelution constraints sufficient to establish compound identity.
- MS-DIAL evidence-score normalization and database candidate ranking use heuristics. Score thresholds require empirical validation; scores are not identification probabilities.
- Empirical groups based on correlation and retention time can include unrelated compounds.
- Accurate-mass matches alone do not establish metabolite identity.
- Local knowledge-base entries require source citations and independent curation. Example entries are not evidence of biological associations.

## Statistics and pathways

Statistics and pathway analysis are experimental components developed separately from annotation. Validation remains necessary for model rank, contrast estimability, sample alignment, false discovery rate control, and ambiguous feature-to-compound assignments.

## Validation status

Full R package checks and end-to-end evaluation against independent synthetic and public datasets remain pending. Comparative performance has not been established.
