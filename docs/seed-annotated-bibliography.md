# Seed Annotated Bibliography: Influenza Distance Metrics and Vaccine-Relevant Antigenic Change

Generated: 2026-05-29

Purpose: seed approximately 100 verified background citations for a manuscript evaluating the working claim that phylogenetic distance is a poor proxy for antigenic change in influenza, especially for vaccine-relevant interpretation of H1N1/H3N2 heterologous human strain panels.

## Verification and Scope Notes

- This is a purposive seed bibliography, not a PRISMA systematic review.
- Metadata were checked against Crossref DOI records on 2026-05-29 after inspecting the current manuscript, existing bibliography, decision log, human-prerequisites document, and publishability review.
- Retraction screening used Crossref DOI metadata relation/update fields for all records. A PubMed ID-conversion/ESummary scan mapped 67 of 100 DOIs to PMIDs and found no PubMed retraction publication-type flags on 2026-05-29; unmapped records still require final PubMed/Retraction Watch audit before submission.
- The one preprint is labeled as a preprint and should not be treated as peer-reviewed evidence.
- Older sources are included only when they are seminal method/background sources or necessary historical anchors for HAI, antigenic cartography, phylogenetics, or matrix-comparison methods.

## Search Strategy Summary

Sources consulted included the existing `products/project-refs.bib`, PubMed/Entrez-style topic searches, Crossref DOI/title searches, and cited method/resource records flagged by the local manuscript review. Search concepts included influenza antigenic drift, antigenic cartography, HAI, antigenic distance, sequence-to-antigenicity prediction, p-epitope, Hamming distance, vaccine mismatch, H1N1, H3N2, immune imprinting/original antigenic sin, human antibody landscapes, maximum-likelihood and neighbor-joining trees, tree-comparison metrics, Mantel tests, permutation tests, phylogenetic reporting, and reproducible workflows.

## Literature Gaps to Track

- Few papers directly test phylogenetic distance as a vaccine-relevant proxy against human HAI/cartographic distances for the same H1N1/H3N2 heterologous panels.
- Neuraminidase antigenic data are much thinner than HA/HAI data for this manuscript's likely panel, so NA should remain a clearly labeled limitation or future-work area.
- Human immune-history effects are well supported, but translating them into a single distance metric for individual heterologous panels remains unsettled.
- Mantel-style matrix tests are more appropriate than treating all pairwise distances as independent, but the methods literature also cautions that Mantel tests have limitations; robustness checks should avoid overclaiming p-values.
- Vaccine strain selection uses antigenic, genetic, surveillance, manufacturing, and epidemiological evidence; this bibliography should not be used to imply that any one distance metric is sufficient for selection.

## Annotated Bibliography

### 1. Smith, D. J., Lapedes, A. S., de Jong, J. C., Bestebroer, T. M., Rimmelzwaan, G. F., Osterhaus, A. D. M. E., et al. (2004). Mapping the Antigenic and Genetic Evolution of Influenza Virus. *Science*, 305(5682), 371-376. https://doi.org/10.1126/science.1097211
- BibTeX key: `smith2004AntigenicGeneticEvolution`
- Identifier: DOI [10.1126/science.1097211](https://doi.org/10.1126/science.1097211)
- Publication type: seminal older source
- Topical category: Influenza antigenic drift and antigenic cartography
- Annotation: This foundational antigenic-cartography paper linked H3N2 antigenic clusters to genetic evolution and remains a baseline citation for influenza antigenic drift. It is useful for explaining why antigenic maps are treated as a distinct empirical phenotype rather than a simple transformation of phylogenetic distance.
- How this citation might be used: introduction; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 2. Russell, C. A., Jones, T. C., Barr, I. G., Cox, N. J., Garten, R. J., Gregory, V., et al. (2008). The Global Circulation of Seasonal Influenza A (H3N2) Viruses. *Science*, 320(5874), 340-346. https://doi.org/10.1126/science.1154137
- BibTeX key: `russell2008GlobalCirculationH3N2`
- Identifier: DOI [10.1126/science.1154137](https://doi.org/10.1126/science.1154137)
- Publication type: seminal older source
- Topical category: Influenza antigenic evolution and global circulation
- Annotation: This study connects antigenic drift with global movement and replacement patterns in seasonal H3N2 viruses. It helps frame why vaccine-relevant interpretation requires antigenic and epidemiological context in addition to tree topology.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 3. Rambaut, A., Pybus, O. G., Nelson, M. I., Viboud, C., Taubenberger, J. K., & Holmes, E. C. (2008). The genomic and epidemiological dynamics of human influenza A virus. *Nature*, 453(7195), 615-619. https://doi.org/10.1038/nature06945
- BibTeX key: `rambaut2008GenomicEpidemiologicalDynamics`
- Identifier: DOI [10.1038/nature06945](https://doi.org/10.1038/nature06945); PMID [18418375](https://pubmed.ncbi.nlm.nih.gov/18418375/)
- Publication type: seminal older source
- Topical category: Influenza molecular evolution and phylogenetics
- Annotation: This paper characterizes the genomic and epidemiological dynamics of human influenza A viruses using evolutionary analysis. It provides broader phylogenetic context for why H3N2 often has ladder-like structure but still requires phenotype-aware interpretation.
- How this citation might be used: background; methods context
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 18418375 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 4. Garten, R. J., Davis, C. T., Russell, C. A., Shu, B., Lindstrom, S., Balish, A., et al. (2009). Antigenic and Genetic Characteristics of Swine-Origin 2009 A(H1N1) Influenza Viruses Circulating in Humans. *Science*, 325(5937), 197-201. https://doi.org/10.1126/science.1176225
- BibTeX key: `garten2009SwineOriginH1N1`
- Identifier: DOI [10.1126/science.1176225](https://doi.org/10.1126/science.1176225); PMID [19465683](https://pubmed.ncbi.nlm.nih.gov/19465683/)
- Publication type: seminal older source
- Topical category: H1N1 evolution and antigenic characterization
- Annotation: This early characterization of the 2009 pandemic H1N1 lineage is relevant for explaining why H1N1 temporal distance can be misleading around lineage replacement events. It supports careful subtype-specific framing rather than assuming chronological distance tracks genetic or antigenic distance.
- How this citation might be used: background; subtype framing
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 19465683 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 5. Wolf, Y. I., Viboud, C., Holmes, E. C., Koonin, E. V., & Lipman, D. J. (2006). Long intervals of stasis punctuated by bursts of positive selection in the seasonal evolution of influenza A virus. *Biology Direct*, 1(1). https://doi.org/10.1186/1745-6150-1-34
- BibTeX key: `wolf2006StasisBurstsInfluenza`
- Identifier: DOI [10.1186/1745-6150-1-34](https://doi.org/10.1186/1745-6150-1-34); PMID [17067369](https://pubmed.ncbi.nlm.nih.gov/17067369/)
- Publication type: seminal older source
- Topical category: Influenza antigenic and molecular evolution
- Annotation: This paper describes punctuated patterns in seasonal influenza evolution rather than smooth continuous change. It is useful for motivating skepticism toward simple temporal or path-length proxies when antigenic change can be episodic.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 17067369 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 6. Sandbulte, M. R., Westgeest, K. B., Gao, J., Xu, X., Klimov, A. I., Russell, C. A., et al. (2011). Discordant antigenic drift of neuraminidase and hemagglutinin in H1N1 and H3N2 influenza viruses. *Proceedings of the National Academy of Sciences*, 108(51), 20748-20753. https://doi.org/10.1073/pnas.1113801108
- BibTeX key: `sandbulte2011DiscordantDriftNAHA`
- Identifier: DOI [10.1073/pnas.1113801108](https://doi.org/10.1073/pnas.1113801108); PMID [22143798](https://pubmed.ncbi.nlm.nih.gov/22143798/)
- Publication type: peer-reviewed article
- Topical category: HA and NA evolution
- Annotation: This paper shows that neuraminidase and hemagglutinin can drift antigenically in discordant ways. It directly supports a limitation/future-work point that HA-only phylogenetic distance may miss vaccine-relevant antigenic change contributed by NA.
- How this citation might be used: discussion; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 22143798 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 7. Koel, B. F., Burke, D. F., Bestebroer, T. M., van der Vliet, S., Zondag, G. C. M., Vervaet, G., et al. (2013). Substitutions Near the Receptor Binding Site Determine Major Antigenic Change During Influenza Virus Evolution. *Science*, 342(6161), 976-979. https://doi.org/10.1126/science.1244730
- BibTeX key: `koel2013ReceptorBindingAntigenicChange`
- Identifier: DOI [10.1126/science.1244730](https://doi.org/10.1126/science.1244730)
- Publication type: peer-reviewed article
- Topical category: HA antigenic drift mechanisms
- Annotation: Koel and colleagues identify substitutions near the HA receptor-binding site that drive major antigenic transitions in H3N2. The paper supports the idea that specific substitutions, not overall phylogenetic distance alone, can dominate antigenic interpretation.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 8. Sun, H., Yang, J., Zhang, T., Long, L. P., Jia, K., Yang, G., et al. (2013). Using Sequence Data To Infer the Antigenicity of Influenza Virus. *mBio*, 4(4). https://doi.org/10.1128/mbio.00230-13
- BibTeX key: `sun2013SequenceInferAntigenicity`
- Identifier: DOI [10.1128/mbio.00230-13](https://doi.org/10.1128/mbio.00230-13); PMID [23820391](https://pubmed.ncbi.nlm.nih.gov/23820391/)
- Publication type: peer-reviewed article
- Topical category: Sequence distance versus antigenic distance
- Annotation: This study evaluates sequence-based inference of influenza antigenicity and is directly relevant to the manuscript's comparison of sequence-derived and antigenic metrics. It can help distinguish amino-acid similarity, antigenic phenotype prediction, and phylogenetic relatedness as related but non-identical targets.
- How this citation might be used: introduction; methods background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 23820391 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 9. Bedford, T., Suchard, M. A., Lemey, P., Dudas, G., Gregory, V., Hay, A. J., et al. (2014). Integrating influenza antigenic dynamics with molecular evolution. *eLife*, 3. https://doi.org/10.7554/elife.01914
- BibTeX key: `bedford2014AntigenicDynamicsMolecular`
- Identifier: DOI [10.7554/elife.01914](https://doi.org/10.7554/elife.01914); PMID [24497547](https://pubmed.ncbi.nlm.nih.gov/24497547/)
- Publication type: peer-reviewed article
- Topical category: Antigenic cartography and molecular evolution
- Annotation: This paper integrates antigenic dynamics with viral phylogenies across seasonal influenza lineages. It is a central citation for explaining how antigenic and molecular evolution can be modeled together while still remaining empirically distinct.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 24497547 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 10. Bedford, T., Riley, S., Barr, I. G., Broor, S., Chadha, M., Cox, N. J., et al. (2015). Global circulation patterns of seasonal influenza viruses vary with antigenic drift. *Nature*, 523(7559), 217-220. https://doi.org/10.1038/nature14460
- BibTeX key: `bedford2015GlobalCirculationAntigenicDrift`
- Identifier: DOI [10.1038/nature14460](https://doi.org/10.1038/nature14460); PMID [26053121](https://pubmed.ncbi.nlm.nih.gov/26053121/)
- Publication type: peer-reviewed article
- Topical category: Antigenic drift and global circulation
- Annotation: This study shows that seasonal influenza lineages differ in global circulation patterns and that antigenic drift affects those patterns. It gives subtype-aware context for why H1N1 and H3N2 should not be collapsed into one undifferentiated analysis.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26053121 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 11. Łuksza, M., & Lässig, M. (2014). A predictive fitness model for influenza. *Nature*, 507(7490), 57-61. https://doi.org/10.1038/nature13087
- BibTeX key: `luksza2014PredictiveFitnessInfluenza`
- Identifier: DOI [10.1038/nature13087](https://doi.org/10.1038/nature13087)
- Publication type: peer-reviewed article
- Topical category: Forecasting antigenic evolution
- Annotation: This paper proposes a predictive fitness model for influenza evolution using genetic and antigenic considerations. It is useful background for explaining why sequence features can be informative while still falling short of directly observed antigenic measurements.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 12. Liu, M., Zhao, X., Hua, S., Du, X., Peng, Y., Li, X., et al. (2015). Antigenic Patterns and Evolution of the Human Influenza A (H1N1) Virus. *Scientific Reports*, 5(1). https://doi.org/10.1038/srep14171
- BibTeX key: `liu2015AntigenicPatternsH1N1`
- Identifier: DOI [10.1038/srep14171](https://doi.org/10.1038/srep14171); PMID [26412348](https://pubmed.ncbi.nlm.nih.gov/26412348/)
- Publication type: peer-reviewed article
- Topical category: H1N1 antigenic evolution
- Annotation: This paper focuses on antigenic patterns in human H1N1 evolution. It can help the manuscript justify discussing H1N1 separately from H3N2, especially when temporal and lineage structure do not behave the same way.
- How this citation might be used: background; subtype framing
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26412348 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 13. Tewawong, N., Prachayangprecha, S., Vichiwattana, P., Korkong, S., Klinfueng, S., Vongpunsawad, S., et al. (2015). Assessing Antigenic Drift of Seasonal Influenza A(H3N2) and A(H1N1)pdm09 Viruses. *PLOS ONE*, 10(10), e0139958. https://doi.org/10.1371/journal.pone.0139958
- BibTeX key: `tewawong2015AssessingAntigenicDrift`
- Identifier: DOI [10.1371/journal.pone.0139958](https://doi.org/10.1371/journal.pone.0139958); PMID [26440103](https://pubmed.ncbi.nlm.nih.gov/26440103/)
- Publication type: peer-reviewed article
- Topical category: H1N1/H3N2 antigenic drift
- Annotation: This study assesses antigenic drift in seasonal A(H3N2) and A(H1N1)pdm09 viruses. It can support subtype-specific discussion of antigenic change in the post-2009 era and provide contrast to sequence-only summaries.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26440103 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 14. Neher, R. A., Bedford, T., Daniels, R. S., Russell, C. A., & Shraiman, B. I. (2016). Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza viruses. *Proceedings of the National Academy of Sciences*, 113(12). https://doi.org/10.1073/pnas.1525578113
- BibTeX key: `neher2016AntigenicPhenotypePrediction`
- Identifier: DOI [10.1073/pnas.1525578113](https://doi.org/10.1073/pnas.1525578113); PMID [26951657](https://pubmed.ncbi.nlm.nih.gov/26951657/)
- Publication type: peer-reviewed article
- Topical category: Sequence-to-antigenic phenotype prediction
- Annotation: Neher and Bedford develop methods to predict and visualize antigenic phenotypes of seasonal influenza viruses from genetic data. The paper is valuable for distinguishing modern antigenic phenotype prediction from simpler metrics such as Hamming distance or tree distance.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26951657 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 15. Anderson, C. S., McCall, P. R., Stern, H. A., Yang, H., & Topham, D. J. (2018). Antigenic cartography of H1N1 influenza viruses using sequence-based antigenic distance calculation. *BMC Bioinformatics*, 19(1). https://doi.org/10.1186/s12859-018-2042-4
- BibTeX key: `anderson2018H1N1SequenceCartography`
- Identifier: DOI [10.1186/s12859-018-2042-4](https://doi.org/10.1186/s12859-018-2042-4); PMID [29433425](https://pubmed.ncbi.nlm.nih.gov/29433425/)
- Publication type: peer-reviewed article
- Topical category: Sequence distance versus antigenic cartography
- Annotation: This paper creates antigenic maps for H1N1 using sequence-based antigenic distance calculations. It is especially relevant because the manuscript compares H1N1 sequence-derived metrics with cartographic distances.
- How this citation might be used: introduction; methods
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29433425 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 16. Petrova, V. N., & Russell, C. A. (2018). The evolution of seasonal influenza viruses. *Nature Reviews Microbiology*, 16(1), 47-60. https://doi.org/10.1038/nrmicro.2017.118
- BibTeX key: `petrova2018SeasonalInfluenzaEvolution`
- Identifier: DOI [10.1038/nrmicro.2017.118](https://doi.org/10.1038/nrmicro.2017.118)
- Publication type: review
- Topical category: Influenza antigenic evolution
- Annotation: This review summarizes how seasonal influenza viruses evolve antigenically and genetically. It is a strong background citation for readers who need the broader drift and vaccine-update context before the manuscript narrows to distance metrics.
- How this citation might be used: introduction; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 17. Morris, D. H., Gostic, K. M., Pompei, S., Bedford, T., Łuksza, M., Neher, R. A., et al. (2018). Predictive Modeling of Influenza Shows the Promise of Applied Evolutionary Biology. *Trends in Microbiology*, 26(2), 102-118. https://doi.org/10.1016/j.tim.2017.09.004
- BibTeX key: `morris2018PredictiveModelingInfluenza`
- Identifier: DOI [10.1016/j.tim.2017.09.004](https://doi.org/10.1016/j.tim.2017.09.004); PMID [29097090](https://pubmed.ncbi.nlm.nih.gov/29097090/)
- Publication type: review
- Topical category: Forecasting influenza evolution
- Annotation: This review explains the promise and limits of predictive modeling for influenza evolution. It can help position the manuscript as a metric-evaluation study rather than as a complete vaccine-strain forecasting framework.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29097090 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 18. Yao, Y., Li, X., Liao, B., Huang, L., He, P., Wang, F., et al. (2017). Predicting influenza antigenicity from Hemagglutintin sequence data based on a joint random forest method. *Scientific Reports*, 7(1). https://doi.org/10.1038/s41598-017-01699-z
- BibTeX key: `yao2017RandomForestAntigenicity`
- Identifier: DOI [10.1038/s41598-017-01699-z](https://doi.org/10.1038/s41598-017-01699-z); PMID [28484283](https://pubmed.ncbi.nlm.nih.gov/28484283/)
- Publication type: peer-reviewed article
- Topical category: Sequence-to-antigenicity prediction
- Annotation: This paper uses a joint random forest approach to predict influenza antigenicity from hemagglutinin sequence data. It is useful as a contrast to simpler sequence-distance measures and shows the field has moved beyond raw pairwise similarity.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28484283 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 19. Doud, M. B., Hensley, S. E., & Bloom, J. D. (2017). Complete mapping of viral escape from neutralizing antibodies. *PLOS Pathogens*, 13(3), e1006271. https://doi.org/10.1371/journal.ppat.1006271
- BibTeX key: `doud2017CompleteEscapeMapping`
- Identifier: DOI [10.1371/journal.ppat.1006271](https://doi.org/10.1371/journal.ppat.1006271); PMID [28288189](https://pubmed.ncbi.nlm.nih.gov/28288189/)
- Publication type: peer-reviewed article
- Topical category: HA antigenic sites and escape mapping
- Annotation: This study maps viral escape from neutralizing antibodies at high resolution. It supports the manuscript's caution that antigenic change can depend on specific epitope-level substitutions rather than aggregate phylogenetic distance.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28288189 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 20. Doud, M. B., Lee, J. M., & Bloom, J. D. (2018). How single mutations affect viral escape from broad and narrow antibodies to H1 influenza hemagglutinin. *Nature Communications*, 9(1). https://doi.org/10.1038/s41467-018-03665-3
- BibTeX key: `doud2018SingleMutationsEscape`
- Identifier: DOI [10.1038/s41467-018-03665-3](https://doi.org/10.1038/s41467-018-03665-3); PMID [29643370](https://pubmed.ncbi.nlm.nih.gov/29643370/)
- Publication type: peer-reviewed article
- Topical category: HA mutation effects and antigenic escape
- Annotation: This paper examines how individual H1 hemagglutinin mutations affect escape from broad and narrow antibodies. It provides mechanistic background for why two viruses close on a tree may still differ in vaccine-relevant antibody recognition.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29643370 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 21. Lee, J. M., Huddleston, J., Doud, M. B., Hooper, K. A., Wu, N. C., Bedford, T., et al. (2018). Deep mutational scanning of hemagglutinin helps predict evolutionary fates of human H3N2 influenza variants. *Proceedings of the National Academy of Sciences*, 115(35). https://doi.org/10.1073/pnas.1806133115
- BibTeX key: `lee2018DMSH3N2EvolutionaryFates`
- Identifier: DOI [10.1073/pnas.1806133115](https://doi.org/10.1073/pnas.1806133115); PMID [30104379](https://pubmed.ncbi.nlm.nih.gov/30104379/)
- Publication type: peer-reviewed article
- Topical category: H3N2 HA evolution and prediction
- Annotation: This paper uses deep mutational scanning to help predict the evolutionary fates of H3N2 variants. It can support discussion that evolutionary success and antigenic phenotype are tied to specific mutational effects rather than generic tree distance.
- How this citation might be used: discussion; robustness justification
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 30104379 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 22. Lee, J. M., Eguia, R., Zost, S. J., Choudhary, S., Wilson, P. C., Bedford, T., et al. (2019). Mapping person-to-person variation in viral mutations that escape polyclonal serum targeting influenza hemagglutinin. *eLife*, 8. https://doi.org/10.7554/elife.49324
- BibTeX key: `lee2019PersonVariationEscape`
- Identifier: DOI [10.7554/elife.49324](https://doi.org/10.7554/elife.49324); PMID [31452511](https://pubmed.ncbi.nlm.nih.gov/31452511/)
- Publication type: peer-reviewed article
- Topical category: Human antibody variation and antigenic escape
- Annotation: This paper maps person-to-person variation in serum escape mutations against influenza HA. It is highly relevant to the manuscript's human heterologous panel context because it emphasizes that antigenic interpretation can vary across immune histories.
- How this citation might be used: discussion; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 31452511 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 23. Huddleston, J., Barnes, J. R., Rowe, T., Xu, X., Kondor, R., Wentworth, D. E., et al. (2020). Integrating genotypes and phenotypes improves long-term forecasts of seasonal influenza A/H3N2 evolution. *eLife*, 9. https://doi.org/10.7554/elife.60067
- BibTeX key: `huddleston2020GenotypePhenotypeForecasts`
- Identifier: DOI [10.7554/elife.60067](https://doi.org/10.7554/elife.60067); PMID [32876050](https://pubmed.ncbi.nlm.nih.gov/32876050/)
- Publication type: peer-reviewed article
- Topical category: Influenza forecasting and antigenic phenotypes
- Annotation: Huddleston and colleagues show that combining genotype and phenotype improves long-term H3N2 forecasts. It supports a balanced argument: sequence and phylogeny are useful, but phenotype-aware data improve vaccine-relevant inference.
- How this citation might be used: discussion; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32876050 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 24. Perofsky, A. C., Huddleston, J., Hansen, C. L., Barnes, J. R., Rowe, T., Xu, X., et al. (2024). Antigenic drift and subtype interference shape A(H3N2) epidemic dynamics in the United States. *eLife*, 13. https://doi.org/10.7554/elife.91849
- BibTeX key: `perofsky2024AntigenicDriftSubtypeInterference`
- Identifier: DOI [10.7554/elife.91849](https://doi.org/10.7554/elife.91849); PMID [39319780](https://pubmed.ncbi.nlm.nih.gov/39319780/)
- Publication type: peer-reviewed article
- Topical category: H3N2 antigenic drift and epidemiology
- Annotation: This recent eLife study links antigenic drift and subtype interference with H3N2 epidemic dynamics. It helps update the background beyond classic cartography papers and reinforces the vaccine relevance of antigenic rather than purely genetic summaries.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 39319780 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 25. Air, G. M. (2015). Influenza virus antigenicity and broadly neutralizing epitopes. *Current Opinion in Virology*, 11, 113-121. https://doi.org/10.1016/j.coviro.2015.03.006
- BibTeX key: `air2015InfluenzaAntigenicityEpitopes`
- Identifier: DOI [10.1016/j.coviro.2015.03.006](https://doi.org/10.1016/j.coviro.2015.03.006); PMID [25846699](https://pubmed.ncbi.nlm.nih.gov/25846699/)
- Publication type: review
- Topical category: Influenza antigenicity and epitopes
- Annotation: This review summarizes influenza antigenicity and broadly neutralizing epitopes. It can help explain why antigenic change is an epitope-level phenotype and why broad phylogenetic distance is an imperfect proxy.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 25846699 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 26. Meyer, A. G., & Wilke, C. O. (2015). Geometric Constraints Dominate the Antigenic Evolution of Influenza H3N2 Hemagglutinin. *PLOS Pathogens*, 11(5), e1004940. https://doi.org/10.1371/journal.ppat.1004940
- BibTeX key: `meyer2015GeometricConstraintsH3`
- Identifier: DOI [10.1371/journal.ppat.1004940](https://doi.org/10.1371/journal.ppat.1004940); PMID [26020774](https://pubmed.ncbi.nlm.nih.gov/26020774/)
- Publication type: peer-reviewed article
- Topical category: H3N2 HA antigenic evolution
- Annotation: This paper argues that geometric constraints shape antigenic evolution of H3N2 hemagglutinin. It is useful for discussion because it suggests antigenic evolution has structural constraints not captured by simple pairwise tree distance.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26020774 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 27. Popova, L., Smith, K., West, A. H., Wilson, P. C., James, J. A., Thompson, L. F., et al. (2012). Immunodominance of Antigenic Site B over Site A of Hemagglutinin of Recent H3N2 Influenza Viruses. *PLoS ONE*, 7(7), e41895. https://doi.org/10.1371/journal.pone.0041895
- BibTeX key: `popova2012ImmunodominanceSiteB`
- Identifier: DOI [10.1371/journal.pone.0041895](https://doi.org/10.1371/journal.pone.0041895); PMID [22848649](https://pubmed.ncbi.nlm.nih.gov/22848649/)
- Publication type: peer-reviewed article
- Topical category: HA antigenic sites and immunodominance
- Annotation: This study highlights immunodominance of antigenic site B in recent H3N2 viruses. It supports careful interpretation of p-epitope and site-focused metrics, especially when antigenic changes are concentrated in particular regions.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 22848649 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 28. Yang, H., Carney, P. J., Chang, J. C., Guo, Z., Villanueva, J. M., & Stevens, J. (2015). Structure and receptor binding preferences of recombinant human A(H3N2) virus hemagglutinins. *Virology*, 477, 18-31. https://doi.org/10.1016/j.virol.2014.12.024
- BibTeX key: `yang2015ReceptorBindingH3N2HA`
- Identifier: DOI [10.1016/j.virol.2014.12.024](https://doi.org/10.1016/j.virol.2014.12.024); PMID [25617824](https://pubmed.ncbi.nlm.nih.gov/25617824/)
- Publication type: peer-reviewed article
- Topical category: H3N2 HA structure and receptor binding
- Annotation: This paper examines structure and receptor-binding preferences of recombinant human H3N2 hemagglutinins. It is lower priority but useful background when discussing how HA sequence changes can affect properties beyond simple antigenic-map position.
- How this citation might be used: background only
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 25617824 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 29. Lee, M. S., & Chen, J. S. E. (2004). Predicting Antigenic Variants of Influenza A/H3N2 Viruses. *Emerging Infectious Diseases*, 10(8), 1385-1390. https://doi.org/10.3201/eid1008.040107
- BibTeX key: `lee2004PredictingH3N2AntigenicVariants`
- Identifier: DOI [10.3201/eid1008.040107](https://doi.org/10.3201/eid1008.040107); PMID [15496238](https://pubmed.ncbi.nlm.nih.gov/15496238/)
- Publication type: seminal older source
- Topical category: Sequence predictors of antigenic variants
- Annotation: This earlier paper predicts antigenic variants of H3N2 viruses from sequence information. It is an older but relevant precursor to modern sequence-to-antigenic prediction work.
- How this citation might be used: background; methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 15496238 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 30. Liao, Y. C., Lee, M. S., Ko, C. Y., & Hsiung, C. A. (2008). Bioinformatics models for predicting antigenic variants of influenza A/H3N2 virus. *Bioinformatics*, 24(4), 505-512. https://doi.org/10.1093/bioinformatics/btm638
- BibTeX key: `liao2008BioinformaticsAntigenicVariants`
- Identifier: DOI [10.1093/bioinformatics/btm638](https://doi.org/10.1093/bioinformatics/btm638)
- Publication type: seminal older source
- Topical category: Sequence predictors of antigenic variants
- Annotation: This paper evaluates bioinformatics models for predicting antigenic variants of H3N2. It provides methodological context for why sequence-derived distances should be evaluated against antigenic data rather than assumed valid.
- How this citation might be used: background; methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 31. Steinbrück, L., Klingen, T. R., & McHardy, A. C. (2014). Computational Prediction of Vaccine Strains for Human Influenza A (H3N2) Viruses. *Journal of Virology*, 88(20), 12123-12132. https://doi.org/10.1128/jvi.01861-14
- BibTeX key: `steinbruck2014ComputationalVaccineStrains`
- Identifier: DOI [10.1128/jvi.01861-14](https://doi.org/10.1128/jvi.01861-14); PMID [25122778](https://pubmed.ncbi.nlm.nih.gov/25122778/)
- Publication type: peer-reviewed article
- Topical category: Vaccine strain prediction
- Annotation: This study presents computational prediction of vaccine strains for human H3N2 viruses. It is relevant because vaccine strain prediction depends on relationships among sequences, antigenic phenotype, and future circulation, not tree distance alone.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 25122778 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 32. Klingen, T. R., Reimering, S., Guzmán, C. A., & McHardy, A. C. (2018). In Silico Vaccine Strain Prediction for Human Influenza Viruses. *Trends in Microbiology*, 26(2), 119-131. https://doi.org/10.1016/j.tim.2017.09.001
- BibTeX key: `klingen2018InSilicoVaccinePrediction`
- Identifier: DOI [10.1016/j.tim.2017.09.001](https://doi.org/10.1016/j.tim.2017.09.001)
- Publication type: review
- Topical category: Vaccine strain selection and prediction
- Annotation: This review covers in silico vaccine strain prediction for human influenza viruses. It helps place the manuscript in the broader literature on computational tools used to inform vaccine decisions.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 33. Cai, Z., Zhang, T., & Wan, X. F. (2012). Antigenic distance measurements for seasonal influenza vaccine selection. *Vaccine*, 30(2), 448-453. https://doi.org/10.1016/j.vaccine.2011.10.051
- BibTeX key: `cai2012AntigenicDistanceVaccineSelection`
- Identifier: DOI [10.1016/j.vaccine.2011.10.051](https://doi.org/10.1016/j.vaccine.2011.10.051); PMID [22063385](https://pubmed.ncbi.nlm.nih.gov/22063385/)
- Publication type: peer-reviewed article
- Topical category: Antigenic distance metrics and vaccine selection
- Annotation: This paper evaluates antigenic distance measurements for seasonal influenza vaccine selection. It is directly relevant to the manuscript's claim that metric choice matters for vaccine-relevant interpretation.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 22063385 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 34. Deem, M. W., & Pan, K. (2009). The epitope regions of H1-subtype influenza A, with application to vaccine efficacy. *Protein Engineering Design and Selection*, 22(9), 543-546. https://doi.org/10.1093/protein/gzp027
- BibTeX key: `deem2009H1EpitopeRegions`
- Identifier: DOI [10.1093/protein/gzp027](https://doi.org/10.1093/protein/gzp027); PMID [19578121](https://pubmed.ncbi.nlm.nih.gov/19578121/)
- Publication type: seminal older source
- Topical category: p-epitope and sequence distance metrics
- Annotation: This paper defines H1 epitope regions and applies them to vaccine efficacy. It supports the manuscript's use and interpretation of epitope-focused sequence metrics for H1N1.
- How this citation might be used: methods; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 19578121 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 35. Pan, K., Subieta, K. C., & Deem, M. W. (2011). A novel sequence-based antigenic distance measure for H1N1, with application to vaccine effectiveness and the selection of vaccine strains. *Protein Engineering Design and Selection*, 24(3), 291-299. https://doi.org/10.1093/protein/gzq105
- BibTeX key: `pan2011H1N1SequenceAntigenicDistance`
- Identifier: DOI [10.1093/protein/gzq105](https://doi.org/10.1093/protein/gzq105); PMID [21123189](https://pubmed.ncbi.nlm.nih.gov/21123189/)
- Publication type: peer-reviewed article
- Topical category: p-epitope and sequence antigenic distance
- Annotation: This paper proposes a sequence-based antigenic distance measure for H1N1 and links it to vaccine effectiveness and strain selection. It is a key citation when comparing p-epitope-like measures with cartographic and phylogenetic distances.
- How this citation might be used: methods; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 21123189 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 36. Bonomo, M. E., & Deem, M. W. (2018). Predicting Influenza H3N2 Vaccine Efficacy From Evolution of the Dominant Epitope. *Clinical Infectious Diseases*, 67(7), 1129-1131. https://doi.org/10.1093/cid/ciy323
- BibTeX key: `bonomo2018DominantEpitopeVE`
- Identifier: DOI [10.1093/cid/ciy323](https://doi.org/10.1093/cid/ciy323)
- Publication type: peer-reviewed article
- Topical category: p-epitope and vaccine effectiveness
- Annotation: This study predicts H3N2 vaccine efficacy from evolution of the dominant epitope. It is useful for connecting epitope-level distance metrics to vaccine-effectiveness interpretation.
- How this citation might be used: discussion; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 37. Li, X., & Deem, M. W. (2016). Influenza evolution and H3N2 vaccine effectiveness, with application to the 2014/2015 season. *Protein Engineering Design and Selection*, 29(8), 309-315. https://doi.org/10.1093/protein/gzw017
- BibTeX key: `li2016H3N2EvolutionVaccineEffectiveness`
- Identifier: DOI [10.1093/protein/gzw017](https://doi.org/10.1093/protein/gzw017); PMID [27313229](https://pubmed.ncbi.nlm.nih.gov/27313229/)
- Publication type: peer-reviewed article
- Topical category: p-epitope and vaccine mismatch
- Annotation: This paper applies influenza evolution and epitope-based reasoning to H3N2 vaccine effectiveness, including the 2014-2015 season. It supports discussion of mismatch and why vaccine interpretation cannot rely on chronological or phylogenetic distance alone.
- How this citation might be used: discussion; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 27313229 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 38. Flannery, B., Zimmerman, R. K., Gubareva, L. V., Garten, R. J., Chung, J. R., Nowalk, M. P., et al. (2016). Enhanced Genetic Characterization of Influenza A(H3N2) Viruses and Vaccine Effectiveness by Genetic Group, 2014–2015. *Journal of Infectious Diseases*, 214(7), 1010-1019. https://doi.org/10.1093/infdis/jiw181
- BibTeX key: `flannery2016GeneticGroupVE`
- Identifier: DOI [10.1093/infdis/jiw181](https://doi.org/10.1093/infdis/jiw181); PMID [27190176](https://pubmed.ncbi.nlm.nih.gov/27190176/)
- Publication type: peer-reviewed article
- Topical category: Vaccine mismatch and H3N2 genetic groups
- Annotation: This study links enhanced genetic characterization of H3N2 viruses with vaccine effectiveness by genetic group. It can support discussion of vaccine mismatch while keeping claims tied to specific seasons and genetic groups.
- How this citation might be used: discussion; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 27190176 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 39. Belongia, E. A., & McLean, H. Q. (2019). Influenza Vaccine Effectiveness: Defining the H3N2 Problem. *Clinical Infectious Diseases*, 69(10), 1817-1823. https://doi.org/10.1093/cid/ciz411
- BibTeX key: `belongia2019H3N2Problem`
- Identifier: DOI [10.1093/cid/ciz411](https://doi.org/10.1093/cid/ciz411)
- Publication type: review
- Topical category: Vaccine effectiveness and H3N2 mismatch
- Annotation: This review frames low and variable H3N2 vaccine effectiveness as a persistent problem. It is valuable for motivating why antigenic-distance metrics matter for vaccine interpretation.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 40. Skowronski, D. M., & De Serres, G. (2018). Role of Egg-adaptation Mutations in Low Influenza A(H3N2) Vaccine Effectiveness During the 2012–2013 Season. *Clinical Infectious Diseases*, 67(9), 1474-1476. https://doi.org/10.1093/cid/ciy350
- BibTeX key: `skowronski2018EggAdaptationH3N2VE`
- Identifier: DOI [10.1093/cid/ciy350](https://doi.org/10.1093/cid/ciy350); PMID [29688295](https://pubmed.ncbi.nlm.nih.gov/29688295/)
- Publication type: peer-reviewed article
- Topical category: Vaccine mismatch and egg adaptation
- Annotation: This paper links egg-adaptation mutations with low H3N2 vaccine effectiveness. It is useful for reminding readers that vaccine mismatch can arise from production-related antigenic changes as well as circulating-virus evolution.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29688295 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 41. de Rooij, A. J. H., Lempers, V. J. C., Park, Y., Vicic, N., Han, A. X., Russell, C. A., et al. (2025). Reproducible and later vaccine strain selection can improve vaccine match to A/H3N2 seasonal influenza viruses. *npj Vaccines*, 10(1). https://doi.org/10.1038/s41541-025-01292-w
- BibTeX key: `deRooij2025ReproducibleVaccineSelection`
- Identifier: DOI [10.1038/s41541-025-01292-w](https://doi.org/10.1038/s41541-025-01292-w); PMID [41271798](https://pubmed.ncbi.nlm.nih.gov/41271798/)
- Publication type: peer-reviewed article
- Topical category: Vaccine strain selection and reproducibility
- Annotation: This recent study argues that reproducible and later vaccine strain selection can improve A/H3N2 vaccine match. It is useful for connecting this manuscript's distance-metric question to practical strain-selection workflows.
- How this citation might be used: discussion; robustness justification
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 41271798 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 42. Meng, J., Liu, J., Song, W., Li, H., Wang, J., Zhang, L., et al. (2024). PREDAC-CNN: predicting antigenic clusters of seasonal influenza A viruses with convolutional neural network. *Briefings in Bioinformatics*, 25(2). https://doi.org/10.1093/bib/bbae033
- BibTeX key: `meng2024PredacCnn`
- Identifier: DOI [10.1093/bib/bbae033](https://doi.org/10.1093/bib/bbae033); PMID [38343322](https://pubmed.ncbi.nlm.nih.gov/38343322/)
- Publication type: methods paper
- Topical category: Sequence-to-antigenic cluster prediction
- Annotation: This paper uses a convolutional neural network to predict antigenic clusters of seasonal influenza A viruses. It is lower priority but useful as a modern computational contrast to simple Hamming, p-epitope, and tree-distance metrics.
- How this citation might be used: background only
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 38343322 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 43. Chen, Y., Xu, Y., Cheng, Y., Qi, X., Bai, T., Yang, J., et al. (2026). A sequence-based proactive intelligence for influenza antigenic profiling improves vaccine strain selection. *openRxiv*. https://doi.org/10.64898/2026.04.18.719333
- BibTeX key: `chen2026SequenceProactiveIntelligencePreprint`
- Identifier: DOI [10.64898/2026.04.18.719333](https://doi.org/10.64898/2026.04.18.719333)
- Publication type: preprint
- Topical category: Sequence-based antigenic profiling
- Annotation: This current preprint proposes sequence-based proactive intelligence for influenza antigenic profiling and vaccine strain selection. It should be cited only with an explicit preprint caveat and rechecked before manuscript submission because it has not been peer reviewed.
- How this citation might be used: background only
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 44. Gupta, V., Earl, D. J., & Deem, M. W. (2006). Quantifying influenza vaccine efficacy and antigenic distance. *Vaccine*, 24(18), 3881-3888. https://doi.org/10.1016/j.vaccine.2006.01.010
- BibTeX key: `gupta2006VaccineEfficacyAntigenicDistance`
- Identifier: DOI [10.1016/j.vaccine.2006.01.010](https://doi.org/10.1016/j.vaccine.2006.01.010); PMID [16460844](https://pubmed.ncbi.nlm.nih.gov/16460844/)
- Publication type: seminal older source
- Topical category: p-epitope and antigenic distance
- Annotation: This paper introduces a p-epitope-style measure linking antigenic distance with influenza vaccine efficacy. It is a core citation for any manuscript comparing p-epitope, sequence distance, and phylogenetic distance.
- How this citation might be used: methods; introduction
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 16460844 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 45. Durviaux, S., Treanor, J., Beran, J., Duval, X., Esen, M., Feldman, G., et al. (2014). Genetic and Antigenic Typing of Seasonal Influenza Virus Breakthrough Cases from a 2008-2009 Vaccine Efficacy Trial. *Clinical and Vaccine Immunology*, 21(3), 271-279. https://doi.org/10.1128/CVI.00544-13
- BibTeX key: `durviaux2014GeneticAntigenicBreakthrough`
- Identifier: DOI [10.1128/CVI.00544-13](https://doi.org/10.1128/CVI.00544-13); PMID [24371255](https://pubmed.ncbi.nlm.nih.gov/24371255/)
- Publication type: peer-reviewed article
- Topical category: Sequence versus HAI antigenic classification
- Annotation: This vaccine trial analysis compares genetic and antigenic typing of breakthrough influenza cases and finds limits to sequence-only antigenic classification. It is directly aligned with the manuscript's warning against using sequence or tree-based proxies too casually.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 24371255 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 46. Hobson, D., Curry, R. L., Beare, A. S., & Ward-Gardner, A. (1972). The role of serum haemagglutination-inhibiting antibody in protection against challenge infection with influenza A2 and B viruses. *Epidemiology and Infection*, 70(4), 767-777. https://doi.org/10.1017/S0022172400022610
- BibTeX key: `hobson1972HaiProtection`
- Identifier: DOI [10.1017/S0022172400022610](https://doi.org/10.1017/S0022172400022610); PMID [4509641](https://pubmed.ncbi.nlm.nih.gov/4509641/)
- Publication type: seminal older source
- Topical category: HAI assay and protection
- Annotation: This classic challenge study links serum hemagglutination-inhibition antibody with protection against influenza. It is older but necessary for grounding why HAI titers remain central in vaccine and antigenic-distance discussions.
- How this citation might be used: introduction; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 4509641 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 47. Coudeville, L., Bailleux, F., Riche, B., Megas, F., Andre, P., & Ecochard, R. (2010). Relationship between haemagglutination-inhibiting antibody titres and clinical protection against influenza: development and application of a bayesian random-effects model. *BMC Medical Research Methodology*, 10(1). https://doi.org/10.1186/1471-2288-10-18
- BibTeX key: `coudeville2010HaiTitersProtection`
- Identifier: DOI [10.1186/1471-2288-10-18](https://doi.org/10.1186/1471-2288-10-18); PMID [20210985](https://pubmed.ncbi.nlm.nih.gov/20210985/)
- Publication type: methods paper
- Topical category: HAI titers and clinical protection
- Annotation: This paper models the relationship between HAI titers and clinical protection. It helps the manuscript explain both the value and the limitations of HAI-derived antigenic data as vaccine-relevant evidence.
- How this citation might be used: background; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 20210985 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 48. Fonville, J. M., Wilks, S. H., James, S. L., Fox, A., Ventresca, M., Aban, M., et al. (2014). Antibody landscapes after influenza virus infection or vaccination. *Science*, 346(6212), 996-1000. https://doi.org/10.1126/science.1256427
- BibTeX key: `fonville2014AntibodyLandscapes`
- Identifier: DOI [10.1126/science.1256427](https://doi.org/10.1126/science.1256427); PMID [25414313](https://pubmed.ncbi.nlm.nih.gov/25414313/)
- Publication type: peer-reviewed article
- Topical category: Human antibody landscapes and antigenic cartography
- Annotation: This paper introduces antibody landscapes after influenza infection or vaccination. It is central background for human heterologous strain panels because it shows immune responses across many strains are structured by exposure history.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 25414313 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 49. Fonville, J. M., Fraaij, P. L. A., de Mutsert, G., Wilks, S. H., van Beek, R., Fouchier, R. A. M., et al. (2016). Antigenic Maps of Influenza A(H3N2) Produced With Human Antisera Obtained After Primary Infection. *Journal of Infectious Diseases*, 213(1), 31-38. https://doi.org/10.1093/infdis/jiv367
- BibTeX key: `fonville2016HumanAntiseraMaps`
- Identifier: DOI [10.1093/infdis/jiv367](https://doi.org/10.1093/infdis/jiv367); PMID [26142433](https://pubmed.ncbi.nlm.nih.gov/26142433/)
- Publication type: peer-reviewed article
- Topical category: Human HAI antigenic maps
- Annotation: This study compares antigenic maps made from human primary-infection sera with ferret sera for H3N2. It is especially relevant because the manuscript uses human immune-response data rather than standard ferret antisera alone.
- How this citation might be used: introduction; methods background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26142433 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 50. Lessler, J., Riley, S., Read, J. M., Wang, S., Zhu, H., Smith, G. J. D., et al. (2012). Evidence for Antigenic Seniority in Influenza A (H3N2) Antibody Responses in Southern China. *PLoS Pathogens*, 8(7), e1002802. https://doi.org/10.1371/journal.ppat.1002802
- BibTeX key: `lessler2012AntigenicSeniority`
- Identifier: DOI [10.1371/journal.ppat.1002802](https://doi.org/10.1371/journal.ppat.1002802); PMID [22829765](https://pubmed.ncbi.nlm.nih.gov/22829765/)
- Publication type: peer-reviewed article
- Topical category: Immune imprinting and antibody profiles
- Annotation: This paper provides evidence for antigenic seniority in H3N2 antibody responses. It supports discussion that observed human HAI distances can reflect immune history, not only current viral sequence relationships.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 22829765 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 51. Andrews, S. F., Huang, Y., Kaur, K., Popova, L. I., Ho, I. Y., Pauli, N. T., et al. (2015). Immune history profoundly affects broadly protective B cell responses to influenza. *Science Translational Medicine*, 7(316). https://doi.org/10.1126/scitranslmed.aad0522
- BibTeX key: `andrews2015ImmuneHistoryBroadBCells`
- Identifier: DOI [10.1126/scitranslmed.aad0522](https://doi.org/10.1126/scitranslmed.aad0522); PMID [26631631](https://pubmed.ncbi.nlm.nih.gov/26631631/)
- Publication type: peer-reviewed article
- Topical category: Immune history and vaccine response
- Annotation: This study shows that immune history profoundly affects broadly protective B-cell responses to influenza. It is useful for framing human vaccine responses as history-dependent rather than as direct readouts of phylogenetic distance.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 26631631 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 52. Gostic, K. M., Ambrose, M., Worobey, M., & Lloyd-Smith, J. O. (2016). Potent protection against H5N1 and H7N9 influenza via childhood hemagglutinin imprinting. *Science*, 354(6313), 722-726. https://doi.org/10.1126/science.aag1322
- BibTeX key: `gostic2016ChildhoodHAImprinting`
- Identifier: DOI [10.1126/science.aag1322](https://doi.org/10.1126/science.aag1322); PMID [27846599](https://pubmed.ncbi.nlm.nih.gov/27846599/)
- Publication type: peer-reviewed article
- Topical category: Immune imprinting/original antigenic sin
- Annotation: This paper demonstrates birth-year and childhood hemagglutinin imprinting effects for avian influenza risk. It provides strong background for why first exposures can shape later heterologous responses.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 27846599 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 53. Gostic, K. M., Bridge, R., Brady, S., Viboud, C., Worobey, M., & Lloyd-Smith, J. O. (2019). Childhood immune imprinting to influenza A shapes birth year-specific risk during seasonal H1N1 and H3N2 epidemics. *PLOS Pathogens*, 15(12), e1008109. https://doi.org/10.1371/journal.ppat.1008109
- BibTeX key: `gostic2019BirthYearSeasonalRisk`
- Identifier: DOI [10.1371/journal.ppat.1008109](https://doi.org/10.1371/journal.ppat.1008109); PMID [31856206](https://pubmed.ncbi.nlm.nih.gov/31856206/)
- Publication type: peer-reviewed article
- Topical category: Immune imprinting and seasonal H1N1/H3N2
- Annotation: This study extends imprinting concepts to seasonal H1N1 and H3N2 epidemic risk. It is directly relevant to human strain-panel interpretation because immune history can vary by birth cohort and subtype.
- How this citation might be used: discussion; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 31856206 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 54. Cobey, S., & Hensley, S. E. (2017). Immune history and influenza virus susceptibility. *Current Opinion in Virology*, 22, 105-111. https://doi.org/10.1016/j.coviro.2016.12.004
- BibTeX key: `cobey2017ImmuneHistorySusceptibility`
- Identifier: DOI [10.1016/j.coviro.2016.12.004](https://doi.org/10.1016/j.coviro.2016.12.004); PMID [28088686](https://pubmed.ncbi.nlm.nih.gov/28088686/)
- Publication type: review
- Topical category: Immune history and influenza susceptibility
- Annotation: This review summarizes how immune history affects influenza susceptibility. It is a strong framing citation for limitations on interpreting human HAI panels as if all sera respond to viruses only by viral distance.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28088686 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 55. Monto, A. S., Malosh, R. E., Petrie, J. G., & Martin, E. T. (2017). The Doctrine of Original Antigenic Sin: Separating Good From Evil. *The Journal of Infectious Diseases*, 215(12), 1782-1788. https://doi.org/10.1093/infdis/jix173
- BibTeX key: `monto2017DoctrineOriginalAntigenicSin`
- Identifier: DOI [10.1093/infdis/jix173](https://doi.org/10.1093/infdis/jix173); PMID [28398521](https://pubmed.ncbi.nlm.nih.gov/28398521/)
- Publication type: review
- Topical category: Original antigenic sin
- Annotation: This review revisits original antigenic sin and separates harmful from useful aspects of the concept. It can help the manuscript use imprinting language precisely without overclaiming.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28398521 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 56. Henry, C., Palm, A. K. E., Krammer, F., & Wilson, P. C. (2018). From Original Antigenic Sin to the Universal Influenza Virus Vaccine. *Trends in Immunology*, 39(1), 70-79. https://doi.org/10.1016/j.it.2017.08.003
- BibTeX key: `henry2018OriginalSinUniversalVaccine`
- Identifier: DOI [10.1016/j.it.2017.08.003](https://doi.org/10.1016/j.it.2017.08.003); PMID [28867526](https://pubmed.ncbi.nlm.nih.gov/28867526/)
- Publication type: review
- Topical category: Original antigenic sin and universal vaccines
- Annotation: This review connects original antigenic sin to universal influenza vaccine design. It is useful for linking distance metrics and heterologous panels to vaccine-development concerns.
- How this citation might be used: introduction; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28867526 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 57. Arevalo, C. P., Le Sage, V., Bolton, M. J., Eilola, T., Jones, J. E., Kormuth, K. A., et al. (2020). Original antigenic sin priming of influenza virus hemagglutinin stalk antibodies. *Proceedings of the National Academy of Sciences*, 117(29), 17221-17227. https://doi.org/10.1073/pnas.1920321117
- BibTeX key: `arevalo2020OriginalSinStalk`
- Identifier: DOI [10.1073/pnas.1920321117](https://doi.org/10.1073/pnas.1920321117); PMID [32631992](https://pubmed.ncbi.nlm.nih.gov/32631992/)
- Publication type: peer-reviewed article
- Topical category: Immune imprinting and HA stalk antibodies
- Annotation: This paper studies original-antigenic-sin priming of HA stalk antibodies. It reinforces that antigenic responses can depend on immune memory and epitope focus, not only on overall viral relatedness.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32631992 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 58. Knight, M., Changrob, S., Li, L., & Wilson, P. C. (2020). Imprinting, immunodominance, and other impediments to generating broad influenza immunity. *Immunological Reviews*, 296(1), 191-204. https://doi.org/10.1111/imr.12900
- BibTeX key: `knight2020ImprintingImmunodominance`
- Identifier: DOI [10.1111/imr.12900](https://doi.org/10.1111/imr.12900)
- Publication type: review
- Topical category: Imprinting and immunodominance
- Annotation: This review summarizes imprinting, immunodominance, and other barriers to broad influenza immunity. It can support a cautious discussion of why human heterologous responses are difficult to reduce to a single antigenic-distance metric.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 59. Krammer, F. (2019). The human antibody response to influenza A virus infection and vaccination. *Nature Reviews Immunology*, 19(6), 383-397. https://doi.org/10.1038/s41577-019-0143-6
- BibTeX key: `krammer2019HumanAntibodyResponse`
- Identifier: DOI [10.1038/s41577-019-0143-6](https://doi.org/10.1038/s41577-019-0143-6)
- Publication type: review
- Topical category: Human influenza antibody response
- Annotation: This review provides broad background on the human antibody response to influenza infection and vaccination. It is a high-value citation for explaining the immunological substrate behind HAI panels and heterologous vaccine responses.
- How this citation might be used: introduction; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 60. Krammer, F., Smith, G. J. D., Fouchier, R. A. M., Peiris, M., Kedzierska, K., Doherty, P. C., et al. (2018). Influenza. *Nature Reviews Disease Primers*, 4(1). https://doi.org/10.1038/s41572-018-0002-y
- BibTeX key: `krammer2018InfluenzaPrimer`
- Identifier: DOI [10.1038/s41572-018-0002-y](https://doi.org/10.1038/s41572-018-0002-y); PMID [29955068](https://pubmed.ncbi.nlm.nih.gov/29955068/)
- Publication type: review
- Topical category: Influenza overview
- Annotation: This Nature Reviews Disease Primers article provides a broad influenza overview, including virology, immunity, and vaccines. It is useful for nontechnical readers when the rebuilt manuscript needs concise foundational context.
- How this citation might be used: introduction; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29955068 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 61. Dugan, H. L., Guthmiller, J. J., Arevalo, P., Huang, M., Chen, Y. Q., Neu, K. E., et al. (2020). Preexisting immunity shapes distinct antibody landscapes after influenza virus infection and vaccination in humans. *Science Translational Medicine*, 12(573). https://doi.org/10.1126/scitranslmed.abd3601
- BibTeX key: `dugan2020PreexistingImmunityLandscapes`
- Identifier: DOI [10.1126/scitranslmed.abd3601](https://doi.org/10.1126/scitranslmed.abd3601); PMID [33298562](https://pubmed.ncbi.nlm.nih.gov/33298562/)
- Publication type: peer-reviewed article
- Topical category: Human preexisting immunity and antibody landscapes
- Annotation: This study shows that preexisting immunity shapes distinct antibody landscapes after infection and vaccination. It supports the manuscript's emphasis that human antigenic data can encode immune history as well as viral antigenic relationships.
- How this citation might be used: discussion; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 33298562 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 62. Yang, B., Lessler, J., Zhu, H., Jiang, C. Q., Read, J. M., Hay, J. A., et al. (2020). Life course exposures continually shape antibody profiles and risk of seroconversion to influenza. *PLOS Pathogens*, 16(7), e1008635. https://doi.org/10.1371/journal.ppat.1008635
- BibTeX key: `yang2020LifeCourseExposures`
- Identifier: DOI [10.1371/journal.ppat.1008635](https://doi.org/10.1371/journal.ppat.1008635); PMID [32702069](https://pubmed.ncbi.nlm.nih.gov/32702069/)
- Publication type: peer-reviewed article
- Topical category: Life-course exposure and antibody profiles
- Annotation: This paper uses longitudinal antibody profiles to show how life-course exposures shape serological risk and response. It is a strong citation for discussing why human HAI-derived distances may not map cleanly onto HA phylogenies.
- How this citation might be used: discussion; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32702069 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 63. Li, Z. N., Liu, F., Gross, F. L., Kim, L., Ferdinands, J., Carney, P., et al. (2021). Antibody Landscape Analysis following Influenza Vaccination and Natural Infection in Humans with a High-Throughput Multiplex Influenza Antibody Detection Assay. *mBio*, 12(1). https://doi.org/10.1128/mBio.02808-20
- BibTeX key: `li2021AntibodyLandscapeMultiplex`
- Identifier: DOI [10.1128/mBio.02808-20](https://doi.org/10.1128/mBio.02808-20); PMID [33531397](https://pubmed.ncbi.nlm.nih.gov/33531397/)
- Publication type: peer-reviewed article
- Topical category: Human antibody landscapes and multiplex assays
- Annotation: This study uses a high-throughput multiplex assay to analyze antibody landscapes after vaccination and natural infection. It is useful background for future larger strain panels and alternatives to traditional HAI-only measurement.
- How this citation might be used: background; future work
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 33531397 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 64. Hinojosa, M., Shepard, S. S., Chung, J. R., King, J. P., McLean, H. Q., Flannery, B., et al. (2021). Impact of Immune Priming, Vaccination, and Infection on Influenza A(H3N2) Antibody Landscapes in Children. *The Journal of Infectious Diseases*, 224(3), 469-480. https://doi.org/10.1093/infdis/jiaa665
- BibTeX key: `hinojosa2021ChildrenH3AntibodyLandscapes`
- Identifier: DOI [10.1093/infdis/jiaa665](https://doi.org/10.1093/infdis/jiaa665); PMID [33090202](https://pubmed.ncbi.nlm.nih.gov/33090202/)
- Publication type: peer-reviewed article
- Topical category: H3N2 antibody landscapes and immune priming
- Annotation: This paper evaluates immune priming, vaccination, and infection effects on H3N2 antibody landscapes in children. It is helpful for discussing how age and priming can shape heterologous panels.
- How this citation might be used: discussion; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 33090202 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 65. Auladell, M., Phuong, H. V. M., Mai, L. T. Q., Tseng, Y. Y., Carolan, L., Wilks, S., et al. (2022). Influenza virus infection history shapes antibody responses to influenza vaccination. *Nature Medicine*, 28(2), 363-372. https://doi.org/10.1038/s41591-022-01690-w
- BibTeX key: `auladell2022InfectionHistoryVaccination`
- Identifier: DOI [10.1038/s41591-022-01690-w](https://doi.org/10.1038/s41591-022-01690-w)
- Publication type: peer-reviewed article
- Topical category: Infection history and vaccine response
- Annotation: This cohort study shows that influenza infection history shapes antibody responses to vaccination. It is directly relevant because the manuscript interprets human HAI reactivity across historical H3N2 strains.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 66. Tsang, T. K., Perera, R. A. P. M., Fang, V. J., Wong, J. Y., Shiu, E. Y., So, H. C., et al. (2022). Reconstructing antibody dynamics to estimate the risk of influenza virus infection. *Nature Communications*, 13(1). https://doi.org/10.1038/s41467-022-29310-8
- BibTeX key: `tsang2022ReconstructingAntibodyDynamics`
- Identifier: DOI [10.1038/s41467-022-29310-8](https://doi.org/10.1038/s41467-022-29310-8); PMID [35322048](https://pubmed.ncbi.nlm.nih.gov/35322048/)
- Publication type: peer-reviewed article
- Topical category: Serology and infection-history inference
- Annotation: This study reconstructs antibody dynamics to estimate influenza infection risk. It provides methodological context for moving from static pairwise distances toward longitudinal immune-history-aware models.
- How this citation might be used: discussion; future work
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 35322048 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 67. Kubo, M., & Miyauchi, K. (2020). Breadth of Antibody Responses during Influenza Virus Infection and Vaccination. *Trends in Immunology*, 41(5), 394-405. https://doi.org/10.1016/j.it.2020.03.005
- BibTeX key: `kubo2020BreadthAntibodyResponses`
- Identifier: DOI [10.1016/j.it.2020.03.005](https://doi.org/10.1016/j.it.2020.03.005)
- Publication type: review
- Topical category: Influenza antibody breadth
- Annotation: This review summarizes breadth of antibody responses during influenza infection and vaccination. It supports the manuscript's use of heterologous strain panels while emphasizing the complexity of interpreting breadth.
- How this citation might be used: background; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 68. Nuñez, I. A., Carlock, M. A., Allen, J. D., Owino, S. O., Moehling, K. K., Nowalk, P., et al. (2017). Impact of age and pre-existing influenza immune responses in humans receiving split inactivated influenza vaccine on the induction of the breadth of antibodies to influenza A strains. *PLOS ONE*, 12(11), e0185666. https://doi.org/10.1371/journal.pone.0185666
- BibTeX key: `nunez2017AgePreexistingBreadth`
- Identifier: DOI [10.1371/journal.pone.0185666](https://doi.org/10.1371/journal.pone.0185666); PMID [29091724](https://pubmed.ncbi.nlm.nih.gov/29091724/)
- Publication type: peer-reviewed article
- Topical category: Human vaccine response and breadth
- Annotation: This paper from the relevant research ecosystem studies age and pre-existing immune responses in split inactivated vaccine recipients. It is important context for the human data source and for interpreting breadth across influenza A strains.
- How this citation might be used: methods context; background
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29091724 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 69. Abreu, R. B., Clutter, E. F., Attari, S., Sautto, G. A., & Ross, T. M. (2020). IgA Responses Following Recurrent Influenza Virus Vaccination. *Frontiers in Immunology*, 11. https://doi.org/10.3389/fimmu.2020.00902
- BibTeX key: `abreu2020IgARecurrentVaccination`
- Identifier: DOI [10.3389/fimmu.2020.00902](https://doi.org/10.3389/fimmu.2020.00902); PMID [32508822](https://pubmed.ncbi.nlm.nih.gov/32508822/)
- Publication type: peer-reviewed article
- Topical category: Human vaccine response
- Annotation: This paper characterizes IgA responses following recurrent influenza vaccination. It broadens the immune-response background beyond HAI IgG endpoints and helps avoid treating HAI as the only immune correlate.
- How this citation might be used: methods context; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32508822 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 70. Allen, J. D., & Ross, T. M. (2021). Evaluation of Next-Generation H3 Influenza Vaccines in Ferrets Pre-Immune to Historical H3N2 Viruses. *Frontiers in Immunology*, 12. https://doi.org/10.3389/fimmu.2021.707339
- BibTeX key: `allen2021H3VaccinesPreImmuneFerrets`
- Identifier: DOI [10.3389/fimmu.2021.707339](https://doi.org/10.3389/fimmu.2021.707339); PMID [34475872](https://pubmed.ncbi.nlm.nih.gov/34475872/)
- Publication type: peer-reviewed article
- Topical category: H3N2 vaccine response and pre-immunity
- Annotation: This ferret study evaluates next-generation H3 vaccines in animals pre-immune to historical H3N2 viruses. It is lower priority for a human-panel manuscript but useful for the broader vaccine-design context.
- How this citation might be used: background only
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 34475872 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 71. Viboud, C., Gostic, K., Nelson, M. I., Price, G. E., Perofsky, A., Sun, K., et al. (2020). Beyond clinical trials: Evolutionary and epidemiological considerations for development of a universal influenza vaccine. *PLOS Pathogens*, 16(9), e1008583. https://doi.org/10.1371/journal.ppat.1008583
- BibTeX key: `viboud2020UniversalVaccineModeling`
- Identifier: DOI [10.1371/journal.ppat.1008583](https://doi.org/10.1371/journal.ppat.1008583); PMID [32970783](https://pubmed.ncbi.nlm.nih.gov/32970783/)
- Publication type: review
- Topical category: Universal vaccine, imprinting, and viral evolution
- Annotation: This review discusses evolutionary and epidemiological issues for universal influenza vaccines, including immune imprinting and viral evolution. It can help keep universal-vaccine framing realistic and not stronger than the evidence.
- How this citation might be used: introduction; discussion
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32970783 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 72. Erbelding, E. J., Post, D. J., Stemmy, E. J., Roberts, P. C., Augustine, A. D., Ferguson, S., et al. (2018). A Universal Influenza Vaccine: The Strategic Plan for the National Institute of Allergy and Infectious Diseases. *The Journal of Infectious Diseases*, 218(3), 347-354. https://doi.org/10.1093/infdis/jiy103
- BibTeX key: `erbelding2018UniversalVaccinePlan`
- Identifier: DOI [10.1093/infdis/jiy103](https://doi.org/10.1093/infdis/jiy103); PMID [29506129](https://pubmed.ncbi.nlm.nih.gov/29506129/)
- Publication type: guideline
- Topical category: Universal influenza vaccine strategy
- Annotation: This strategic plan outlines NIAID priorities for universal influenza vaccine development. It can justify why better distance metrics matter without implying this manuscript itself tests a universal vaccine.
- How this citation might be used: introduction; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29506129 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 73. Belongia, E. A., Simpson, M. D., King, J. P., Sundaram, M. E., Kelley, N. S., Osterholm, M. T., et al. (2016). Variable influenza vaccine effectiveness by subtype: a systematic review and meta-analysis of test-negative design studies. *The Lancet Infectious Diseases*, 16(8), 942-951. https://doi.org/10.1016/S1473-3099(16)00129-8
- BibTeX key: `belongia2016VEMetaAnalysis`
- Identifier: DOI [10.1016/S1473-3099(16)00129-8](https://doi.org/10.1016/S1473-3099(16)00129-8)
- Publication type: review
- Topical category: Vaccine effectiveness and subtype differences
- Annotation: This systematic review/meta-analysis summarizes variable vaccine effectiveness by influenza subtype. It supports the manuscript's subtype-specific framing and the practical importance of H3N2 mismatch.
- How this citation might be used: introduction; discussion
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 74. Kim, H., Webster, R. G., & Webby, R. J. (2018). Influenza Virus: Dealing with a Drifting and Shifting Pathogen. *Viral Immunology*, 31(2), 174-183. https://doi.org/10.1089/vim.2017.0141
- BibTeX key: `kim2018DriftingShiftingPathogen`
- Identifier: DOI [10.1089/vim.2017.0141](https://doi.org/10.1089/vim.2017.0141)
- Publication type: review
- Topical category: Influenza drift and vaccine challenge
- Annotation: This review summarizes influenza as a drifting and shifting pathogen. It is useful for a concise background paragraph aimed at vaccine scientists who do not need a full molecular-evolution review.
- How this citation might be used: introduction; background
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 75. Dang, C. C., Le, Q. S., Gascuel, O., & Le, V. S. (2010). FLU, an amino acid substitution model for influenza proteins. *BMC Evolutionary Biology*, 10(1). https://doi.org/10.1186/1471-2148-10-99
- BibTeX key: `dang2010FluSubstitutionModel`
- Identifier: DOI [10.1186/1471-2148-10-99](https://doi.org/10.1186/1471-2148-10-99); PMID [20384985](https://pubmed.ncbi.nlm.nih.gov/20384985/)
- Publication type: methods paper
- Topical category: Phylogenetic model choice
- Annotation: This paper introduces the FLU amino-acid substitution model for influenza proteins. It is a required methods citation if the manuscript uses the FLU model for HA phylogenetic inference.
- How this citation might be used: methods
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 20384985 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 76. Edgar, R. C. (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity. *BMC Bioinformatics*, 5(1). https://doi.org/10.1186/1471-2105-5-113
- BibTeX key: `edgar2004Muscle`
- Identifier: DOI [10.1186/1471-2105-5-113](https://doi.org/10.1186/1471-2105-5-113); PMID [15318951](https://pubmed.ncbi.nlm.nih.gov/15318951/)
- Publication type: methods paper
- Topical category: Multiple sequence alignment
- Annotation: This paper introduces MUSCLE for multiple sequence alignment. It supports the manuscript's alignment step and should be cited wherever the HA MSA is described.
- How this citation might be used: methods
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 15318951 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 77. Bodenhofer, U., Bonatesta, E., Horejš-Kainrath, C., & Hochreiter, S. (2015). msa: an R package for multiple sequence alignment. *Bioinformatics*, 31(24), 3997-3999. https://doi.org/10.1093/bioinformatics/btv494
- BibTeX key: `bodenhofer2015MsaRPackage`
- Identifier: DOI [10.1093/bioinformatics/btv494](https://doi.org/10.1093/bioinformatics/btv494)
- Publication type: methods paper
- Topical category: Multiple sequence alignment software
- Annotation: This paper documents the `msa` R package used to access multiple alignment algorithms. It is useful for software provenance and reproducibility in the methods section.
- How this citation might be used: methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 78. Schliep, K. P. (2011). phangorn: phylogenetic analysis in R. *Bioinformatics*, 27(4), 592-593. https://doi.org/10.1093/bioinformatics/btq706
- BibTeX key: `schliep2011Phangorn`
- Identifier: DOI [10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706); PMID [21169378](https://pubmed.ncbi.nlm.nih.gov/21169378/)
- Publication type: methods paper
- Topical category: Phylogenetic inference software
- Annotation: This paper documents `phangorn`, the R package currently used for phylogenetic reconstruction and tree comparison. It is a necessary software citation for the manuscript's ML, NJ, and tree-distance workflow.
- How this citation might be used: methods
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 21169378 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 79. Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3), 526-528. https://doi.org/10.1093/bioinformatics/bty633
- BibTeX key: `paradis2019Ape5`
- Identifier: DOI [10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633)
- Publication type: methods paper
- Topical category: Phylogenetic methods software
- Annotation: This paper describes the modern `ape` ecosystem for phylogenetics in R. It is useful if the rebuilt pipeline relies on `ape` data structures or tree utilities alongside `phangorn`.
- How this citation might be used: methods; reproducibility
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 80. Unknown author (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees.. *Molecular Biology and Evolution*. https://doi.org/10.1093/oxfordjournals.molbev.a040454
- BibTeX key: `saitou1987NeighborJoining`
- Identifier: DOI [10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- Publication type: seminal older source
- Topical category: Phylogenetic tree construction
- Annotation: This is the original neighbor-joining method paper. It should be cited when distance-based NJ trees are built from temporal, sequence, p-epitope, or cartographic matrices.
- How this citation might be used: methods
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 81. Felsenstein, J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368-376. https://doi.org/10.1007/BF01734359
- BibTeX key: `felsenstein1981MaximumLikelihoodTrees`
- Identifier: DOI [10.1007/BF01734359](https://doi.org/10.1007/BF01734359)
- Publication type: seminal older source
- Topical category: Maximum-likelihood phylogenetics
- Annotation: This foundational paper introduced maximum-likelihood approaches for phylogenetic inference from sequence data. It provides historical methods context for the ML tree used as a comparator.
- How this citation might be used: methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 82. Felsenstein, J. (1985). Confidence Limits on Phylogenies: An Approach Using the Bootstrap. *Evolution*, 39(4), 783. https://doi.org/10.2307/2408678
- BibTeX key: `felsenstein1985BootstrapPhylogenies`
- Identifier: DOI [10.2307/2408678](https://doi.org/10.2307/2408678)
- Publication type: seminal older source
- Topical category: Phylogenetic support and uncertainty
- Annotation: This paper introduced bootstrap confidence limits for phylogenies. It should support any branch-support or topology-stability checks added in the rebuilt analysis.
- How this citation might be used: methods; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 83. Shimodaira, H., & Hasegawa, M. (1999). Multiple Comparisons of Log-Likelihoods with Applications to Phylogenetic Inference. *Molecular Biology and Evolution*, 16(8), 1114-1116. https://doi.org/10.1093/oxfordjournals.molbev.a026201
- BibTeX key: `shimodaira1999ShTest`
- Identifier: DOI [10.1093/oxfordjournals.molbev.a026201](https://doi.org/10.1093/oxfordjournals.molbev.a026201)
- Publication type: methods paper
- Topical category: Tree comparison and likelihood tests
- Annotation: This paper introduces the Shimodaira-Hasegawa test for comparing phylogenetic trees by log-likelihood. It is required if the manuscript reports SH-test results or explains why likelihood comparisons are preferred.
- How this citation might be used: methods; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 84. Robinson, D., & Foulds, L. (1981). Comparison of phylogenetic trees. *Mathematical Biosciences*, 53(1-2), 131-147. https://doi.org/10.1016/0025-5564(81)90043-2
- BibTeX key: `robinson1981ComparisonPhylogeneticTrees`
- Identifier: DOI [10.1016/0025-5564(81)90043-2](https://doi.org/10.1016/0025-5564(81)90043-2)
- Publication type: seminal older source
- Topical category: Tree comparison metrics
- Annotation: This paper defines the Robinson-Foulds tree comparison distance. It supports reporting RF distances as topology summaries distinct from likelihood-based comparisons.
- How this citation might be used: methods; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 85. Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. *Molecular Biology and Evolution*, 32(1), 268-274. https://doi.org/10.1093/molbev/msu300
- BibTeX key: `nguyen2015Iqtree`
- Identifier: DOI [10.1093/molbev/msu300](https://doi.org/10.1093/molbev/msu300); PMID [25371430](https://pubmed.ncbi.nlm.nih.gov/25371430/)
- Publication type: methods paper
- Topical category: Maximum-likelihood phylogenetics
- Annotation: This paper introduces IQ-TREE for efficient maximum-likelihood phylogenetic inference. It is optional for the current R workflow but useful if future robustness checks compare tree-fitting engines.
- How this citation might be used: methods background; future workflow
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 25371430 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 86. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., et al. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*, 37(5), 1530-1534. https://doi.org/10.1093/molbev/msaa015
- BibTeX key: `minh2020Iqtree2`
- Identifier: DOI [10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015); PMID [32011700](https://pubmed.ncbi.nlm.nih.gov/32011700/)
- Publication type: methods paper
- Topical category: Maximum-likelihood phylogenetics
- Annotation: This paper describes IQ-TREE 2 and expanded model/inference capabilities. It can be used as a methods alternative citation if future pipeline versions move beyond the current R-only implementation.
- How this citation might be used: methods background; future workflow
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 32011700 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 87. Le, S. Q., & Gascuel, O. (2008). An Improved General Amino Acid Replacement Matrix. *Molecular Biology and Evolution*, 25(7), 1307-1320. https://doi.org/10.1093/molbev/msn067
- BibTeX key: `le2008LgMatrix`
- Identifier: DOI [10.1093/molbev/msn067](https://doi.org/10.1093/molbev/msn067)
- Publication type: methods paper
- Topical category: Amino-acid substitution models
- Annotation: This paper introduces the LG amino-acid replacement matrix. It is relevant to model-choice discussion because FLU should be justified against general protein models such as LG.
- How this citation might be used: methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 88. Posada, D., & Crandall, K. A. (1998). MODELTEST: testing the model of DNA substitution.. *Bioinformatics*, 14(9), 817-818. https://doi.org/10.1093/bioinformatics/14.9.817
- BibTeX key: `posada1998Modeltest`
- Identifier: DOI [10.1093/bioinformatics/14.9.817](https://doi.org/10.1093/bioinformatics/14.9.817)
- Publication type: methods paper
- Topical category: Phylogenetic model choice
- Annotation: MODELTEST is an older but foundational citation for formal model selection in phylogenetics. It supports the general principle that substitution models should be chosen and reported rather than assumed.
- How this citation might be used: methods; robustness justification
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 89. Darriba, D., Taboada, G. L., Doallo, R., & Posada, D. (2011). ProtTest 3: fast selection of best-fit models of protein evolution. *Bioinformatics*, 27(8), 1164-1165. https://doi.org/10.1093/bioinformatics/btr088
- BibTeX key: `darriba2011Prottest3`
- Identifier: DOI [10.1093/bioinformatics/btr088](https://doi.org/10.1093/bioinformatics/btr088); PMID [21335321](https://pubmed.ncbi.nlm.nih.gov/21335321/)
- Publication type: methods paper
- Topical category: Protein model choice
- Annotation: ProtTest 3 is directly relevant to selecting protein evolution models for HA amino-acid alignments. It can support model-selection reporting in the rebuilt methods section.
- How this citation might be used: methods; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 21335321 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 90. LEGENDRE, P., & FORTIN, M. (2010). Comparison of the Mantel test and alternative approaches for detecting complex multivariate relationships in the spatial analysis of genetic data. *Molecular Ecology Resources*, 10(5), 831-844. https://doi.org/10.1111/j.1755-0998.2010.02866.x
- BibTeX key: `legendre2010MantelAlternatives`
- Identifier: DOI [10.1111/j.1755-0998.2010.02866.x](https://doi.org/10.1111/j.1755-0998.2010.02866.x)
- Publication type: methods paper
- Topical category: Distance-matrix comparison
- Annotation: This paper compares the Mantel test with alternative approaches for matrix relationships. It is a key citation for replacing ordinary independent-observation correlations with permutation-aware matrix comparisons.
- How this citation might be used: methods; limitations; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 91. Guillot, G., & Rousset, F. (2013). Dismantling the Mantel tests. *Methods in Ecology and Evolution*, 4(4), 336-344. https://doi.org/10.1111/2041-210X.12018
- BibTeX key: `guillot2013DismantlingMantel`
- Identifier: DOI [10.1111/2041-210X.12018](https://doi.org/10.1111/2041-210X.12018)
- Publication type: methods paper
- Topical category: Mantel-test limitations
- Annotation: This critique highlights limitations of Mantel tests under spatial/genetic data structures. It should be cited to make the manuscript's matrix-comparison inference appropriately cautious rather than presenting Mantel as a magic fix.
- How this citation might be used: limitations; robustness justification
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 92. Legendre, P., Fortin, M., & Borcard, D. (2015). Should the Mantel test be used in spatial analysis?. *Methods in Ecology and Evolution*, 6(11), 1239-1247. https://doi.org/10.1111/2041-210X.12425
- BibTeX key: `legendre2015ShouldMantel`
- Identifier: DOI [10.1111/2041-210X.12425](https://doi.org/10.1111/2041-210X.12425)
- Publication type: methods paper
- Topical category: Mantel-test use and limitations
- Annotation: This paper discusses when Mantel tests should and should not be used in spatial analysis. It provides nuance for justifying Mantel-style tests while also acknowledging limitations and possible alternatives.
- How this citation might be used: methods; limitations
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 93. Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple Regression and Correlation Extensions of the Mantel Test of Matrix Correspondence. *Systematic Zoology*, 35(4), 627. https://doi.org/10.2307/2413122
- BibTeX key: `smouse1986MantelExtensions`
- Identifier: DOI [10.2307/2413122](https://doi.org/10.2307/2413122)
- Publication type: seminal older source
- Topical category: Matrix correlation and permutation tests
- Annotation: This paper extends Mantel-style matrix correspondence tests to multiple regression and correlation settings. It is useful if the project considers multivariable or partial matrix-comparison approaches.
- How this citation might be used: methods; robustness justification
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 94. Harmon, L. J., & Glor, R. E. (2010). POOR STATISTICAL PERFORMANCE OF THE MANTEL TEST IN PHYLOGENETIC COMPARATIVE ANALYSES. *Evolution*. https://doi.org/10.1111/j.1558-5646.2010.00973.x
- BibTeX key: `harmon2010MantelPhylogeneticComparative`
- Identifier: DOI [10.1111/j.1558-5646.2010.00973.x](https://doi.org/10.1111/j.1558-5646.2010.00973.x)
- Publication type: methods paper
- Topical category: Mantel limitations in phylogenetics
- Annotation: This paper shows poor statistical performance of Mantel tests in phylogenetic comparative analyses. It is a useful cautionary citation to avoid over-interpreting matrix correlations as definitive evidence.
- How this citation might be used: limitations; robustness justification
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 95. Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test. *Oecologia*, 129(2), 169-178. https://doi.org/10.1007/s004420100720
- BibTeX key: `peresNeto2001ProcrustesMantel`
- Identifier: DOI [10.1007/s004420100720](https://doi.org/10.1007/s004420100720)
- Publication type: methods paper
- Topical category: Matrix comparison alternatives
- Annotation: This paper compares Procrustean superimposition with Mantel approaches for matching multivariate datasets. It is relevant if the analysis later compares antigenic maps, ordinations, or distance matrices with methods beyond simple correlations.
- How this citation might be used: robustness justification; limitations
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 96. Anderson, M. J. (2001). A new method for non‐parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32-46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x
- BibTeX key: `anderson2001Permanova`
- Identifier: DOI [10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x)
- Publication type: methods paper
- Topical category: Permutation-based multivariate inference
- Annotation: This paper introduces permutation-based multivariate analysis of variance. It is lower priority but useful background for permutation logic and dependence-aware inference in matrix-derived data.
- How this citation might be used: methods background; robustness justification
- Priority: low
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 97. Hadfield, J., Megill, C., Bell, S. M., Huddleston, J., Potter, B., Callender, C., et al. (2018). Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics*, 34(23), 4121-4123. https://doi.org/10.1093/bioinformatics/bty407
- BibTeX key: `hadfield2018Nextstrain`
- Identifier: DOI [10.1093/bioinformatics/bty407](https://doi.org/10.1093/bioinformatics/bty407); PMID [29790939](https://pubmed.ncbi.nlm.nih.gov/29790939/)
- Publication type: dataset/resource
- Topical category: Influenza phylogenetic surveillance resources
- Annotation: Nextstrain is a real-time pathogen-evolution platform that integrates sequence data, phylogenies, and metadata. It provides resource context for modern influenza phylogenetics and for transparent, reproducible visualization of viral evolution.
- How this citation might be used: background; reproducibility context
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 29790939 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 98. Shu, Y., & McCauley, J. (2017). GISAID: Global initiative on sharing all influenza data – from vision to reality. *Eurosurveillance*, 22(13). https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494
- BibTeX key: `shu2017Gisaid`
- Identifier: DOI [10.2807/1560-7917.ES.2017.22.13.30494](https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494); PMID [28382917](https://pubmed.ncbi.nlm.nih.gov/28382917/)
- Publication type: dataset/resource
- Topical category: Influenza sequence data sharing
- Annotation: This paper describes GISAID's role in sharing influenza data. It is important for provenance and data-availability discussions if the manuscript's sequences or accession history involve GISAID-derived records.
- How this citation might be used: methods context; data provenance
- Priority: high
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 28382917 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

### 99. Landau, W. (2021). The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing. *Journal of Open Source Software*, 6(57), 2959. https://doi.org/10.21105/joss.02959
- BibTeX key: `landau2021Targets`
- Identifier: DOI [10.21105/joss.02959](https://doi.org/10.21105/joss.02959)
- Publication type: methods paper
- Topical category: Reproducible workflow
- Annotation: This paper documents the `targets` R package for reproducible, pipeline-oriented analysis. It supports the decision-log direction to rebuild the project as an auditable pipeline from inputs to manuscript-ready outputs.
- How this citation might be used: reproducibility; methods
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found on 2026-05-29; no PubMed ID-conversion match in batch scan, so recheck PubMed/Retraction Watch before final submission.

### 100. Leebens-Mack, J., Vision, T., Brenner, E., Bowers, J. E., Cannon, S., Clement, M. J., et al. (2006). Taking the First Steps towards a Standard for Reporting on Phylogenies: Minimum Information about a Phylogenetic Analysis (MIAPA). *OMICS: A Journal of Integrative Biology*, 10(2), 231-237. https://doi.org/10.1089/omi.2006.10.231
- BibTeX key: `leebensMack2006Miapa`
- Identifier: DOI [10.1089/omi.2006.10.231](https://doi.org/10.1089/omi.2006.10.231); PMID [16901231](https://pubmed.ncbi.nlm.nih.gov/16901231/)
- Publication type: guideline
- Topical category: Phylogenetic reporting
- Annotation: MIAPA provides minimum-information guidance for reporting phylogenetic analyses. It is useful for ensuring the rebuilt manuscript reports alignment, model choice, support, and tree-comparison details transparently.
- How this citation might be used: methods; reporting
- Priority: medium
- Verification status: Verified against Crossref DOI metadata on 2026-05-29.
- Retraction check status: No Crossref retraction/update relation found and no PubMed ESummary retraction-type flag found for PMID 16901231 on 2026-05-29; recheck PubMed/Retraction Watch before final submission.

## Notes for Later Citation Work

- Before integrating these entries into the manuscript bibliography, deduplicate against `products/project-refs.bib` and preserve any existing citation keys already used in `products/manuscript.qmd`.
- Re-run citation verification after the manuscript scope is finalized, especially for the preprint and for any claims about current WHO vaccine strain selection policy.
- Convert any broad background citation into a claim-specific citation only after reading the relevant full text or abstract in context.
