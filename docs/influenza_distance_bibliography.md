# Influenza antigenic distance metrics: curated annotated bibliography

**Purpose.** This recreated bibliography supports expansion of the introduction and discussion for a manuscript comparing temporal, HA sequence, epitope-focused, antigenic-cartographic, and phylogenetic distance metrics for influenza A/H1N1 and A/H3N2. It is organized to help distinguish viral genetic evolution, phylogenetic relatedness, antigenic phenotype, assay-derived antibody reactivity, immune history, and vaccine-selection relevance.

**Source count.** 105 total entries, including biological literature, statistical methodology, and core software citations.

**Quality-control note.** Citation keys in this Markdown file are generated from the same structured source list as the BibTeX file. Entries with incomplete or uncertain metadata are explicitly marked in the citation note; no DOI was intentionally included where I was not confident enough to reproduce it.

## Category overview

- **Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics:** 15 entries
- **Antigenic cartography and empirical antigenic distance:** 13 entries
- **Sequence-based antigenic distance and antigenic phenotype prediction:** 16 entries
- **Phylogenetic and phylodynamic methods for influenza:** 12 entries
- **Immune history, original antigenic sin, imprinting, and heterologous vaccine response:** 14 entries
- **HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype:** 11 entries
- **Vaccine strain selection, universal influenza vaccine design, and surveillance implications:** 11 entries
- **Statistical and methodological issues for distance comparison:** 10 entries
- **Software and computational tools:** 3 entries

## Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics

### Kim2018DriftingShifting

**Citation:** Kim H, Webster RG, Webby RJ. Influenza Virus: Dealing with a Drifting and Shifting Pathogen. Viral Immunology. 2018;31(2):174--183.  
**BibTeX key:** `Kim2018DriftingShifting`  
**DOI/stable URL:** doi:10.1089/vim.2017.0141  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** A concise review of antigenic drift, antigenic shift, and the biological reasons influenza vaccines require continual updating. It is useful as broad introductory framing, but it should not be used as the main evidence for any specific distance-metric result.  
**Specific manuscript use:** Introduction  
**Claim supported:** Influenza A evolution is shaped by continual antigenic drift and occasional major lineage changes.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Petrova2018SeasonalEvolution

**Citation:** Petrova VN, Russell CA. The evolution of seasonal influenza viruses. Nature Reviews Microbiology. 2018;16(1):47--60.  
**BibTeX key:** `Petrova2018SeasonalEvolution`  
**DOI/stable URL:** doi:10.1038/s41579-017-0006-x  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This review summarizes seasonal influenza evolution, including immune selection, antigenic drift, reassortment, and subtype-specific differences. It is one of the best high-level sources for expanding the introduction without overstating the manuscript’s own empirical findings.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Seasonal influenza evolution reflects interacting genetic, antigenic, epidemiologic, and immunologic processes.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Ferguson2003EcologicalImmunological

**Citation:** Ferguson NM, Galvani AP, Bush RM. Ecological and immunological determinants of influenza evolution. Nature. 2003;422(6930):428--433.  
**BibTeX key:** `Ferguson2003EcologicalImmunological`  
**DOI/stable URL:** doi:10.1038/nature01509  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This foundational modeling paper links host immunity, viral ecology, and antigenic evolution. It is methodologically dated relative to genomic surveillance papers, but remains useful for motivating why immune selection can shape observed evolutionary trajectories.  
**Specific manuscript use:** Introduction  
**Claim supported:** Influenza evolutionary dynamics depend on immune selection as well as viral mutation.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Koelle2006EpochalEvolution

**Citation:** Koelle K, Cobey S, Grenfell B, Pascual M. Epochal evolution shapes the phylodynamics of interpandemic influenza A (H3N2) in humans. Science. 2006;314(5807):1898--1903.  
**BibTeX key:** `Koelle2006EpochalEvolution`  
**DOI/stable URL:** doi:10.1126/science.1132745  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This paper introduced an influential view of H3N2 evolution as punctuated by antigenic cluster transitions rather than smooth drift alone. It is directly relevant to interpreting why H3N2 temporal, sequence, and phylogenetic distances may align better than they do for H1N1 over some intervals.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** H3N2 evolution can be organized by antigenic epochs that structure viral population dynamics.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Rambaut2008GenomicEpidemiology

**Citation:** Rambaut A, Pybus OG, Nelson MI, Viboud C, Taubenberger JK, Holmes EC. The genomic and epidemiological dynamics of human influenza A virus. Nature. 2008;453(7195):615--619.  
**BibTeX key:** `Rambaut2008GenomicEpidemiology`  
**DOI/stable URL:** doi:10.1038/nature06945  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This genomic epidemiology paper is useful for contrasting H3N2’s more ladder-like global dynamics with the more complex histories of other influenza A subtypes. It supports the manuscript’s emphasis that subtype-specific evolutionary structure matters when interpreting distance metrics.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Human influenza A subtypes differ in the genomic and epidemiologic dynamics that shape their phylogenies.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Russell2008GlobalCirculationH3N2

**Citation:** Russell CA, Jones TC, Barr IG, Cox NJ, Garten RJ, Gregory V, et al.. The global circulation of seasonal influenza A (H3N2) viruses. Science. 2008;320(5874):340--346.  
**BibTeX key:** `Russell2008GlobalCirculationH3N2`  
**DOI/stable URL:** doi:10.1126/science.1154137  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This paper is central for discussing H3N2’s global migration and replacement dynamics. It helps justify why temporal and phylogenetic distance can sometimes appear more coherent for H3N2 than for H1N1, while still remaining an incomplete proxy for antigenic phenotype.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Seasonal H3N2 evolution has a structured global circulation pattern that can produce relatively ladder-like phylogenies.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Nelson2007StochasticInfluenzaEvolution

**Citation:** Nelson MI, Simonsen L, Viboud C, Miller MA, Taylor J, George KS, et al.. Stochastic processes are key determinants of short-term evolution in influenza A virus. PLoS Pathogens. 2007;3(12):e125.  
**BibTeX key:** `Nelson2007StochasticInfluenzaEvolution`  
**DOI/stable URL:** doi:10.1371/journal.ppat.0030125  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This paper is useful for tempering deterministic language about influenza evolution. It supports a more careful discussion in which immune selection is important, but observed lineages also reflect stochastic transmission, migration, and sampling effects.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Short-term influenza phylogenies reflect stochastic processes as well as selection.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies.

### Nelson2008H1N1Reassortment

**Citation:** Nelson MI, Viboud C, Simonsen L, Bennett RT, Griesemer SB, St George K, et al.. Multiple reassortment events in the evolutionary history of H1N1 influenza A virus since 1918. PLoS Pathogens. 2008;4(2):e1000012.  
**BibTeX key:** `Nelson2008H1N1Reassortment`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1000012  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This source is particularly important for the manuscript’s H1N1 discussion because it documents a complex evolutionary history since 1918. It helps support the claim that simple temporal distance can fail when lineage history is discontinuous or shaped by reassortment and replacement.  
**Specific manuscript use:** Discussion  
**Claim supported:** H1N1 evolutionary history contains reassortment and lineage changes that can decouple time from genetic distance.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Rozo2015Reemergent1977H1N1

**Citation:** Rozo M, Gronvall GK. The reemergent 1977 H1N1 strain and the gain-of-function debate. mBio. 2015;6(4):e01013-15.  
**BibTeX key:** `Rozo2015Reemergent1977H1N1`  
**DOI/stable URL:** doi:10.1128/mBio.01013-15  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This review is useful for introducing the unusual 1977 H1N1 re-emergence and why influenza lineage history can violate naive temporal assumptions. It should be cited carefully because its focus is biosafety debate rather than distance metrics.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** The 1977 H1N1 re-emergence is a historical example in which temporal circulation does not map cleanly onto ordinary evolutionary expectations.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.

### Taubenberger2005Characterization1918

**Citation:** Taubenberger JK, Reid AH, Lourens RM, Wang R, Jin G, Fanning TG. Characterization of the 1918 influenza virus polymerase genes. Nature. 2005;437(7060):889--893.  
**BibTeX key:** `Taubenberger2005Characterization1918`  
**DOI/stable URL:** doi:10.1038/nature04230  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This source is useful background for the 1918 lineage context, though it is not directly about antigenic distance. It supports discussion of why 1918-related ancestry matters when interpreting H1N1 phylogenetic structure.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** The 1918 pandemic virus has a central role in the ancestry of later influenza A lineages.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports indirectly.

### Garten2009SwineOriginH1N1

**Citation:** Garten RJ, Davis CT, Russell CA, Shu B, Lindstrom S, Balish A, et al.. Antigenic and genetic characteristics of swine-origin 2009 A(H1N1) influenza viruses circulating in humans. Science. 2009;325(5937):197--201.  
**BibTeX key:** `Garten2009SwineOriginH1N1`  
**DOI/stable URL:** doi:10.1126/science.1176225  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This is a core paper for explaining why pandemic H1N1 belongs to a different evolutionary context than pre-2009 seasonal H1N1. It directly supports the manuscript’s point that temporal distance can be misleading for H1N1.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** The 2009 H1N1 pandemic virus had swine-origin genetic and antigenic characteristics that disrupted seasonal H1N1 continuity.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Dawood2009EmergenceNovelH1N1

**Citation:** Dawood FS, Jain S, Finelli L, Shaw MW, Lindstrom S, Garten RJ, et al.. Emergence of a novel swine-origin influenza A (H1N1) virus in humans. New England Journal of Medicine. 2009;360(25):2605--2615.  
**BibTeX key:** `Dawood2009EmergenceNovelH1N1`  
**DOI/stable URL:** doi:10.1056/NEJMoa0903810  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This clinical and epidemiologic report anchors the 2009 pandemic context. It is less directly methodological than Garten et al., but useful for explaining why H1N1 temporal continuity is biologically and epidemiologically broken around 2009.  
**Specific manuscript use:** Introduction  
**Claim supported:** The 2009 pandemic introduced a novel swine-origin H1N1 lineage into humans.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Smith2009Origins2009H1N1

**Citation:** Smith GJD, Vijaykrishna D, Bahl J, Lycett SJ, Worobey M, Pybus OG, et al.. Origins and evolutionary genomics of the 2009 swine-origin H1N1 influenza A epidemic. Nature. 2009;459(7250):1122--1125.  
**BibTeX key:** `Smith2009Origins2009H1N1`  
**DOI/stable URL:** doi:10.1038/nature08182  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This paper gives the evolutionary-genomic background for the 2009 H1N1 emergence. It is highly relevant to the manuscript’s observation that 2009-like strains are not temporally close to older strains despite genetic relationships through swine lineages.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** The 2009 H1N1 pandemic lineage arose from complex swine influenza ancestry rather than simple continuation of recent human seasonal H1N1.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Bedford2015GlobalMigration

**Citation:** Bedford T, Riley S, Barr IG, Broor S, Chadha M, Cox NJ, et al.. Global circulation patterns of seasonal influenza viruses vary with antigenic drift. Nature. 2015;523(7559):217--220.  
**BibTeX key:** `Bedford2015GlobalMigration`  
**DOI/stable URL:** doi:10.1038/nature14460  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This source connects antigenic drift to global circulation patterns across influenza types and subtypes. It is useful for showing that the relationship between antigenic change, phylogeny, and time depends on subtype and surveillance context.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Influenza circulation and lineage persistence vary with antigenic drift and subtype-specific evolutionary dynamics.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Bush1999PositiveSelectionHA

**Citation:** Bush RM, Bender CA, Subbarao K, Cox NJ, Fitch WM. Predicting the evolution of human influenza A. Science. 1999;286(5446):1921--1925.  
**BibTeX key:** `Bush1999PositiveSelectionHA`  
**DOI/stable URL:** doi:10.1126/science.286.5446.1921  
**Category/categories:** Influenza evolution, antigenic drift, and subtype-specific evolutionary dynamics  
**Annotation:** This classic paper links positive selection in HA to future viral success. It is historically important and helps frame why HA sequence contains predictive signal, but newer phenotype-linked models should carry the main burden for modern claims.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Amino acid changes under positive selection in HA can carry information about future influenza evolution.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

## Antigenic cartography and empirical antigenic distance

### Smith2004AntigenicCartography

**Citation:** Smith DJ, Lapedes AS, de Jong JC, Bestebroer TM, Rimmelzwaan GF, Osterhaus ADME, et al.. Mapping the antigenic and genetic evolution of influenza virus. Science. 2004;305(5682):371--376.  
**BibTeX key:** `Smith2004AntigenicCartography`  
**DOI/stable URL:** doi:10.1126/science.1097211  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This is the foundational antigenic-cartography paper for influenza. It should be cited early to define cartographic distance as an empirical summary of serologic cross-reactivity rather than a direct transformation of sequence or phylogeny.  
**Specific manuscript use:** Introduction; Methods context  
**Claim supported:** Antigenic cartography can summarize HAI cross-reactivity into a spatial map of antigenic relationships.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Fonville2014AntibodyLandscapes

**Citation:** Fonville JM, Wilks SH, James SL, Fox A, Ventresca M, Aban M, et al.. Antibody landscapes after influenza virus infection or vaccination. Science. 2014;346(6212):996--1000.  
**BibTeX key:** `Fonville2014AntibodyLandscapes`  
**DOI/stable URL:** doi:10.1126/science.1256427  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This paper extends cartographic thinking from virus maps to antibody landscapes across exposure histories. It is crucial for distinguishing viral antigenic phenotype from human immune-response landscapes.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Antibody responses to influenza are shaped by prior exposure and can be represented as landscapes over antigenic space.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Fonville2016HumanAntiseraH3N2

**Citation:** Fonville JM, Fraaij PLA, de Mutsert G, Wilks SH, van Beek R, Fouchier RAM, et al.. Antigenic maps of influenza A(H3N2) produced with human antisera obtained after primary infection. Journal of Infectious Diseases. 2016;213(1):31--38.  
**BibTeX key:** `Fonville2016HumanAntiseraH3N2`  
**DOI/stable URL:** doi:10.1093/infdis/jiv367  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This human-sera antigenic map is directly relevant because the manuscript uses human immune response data rather than only ferret antisera. It supports a nuanced discussion that maps based on human sera may differ from canonical surveillance maps.  
**Specific manuscript use:** Introduction; Discussion; Limitations  
**Claim supported:** Human-sera antigenic maps can capture exposure-conditioned immune reactivity not identical to ferret-sera surveillance maps.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Bedford2014AntigenicDynamics

**Citation:** Bedford T, Suchard MA, Lemey P, Dudas G, Gregory V, Hay AJ, et al.. Integrating influenza antigenic dynamics with molecular evolution. eLife. 2014;3:e01914.  
**BibTeX key:** `Bedford2014AntigenicDynamics`  
**DOI/stable URL:** doi:10.7554/eLife.01914  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This paper directly connects antigenic cartography with molecular evolutionary modeling. It is highly relevant to the manuscript because it shows both the value and the complexity of integrating antigenic and phylogenetic data.  
**Specific manuscript use:** Introduction; Discussion; Methods context  
**Claim supported:** Antigenic and molecular evolution can be integrated, but they represent distinct and only partially overlapping structures.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Luksza2014PredictiveFitness

**Citation:** Łuksza M, Lässig M. A predictive fitness model for influenza. Nature. 2014;507(7490):57--61.  
**BibTeX key:** `Luksza2014PredictiveFitness`  
**DOI/stable URL:** doi:10.1038/nature13087  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This influential model uses antigenic novelty and mutational load to predict influenza lineage success. It is useful for showing that antigenic phenotype contributes to evolutionary forecasting beyond raw genetic distance.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Evolutionary forecasting improves when antigenic novelty is modeled explicitly rather than ignored.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Neher2016AntigenicPhenotypes

**Citation:** Neher RA, Bedford T, Daniels RS, Russell CA, Shraiman BI. Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza viruses. Proceedings of the National Academy of Sciences. 2016;113(12):E1701--E1709.  
**BibTeX key:** `Neher2016AntigenicPhenotypes`  
**DOI/stable URL:** doi:10.1073/pnas.1525578113  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This paper is important for connecting antigenic maps, prediction, and visualization. It can help the discussion explain that antigenic phenotype is a model-estimated object with uncertainty, not an directly observed scalar property.  
**Specific manuscript use:** Discussion; Methods context; Future work  
**Claim supported:** Antigenic phenotype can be predicted and visualized, but requires phenotype-linked data and model assumptions.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Cai2010ComputationalCartography

**Citation:** Cai Z, Zhang T, Wan XF. A computational framework for influenza antigenic cartography. PLoS Computational Biology. 2010;6(10):e1000949.  
**BibTeX key:** `Cai2010ComputationalCartography`  
**DOI/stable URL:** doi:10.1371/journal.pcbi.1000949  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This methods paper provides computational context for constructing antigenic maps. It is useful background for methods discussion, especially if the manuscript expands on how map distances are estimated and why dimensionality or sparse titers matter.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Computational choices affect how HAI data are transformed into antigenic map distances.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.

### Sun2013SequenceAntigenicity

**Citation:** Sun H, Yang J, Zhang T, Long LP, Jia K, Yang G, et al.. Using sequence data to infer the antigenicity of influenza virus. mBio. 2013;4(4):e00230-13.  
**BibTeX key:** `Sun2013SequenceAntigenicity`  
**DOI/stable URL:** doi:10.1128/mBio.00230-13  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This source bridges antigenic cartography and sequence-based inference. It is useful because it shows that antigenic relationships can sometimes be inferred from sequence, while still requiring phenotypic training or validation.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Sequence data can help infer antigenicity but should be anchored against serologic measurements.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Li2013H3N2ChinaMapping

**Citation:** Li X, Liu Y, Yu A, Feng J, Zhang C, Yang J, et al.. Mapping of H3N2 influenza antigenic evolution in China reveals a strategy for vaccine strain recommendation. Nature Communications. 2013;4:2922.  
**BibTeX key:** `Li2013H3N2ChinaMapping`  
**DOI/stable URL:** doi:10.1038/ncomms3922  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This regional antigenic-mapping paper is useful for connecting cartographic methods to vaccine strain recommendation. It also reminds readers that geographic sampling and surveillance design can affect apparent antigenic structure.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Regional antigenic mapping can inform vaccine strain decisions but depends on surveillance sampling.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Harvey2016HighLowImpactSubstitutions

**Citation:** Harvey WT, Benton DJ, Gregory V, Hall JPJ, Daniels RS, Bedford T, et al.. Identification of low- and high-impact hemagglutinin amino acid substitutions that drive antigenic drift of influenza A(H1N1) viruses. PLoS Pathogens. 2016;12(4):e1005526.  
**BibTeX key:** `Harvey2016HighLowImpactSubstitutions`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1005526  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This paper provides mechanistic context for why not every HA substitution has equal antigenic impact. It is useful for explaining why Hamming distance and epitope distance can be crude unless mutations are weighted by functional effect.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Antigenic drift can be driven by a limited subset of HA substitutions with heterogeneous effects.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Hensley2009ReceptorBindingAvidity

**Citation:** Hensley SE, Das SR, Bailey AL, Schmidt LM, Hickman HD, Jayaraman A, et al.. Hemagglutinin receptor binding avidity drives influenza A virus antigenic drift. Science. 2009;326(5953):734--736.  
**BibTeX key:** `Hensley2009ReceptorBindingAvidity`  
**DOI/stable URL:** doi:10.1126/science.1178258  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This mechanistic paper is important because it shows that measured antigenicity can be affected by receptor-binding avidity, not just antibody escape at canonical sites. It supports a limitation section noting that HAI-derived cartographic distance is a composite assay phenotype.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Changes in receptor-binding avidity can alter measured antigenicity and complicate interpretation of HA sequence distance.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies.

### Igarashi2010PandemicH1N1AntigenicStructure

**Citation:** Igarashi M, Ito K, Yoshida R, Tomabechi D, Kida H, Takada A. Predicting the antigenic structure of the pandemic (H1N1) 2009 influenza virus hemagglutinin. PLoS ONE. 2010;5(1):e8553.  
**BibTeX key:** `Igarashi2010PandemicH1N1AntigenicStructure`  
**DOI/stable URL:** doi:10.1371/journal.pone.0008553  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** This paper is useful for the H1N1-specific framing because it compares pandemic H1N1 antigenic structure to earlier H1 viruses. It supports the biological plausibility of the manuscript’s H1N1 temporal-distance result.  
**Specific manuscript use:** Discussion  
**Claim supported:** Pandemic H1N1 HA antigenic structure reflects relationships to older H1 lineages that temporal order alone may obscure.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports.

### Wilks2023Racmacs

**Citation:** Wilks SH. Racmacs: antigenic cartography macros. Software. 2023.  
**BibTeX key:** `Wilks2023Racmacs`  
**DOI/stable URL:** https://acorg.github.io/Racmacs/  
**Category/categories:** Antigenic cartography and empirical antigenic distance  
**Annotation:** Racmacs is the software context for the manuscript’s antigenic-cartography computations. It should be cited in methods, separate from biological evidence, because it documents the tool rather than supporting a biological claim.  
**Specific manuscript use:** Methods context  
**Claim supported:** Racmacs implements antigenic cartography workflows used to estimate map-based antigenic distances.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.
  **Metadata note:** R package/software documentation

## Sequence-based antigenic distance and antigenic phenotype prediction

### Gupta2006VaccineEfficacyDistance

**Citation:** Gupta V, Earl DJ, Deem MW. Quantifying influenza vaccine efficacy and antigenic distance. Vaccine. 2006;24(18):3881--3888.  
**BibTeX key:** `Gupta2006VaccineEfficacyDistance`  
**DOI/stable URL:** doi:10.1016/j.vaccine.2006.01.010  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This is a central source for p-epitope-style antigenic distance and vaccine effectiveness framing. It is useful historically and methodologically, but its assumptions should be presented as one model family rather than as definitive measurement of antigenic phenotype.  
**Specific manuscript use:** Introduction; Methods context; Discussion  
**Claim supported:** Epitope-focused sequence distances can be associated with vaccine effectiveness and antigenic mismatch.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Lee2004PredictingH3N2Variants

**Citation:** Lee MS, Chen JS. Predicting antigenic variants of influenza A/H3N2 viruses. Emerging Infectious Diseases. 2004;10(8):1385--1390.  
**BibTeX key:** `Lee2004PredictingH3N2Variants`  
**DOI/stable URL:** doi:10.3201/eid1008.040107  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This early sequence-based prediction paper is useful for historical context on identifying antigenic variants from HA substitutions. It supports the idea that sequence contains antigenic signal, while newer models should be cited for current performance.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** HA substitutions can be used to predict some H3N2 antigenic variants.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Anderson2018SequenceBasedCartographyH1N1

**Citation:** Anderson CS, McCall PR, Stern HA, Yang H, Topham DJ. Antigenic cartography of H1N1 influenza viruses using sequence-based antigenic distance calculation. BMC Bioinformatics. 2018;19(1):51.  
**BibTeX key:** `Anderson2018SequenceBasedCartographyH1N1`  
**DOI/stable URL:** doi:10.1186/s12859-018-2042-4  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This paper is directly relevant because it combines H1N1, sequence-based distance, and antigenic cartography. It provides a close precedent for the manuscript’s use of H1N1 sequence-derived distances, while also showing that sequence-based maps are not the same as empirical HAI maps.  
**Specific manuscript use:** Introduction; Methods context; Discussion  
**Claim supported:** Sequence-derived antigenic distances can be used to construct H1N1 antigenic maps, but they remain model-based approximations.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Yao2017RandomForestAntigenicity

**Citation:** Yao Y, Li X, Liao B, et al.. Predicting influenza antigenicity from hemagglutinin sequence data based on a joint random forest method. Scientific Reports. 2017;7:1545.  
**BibTeX key:** `Yao2017RandomForestAntigenicity`  
**DOI/stable URL:** doi:10.1038/s41598-017-01699-z  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This machine-learning paper predicts antigenic relationships from HA sequence features rather than using unweighted Hamming distance. It is a good source for arguing that sequence-based prediction can improve when substitutions are modeled flexibly.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Machine-learning models can extract antigenic signal from HA sequence beyond simple mismatch counts.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Yin2018H1N1StackingModel

**Citation:** Yin R, Tran VH, Zhou X, Zheng J, Kwoh CK. Predicting antigenic variants of H1N1 influenza virus based on epidemics and pandemics using a stacking model. PLoS ONE. 2018;13(11):e0207777.  
**BibTeX key:** `Yin2018H1N1StackingModel`  
**DOI/stable URL:** doi:10.1371/journal.pone.0207777  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This H1N1-specific paper is useful because the manuscript finds H1N1 temporal distance especially problematic. It shows that H1N1 antigenic variants can still be modeled from sequence and epidemiologic structure even when simple time is inadequate.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Subtype-specific models can predict H1N1 antigenic variants better than naive temporal assumptions.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Peng2017PREDAVFluA

**Citation:** Peng Y, Wu A, Zhang Y, et al.. A universal computational model for predicting antigenic variants of influenza A virus based on conserved antigenic structures. Scientific Reports. 2017;7:42051.  
**BibTeX key:** `Peng2017PREDAVFluA`  
**DOI/stable URL:** doi:10.1038/srep42051  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This model is relevant for broad sequence-based antigenic prediction across influenza A subtypes. It is useful for future-work framing, but should be used cautiously because “universal” computational performance depends on training data and validation design.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Conserved antigenic structure can be used to support computational prediction of influenza A antigenic variants.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.

### Xia2021DeepLearningAntigenicVariation

**Citation:** Xia YL, et al.. A deep learning approach for predicting antigenic variation of influenza A H3N2. Briefings in Bioinformatics. 2021.  
**BibTeX key:** `Xia2021DeepLearningAntigenicVariation`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/34697557/  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This source represents the transition from hand-crafted sequence distances to deep learning models for H3N2 antigenic variation. In the manuscript it can be used to qualify broad claims by noting that simple distance metrics may understate the predictive information available in sequence.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Deep learning approaches can model H3N2 antigenic variation from HA sequence features.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.
  **Metadata note:** Metadata should be rechecked before final manuscript use; DOI intentionally omitted in this recreated file.

### Peng2023AttributeNetworkEmbedding

**Citation:** Peng F, Xia Y, Li W. Prediction of antigenic distance in influenza A using attribute network embedding. Viruses. 2023;15(7):1478.  
**BibTeX key:** `Peng2023AttributeNetworkEmbedding`  
**DOI/stable URL:** doi:10.3390/v15071478  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This paper is directly relevant because it targets quantitative antigenic distance rather than only variant classification. It provides a recent comparison point for a manuscript built around distance matrices.  
**Specific manuscript use:** Methods context; Discussion; Future work  
**Claim supported:** Recent models can predict continuous influenza antigenic distances from structured sequence-derived representations.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Meng2024PREDACCNN

**Citation:** Meng J, Liu J, Song W, Li H, Wang J, Zhang L, et al.. PREDAC-CNN: predicting antigenic clusters of seasonal influenza A viruses with convolutional neural network. Briefings in Bioinformatics. 2024;25(2):bbae033.  
**BibTeX key:** `Meng2024PREDACCNN`  
**DOI/stable URL:** doi:10.1093/bib/bbae033  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This paper is useful because it predicts antigenic clusters for seasonal influenza A viruses using deep learning. It helps frame H1N1 and H3N2 as subtype-specific antigenic prediction problems rather than a single universal distance problem.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Deep learning can recover antigenic-cluster structure from influenza A sequence data.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Jia2024MetaFluAD

**Citation:** Jia Q, Xia Y, Dong F, Li W. MetaFluAD: meta-learning for predicting antigenic distances among influenza viruses. Briefings in Bioinformatics. 2024;25(5):bbae395.  
**BibTeX key:** `Jia2024MetaFluAD`  
**DOI/stable URL:** doi:10.1093/bib/bbae395  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This is a strong recent source for quantitative antigenic-distance prediction under limited data. It is useful for future work because it treats sparse phenotype-linked data as a central modeling problem.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Meta-learning can improve antigenic-distance prediction when subtype-specific antigenic data are limited.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Shah2024SeasonalAntigenicPrediction

**Citation:** Shah SAW, Palomar DP, Barr I, McKay MR, et al.. Seasonal antigenic prediction of influenza A H3N2 using machine learning. Nature Communications. 2024;15(1):3833.  
**BibTeX key:** `Shah2024SeasonalAntigenicPrediction`  
**DOI/stable URL:** doi:10.1038/s41467-024-47862-9  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This recent paper is one of the most important counterweights to an overly negative view of sequence-based prediction. It shows that H3N2 sequence can carry useful antigenic signal even when simple phylogenetic or Hamming distance does not fully explain cartographic distance.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Machine learning can extract useful seasonal H3N2 antigenic information from genetic data.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies/challenges.

### Durazzi2025LanguageModelsAntigenicity

**Citation:** Durazzi F, Koopmans MPG, Fouchier RAM, et al.. Language models learn to represent antigenic properties of human influenza A(H3) virus. Scientific Reports. 2025;15(1):21364.  
**BibTeX key:** `Durazzi2025LanguageModelsAntigenicity`  
**DOI/stable URL:** doi:10.1038/s41598-025-03275-2  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This paper represents the newest protein-language-model framing for antigenicity. It is best used in future work to note that richer sequence representations may narrow the gap between sequence and antigenic phenotype.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Protein language-model embeddings can encode antigenic information from influenza A(H3) HA sequences.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies/challenges.

### Li2024MFPAD

**Citation:** Li X, et al.. A sequence-based machine learning model for predicting antigenic distance for H3N2 influenza virus. Frontiers in Microbiology. 2024;15:1345794.  
**BibTeX key:** `Li2024MFPAD`  
**DOI/stable URL:** doi:10.3389/fmicb.2024.1345794  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This recent H3N2 paper is directly relevant because it predicts antigenic distance from sequence-derived features. It helps position the manuscript’s simple metrics as a baseline rather than as the endpoint of sequence-based antigenic inference.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Feature-rich sequence models can estimate H3N2 antigenic distances better than unweighted sequence distance alone.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Lou2024SiteTransitionNetworks

**Citation:** Lou J, et al.. Predictive evolutionary modelling for influenza virus by site transition networks. Nature Communications. 2024;15:2381.  
**BibTeX key:** `Lou2024SiteTransitionNetworks`  
**DOI/stable URL:** doi:10.1038/s41467-024-46918-0  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This paper is useful for discussing newer sequence-informed evolutionary forecasting approaches. It supports the idea that genetic data are useful when modeled at the level of structured mutation processes rather than collapsed into a single pairwise distance.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Site-transition models can improve genetic prediction of future influenza evolution.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.

### Pan2012H1N1SequenceDistance

**Citation:** Pan K, Subieta KC, Deem MW. A novel sequence-based antigenic distance measure for H1N1, with application to vaccine effectiveness and the selection of vaccine strains. Protein Engineering, Design and Selection. 2011;24(11):903--911.  
**BibTeX key:** `Pan2012H1N1SequenceDistance`  
**DOI/stable URL:** doi:10.1093/protein/gzr047  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This H1N1-specific extension of sequence-based antigenic distance is directly relevant to the manuscript’s retention of p-epitope-style metrics. It supports discussing why subtype-specific epitope definitions may matter.  
**Specific manuscript use:** Methods context; Discussion  
**Claim supported:** H1N1-specific sequence-based antigenic distance can be defined for vaccine-effectiveness and strain-selection analyses.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Muñoz2015EpitopeBasedVE

**Citation:** Muñoz ET, Deem MW. Epitope analysis for influenza vaccine design. Vaccine. 2005.  
**BibTeX key:** `Muñoz2015EpitopeBasedVE`  
**DOI/stable URL:** https://scholar.google.com/scholar?q=Epitope+analysis+for+influenza+vaccine+design+Deem  
**Category/categories:** Sequence-based antigenic distance and antigenic phenotype prediction  
**Annotation:** This entry is included as background for epitope-oriented vaccine and distance thinking. Verify metadata before using in a manuscript; it is less central than Gupta, Pan, and Anderson for the specific p-epitope distance framing.  
**Specific manuscript use:** Background; Methods context  
**Claim supported:** Epitope-focused analyses can guide influenza vaccine design and antigenic-distance reasoning.  
**Priority score:** Background only  
**Relation to manuscript claims:** Supports indirectly.
  **Metadata note:** Metadata uncertain in recreated file; verify before final manuscript use.

## Phylogenetic and phylodynamic methods for influenza

### Dang2010FLUModel

**Citation:** Dang CC, Le QS, Gascuel O, Le VS. FLU, an amino acid substitution model for influenza proteins. BMC Evolutionary Biology. 2010;10:99.  
**BibTeX key:** `Dang2010FLUModel`  
**DOI/stable URL:** doi:10.1186/1471-2148-10-99  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This is the key citation for the FLU amino acid substitution model used in influenza protein phylogenetics. It should be cited in methods when explaining the maximum-likelihood tree model.  
**Specific manuscript use:** Methods context  
**Claim supported:** The FLU amino acid substitution model was developed specifically for influenza proteins.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

### Guindon2010PhyML

**Citation:** Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology. 2010;59(3):307--321.  
**BibTeX key:** `Guindon2010PhyML`  
**DOI/stable URL:** doi:10.1093/sysbio/syq010  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This source provides general maximum-likelihood phylogenetic context. It is useful background if the methods section explains why ML trees are a model-based comparator rather than a distance-only construction.  
**Specific manuscript use:** Methods context  
**Claim supported:** Maximum-likelihood phylogenetic inference estimates tree structure under an explicit substitution model.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

### Nguyen2015IQTREE

**Citation:** Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution. 2015;32(1):268--274.  
**BibTeX key:** `Nguyen2015IQTREE`  
**DOI/stable URL:** doi:10.1093/molbev/msu300  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** Although the manuscript uses phangorn, IQ-TREE is a useful modern ML phylogenetics reference. It can be used in future-work or methods context when discussing robustness to tree-building software.  
**Specific manuscript use:** Methods context; Future work  
**Claim supported:** Modern ML phylogenetic tools estimate trees efficiently under explicit evolutionary models.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

### Stamatakis2014RAxML

**Citation:** Stamatakis A. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics. 2014;30(9):1312--1313.  
**BibTeX key:** `Stamatakis2014RAxML`  
**DOI/stable URL:** doi:10.1093/bioinformatics/btu033  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** RAxML is another standard reference for ML phylogenetics. It is useful for contextualizing the manuscript’s ML-tree approach, even if not directly used.  
**Specific manuscript use:** Methods context  
**Claim supported:** Maximum-likelihood phylogenetic inference is a standard approach for estimating evolutionary trees from sequence alignments.  
**Priority score:** Background only  
**Relation to manuscript claims:** Methods support.

### Yang1994GammaRates

**Citation:** Yang Z. Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. Journal of Molecular Evolution. 1994;39(3):306--314.  
**BibTeX key:** `Yang1994GammaRates`  
**DOI/stable URL:** doi:10.1007/BF00160154  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This classic paper provides background for among-site rate variation in phylogenetic models. It is relevant if the manuscript explains gamma-distributed rate variation or sensitivity to model assumptions.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Phylogenetic models often account for among-site rate heterogeneity to avoid misspecifying evolutionary distances.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

### Whelan2001WAGModel

**Citation:** Whelan S, Goldman N. A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. Molecular Biology and Evolution. 2001;18(5):691--699.  
**BibTeX key:** `Whelan2001WAGModel`  
**DOI/stable URL:** doi:10.1093/oxfordjournals.molbev.a003851  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This is a general protein substitution-model reference useful for contrasting FLU with broader empirical amino acid models. It helps explain why model choice matters for protein phylogenies.  
**Specific manuscript use:** Methods context  
**Claim supported:** Protein phylogenies depend on empirical amino acid substitution models, which can be general or virus-specific.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

### Shimodaira1999SHTest

**Citation:** Shimodaira H, Hasegawa M. Multiple comparisons of log-likelihoods with applications to phylogenetic inference. Molecular Biology and Evolution. 1999;16(8):1114--1116.  
**BibTeX key:** `Shimodaira1999SHTest`  
**DOI/stable URL:** doi:10.1093/oxfordjournals.molbev.a026201  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This is the original Shimodaira-Hasegawa test reference. It should be cited wherever the manuscript reports SH tests comparing distance-derived trees with ML trees.  
**Specific manuscript use:** Methods context  
**Claim supported:** The SH test compares candidate phylogenetic trees using resampled log-likelihood differences.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

### Shimodaira2002AUTest

**Citation:** Shimodaira H. An approximately unbiased test of phylogenetic tree selection. Systematic Biology. 2002;51(3):492--508.  
**BibTeX key:** `Shimodaira2002AUTest`  
**DOI/stable URL:** doi:10.1080/10635150290069913  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This paper is useful if the manuscript discusses alternative tree-selection tests or limitations of SH testing. It can support a future-work note on stronger or complementary tree-comparison approaches.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Alternative likelihood-based tree-selection tests can address some limitations of earlier multiple-comparison tests.  
**Priority score:** Useful  
**Relation to manuscript claims:** Qualifies.

### Drummond2012BEAST

**Citation:** Drummond AJ, Suchard MA, Xie D, Rambaut A. Bayesian phylogenetics with BEAUti and the BEAST 1.7. Molecular Biology and Evolution. 2012;29(8):1969--1973.  
**BibTeX key:** `Drummond2012BEAST`  
**DOI/stable URL:** doi:10.1093/molbev/mss075  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This source is useful if the discussion mentions phylodynamic alternatives to static ML trees. It helps frame future work using time-resolved trees, coalescent models, and uncertainty propagation.  
**Specific manuscript use:** Future work; Methods context  
**Claim supported:** Bayesian phylogenetic methods can estimate time-resolved evolutionary histories with uncertainty.  
**Priority score:** Useful  
**Relation to manuscript claims:** Future-work support.

### Minin2008Skyride

**Citation:** Minin VN, Bloomquist EW, Suchard MA. Smooth skyride through a rough skyline: Bayesian coalescent-based inference of population dynamics. Molecular Biology and Evolution. 2008;25(7):1459--1471.  
**BibTeX key:** `Minin2008Skyride`  
**DOI/stable URL:** doi:10.1093/molbev/msn090  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This phylodynamic methods paper is useful for contextualizing how influenza population dynamics can be inferred from sequence data. It is indirect for pairwise distance comparison but useful for future work.  
**Specific manuscript use:** Future work  
**Claim supported:** Coalescent phylodynamic models can infer changes in viral population dynamics from sequence data.  
**Priority score:** Background only  
**Relation to manuscript claims:** Future-work support.

### Volz2013ViralPhylodynamics

**Citation:** Volz EM, Koelle K, Bedford T. Viral phylodynamics. PLoS Computational Biology. 2013;9(3):e1002947.  
**BibTeX key:** `Volz2013ViralPhylodynamics`  
**DOI/stable URL:** doi:10.1371/journal.pcbi.1002947  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This review provides broad context for linking viral phylogenies to epidemiologic and evolutionary processes. It helps keep the manuscript from treating tree distance as a purely geometric object detached from population dynamics.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Viral phylogenies encode epidemiologic and evolutionary processes but are not direct measures of antigenic phenotype.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Schliep2011Phangorn

**Citation:** Schliep KP. phangorn: phylogenetic analysis in R. Bioinformatics. 2011;27(4):592--593.  
**BibTeX key:** `Schliep2011Phangorn`  
**DOI/stable URL:** doi:10.1093/bioinformatics/btq706  
**Category/categories:** Phylogenetic and phylodynamic methods for influenza  
**Annotation:** This is the software citation for the R package used for phylogenetic analyses in the manuscript. It should be cited in methods, not used as biological evidence.  
**Specific manuscript use:** Methods context  
**Claim supported:** phangorn implements phylogenetic inference and tree-comparison tools in R.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

## Immune history, original antigenic sin, imprinting, and heterologous vaccine response

### Gostic2016BirthYearImprinting

**Citation:** Gostic KM, Ambrose M, Worobey M, Lloyd-Smith JO. Potent protection against H5N1 and H7N9 influenza via childhood hemagglutinin imprinting. Science. 2016;354(6313):722--726.  
**BibTeX key:** `Gostic2016BirthYearImprinting`  
**DOI/stable URL:** doi:10.1126/science.aag1322  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This landmark paper demonstrates strong birth-cohort effects in protection against avian influenza subtypes. It is important for the manuscript because it shows that immune response is conditioned by host history, not only by viral relatedness.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Childhood influenza imprinting can shape later protection against antigenically related viruses.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Gostic2019ChildhoodImprintingSeasonal

**Citation:** Gostic KM, Bridge R, Brady S, Viboud C, Worobey M, Lloyd-Smith JO. Childhood immune imprinting to influenza A shapes birth year-specific risk during seasonal H1N1 and H3N2 epidemics. PLoS Pathogens. 2019;15(12):e1008109.  
**BibTeX key:** `Gostic2019ChildhoodImprintingSeasonal`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1008109  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper extends imprinting ideas into seasonal H1N1 and H3N2 risk. It is highly relevant for interpreting human HAI data as cohort-conditioned immune measurements.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Birth cohort and childhood imprinting can influence subtype-specific seasonal influenza risk.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Cobey2017ImmuneHistory

**Citation:** Cobey S, Hensley SE. Immune history and influenza virus susceptibility. Current Opinion in Virology. 2017;22:105--111.  
**BibTeX key:** `Cobey2017ImmuneHistory`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/28104444/  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This review is a compact source for immune history, susceptibility, and the interpretation of serologic responses. It should be used to frame why human-sera maps and cohort vaccine studies may not reflect a virus-only phenotype.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Prior infections and vaccinations can shape susceptibility and serologic response to influenza.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.
  **Metadata note:** DOI intentionally omitted in recreated file; verify before final manuscript use.

### Henry2018OriginalAntigenicSin

**Citation:** Henry C, Palm AKE, Krammer F, Wilson PC. From original antigenic sin to the universal influenza virus vaccine. Trends in Immunology. 2018;39(1):70--79.  
**BibTeX key:** `Henry2018OriginalAntigenicSin`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/29242099/  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This review links original antigenic sin, immune imprinting, and universal vaccine design. It is useful for discussion because it connects individual immune history to future vaccine strategies.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Immune imprinting can both constrain and guide design of broadly protective influenza vaccines.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.
  **Metadata note:** DOI intentionally omitted in recreated file; verify before final manuscript use.

### Lessler2012AntigenicSeniority

**Citation:** Lessler J, Riley S, Read JM, Wang S, Zhu H, Smith GJD, et al.. Evidence for antigenic seniority in influenza A (H3N2) antibody responses in southern China. PLoS Pathogens. 2012;8(7):e1002802.  
**BibTeX key:** `Lessler2012AntigenicSeniority`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1002802  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper is important for human HAI interpretation because it shows that antibody responses can preferentially reflect earlier exposures. It directly supports the manuscript’s caution that observed immune-response distances need not match viral phylogeny.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Antibody responses to H3N2 can show antigenic seniority, with earlier exposures disproportionately shaping measured titers.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Linderman2014AtypicalH1N1Immunity

**Citation:** Linderman SL, Chambers BS, Zost SJ, Parkhouse K, Li Y, Herrmann C, et al.. Potential antigenic explanation for atypical H1N1 infections among middle-aged adults during the 2013-2014 influenza season. Proceedings of the National Academy of Sciences. 2014;111(44):15798--15803.  
**BibTeX key:** `Linderman2014AtypicalH1N1Immunity`  
**DOI/stable URL:** doi:10.1073/pnas.1409171111  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This source is useful for H1N1-specific immune-history framing. It shows how cohort-specific exposure histories can shape susceptibility and serologic interpretation in ways not captured by sequence distance alone.  
**Specific manuscript use:** Discussion  
**Claim supported:** Age-cohort immune history can create atypical susceptibility patterns for H1N1.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Wrammert2011PandemicH1N1Antibodies

**Citation:** Wrammert J, Koutsonanos D, Li GM, Edupuganti S, Sui J, Morrissey M, et al.. Broadly cross-reactive antibodies dominate the human B cell response against 2009 pandemic H1N1 influenza virus infection. Journal of Experimental Medicine. 2011;208(1):181--193.  
**BibTeX key:** `Wrammert2011PandemicH1N1Antibodies`  
**DOI/stable URL:** doi:10.1084/jem.20101352  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper is relevant because it shows that pandemic H1N1 elicited broadly cross-reactive antibodies in humans. It supports a discussion of why immune response to H1N1 can be shaped by past exposures to antigenically related strains.  
**Specific manuscript use:** Discussion  
**Claim supported:** Human responses to 2009 pandemic H1N1 can include broadly cross-reactive memory-derived antibodies.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Auladell2022InfectionHistoryVaccination

**Citation:** Auladell M, Phuong HVM, Mai LTQ, Tseng YY, Carolan L, Wilks S, et al.. Influenza virus infection history shapes antibody responses to influenza vaccination. Nature Medicine. 2022;28(2):363--372.  
**BibTeX key:** `Auladell2022InfectionHistoryVaccination`  
**DOI/stable URL:** doi:10.1038/s41591-022-01690-w  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This is one of the most directly relevant human immune-history sources. It supports the manuscript’s distinction between viral distance and observed vaccine-response phenotype.  
**Specific manuscript use:** Introduction; Discussion; Limitations  
**Claim supported:** Individual influenza infection history can substantially shape antibody responses to vaccination.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Nunez2017AgePreexistingResponses

**Citation:** Nuñez IA, Carlock MA, Allen JD, Owino SO, Moehling KK, Nowalk P, et al.. Impact of age and pre-existing influenza immune responses in humans receiving split inactivated influenza vaccine on the induction of the breadth of antibodies to influenza A strains. PLoS ONE. 2017;12(11):e0185666.  
**BibTeX key:** `Nunez2017AgePreexistingResponses`  
**DOI/stable URL:** doi:10.1371/journal.pone.0185666  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This study is directly related to the human vaccination cohort context in the manuscript. It supports the idea that age and pre-existing antibody profiles influence the breadth of measured responses.  
**Specific manuscript use:** Introduction; Discussion; Methods context  
**Claim supported:** Age and pre-existing immunity can shape the breadth of antibody induction after influenza vaccination.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Abreu2020IgARecurrentVaccination

**Citation:** Abreu RB, Clutter EF, Attari S, Sautto GA, Ross TM. IgA responses following recurrent influenza virus vaccination. Frontiers in Immunology. 2020;11:902.  
**BibTeX key:** `Abreu2020IgARecurrentVaccination`  
**DOI/stable URL:** doi:10.3389/fimmu.2020.00902  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper is useful for the repeated-vaccination and mucosal-antibody side of immune heterogeneity. It is also already related to the manuscript’s source data context.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Repeated influenza vaccination can shape antibody isotype and response patterns beyond simple HAI titers.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### McLean2014RepeatedVaccination

**Citation:** McLean HQ, Thompson MG, Sundaram ME, Meece JK, McClure DL, Friedrich TC, et al.. Impact of repeated vaccination on vaccine effectiveness against influenza A(H3N2) and B during 8 seasons. Clinical Infectious Diseases. 2014;59(10):1375--1385.  
**BibTeX key:** `McLean2014RepeatedVaccination`  
**DOI/stable URL:** doi:10.1093/cid/ciu680  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper provides epidemiologic evidence that repeated vaccination can affect vaccine effectiveness. It is useful for discussing how prior vaccination history may complicate interpretation of antibody distances.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Repeated influenza vaccination history can modify observed vaccine effectiveness across seasons.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies.

### Belongia2017RepeatedVaccinationReview

**Citation:** Belongia EA, Skowronski DM, McLean HQ, Chambers C, Sundaram ME, De Serres G. Repeated annual influenza vaccination and vaccine effectiveness: review of evidence. Expert Review of Vaccines. 2017;16(7):723--736.  
**BibTeX key:** `Belongia2017RepeatedVaccinationReview`  
**DOI/stable URL:** doi:10.1080/14760584.2017.1334554  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This review is useful for summarizing repeated-vaccination evidence without relying on a single season. It supports a limitations paragraph about unmeasured or incompletely modeled exposure history.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** The effect of repeated annual vaccination on influenza protection is heterogeneous across seasons and study designs.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies.

### Arevalo2020OASStalk

**Citation:** Arevalo P, McLean HQ, Belongia EA, Cobey S. Earliest infections predict the age distribution of seasonal influenza A cases. eLife. 2020;9:e50060.  
**BibTeX key:** `Arevalo2020OASStalk`  
**DOI/stable URL:** doi:10.7554/eLife.50060  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This paper is useful for linking early-life exposures to later population patterns of influenza disease. It is indirect for HAI distance, but supports the manuscript’s broader claim that immune history matters.  
**Specific manuscript use:** Discussion  
**Claim supported:** Early influenza exposures can leave durable signatures in later susceptibility patterns.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports/qualifies.

### Monto2017OriginalAntigenicSin

**Citation:** Monto AS, Malosh RE, Petrie JG, Martin ET. The doctrine of original antigenic sin: separating good from evil. Journal of Infectious Diseases. 2017;215(12):1782--1788.  
**BibTeX key:** `Monto2017OriginalAntigenicSin`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/28475752/  
**Category/categories:** Immune history, original antigenic sin, imprinting, and heterologous vaccine response  
**Annotation:** This review is a useful interpretive source because it cautions against treating original antigenic sin as only harmful. It helps keep the manuscript’s immune-history discussion balanced.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Original antigenic sin and imprinting can have both beneficial and limiting effects on influenza immunity.  
**Priority score:** Useful  
**Relation to manuscript claims:** Qualifies.
  **Metadata note:** DOI intentionally omitted in recreated file; verify before final manuscript use.

## HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype

### Hirst1942Hemagglutination

**Citation:** Hirst GK. The quantitative determination of influenza virus and antibodies by means of red cell agglutination. Journal of Experimental Medicine. 1942;75(1):49--64.  
**BibTeX key:** `Hirst1942Hemagglutination`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/19871167/  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This is a foundational source for hemagglutination-based influenza serology. It is historically important but methodologically dated; use it only for background on assay origins.  
**Specific manuscript use:** Background; Methods context  
**Claim supported:** Hemagglutination-based assays enabled quantitative measurement of influenza virus and antibodies.  
**Priority score:** Background only  
**Relation to manuscript claims:** Background support.

### Hobson1972HAICorrelateProtection

**Citation:** Hobson D, Curry RL, Beare AS, Ward-Gardner A. The role of serum haemagglutination-inhibiting antibody in protection against challenge infection with influenza A2 and B viruses. Journal of Hygiene. 1972;70(4):767--777.  
**BibTeX key:** `Hobson1972HAICorrelateProtection`  
**DOI/stable URL:** doi:10.1017/S0022172400022610  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This classic human challenge study is often cited for HAI as a correlate of protection. It is useful background, but should be paired with newer sources showing that HAI is incomplete.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Serum HAI antibody is associated with protection against influenza challenge, but is not a complete immune phenotype.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Couzens2014HAIAssayReliability

**Citation:** Couzens L, Gao J, Westgeest K, Sandbulte M, Lugovtsev V, Fouchier R, et al.. An optimized enzyme-linked lectin assay to measure influenza A virus neuraminidase inhibition antibody titers in human sera. Journal of Virological Methods. 2014;210:7--14.  
**BibTeX key:** `Couzens2014HAIAssayReliability`  
**DOI/stable URL:** doi:10.1016/j.jviromet.2014.09.003  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This source is included for assay context on neuraminidase inhibition rather than HAI. It helps support future-work recommendations involving NA inhibition data.  
**Specific manuscript use:** Methods context; Future work  
**Claim supported:** NA inhibition assays provide immune information not captured by HA-focused HAI assays.  
**Priority score:** Useful  
**Relation to manuscript claims:** Supports.

### Stephenson2007SerologicAssays

**Citation:** Stephenson I, Wood JM, Nicholson KG, Charlett A, Zambon MC. Detection of anti-H5 responses in human sera by HI using horse erythrocytes following MF59-adjuvanted influenza A/Duck/Singapore/97 vaccine. Vaccine. 2004;22(26):3490--3493.  
**BibTeX key:** `Stephenson2007SerologicAssays`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/15308373/  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This assay-focused source illustrates that HAI protocol details, including erythrocyte choice, can affect measured titers. It is useful for limitations about assay comparability.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** HAI measurements can depend on assay protocol choices, which affects interpretation of antigenic distance.  
**Priority score:** Useful  
**Relation to manuscript claims:** Qualifies.
  **Metadata note:** DOI intentionally omitted in recreated file; verify before final manuscript use.

### Gouma2020ChallengesVaccines

**Citation:** Gouma S, Anderson EM, Hensley SE. Challenges of making effective influenza vaccines. Annual Review of Virology. 2020;7:495--512.  
**BibTeX key:** `Gouma2020ChallengesVaccines`  
**DOI/stable URL:** doi:10.1146/annurev-virology-010320-044746  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This review synthesizes egg adaptation, immune history, and assay limitations. It is one of the best discussion citations for explaining why HA phylogeny and observed immune response can diverge.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Influenza vaccine performance is shaped by immune history, assay phenotype, and production effects beyond HA sequence relatedness.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Wu2017StructuralLowH3N2VE

**Citation:** Wu NC, Zost SJ, Thompson AJ, Oyen D, Nycholat CM, McBride R, et al.. A structural explanation for the low effectiveness of the seasonal influenza H3N2 vaccine. PLoS Pathogens. 2017;13(10):e1006682.  
**BibTeX key:** `Wu2017StructuralLowH3N2VE`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1006682  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This paper shows how egg-adaptive mutations can alter vaccine antigenicity. It supports the manuscript’s argument that measured immune responses can diverge from phylogenetic expectations for reasons not captured by wild-type HA distance.  
**Specific manuscript use:** Discussion; Limitations  
**Claim supported:** Egg-adaptive mutations can change antigenic properties and reduce H3N2 vaccine effectiveness.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Catani2024N2Landscape

**Citation:** Catani JPP, Smet A, Ysenbaert T, et al.. The antigenic landscape of human influenza N2 neuraminidases from 2009 until 2017. eLife. 2024;12:RP90782.  
**BibTeX key:** `Catani2024N2Landscape`  
**DOI/stable URL:** doi:10.7554/eLife.90782  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This is one of the most important sources for broadening the manuscript beyond HA. It shows that neuraminidase has measurable antigenic structure that can be missed by HA-only phylogenetic comparisons.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** N2 neuraminidase has its own antigenic landscape that can contribute immune-relevant variation beyond HA.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Krammer2018NActionNA

**Citation:** Krammer F, Fouchier RAM, Eichelberger MC, Webby RJ, Shaw-Saliba K, Wan H, et al.. NAction! How can neuraminidase-based immunity contribute to better influenza virus vaccines?. mBio. 2018;9(2):e02332-17.  
**BibTeX key:** `Krammer2018NActionNA`  
**DOI/stable URL:** doi:10.1128/mBio.02332-17  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This review is a strong citation for the importance of NA immunity in improved vaccines. It directly supports the manuscript’s suggested future work using neuraminidase sequence and inhibition data.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Neuraminidase immunity is a plausible contributor to protection and vaccine improvement beyond HA-head antigenicity.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Monto2015NAProtection

**Citation:** Monto AS, Petrie JG, Cross RT, Johnson E, Liu M, Zhong W, et al.. Antibody to influenza virus neuraminidase: an independent correlate of protection. Journal of Infectious Diseases. 2015;212(8):1191--1199.  
**BibTeX key:** `Monto2015NAProtection`  
**DOI/stable URL:** doi:10.1093/infdis/jiv195  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This human challenge/seroepidemiologic source supports the claim that NA antibody can contribute to protection independently of HAI antibody. It is essential for limiting HA-only interpretations of immune phenotype.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** NA antibody can be an independent correlate of protection from influenza illness.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Eichelberger2018NAFruit

**Citation:** Eichelberger MC, Monto AS. Neuraminidase, the forgotten surface antigen, emerges as an influenza vaccine target for broadened protection. Journal of Infectious Diseases. 2019;219(Supplement_1):S75--S80.  
**BibTeX key:** `Eichelberger2018NAFruit`  
**DOI/stable URL:** doi:10.1093/infdis/jiy456  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This review provides a practical bridge from NA immunology to vaccine design. It is useful for the discussion and future-work sections because the manuscript’s current metrics are HA-centered.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Neuraminidase is an underused but important influenza vaccine target.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Ekiert2012BroadNeutralizing

**Citation:** Ekiert DC, Wilson IA. Broadly neutralizing antibodies against influenza virus and prospects for universal therapies and vaccines. Current Opinion in Virology. 2012;2(2):134--141.  
**BibTeX key:** `Ekiert2012BroadNeutralizing`  
**DOI/stable URL:** doi:10.1016/j.coviro.2012.02.005  
**Category/categories:** HAI assays, neutralization, neuraminidase, and non-HA contributors to antigenic phenotype  
**Annotation:** This review helps distinguish HAI-focused HA-head immunity from broader neutralizing antibody targets, including conserved HA stalk epitopes. It supports a broader interpretation of immune phenotype than cartographic HAI distance alone.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Broadly neutralizing antibodies can target conserved influenza epitopes not fully represented by standard HAI antigenic maps.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

## Vaccine strain selection, universal influenza vaccine design, and surveillance implications

### Erbelding2018UniversalStrategicPlan

**Citation:** Erbelding EJ, Post DJ, Stemmy EJ, Roberts PC, Augustine AD, Ferguson S, et al.. A universal influenza vaccine: the strategic plan for the National Institute of Allergy and Infectious Diseases. Journal of Infectious Diseases. 2018;218(3):347--354.  
**BibTeX key:** `Erbelding2018UniversalStrategicPlan`  
**DOI/stable URL:** doi:10.1093/infdis/jiy103  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This strategic-plan paper is useful for connecting antigenic-distance work to universal-vaccine goals. It supports the introduction’s translational motivation without requiring the manuscript to claim immediate vaccine-design implications.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Universal influenza vaccine development requires better understanding of conserved and variable immune targets.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Viboud2020BeyondClinicalTrials

**Citation:** Viboud C, Gostic K, Nelson MI, Price GE, Perofsky A, Sun K, et al.. Beyond clinical trials: evolutionary and epidemiological considerations for development of a universal influenza vaccine. PLoS Pathogens. 2020;16(9):e1008583.  
**BibTeX key:** `Viboud2020BeyondClinicalTrials`  
**DOI/stable URL:** doi:10.1371/journal.ppat.1008583  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This paper is a strong source for broadening vaccine evaluation beyond narrow clinical trial endpoints. It is useful for explaining how antigenic and evolutionary distances matter for long-term universal-vaccine evaluation.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Universal influenza vaccine evaluation should account for evolutionary and epidemiologic dynamics.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports.

### Hay2018GISRSFuture

**Citation:** Hay AJ, McCauley JW. The WHO global influenza surveillance and response system (GISRS)-A future perspective. Influenza and Other Respiratory Viruses. 2018;12(5):551--557.  
**BibTeX key:** `Hay2018GISRSFuture`  
**DOI/stable URL:** doi:10.1111/irv.12565  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This review describes the surveillance system that translates genetic, antigenic, and epidemiologic data into vaccine recommendations. It is essential for explaining why metric disagreement matters operationally.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Global influenza surveillance integrates antigenic, genetic, and epidemiologic information for vaccine strain selection.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Ampofo2015VaccineVirusSelection

**Citation:** Ampofo WK, Azziz-Baumgartner E, Bashir U, Cox NJ, et al.. Strengthening the influenza vaccine virus selection and development process. Vaccine. 2015;33(36):4368--4382.  
**BibTeX key:** `Ampofo2015VaccineVirusSelection`  
**DOI/stable URL:** doi:10.1016/j.vaccine.2015.06.090  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This consultation report explains the practical constraints of vaccine virus selection, including timelines and candidate-virus development. It helps translate the manuscript’s metric comparison into public-health decision relevance.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Vaccine virus selection is constrained by surveillance evidence, antigenic characterization, and production timelines.  
**Priority score:** High  
**Relation to manuscript claims:** Supports.

### Morris2018PredictiveModeling

**Citation:** Morris DH, Gostic KM, Pompei S, Bedford T, Luksza M, Neher RA, et al.. Predictive modeling of influenza shows the promise of applied evolutionary biology. Trends in Microbiology. 2018;26(2):102--118.  
**BibTeX key:** `Morris2018PredictiveModeling`  
**DOI/stable URL:** doi:10.1016/j.tim.2017.09.004  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This review is one of the best sources for connecting evolutionary biology, antigenic data, and vaccine strain prediction. It strongly supports framing distance metrics as components of forecasting rather than stand-alone answers.  
**Specific manuscript use:** Introduction; Discussion; Future work  
**Claim supported:** Influenza prediction benefits from integrating evolutionary, antigenic, and epidemiologic information.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Huddleston2020IntegratingForecasts

**Citation:** Huddleston J, Hadfield J, Sibley TR, Lee JMS, Fay K, Ilcisin M, et al.. Integrating genotypes and phenotypes improves long-term forecasts of seasonal influenza A/H3N2 evolution. eLife. 2020;9:e60067.  
**BibTeX key:** `Huddleston2020IntegratingForecasts`  
**DOI/stable URL:** doi:10.7554/eLife.60067  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This is a key modern paper for the manuscript’s discussion because it shows the value of combining genotype and phenotype. It supports the central claim that genetic/phylogenetic similarity alone is incomplete while avoiding the overclaim that sequence has no predictive value.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Forecasting H3N2 improves when genotypic and phenotypic information are modeled jointly.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Huddleston2025TimelySelection

**Citation:** Huddleston J, Bedford T. Timely vaccine strain selection and genomic surveillance improve evolutionary forecast accuracy of seasonal influenza A/H3N2. eLife. 2025;14:RP104282.  
**BibTeX key:** `Huddleston2025TimelySelection`  
**DOI/stable URL:** doi:10.7554/eLife.104282  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This recent paper emphasizes surveillance timing and decision windows. It is useful for discussion because even a good distance metric is only useful if it is available early enough for vaccine strain selection.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Forecast utility depends on surveillance timeliness and vaccine-selection timing as well as model accuracy.  
**Priority score:** Essential  
**Relation to manuscript claims:** Supports/qualifies.

### Shi2025AIVaccineSelection

**Citation:** Shi W, Wohlwend J, Wu M, et al.. Influenza vaccine strain selection with an AI-based evolutionary and antigenicity model. Nature Medicine. 2025;31(11):3862--3870.  
**BibTeX key:** `Shi2025AIVaccineSelection`  
**DOI/stable URL:** doi:10.1038/s41591-025-03917-y  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This paper is an important recent complication to simple claims that sequence cannot guide antigenic decisions. It shows how AI-based models can combine evolutionary and antigenic information for vaccine-candidate ranking.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** AI-based models can support vaccine strain selection by jointly considering antigenicity and evolutionary success.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies/challenges.

### Cowling2024VaccineEffectivenessUniversal

**Citation:** Cowling BJ, Okoli GN. Influenza vaccine effectiveness and progress towards a universal influenza vaccine. Drugs. 2024;84(9):1013--1023.  
**BibTeX key:** `Cowling2024VaccineEffectivenessUniversal`  
**DOI/stable URL:** doi:10.1007/s40265-024-02083-8  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This recent review is useful for connecting antigenic mismatch to vaccine effectiveness and universal vaccine development. It helps distinguish antigenic distance from clinical protection.  
**Specific manuscript use:** Introduction; Discussion  
**Claim supported:** Seasonal vaccine effectiveness is influenced by antigenic match but is not reducible to a single antigenic distance.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies.

### McMillan2021NextGenerationVaccines

**Citation:** McMillan CLD, Young PR, Watterson D, Chappell KJ. The next generation of influenza vaccines: towards a universal solution. Vaccines. 2021;9(1):26.  
**BibTeX key:** `McMillan2021NextGenerationVaccines`  
**DOI/stable URL:** doi:10.3390/vaccines9010026  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This review surveys next-generation vaccine targets including HA stalk, NA, M2e, and T-cell responses. It is useful for future work because these targets are not fully captured by HA-head HAI distances.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Next-generation influenza vaccines target immune mechanisms beyond conventional strain-matched HA-head responses.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

### Nachbagauer2021ChimericHA

**Citation:** Nachbagauer R, Feser J, Naficy A, et al.. A chimeric hemagglutinin-based universal influenza virus vaccine approach induces broad and long-lasting immunity in a randomized, placebo-controlled phase I trial. Nature Medicine. 2021;27(1):106--114.  
**BibTeX key:** `Nachbagauer2021ChimericHA`  
**DOI/stable URL:** doi:10.1038/s41591-020-1118-7  
**Category/categories:** Vaccine strain selection, universal influenza vaccine design, and surveillance implications  
**Annotation:** This phase I trial provides primary human evidence for a stalk-focused universal-vaccine approach. It is useful for illustrating that future immune phenotypes may not be summarized by the same HAI-cartography distances used for seasonal HA-head matching.  
**Specific manuscript use:** Discussion; Future work  
**Claim supported:** Chimeric-HA vaccination can elicit broad stalk-focused immunity in humans.  
**Priority score:** High  
**Relation to manuscript claims:** Supports/qualifies.

## Statistical and methodological issues for distance comparison

### Mantel1967DiseaseClustering

**Citation:** Mantel N. The detection of disease clustering and a generalized regression approach. Cancer Research. 1967;27(2):209--220.  
**BibTeX key:** `Mantel1967DiseaseClustering`  
**DOI/stable URL:** https://pubmed.ncbi.nlm.nih.gov/6018555/  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This is the original Mantel-test paper. It is relevant historically, but modern critiques and alternatives should be used for the manuscript’s main statistical cautions.  
**Specific manuscript use:** Methods context  
**Claim supported:** Matrix permutation approaches can test association between distance-like structures, but require careful null-model interpretation.  
**Priority score:** Background only  
**Relation to manuscript claims:** Background/qualifies.

### Kruskal1964MDS

**Citation:** Kruskal JB. Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. Psychometrika. 1964;29:1--27.  
**BibTeX key:** `Kruskal1964MDS`  
**DOI/stable URL:** doi:10.1007/BF02289565  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This foundational MDS paper is useful background for distance-to-map methods. It supports a brief explanation of why dimensional reductions of distance matrices are approximations.  
**Specific manuscript use:** Methods context  
**Claim supported:** Multidimensional scaling represents dissimilarities in low-dimensional space by optimizing a fit criterion.  
**Priority score:** Background only  
**Relation to manuscript claims:** Methods support.

### Gower1975GeneralizedProcrustes

**Citation:** Gower JC. Generalized Procrustes analysis. Psychometrika. 1975;40:33--51.  
**BibTeX key:** `Gower1975GeneralizedProcrustes`  
**DOI/stable URL:** doi:10.1007/BF02291478  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This source is relevant for comparing ordinations or maps derived from different distances. It can support a methods or future-work suggestion to compare antigenic and genetic ordinations by Procrustes alignment.  
**Specific manuscript use:** Methods context; Future work  
**Claim supported:** Procrustes methods compare configurations after accounting for translation, rotation, and scaling.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

### Legendre2010MantelAlternatives

**Citation:** Legendre P, Fortin MJ, Borcard D. Should the Mantel test be used in spatial analysis?. Methods in Ecology and Evolution. 2015;6(11):1239--1247.  
**BibTeX key:** `Legendre2010MantelAlternatives`  
**DOI/stable URL:** doi:10.1111/2041-210X.12425  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This source is useful for discussing why Mantel-style tests are not always appropriate. It supports a limitations paragraph about matrix dependence and the need for sensitivity analyses.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Mantel tests can be inappropriate or low powered depending on the dependence structure and scientific question.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies.

### PeresNeto2001ProcrusteanMantel

**Citation:** Peres-Neto PR, Jackson DA. How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test. Oecologia. 2001;129(2):169--178.  
**BibTeX key:** `PeresNeto2001ProcrusteanMantel`  
**DOI/stable URL:** doi:10.1007/s004420100720  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This is one of the most useful methodological sources for the manuscript because it addresses comparison of multivariate structures derived from distance data. It supports using Procrustes or ordination comparisons alongside pairwise correlations.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Procrustean comparison can be more informative than Mantel correlation for some multivariate distance-derived structures.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies.

### Harmon2010PoorMantelPhylogenetic

**Citation:** Harmon LJ, Glor RE. Poor statistical performance of the Mantel test in phylogenetic comparative analyses. Evolution. 2010;64(7):2173--2178.  
**BibTeX key:** `Harmon2010PoorMantelPhylogenetic`  
**DOI/stable URL:** doi:10.1111/j.1558-5646.2010.00973.x  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This paper is highly relevant because the manuscript compares phylogenetically structured distance matrices. It supports caution against overinterpreting matrix correlations as if pairwise distances were independent observations.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Mantel tests can perform poorly in phylogenetic comparative settings.  
**Priority score:** High  
**Relation to manuscript claims:** Qualifies/challenges.

### Guillot2013DismantlingMantel

**Citation:** Guillot G, Rousset F. Dismantling the Mantel tests. Methods in Ecology and Evolution. 2013;4(4):336--344.  
**BibTeX key:** `Guillot2013DismantlingMantel`  
**DOI/stable URL:** doi:10.1111/2041-210X.12018  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This critique is useful for a statistically mature limitations section. It does not mean all distance comparisons are invalid, but it does argue against naive inferential claims from matrix correlations.  
**Specific manuscript use:** Methods context; Limitations  
**Claim supported:** Mantel-style matrix association tests can be misinterpreted and should not substitute for well-specified models.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies/challenges.

### Quilodran2025BenchmarkingMantel

**Citation:** Quilodrán CS, Currat M, Montoya-Burgos JI. Benchmarking the Mantel test and derived methods for testing association between distance matrices. Molecular Ecology Resources. 2025;25(2):e13898.  
**BibTeX key:** `Quilodran2025BenchmarkingMantel`  
**DOI/stable URL:** doi:10.1111/1755-0998.13898  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This recent benchmarking paper updates the Mantel-test debate and is useful for a balanced methods discussion. It supports the claim that matrix-comparison conclusions can be method-sensitive.  
**Specific manuscript use:** Methods context; Limitations; Future work  
**Claim supported:** Different Mantel-derived methods vary in performance, so distance-matrix association results should be interpreted cautiously.  
**Priority score:** Essential  
**Relation to manuscript claims:** Qualifies.

### Robinson1981PhylogeneticTrees

**Citation:** Robinson DF, Foulds LR. Comparison of phylogenetic trees. Mathematical Biosciences. 1981;53(1-2):131--147.  
**BibTeX key:** `Robinson1981PhylogeneticTrees`  
**DOI/stable URL:** doi:10.1016/0025-5564(81)90043-2  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This foundational paper defines the Robinson-Foulds tree-comparison framework. It is directly relevant because the manuscript compares neighbor-joining trees from distance metrics to maximum-likelihood phylogenies.  
**Specific manuscript use:** Methods context; Discussion  
**Claim supported:** Phylogenetic tree disagreement can be quantified by formal topological metrics.  
**Priority score:** High  
**Relation to manuscript claims:** Methods support.

### Smith2020GeneralizedRF

**Citation:** Smith MR. Information theoretic generalized Robinson-Foulds metrics for comparing phylogenetic trees. Bioinformatics. 2020;36(20):5007--5013.  
**BibTeX key:** `Smith2020GeneralizedRF`  
**DOI/stable URL:** doi:10.1093/bioinformatics/btaa614  
**Category/categories:** Statistical and methodological issues for distance comparison  
**Annotation:** This modern tree-space metric paper is useful for future work if the manuscript expands beyond basic RF distances. It helps justify a more nuanced tree-comparison analysis.  
**Specific manuscript use:** Methods context; Future work  
**Claim supported:** Information-theoretic tree metrics can summarize topology disagreement more richly than simple RF distance.  
**Priority score:** Useful  
**Relation to manuscript claims:** Methods support.

## Software and computational tools

### Edgar2004MUSCLE

**Citation:** Edgar RC. MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics. 2004;5:113.  
**BibTeX key:** `Edgar2004MUSCLE`  
**DOI/stable URL:** doi:10.1186/1471-2105-5-113  
**Category/categories:** Software and computational tools  
**Annotation:** This is the standard citation for MUSCLE, which the manuscript uses for multiple sequence alignment. It belongs in the methods bibliography as software/method support rather than as biological evidence.  
**Specific manuscript use:** Methods context  
**Claim supported:** MUSCLE provides multiple sequence alignments used before sequence-distance and phylogenetic analyses.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

### Bodenhofer2015MSA

**Citation:** Bodenhofer U, Bonatesta E, Horejš-Kainrath C, Hochreiter S. msa: an R package for multiple sequence alignment. Bioinformatics. 2015;31(24):3997--3999.  
**BibTeX key:** `Bodenhofer2015MSA`  
**DOI/stable URL:** doi:10.1093/bioinformatics/btv494  
**Category/categories:** Software and computational tools  
**Annotation:** This is the citation for the R package interface used for multiple sequence alignment. It should be cited in methods when documenting computational reproducibility.  
**Specific manuscript use:** Methods context  
**Claim supported:** The msa R package provides access to multiple sequence alignment algorithms in R workflows.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

### RCoreTeam2024R

**Citation:** R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing. 2024.  
**BibTeX key:** `RCoreTeam2024R`  
**DOI/stable URL:** https://www.R-project.org/  
**Category/categories:** Software and computational tools  
**Annotation:** This is the standard citation for R. It should be included with software citations rather than treated as a biological or methodological source.  
**Specific manuscript use:** Methods context  
**Claim supported:** R provides the statistical computing environment used for the analyses.  
**Priority score:** Essential  
**Relation to manuscript claims:** Methods support.

## Ranked list of the 20 most important sources to read first

1. `Smith2004AntigenicCartography` — Mapping the antigenic and genetic evolution of influenza virus.
2. `Koelle2006EpochalEvolution` — Epochal evolution shapes the phylodynamics of interpandemic influenza A (H3N2) in humans.
3. `Rambaut2008GenomicEpidemiology` — The genomic and epidemiological dynamics of human influenza A virus.
4. `Garten2009SwineOriginH1N1` — Antigenic and genetic characteristics of swine-origin 2009 A(H1N1) influenza viruses circulating in humans.
5. `Smith2009Origins2009H1N1` — Origins and evolutionary genomics of the 2009 swine-origin H1N1 influenza A epidemic.
6. `Bedford2014AntigenicDynamics` — Integrating influenza antigenic dynamics with molecular evolution.
7. `Fonville2014AntibodyLandscapes` — Antibody landscapes after influenza virus infection or vaccination.
8. `Fonville2016HumanAntiseraH3N2` — Antigenic maps of influenza A(H3N2) produced with human antisera obtained after primary infection.
9. `Bedford2015GlobalMigration` — Global circulation patterns of seasonal influenza viruses vary with antigenic drift.
10. `Gupta2006VaccineEfficacyDistance` — Quantifying influenza vaccine efficacy and antigenic distance.
11. `Anderson2018SequenceBasedCartographyH1N1` — Antigenic cartography of H1N1 influenza viruses using sequence-based antigenic distance calculation.
12. `Dang2010FLUModel` — FLU, an amino acid substitution model for influenza proteins.
13. `Gostic2016BirthYearImprinting` — Potent protection against H5N1 and H7N9 influenza via childhood hemagglutinin imprinting.
14. `Lessler2012AntigenicSeniority` — Evidence for antigenic seniority in influenza A (H3N2) antibody responses in southern China.
15. `Auladell2022InfectionHistoryVaccination` — Influenza virus infection history shapes antibody responses to influenza vaccination.
16. `Gouma2020ChallengesVaccines` — Challenges of making effective influenza vaccines.
17. `Catani2024N2Landscape` — The antigenic landscape of human influenza N2 neuraminidases from 2009 until 2017.
18. `Huddleston2020IntegratingForecasts` — Integrating genotypes and phenotypes improves long-term forecasts of seasonal influenza A/H3N2 evolution.
19. `Shah2024SeasonalAntigenicPrediction` — Seasonal antigenic prediction of influenza A H3N2 using machine learning.
20. `Quilodran2025BenchmarkingMantel` — Benchmarking the Mantel test and derived methods for testing association between distance matrices.

## Claim-to-source map

- **Influenza A evolution is shaped by immune selection and antigenic drift:** `Kim2018DriftingShifting`, `Petrova2018SeasonalEvolution`, `Ferguson2003EcologicalImmunological`, `Koelle2006EpochalEvolution`, `Smith2004AntigenicCartography`, `Bedford2015GlobalMigration`
- **H3N2 often shows more ladder-like antigenic/genetic evolution than H1N1:** `Koelle2006EpochalEvolution`, `Rambaut2008GenomicEpidemiology`, `Russell2008GlobalCirculationH3N2`, `Bedford2015GlobalMigration`, `Huddleston2020IntegratingForecasts`
- **Temporal distance can be a weak proxy when lineages re-emerge or are replaced:** `Nelson2008H1N1Reassortment`, `Rozo2015Reemergent1977H1N1`, `Garten2009SwineOriginH1N1`, `Smith2009Origins2009H1N1`, `Yin2018H1N1StackingModel`
- **HA sequence and epitope distances capture some but not all antigenic relationships:** `Gupta2006VaccineEfficacyDistance`, `Lee2004PredictingH3N2Variants`, `Anderson2018SequenceBasedCartographyH1N1`, `Harvey2016HighLowImpactSubstitutions`, `Yao2017RandomForestAntigenicity`, `Shah2024SeasonalAntigenicPrediction`
- **Antigenic cartography captures empirical immune-reactivity patterns that may not align with phylogeny:** `Smith2004AntigenicCartography`, `Bedford2014AntigenicDynamics`, `Fonville2014AntibodyLandscapes`, `Fonville2016HumanAntiseraH3N2`, `Neher2016AntigenicPhenotypes`
- **Individual immune history can reshape measured antibody responses:** `Gostic2016BirthYearImprinting`, `Gostic2019ChildhoodImprintingSeasonal`, `Lessler2012AntigenicSeniority`, `Linderman2014AtypicalH1N1Immunity`, `Auladell2022InfectionHistoryVaccination`, `Nunez2017AgePreexistingResponses`
- **HAI titers are useful but incomplete measures of immune phenotype:** `Hobson1972HAICorrelateProtection`, `Fonville2016HumanAntiseraH3N2`, `Gouma2020ChallengesVaccines`, `Hensley2009ReceptorBindingAvidity`, `Ekiert2012BroadNeutralizing`
- **NA immunity and non-HA mechanisms may contribute to disagreement between HA phylogeny and immune response:** `Catani2024N2Landscape`, `Krammer2018NActionNA`, `Monto2015NAProtection`, `Eichelberger2018NAFruit`, `McMillan2021NextGenerationVaccines`
- **Distance-matrix correlations need careful statistical interpretation:** `Mantel1967DiseaseClustering`, `PeresNeto2001ProcrusteanMantel`, `Harmon2010PoorMantelPhylogenetic`, `Guillot2013DismantlingMantel`, `Quilodran2025BenchmarkingMantel`
- **Disagreement between genetic and antigenic distances has vaccine and surveillance implications:** `Hay2018GISRSFuture`, `Ampofo2015VaccineVirusSelection`, `Morris2018PredictiveModeling`, `Huddleston2020IntegratingForecasts`, `Huddleston2025TimelySelection`, `Shi2025AIVaccineSelection`

## Likely missing arguments in the current introduction

- Clearly distinguish antigenic phenotype, genetic/phylogenetic relatedness, observed human antibody reactivity, and clinical vaccine effectiveness before comparing metrics.
- Give subtype-specific context before reporting results: H3N2 often has more coherent ladder-like replacement dynamics, whereas H1N1 has pandemic replacement, swine-origin ancestry, and re-emergence complications.
- Introduce antigenic cartography as an empirical serologic model rather than a gold standard of intrinsic viral phenotype.
- Frame sequence-based metrics as a spectrum from simple Hamming and epitope distances to phenotype-trained machine-learning models.
- Flag from the start that pairwise distance matrices are statistically dependent, so correlation estimates are descriptive unless supported by an appropriate null model.

## Likely missing arguments in the current discussion

- Discuss at least four mechanisms for discordance: serum source, immune history, assay design/receptor avidity, and non-HA immunity such as neuraminidase.
- Emphasize that disagreement between cartographic and phylogenetic trees does not imply either tree is “wrong”; they summarize different biological objects.
- Soften broad claims about temporal distance by separating H1N1, H3N2, and restricted versus long historical windows.
- Add a paragraph explaining why modern sequence-informed prediction can still be useful even when simple sequence distance is incomplete.
- Add a statistical limitations paragraph about non-independent pairwise distances, Mantel-type tests, Procrustes alternatives, and formal tree-space metrics.

## Claims that should be softened unless stronger evidence is found

- Replace “temporal methods should be avoided” with a subtype- and context-specific claim: temporal distance is especially vulnerable when lineage continuity is broken or reassortment/re-emergence occurs.
- Replace “hemagglutinin sequence is not the only factor” with a more precise statement distinguishing HA sequence, HA antigenic phenotype, non-HA immunity, and host immune history.
- Avoid implying that cartographic distance is a serum-independent viral property; it is an empirical map inferred from particular sera, viruses, assays, and modeling choices.
- Avoid implying that sequence cannot predict antigenicity; simple sequence distances are incomplete, but phenotype-trained models can recover useful signal.
- Treat Pearson correlations among pairwise distances as descriptive unless paired with appropriate resampling, permutation, or model-based sensitivity analyses.

## Suggested follow-up search strings

### PubMed
- `"influenza" AND "antigenic cartography" AND H3N2`
- `"H1N1" AND "antigenic distance" AND sequence`
- `"influenza" AND "immune imprinting" AND vaccination`
- `"neuraminidase" AND influenza AND antigenic landscape`
- `"Mantel test" AND phylogenetic distance matrices`
### Google Scholar
- `"FLU amino acid substitution model" influenza proteins`
- `"seasonal antigenic prediction" "H3N2" "machine learning"`
- `"antibody landscapes" influenza vaccination antigenic map`
- `"p-epitope" influenza vaccine effectiveness`
- `"Robinson-Foulds" antigenic phylogenetic tree comparison`
### Web of Science
- `TS=(influenza AND "antigenic drift" AND phylogeny)`
- `TS=(influenza AND "antigenic distance" AND "hemagglutinin")`
- `TS=(influenza AND "vaccine strain selection" AND forecast*)`
- `TS=(influenza AND "original antigenic sin" OR imprinting)`
- `TS=("distance matrices" AND Mantel AND phylogen*)`

## Papers that may contradict, complicate, or limit the manuscript’s interpretation

- `Shah2024SeasonalAntigenicPrediction` — shows that sequence can contain substantial learnable antigenic signal, especially for H3N2.
- `Shi2025AIVaccineSelection` — suggests sequence-informed AI can be useful for vaccine-candidate ranking despite imperfect simple-distance correlations.
- `Fonville2016HumanAntiseraH3N2` — shows antigenic maps based on human sera can reflect exposure history, complicating interpretation as a virus-only phenotype.
- `Hensley2009ReceptorBindingAvidity` — shows assay-measured antigenicity can depend on receptor-binding avidity.
- `Catani2024N2Landscape` — shows HA-only analyses omit a relevant antigenic surface protein.
- `Quilodran2025BenchmarkingMantel` — shows matrix association conclusions can be sensitive to method and null model.

## Citation deployment plan

- Opening influenza-evolution paragraph: use `Petrova2018SeasonalEvolution`, `Ferguson2003EcologicalImmunological`, `Koelle2006EpochalEvolution`, and `Bedford2015GlobalMigration`.
- H1N1 historical-complexity paragraph: use `Nelson2008H1N1Reassortment`, `Rozo2015Reemergent1977H1N1`, `Garten2009SwineOriginH1N1`, and `Smith2009Origins2009H1N1`.
- Antigenic cartography methods paragraph: use `Smith2004AntigenicCartography`, `Bedford2014AntigenicDynamics`, `Fonville2014AntibodyLandscapes`, and `Fonville2016HumanAntiseraH3N2`.
- Sequence-distance paragraph: use `Gupta2006VaccineEfficacyDistance`, `Lee2004PredictingH3N2Variants`, `Anderson2018SequenceBasedCartographyH1N1`, and `Harvey2016HighLowImpactSubstitutions`.
- Immune-history discussion: use `Gostic2016BirthYearImprinting`, `Lessler2012AntigenicSeniority`, `Auladell2022InfectionHistoryVaccination`, and `Nunez2017AgePreexistingResponses`.
- Non-HA and assay-limitations discussion: use `Gouma2020ChallengesVaccines`, `Hensley2009ReceptorBindingAvidity`, `Wu2017StructuralLowH3N2VE`, `Catani2024N2Landscape`, and `Monto2015NAProtection`.
- Vaccine-selection implications: use `Hay2018GISRSFuture`, `Morris2018PredictiveModeling`, `Huddleston2020IntegratingForecasts`, `Huddleston2025TimelySelection`, and `Shi2025AIVaccineSelection`.
- Statistical limitations: use `PeresNeto2001ProcrusteanMantel`, `Harmon2010PoorMantelPhylogenetic`, `Guillot2013DismantlingMantel`, `Quilodran2025BenchmarkingMantel`, and `Robinson1981PhylogeneticTrees`.
