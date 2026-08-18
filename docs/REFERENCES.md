# References

This file records external sources that materially informed liftAssess's scientific motivation,
evidence model, provenance rules, validation strategy, provider/resource handling, or documented
real-world problem cases.

References are formatted in a Vancouver-style numbered bibliography. For living web resources,
the access date is included because provider documentation, directory contents, and online support
threads can change over time.

This audit covers sources recoverable from the liftAssess repository/code history and prior
liftAssess project research discussions through **2026-08-17**. Sources mentioned only as possible
future additions (for example, an eventual Ensembl Compara evidence source) are intentionally not
listed unless they materially informed a current design or implementation decision.

## Peer-reviewed literature

1. Luu PL, Ong PT, Dinh TP, Clark SJ. Benchmark study comparing liftover tools for genome
   conversion of epigenome sequencing data. *NAR Genom Bioinform.* 2020;2(3):lqaa054.
   doi:10.1093/nargab/lqaa054. Available from:
   https://doi.org/10.1093/nargab/lqaa054

2. Genovese G, Rockweiler NB, Gorman BR, Bigdeli TB, Pato MT, Pato CN, et al.
   BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across
   genome assemblies. *Bioinformatics.* 2024;40(2):btae038.
   doi:10.1093/bioinformatics/btae038. Available from:
   https://doi.org/10.1093/bioinformatics/btae038

3. Shumate A, Salzberg SL. Liftoff: accurate mapping of gene annotations.
   *Bioinformatics.* 2021;37(12):1639-1643. doi:10.1093/bioinformatics/btaa1016.
   Available from: https://doi.org/10.1093/bioinformatics/btaa1016

4. Jagannathan V, Hitte C, Kidd JM, Masterson P, Murphy TD, Emery S, et al.
   Dog10K_Boxer_Tasha_1.0: A Long-Read Assembly of the Dog Reference Genome.
   *Genes (Basel).* 2021;12(6):847. doi:10.3390/genes12060847.
   Available from: https://doi.org/10.3390/genes12060847

5. Wang C, Wallerman O, Arendt ML, Sundström E, Karlsson Å, Nordin J, et al.
   A novel canine reference genome resolves genomic architecture and uncovers transcript
   complexity. *Commun Biol.* 2021;4:185. doi:10.1038/s42003-021-01698-x.
   Available from: https://doi.org/10.1038/s42003-021-01698-x

## Provider, format, standards, and implementation documentation

6. UCSC Genome Browser Project. Chain Format [Internet]. Santa Cruz (CA): University of
   California, Santa Cruz; [cited 2026 Aug 17]. Available from:
   https://genome.ucsc.edu/goldenPath/help/chain.html

7. UCSC Genome Browser Project. Net Format [Internet]. Santa Cruz (CA): University of
   California, Santa Cruz; [cited 2026 Aug 17]. Available from:
   https://genome.ucsc.edu/goldenPath/help/net.html

8. UCSC Genome Browser Project. Genome Browser User's Guide [Internet]. Santa Cruz (CA):
   University of California, Santa Cruz; [cited 2026 Aug 17]. Available from:
   https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html

9. UCSC Genome Browser Project. Genome Browser Licensing [Internet]. Santa Cruz (CA):
   University of California, Santa Cruz; [cited 2026 Aug 17]. Available from:
   https://genome.ucsc.edu/license/

10. UCSC Genome Browser Project. canFam3 liftOver download directory and README/terms
    [Internet]. Santa Cruz (CA): University of California, Santa Cruz;
    [cited 2026 Aug 17]. Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/

11. UCSC Genome Browser Project. canFam3 versus canFam4 comparative alignment directory and
    README [Internet]. Santa Cruz (CA): University of California, Santa Cruz;
    [cited 2026 Aug 17]. Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/

12. UCSC Genome Browser Project. canFam3 versus canFam4 provider MD5 metadata
    [Internet]. Santa Cruz (CA): University of California, Santa Cruz;
    [cited 2026 Aug 17]. Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt

13. UCSC Genome Browser Project. canFam4 versus canFam3 reciprocal-best alignment directory
    [Internet]. Santa Cruz (CA): University of California, Santa Cruz;
    [cited 2026 Aug 17]. Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/

14. UCSC Genome Browser Project. canFam4 versus canFam3 reciprocal-best provider MD5 metadata
    [Internet]. Santa Cruz (CA): University of California, Santa Cruz;
    [cited 2026 Aug 17]. Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/md5sum.txt

15. UCSC Genome Browser Project. canFam6 liftOver download directory [Internet].
    Santa Cruz (CA): University of California, Santa Cruz; [cited 2026 Aug 17].
    Available from:
    https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/liftOver/

16. UCSC Genome Browser Project. kent source tree license [Internet]. GitHub;
    [cited 2026 Aug 17]. Available from:
    https://github.com/ucscGenomeBrowser/kent/blob/master/LICENSE

17. UCSC Genome Browser Project. doBlastzChainNet.pl [source code on Internet]. GitHub;
    [cited 2026 Aug 17]. Available from:
    https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doBlastzChainNet.pl

18. UCSC Genome Browser Project. doRecipBest.pl [source code on Internet]. GitHub;
    [cited 2026 Aug 17]. Available from:
    https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doRecipBest.pl

19. UCSC Genome Browser Project. Genome Browser Data Downloads [Internet]. Santa Cruz (CA):
    University of California, Santa Cruz; [cited 2026 Aug 17]. Available from:
    https://genome.ucsc.edu/goldenPath/help/ftp.html

20. Fielding R, Nottingham M, Reschke J. HTTP Semantics. RFC 9110 [Internet].
    Internet Engineering Task Force; 2022 Jun [cited 2026 Aug 17].
    Available from: https://www.rfc-editor.org/rfc/rfc9110

21. National Center for Biotechnology Information. Genome sequence report
    [Internet]. Bethesda (MD): National Library of Medicine (US);
    [cited 2026 Aug 17]. Available from:
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/

22. National Center for Biotechnology Information. Create a table from the genome data reports:
    genome sequence report fields [Internet]. Bethesda (MD): National Library of Medicine (US);
    [cited 2026 Aug 17]. Available from:
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/using-dataformat/genome-data-reports/

23. Global Alliance for Genomics and Health. Refget Sequences v2.0.0
    [Internet]. GA4GH; [cited 2026 Aug 17]. Available from:
    https://ga4gh.github.io/refget/sequences/

24. Global Alliance for Genomics and Health. Refget Sequence Collections v1.0.0
    [Internet]. GA4GH; [cited 2026 Aug 17]. Available from:
    https://ga4gh.github.io/refget/seqcols/

### Validation-data note

The real `canFam3` to `canFam4` mechanical fixture consumes five exact UCSC comparative
artifacts: the all-chain, ordinary net, syntenic net, reciprocal-best chain, and reciprocal-best
net. The provider directories and checksum manifests above are the bibliographic/provider
references for those artifacts. Exact per-file URLs, byte sizes, provider checksums, liftAssess
SHA-256 identities, retrieval metadata, and consumption status belong in the verifier/cache/report
provenance rather than being duplicated as static bibliography entries here.

## Support reports, issues, and other technical evidence

The sources in this section are evidence that researchers encounter recurring liftover ambiguity
and interpretation problems. They are **not** treated as biological ground truth or as validation
of a particular liftAssess verdict.

25. Biostars. Which locus should be considered the true ortholog after liftOver disagreement
    between chrUn and chr16 in canFam4 miRNA mapping? [Internet].
    [cited 2026 Aug 17]. Available from:
    https://www.biostars.org/p/9619818/

26. Biostars. UCSC liftover [Internet]. 2021 Jul 15 [cited 2026 Aug 17].
    Available from: https://www.biostars.org/p/9480565/

27. Bioconductor Support. Bug in rtracklayer liftOver function [Internet].
    [cited 2026 Aug 17]. Available from:
    https://support.bioconductor.org/p/9136646/

28. Bioconductor Support. rtracklayer::lifOver one-to-many regions and not match UCSC
    [Internet]. [cited 2026 Aug 17]. Available from:
    https://support.bioconductor.org/p/99306/

29. hbc/giab_remap_38. Crossmap liftover issues. Issue #1 [Internet]. GitHub;
    2016 Apr 20 [cited 2026 Aug 17]. Available from:
    https://github.com/hbc/giab_remap_38/issues/1
