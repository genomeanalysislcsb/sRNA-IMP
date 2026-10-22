"""Reference library for srna_methods.py: key -> (authors, year, title, journal, volume, pages, doi).
Authors as 'Surname, I.I.' strings; software without a paper: journal = '' and doi = URL."""

REFS = {
    # ---- read processing ----
    "cutadapt": (["Martin, M."], 2011, "Cutadapt removes adapter sequences from high-throughput sequencing reads",
                 "EMBnet.journal", "17", "10-12", "10.14806/ej.17.1.200"),
    "fastqc": (["Andrews, S."], 2010, "FastQC: a quality control tool for high throughput sequence data",
               "", "", "", "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"),
    "multiqc": (["Ewels, P.", "Magnusson, M.", "Lundin, S.", "Käller, M."], 2016,
                "MultiQC: summarize analysis results for multiple tools and samples in a single report",
                "Bioinformatics", "32", "3047-3048", "10.1093/bioinformatics/btw354"),
    "bowtie": (["Langmead, B.", "Trapnell, C.", "Pop, M.", "Salzberg, S.L."], 2009,
               "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome",
               "Genome Biol.", "10", "R25", "10.1186/gb-2009-10-3-r25"),
    "samtools": (["Danecek, P.", "Bonfield, J.K.", "Liddle, J.", "Marshall, J.", "Ohan, V.", "Pollard, M.O.",
                  "Whitwham, A.", "Keane, T.", "McCarthy, S.A.", "Davies, R.M.", "Li, H."], 2021,
                 "Twelve years of SAMtools and BCFtools", "GigaScience", "10", "giab008",
                 "10.1093/gigascience/giab008"),
    "bedtools": (["Quinlan, A.R.", "Hall, I.M."], 2010, "BEDTools: a flexible suite of utilities for comparing "
                 "genomic features", "Bioinformatics", "26", "841-842", "10.1093/bioinformatics/btq033"),
    "bedops": (["Neph, S.", "Kuehn, M.S.", "Reynolds, A.P.", "Haugen, E.", "Thurman, R.E.", "Johnson, A.K.",
                "Rynes, E.", "Maurano, M.T.", "Vierstra, J.", "Thomas, S.", "et al."], 2012,
               "BEDOPS: high-performance genomic feature operations", "Bioinformatics", "28", "1919-1920",
               "10.1093/bioinformatics/bts277"),
    "seqkit": (["Shen, W.", "Sipos, B.", "Zhao, L."], 2024, "SeqKit2: A Swiss army knife for sequence and "
               "alignment processing", "iMeta", "3", "e191", "10.1002/imt2.191"),
    "ribodetector": (["Deng, Z.-L.", "Münch, P.C.", "Mreches, R.", "McHardy, A.C."], 2022,
                     "Rapid and accurate identification of ribosomal RNA sequences via deep learning",
                     "Nucleic Acids Res.", "50", "e60", "10.1093/nar/gkac112"),
    "kraken2": (["Wood, D.E.", "Lu, J.", "Langmead, B."], 2019, "Improved metagenomic analysis with Kraken 2",
                "Genome Biol.", "20", "257", "10.1186/s13059-019-1891-0"),
    "bracken": (["Lu, J.", "Breitwieser, F.P.", "Thielen, P.", "Salzberg, S.L."], 2017,
                "Bracken: estimating species abundance in metagenomics data", "PeerJ Comput. Sci.", "3", "e104",
                "10.7717/peerj-cs.104"),
    "gtdb": (["Parks, D.H.", "Chuvochina, M.", "Rinke, C.", "Mussig, A.J.", "Chaumeil, P.-A.", "Hugenholtz, P."],
             2022, "GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically "
             "consistent, rank normalized and complete genome-based taxonomy", "Nucleic Acids Res.", "50",
             "D785-D794", "10.1093/nar/gkab776"),
    # ---- references / databases ----
    "t2t": (["Nurk, S.", "Koren, S.", "Rhie, A.", "Rautiainen, M.", "Bzikadze, A.V.", "Mikheenko, A.",
             "Vollger, M.R.", "Altemose, N.", "Uralsky, L.", "Gershman, A.", "et al."], 2022,
            "The complete sequence of a human genome", "Science", "376", "44-53", "10.1126/science.abj6987"),
    "mirbase": (["Kozomara, A.", "Birgaoanu, M.", "Griffiths-Jones, S."], 2019,
                "miRBase: from microRNA sequences to function", "Nucleic Acids Res.", "47", "D155-D162",
                "10.1093/nar/gky1141"),
    "gtrnadb": (["Chan, P.P.", "Lowe, T.M."], 2016, "GtRNAdb 2.0: an expanded database of transfer RNA genes "
                "identified in complete and draft genomes", "Nucleic Acids Res.", "44", "D184-D189",
                "10.1093/nar/gkv1309"),
    "rfam": (["Kalvari, I.", "Nawrocki, E.P.", "Ontiveros-Palacios, N.", "Argasinska, J.", "Lamkiewicz, K.",
              "Marz, M.", "Griffiths-Jones, S.", "Toffano-Nioche, C.", "Gautheret, D.", "Weinberg, Z.", "et al."],
             2021, "Rfam 14: expanded coverage of metagenomic, viral and microRNA families",
             "Nucleic Acids Res.", "49", "D192-D200", "10.1093/nar/gkaa1047"),
    # ---- multi-omics context ----
    "imp": (["Narayanasamy, S.", "Jarosz, Y.", "Muller, E.E.L.", "Heintz-Buschart, A.", "Herold, M.",
             "Kaysen, A.", "Laczny, C.C.", "Pinel, N.", "May, P.", "Wilmes, P."], 2016,
            "IMP: a pipeline for reproducible reference-independent integrated metagenomic and "
            "metatranscriptomic analyses", "Genome Biol.", "17", "260", "10.1186/s13059-016-1116-8"),
    "bakta": (["Schwengers, O.", "Jelonek, L.", "Dieckmann, M.A.", "Beyvers, S.", "Blom, J.", "Goesmann, A."],
              2021, "Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence "
              "identification", "Microb. Genom.", "7", "000685", "10.1099/mgen.0.000685"),
    "mantis": (["Queirós, P.", "Delogu, F.", "Hickl, O.", "May, P.", "Wilmes, P."], 2021,
               "Mantis: flexible and consensus-driven genome annotation", "GigaScience", "10", "giab042",
               "10.1093/gigascience/giab042"),
    "kegg": (["Kanehisa, M.", "Goto, S."], 2000, "KEGG: Kyoto Encyclopedia of Genes and Genomes",
             "Nucleic Acids Res.", "28", "27-30", "10.1093/nar/28.1.27"),
    "pfam": (["Mistry, J.", "Chuguransky, S.", "Williams, L.", "Qureshi, M.", "Salazar, G.A.",
              "Sonnhammer, E.L.L.", "Tosatto, S.C.E.", "Paladin, L.", "Raj, S.", "Richardson, L.J.", "et al."],
             2021, "Pfam: The protein families database in 2021", "Nucleic Acids Res.", "49", "D412-D419",
             "10.1093/nar/gkaa913"),
    "cog": (["Galperin, M.Y.", "Wolf, Y.I.", "Makarova, K.S.", "Vera Alvarez, R.", "Landsman, D.",
             "Koonin, E.V."], 2021, "COG database update: focus on microbial diversity, model organisms, and "
             "widespread pathogens", "Nucleic Acids Res.", "49", "D274-D281", "10.1093/nar/gkaa1018"),
    "metawrap": (["Uritskiy, G.V.", "DiRuggiero, J.", "Taylor, J."], 2018, "MetaWRAP—a flexible pipeline for "
                 "genome-resolved metagenomic data analysis", "Microbiome", "6", "158",
                 "10.1186/s40168-018-0541-1"),
    "gtdbtk": (["Chaumeil, P.-A.", "Mussig, A.J.", "Hugenholtz, P.", "Parks, D.H."], 2022,
               "GTDB-Tk v2: memory friendly classification with the genome taxonomy database",
               "Bioinformatics", "38", "5315-5316", "10.1093/bioinformatics/btac672"),
    # ---- ncRNA discovery ----
    "infernal": (["Nawrocki, E.P.", "Eddy, S.R."], 2013, "Infernal 1.1: 100-fold faster RNA homology searches",
                 "Bioinformatics", "29", "2933-2935", "10.1093/bioinformatics/btt509"),
    "aragorn": (["Laslett, D.", "Canback, B."], 2004, "ARAGORN, a program to detect tRNA genes and tmRNA genes "
                "in nucleotide sequences", "Nucleic Acids Res.", "32", "11-16", "10.1093/nar/gkh152"),
    "trnascan": (["Chan, P.P.", "Lin, B.Y.", "Mak, A.J.", "Lowe, T.M."], 2021, "tRNAscan-SE 2.0: improved "
                 "detection and functional classification of transfer RNA genes", "Nucleic Acids Res.", "49",
                 "9077-9096", "10.1093/nar/gkab688"),
    "viennarna": (["Lorenz, R.", "Bernhart, S.H.", "Höner zu Siederdissen, C.", "Tafer, H.", "Flamm, C.",
                   "Stadler, P.F.", "Hofacker, I.L."], 2011, "ViennaRNA Package 2.0",
                  "Algorithms Mol. Biol.", "6", "26", "10.1186/1748-7188-6-26"),
    "trf": (["Kumar, P.", "Kuscu, C.", "Dutta, A."], 2016, "Biogenesis and function of transfer RNA-related "
            "fragments (tRFs)", "Trends Biochem. Sci.", "41", "679-689", "10.1016/j.tibs.2016.05.004"),
    "rin": (["Schroeder, A.", "Mueller, O.", "Stocker, S.", "Salowsky, R.", "Leiber, M.", "Gassmann, M.",
             "Lightfoot, S.", "Menzel, W.", "Granzow, M.", "Ragg, T."], 2006, "The RIN: an RNA integrity number "
            "for assigning integrity values to RNA measurements", "BMC Mol. Biol.", "7", "3",
            "10.1186/1471-2199-7-3"),
    # ---- software environment ----
    "numpy": (["Harris, C.R.", "Millman, K.J.", "van der Walt, S.J.", "Gommers, R.", "Virtanen, P.",
               "Cournapeau, D.", "Wieser, E.", "Taylor, J.", "Berg, S.", "Smith, N.J.", "et al."], 2020,
              "Array programming with NumPy", "Nature", "585", "357-362", "10.1038/s41586-020-2649-2"),
    "matplotlib": (["Hunter, J.D."], 2007, "Matplotlib: a 2D graphics environment", "Comput. Sci. Eng.", "9",
                   "90-95", "10.1109/MCSE.2007.55"),
    # ---- statistics ----
    "kruskal": (["Kruskal, W.H.", "Wallis, W.A."], 1952, "Use of ranks in one-criterion variance analysis",
                "J. Am. Stat. Assoc.", "47", "583-621", "10.1080/01621459.1952.10483441"),
    "mannwhitney": (["Mann, H.B.", "Whitney, D.R."], 1947, "On a test of whether one of two random variables is "
                    "stochastically larger than the other", "Ann. Math. Stat.", "18", "50-60",
                    "10.1214/aoms/1177730491"),
    "wilcoxon": (["Wilcoxon, F."], 1945, "Individual comparisons by ranking methods", "Biometrics Bull.", "1",
                 "80-83", "10.2307/3001968"),
    "spearman": (["Spearman, C."], 1904, "The proof and measurement of association between two things",
                 "Am. J. Psychol.", "15", "72-101", "10.2307/1412159"),
    "bh": (["Benjamini, Y.", "Hochberg, Y."], 1995, "Controlling the false discovery rate: a practical and "
           "powerful approach to multiple testing", "J. R. Stat. Soc. B", "57", "289-300",
           "10.1111/j.2517-6161.1995.tb02031.x"),
    "braycurtis": (["Bray, J.R.", "Curtis, J.T."], 1957, "An ordination of the upland forest communities of "
                   "southern Wisconsin", "Ecol. Monogr.", "27", "325-349", "10.2307/1942268"),
    "permanova": (["Anderson, M.J."], 2001, "A new method for non-parametric multivariate analysis of variance",
                  "Austral Ecol.", "26", "32-46", "10.1111/j.1442-9993.2001.01070.pp.x"),
    "pcoa": (["Gower, J.C."], 1966, "Some distance properties of latent root and vector methods used in "
             "multivariate analysis", "Biometrika", "53", "325-338", "10.1093/biomet/53.3-4.325"),
    "nmds": (["Kruskal, J.B."], 1964, "Multidimensional scaling by optimizing goodness of fit to a nonmetric "
             "hypothesis", "Psychometrika", "29", "1-27", "10.1007/BF02289565"),
    "smacof": (["de Leeuw, J.", "Mair, P."], 2009, "Multidimensional scaling using majorization: SMACOF in R",
               "J. Stat. Softw.", "31", "1-30", "10.18637/jss.v031.i03"),
    "clr": (["Aitchison, J."], 1982, "The statistical analysis of compositional data", "J. R. Stat. Soc. B",
            "44", "139-177", "10.1111/j.2517-6161.1982.tb01195.x"),
    "gloor": (["Gloor, G.B.", "Macklaim, J.M.", "Pawlowsky-Glahn, V.", "Egozcue, J.J."], 2017,
              "Microbiome datasets are compositional: and this is not optional", "Front. Microbiol.", "8", "2224",
              "10.3389/fmicb.2017.02224"),
    "maaslin2": (["Mallick, H.", "Rahnavard, A.", "McIver, L.J.", "Ma, S.", "Zhang, Y.", "Nguyen, L.H.",
                  "Tickle, T.L.", "Weingart, G.", "Ren, B.", "Schwager, E.H.", "et al."], 2021,
                 "Multivariable association discovery in population-scale meta-omics studies",
                 "PLoS Comput. Biol.", "17", "e1009442", "10.1371/journal.pcbi.1009442"),
    "hc3": (["MacKinnon, J.G.", "White, H."], 1985, "Some heteroskedasticity-consistent covariance matrix "
            "estimators with improved finite sample properties", "J. Econom.", "29", "305-325",
            "10.1016/0304-4076(85)90158-7"),
    "cluster": (["Cameron, A.C.", "Miller, D.L."], 2015, "A practitioner's guide to cluster-robust inference",
                "J. Hum. Resour.", "50", "317-372", "10.3368/jhr.50.2.317"),
    "vanelteren": (["van Elteren, P.H."], 1960, "On the combination of independent two-sample tests of "
                   "Wilcoxon", "Bull. Int. Stat. Inst.", "37", "351-361", ""),
    "iglewicz": (["Iglewicz, B.", "Hoaglin, D.C."], 1993, "How to Detect and Handle Outliers",
                 "ASQC Quality Press", "", "", ""),
}
