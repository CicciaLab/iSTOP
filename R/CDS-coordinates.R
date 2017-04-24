# ---- CDS from UCSC ----
#' Generate a CDS coordinates table
#'
#' `CDS` takes transcript annotation tables in [UCSC format](https://genome.ucsc.edu/cgi-bin/hgTables)
#' and reshapes them to have coordinates for each exon represented on a single
#' row, rather than collapsed into a comma separated string in a single cell.
#'
#' @param tx A URL to a genome's transcript reference file. This table must have
#' tab separated fields and contain, identifiers for each transcript,
#' chromosome, strand, CDS start/end, and exon start/end information.
#' There should be only one row per transcript and exon start/end columns should
#' contain comma separated cooordinates for each exon in the transcript.
#' An example file can be found
#' [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz).
#'
#' @param gene A URL to a tab separated file that maps transcript identifiers
#' to common gene names. An example file can be found
#' [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/kgXref.txt.gz).
#'
#' @param tx_cols A character vector of expected column names for the
#' known-gene reference file. Required columns: "tx", "chr", "strand",
#' "cds_start", "cds_end", "exon_start" and "exon_end". All other columns will
#' be ignored.
#'
#' @param gene_cols A character vector of expected column names for the
#' cross-reference file. Required columns: "tx" and "gene". All other columns
#' will be ignored.
#'
#' @param shift_start Number of bases to shift the start positions. Defaults
#' to 1 as this is necessary for compatibility with [Biostrings::getSeq] which
#' includes the start position in the returned sequence and begins counting bases
#' at 1.
#'
#' @param shift_end Number of bases to shift the end positions. Defaults to 0.
#'
#' @details The output of `CDS` should meet the following standards (1) each row should
#' represent the coordinates of a single exon, (2) exons should be numbered in
#' order with reference to the transcript's strand (e.g. the first exon should
#' include the start codon). The absolute numbering is unimportant so long as they
#' are numbered in the correct order. (3) The first and last exon coordinates should
#' begin with the start codon and end with the stop codon.
#'
#' To save the trouble of looking up URLs, pre-defined CDS builders
#' are provided. They are named `CDS_<Species>_<data-source>_<genome-assembly-ID>()`.
#'
#' @return A data.frame with the following columns where each row represents
#' a single exon:
#'
#'   - __COLUMN-NAME__ _`DATA-TYPE`_ DESCRIPTION
#'   - __tx__      _`chr`_ Transcript symbol
#'   - __gene__    _`chr`_ Gene symbol
#'   - __exon__    _`int`_ Exon rank in gene (lowest contains ATG, highest contains native Stop)
#'   - __chr__     _`chr`_ Chromosome
#'   - __strand__  _`chr`_ Strand (+/-)
#'   - __start__   _`int`_ CDS coordinate start (always <= __end__)
#'   - __end__     _`int`_ CDS coordinate end (always >= __start__)
#'
#' @rdname CDS
#' @importFrom tidyr separate_rows
#' @export
#' @md

CDS <- function(tx, gene, tx_cols, gene_cols, shift_start = 1L, shift_end = 0L) {

  Validate_CDS_args(tx, gene, tx_cols, gene_cols, shift_start, shift_end)

  message('Downloading CDS coordinates from:\n    ', tx, '\n    ', gene)

  # Download known gene dataset from UCSC
  txs <-
    tx %>%
    read_tsv(
      col_names = tx_cols,
      col_types = cols_only(  # keep only the following columns
        tx         = col_character(),
        chr        = col_character(),
        strand     = col_character(),
        cds_start  = col_integer(),
        cds_end    = col_integer(),
        exon_start = col_character(),
        exon_end   = col_character()))

  # Download known gene cross reference dataset from UCSC
  genes <-
    gene %>%
    read_tsv(
      col_names = gene_cols,
      col_types = cols_only(  # keep only the following columns
        tx   = col_character(),
        gene = col_character()))

  # Expand coordinates and keep CDS coordinates only
  expanded <-
    txs %>%
    separate_rows(exon_start, exon_end, sep = ',', convert = T) %>%
    filter(
      !is.na(exon_start),    # trailing commas result in NA after separation
      cds_start != cds_end   # remove non-coding transcripts
    ) %>%
    group_by(tx) %>%
    mutate(
      exon = switch(unique(strand), '+' = 1:n(), '-' = n():1),
      # Identify first and last exon based on location of CDS start and end
      first = (cds_start >= exon_start & cds_start <= exon_end),
      last  = (cds_end   >= exon_start & cds_end   <= exon_end),
      start = if_else(first, cds_start, exon_start) + shift_start,
      end   = if_else(last,  cds_end,   exon_end)   + shift_end,
      # Mark first and last exons of CDS
      max_exon = max(exon[first | last]),
      min_exon = min(exon[first | last])
    ) %>%
    # Remove exons that are not part of CDS
    filter(exon >= min_exon, exon <= max_exon) %>%
    select(tx, exon, chr, strand, start, end)

  # Add gene names to table
  right_join(genes, expanded, by = 'tx') %>% arrange(tx, exon)
}

# ---- Pre-defined CDS constructors ----
#' @rdname CDS
#' @export
CDS_example <- function() {
  read_csv(system.file('db/CDS-ATM-ATR-example.csv', package = 'iSTOP'), col_types = cols())
}

#' @rdname CDS
#' @export
CDS_Celegans_UCSC_ce11 <- function() {
  message('Downloading CDS coordinates for all transcripts in Celegans from UCSC ce11...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/ensGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/ensemblToGeneName.txt.gz',
      tx_cols = c(
        'bin', 'tx', 'chr', 'strand', 'txStart', 'txEnd', 'cds_start', 'cds_end',
        'exonCount', 'exon_start', 'exon_end', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'),
      gene_cols = c('tx', 'gene'))
}

#' @rdname CDS
#' @export
CDS_Dmelanogaster_UCSC_dm6 <- function() {
  message('Gathering CDS coordinates for all transcripts in Dmelanogaster from UCSC dm6...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/ensGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/ensemblToGeneName.txt.gz',
      tx_cols = c(
        'bin', 'tx', 'chr', 'strand', 'txStart', 'txEnd', 'cds_start', 'cds_end',
        'exonCount', 'exon_start', 'exon_end', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'),
      gene_cols = c('tx', 'gene'))
}

#' @rdname CDS
#' @export
CDS_Drerio_UCSC_danRer10 <- function() {
  message('Gathering CDS coordinates for all transcripts in Drerio from UCSC danRer10...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/danRer10/database/ensGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/danRer10/database/ensemblToGeneName.txt.gz',
      tx_cols = c(
        'bin', 'tx', 'chr', 'strand', 'txStart', 'txEnd', 'cds_start', 'cds_end',
        'exonCount', 'exon_start', 'exon_end', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'),
      gene_cols = c('tx', 'gene'))
}

#' @rdname CDS
#' @export
CDS_Hsapiens_UCSC_hg38 <- function() {
  message('Gathering CDS coordinates for all transcripts in Hsapiens from UCSC hg38...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/kgXref.txt.gz',
      tx_cols = c(
        'tx', 'chr', 'strand', 'tx_start', 'tx_end', 'cds_start', 'cds_end',
        'exon_count', 'exon_start', 'exon_end', 'protein_id', 'align_id'),
      gene_cols = c(
        'tx', 'mRNA', 'spID', 'spDisplayID', 'gene', 'refseq', 'protAcc',
        'description', 'rfamAcc', 'tRnaName'))
}

#' @rdname CDS
#' @export
CDS_Mmusculus_UCSC_mm10 <- function() {
  message('Gathering CDS coordinates for all transcripts in Mmusculus from UCSC mm10...')
  CDS(tx   = 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz',
      gene = 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/kgXref.txt.gz',
      tx_cols = c(
        'tx', 'chr', 'strand', 'tx_start', 'tx_end', 'cds_start', 'cds_end',
        'exon_count', 'exon_start', 'exon_end', 'protein_id', 'align_id'),
      gene_cols = c(
        'tx', 'mRNA', 'spID', 'spDisplayID', 'gene', 'refseq', 'protAcc',
        'description', 'rfamAcc', 'tRnaName'))
}

#' @rdname CDS
#' @export
CDS_Rnorvegicus_UCSC_rn6 <- function() {
  message('Gathering CDS coordinates for all transcripts in Rnorvegicus from UCSC rn6...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/rn6/database/ensGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/rn6/database/ensemblToGeneName.txt.gz',
      tx_cols = c(
        'bin', 'tx', 'chr', 'strand', 'txStart', 'txEnd', 'cds_start', 'cds_end',
        'exonCount', 'exon_start', 'exon_end', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'),
      gene_cols = c('tx', 'gene'))
}

#' @rdname CDS
#' @export
CDS_Scerevisiae_UCSC_sacCer3 <- function() {
  message('Gathering CDS coordinates for all transcripts in Scerevisiae UCSC sacCer3...')
  CDS(tx      = 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/sgdGene.txt.gz',
      gene    = 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/sgdToName.txt.gz',
      tx_cols = c(
        'bin', 'tx', 'chr', 'strand', 'txStart', 'txEnd', 'cds_start', 'cds_end',
        'exonCount', 'exon_start', 'exon_end', 'proteinID'),
      gene_cols = c('tx', 'gene'))
}

#' @rdname CDS
#' @export
CDS_Athaliana_BioMart_plantsmart28 <- function() {
  require_bioc_install(
    require   = c('TxDb.Athaliana.BioMart.plantsmart28', 'org.At.tair.db'),
    recommend = 'BSgenome.Athaliana.TAIR.TAIR9'
  )

  message('Gathering CDS coordinates for all transcripts in Athaliana BioMart TAIR9 plantsmart28...')

  CDS_from_Bioconductor(
    txdb    = TxDb.Athaliana.BioMart.plantsmart28::TxDb.Athaliana.BioMart.plantsmart28,
    orgdb   = org.At.tair.db::org.At.tair.db,
    gene_id = 'TAIR',
    fix_chr = function(.) { # Chromosome names in TxDb do not match BSGenome
      c('1'  = 'Chr1',
        '2'  = 'Chr2',
        '3'  = 'Chr3',
        '4'  = 'Chr4',
        '5'  = 'Chr5',
        'Mt' = 'ChrM',
        'Pt' = 'ChrC')[.]
    }
  )
}


# ---- CDS coordinates from Bioconductor ----
# CDS coordinates from Bioconductor packages
#
# Compiles a data frame of coding sequence coordinates for all exons in a
# [TxDb][GenomicFeatures::TxDb] object. Common gene names are also included
# from an [OrgDb][AnnotationDbi::AnnotationDb] object. Note that CDS coordinates are
# nearly equivalent to exon coordinates except the first and last exons that
# contain coding sequence should have coordinates that begin with the
# translation start codon and end with the stop codon.
#
# @param txdb    [TxDb][GenomicFeatures::TxDb]
#                Transcript database object.
# @param orgdb   [OrgDb][AnnotationDbi::AnnotationDb]
#                Organism annotation database object.
# @param gene_id A string.
#                Name of gene ID column in `orgdb`. See [columns][AnnotationDbi::columns].
# @param gene    A string.
#                Name of gene symbol column in `orgdb`. See [columns][AnnotationDbi::columns]
# @param fix_chr **Expert use only** A function to fix the `chr` column.
#                This function should accept a vector and return a
#                vector of the same length. See 'Details'.
# @param fix_ids **Expert use only** A function to fix the gene_id column.
#                This function should accept a vector and return a
#                vector of the same length. See 'Details'.

CDS_from_Bioconductor <- function(txdb, orgdb,
                                  gene_id = 'ENTREZID', gene = 'SYMBOL',
                                  fix_chr = NULL, fix_ids = NULL) {

  # Join CDS coordinates and common gene symbols
  left_join(
    get_cds_coords(txdb, fix_chr, fix_ids), # from a TxDb object
    get_symbols(orgdb, gene_id, gene),      # from an OrgDb object
    by = 'GENEID'
  ) %>%
    # Arrange variables as desired
    select(
      tx  = TXNAME,   gene   = GENE,      exon  = EXONRANK,
      chr = CDSCHROM, strand = CDSSTRAND, start = CDSSTART, end = CDSEND
    ) %>%
    arrange(gene, tx, exon)
}

# ---- Utilities for CDS_from_Bioconductor ----

# Get a data frame of gene IDs and gene symbols. If gene ID maps to multiple
# symbols, symbols will be collapsed into a single string with each symbol
# separated with ' | '.
get_symbols <- function(orgdb, gene_id, gene) {

  # Mapping between entrezgene IDs and common gene symbols
  suppressMessages({
    orgdb %>%
      AnnotationDbi::select(AnnotationDbi::keys(.), c(gene_id, gene)) %>%
      select_(GENEID = gene_id, GENE = gene) %>%
      # Gene ID MUST be present
      filter(!is.na(GENEID)) %>%
      # Guarantee one row per ID
      group_by(GENEID) %>%
      summarise(GENE = stringr::str_c(na.omit(GENE) %||% '', collapse = ' | ')) %>%
      ungroup
  })
}

# Get a data frame of CDS coordinates. Allow user to fix Chromosome names and
# Gene IDs with custom functions.
get_cds_coords <- function(txdb, fix_chr = NULL, fix_ids = NULL) {
  # Gather CDS, exon, transcript and gene information
  cds <- suppressMessages({
    txdb %>%
      AnnotationDbi::select(
        txdb,
        keytype = 'CDSID',
        keys    = AnnotationDbi::keys(., 'CDSID'),
        columns = c(
          'CDSID', 'GENEID', 'TXID', 'EXONID', 'TXNAME', 'EXONRANK',
          'CDSCHROM', 'CDSSTRAND', 'CDSSTART', 'CDSEND')
      ) %>%
      filter(!is.na(GENEID))           # Some CDS IDs have no Gene ID
  })

  # User accessible fixes to CDS data
  if (!is.null(fix_chr)) cds$CDSCHROM <- fix_chr(cds$CDSCHROM)
  if (!is.null(fix_ids)) cds$GENEID   <- fix_ids(cds$GENEID)
  return(cds)
}
