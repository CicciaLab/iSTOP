# ---- CDS coordinates ----
#' CDS coordinates
#'
#' Compiles a data frame of coding sequence coordinates for all exons in a
#' [TxDb][GenomicFeatures::TxDb] object. Common gene names are also included
#' from an [OrgDb][AnnotationDbi::AnnotationDb] object. Note that CDS coordinates are
#' nearly equivalent to exon coordinates except the first and last exons that
#' contain coding sequence should have coordinates that begin with the
#' translation start codon and end with the stop codon.
#'
#' @param txdb    [TxDb][GenomicFeatures::TxDb]
#'                Transcript database object.
#' @param orgdb   [OrgDb][AnnotationDbi::AnnotationDb]
#'                Organism annotation database object.
#' @param gene_id [String][assertthat::is.string]
#'                Name of gene ID column in `orgdb`. See [columns][AnnotationDbi::columns].
#' @param gene    [String][assertthat::is.string]
#'                Name of gene symbol column in `orgdb`. See [columns][AnnotationDbi::columns]
#' @param fix_chr A function to fix the `chr` column.
#'                This function should accept a vector and return a
#'                vector of the same length. See 'Details'.
#' @param fix_ids A function to fix the gene_id column.
#'                This function should accept a vector and return a
#'                vector of the same length. See 'Details'.
#'
#' @details If constructing a CDS data.frame manually (1) each row should
#' represent the coordinates of a single exon, (2) exons should be numbered in
#' order with reference to the transcript's strand (i.e. the first exon should
#' include the start codon), (3) and the first and last exon coordinates should
#' beging with the start codon and end with the stop codon.
#'
#' Pre-defined CDS data.frame builders are provided. They
#' are named `CDS_<Species>_<data-source>_<genome-assembly-ID>()`
#'
#' The arguments `fix_chr` and `fix_ids` are provided for convenience
#' to fix chromosome names or gene IDs in the TxDb to match those in the
#' OrgDb. The fix functions are applied as a transformation of the TxDb `chr` and
#' `gene_id` columns. For example, this functionality is used
#' in `CDS_Athaliana_BioMart_plantsmart28()` to provide a mapping of chromosome
#' names from the TxDb format to the OrgDb format (1 = Chr1, Mt = ChrM ...). The
#' `fix_ids` functionality is used in `CDS_Dmelanogaster_UCSC_dm6()` to remove a
#' trailing '.1' that was added to every ID in the TxDb.
#'
#' @return A data.frame with the following columns where each row represents
#' a single exon:
#'
#'   - __COLUMN-NAME__ _`DATA-TYPE`_ DESCRIPTION
#'   - __gene__    _`chr`_ Gene symbol
#'   - __tx__      _`chr`_ Transcript symbol
#'   - __exon__    _`int`_ Exon rank in gene
#'   - __ncds__    _`int`_ Exon rank in CDS
#'   - __gene_id__ _`chr`_ Gene identifier
#'   - __tx_id__   _`int`_ Transcript identifier
#'   - __exon_id__ _`int`_ Exon identifier
#'   - __cds_id__  _`int`_ CDS identifier
#'   - __chr__     _`chr`_ Chromosome
#'   - __strand__  _`chr`_ Strand (+/-)
#'   - __start__   _`int`_ CDS coordinate start (always < __end__)
#'   - __end__     _`int`_ CDS coordinate end (always > __start__)
#'
#' @rdname CDS
#' @import dplyr
#' @import stringr
#' @importFrom AnnotationDbi keys
#' @export
#' @md

CDS <- function(txdb, orgdb,
                gene_id = 'ENTREZID', gene = 'SYMBOL',
                fix_chr = NULL, fix_ids = NULL) {
  # Join CDS coordinates and common gene symbols
  left_join(
      get_cds_coords(txdb, fix_chr, fix_ids), # from a TxDb object
      get_symbols(orgdb, gene_id, gene),      # from an OrgDb object
      by = 'GENEID'
    ) %>%
    # Arrange variables as desired and write to file
    select(
      gene,               tx     = TXNAME,    exon    = EXONRANK, ncds,
      gene_id = GENEID,   tx_id  = TXID,      exon_id = EXONID,   cds_id = CDSID,
      chr     = CDSCHROM, strand = CDSSTRAND, start   = CDSSTART, end    = CDSEND
    ) %>%
    arrange(gene_id, tx_id, ncds)
}

# ---- Utilities for write_cds_coordinates ----

'%||%' <- function(lhs, rhs) if (length(lhs)) lhs else rhs

# Get a data frame of gene IDs and gene symbols. If gene ID maps to multiple
# symbols, symbols will be collapsed into a single string with each symbol
# separated with ' | '.
get_symbols <- function(orgdb, gene_id, gene) {
  # Mapping between entrezgene IDs and common gene symbols
  suppressMessages({
    orgdb %>%
      AnnotationDbi::select(AnnotationDbi::keys(.), c(gene_id, gene)) %>%
      select_(GENEID = gene_id, gene = gene) %>%
      # Gene ID MUST be present
      filter(!is.na(GENEID)) %>%
      # Guarantee one row per ID
      group_by(GENEID) %>%
      summarise(gene = stringr::str_c(na.omit(gene) %||% '', collapse = ' | ')) %>%
      ungroup
  })
}

# Get a data frame of CDS coordinates. Allow user to fix Chromosome names and
# Gene IDs with custom functions.
get_cds_coords <- function(txdb, fix_chr, fix_ids) {
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
      filter(!is.na(GENEID)) %>%           # Some CDS IDs have no Gene ID
      group_by(TXID) %>%
      arrange(TXID, GENEID, EXONRANK) %>%  # Critical to arrange for CDS rank
      mutate(ncds = 1:n()) %>%             # ncds is the rank in the CDS
      ungroup
  })

  # User accessible fixes to CDS data
  if (!is.null(fix_chr)) cds$CDSCHROM <- fix_chr(cds$CDSCHROM)
  if (!is.null(fix_ids)) cds$GENEID   <- fix_ids(cds$GENEID)
  return(cds)
}

require_bioc_install <- function(require, recommend) {
  installed <- rownames(installed.packages())

  if (!all(recommend %in% installed)) {
    message(
      'Recommended package installation:\n',
      '  source("https://bioconductor.org/biocLite.R")\n',
      '  biocLite(c("', stringr::str_c(recommend, collapse = '", "'), '")'
    )
  }

  if (!all(require %in% installed)) {
    stop(
      'Required package installation:\n',
      '  source("https://bioconductor.org/biocLite.R")\n',
      '  biocLite(c("', stringr::str_c(require, collapse = '", "'), '")',
      call. = FALSE
    )
  }
}

# ---- Pre-defined CDS constructors ----
#' @rdname CDS
#' @export
CDS_Athaliana_BioMart_plantsmart28 <- function() {
  require_bioc_install(
    require   = c('TxDb.Athaliana.BioMart.plantsmart28', 'org.At.tair.db'),
    recommend = 'BSgenome.Athaliana.TAIR.TAIR9'
  )

  message('Gathering CDS coordinates for all transcripts in Athaliana BioMart TAIR9 plantsmart28...')

  CDS(
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

#' @rdname CDS
#' @export
CDS_Celegans_UCSC_ce11 <- function() {
  require_bioc_install(
    require   = c('TxDb.Celegans.UCSC.ce11.refGene', 'org.Ce.eg.db'),
    recommend = 'BSgenome.Celegans.UCSC.ce11'
  )

  message('Gathering CDS coordinates for all transcripts in Celegans UCSC ce11 entrezgene...')

  CDS(
    txdb    = TxDb.Celegans.UCSC.ce11.refGene::TxDb.Celegans.UCSC.ce11.refGene,
    orgdb   = org.Ce.eg.db::org.Ce.eg.db,
    gene_id = 'ENTREZID',
    gene    = 'SYMBOL'
  )
}

#' @rdname CDS
#' @export
CDS_Dmelanogaster_UCSC_dm6 <- function() {
  require_bioc_install(
    require   = c('TxDb.Dmelanogaster.UCSC.dm6.ensGene', 'org.Dm.eg.db'),
    recommend = 'BSgenome.Dmelanogaster.UCSC.dm6'
  )

  message('Gathering CDS coordinates for all transcripts in Dmelanogaster UCSC dm6 ensembl...')

  CDS(
    txdb    = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
    orgdb   = org.Dm.eg.db::org.Dm.eg.db,
    gene_id = 'ENSEMBL',
    gene    = 'SYMBOL',
    fix_ids = function(.) stringr::str_replace(., '\\.[0-9]*?$', '') # gene IDs have trailing '.1'
  )
}

#' @rdname CDS
#' @export
CDS_Drerio_UCSC_danRer10 <- function() {
  require_bioc_install(
    require   = c('TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db'),
    recommend = 'BSgenome.Drerio.UCSC.danRer10'
  )

  message('Gathering CDS coordinates for all transcripts in Drerio UCSC danRer10 entrezgene...')

  CDS(
    txdb    = TxDb.Drerio.UCSC.danRer10.refGene::TxDb.Drerio.UCSC.danRer10.refGene,
    orgdb   = org.Dr.eg.db::org.Dr.eg.db,
    gene_id = 'ENTREZID',
    gene    = 'SYMBOL'
  )
}

#' @rdname CDS
#' @export
CDS_Hsapiens_UCSC_hg38 <- function() {
  require_bioc_install(
    require   = c('TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'),
    recommend = 'BSgenome.Hsapiens.UCSC.hg38'
  )

  message('Gathering CDS coordinates for all transcripts in Hsapiens UCSC hg38 entrezgene...')

  CDS(
    txdb    = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    orgdb   = org.Hs.eg.db::org.Hs.eg.db,
    gene_id = 'ENTREZID',
    gene    = 'SYMBOL'
  )
}

#' @rdname CDS
#' @export
CDS_Mmusculus_UCSC_mm10 <- function() {
  require_bioc_install(
    require   = c('TxDb.Mmusculus.UCSC.mm10.knownGene', 'org.Mm.eg.db'),
    recommend = 'BSgenome.Mmusculus.UCSC.mm10'
  )

  message('Gathering CDS coordinates for all transcripts in Mmusculus UCSC mm10 entrezgene...')

  CDS(
    txdb    = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
    orgdb   = org.Mm.eg.db::org.Mm.eg.db,
    gene_id = 'ENTREZID',
    gene    = 'SYMBOL'
  )
}

#' @rdname CDS
#' @export
CDS_Rnorvegicus_UCSC_rn6 <- function() {
  require_bioc_install(
    require   = c('TxDb.Rnorvegicus.UCSC.rn6.refGene', 'org.Rn.eg.db'),
    recommend = 'BSgenome.Rnorvegicus.UCSC.rn6'
  )

  message('Gathering CDS coordinates for all transcripts in Rnorvegicus UCSC rn6 entrezgene...')

  CDS(
    txdb    = TxDb.Rnorvegicus.UCSC.rn6.refGene::TxDb.Rnorvegicus.UCSC.rn6.refGene,
    orgdb   = org.Rn.eg.db::org.Rn.eg.db,
    gene_id = 'ENTREZID',
    gene    = 'SYMBOL'
  )
}

#' @rdname CDS
#' @export
CDS_Scerevisiae_UCSC_sacCer3 <- function() {
  require_bioc_install(
    require   = c('TxDb.Scerevisiae.UCSC.sacCer3.sgdGene', 'org.Sc.sgd.db'),
    recommend = 'BSgenome.Scerevisiae.UCSC.sacCer3'
  )

  message('Gathering CDS coordinates for all transcripts in Scerevisiae UCSC sacCer3 SGD gene...')

  CDS(
    txdb    = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,
    orgdb   = org.Sc.sgd.db::org.Sc.sgd.db,
    gene_id = 'ORF',    # ORF is used instead of ENTREZID
    gene    = 'COMMON'  # COMMON is used instead of SYMBOL
  )
}
