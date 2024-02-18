library(dplyr)
library(tibble)
library(rtracklayer)
library(limma)

EnsemblID_to_symbol <- function(tpm_matrix, gtf_file_path) {
  # Convert TPM matrix to data frame and add row names as new column
  tpm_df <- as.data.frame(tpm_matrix) %>%
    tibble::rownames_to_column(var = "Ensembl_ID")
  
  # Import GTF file and filter for relevant columns
  gtf_data <- rtracklayer::import(gtf_file_path) %>%
    as.data.frame() %>%
    dplyr::select(gene_id, type,gene_type, gene_name, seqnames, start, end, strand) %>%
    dplyr::filter(gene_type == "protein_coding" & type == "transcript") %>%
    distinct()
  
  # Join TPM data with GTF data and average replicates
  expr_tpm_averaged <- tpm_df %>%
    dplyr::inner_join(gtf_data, by = c("Ensembl_ID" = "gene_id")) %>%
    dplyr::select(gene_name, starts_with("ENCS")) %>%
    {avereps(.[,-1], ID = .$gene_name)} %>% 
    as.data.frame()
  
  # Get a vector of unique gene names that are present in the averaged expression data
  genes_in_expr_tpm_averaged <- unique(rownames(expr_tpm_averaged))
  
  # Filter the GTF data to only include genes present in the averaged expression data
  gene_location_info <- gtf_data %>%
    dplyr::filter(gene_name %in% genes_in_expr_tpm_averaged) %>%
    dplyr::select(gene_name, seqnames, start, end, strand) %>%
    dplyr::distinct()

  # Set the row names of the location info data frame to be the gene_id
  colnames(gene_location_info) <- c("mgi_symbol","chromosome_name","start_position","end_position","strand")
  # Return a list containing both the averaged expression and location information
  return(list(expr_tpm_averaged = expr_tpm_averaged,
              gene_location_info = gene_location_info))
}




# Example usage:
# result_expr_tpm <- process_tpm_data(tpm_matrix = tpm.m, gtf_file_path = '~/chain/gencode.vM21.basic.annotation.gtf.gz')
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
retry_useMart <- function(mart, dataset, host, max_attempts = 5) {
  attempt <- 1
  while (attempt <= max_attempts) {
    tryCatch({
      ensembl <- useMart(mart, dataset = dataset, host = host)
      return(ensembl)
    }, error = function(e) {
      cat(paste("Attempt", attempt, "failed:", conditionMessage(e), "\n"))
      attempt <- attempt + 1
      Sys.sleep(5)  # Wait for 5 seconds before retrying
    })
  }
  stop("Exceeded maximum attempts")
}
