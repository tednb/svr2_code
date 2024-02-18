library(GenomicRanges)
library(data.table)
# Function to read and parse the data
read_and_parse_data <- function(filename) {
  # 使用 fread 读取文件
  dt <- fread(filename, header = FALSE, col.names = c("info", "score"))
  
  # 拆分 info 列的数据
  info_cols <- tstrsplit(dt$info, '[:_\\$\\-]', type.convert = TRUE)    # 使用正则表达式
  
  # 将列表转换为数据框
  df <- data.frame(
    chromosome = info_cols[[1]],
    start = as.integer(info_cols[[2]]),
    end = as.integer(info_cols[[3]]),
    gene_id = info_cols[[4]],
    gene_name = info_cols[[5]],
    gene_chromosome = info_cols[[6]],
    gene_tss = as.integer(info_cols[[7]]),
    score = as.numeric(dt$score)  # 将得分列添加到数据框中
  )
  
  return(na.omit(df))
}

# Function to create a GRanges object from the parsed data
create_granges <- function(parsed_data) {
  # Create a data frame from the parsed data
  df <- do.call(rbind, lapply(parsed_data, as.data.frame))
  
  # Convert it to a GRanges object
  gr <- GRanges(
    seqnames = Rle(df$chromosome),
    ranges = IRanges(start = df$start, end = df$end),
    gene_id = df$gene_id,
    gene_name = df$gene_name,
    gene_chromosome = df$gene_chromosome,
    gene_tss = df$gene_tss,
    score = df$score
  )
  
  # Sort and reduce to merge overlapping ranges
  gr <- reduce(gr)
  
  # Return the GRanges object
  return(gr)
}

# Function to filter GRanges by highest score for identical enhancers
filter_highest_score <- function(gr) {
  # Order by score and then by enhancer identifier (chromosome, start, end)
  gr <- gr[order(gr$score, decreasing = TRUE)]
  gr <- gr[!duplicated(mcols(gr)[, c("gene_id", "gene_name", "gene_chromosome", "gene_tss")])]
  
  # Return the filtered GRanges object
  return(gr)
}

# Main execution function
process_data_to_granges <- function(filename) {
  # Read and parse the data
  parsed_data <- read_and_parse_data(filename)
  
  # Create a GRanges object
  gr <- create_granges(parsed_data)
  
  # Filter by highest score
  gr <- filter_highest_score(gr)
  
  # Return the final GRanges object
  return(gr)
}

# Replace with your actual file path
filename <- "~/chain/Enhancer/Uterus_EP.txt"

# Process the data and get the GRanges object
granges_object <- process_data_to_granges(filename)

# Print the GRanges object
print(granges_object)
