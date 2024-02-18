# merging pdfs to pdf
library(grid)
library(png)
library(gridExtra)

dev.new()
pdf_dir <- "/mnt/local-disk/data/guoxiaolong/Renv" # where your pdf existing
pdf_files <- list.files(pdf_dir, pattern = "*.pdf", full.names = TRUE)
png_files <- list()
# convert each PDF file to a PNG file and add it to the list
for (i in 1:length(pdf_files)) {
  pdf_path <- pdf_files[i]
  png_path <- gsub(".pdf", ".png", pdf_path)
  system(paste("convert -density 300", shQuote(pdf_path), shQuote(png_path)))
  png_files[[length(png_files) + 1]] <- png_path
}
# create a list to store the grid objects
grid_list <- list()
# import each PNG file as a grid object and add it to the list
for (i in 1:length(png_files)) {
  grid_list[[i]] <- rasterGrob(readPNG(png_files[[i]]))
}
# arrange the grid objects in a grid
grid_arrange = function(grid_list, nrow, ncol) {
  grid_matrix = matrix(list(), nrow=nrow, ncol=ncol)
  for (i in 1:length(grid_list)) {
    row = ((i-1) %/% ncol) + 1
    col = ((i-1) %% ncol) + 1
    grid_matrix[[row, col]] = grid_list[[i]]
  }
  do.call("grid.arrange", c(grid_matrix, ncol=ncol, nrow=nrow))
}
pdf("merged.pdf", width = 11.69, height = 8.27, paper = "a4r") # set paper's size
grid_arrange(grid_list, nrow=1, ncol=2) # set arrangement
dev.off()
sapply(png_files, file.remove)