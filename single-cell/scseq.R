###导入包
library(ktplots)
library(beeswarm)
library(WGCNA)
library(edgeR)
library(feather)
library(dendextend)
library(monocle)
library(ggplot2)
library(dplyr)
library(Seurat)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(knitr)
library(parallel)
library(R.utils)
library(biomaRt)
#1#读取注释列表
info_file_all <- read.csv(
  "/picb/neurosys/chiyuhao/0429/qiricheng/code/info_thg_all_0612_final.csv"
) # 注释文件，每个pmid对应一组测序数据
#for(total_index in 1:nrow(info_file_all)){
#for (index in c(549, 278, 280)) {#277,279,281
  total_index <- 276        #34789849: 276~281,549
  x <- info_file_all[total_index, ] # 按行提取每行pmid注释
  source("/picb/neurosys/chiyuhao/0429/onestep_cluster1.R")
  print("program start")
  print(as.character(x["pmid"])) #以字符型打印出pmid
  file_string <- paste0(paste0(
    "/picb/neurosys/chiyuhao/0429/tonghuige/Data/",
    as.character(x["pmid"])
  ), "/") # data下是各组数据位置，此步列出每组数据的绝对路径

  pubmed_id <- as.character(x["pmid"])
  #2#获取数据所在路径
  if (as.character(x["files"]) == "T1") { # files列：F1 T1
    data_index <- as.character(x["description"]) #数据pmid号下子文件名

    temp_pubmed_data_index <- paste0(
      paste0(pubmed_id, "_"),
      data_index
    ) # pmid_子文件名
    file_path_root <- paste0(paste0(file_string, data_index), "/")#数据根目录
  } else {
    temp_pubmed_data_index <- paste0(paste0(pubmed_id, "_"), "1")
    file_path_root <- file_string
  }

    #tryCatch({
#3#获取数据信息，包括数据的物种、测序方式和标准化方式
      species_name <- as.character(x["species"])  ###数据来源物种：小鼠、人类、斑马鱼、猴子
      read_method <- "common"
      data_type <- as.character(x["seq_method"])  ###10X测序
      if (species_name == "Human") { ### 人类数据和其他物种数据分开
        data_state <- "n"
      } else {
        data_state <- "c"
      }
      log_flag <- as.character(x["log"]) # log：F1 T1
      norm_flag <- as.character(x["normalizition"]) # normalizition：F1 T1
 #4#根据物种选择一个参考数据集作为第二步有监督分类器的输入
        ref_data <- read.csv(
          "/picb/neurosys/chiyuhao/0429/ref_data/human.csv", #转录组测序参考数据
          row.names = 1
        ) #
        ref_data <- t(ref_data)
        ref_meta <- read.csv(
          "/picb/neurosys/chiyuhao/0429/ref_data/human_meta.csv", #参考元数据，统计每个细胞的UMI数等各种参数合集
          row.names = 1
        )
        ref_meta$class_label <- sub("human_", "", ref_meta$class_label)
        ref_meta$subclass_label <- sub("human_", "", ref_meta$subclass_label)


################ 基于run_consensus_clust的聚类 ####################

#1#读取文件并对数据进行预处理
      rowname_flag <- FALSE
      result <- tryCatch(
        {
          data <- read.csv(paste0(file_path_root, "data.csv"),
            row.names = 1
          ) # 读取表达数据
          rowname_flag <- TURE
        },
        error = function(e) {
        }
      )




data <- as.matrix(data) #表达矩阵
shinyFail <- (grepl("\\-",rownames(data)))|is.element(substr(rownames(data),1,1),0:9)  # rownames(data)为ENSG（Ensembl ID，Ensembl数据库的基因号）
      excludeGenes = sort(rownames(data)[grepl("LOC",rownames(data))|grepl("LINC",rownames(data))|
                                         grepl("FAM",rownames(data))|grepl("ORF",rownames(data))|grepl("KIAA",rownames(data))|
                                         grepl("FLJ",rownames(data))|grepl("DKFZ",rownames(data))|
                                         grepl("RPS",rownames(data))|grepl("RPL",rownames(data))|shinyFail])  # 排除其他类型包括核糖体基因的其他基因
      
      keepGenes    = setdiff(rownames(data),excludeGenes) 
      data = data[keepGenes,] #提取所有为ENSG的行
      data = data[,!(colSums(data)==0)] #去掉所有没有表达任何一个基因的细胞
      
      if(norm_flag=="T1"){
        norm.dat <- Matrix(data, sparse = TRUE)
      }else{
        norm.dat <- Matrix(data, sparse = TRUE) #转换为稀疏矩阵形式节省内存
      }
      
      rm(data)
      rm(excludeGenes)
      rm(shinyFail)
      rm(result)

      gc()
      if(log_flag=="T1"){
        
      }else{
        if(min(norm.dat@x)==0){
          norm.dat@x <- log2(norm.dat@x+1) #log2(0)无效
        }else{
          norm.dat@x <- log2(norm.dat@x)# 将表达数量对数化，数量为1的都以0计算，放大数据间差异
        }
      }

 #2#根据测序方式选择合适的聚类参数
 
 de.param <- de_param(
   padj.th = 0.01, #差异表达检测的p值上限
   q1.th = 0.4, #q1是每对簇中foreground中每个被检测到的基因表达的细胞比例占fore簇总细胞数的比例
   q.diff.th = 0.5, #值接近1说明某个在前景集表达的基因在背景集基本不表达，用于区分分离的细胞类型，值接近0说明某个基因在两个簇中表达的细胞比例相近，用于研究基因表达等级差异
   de.score.th = 100, #聚类比较的总差异表达分数的下限，de.score:the sum of -log10(adjusted Pvalue) for all DE genes
   min.cells = 10 #每个簇的最小细胞数
 )
 #3#进行聚类
my_dir <- "/picb/neurosys/guoxiaolong/result/"
result_file <- paste(my_dir, pubmed_id, "/", data_index, "/", sep = "")
result_dir <- paste(my_dir, pubmed_id, "/", data_index, "/", "run_cluster", sep = "")
dir.create(result_dir, recursive = TRUE)
result_1 <- run_consensus_cluster1(norm.dat, #结果为一个列表，列表第一项是一个字典，cell条码—>value,第二项是markers
  niter = 1, #迭代次数
  de.param = de.param,
  dim.method = "pca", #pca降维然后无监督分类：louvain
  output_dir = result_dir,
  mc.cores = 1 #并行计算核心数
)
group_meta <- data.frame(names(result_1$cl), result_1$cl)
colnames(group_meta) <- c("cell", "group")
new_data <- norm.dat


###########################Seurat监督分类#################################
#1#参考数据symbol转ensg
sdata <- as.data.frame(ref_data)
symbol <- rownames(sdata)

 while(TRUE){
  flag <- TRUE
  tryCatch(
    {withTimeout({ensemble <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org/")},timeout = 4,onTimeout = "warning") # 连接不稳定选择多次尝试
    }, error = function(e) {
      flag <<- FALSE
      print("ERROR")
    } ,warning = function(w){
      flag <<- FALSE
    } ,finally = {
        if (!flag) {
    next
  }
  print(flag)
  if (flag) {
    break
  }
    }
  )
 }
ensg <- getBM(attributes=c("external_gene_name","ensembl_gene_id"), filters = "external_gene_name",values = symbol, mart = ensemble)#参考数据基因名对应的ensg号
dectg <- intersect(rownames(new_data), ensg$ensembl_gene_id) #因为gene name一对多ensg，所以根据多的ensg在new_data中筛选出共有基因
new_data <- new_data[dectg, ]#筛选原数据

rownames(ensg) <- ensg$ensembl_gene_id
ensg <- ensg[dectg, ]
colnames(ensg)[1] <- 'genename'
#ensg1 <- subset(ensg, select = -ensembl_gene_id)
#ens <- unique(ensg1)
sdata$genename <- rownames(sdata)
sdata <- sdata[ensg$genename, ] ##筛选参考数据,重复的名字会被R在后面加了个'.1'，并且在提取过程中R会根据括号里的列来排序
sdata$ENSG <- ensg$ensembl_gene_id #导入ensg列
rownames(sdata) <- sdata$ENSG
sdata <- subset(sdata, select = -c(ENSG, genename))
sdata <- as.matrix(sdata) 
reference_count <- sdata
reference_meta <- ref_meta
reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example") #创建参考seurat对象
reference_seurat <- AddMetaData(reference_seurat, reference_meta) #导入参考元数据
rm(ref_data)
gc()
#2#对于大类的分类
#load query
query_count <- new_data
query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

# standard pipeline
reference_seurat <- NormalizeData(object = reference_seurat)  #标准化
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = "vst", nfeatures = 2000) 
reference_seurat <- ScaleData(reference_seurat)

query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = "vst", nfeatures = 2000)
query_seurat <- ScaleData(query_seurat)

 ##prediction
k.score <- min(30, floor(ncol(query_seurat) / 2))
k.anchor <- min(5, floor(ncol(query_seurat) / 2))
k.weight <- min(10, floor(ncol(query_seurat) / 2))
k.filter <- min(1000,floor(ncol(query_seurat) / 2))
error_flag <- F

result <- tryCatch(
  {
    sim.anchors <- FindTransferAnchors(  
      reference = reference_seurat, query = query_seurat,
      dims = 1:30, k.score = k.score, k.anchor = 10, k.filter = k.filter
    )
  },
  error = function(e) {
    error_flag <- T
  }
)

if (!exists("sim.anchors")) {
  sim.anchors <- FindTransferAnchors(
    reference = reference_seurat, query = query_seurat,
    dims = 1:30, k.score = k.score, k.anchor = k.anchor, project.query = T
  )
}

##replace Group with the actual column name from meta
k.weight <- min(floor(nrow(sim.anchors@anchors) / 2), k.weight)
k.weight <- max(k.weight, 3)  #3,2,1
predictions <- TransferData(
  anchorset = sim.anchors, refdata = reference_seurat$class_label,
  dims = 1:30, k.weight = k.weight
)
query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)
predict_meta <- data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)
rownames(predict_meta) <- predict_meta$rownames.query_seurat.meta.data.
colnames(predict_meta) <- c("group", "class_label")
for (i in unique(group_meta$group)) {
  temp <- predict_meta[rownames(group_meta[group_meta$group == i, ]), ]
  type <- names(table(temp$class_label)[table(temp$class_label) == max(table(temp$class_label))])
  if (max(table(temp$class_label)) <= (nrow(temp) / 2)) {
    predict_meta[rownames(group_meta[group_meta$group == i, ]), ]$class_label <- NA
  } else {
    predict_meta[rownames(group_meta[group_meta$group == i, ]), ]$class_label <- type
  }
}
predict_meta <- predict_meta[!is.na(predict_meta$class_label), ]
new_data <- new_data[, rownames(predict_meta)]
all_label <- unique(predict_meta$class_label)
all_meta <- predict_meta[1, ]

#3#对于亚类的分类
for (i in all_label) {
  reference_count <- sdata[, ref_meta$class_label == i]
  reference_meta <- ref_meta[ref_meta$class_label == i, ]
  reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
  reference_seurat <- AddMetaData(reference_seurat, reference_meta)

  # load query
  query_count <- new_data[, predict_meta$class_label == i]
  query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

  # standard pipeline
  reference_seurat <- NormalizeData(object = reference_seurat)
  reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = "vst", nfeatures = 2000)
  reference_seurat <- ScaleData(reference_seurat)

  query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = "vst", nfeatures = 2000)
  query_seurat <- ScaleData(query_seurat)

  ## prediction###
  k.score <- min(30, floor(ncol(query_seurat) / 2))
  k.anchor <- min(5, floor(ncol(query_seurat) / 2))
  k.weight <- min(10, floor(ncol(query_seurat) / 2))
  k.filter <- min(1000,floor(ncol(query_seurat) / 2))
  error_flag <- F
  result <- tryCatch(
    {
      sim.anchors <- FindTransferAnchors(
        reference = reference_seurat, query = query_seurat,
        dims = 1:30, k.score = k.score, k.anchor = 10, k.filter = k.filter
      )
    },
    error = function(e) {
      error_flag <- T
    }
  )
  if (!exists("sim.anchors")) {
    sim.anchors <- FindTransferAnchors(
      reference = reference_seurat, query = query_seurat,
      dims = 1:30, k.score = k.score, k.anchor = k.anchor, project.query = T
    )
  }
  ## replace Group with the actual column name from meta
  k.weight <- min(floor(nrow(sim.anchors@anchors) / 2), k.weight)
  k.weight <- max(k.weight, 3)
  predictions <- TransferData(
    anchorset = sim.anchors, refdata = reference_seurat$subclass_label,
    dims = 1:30, k.weight = k.weight
  )
  query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

  predict_meta1 <- data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)

  rownames(predict_meta1) <- predict_meta1$rownames.query_seurat.meta.data.
  colnames(predict_meta1) <- c("group", "class_label")
  all_meta <- rbind(all_meta, predict_meta1)
}

####
all_meta <- all_meta[-1, ]
rownames(all_meta) <- all_meta$group
colnames(all_meta) <- c("cell", "subclass")
all_meta_temp <- all_meta
for (i in unique(group_meta$group)) {
  temp <- all_meta[rownames(group_meta[group_meta$group == i, ]), ]
  type <- names(table(temp$subclass)[table(temp$subclass) == max(table(temp$subclass))])
  if (max(table(temp$subclass)) <= (nrow(temp) / 2)) {
    all_meta[rownames(group_meta[group_meta$group == i, ]), ]$subclass <- NA
  } else {
    all_meta[rownames(group_meta[group_meta$group == i, ]), ]$subclass <- type
  }
}
all_meta <- all_meta[!is.na(all_meta$subclass), ]

data_meta <- data.frame(all_meta$cell, predict_meta[rownames(all_meta), ]$class_label, all_meta$subclass, group_meta[rownames(all_meta), ]$group)
rownames(data_meta) <- data_meta$all_meta.cell
colnames(data_meta) <- c("cell", "class", "subclass_label", "cluster_label")
data_meta$cluster_label <- paste0("a", data_meta$cluster_label)
data_meta <- data_meta[data_meta$cluster_label %in% names(table(data_meta$cluster_label))[table(data_meta$cluster_label) > 1], ]
new_data <- new_data[, rownames(data_meta)]

###################分类结束##################

      ######################################
      ######3. 计算差异表达基因并且寻找marker
      ######################################
      
      mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
      mouse_data <- AddMetaData(mouse_data, data_meta)
      Idents(mouse_data) <- mouse_data$subclass_label
      all_marker_list = list()
      num = 1
      for(i in unique(mouse_data$subclass_label)){
        temp_mouse_data = subset(mouse_data, idents = i)
        Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
        #plan(workers = 6)
        mouse_cells_markers <- FindAllMarkers(temp_mouse_data, test.use = "roc",densify=T)
        mouse_cells_markers = mouse_cells_markers[mouse_cells_markers$avg_log2FC>0,]
        all_marker_list[[num]] = mouse_cells_markers
        num = num + 1
      }
      names(all_marker_list) = unique(mouse_data$subclass_label)
      saveRDS(all_marker_list, paste0(result_file,"de_gene.rds"))
      
      
      cl = mouse_data@meta.data$cluster_label
      names(cl) = rownames(mouse_data@meta.data)
      temp_exp = as.matrix(mouse_data@assays$RNA@data)
      medianExpr = do.call("cbind", tapply(names(cl), cl, function(x) rowMedians(temp_exp[,x])))
      rm(temp_exp)
      gc()
      rownames(medianExpr) = rownames(mouse_data@assays$RNA@data)
      
      
      median_marker_list = list()
      num=1
      yes_num = 0
      cluster_names = c()
      for(i in unique(mouse_data@meta.data$subclass_label)){
        temp_cluster = unique(mouse_data@meta.data[mouse_data@meta.data$subclass_label==i,]$cluster_label)
        temp_expr = medianExpr[,temp_cluster]
        cluster_names <- c(cluster_names, temp_cluster)
        if(length(temp_cluster)==1){
          median_marker_list[[num]] = NA
          num = num + 1
        }else{
          for(j in temp_cluster){
            current = temp_expr[,j]
            if(sum(colnames(temp_expr)!=j)==1){
              other = temp_expr[,colnames(temp_expr)!=j]
            }else{
              other = rowMax(temp_expr[,colnames(temp_expr)!=j])
            }
            names(other) = rownames(temp_expr)
            temp_gene = rownames(temp_expr)[current>=other]
            if(length(temp_gene)==0){
              temp_gene = NA
            }else{
              yes_num = yes_num + 1
            }
            median_marker_list[[num]] = temp_gene
            
            num = num + 1
          }
        }
        
      }
      for(i in 1:length(median_marker_list)){
        print(length(median_marker_list[[i]]))
      }
      
      temp_all_marker_list = all_marker_list
      names(median_marker_list) = cluster_names
      for(i in 1:length(all_marker_list)){
        temp = temp_all_marker_list[[i]][1,]
        for(j in unique(all_marker_list[[i]]$cluster)){
          temp_median_gene = median_marker_list[[j]]
          temp = rbind(temp, all_marker_list[[i]][all_marker_list[[i]]$cluster==j,][all_marker_list[[i]][all_marker_list[[i]]$cluster==j,]$gene %in%intersect(all_marker_list[[i]]$gene,temp_median_gene),])
          print(length(intersect(all_marker_list[[i]]$gene,temp_median_gene)))
        }
        temp_all_marker_list[[i]] = temp[-1,]
      }
      
      all_marker_list = temp_all_marker_list
      
      
      
      Idents(mouse_data) <- mouse_data$subclass_label
      num = 1
      binary_data_frame = data.frame(0,0,0)
      for(i in unique(mouse_data$subclass_label)){
        print(num)
        temp_mouse_data = subset(mouse_data, idents = i)
        Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
        temp_mouse_data_exp = as.matrix(temp_mouse_data@assays$RNA@data)
        temp_cluster = as.character(unique(temp_mouse_data@meta.data$cluster_label))
        
        for(j in temp_cluster){
          temp_score = 0
          temp_current_marker_list_temp = all_marker_list[[num]][all_marker_list[[num]]$cluster==j,]
          temp_current_marker_list = temp_current_marker_list_temp[temp_current_marker_list_temp$avg_log2FC>0,]$gene
          for(g in temp_current_marker_list){
            temp_medianExpr = medianExpr[,temp_cluster]
            # median_sec_median_score = max(temp_medianExpr[g,colnames(temp_medianExpr)!=j]) / temp_medianExpr[g,j]
            # if(median_sec_median_score == Inf){
            #   median_sec_median_score = 1
            # }
            # median_sec_median_score = 1 - median_sec_median_score
            
            temp_score = 0
            other_temp_cluster = temp_cluster[temp_cluster!=j]
            temp_median = median(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])])
            temp_percentage = sum(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])]!=0) / length(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])])
            all_temp_other_percentage = c()
            for(k in other_temp_cluster){
              temp_other_median = median(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])])
              temp_other_percentage = sum(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])]==0)
              all_temp_other_percentage = c(all_temp_other_percentage, temp_other_percentage)
              temp222 = (1 - temp_other_median / temp_median) * temp_other_percentage/length(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])])
              if(!is.na(temp222)){
                if(temp222 < 0){
                  temp222 = 0
                }
              }
              
              temp_score = temp_score + temp222
            }
            temp_score = temp_score / (length(temp_cluster) - 1)
            binary_data_frame = rbind(binary_data_frame, c(j, temp_score,g))
          }
          
        }
        
        num = num + 1
      }
      
      binary_data_frame$X0.1 = as.numeric(binary_data_frame$X0.1)
      binary_data_frame = binary_data_frame[!is.nan(binary_data_frame$X0.1),]
      binary_data_frame = binary_data_frame[!is.infinite(binary_data_frame$X0.1),]
      
      binary_data_frame = binary_data_frame[-1,]
      
      all_fc_list <- c()
      for(i in unique(binary_data_frame$X0)){
        temp_binary_data_frame = binary_data_frame[binary_data_frame$X0==i,]
        for(j in 1:length(all_marker_list)){
          if(i %in% all_marker_list[[j]]$cluster){
            temp_marker_list = all_marker_list[[j]][all_marker_list[[j]]$cluster==i,]
            temp_fc_list = temp_marker_list[temp_binary_data_frame$X0.2,]$avg_log2FC
            all_fc_list <- c(all_fc_list, temp_fc_list)
          }
        }
      }
      binary_data_frame <- cbind(binary_data_frame, all_fc_list)
      binary_gene = data.frame(0,0,0,0,0)
      for(i in unique(binary_data_frame$X0)){
        temp = binary_data_frame[binary_data_frame$X0==i,]
        binary_gene = rbind(binary_gene, temp$X0.2[1:5])
      }
      binary_gene = binary_gene[-1,]
      rownames(binary_gene) = unique(binary_data_frame$X0)
      
      data_meta_cluster_level = unique(data.frame(data_meta$class,data_meta$subclass_label,data_meta$cluster_label))
      rownames(data_meta_cluster_level) = data_meta_cluster_level$data_meta.cluster_label
      colnames(data_meta_cluster_level) = c("class","subclass","group")
      data_meta_final = data.frame(0,0,0,0,0,0,0,0)
      for(i in 1:nrow(data_meta_cluster_level)){
        if(rownames(data_meta_cluster_level)[i] %in% rownames(binary_gene)){
          temp = c(as.character(data_meta_cluster_level[i,]), as.character(binary_gene[rownames(data_meta_cluster_level)[i],]))
        }else{
          temp = c(as.character(data_meta_cluster_level[i,]), NA,NA,NA,NA,NA)
        }
        data_meta_final = rbind(data_meta_final, temp)
      }
      data_meta_final = data_meta_final[-1,]
      
      colnames(data_meta_final) = c("class","subclass","cluster","gene1","gene2","gene3","gene4","gene5")
      #write.csv(data_meta_final, paste0(result_file,"cluster_table.csv"))
      
      
      
      temp <- data.frame(0,0)
      for(i in unique(data_meta_final$subclass)){
        subclass_num = 1
        temp_data = data_meta_final[data_meta_final$subclass==i,]
        for(j in 1:length(temp_data[,1])){
          if(!is.na(temp_data[j,"gene1"])){
            temp_name = paste0(paste0(paste0(paste0(temp_data[j,"class"],"_"),temp_data[j,"subclass"]),"_"),temp_data[j,"gene1"])
            temp <- rbind(temp, c(temp_data[j,"cluster"], temp_name))
          }else{
            temp <- rbind(temp, c(temp_data[j,"cluster"], paste0(paste0(paste0(paste0(paste0(paste0(temp_data[j,"class"],"_"),temp_data[j,"subclass"]),"_"),temp_pubmed_data_index),"_"),subclass_num)))
            subclass_num = subclass_num + 1
          }
        }
      }
      temp = temp[-1,]
      rownames(temp) = temp$X0
      data_meta_final$cluster_new = temp[data_meta_final$cluster,2]
      write.csv(data_meta_final, paste0(result_file,"cluster_table.csv"))
      for(i in unique(data_meta$cluster_label)){
        data_meta[data_meta$cluster_label==i,"cluster_label"] = temp[i,2]
      }
      write.csv(data_meta, paste0(result_file,"new_anno.csv"))
      
      write.csv(new_data, paste0(result_file,"new_data.csv"))
      
      
      ######################################
      ######3. 计算差异表达基因并且寻找marker结束
      ######################################
            ######################################
      ######4. 画umap图
      ######################################
      
    
      
      query_seurat <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
      query_seurat <- AddMetaData(object = query_seurat, metadata = data_meta)
      Idents(query_seurat) = query_seurat$cluster_label
      
      query_seurat <- NormalizeData(object = query_seurat)
      query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
      query_seurat <- ScaleData(query_seurat)
      
      query_seurat <- RunPCA(query_seurat, features = VariableFeatures(object = query_seurat),npcs=100)
      query_seurat <- RunUMAP(query_seurat, dims = 1:100)
      write.csv(query_seurat@reductions$umap@cell.embeddings, paste0(result_file,"umap_index.csv"))
      saveRDS(query_seurat, paste0(result_file,"umap.rds"))
      ######################################
      ######4. 画umap图结束
      ######################################
      print("program end")
}