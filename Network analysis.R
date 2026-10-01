##加载R包
library(phyloseq)
library(microeco)#Microbial Community Ecology Data Analysis
library(WGCNA) # Weighted Correlation Network Analysis
library(magrittr) # A Forward-Pipe Operator for R
#设置工作路径
setwd("F:/编程绘图学习/work/网络图7.31/11.7")
#读入OTU数据、metadata、tax表、tree文件、环境因子表
#OTU<-read.csv("OTUtable_A.csv",row.names = 1)
#metadata<-read.csv("group.csv",row.names = 1)
tax<-read.csv("Group1.csv",row.names = 1)
#tree<- read_tree("TreeA.nwk")
#env<-read.csv("理化指标A.csv",row.names = 1)

#++++++++++++++++++【方法一】将样本发现率低于20%的ASV过滤掉+++++++++++++++++++++#
#+#读入样本ASVs绝对丰度表；
ASV <- read.csv("ckotu.csv",row.names=1) 
#查看数据
ASV[1:5,1:6]
#方法二：将丰度值大于1的值替换为1,便于计算不同ASV的检测率；
ASVs_detection <- ASV
ASVs_detection[ASVs_detection >0] <- 1
#将样本发现率(检出率)低于20%的ASV过滤掉；
Reserve_ASVs <- which(rowSums(ASVs_detection)/ncol(ASVs_detection)>0.2)
length(Reserve_ASVs)
ASV3 <- ASV[Reserve_ASVs,]
#完成由绝对丰度转化为相对丰度
row = as.numeric(length(row.names(ASV3)))#行数
col = as.numeric(length(colnames(ASV3)))#列数
col_sum = rep(colSums(ASV3), row)#求和
col_sum = matrix(col_sum, nrow = col, ncol = row)#转化矩阵
#OTU/ASV相对丰度表
ASV_relative=ASV/t(col_sum)
#检查每列结果是否为1
colSums(ASV_relative)
#计算相对丰度总和
Total <-rowSums(ASV_relative)
#保留相对丰度总和高于0.1%(0.001)的ASV;即绘制≥0.1%的物种网络图
ASV_Shift <- which(Total >= 0.002)
ASV1 <- ASV[ASV_Shift,]
ASV2 = ASV1
ASV2 <- as.matrix(ASV2)#转成矩阵
dim(ASV2)#展示表行数和列数



#根据两个数据的特点，EMY使用了筛选检出率大于0.1的筛选方法，ZMY使用了筛选小于0.1%的方法
OTU<-ASV2


#设置随机变量
set.seed(20240812)

#设置ggplot主题
library(ggplot2)
theme_set(theme_bw())

#确保OTU表是data.frame
OTU<-data.frame(OTU)
class(OTU)

#查看OTU表
OTU[1:5,1:5]
class(tax)
taxonomy_table_16S[1:5,1:3]


#然后，我们创建一个 microtable 类的对象。

#添加进化树（可选）
dataset<-microtable$new(otu_table=OTU,tax_table=tax)
dataset


#为了使 OTU 和示例信息在数据集对象中的所有文件中保持一致，我们使用 function 来修剪数据集。
dataset$tidy_dataset()
print(dataset)

network<-trans_network$new(dataset=dataset,
                             cor_method="spearman",
                             use_WGCNA_pearson_spearman=TRUE
                            )

network$cal_network(COR_p_thres=0.05,#设置p阈值
                    COR_cut=0.6,COR_p_adjust = "fdr",add_taxa_name=c("phylum"),
                    delete_unlinked_nodes=TRUE,
                    COR_optimization=FALSE)#设置R值阈值


###将网络划分模块
network$cal_module(method="cluster_fast_greedy")

###计算模块的特征基因，即基于PCA分析的第一主成分，以及方差百分比
network$cal_eigen()

###计算网络属性，即拓扑机构参数
network$cal_network_attr()


###获取节点属性
network$get_node_table(node_roles=TRUE)

###获取边属性
network$get_edge_table()
network$get_adjacency_matrix()

#install.packages("rgexf")
library(rgexf)
network$save_network(filepath="CK2.gexf")
write.csv(network$res_node_table,'原产地真菌.csv')
##看枢纽节点数据
network$res_node_table
##可视化
network$plot_taxa_roles(use_type=1)
write.csv(network$res_node_table,'中心物种原产地真菌.csv')
#提取网络图相关性矩阵
edge_labels <- edge_attr(network$res_network, "label")
current_weights <- edge_attr(network$res_network, "weight")
# 2. 从 label 中提取符号信息
# 方法：检查 label 是否包含负号（"-"）
has_negative_sign <- grepl("-", edge_labels, fixed = TRUE)
# 3. 应用符号到 weight
# 如果 label 包含负号，将 weight 设为负值
signed_weights <- ifelse(has_negative_sign, -abs(current_weights), abs(current_weights))
network$res_network <- set_edge_attr(
  graph = network$res_network,
  name = "weight",
  value = signed_weights
)
cor_matrix_signed <- as_adjacency_matrix(
  network$res_network,
  attr = "weight",  # 现在包含符号
  sparse = FALSE
)
relative_abundance <- vertex_attr(network$res_network, "RelativeAbundance")

# 添加节点名称
abundance_df <- data.frame(
  Node = V(network$res_network)$name,
  RelativeAbundance = relative_abundance,
  row.names=1)

#读取网络图中的相对丰度表
comm2<-abundance_df
#计算相对丰度
sp.ra<-rowMeans(comm2)/sum(comm2)   #relative abundance of each species

#提取网络图中相关性矩阵
cormatrix2<-cor_matrix_signed
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
#检查相对丰度表和相关性矩阵是否匹配
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

#1.鲁棒性计算
## robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))

currentdat<-dat1

write.csv(currentdat,"原产地真菌鲁棒性.csv")

#计算脆弱性
cormatrix2[abs(cormatrix2)>0]<-1 # adjacency matrix
g = graph_from_adjacency_matrix(as.matrix(cormatrix2), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL) # note: this graph contains isolated nodes.
## End 2) construct a graph from OTU table

## 3) read in the correlation matrix downloaded from MENAP to construct the graph
#cormatrix.input  <- matrix(0, 91, 91)
#cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("input_files/CorrelationMatrix_Y14_W.txt")
#cormatrix <- t(cormatrix.input)
#cormatrix[abs(cormatrix)<0.8]<-0
#cormatrix[abs(cormatrix)>0]<-1 # adjacency matrix
#g = graph_from_adjacency_matrix(as.matrix(cormatrix), mode="upper", weighted = NULL, diag = FALSE, add.colnames = NULL) # note: this graph contains isolated nodes.
## End 3) read in the correlation matrix downloaded from MENAP to construct the graph
### End get graph

# remove isolated nodes
iso_node_id = which(degree(g)==0)
g2 = delete.vertices(g, iso_node_id) # graph without isolated nodes

#check node number and links
length(V(g2));length(E(g2))   
source("F:/文献阅读/通才和专才/数据/网络图/info.centrality.R")#自定义函数
# calculate vulnerability of each node
node.vul<-info.centrality.vertex(g2)
max(node.vul)

