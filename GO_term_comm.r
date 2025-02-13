# 2022-08-19
# marieke kuijjer
# GO term enrichment analysis of Genis' test network

rm(list=ls())
library(topGO)
#BiocManager::install("biomaRt")
#library(biomaRt)

#home <- "/home/mariekek/Dropbox/work/group/research_assistants/genis_calderer_i_garcia_(2021-2022, research assistant)/BiHiDeF_data"
home <- "/storage/kuijjerarea/romana/eland/ELAND/"

setwd(home)

sign <- utils::read.table("results/pvg.nodes", header = FALSE, sep = "\t")

# column 1: cluster and its community
clusters <- as.character(sign[,1])
clusters <- substr(clusters, 8, nchar(clusters))
clusters <- t(matrix(as.numeric(unlist(c(strsplit(clusters, "-")))),2))
table(clusters[,1]) # from 1 to 82 communities in the different layers
sum(table(clusters[,1]))==nrow(clusters) # TRUE, so all communities are listed

# column 2: number of genes in community
summary(sign[,2]) # from 1 to 756 genes

allgenes <- vector()
for(i in 1:nrow(sign)){
	genesincom <- sign[i,3]
	genesincom <- unlist(c(strsplit(as.character(genesincom), " ")))
	#genesincom <- substr(genesincom, 5, nchar(genesincom))
	allgenes <- c(allgenes,genesincom)
	print(i)
}
length(allgenes) # 3963
length(unique(allgenes)) # 756
# so there are only 756 genes in total???

# i should used our gtex annotation as that is what these networks are based on
anno <- read.delim("data/GTEx_gene_names.txt") #  from gtex code review NEW2 folder
annosub <- anno[which(row.names(anno) %in% allgenes),]
allgenes <- substr(allgenes,1,15)
table(nchar(allgenes)) # 16 or 17 characters
table(nchar(row.names(anno))) # all 15 characters


# looking at these genes, it seems genis used networks with cancer risk snps (e.g. there are many HLAs in here)


# ### below, repeated on snps just to see how many there are
# sign <- utils::read.table("pvg.nodes", header = FALSE, sep = "\t")
# clusters <- as.character(sign[,1])
# clusters <- substr(clusters, 8, nchar(clusters))
# clusters <- t(matrix(as.numeric(unlist(c(strsplit(clusters, "-")))),2))
# table(clusters[,1]) # from 1 to 82 communities in the different layers
# sum(table(clusters[,1]))==nrow(clusters) # TRUE, so all communities are listed
# summary(sign[,2]) # from 1 to 35870 snps
# allgenes <- vector()
# for(i in 1:nrow(sign)){
# 	genesincom <- sign[i,3]
# 	genesincom <- unlist(c(strsplit(as.character(genesincom), " ")))
# 	genesincom <- substr(genesincom, 5, nchar(genesincom))
# 	allgenes <- c(allgenes,genesincom)
# 	print(i)
# }
# length(allgenes) # 
# length(unique(allgenes)) # 35867 snps in total


### this is how you can read in a network from networkmedicine.org

	#dat <- readRDS("breast_mammary_tissue.rds")
	#dat <- dat$qtls
	#length(unique(dat[,4])) # indeed 756 genes only, and indeed the ensembl gene ids have a dot and this was removed in genis script

	#dat <- readRDS("skin.rds")
	#dat <- dat$qtls
	#length(unique(dat[,4])) # 9460 genes
	#length(unique(dat[,1])) # 37k snps


### GO term enrichment of breast network (only 756 genes total, will compare all with entire gtex gene background)

	# 1. GO term enrichment with classical test
		background <- anno[,1]
		GO2Symbol.bp <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol") # Human genes, attach GO to it
		symbol2GO.bp <- inverseList(GO2Symbol.bp)
		dir.create("GOresults_classic")
		res <- matrix(NA, nrow(sign), 5)
		colnames(res) <- c("community","sign", "minadp", "topGOID", "topTerm")
		res[,1] <- as.character(sign[,1])
		for(i in 1:nrow(sign)){
			genesincom <- sign[i,3]
			genesincom <- unlist(c(strsplit(as.character(genesincom), " ")))
			#genesincom <- substr(genesincom, 5, nchar(genesincom))
			genesincom <- unique(substr(genesincom, 1, 15))
			interesting <- as.character(anno[which(anno[,1] %in% genesincom),1])

			if(length(interesting) == 0){
				next
			} else {

		    selection = factor(as.integer(background %in% interesting))
		    geneList <- factor(as.integer(background %in% interesting))
		    names(geneList) = toupper(background)
		    GOdata.bp <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = selection, annot = annFUN.gene2GO, gene2GO = symbol2GO.bp)
		    
		    resultFisher = runTest(GOdata.bp, algorithm="classic", statistic="fisher")  
		    enrichRes.bp <- GenTable(GOdata.bp, classic = resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = length(score(resultFisher)))
		    adjP_Fisher <- p.adjust(enrichRes.bp$classic,"BH")
		    enrichRes.bp <- cbind(enrichRes.bp, adjP_Fisher)
		    enrichRes.bp <- enrichRes.bp[order(enrichRes.bp$adjP_Fisher),]
		    write.table(enrichRes.bp,file=paste("GOresults_classic/", as.character(sign[i,1]), ".txt", sep="") ,row.names = F,col.names = T, quote = F,sep="\t")

		   	res[i,2] <- length(which(enrichRes.bp$adjP_Fisher<0.05))
		   	res[i,3] <- min(enrichRes.bp$adjP_Fisher)
		   	res[i,4] <- enrichRes.bp[1,1] # this makes no sense. it should always take the first row as that's the lowest pval
		   	res[i,5] <- enrichRes.bp[1,2]

				print(res[i,])
		write.table(res, "GOresults_summary.txt", sep="\t", quote=F)
			}
		#write.table(res, "GOresults_summary.txt", sep="\t", quote=F)
			}


	# 2. with elim test (likely most interesting as it finds the best fitting pathways)
		dir.create("GOresults_elim")
		res <- matrix(NA, nrow(sign), 5)
		colnames(res) <- c("community","sign", "minadp", "topGOID", "topTerm")
		res[,1] <- as.character(sign[,1])
		for(i in 1:nrow(sign)){
			genesincom <- sign[i,3]
			genesincom <- unlist(c(strsplit(as.character(genesincom), " ")))
			#genesincom <- substr(genesincom, 5, nchar(genesincom))
			genesincom <- unique(substr(genesincom, 1, 15))
			interesting <- as.character(anno[which(anno[,1] %in% genesincom),1])
			if(length(interesting) == 0){
				next
			} else {
		    selection = factor(as.integer(background %in% interesting))
		    geneList <- factor(as.integer(background %in% interesting))
		    names(geneList) = toupper(background)
		    GOdata.bp <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = selection, annot = annFUN.gene2GO, gene2GO = symbol2GO.bp)
		    
		    resultFisher = runTest(GOdata.bp, algorithm="elim", statistic="fisher")
		    enrichRes.bp <- GenTable(GOdata.bp, elim = resultFisher, orderBy = "elim", ranksOf = "elim", topNodes = length(score(resultFisher)))
		    adjP_Fisher <- p.adjust(enrichRes.bp$elim,"BH")
		    enrichRes.bp <- cbind(enrichRes.bp, adjP_Fisher)
		    enrichRes.bp <- enrichRes.bp[order(enrichRes.bp$adjP_Fisher),]
		    write.table(enrichRes.bp,file=paste("GOresults_elim/", as.character(sign[i,1]), ".txt", sep="") ,row.names = F,col.names = T, quote = F,sep="\t")

		   	res[i,2] <- length(which(enrichRes.bp$adjP_Fisher<0.05))
		   	res[i,3] <- min(enrichRes.bp$adjP_Fisher)
		   	res[i,4] <- enrichRes.bp[1,1]  # this makes no sense. it should always take the first row as that's the lowest pval
		   	res[i,5] <- enrichRes.bp[1,2]

				print(res[i,])
				write.table(res, "GOresults_summary_elim.txt", sep="\t", quote=F)
			}
		}


### cytoscape
# i imported pvr.edges into cytoscape, together with the GO term summary results (top elim GO term per pathway with FDR<0.25)
# what i could do now is include the amount of overlapping pathways between the communities
	# however, it would be best to base that on the classic test
# workflow
	# 1. loop over all of the cluster files in classic, select pathways with FDR<0.25 (for now) and combine all into one file
		setwd("GOresults_classic")
		filenames <- list.files()
		signgoterms <- matrix(NA,0,2)

		for (f in 1:length(filenames)){
			gores <- read.delim(filenames[f])
			gores <- as.character(gores[which(gores$adjP_Fisher<0.25),,drop=F]$GO.ID)
			gosave <- cbind(substr(filenames[f], 1, nchar(filenames[f])-4), gores)
			if(nrow(gosave)>1){
				signgoterms <- rbind(signgoterms, gosave)
			}
		print(f)
		} # end loop over f filenames

	# some stats
		nrow(signgoterms) # 1469 in total
		table(signgoterms[,1]) # 118 36 187 etc
		length(table(signgoterms[,1])) # 29 so make pairwise comparisons for all parent-child combinations of these clusters

			# OBS! i am now using cytoscape to check the dependencies. this is probably not feasible if this tree is going to be much larger
			# parent child
			# Cluster0-0 Cluster1-0
			# Cluster0-0 Cluster1-2
			# Cluster0-0 Cluster1-4
			# Cluster1-0 Cluster2-0
			# Cluster1-2 Cluster2-1
			# Cluster1-4 Cluster2-10
			# Cluster1-4 Cluster2-14
			# Cluster1-0 Cluster2-2
			# Cluster1-2 Cluster2-6
			# Cluster2-2 Cluster3-1
			# Cluster2-1 Cluster3-18
			# Cluster2-0 Cluster3-2
			# Cluster2-6 Cluster3-20
			# Cluster2-1 Cluster3-25
			# Cluster2-0 Cluster3-3
			# Cluster2-1 Cluster3-4
			# Cluster2-1 Cluster3-6
			# Cluster3-6 Cluster4-15
			# Cluster3-3 Cluster4-2
			# Cluster3-1 Cluster4-22
			# Cluster3-4 Cluster4-4
			# Cluster3-3 Cluster4-5
			# Cluster4-15 Cluster5-10
			# Cluster4-21 Cluster5-23 -->> OBS! Cluster4-21 is not enriched! (but that is OK, i do not have to look into all possible dependencies, as none of the other nodes are enriched)
			# Cluster4-15 Cluster5-5
			# Cluster4-5 Cluster5-6
			# Cluster5-10 Cluster6-14
			# Cluster5-10 Cluster6-5

# use object "clusters" from above to screen though all
edges <- read.delim("../pvr.edges")
jacres <- cbind(edges[,1:2], 0) # default is 0
for(e in 1:nrow(edges)){
	f_parent <- paste(as.character(edges[e,1]), ".txt", sep="")
	f_child <- paste(as.character(edges[e,2]), ".txt", sep="")
	go_parent <- read.delim(f_parent)
	go_child <- read.delim(f_child)	
	go_parent <- as.character(go_parent[which(go_parent$adjP_Fisher<0.25),,drop=F]$GO.ID)
	go_child <- as.character(go_child[which(go_child$adjP_Fisher<0.25),,drop=F]$GO.ID)
	jacres[e,3] <- length(intersect(go_parent, go_child))/length(union(go_parent, go_child))
	print(e)
}
colnames(jacres) <- c("parent","child","jaccard")
write.table(jacres, "../GOclassic_edges_jaccard_FDR025.txt", sep="\t", quote=F, row.names=F)




gosave <- cbind(substr(filenames[f], 1, nchar(filenames[f])-4), gores)

	# 4. count the jaccard index for overlapping pathways for each combination of parent-child clusters
	# 5. that will be an edge attribute to overlay in cytoscape