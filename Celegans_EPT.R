#Setting working directory to my workspace
setwd("C:/Users/mahaz/OneDrive/Documents/Class Documents/Fall 2020/Computational Biology/Final Project/Biol_536_final_project_code")

#Importing data from the GLMPCR code
source("./Celegans_GLMPCA.R")

#Importing libraries necessary for analysis
library("ElPiGraph.R")
library(igraph)
library(magrittr)

#Creating object based off of the reduced GLMPCA data obrained from Celegans_GLMPCR.R
tree_data = data.matrix(pd, rownames.force = NA)

#Creating and visualizing elastic principal tree
TreeEPG<- computeElasticPrincipalTree(X = tree_data, Do_PCA = FALSE, NumNodes = 80)

# !!! Special Note: From here on, the code exactly mimics the procedure described in Albluca/ElPiGraph.R for visualization of graph substructures

#Construct Tree_Graph object
Tree_Graph <- ConstructGraph(PrintGraph = TreeEPG[[1]])

#Pull out substructures of the graph object
Tree_e2e <- GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
Tree_Brches <- GetSubGraph(Net = Tree_Graph, Structure = 'branches')
Tree_BrBrPt <- GetSubGraph(Net = Tree_Graph, Structure = 'branches&bpoints')
Tree_SubTrees <- GetSubGraph(Net = Tree_Graph, Structure = 'branching')
Tree_Brches_NoEnds <- GetSubGraph(Net = Tree_Graph, Structure = 'branches', KeepEnds = FALSE)

#Assign branches a certain length and generate their data
BrID <- sapply(1:length(Tree_Brches_NoEnds), function(i){
  rep(i, length(Tree_Brches_NoEnds[[i]]))}) %>%
  unlist()

#Grab the ID of the nodes in the graph and assign them to branches
NodesID <- rep(0, vcount(Tree_Graph))
NodesID[unlist(Tree_Brches_NoEnds)] <- BrID

#Associate graph points to nodes
PartStruct <- PartitionData(X = tree_data, NodePositions = TreeEPG[[1]]$NodePositions)

#Obtaining points associated with each substructure
PtInBr <- lapply(Tree_Brches, function(x){which(PartStruct$Partition %in% x)})

#Partitioning the data into specific branches and branching points
PointLabel = rep("", length(PartStruct$Partition))

for(i in 1:length(Tree_BrBrPt)){
  PointLabel[PartStruct$Partition %in% Tree_BrBrPt[[i]]] <- names(Tree_BrBrPt)[i]
}

#Visualizing the branches and branching points in the data
PlotPG(X = tree_data, TargetPG = TreeEPG[[1]], GroupsLab = PointLabel)


#Again, credit for this procedure belongs to the author of Albluca/ElPiGraph.R.
#Exact procedure followed is at https://github.com/Albluca/ElPiGraph.R/blob/master/guides/struct.md