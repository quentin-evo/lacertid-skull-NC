library("bayestraitr")
library("coda")
library("ape")
library("phytools")
library("tidytree")
library("tidyverse")

#load("G:/macroevo/lacertids.Rdata")
log <- bt_read.log("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PCAscores.txt.Log.txt")

log2 <- read.table("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PCAscores.txt.Log.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Convert to MCMC object
log.mcmc <- as.mcmc(log)


plot(log.mcmc)

output <-  read.csv("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PCAscores.txt.VarRates.txt.csv", header=TRUE)

t1 <- read.nexus("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/software-files/Artiodactyl.trees")

tree <- red.data$phy 
output$ID <- 1:347  # Shift IDs to match R's 1-based index
output <- output[-1,]

# Identify branches meeting the conditions
red_branches <- tree$edge[output$P...1.1 > 0.79, 2]  # Extract node numbers for red dots
blue_branches <- tree$edge[output$P...1 > 0.79, 2]   # Extract node numbers for blue dots

# Plot the phylogenetic tree as a fan
plot(tree, type = "fan", show.tip.label = T, edge.width = 2.5, cex = 0.5)

# Add red and blue dots at branch midpoints
nodelabels(pch = 16, col = "orange", cex = 1.2, node = red_branches)
nodelabels(pch = 16, col = "blue", cex = 1.2, node = blue_branches)



#### PCA - NC ####

# logNC <- bt_read.log("G:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PC-scoresNC.txt.Log.txt")
# logNC.mcmc <- as.mcmc(logNC)
# plot(logNC.mcmc)

outputNC <-  read.csv("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PC-scoresNC.txt.VarRates.txt.csv", header=TRUE)

treeNC <- read.nexus("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/lacertid-nigri.nex") 
outputNC$ID <- 1:349  # Shift IDs to match R's 1-based index
outputNC <- outputNC[-1,]

# Identify branches meeting the conditions
red_branchesNC <- treeNC$edge[outputNC$P...1.1 > 0.80, 2]  # Extract node numbers for red dots
blue_branchesNC <- treeNC$edge[outputNC$P...1 > 0.80, 2]   # Extract node numbers for blue dots

# Initialize edge colors as black (default)
edge_colors <- rep("black", nrow(treeNC$edge))

# Function to get all descendant edges of a given node
get_descendants <- function(tree, node) {
  descendants <- c()
  to_visit <- node
  while (length(to_visit) > 0) {
    current <- to_visit[1]
    to_visit <- to_visit[-1]
    
    child_edges <- which(tree$edge[, 1] == current)  # Find edges where current is the parent
    if (length(child_edges) > 0) {
      child_nodes <- tree$edge[child_edges, 2]  # Get descendant nodes
      descendants <- c(descendants, child_edges)  # Store descendant edges
      to_visit <- c(to_visit, child_nodes)  # Add to visit queue
    }
  }
  return(descendants)
}

# Identify edges that meet the conditions
red_edges <- which(outputNC$P...0.1 > 0.8)  # Edges where P...0.1 > 0.8
blue_edges <- which(outputNC$P...0 > 0.8)   # Edges where P...0 > 0.8

# Function to propagate color without overwriting an already colored edge
propagate_color <- function(tree, edges, color) {
  for (edge in edges) {
    descendants <- get_descendants(tree, tree$edge[edge, 2])
    for (desc in descendants) {
      if (edge_colors[desc] == "black") {  # Only change if it's still black
        edge_colors[desc] <- color
      }
    }
  }
}

# Assign initial colors
edge_colors[red_edges] <- "orange2"
edge_colors[blue_edges] <- "blue"

# First propagation without overwriting
propagate_color(treeNC, red_edges, "orange2")
propagate_color(treeNC, blue_edges, "blue")

# Function to color remaining black edges based on closest ancestor's color
color_remaining_black_edges <- function(tree, edge_colors) {
  for (i in seq_along(edge_colors)) {
    if (edge_colors[i] == "black") {  # Only process black edges
      parent <- tree$edge[i, 1]  # Find parent node
      while (parent > 0) {  # Traverse up until we find a colored ancestor
        ancestor_edge <- which(tree$edge[, 2] == parent)
        if (length(ancestor_edge) > 0) {
          ancestor_color <- edge_colors[ancestor_edge]
          if (ancestor_color != "black") {
            edge_colors[i] <- ancestor_color  # Assign ancestor's color
            break
          }
          parent <- tree$edge[ancestor_edge, 1]  # Move up to next ancestor
        } else {
          break  # Stop if no ancestor found
        }
      }
    }
  }
  return(edge_colors)
}

# Apply the function to color remaining black edges
edge_colors <- color_remaining_black_edges(treeNC, edge_colors)

# Plot the phylogenetic tree as a fan
par(mar = c(0.1,0.1,0.1,0.1))
plot(treeNC, type = "fan", show.tip.label = T, edge.width = 2.5, cex = 0.7 , edge.color = edge_colors)

# Add red and blue dots at branch midpoints
nodelabels(pch = 16, col = "orange", cex = 1.2, node = red_branchesNC)
nodelabels(pch = 16, col = "blue", cex = 1.2, node = blue_branchesNC)

##### PC MES ####

# logMES <- bt_read.log("G:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PC-scoresMES.txt.Log.txt")
# log.mcmcMES <- as.mcmc(logMES)
# plot(log.mcmcMES)

outputMES <-  read.csv("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PC-scoresMES.txt.VarRates.txt.csv", header=TRUE)

treeMES <- read.nexus("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/lacertid-nigri.nex") 
outputMES$ID <- 1:349  # Shift IDs to match R's 1-based index
outputMES <- outputMES[-1,]

# Identify branches meeting the conditions
red_branchesMES <- treeMES$edge[outputMES$P...1.1 > 0.80, 2]  # Extract node numbers for red dots
blue_branchesMES <- treeMES$edge[outputMES$P...1 > 0.80, 2]   # Extract node numbers for blue dots

# Initialize edge colors as black (default)
edge_colors <- rep("black", nrow(treeMES$edge))

# Function to get all descendant edges of a given node
get_descendants <- function(tree, node) {
  descendants <- c()
  to_visit <- node
  while (length(to_visit) > 0) {
    current <- to_visit[1]
    to_visit <- to_visit[-1]
    
    child_edges <- which(tree$edge[, 1] == current)  # Find edges where current is the parent
    if (length(child_edges) > 0) {
      child_nodes <- tree$edge[child_edges, 2]  # Get descendant nodes
      descendants <- c(descendants, child_edges)  # Store descendant edges
      to_visit <- c(to_visit, child_nodes)  # Add to visit queue
    }
  }
  return(descendants)
}

# Identify edges that meet the conditions
red_edges <- which(outputMES$P...0.1 > 0.8)  # Edges where P...0.1 > 0.8
blue_edges <- which(outputMES$P...0 > 0.8)   # Edges where P...0 > 0.8

# Function to propagate color without overwriting an already colored edge
propagate_color <- function(tree, edges, color) {
  for (edge in edges) {
    descendants <- get_descendants(tree, tree$edge[edge, 2])
    for (desc in descendants) {
      if (edge_colors[desc] == "black") {  # Only change if it's still black
        edge_colors[desc] <- color
      }
    }
  }
}

# Assign initial colors
edge_colors[red_edges] <- "orange2"
edge_colors[blue_edges] <- "blue"

# First propagation without overwriting
propagate_color(treeMES, red_edges, "orange2")
propagate_color(treeMES, blue_edges, "blue")

# Function to color remaining black edges based on closest ancestor's color
color_remaining_black_edges <- function(tree, edge_colors) {
  for (i in seq_along(edge_colors)) {
    if (edge_colors[i] == "black") {  # Only process black edges
      parent <- tree$edge[i, 1]  # Find parent node
      while (parent > 0) {  # Traverse up until we find a colored ancestor
        ancestor_edge <- which(tree$edge[, 2] == parent)
        if (length(ancestor_edge) > 0) {
          ancestor_color <- edge_colors[ancestor_edge]
          if (ancestor_color != "black") {
            edge_colors[i] <- ancestor_color  # Assign ancestor's color
            break
          }
          parent <- tree$edge[ancestor_edge, 1]  # Move up to next ancestor
        } else {
          break  # Stop if no ancestor found
        }
      }
    }
  }
  return(edge_colors)
}

# Apply the function to color remaining black edges
edge_colors <- color_remaining_black_edges(treeMES, edge_colors)

# Plot the phylogenetic tree as a fan
par(mar = c(0.1,0.1,0.1,0.1))
plot(treeMES, type = "fan", show.tip.label = T, edge.width = 2.5, cex = 0.5 , edge.color = edge_colors)

# Add red and blue dots at branch midpoints
nodelabels(pch = 16, col = "orange", cex = 1.2, node = red_branchesMES)
nodelabels(pch = 16, col = "blue", cex = 1.2, node = blue_branchesMES)

#####  Rotated axis #### 

#load("G:/macroevo/lacertids.Rdata")
logR <- bt_read.log("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/rotated-morphospace-green-brown.txt.Log.txt")
logR.mcmc <- as.mcmc(logR)
plot(logR.mcmc)

outputR <- read.csv("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/rotated-morphospace-green-brown.txt.VarRates.txt.csv", header=TRUE)
treeR <- read.nexus("E:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/lacertid-nigri.nex") 
outputR$ID <- 1:349  # Shift IDs to match R's 1-based index
outputR <- outputR[-1,]

# Identify branches meeting the conditions
red_branchesR <- treeR$edge[outputR$P...1.1 > 0.80, 2]  # Extract node numbers for red dots
blue_branchesR <- treeR$edge[outputR$P...1 > 0.80, 2]   # Extract node numbers for blue dots

# Initialize edge colors as black (default)
edge_colors <- rep("black", nrow(treeR$edge))

# Function to get all descendant edges of a given node
get_descendants <- function(tree, node) {
  descendants <- c()
  to_visit <- node
  while (length(to_visit) > 0) {
    current <- to_visit[1]
    to_visit <- to_visit[-1]
    
    child_edges <- which(tree$edge[, 1] == current)  # Find edges where current is the parent
    if (length(child_edges) > 0) {
      child_nodes <- tree$edge[child_edges, 2]  # Get descendant nodes
      descendants <- c(descendants, child_edges)  # Store descendant edges
      to_visit <- c(to_visit, child_nodes)  # Add to visit queue
    }
  }
  return(descendants)
}

# Identify edges that meet the conditions
red_edges <- which(outputR$P...0.1 > 0.8)  # Edges where P...0.1 > 0.8
blue_edges <- which(outputR$P...0 > 0.8)   # Edges where P...0 > 0.8

# Function to propagate color without overwriting an already colored edge
propagate_color <- function(tree, edges, color) {
  for (edge in edges) {
    descendants <- get_descendants(tree, tree$edge[edge, 2])
    for (desc in descendants) {
      if (edge_colors[desc] == "black") {  # Only change if it's still black
        edge_colors[desc] <- color
      }
    }
  }
}

# Assign initial colors
edge_colors[red_edges] <- "orange"
edge_colors[blue_edges] <- "blue"

# First propagation without overwriting
propagate_color(treeR, red_edges, "orange")
propagate_color(treeR, blue_edges, "blue")

# Function to color remaining black edges based on closest ancestor's color
color_remaining_black_edges <- function(tree, edge_colors) {
  for (i in seq_along(edge_colors)) {
    if (edge_colors[i] == "black") {  # Only process black edges
      parent <- tree$edge[i, 1]  # Find parent node
      while (parent > 0) {  # Traverse up until we find a colored ancestor
        ancestor_edge <- which(tree$edge[, 2] == parent)
        if (length(ancestor_edge) > 0) {
          ancestor_color <- edge_colors[ancestor_edge]
          if (ancestor_color != "black") {
            edge_colors[i] <- ancestor_color  # Assign ancestor's color
            break
          }
          parent <- tree$edge[ancestor_edge, 1]  # Move up to next ancestor
        } else {
          break  # Stop if no ancestor found
        }
      }
    }
  }
  return(edge_colors)
}

# Apply the function to color remaining black edges
edge_colors <- color_remaining_black_edges(treeR, edge_colors)
# Plot the phylogenetic tree as a fan
par(mar = c(0.1,0.1,0.1,0.1))
plot(treeR, type = "fan", show.tip.label = T, edge.width = 2.5, cex = 0.5 , edge.color = edge_colors)

# Add red and blue dots at branch midpoints
nodelabels(pch = 16, col = "red", cex = 1.2, node = red_branchesR)
nodelabels(pch = 16, col = "blue", cex = 1.2, node = blue_branchesR)




log <- bt_read.log("G:/macroevo/BayesTraitsV4.0.0-Win64/BayesTraitsV4.0.0-Win64/PCAscores.txt.Log.txt")


