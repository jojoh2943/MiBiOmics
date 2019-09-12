#-------------------------#
#### GENERAL FUNCTIONS ####
#-------------------------#

'%!in%' <- function(x,y)!('%in%'(x,y))

rename_axis_drivers_simplified <- function(axis_drivers){
  right_names <- substr(rownames(axis_drivers),3, nchar(rownames(axis_drivers)))
  return(right_names)
}

rename_axis_drivers <- function(axis_drivers, colnames_df1, colnames_df2, colnames_df3, n_df =3){
  right_names <- c()
  if (n_df == 3){
    for (row in 1:nrow(axis_drivers)){
      if (substr(rownames(axis_drivers)[row], nchar(rownames(axis_drivers)[row])- 3, nchar(rownames(axis_drivers)[row])) == ".df1"){
        if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
          if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df1){
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
          }else{
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
          }
        }else{
          if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %in% colnames_df1){
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
          }else{
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
          }
          
        }
      }else{
        if (substr(rownames(axis_drivers)[row], nchar(rownames(axis_drivers)[row])- 3, nchar(rownames(axis_drivers)[row])) == ".df2"){
          if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df2){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
            }
          }else{
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df2){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
            }
          }
        }else{
          if (substr(rownames(axis_drivers)[row], nchar(rownames(axis_drivers)[row])- 3, nchar(rownames(axis_drivers)[row])) == ".df3"){
            if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
              if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df3){
                right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
              }else{
                right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
              }
            }else{
              if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df3){
                right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
              }else{
                right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
              }            }
          }
        }
      }
    }
  }else{
    for (row in 1:nrow(axis_drivers)){
      if (substr(rownames(axis_drivers)[row], nchar(rownames(axis_drivers)[row])- 3, nchar(rownames(axis_drivers)[row])) == ".df1"){
        if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
          if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df1){
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
          }else{
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
          }
        }else{
          if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df1){
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
          }else{
            right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
          }        
          }
      }else{
        if (substr(rownames(axis_drivers)[row], nchar(rownames(axis_drivers)[row])- 3, nchar(rownames(axis_drivers)[row])) == ".df2"){
          if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df2){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
            }
          }else{
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4) %!in% colnames_df1){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])- 4))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])- 4))
            }          
            }
        }else{
          if (substr(rownames(axis_drivers)[row], 0, 1)== "X"){
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])) %!in% colnames_df2){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])))
            }
          }else{
            if (substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])) %!in% colnames_df1){
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 2, nchar(rownames(axis_drivers)[row])))
            }else{
              right_names <- c(right_names, substr(rownames(axis_drivers)[row], 0, nchar(rownames(axis_drivers)[row])))
            }          
          }
        }
      }
    }
  }
  
  return(right_names)
}


TSS.divide = function(x){
  x/sum(x)
}

'%!in%' <- function(x,y)!('%in%'(x,y)) 

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


## Specific functions: PLS and VIP Scores
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}
#Function for hive plot
vector_geom_curve <- function(edge_df, df_ggplot){
  vector_x <- c()
  vector_y <- c()
  label_fromNode <- c()
  vector_xend <- c()
  vector_yend <-c()
  label_toNode <- c()
  for (i in 1:nrow(edge_df)){
    x <- df_ggplot$x1[which(df_ggplot$label == edge_df$fromNode[i])]
    label_fromNode_ind <- edge_df$fromNode[i]
    y = 0
    xend = 0
    yend <-df_ggplot$y2[which(df_ggplot$label == edge_df$toNode[i])]
    label_toNode_ind <- edge_df$toNode[i]
    vector_x <- c(vector_x, x)
    vector_y <- c(vector_y, y)
    label_fromNode <- c(label_fromNode, label_fromNode_ind)
    label_toNode <- c(label_toNode, label_toNode_ind)
    vector_xend <- c(vector_xend, xend)
    vector_yend <-c(vector_yend, yend)
  }
  df_geom_curve <- data.frame(label_fromNode, vector_x, vector_y, label_toNode, vector_xend, vector_yend)
  return(df_geom_curve)
}

#### HIVE FUNCTIONS ####

hive_myLayers <- function(l_WGCNA_D1, l_WGCNA_D2, l_WGCNA_D3, myCorrs, myAnnots, correlation= "spearman", trait, exportPDF = FALSE, exportSVG = FALSE, sizePlot = 10, nameFile = "hive_D1_D2_D3", cureCorr = FALSE, exportCSV = FALSE){
  require(leaflet)
  
  annot_D1 <- myAnnots[[1]]
  annot_D2 <- myAnnots[[2]]
  annot_D3 <- myAnnots[[3]]
  WGCNA_list <- list()
  WGCNA_list[[1]] <- l_WGCNA_D1
  WGCNA_list[[2]] <- l_WGCNA_D2
  WGCNA_list[[3]] <- l_WGCNA_D3
  pal <- colorNumeric(palette = "RdBu", 1:-1)
  
  
  # Create nodes dataframe:
  id <- c(seq(from = 1, to = (ncol(l_WGCNA_D1[[5]])+ ncol(l_WGCNA_D2[[5]]) + ncol(l_WGCNA_D3[[5]]) +6 )))
  label <- c("min_DF1", "min_DF2", "min_DF3", as.vector(paste(colnames(l_WGCNA_D1[[5]]), "DF1", sep = "_")), as.vector(paste( colnames(l_WGCNA_D2[[5]]), "DF2", sep = "_")), as.vector(paste( colnames(l_WGCNA_D3[[5]]), "DF3", sep = "_")), "extreme_DF1", "extreme_DF2", "extreme_DF3" )
  color <- c("white", "white", "white", as.vector(substr(colnames(l_WGCNA_D1[[5]]), 3, 30)), as.vector(substr( colnames(l_WGCNA_D2[[5]]), 3, 30)), as.vector(substr(colnames(l_WGCNA_D3[[5]]), 3, 30)), "white", "white", "white" )
  axis <- c(1, 2, 3, as.vector(rep(1, ncol(l_WGCNA_D1[[5]]))), as.vector(rep(2, ncol(l_WGCNA_D2[[5]]))), as.vector(rep(3, ncol(l_WGCNA_D3[[5]]))), 1, 2, 3 )
  size <- c(0, 0, 0, as.vector(rep(1, sum(ncol(l_WGCNA_D1[[5]]), ncol(l_WGCNA_D2[[5]]), ncol(l_WGCNA_D3[[5]])))), 0, 0, 0)
  radius <- c()
  v_pv <- c()
  for (i in 1:3){
    annot <- myAnnots[[i]]
    myWGCNA <- WGCNA_list[[i]]
    if (is.numeric(annot[,trait])){
      moduleTraitCor = cor(myWGCNA[[5]], annot[,trait], use = "p", method = correlation)
    }else{
      annot2<- annot
      annot2[,trait] <- as.factor(annot2[, trait])
      annot2[,trait] <- as.numeric(annot2[,trait])
      moduleTraitCor = cor(myWGCNA[[5]], annot2[,trait], use = "p", method = correlation)
    }
    pv <- corPvalueStudent(moduleTraitCor, nrow(annot))
    
    radius <- c(radius, abs(moduleTraitCor)*100)
    v_pv <- c(v_pv, pv)
  }
  v_pv <- c(0, 0, 0 ,v_pv, 0, 0, 0)
  df_size <- data.frame("size" = size, "pv" = v_pv)
  df_size[df_size$pv > 0.05,] <- 0
  radius <- c(0, 0, 0, radius, 100, 100, 100)
  nodes <- data.frame(id = id, lab= as.character(label), axis = as.integer(axis), radius= radius, size= df_size$size, color= as.character(color))
  # Do we set the size of the node to zero if the p-value is not significative ?
  
  # create edge dataframe:
  id1 <- c()
  id2 <- c()
  weight <- c()
  color <- c()
  for (k in 1:3){
    myCorr <- myCorrs[[k]]
    for (i in 1:nrow(myCorr[[2]])){
      for (j in 1:ncol(myCorr[[2]])){
        if (myCorr[[2]][i, j] < 0.05 && abs(myCorr[[1]][i, j]) > 0.3){
          myNode1 <-nodes[which(nodes$lab==rownames(myCorr[[2]])[i]),]
          myNode2 <-nodes[which(nodes$lab==colnames(myCorr[[2]])[j]),]
          id1 <- c(id1, myNode1$id)
          id2 <- c(id2, myNode2$id)
          weight <- c(weight, round(exp(abs(myCorr[[1]][i, j]))^3))
          color <- c(color, pal(c(myCorr[[1]][i, j])))
        }
      }
    }
  }
  if (length(id1) == 0){
    edges <- 0
  }else{
    edges <- data.frame(id1 = id1, id2 = id2, weight = weight, color = as.character(color))
    id_to_keep <- nodes[which(nodes$size != 0),]
    id_to_keep <- id_to_keep$id
    id_to_keep <- c(1, 2, 3, id_to_keep, nodes$id[(length(nodes$id)-2):length(nodes$id)])
    nodes <- nodes[which(nodes$id %in% id_to_keep),]
    
    
    edges <- edges[which(edges$id1 %in% id_to_keep),]
    edges <- edges[which(edges$id2 %in% id_to_keep),]

    if (cureCorr){
      nodes <- nodes[which(nodes$radius > 40),]
      edges <- edges[which(edges$id1 %in% nodes$id),]
      edges <- edges[which(edges$id2 %in% nodes$id),]
    }
    type <- "2D"
    desc <- "Hive Plot 3 Layers"
    axis.cols <- c("#636363", "#636363", "#636363")
    if (exportCSV){
      write.csv(nodes, file = paste("nodes_", trait, "_", nameFile, ".csv", sep = ""))
      write.csv(edges, file = paste("edges_", trait, "_", nameFile, ".csv", sep = ""))
      
    }
  }
  myHive <- list(nodes = nodes, edges = edges, type = type, desc = desc, axis.cols= axis.cols)
  
  return(myHive)
}




hive_my2Layers <- function(l_WGCNA_D1, l_WGCNA_D2, myCorr, myAnnots, correlation= "spearman", trait, exportPDF = FALSE, exportSVG = FALSE, sizePlot = 10, nameFile = "hive_D1_D2_D3", cureCorr = FALSE, exportCSV = FALSE){
  require(leaflet)
  
  annot_D1 <- myAnnots[[1]]
  annot_D2 <- myAnnots[[2]]
  WGCNA_list <- list()
  WGCNA_list[[1]] <- l_WGCNA_D1 
  WGCNA_list[[2]] <- l_WGCNA_D2
  pal <- colorNumeric(palette = "RdBu", 1:-1)
  
  
  # Create nodes dataframe:
  id <- c(seq(from = 1, to = (ncol(l_WGCNA_D1[[5]])+ ncol(l_WGCNA_D2[[5]]) + 4 )))
  label <- c("min_DF1", "min_DF2", as.vector(paste(colnames(l_WGCNA_D1[[5]]), "DF1", sep = "_")), as.vector(paste( colnames(l_WGCNA_D2[[5]]), "DF2", sep = "_")), "extreme_DF1", "extreme_DF2")
  color <- c("white", "white", as.vector(substr(colnames(l_WGCNA_D1[[5]]), 3, 30)), as.vector(substr( colnames(l_WGCNA_D2[[5]]), 3, 30)), "white", "white")
  axis <- c(1, 2, as.vector(rep(1, ncol(l_WGCNA_D1[[5]]))), as.vector(rep(2, ncol(l_WGCNA_D2[[5]]))), 1, 2 )
  size <- c(0, 0, as.vector(rep(1, sum(ncol(l_WGCNA_D1[[5]]), ncol(l_WGCNA_D2[[5]])))), 0, 0)
  radius <- c()
  v_pv <- c()
  for (i in 1:2){
    annot <- myAnnots[[i]]
    myWGCNA <- WGCNA_list[[i]]

    if (is.numeric(annot[,trait])){

      moduleTraitCor = cor(myWGCNA[[5]], annot[,trait], use = "p", method = correlation)
    }else{
      annot2<- annot
      annot2[,trait] <- as.factor(annot2[, trait])
      annot2[,trait] <- as.numeric(annot2[,trait])
      moduleTraitCor = cor(myWGCNA[[5]], annot2[,trait], use = "p", method = correlation)
    }
    pv <- corPvalueStudent(moduleTraitCor, nrow(annot))
    
    radius <- c(radius, abs(moduleTraitCor)*100)
    v_pv <- c(v_pv, pv)
  }
  v_pv <- c(0, 0,v_pv, 0, 0)
  df_size <- data.frame("size" = size, "pv" = v_pv)
  df_size[df_size$pv > 0.05,] <- 0
  radius <- c(0, 0, radius, 100, 100)
  nodes <- data.frame(id = id, lab= as.character(label), axis = as.integer(axis), radius= radius, size= df_size$size, color= as.character(color))
  # Do we set the size of the node to zero if the p-value is not significative ?
  
  # create edge dataframe:
  id1 <- c()
  id2 <- c()
  weight <- c()
  color <- c()
  
  for (i in 1:nrow(myCorr[[2]])){
    for (j in 1:ncol(myCorr[[2]])){
      if (myCorr[[2]][i, j] < 0.05 && abs(myCorr[[1]][i, j]) > 0.4){
        myNode1 <-nodes[which(nodes$lab==rownames(myCorr[[2]])[i]),]
        myNode2 <-nodes[which(nodes$lab==colnames(myCorr[[2]])[j]),]
        id1 <- c(id1, myNode1$id)
        id2 <- c(id2, myNode2$id)
        weight <- c(weight, round(exp(abs(myCorr[[1]][i, j]))^3))
        color <- c(color, pal(c(myCorr[[1]][i, j])))
      }
    }
  }
  if (length(id1) == 0){
    myHive <- 0
  }else{
    edges <- data.frame(id1 = id1, id2 = id2, weight = weight, color = as.character(color))
    id_to_keep <- nodes[which(nodes$size != 0),]
    id_to_keep <- id_to_keep$id
    id_to_keep <- c(1, 2, id_to_keep, nodes$id[(length(nodes$id)-1):length(nodes$id)])
    nodes <- nodes[which(nodes$id %in% id_to_keep),]
    edges <- edges[which(edges$id1 %in% id_to_keep),]
    if (nrow(edges) != 0){
      edges <- edges[which(edges$id2 %in% id_to_keep),]
    }
    
    
    if (cureCorr){
      nodes <- nodes[which(nodes$radius > 40),]
      edges <- edges[which(edges$id1 %in% nodes$id),]
      edges <- edges[which(edges$id2 %in% nodes$id),]
    }
    
    if (exportCSV){
      write.csv(nodes, file = paste("nodes_", trait, "_", nameFile, ".csv", sep = ""))
      write.csv(edges, file = paste("edges_", trait, "_", nameFile, ".csv", sep = ""))
      
    }
    type <- "2D"
    desc <- "Hive Plot 2 Layers"
    axis.cols <- c("#636363", "#636363")
    myHive <- list(nodes = nodes, edges = edges, type = type, desc = desc, axis.cols= axis.cols)
  }
  
  return(myHive)
  
}




plotMyHive <- function(HPD, ch = 1, method = "abs",
                       dr.nodes = TRUE, bkgnd = "black",
                       axLabs = NULL, axLab.pos = NULL, axLab.gpar = NULL,
                       anNodes = NULL, anNode.gpar = NULL, grInfo = NULL,
                       arrow = NULL, np = TRUE, anCoord = "local", ...) {
  
  # Function to plot hive plots using grid graphics
  # Inspired by the work of Martin Kryzwinski
  # Bryan Hanson, DePauw Univ, Feb 2011 onward
  
  # This function is intended to draw in 2D for nx from 2 to 6
  # The results will be similar to the original hive plot concept
  
  ##### Set up some common parameters
  
  if (!HPD$type == "2D") stop("This is not a 2D hive data set: use plot3dHive instead")		
  # chkHPD(HPD)
  nx <- length(unique(HPD$nodes$axis))
  
  if (nx == 1) stop("Something is wrong: only one axis seems to be present")
  
  # Send out for ranking/norming/pruning/inverting if requested
  
  if (!method == "abs") HPD <- manipAxis(HPD, method, ...)
  
  nodes <- HPD$nodes
  edges <- HPD$edges
  axis.cols <- HPD$axis.cols
  
  # Fix up center hole
  
  nodes$radius <- nodes$radius + ch
  HPD$nodes$radius <- nodes$radius
  
  ##### Some convenience functions, only defined in this function environ.
  ##### The two long functions need to stay here for simplicity, since
  ##### all of the radius checking etc is here and if moved elsewhere,
  ##### these calculations would have to be redone or results passed.
  
  p2cX <- function(r, theta) { x <- r*cos(theta*2*pi/360) }
  p2cY <- function(r, theta) { y <- r*sin(theta*2*pi/360) }
  
  addArrow <- function(arrow, nx) {
    if (!length(arrow) >= 5) stop("Too few arrow components")
    if (is.null(axLab.gpar)) {
      if (bkgnd == "black") axLab.gpar <- gpar(fontsize = 12, col = "white", lwd = 2)
      if (!bkgnd == "black") axLab.gpar <- gpar(fontsize = 12, col = "black", lwd = 2)
    }
    a <- as.numeric(arrow[2])
    rs <- as.numeric(arrow[3])
    re <- as.numeric(arrow[4])
    b <- as.numeric(arrow[5]) # label offset from end of arrow
    
    x.st <- p2cX(rs, a)
    y.st <- p2cY(rs, a)
    x.end <- p2cX(re, a)
    y.end <- p2cY(re, a)
    
    x.lab <- p2cX(re + b, a) # figure arrow label position
    y.lab <- p2cY(re + b, a)
    al <- 0.2*(re-rs) # arrow head length
    
    # for nx = 2 only, offset the arrow
    # in the y direction to save space overall
    
    if (nx == 2) {
      if (is.na(arrow[6])) {
        arrow[6] <- 0
        cat("\tThe arrow can be offset vertically; see ?plotHive\n")
      }
      y.st <- y.st + as.numeric(arrow[6])
      y.end <- y.end + as.numeric(arrow[6])
      y.lab <- y.lab + as.numeric(arrow[6])		
    }
    
    grid.lines(x = c(x.st, x.end), y = c(y.st, y.end),
               arrow = arrow(length = unit(al, "native")),
               default.units = "native", gp = axLab.gpar)
    grid.text(arrow[1], x.lab, y.lab, default.units = "native", gp = axLab.gpar)
  }
  
  annotateNodes <- function(anNodes, nodes, nx, anCoord) {
    
    if (is.null(anNode.gpar)) {
      if (bkgnd == "black") anNode.gpar <- gpar(fontsize = 10, col = "white", lwd = 0.5)
      if (!bkgnd == "black") anNode.gpar <- gpar(fontsize = 10, col = "black", lwd = 0.5)
    }
    
    ann <- utils::read.csv(anNodes, header = TRUE, colClasses = c(rep("character", 2), rep("numeric", 5)))
    cds <- getCoords(anNodes, anCoord, nodes)
    
    grid.segments(x0 = cds$x.st, x1 = cds$x.end, y0 = cds$y.st, y1 = cds$y.end,
                  default.units = "native", gp = anNode.gpar)
    grid.text(ann$node.text, cds$x.lab, cds$y.lab, hjust = ann$hjust, vjust = ann$vjust,
              default.units = "native", gp = anNode.gpar, ...)
  }
  
  addGraphic <- function(grInfo, nodes, nx, anCoord) {
    
    gr <- utils::read.csv(grInfo, header = TRUE, stringsAsFactors = FALSE)
    cds <- getCoords(grInfo, anCoord, nodes)
    
    grid.segments(x0 = cds$x.st, x1 = cds$x.end, y0 = cds$y.st, y1 = cds$y.end,
                  default.units = "native", gp = anNode.gpar)
    
    # readJPEG and readPNG are not vectorized, grab each graphic in turn
    # Figure out if we are using jpg or png files
    
    ext <- substr(gr$path[1], nchar(gr$path[1])-2, nchar(gr$path[1]))
    if ((ext == "png") | (ext == "PNG")) ext <- "png"
    if ((ext == "jpg") | (ext == "JPG") | (ext == "peg") | (ext =="PEG")) ext <- "jpg"
    
    # Now draw the images
    
    if (ext == "jpg") {
      for (n in 1:nrow(gr)) {
        grid.raster(readJPEG(gr$path[n]),
                    x = cds$x.lab[n], y = cds$y.lab[n], default.units = "native", width = gr$width[n])			
      }
    }
    
    if (ext == "png") {
      for (n in 1:nrow(gr)) {
        grid.raster(readPNG(gr$path[n]),
                    x = cds$x.lab[n], y = cds$y.lab[n], default.units = "native", width = gr$width[n])			
      }
    }
    
  }
  
  getCoords <- function(file, anCoord, nodes) {
    
    # Figure out the coordinates of the line segments and labels/graphics
    # anNodes and grInfo both contains certain columns which are used here
    
    df <- utils::read.csv(file, header = TRUE)
    
    id <- rep(NA, nrow(df))	
    for (n in 1:nrow(df)) {			
      pat <- paste("\\b", df$node.lab[n], "\\b", sep = "") 
      id[n] <- grep(pat, nodes$lab)
    }
    
    N <- matrix(data = c(
      0, 180, NA, NA, NA, NA,
      90, 210, 330, NA, NA, NA,
      90, 180, 270, 0, NA, NA,
      90, 162, 234, 306, 18, NA,
      90, 150, 210, 270, 330, 390),
      byrow = TRUE, nrow = 5)
    
    ax <- nodes$axis[id] # axis number
    for (n in 1:length(ax)) {
      ax[n] <- N[nx-1,ax[n]]			
    }
    
    # Figure coords in requested reference frame
    
    x.st <- p2cX(nodes$radius[id], ax)
    y.st <- p2cY(nodes$radius[id], ax)
    
    if (anCoord == "local") {
      x.end <- x.st + p2cX(df$radius, df$angle)
      y.end <- y.st + p2cY(df$radius, df$angle)
      
      x.lab <- x.st + p2cX(df$radius + df$offset, df$angle)
      y.lab <- y.st + p2cY(df$radius + df$offset, df$angle)		
    }
    
    if (anCoord == "global") {					
      x.end <- p2cX(df$radius, df$angle)
      y.end <- p2cY(df$radius, df$angle)
      
      x.lab <- p2cX(df$radius + df$offset, df$angle)
      y.lab <- p2cY(df$radius + df$offset, df$angle)		
    }
    
    retval <- data.frame(x.st, y.st, x.end, y.end, x.lab, y.lab)
    retval
  }
  
  ###############
  
  # Figure out which nodes to draw for each edge
  # Since they are in random order
  # Do this once/early to save time
  
  id1 <- id2 <- c()
  
  for (n in 1:nrow(edges)) {
    
    pat1 <- paste("\\b", edges$id1[n], "\\b", sep = "") 
    pat2 <- paste("\\b", edges$id2[n], "\\b", sep = "")
    id1 <- c(id1, grep(pat1, nodes$id))
    id2 <- c(id2, grep(pat2, nodes$id))
    
  }
  
  ##### Two dimensional case  (using grid graphics)
  
  # Prep axes first
  
  if (nx == 2) {
    
    # n1 <- subset(nodes, axis == 1)
    # n2 <- subset(nodes, axis == 2)
    
    n1 <- nodes[nodes[,"axis"] == 1,]
    n2 <- nodes[nodes[,"axis"] == 2,]
    
    max1 <- max(n1$radius)
    max2 <- max(n2$radius)
    min1 <- min(n1$radius)
    min2 <- min(n2$radius)
    
    r.st <- c(min1, min2) # in polar coordinates
    axst <- c(0, 180)
    x0a = p2cX(r.st, axst)
    y0a = p2cY(r.st, axst)
    
    r.end <- c(max1, max2)
    axend <- c(0, 180)
    x1a = p2cX(r.end, axend)
    y1a = p2cY(r.end, axend)
    
    # Set up grid graphics viewport
    
    md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.5 # max dimension
    # 1.5 is used in case of labels
    
    if (np) grid.newpage()
    grid.rect(gp = gpar(col = NA, fill = bkgnd))
    vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                   xscale = c(-md, md), yscale = c(-md, md),
                   name = "3DHivePlot")
    
    pushViewport(vp)
    
    # Now draw edges
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if (nodes$axis[id1[n]] == 1) { # set up edge start params 1st
        th.st <- c(th.st, 0)
        r.st <- c(r.st, nodes$radius[id1[n]])
      }
      if (nodes$axis[id1[n]] == 2) {
        th.st <- c(th.st, 180)
        r.st <- c(r.st, nodes$radius[id1[n]])
      }
      
      if (nodes$axis[id2[n]] == 1) { # now edge end params
        th.end <- c(th.end, 0)
        r.end <- c(r.end, nodes$radius[id2[n]])
      }
      if (nodes$axis[id2[n]] == 2) {
        th.end <- c(th.end, 180)
        r.end <- c(r.end, nodes$radius[id2[n]])
      }
      
      ecol <- c(ecol, edges$color[n])
      ewt <- c(ewt, edges$weight[n])
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Draw axes
    
    grid.segments(x0a, y0a, x1a, y1a,
                  gp = gpar(col = HPD$axis.cols, lwd = 8),
                  default.units = "native")
    
    # Now add nodes
    
    if (dr.nodes) {
      r <- c(n1$radius, n2$radius) 
      theta <- c(rep(0, length(n1$radius)),
                 rep(180, length(n2$radius)))
      x = p2cX(r, theta)
      y = p2cY(r, theta)
      grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size), col = c(n1$color, n2$color)))
    }
    
    # Now label axes
    
    if (!is.null(axLabs)) {
      if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
      if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
      r <- c(max1, max2)
      if (is.null(axLab.pos)) axLab.pos <- r*0.1
      r <- r + axLab.pos
      th <- c(0, 180)
      x <- p2cX(r, th)
      y <- p2cY(r, th)
      grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
    }
    
    # Add a legend arrow & any annotations
    
    if (!is.null(arrow)) addArrow(arrow, nx)
    if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx, anCoord)
    if (!is.null(grInfo)) addGraphic(grInfo, nodes, nx, anCoord)
    
  } # end of 2D
  
  ##### Three dimensional case (using grid graphics)
  
  # Prep axes first
  
  if (nx == 3) {
    
    # n1 <- subset(nodes, axis == 1)
    # n2 <- subset(nodes, axis == 2)
    # n3 <- subset(nodes, axis == 3)
    
    n1 <- nodes[nodes[,"axis"] == 1,]
    n2 <- nodes[nodes[,"axis"] == 2,]
    n3 <- nodes[nodes[,"axis"] == 3,]
    
    max1 <- max(n1$radius)
    max2 <- max(n2$radius)
    max3 <- max(n3$radius)
    min1 <- min(n1$radius)
    min2 <- min(n2$radius)
    min3 <- min(n3$radius)
    
    r.st <- c(min1, min2, min3) # in polar coordinates
    axst <- c(90, 210, 330)
    x0a = p2cX(r.st, axst)
    y0a = p2cY(r.st, axst)
    
    r.end <- c(max1, max2, max3)
    axend <- c(90, 210, 330)
    x1a = p2cX(r.end, axend)
    y1a = p2cY(r.end, axend)
    
    # Set up grid graphics viewport
    
    md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
    if (np) grid.newpage()
    grid.rect(gp = gpar(col = NA, fill = bkgnd))
    vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                   xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
    pushViewport(vp)
    
    # Now draw edges (must do in sets as curvature is not vectorized)
    
    # Axis 1 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 2 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 3 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 1 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 3 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 2 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Draw axes
    
    grid.segments(x0a, y0a, x1a, y1a,
                  gp = gpar(col = HPD$axis.cols, lwd = 3),
                  default.units = "native")
    
    # Now add nodes
    
    if (dr.nodes) {
      r <- c(n1$radius, n2$radius, n3$radius) 
      theta <- c(rep(90, length(n1$radius)),
                 rep(210, length(n2$radius)),
                 rep(330, length(n3$radius)))
      x = p2cX(r, theta)
      y = p2cY(r, theta)
      grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size),
                                            col = c(n1$color, n2$color, n3$color)))
    }
    
    # Now label axes
    
    if (!is.null(axLabs)) {
      if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
      if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
      r <- c(max1, max2, max3)
      if (is.null(axLab.pos)) axLab.pos <- r*0.1
      r <- r + axLab.pos
      th <- c(90, 210, 330)
      x <- p2cX(r, th)
      y <- p2cY(r, th)
      grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
    }
    
    # Add a legend arrow & any annotations
    
    if (!is.null(arrow)) addArrow(arrow, nx)
    if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx, anCoord)
    if (!is.null(grInfo)) addGraphic(grInfo, nodes, nx, anCoord)
    
  } # end of 3D
  
  
  ##### Four dimensional case (using grid graphics)
  
  # Prep axes first
  
  if (nx == 4) {
    
    # n1 <- subset(nodes, axis == 1)
    # n2 <- subset(nodes, axis == 2)
    # n3 <- subset(nodes, axis == 3)
    # n4 <- subset(nodes, axis == 4)
    
    n1 <- nodes[nodes[,"axis"] == 1,]
    n2 <- nodes[nodes[,"axis"] == 2,]
    n3 <- nodes[nodes[,"axis"] == 3,]
    n4 <- nodes[nodes[,"axis"] == 4,]
    
    max1 <- max(n1$radius)
    max2 <- max(n2$radius)
    max3 <- max(n3$radius)
    max4 <- max(n4$radius)
    min1 <- min(n1$radius)
    min2 <- min(n2$radius)
    min3 <- min(n3$radius)
    min4 <- min(n4$radius)
    
    r.st <- c(min1, min2, min3, min4) # in polar coordinates
    axst <- c(90, 180, 270, 0)
    x0a = p2cX(r.st, axst)
    y0a = p2cY(r.st, axst)
    
    r.end <- c(max1, max2, max3, max4)
    axend <- c(90, 180, 270, 0)
    x1a = p2cX(r.end, axend)
    y1a = p2cY(r.end, axend)
    
    # Set up grid graphics viewport
    
    md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.5 # max dimension
    if (np) grid.newpage()
    grid.rect(gp = gpar(col = NA, fill = bkgnd))
    vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                   xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
    pushViewport(vp)
    
    # Now draw edges (must do in sets as curvature is not vectorized)
    
    # Axis 1 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 180)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 2 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 180)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 3 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 0)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 4 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 0)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 1 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 0)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 4 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 0)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 3 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 180)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }				
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 2 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 180)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }				
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 180)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 180)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 0)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 0)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Draw axes
    
    grid.segments(x0a, y0a, x1a, y1a,
                  gp = gpar(col = HPD$axis.cols, lwd = 3),
                  default.units = "native")
    
    # Now add nodes
    
    if (dr.nodes) {
      r <- c(n1$radius, n2$radius, n3$radius, n4$radius) 
      theta <- c(rep(90, length(n1$radius)),
                 rep(180, length(n2$radius)),
                 rep(270, length(n3$radius)),
                 rep(0, length(n4$radius)))
      x = p2cX(r, theta)
      y = p2cY(r, theta)
      grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size),
                                            col = c(n1$color, n2$color, n3$color, n4$color)))
    }
    
    # Now label axes
    
    if (!is.null(axLabs)) {
      if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
      if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
      r <- c(max1, max2, max3, max4)
      if (is.null(axLab.pos)) axLab.pos <- r*0.1
      r <- r + axLab.pos
      th <- c(90, 180, 270, 0)
      x <- p2cX(r, th)
      y <- p2cY(r, th)
      grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
    }
    
    # Add a legend arrow & any annotations
    
    if (!is.null(arrow)) addArrow(arrow, nx)
    if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx, anCoord)
    if (!is.null(grInfo)) addGraphic(grInfo, nodes, nx, anCoord)
    
  } # end of 4D
  
  ##### Five dimensional case (using grid graphics)
  
  # Prep axes first
  
  if (nx == 5) {
    
    # n1 <- subset(nodes, axis == 1)
    # n2 <- subset(nodes, axis == 2)
    # n3 <- subset(nodes, axis == 3)
    # n4 <- subset(nodes, axis == 4)
    # n5 <- subset(nodes, axis == 5)
    
    n1 <- nodes[nodes[,"axis"] == 1,]
    n2 <- nodes[nodes[,"axis"] == 2,]
    n3 <- nodes[nodes[,"axis"] == 3,]
    n4 <- nodes[nodes[,"axis"] == 4,]
    n5 <- nodes[nodes[,"axis"] == 5,]
    
    max1 <- max(n1$radius)
    max2 <- max(n2$radius)
    max3 <- max(n3$radius)
    max4 <- max(n4$radius)
    max5 <- max(n5$radius)
    min1 <- min(n1$radius)
    min2 <- min(n2$radius)
    min3 <- min(n3$radius)
    min4 <- min(n4$radius)
    min5 <- min(n5$radius)
    
    r.st <- c(min1, min2, min3, min4, min5) # in polar coordinates
    axst <- c(90, 162, 234, 306, 18)
    x0a = p2cX(r.st, axst)
    y0a = p2cY(r.st, axst)
    
    r.end <- c(max1, max2, max3, max4, max5)
    axend <- c(90, 162, 234, 306, 18)
    x1a = p2cX(r.end, axend)
    y1a = p2cY(r.end, axend)
    
    # Set up grid graphics viewport
    
    md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
    if (np) grid.newpage()
    grid.rect(gp = gpar(col = NA, fill = bkgnd))
    vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                   xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
    pushViewport(vp)
    
    # Now draw edges (must do in sets as curvature is not vectorized)
    
    # Axis 1 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 162)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }			
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 2 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 162)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 234)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 3 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 234)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 306)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 4 -> 5
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 306)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 18)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 5 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 18)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    
    # Axis 1 -> 5
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 18)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 5 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 18)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 306)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 4 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 306)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 234)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 3 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 234)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 162)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 2 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 162)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 162)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 162)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 234)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 234)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 306)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 306)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 18)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 18)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Draw axes
    
    grid.segments(x0a, y0a, x1a, y1a,
                  gp = gpar(col = HPD$axis.cols, lwd = 3),
                  default.units = "native")
    
    # Now add nodes
    
    if (dr.nodes) {
      r <- c(n1$radius, n2$radius, n3$radius, n4$radius, n5$radius) 
      theta <- c(rep(90, length(n1$radius)),
                 rep(162, length(n2$radius)),
                 rep(234, length(n3$radius)),
                 rep(306, length(n4$radius)),
                 rep(18, length(n5$radius)))
      x = p2cX(r, theta)
      y = p2cY(r, theta)
      grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size, n5$size),
                                            col = c(n1$color, n2$color, n3$color, n4$color, n5$color)))
    }
    
    # Now label axes
    
    if (!is.null(axLabs)) {
      if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
      if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
      r <- c(max1, max2, max3, max4, max5)
      if (is.null(axLab.pos)) axLab.pos <- r*0.1
      r <- r + axLab.pos
      th <- c(90, 162, 234, 306, 18)
      x <- p2cX(r, th)
      y <- p2cY(r, th)
      grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
    }
    
    # Add a legend arrow & any annotations
    
    if (!is.null(arrow)) addArrow(arrow, nx)
    if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx, anCoord)
    if (!is.null(grInfo)) addGraphic(grInfo, nodes, nx, anCoord)
    
  } # end of 5D
  
  ##### Six dimensional case (using grid graphics)
  
  # Prep axes first
  
  if (nx == 6) {
    
    # n1 <- subset(nodes, axis == 1)
    # n2 <- subset(nodes, axis == 2)
    # n3 <- subset(nodes, axis == 3)
    # n4 <- subset(nodes, axis == 4)
    # n5 <- subset(nodes, axis == 5)
    # n6 <- subset(nodes, axis == 6)
    
    n1 <- nodes[nodes[,"axis"] == 1,]
    n2 <- nodes[nodes[,"axis"] == 2,]
    n3 <- nodes[nodes[,"axis"] == 3,]
    n4 <- nodes[nodes[,"axis"] == 4,]
    n5 <- nodes[nodes[,"axis"] == 5,]
    n6 <- nodes[nodes[,"axis"] == 6,]
    
    max1 <- max(n1$radius)
    max2 <- max(n2$radius)
    max3 <- max(n3$radius)
    max4 <- max(n4$radius)
    max5 <- max(n5$radius)
    max6 <- max(n6$radius)
    min1 <- min(n1$radius)
    min2 <- min(n2$radius)
    min3 <- min(n3$radius)
    min4 <- min(n4$radius)
    min5 <- min(n5$radius)
    min6 <- min(n6$radius)
    
    r.st <- c(min1, min2, min3, min4, min5, min6) # in polar coordinates
    axst <- c(90, 150, 210, 270, 330, 390)
    x0a = p2cX(r.st, axst)
    y0a = p2cY(r.st, axst)
    
    r.end <- c(max1, max2, max3, max4, max5, max6)
    axend <- c(90, 150, 210, 270, 330, 390)
    x1a = p2cX(r.end, axend)
    y1a = p2cY(r.end, axend)
    
    # Set up grid graphics viewport
    
    md <- max(abs(c(x0a, y0a, x1a, y1a)))*1.3 # max dimension
    if (np) grid.newpage()
    grid.rect(gp = gpar(col = NA, fill = bkgnd))
    vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                   xscale = c(-md, md), yscale = c(-md, md), name = "3DHivePlot")
    pushViewport(vp)
    
    # Now draw edges (must do in sets as curvature is not vectorized)
    
    # Axis 1 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 150)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 2 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 150)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 3 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 4 -> 5
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 5 -> 6
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 6)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 390)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 6 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 390)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Axis 1 -> 6
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 6)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 390)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 6 -> 5
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 390)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 5 -> 4
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 4 -> 3
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 3 -> 2
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 150)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 2 -> 1
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 150)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])				}
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = -0.5)
    }
    
    # Axis 1 -> 1, 2 -> 2 etc (can be done as a group since curvature can be fixed)
    
    r.st <- r.end <- th.st <- th.end <- ecol <- ewt <- c()
    
    for (n in 1:nrow(edges)) {
      
      if ((nodes$axis[id1[n]] == 1) & (nodes$axis[id2[n]] == 1)) {
        th.st <- c(th.st, 90)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 90)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 2) & (nodes$axis[id2[n]] == 2)) {
        th.st <- c(th.st, 150)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 150)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 3) & (nodes$axis[id2[n]] == 3)) {
        th.st <- c(th.st, 210)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 210)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 4) & (nodes$axis[id2[n]] == 4)) {
        th.st <- c(th.st, 270)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 270)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 5) & (nodes$axis[id2[n]] == 5)) {
        th.st <- c(th.st, 330)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 330)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
      if ((nodes$axis[id1[n]] == 6) & (nodes$axis[id2[n]] == 6)) {
        th.st <- c(th.st, 390)
        r.st <- c(r.st, nodes$radius[id1[n]])
        th.end <- c(th.end, 390)
        r.end <- c(r.end, nodes$radius[id2[n]])
        ecol <- c(ecol, edges$color[n])
        ewt <- c(ewt, edges$weight[n])
      }
      
    }
    
    x0 = p2cX(r.st, th.st)
    y0 = p2cY(r.st, th.st)
    x1 = p2cX(r.end, th.end)
    y1 = p2cY(r.end, th.end)
    
    if (!length(x0) == 0) {
      grid.curve(x0, y0, x1, y1,
                 default.units = "native", ncp = 5, square = FALSE,
                 gp = gpar(col = ecol, lwd = ewt), curvature = 0.5)
    }
    
    # Draw axes
    
    # grid.segments(x0a, y0a, x1a, y1a,
    # gp = gpar(col = "black", lwd = 7),
    # default.units = "native") # more like linnet
    
    grid.segments(x0a, y0a, x1a, y1a,
                  gp = gpar(col = HPD$axis.cols, lwd = 3),
                  default.units = "native")
    
    # Now add nodes
    
    if (dr.nodes) {
      r <- c(n1$radius, n2$radius, n3$radius, n4$radius, n5$radius, n6$radius) 
      theta <- c(rep(90, length(n1$radius)),
                 rep(150, length(n2$radius)),
                 rep(210, length(n3$radius)),
                 rep(270, length(n4$radius)),
                 rep(330, length(n5$radius)),
                 rep(390, length(n6$radius)))
      x = p2cX(r, theta)
      y = p2cY(r, theta)
      grid.points(x, y, pch = 20, gp = gpar(cex = c(n1$size, n2$size, n3$size, n4$size, n5$size, n6$size),
                                            col = c(n1$color, n2$color, n3$color, n4$color, n5$color, n6$color)))
    }
    
    # Now label axes
    
    if (!is.null(axLabs)) {
      if (!length(axLabs) == nx) stop("Incorrect number of axis labels")
      if (is.null(axLab.gpar)) axLab.gpar <- gpar(fontsize = 12, col = "white")
      r <- c(max1, max2, max3, max4, max5, max6)
      if (is.null(axLab.pos)) axLab.pos <- r*0.1
      r <- r + axLab.pos
      th <- c(90, 150, 210, 270, 330, 390)
      x <- p2cX(r, th)
      y <- p2cY(r, th)
      grid.text(axLabs, x, y, gp = axLab.gpar,  default.units = "native", ...)
    }
    
    # Add a legend arrow & any annotations
    
    if (!is.null(arrow)) addArrow(arrow, nx)
    if (!is.null(anNodes)) annotateNodes(anNodes, nodes, nx, anCoord)
    if (!is.null(grInfo)) addGraphic(grInfo, nodes, nx, anCoord)
    
  } # end of 6D
  
  
} # closing brace, this is the end


