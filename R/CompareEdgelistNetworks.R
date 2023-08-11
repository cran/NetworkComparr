#' @title Statistical Comparison of Networks based on Edgelists
#'
#' @description
#' A permutation-based hypothesis test for statistical comparison of two networks based on the invariance measures of the R package 'NetworkComparisonTest' by van Borkulo et al. (2022): global strength invariance, network structure invariance, edge invariance, and various centrality measures. Edgelists of two comparable group networks from dependent or independent samples are used as input. These edgelists can be generated from individual graphical network representations, such as concept maps. The networks can be directed or undirected.
#'
#' Keywords: concept-maps, edgelists, comparing networks, permutation-test
#' @author Lara Trani <lara.trani@rptu.de>
#' Maike Sauer <maike.sauer@rptu.de>
#'
#' @references
#'Salmaso, L., & Pesarin, F. (2010). Permutation tests for complex data: theory, applications and software. John Wiley & Sons. https://doi.org/10.1002/9780470689516
#'
#'van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022). Comparing network structures on three aspects: A permutation test. Psychological methods. https://doi.org/10.1037/met0000476.
#'
#' @keywords concept-maps edgelists comparing-networks permutation-test
#'
#' @name NetworkComparr
NULL


#' @title Statistical Comparison of Networks based on Edgelists
#' @name CompareEdgelistNetworks
#' @description This permutation-based hypothesis test assesses the differences between two networks based on the invariance measures of the R package 'NetworkComparisonTest' by van Borkulo et al. (2022): global strength invariance, network structure invariance, edge invariance, and various centrality measures. The global strength invariance represents the differences in the total number of all edges. The network structure invariance involves only the edge that differs the most in its count. The edge invariance covers the differences of the specified individual edges. The centrality measures covers the differences regarding the importance of the individual nodes. Edgelists of two comparable group networks from dependent or independent samples are used as input. These edgelists can be generated from individual graphical network representations, such as concept maps. The networks can be directed or undirected.
#' @usage CEN(EL1, EL2, noc, it = 10000,
#'              paired = FALSE, directed = FALSE, abs = TRUE, plot = TRUE,
#'              test.edges = TRUE, edges = "all", progressbar = TRUE,
#'              p.adjust.methods = c("none", "holm", "hochberg","hommel",
#'                                    "bonferroni", "BH", "BY", "fdr"),
#'              test.centrality = TRUE, centrality = c("all"), cen.nodes = "all",
#'              test.bridge.centrality = FALSE, bridge.centrality = c("all"),
#'              brg.nodes= "all", communities = NULL, useCommunities = "all",
#'              verbose = TRUE)
#' @param EL1 Edgelist 1. One of two edgelists containing the subject number and the coding of the network. The object type is "data.frame", the columns must be labeled "person", "from" and "to". "Person" contains the subject number. For dependent samples, the same person must have the same label in both edgelists. For independent samples, the same subject number cannot occur twice in the two edgelists. "From" indicates where an edge begins, "to" indicates where it ends. For undirected networks, an edge can be specified in one direction only. The direction is irrelevant (e.g., from 1 to 2, or from 2 to 1). The nodes must be numbered and must not have labels of the class character. The example serves as an orientation for the structure of the data frame.
#' @param EL2 Edgelist 2. The other of two edgelists containing the subject number and the coding of the network. The object type is "data.frame", the columns must be labeled "person", "from" and "to". "Person" contains the subject number. For dependent samples, the same person must have the same label in both edgelists. For independent samples, the same subject number cannot occur twice in the two edgelists. "From" indicates where an edge begins, "to" indicates where it ends. For undirected networks, an edge can be specified in one direction only. The direction is irrelevant (e.g., from 1 to 2, or from 2 to 1). The nodes must be numbered and must not have labels of the class character. The example serves as an orientation for the structure of the data frame.
#' @param noc Number of nodes. Total number of nodes that occur in the networks.
#' @param it  Iterations. The number of iterations (permutations). Defaults to 10.000. The ideal number of permutations can be calculated in advance using the R package 'coin'.
#' @param paired Logical. Can be TRUE or FALSE to indicate whether the samples are dependent or not. If paired is TRUE, relabeling is performed within each pair of observations. If paired is FALSE, relabeling is not restricted to pairs of observations. Note that, currently, dependent data is assumed to entail one group measured twice. Defaults to FALSE.
#' @param directed Logical. Can be TRUE or FALSE to indicate whether the networks are directed or not. Undirected networks have edges without directions. For directed networks, the edges are often represented with arrows in one direction, however, the arrow can also point in both directions. Defaults to FALSE.
#' @param abs Logical. Should global strength consider the absolute value of edge weights, or the raw value (i.e., global expected influence)? Important if the edgelists contain negative values. Defaults to TRUE.
#' @param plot Logical. Should the networks be displayed graphically? Defaults to TRUE.
#' @param test.edges Logical. Can be TRUE or FALSE to indicate whether or not differences in individual edges should be tested. Defaults to TRUE.
#' @param edges Character or list. When 'all', differences between all individual edges are tested. When provided a list with one or more pairs of nodes, the provided edges are tested. The list must contain the two nodes that delimit the edge, e.g., edges = list(c(1,2), c(3,2)).
#' @param progressbar Logical. Should the progressbar be plotted in order to see the progress of the permutation procedure? Defaults to TRUE.
#' @param p.adjust.methods Character. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none". To control (or not) for testing of multiple edges. Defaults to "none".
#' @param test.centrality Logical. Should centrality metrics be compared across networks? Defaults to FALSE.
#' @param centrality Type of centrality metrics to test. Can be either "all" or individually selected. For undirected networks it can be any of c("Closeness", "Betweenness", "Strength", "ExpectedInfluence"). For directed networks it can be any of c("Closeness", "Betweenness", "InStrength", "OutStrength", "InExpectedInfluence", "OutExpectedInfluence"). For more information see '\link[qgraph]{centrality}'.
#' @param cen.nodes Specific nodes for centrality tests. Indicated by the node numbers, e.g., c(1,3), to test node 1 and 3. Only used if test.centrality = TRUE.
#' @param test.bridge.centrality Logical. Should bridge centrality metrics be compared across networks? Defaults to FALSE.
#' @param bridge.centrality Type of bridge centrality metrics to test. Currently, only "all" can be selected.
#' @param brg.nodes Specific nodes for bridge centrality tests. Indicated by the node numbers, e.g., c(1,3), to test node 1 and 3. Only used if test.bridge.centrality = TRUE.
#' @param communities A character vector of community assignments for each node (e.g., c("Comm1", "Comm1", "Comm2", "Comm2)). The ordering of this vector should correspond to the order of the nodes. Can also be in list format (e.g., list("Comm1"=c(1:3), "Comm2"=c(4:6))). For more information see '\link[networktools]{bridge}'.
#' @param useCommunities A character vector specifying which communities should be included. Default set to "all". For more information see '\link[networktools]{bridge}'.
#' @param verbose Logical: Should warnings and notes be printed?
#'
#' @return CompareEdgelistNetworks returns an object of class "CEN" containing the following items:
#' @return \strong{glstrinv.real} The difference in global strength between the networks of the observed data sets.
#' @return \strong{glstrinv.sep} The global strength values of the individual networks.
#' @return \strong{glstrinv.pval} The p value resulting from the permutation test concerning difference in global strength.
#' @return \strong{glstrinv.perm}	The difference in global strength between the networks of the permutated data sets.
#' @return \strong{nwinv.real} The value of the maximum difference in edge weights of the observed networks.
#' @return \strong{nwinv.pval} The p value resulting from the permutation test concerning the maximum difference in edge weights.
#' @return \strong{nwinv.perm} The values of the maximum difference in edge weights of the permuted networks.
#' @return \strong{edges.tested} The edges, that are called to be tested. Only if test.edges = TRUE.
#' @return \strong{einv.real} The values of the differences in edge count of the observed networks. This can be either of all edges or of the selected edge(s). Only if test.edges = TRUE.
#' @return \strong{einv.pvals} p-values (corrected for multiple testing or not according to 'p.adjust.methods') per edge from the permutation test concerning differences in edges count. Only returned if test.edges = TRUE.
#' @return \strong{einv.perm}	The values of the differences in edge count of the permuted networks. Only if test.edges = TRUE.
#' @return \strong{diffcen.real} The values of the difference in centralities of the observed networks. Only if test.centrality = TRUE.
#' @return \strong{diffcen.pval} p-values (corrected for multiple testing or not according to 'p.adjust.methods') per node from the permutation test concerning differences in centralities. Only if test.centrality = TRUE.
#' @return \strong{diffcen.perm} The values of the difference in centralities of the permuted networks. Only if test.centrality = TRUE.
#' @return \strong{cen.nw1}	Values of the selected centrality measures of network 1.
#' @return \strong{cen.nw2}	Values of the selected centrality measures of network 2.
#' @return \strong{cen.plot1}	Visualization of the selected centrality measures of network 1.
#' @return \strong{cen.plot2}	Visualization of the selected centrality measures of network 2.
#' @return \strong{diffbcen.real}	The values of the difference in bridge centralities of the observed networks. Only if test.bridge.centrality = TRUE.
#' @return \strong{diffbcen.pval}	p-values (corrected for multiple testing or not according to 'p.adjust.methods') per node from the permutation test concerning differences in bridge centralities. Only if test.bridge.centrality = TRUE.
#' @return \strong{diffbcen.perm}	The values of the difference in bridge centralities of the permuted networks. Only if test.bridge.centrality = TRUE.
#' @return \strong{bridgecen.nw1}	Values of the selected bridge centrality measures of network 1.
#' @return \strong{bridgecen.nw2} Values of the selected bridge centrality measures of network 2.
#' @return \strong{bridge.plot1} Visualization of the selected bridge centrality measures of network 1.
#' @return \strong{bridge.plot2} Visualization of the selected bridge centrality measures of network 2.
#' @return \strong{info} Returns the selected arguments of the CEN usage.
#'
#' @details This function is predominantly based on the R package '\pkg{NetworkComparisonTest}' by van Borkulo et al. (2022). By changing the input in the form of existing edgelists, the usage has been changed.
#' @importFrom stats p.adjust runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' ###Simulate datasets with a dependent sample and directed networks
#' #For an example with an undependent sample and undirected networks see README
#' data1 <- dplyr::tibble(person = c("1","1","2","2","2","3","3","3","4","4","4","4","4","5","5"),
#'                       from =    c("1","3","3","1","4","3","1","3","1","4","3","3","3","2","1"),
#'                       to =      c("3","4","4","3","2","1","3","4","3","3","4","2","4","4","3"))
#' data2 <- dplyr::tibble(person = c("1","2","2","3","3","3","3","3","3","3","4","4","4","5","5","5"),
#'                       from =    c("1","1","1","1","2","1","2","4","4","3","1","3","2","1","1","3"),
#'                       to =      c("2","2","4","2","4","4","3","3","1","4","2","2","1","4","2","4"))
#'
#' ### Compare networks of data sets using CEN
#' Res <- CEN(data1, data2, noc=4, it=50, paired=TRUE, directed=TRUE, abs=TRUE,
#'             test.edges=TRUE, edge=list(c(1,3),c(4,2),c(3,2)), p.adjust.methods= "none",
#'             test.centrality=TRUE, centrality=c("Closeness", "Betweenness"), cen.nodes="all",
#'             test.bridge.centrality=FALSE, bridge.centrality="all", brg.nodes=c(1,3),
#'             communities=c("1","1","2","2"), useCommunities="all")
#'
#' ###See results
#' summary(Res)
#' print(Res)
#' Res$glstrinv.sep
#' Res$glstrinv.pval
#' Res$nwinv.real
#' Res$nwinv.pval
#' Res$einv.real
#' Res$einv.pvals
#' Res$diffcen.real
#' Res$diffcen.pval
#'
#' ###Plot results
#' plot(Res, what="network")
#' plot(Res, what="strength")
#' plot(Res, what="edge")
#' plot(Res, what="centrality")
#'
#'
#' @export
CEN <- function(EL1, EL2, noc, it = 10000,
                paired = FALSE, directed=FALSE, abs = TRUE, plot = TRUE,
                test.edges = TRUE, edges = "all", progressbar =TRUE,
                p.adjust.methods = c("none", "holm", "hochberg", "hommel",
                                     "bonferroni", "BH", "BY", "fdr"),
                test.centrality = TRUE, centrality=c("all"), cen.nodes="all",
                test.bridge.centrality = FALSE, bridge.centrality=c("all"), brg.nodes="all",
                communities=NULL, useCommunities="all", verbose = TRUE){

  #Struktur####
  # store function call, including default arguments not explicitly set - credit to Neal Fultz
  match.call.defaults <- function(...) {
    # Extract explicit call arguments
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    # Extract default call arguments
    formals <- evalq(formals(), parent.frame(1))

    # When extracting calls, the symbols T or F not properly mapped to TRUE and FALSE
    clean.call.arg <- function(arg){
      if(is.null(arg)){
        return(list(arg))
      }
      if(arg == "T"){
        return(list(TRUE))
      }
      if(arg == "F"){
        return(list(FALSE))
      }
      return(list(arg))
    }
    for(i in 1:length(names(call))){
      call[i] <- clean.call.arg(call[[i]])
    }

    # if default argument not explicitly written, add the default value to the saved call
    for(i in setdiff(names(formals), names(call))){
      call[i] <- list(formals[[i]])
    }

    match.call(sys.function(sys.parent()), call)
  }
  cl <- match.call.defaults()

  p.adjust.methods <- match.arg(p.adjust.methods)

  # Fix for networktools example:
  if (missing(edges)) edges <- "all" #Wenn test.edges = TRUE und keine bestimmten edges ausgewählt werden, werden alle edges überprüft
  if (progressbar==TRUE) pb <- txtProgressBar(max=it, style = 3) #Definiert die Progressbar, die während der Durchführung der Funktion angezeigt wird, falls sie angezeigt werden soll



  #Objekte erstellen####
  ####Gesamt-Edgelist####
  ELAll <- rbind(EL1,EL2) #Gemeinsame Edgeslist (EL1 + EL2) erstellen für die Permutationen, bevor VPN geöscht werden
  data1 <- EL1 #Einzeldatensätze mit Personen für die Permutationen später
  data2 <- EL2 #Einzeldatensätze mit Personen für die Permutationen später
  VPN <- nrow(table(data1$person)) #Versuchspersonen Data1 für die Permutationen später

  #### Netzwerke ####
  #Netzwerk 1 aus der Edgelist 1 erstellen
  EL1[ , 'person'] <- list(NULL) #Lösche die Zeile Person" aus den Edgelisten
  EL_M1 <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(EL1))) #Edgelist wird zu Matrix
  nw1 <-  matrix(0, noc, noc) #Leere Matrix mit allen Knoten erstellen, um Matrix EL_M1 hier einzufügen, damit am Ende alle Knoten vorhanden sind und nicht-verbundene Knoten nicht weggelassen werden
  rownames(nw1) <- c(1:noc) #Format der Zeilen-Knotennamen ändern, damit beide Matrixen verknüpft werden können
  colnames(nw1) <- c(1:noc) #Format der Spalten-Knotennamen ändern, damit beide Matrixen verknüpft werden können
  cols <- colnames(nw1)[colnames(nw1) %in% colnames(EL_M1)] #Abspeichern der Spalten-Namen
  rows <- rownames(nw1)[rownames(nw1) %in% rownames(EL_M1)] #Abspeichern der Zeilen-Namen
  nw1[rows, cols] <- EL_M1[rows, cols] #Überführung der Matrix EL_M1 in die leere Matrix mit allen Knoten

  if(directed==FALSE){ #dient dem Zweck, dass bei ungerichteten Netzwerken die Edgelists die Verbindungen nur in eine Richtung beinhalten müssen
    if(isSymmetric(nw1) == FALSE){
    m <- matrix(nrow = noc, ncol = noc) #erstelle eine leere Matrix mit der anzahl der Knoten
    a <- gdata::lowerTriangle(nw1, diag=TRUE, byrow=FALSE) + gdata::upperTriangle(nw1, diag=TRUE, byrow=TRUE) #addiere alle zahlen im unteren dreieck mit dem oberen dreieck
    m[lower.tri(m, diag=TRUE)] <- a #füge die addierten zahlen in die leere Matrix
    diag(m) <- diag(m/2) #Teile die Mittelspalte durch zwei, da diese mit sich selbst addiert wurde
    gdata::upperTriangle(m) = gdata::lowerTriangle(m, byrow=TRUE) #spiegel die matrix, sodass auf beiden seiten die selben zahlen stehen, die nun eben addiert wurden
    nw1 <- m #überschreibe die ursprüngliche matrix nw1 mit der neu erstellten, symmetrischen matrix
    rownames(nw1) <- c(1:noc) #Benenne die Zeilen nach den Knoten
    colnames(nw1) <- c(1:noc) #Benenne die Spalten nach den Knoten
    }
 }

  #Netzwerk 2 aus der Edgelist 2 erstellen
  EL2[ , 'person'] <- list(NULL) #Selbes Prozedere für die zweite Edgelist
  EL_M2 <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(EL2)))
  nw2 <-  matrix(0, noc, noc)
  rownames(nw2) <- c(1:noc)
  colnames(nw2) <- c(1:noc)
  cols <- colnames(nw2)[colnames(nw2) %in% colnames(EL_M2)]
  rows <- rownames(nw2)[rownames(nw2) %in% rownames(EL_M2)]
  nw2[rows, cols] <- EL_M2[rows, cols]

  if(directed==FALSE){ #dient dem Zweck, dass bei ungerichteten Netzwerken die Edgelists die Verbindungen nur in eine Richtung beinhalten müssen
    if(isSymmetric(nw2) == FALSE){
    m <- matrix(nrow = noc, ncol = noc) #erstelle eine leere Matrix mit der anzahl der Knoten
    a <- gdata::lowerTriangle(nw2, diag=TRUE, byrow=FALSE) + gdata::upperTriangle(nw2, diag=TRUE, byrow=TRUE) #addiere alle zahlen im unteren dreieck mit dem oberen dreieck
    m[lower.tri(m, diag=TRUE)] <- a #füge die addierten zahlen in die leere Matrix
    diag(m) <- diag(m/2) #Teile die Mittelspalte durch zwei, da diese mit sich selbst addiert wurde
    gdata::upperTriangle(m) = gdata::lowerTriangle(m, byrow=TRUE) #spiegel die matrix, sodass auf beiden seiten die selben zahlen stehen, die nun eben addiert wurden
    nw2 <- m #überschreibe die ursprüngliche matrix nw2 mit der neu erstellten, symmetrischen matrix
    rownames(nw2) <- c(1:noc) #Benenne die Zeilen nach den Knoten
    colnames(nw2) <- c(1:noc) #Benenne die Spalten nach den Knoten
    }
  }

  #### Weitere Objekte ####
  #Maximale Knotenanzahl
  if(directed == TRUE){
    nedges <- (noc*noc)  #Maximale Kantenanzahl bei noc Knoten, if directed = TRUE
  } else {nedges <- (noc*(noc-1)/2)} #Maximale Kantenanzahl bei noc Knoten, if directed = FALSE

  #Spezifische Knoten für Zentralitätstestung
  valid.cen.nodes <- 1:noc
  cen.nodes.all <- if(cen.nodes[1]=="all") { #Wenn alle Knoten überprüft werden sollen, nimm alle Werte die im Objekt cen.nodes gespeichert sind, ...
  valid.cen.nodes } else { cen.nodes  }

  c.nodes <- ifelse(cen.nodes[1]=="all",noc,length(cen.nodes)) #Hier kannn angegeben werden, welche Knoten auf Zentralität überprüft werden, wird kein bestimmter Knoten ausgewählt, dann alle
  #cen.nodes <- if(is.numeric(cen.nodes)){colnames(nw1)[cen.nodes]} else{cen.nodes}#Wenn Knotenanzahl einen nummerischen Namen haben (bei uns ein muss), nutze diesen Namen, ansonsten werden die nicht-numerischen Knotennamen übernommen

  #Spezifische Knoten für Bridge-Zentralitätstestung
  valid.brg.nodes <- 1:noc
  brg.nodes.all <- if(brg.nodes[1]=="all") { #Wenn alle Knoten überprüft werden sollen, nimm alle Werte die im Objekt brg.nodes gespeichert sind, ...
  valid.brg.nodes } else { brg.nodes  }

  b.nodes <- ifelse(brg.nodes[1]=="all",noc,length(brg.nodes)) #Hier kannn angegeben werden, welche Knoten auf Zentralität überprüft werden, wird kein bestimmter Knoten ausgewählt, dann alle
  #brg.nodes <- if(is.numeric(brg.nodes)){colnames(nw1)[brg.nodes]} else{brg.nodes}#Wenn Knotenanzahl einen nummerischen Namen haben (bei uns ein muss), nutze diesen Namen, ansonsten werden die nicht-numerischen Knotennamen übernommen

  #Knoten, die bei test.edges überprüft werden sollen
  if(is.list(edges)){ #wenn die Knoten die wir auf Zentralitäten überprüfen wollen eine Liste sind, ...
    edges.tested <- edges # ... verwende dies ...
    if(is.character(edges[[1]])){ #... wenn sie ein Character-Objekt sind, verwende das
      whichfun <- function(x){which(colnames(nw1)%in%x)} #??
      edges <- lapply(edges,whichfun)} #?? gehört zu Zeile 71
  }

  #Leere Objekte, in denen später die Zwischenergebnnisse oder Ergebnisse gespeichert werden
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c() #4 Objekte werden erstellt
  diffedges.perm <- matrix(0,it,nedges) #nedges = was ist mit Verbindungen zu sich selbst????
  einv.perm.all <- array(NA,dim=c(noc, noc, it)) #Individual edge invariance
  corrpvals.all <- matrix(NA,noc,noc) #Objekt in dem die p-Werte für ? gespeichert werden
  edges.pvalmattemp <- matrix(0,noc,noc) #Objekt in dem die p-Werte für ? gespeichert werden

  #Leere Objekte, in denen später die Zentraliätsmaße und die Bridge-Zentralitätsmaße gespeichert werden
  #Festlegen der Zentralitätsmaße, wenn das Netzwerk gerichtet bzw. ungerichtet ist
  if(directed == TRUE){
    validCentrality <- c("Closeness", "Betweenness", "InStrength",
                         "OutStrength", "InExpectedInfluence", "OutExpectedInfluence")
  } else {validCentrality <- c("Closeness", "Betweenness", "Strength", "ExpectedInfluence")}

  centrality <- if(centrality[1]=="all") { #Wenn alle Zentralitäten überprüft werden sollen, nimm alle Werte die im Objekt validCentrality gespeichert sind, ...
    validCentrality } else { centrality  } #... ansonsten nimm eben die Zentralitätsmaße, die händisch eingegeben wurden

  diffcen.perm <- matrix(NA, it, c.nodes*length(centrality)) #Leeres Objekt für die Ergebnisse der Zentralitätsmaße

  #Festlegen der Bridge-Zentralitätsmaße, wenn das Netzwerk gerichtet bzw. ungerichtet ist
  if(directed == TRUE){
    bridgecen <- c("BridgeStrengthIn", "BridgeStrengthOut", "BridgeStrength",
                   "BridgeBetweenness", "BridgeCloseness", "BridgeExpectedInfluence")
  } else {bridgecen <- c("BridgeStrength", "BridgeBetweenness", "BridgeCloseness", "BridgeExpectedInfluence")}

  bridge.centrality <- if(bridge.centrality[1]=="all") { #Wenn alle Zentralitäten überprüft werden sollen, nimm alle Werte die im Objekt validCentrality gespeichert sind, ...
    bridgecen } else { bridge.centrality } #... ansonsten nimm eben die Zentralitätsmaße, die händisch eingegeben wurden

  diffbcen.perm <- matrix(NA, it, b.nodes*length(bridgecen)) #Leeres Objekt für die Ergebnisse der Bridge-Zentralitätsmaße







  #Nächstes Kapitel#
  #Testen der Invarianzmaße + Zentralitäsmaße der realen Daten berechnen ####

  ##### Global strength invariance #####

  ifelse(directed == TRUE,          #Wenn das Netzwerk gerichtet ist

         yes = ifelse(abs == TRUE,

                      #Wenn wir die absoluten Werte wollen
                      yes =  c(glstrinv.real <- abs(sum(abs(nw1))-sum(abs(nw2))),     # Global strength invariance: Differenz der Kantenanzahl der beiden Netzwerke als Summe
                               glstrinv.sep <- c(sum(abs(nw1)), sum(abs(nw2)))),   # Global strength of individual networks: Anzahl der Kanten der einzelnen Netzwerke

                      #Wenn wir nicht die absoluten Werte wollen
                      no =   c(glstrinv.real <- abs(sum(nw1)-sum(nw2)),     # Global strength invariance: Differenz der Kantenanzahl der beiden Netzwerke als Summe
                               glstrinv.sep <- c(sum(nw1), sum(nw2)))),   # Global strength of individual networks: Anzahl der Kanten der einzelnen Netzwerke

         no = ifelse(abs == TRUE,   #Wenn das Netzwerk nicht gerichtet ist

                     #Wenn wir die absoluten Werte wollen
                     yes =  c(glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)]))),    # Global strength invariance: Differenz der Kantenanzahl der beiden Netzwerke als Summe
                              glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))),   #Global strength of individual networks: Anzahl der Kanten der einzelnen Netzwerke

                     #Wenn wir nicht die absoluten Werte wollen
                     no =   c(glstrinv.real <- abs(sum(nw1[upper.tri(nw1)])-sum(nw2[upper.tri(nw2)])),    # Global strength invariance: Differenz der Kantenanzahl der beiden Netzwerke als Summe
                              glstrinv.sep <- c(sum(nw1[upper.tri(nw1)]), sum(nw2[upper.tri(nw2)])))))  # Global strength of individual networks: Anzahl der Kanten der einzelnen Netzwerke

  ##### Individual edge invariance #####
  if (directed == TRUE){
    diffedges.real <- abs(nw1-nw2) #Das eine Netzwerk wird vom anderen Netzwerk abgezogen, Zahlen werden als Betrag genommen
    diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE) #Matrix, bei der jede Spalte für eine Kante steht, in der die Unterschiede zwischen NW1 und NW2 stehen, die erste Zeile wird "it" mal wiederholt
    diffedges.realoutput <- abs(nw1-nw2) #Das eine Netzwerk wird vom anderen Netzwerk abgezogen, Zahlen werden als Betrag genommen
  } else {
    diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] #Die obere Hälfte des eine Netzwerk wird von der oberen Hälfte des anderen Netzwerk abgezogen (für ungerichtete Netzwerke), Zahlen werden als Betrag genommen
    diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE) #Matrix, bei der jede Spalte für eine Kante steht, in der die Unterschiede zwischen NW1 und NW2 stehen, die erste Zeile wird "it" mal wiederholt, für ungerichtete Netzwerke, daher nur die halbe Matrix
    diffedges.realoutput <- abs(nw1-nw2) #Das eine Netzwerk wird vom anderen Netzwerk abgezogen, Zahlen werden als Betrag genommen
  }

  ##### Network structure invariance #####
  nwinv.real <- max(diffedges.real) #Maximale Unterschiedlichkeit um zu bestimmen, ob sich die Netzwerke in mindestens einer Kante unterscheiden


  ##### Centrality invariance #####
  if(test.centrality==TRUE){ #wenn Zentralitäten überprüft werden sollen, ...
    if (!all(centrality %in% validCentrality)) { #... und dabei nicht alle Zentralitäten verwendet werden, dann packe in das Objekt "validCentrality" nur die Zentraliäten, die ausgewählt wurden
      stop(paste0("'centrality' must be one of: ", paste0("'", validCentrality, "'", collapse = ", "))) #Sollte dabei eine unbekannte Zentralität verlangt werden, brech ab und gibt eine Warnmeldung raus
    }
    cen1 <- qgraph::centrality_auto(nw1)$node.centrality #Zentralitäten für NW1 werden berechnet
    cen2 <- qgraph::centrality_auto(nw2)$node.centrality #Zentralitäten für NW2 werden berechnet

    if(directed == TRUE){
      names(cen1) <- names(cen2) <- c("Betweenness","Closeness","InStrength","OutStrength", "OutExpectedInfluence", "InExpectedInfluence") #Zentralitäten werden umgenannt, entweder so, beo gerichteten Netzwerken...
    } else {names(cen1) <- names(cen2) <- c("Betweenness","Closeness","Strength", "ExpectedInfluence")} #... und bei ungerichteten Netzwerken eben so
    diffcen.real <- data.matrix(cen1) - data.matrix(cen2) #Die Zentralitätsmaße von NW1 und NW2 werden voneinander abgezogen, dass ein Dataframe mit den Differenzen der Zentralitäten entsteht
  }

  #Centrality Plot
  plot1 <- qgraph::centralityPlot(nw1, include = centrality, print = FALSE, scale = "relative", verbose = FALSE) #Zum Plotten der Zentralitäten ! Gerade werden die relativen Daten geplottet und nicht die Rohdaten - Ändern?
  plot2 <- qgraph::centralityPlot(nw2, include = centrality, print = FALSE, scale = "relative", verbose = FALSE) #Zum Plotten der Zentralitäten ! Gerade werden die relativen Daten geplottet und nicht die Rohdaten - Ändern?

  #Bridge Centrality
  if(test.bridge.centrality==TRUE){#wenn Bridge-Zentralitäten überprüft werden sollen, ...
    if (!all(bridge.centrality %in% bridgecen)) { #... und dabei nicht alle Zentralitäten verwendet werden, dann packe in das Objekt "bridgecen" nur die Bridge-Zentraliäten, die ausgewählt wurden
      stop(paste0("'bridge-centrality' must be one of: ", paste0("'", bridgecen, "'", collapse = ", "))) #Sollte dabei eine unbekannte Zentralität verlangt werden, brech ab und gibt eine Warnmeldung raus
    }
    b1 <- networktools::bridge(nw1, communities=communities, useCommunities=useCommunities, directed = directed) #Bridge-Zentralitäten von NW1 werden berechnet
    b2 <- networktools::bridge(nw2, communities=communities, useCommunities=useCommunities, directed = directed) #Bridge-Zentralitäten von NW2 werden berechnet
    names(b1) <- names(b2) <- c(bridgecen, "bridgeExpectedInfluence2step", "communities")  #Datenobjekt mir Bridge-Zentralitäten von NW1 und NW2 werden umbenannt
    suppressWarnings(plot3 <- plot(b1,include = bridge.centrality, zscore=TRUE)) #Zum Plotten der Zentralitäten ! Gerade werden die relativen Daten geplottet und nicht die Rohdaten - Ändern?
    suppressWarnings(plot4 <- plot(b2,include = bridge.centrality, zscore=TRUE)) #Zum Plotten der Zentralitäten ! Gerade werden die relativen Daten geplottet und nicht die Rohdaten - Ändern?
    b1$communities <- b2$communities <- NULL #Spalte mit den Communities wird gelöscht
    dfb1 <- data.frame(matrix(ncol = 0, nrow = b.nodes)) #Leerer Dataframe mit so vielen Zeilen, wie viele Knoten auf bridge-Zentralität überprüft werden sollen
    dfb2 <- data.frame(matrix(ncol = 0, nrow = b.nodes)) #-||-
    b1 <- data.frame(c(dfb1,b1)) #Bridge-Zentralitäten werden mit den leeren Dataframes kombiniert und dadurch von einem Listen-Objekt in ein data.frame umgewandelt
    b2 <- data.frame(c(dfb2,b2)) #-||-

    diffbcen.real <- data.matrix(b1) - data.matrix(b2) #Dataframes werden in Matrixen umgewandelt und die Bridge-Zentralitäten von Netzwerk 2 werden von Netzwerk 1 abgezogen, sodass eine Differenz-Matrix über bleibt

  }



  #Nächstes Kapitel#
  #Start der Permutationen #####
  ##### Resampling #####
  for (i in 1:it) #Start der Permutationen, das prozedere was jetzt kommt wird "it" mal wiederholt
  {diffedges.permtemp <- matrix(0, noc, noc) #Leere Hülle für die Zwischenergebnisse der edge-invariance der Permutationen

  # If not paired data
  if(paired==FALSE)
  {
    {
      sp <- split(ELAll, ELAll$person) #Erstelle ein Listen-Objekt, dass nach Personen sortiert ist (alle Kanten, die von einer Person entdeckt wurden sind quasi in einem "Kapitel")
      x1perm <- sample(sp,VPN,replace=FALSE) #Erstelle eine Liste, in dem zufällig "Anzahl VPN"-Personen mit ihren stecken
      x2perm <- setdiff(sp, x1perm) #Alle übrigen Personen und ihre Kanten, die nicht in der obrigen Liste sind, sollen in diese Liste
    }


    # Netzwerke aus den eben erstellten Datensätzen erstellen:
    DF <- as.data.frame(do.call(rbind, x1perm)) #Wandel die eben erstelle Liste in einen Datensatz um
    DF[ , 'person'] <- list(NULL) #Lösche die Zeile Personen
    rownames(DF) <- 1:nrow(DF) #Nummeriere den Datensatz
    r1perm <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(DF))) #Erstelle eine Matrix aus diesem Dataframe
    r1leer <-  matrix(0, noc, noc) #Erstelle eine leere Matrix mit Länge noc x noc
    rownames(r1leer) <- c(1:noc) #Nummere ihre Zeilen von 1 bis noc
    colnames(r1leer) <- c(1:noc) #Nummere ihre Spalten von 1 bis noc
    cols <- colnames(r1leer)[colnames(r1leer) %in% colnames(r1perm)] #Die Spaltennnamen von der leeren Matrix und der richtigen Matrix müssen übereinstimmen und werden abgespeichert
    rows <- rownames(r1leer)[rownames(r1leer) %in% rownames(r1perm)] #Ebenso mit Zeilennamen
    r1leer[rows, cols] <- r1perm[rows, cols] #Zeilennamen und Soaltennamen müssen wieder übereinstimmen
    r1perm <- r1leer #Die leere Matrix wird mit der richtigen Matrix (von NW1) befüllt

    DF <- as.data.frame(do.call(rbind, x2perm)) #selbes Prozedere für den 2. permutierten Datensatz
    DF[ , 'person'] <- list(NULL)
    rownames(DF) <- 1:nrow(DF)
    r2perm <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(DF)))
    r1leer <-  matrix(0, noc, noc)
    rownames(r1leer) <- c(1:noc)
    colnames(r1leer) <- c(1:noc)
    cols <- colnames(r1leer)[colnames(r1leer) %in% colnames(r2perm)]
    rows <- rownames(r1leer)[rownames(r1leer) %in% rownames(r2perm)]
    r1leer[rows, cols] <- r2perm[rows, cols]
    r2perm <- r1leer
  }

  # If paired data
  if(paired==TRUE)
  {
    {
      sp <- split(data1, data1$person) #Erstelle ein Listen-Objekt des ersten Datensatzes, dass nach Personen sortiert ist
      VPNX =  floor(runif(1, min=1, max=VPN)) #Erstelle eine Zufallswert die zwischen 1 und der maximalen Versuchspersonenanzahl von Data 1 liegt
      x1permx <- sample(sp,VPNX,replace=FALSE) #Erstelle einen neuen Datensatz der aus dem Listenobjekt zufällig so viele Personen rauszieht, wie zuvor zufällig entschieden wurde
      x2permx <- dplyr::setdiff(sp, x1permx) #Packe in ein zweites Objekt alle Personen aus Data1, die nicht in dem eben erstellen Objekt liegen -> Data1 ist jetzt zufällig auf 2 neue Datensätze aufgeteilt
      DF1 <- as.data.frame(do.call(rbind, x1permx)) #Verwandle das eben erstelle Listenobjekt (DF1) in einen Datensatz
      DF2 <- as.data.frame(do.call(rbind, x2permx)) #Verwandle das eben erstelle Listenobjekt (DF2) in einen Datensatz
      x1permxx <- dplyr::setdiff(data2$person, DF1$person) #Nenne alle Personen, die in Data2 vorkommen nicht aber in DF1
      x2permxx <- dplyr::setdiff(data2$person, DF2$person) #Nenne alle Personen, die in Data2 vorkommen nicht aber in DF2
      DF3 <- data2[data2$person %in% c(x1permxx), ] #Erstelle einen Datensatz mit den Daten aller Personen, die zuvor ermittelt wurden
      DF4 <- data2[data2$person %in% c(x2permxx), ] #Erstelle einen Datensatz mit den Daten aller Personen, die zuvor ermittelt wurden
      x1perm <- rbind(DF1, DF3) #Verbinde den eben neu erstellten Datensatz mit den übrigen Personen von Data2 (DF3) mit den Daten der ergänzenden Personen von Data1 (DF1)
      x2perm <- rbind(DF2, DF4) #Verbinde den eben neu erstellten Datensatz mit den übrigen Personen von Data2 (DF4) mit den Daten der ergänzenden Personen von Data1 (DF2)
    }

    # Netzwerke aus Edgelisten erstellen:
    DF <- x1perm #Nenne den Datensatz zur Weiterverarbeitung kurzzeitig um
    DF[ , 'person'] <- list(NULL) #Lösche die Zeile Personen
    rownames(DF) <- 1:nrow(DF) #Nummeriere den Datensatz
    r1perm <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(DF))) #Erstelle eine Matrix aus diesem Dataframe
    r1leer <-  matrix(0, noc, noc) #Erstelle eine leere Matrix mit Länge noc x noc
    rownames(r1leer) <- c(1:noc) #Nummere ihre Zeilen von 1 bis noc
    colnames(r1leer) <- c(1:noc) #Nummere ihre Spalten von 1 bis noc
    cols <- colnames(r1leer)[colnames(r1leer) %in% colnames(r1perm)] #Die Spaltennnamen von der leeren Matrix und der richtigen Matrix müssen übereinstimmen und werden abgespeichert
    rows <- rownames(r1leer)[rownames(r1leer) %in% rownames(r1perm)] #Ebenso mit Zeilennamen
    r1leer[rows, cols] <- r1perm[rows, cols] #Zeilennamen und Soaltennamen müssen wieder übereinstimmen
    r1perm <- r1leer  #Die leere Matrix wird mit der richtigen Matrix (von NW1) befüllt

    DF <- x2perm #selbes Prozedere für den 2. permutierten Datensatz
    DF[ , 'person'] <- list(NULL)
    rownames(DF) <- 1:nrow(DF)
    r2perm <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(DF)))
    r1leer <-  matrix(0, noc, noc)
    rownames(r1leer) <- c(1:noc)
    colnames(r1leer) <- c(1:noc)
    cols <- colnames(r1leer)[colnames(r1leer) %in% colnames(r2perm)]
    rows <- rownames(r1leer)[rownames(r1leer) %in% rownames(r2perm)]
    r1leer[rows, cols] <- r2perm[rows, cols]
    r2perm <- r1leer
  }


  ## Invariance measures for permuted data
  ##### Invariance measures #####
  #Die Invarianz-Maße der permutierten Daten werden berechnet, nach dem selben Prinzip wie oben, nur Unterschiede werden beschriftet, das i steht immer für die Permutation, die gerade durchläuft

  ## Global strength invariance:
  #Global strength invariance: Differenz der Kantenanzahl der beiden Netzwerke als Summe
  ifelse(directed == TRUE,

         yes = ifelse(abs == TRUE,
                      yes =  glstrinv.perm[i] <- abs(sum(abs(r1perm))-sum(abs(r2perm))),
                      no =   glstrinv.perm[i] <- abs(sum(r1perm)-sum(r2perm))),

         no = ifelse(abs == TRUE,
                     yes =  glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)]))),
                     no =   glstrinv.perm[i] <- abs(sum(r1perm[upper.tri(r1perm)])-sum(r2perm[upper.tri(r2perm)]))))


  ## Individual edge invariance
  if (directed == TRUE){
    diffedges.perm[i,] <- c(abs(r1perm-r2perm)) #Warum das c? Unterschied zu oben, den wir nicht mehr nachvollziehen können, es geht wohl beides, man muss dann nur in der nnächsten Reihe byrow = FALSE setzen (oben TRUE), damit sich die Matrix in der richtigen Reihenfolge auffädelt
    diffedges.permtemp <- matrix(diffedges.perm[i,], byrow=FALSE, ncol = noc) #Matrix der Permutation i wird erstellt
    einv.perm.all[,,i] <- diffedges.permtemp # Alle Permutations-Matrixen werden hier abgelegt

  } else {
    diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
    diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,] #Erstellen der Matrix in zwei Schritten
    diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp) #Schritt 2: die Matrix von Schritt eins eins wird gespiegelt, da bei gerichteten Netzwerken ja nur mit "halben Matrixen" gearbeitet wird
    einv.perm.all[,,i] <- diffedges.permtemp #Matrix der Permutation i wird hier abgelegt, zusammen mit allen Matrixen der Permutationen
  }

  ## Network structure invariance
  nwinv.perm[i] <- max(diffedges.perm[i,])


  ## Zentralitäten
  if(test.centrality==TRUE){
    cen1permtemp <- qgraph::centrality_auto(r1perm)$node.centrality
    cen2permtemp <- qgraph::centrality_auto(r2perm)$node.centrality
    if(directed == TRUE){ names(cen1permtemp) <- names(cen2permtemp) <- c("Betweenness","Closeness","InStrength","OutStrength", "OutExpectedInfluence", "InExpectedInfluence")
    } else {names(cen1permtemp) <- names(cen2permtemp) <- c("Betweenness","Closeness","Strength", "ExpectedInfluence")}

    diffcen.permtemp <- as.matrix(cen1permtemp) - as.matrix(cen2permtemp)
    if(cen.nodes[1]=="all"){diffcen.perm[i,] <- reshape2::melt(diffcen.permtemp[,centrality])$value #Wenn die Zentralitäten aller Knoten berachtet wird, verschmelze alle Ergebnisse in einem Objekt und verwende als Spaltennamen die Namen der Zentralitäten...
    } else {diffcen.perm[i,] <- reshape2::melt(diffcen.permtemp[cen.nodes,centrality])$value} #... oder bei nicht-allen Knoten erstelle das obrige Objekt eben nur mit den ausgewählten Knoten #funktioniert noch nicht
  }

  ## Bridge-Zentralitäten #funktioniert noch nicht
  if(test.bridge.centrality==TRUE){
    b1permtemp <- networktools::bridge(r1perm, communities=communities, useCommunities=useCommunities, directed = directed)
    b2permtemp <- networktools::bridge(r2perm, communities=communities, useCommunities=useCommunities, directed = directed)
    names(b1permtemp) <- names(b2permtemp) <- c(bridgecen, "bridgeExpectedInfluence2step", "communities")
    b1permtemp$communities <- b2permtemp$communities <- NULL
    dfpermb1 <- data.frame(matrix(ncol = 0, nrow = b.nodes))
    dfpermb2 <- data.frame(matrix(ncol = 0, nrow = b.nodes))
    b1permtemp <- data.frame(c(dfpermb1,b1permtemp))
    b2permtemp <- data.frame(c(dfpermb2,b2permtemp))

    diffbcen.permtemp <- data.matrix(b1permtemp) - data.matrix(b2permtemp)
    if(brg.nodes[1]=="all"){diffbcen.perm[i,] <- reshape2::melt(diffbcen.permtemp[,bridgecen])$value #Wenn die Bridge-Zentralitäten aller Knoten berachtet wird, verschmelze alle Ergebnisse in einem Objekt und verwende als Spaltennamen die Namen der Bridge-Zentralitäten...
    } else {diffbcen.perm[i,] <- reshape2::melt(diffbcen.permtemp[brg.nodes,bridgecen])$value}  #... oder bei nicht-allen Knoten erstelle das obrige Objekt eben nur mit den ausgewählten Knoten
  }


  if (progressbar==TRUE) setTxtProgressBar(pb, i) #die Progressbar soll so für jede Permutation (i) weiterlaufen, bis alle durch sind
  }

  #End permutations#







  #Calculate p-values #####
  #### Mit Individual Edge Invariance ####
  if(test.edges==TRUE) #Wenn die Individual Edge Invariance berechnet werden soll, dann berechne die folgenden Werte dafür
  {
    # vector with uncorrected p values
    edges.pvaltemp <- (colSums(diffedges.perm >= diffedges.realmat) + 1) / (it + 1) #p-Werte-Berechnung mit den Kantenunterschieden der realen Daten und den der Permuationen

    ## If all edges should be tested: falls nicht spezifische Kanten ausgewählt wurden, die überprüft werden sollen
    if(is.character(edges)) #Wenn "edges" ein Charakterobkejt ist, "alle", dann berechne das Folgende
    {
      # corrected p-values (or not if p.adjust.methods='none')
      corrpvals.all.temp <- p.adjust(edges.pvaltemp, method=p.adjust.methods) #Korrigiere die berechneten p-Werte mit der ausgewählten Methode #Muss die Adjustierung an die Mittelzeile angepasst werden?

      # matrix with corrected p values
      if (directed == TRUE){
        corrpvals.all <- matrix(corrpvals.all.temp, byrow=FALSE, ncol = noc) #Erstelle eine Matrix (Kn x Kn), in denen die p-Werte erhateln sind, für gerichtete Netzwerke...
      } else {
        corrpvals.all[upper.tri(corrpvals.all,diag=FALSE)] <- corrpvals.all.temp #oder ungerichtete Netzwerke
      }

      rownames(corrpvals.all) <- colnames(corrpvals.all) <- colnames(nw1) #Bennene die Zeilen/Spalten nach den Knotennamen
      einv.pvals <- reshape2::melt(corrpvals.all, na.rm=TRUE, value.name = 'p-value') #Erstelle einen Dataframe mit 3 Spalten, für jeglische Kanten die es gibt und den p-Wert dazu
      einv.perm <- einv.perm.all #Objekt umbenennen als Endergebnis
      einv.real <- diffedges.realoutput #Objekt umbenennen als Endergebnis
      if (directed == TRUE){
        einv.pvals <- cbind(einv.pvals, c(einv.real)) #p-Werte mit Differenzmatrix kombinieren
      } else {einv.pvals <- cbind(einv.pvals, (einv.real[upper.tri(einv.real)]))} #p-Werte mit Differenzmatrix kombinieren
      colnames(einv.pvals) <- c('from', 'to', 'p-value', "differences") #benenne die Spalten wie angegeben

      edges.tested <- "all" #Hier wird hinterlegt, welche Edges getestet wurden, in diesem Fall ja alle

    }

    ## If a selection of edges should be tested: falls  spezifische Kanten ausgewählt wurden, die überprüft werden sollen
    if(is.list(edges)) #Wenn "edges" eine Liste ist, "Liste mit spezifischen Knoten die untersucht werden sollen", dann berechne das Folgende
    {
      einv.perm <- matrix(NA,it,length(edges)) #Erstelle eine leere Matrix mit der Länge der edges, die getestet werden sollen (Spalten) x der Anzahl der Permutationen (it, Zeilen)
      colnames(einv.perm) <- edges #Benenne die Spalten so, wie die edges heißen, die Untersucht werden sollen c(1,3): Kante von Knoten 1 zu Knoten 3
      uncorrpvals <- einv.real <- pairs <- c() #Erstelle leere Objekte
      edges.tested <- as.data.frame(do.call(rbind, edges)) #Verwandel die Liste in einen Dataframe
      colnames(edges.tested) <- c('from', 'to') #Bennene die Zeilen des Dataframes

      # matrix with uncorrected p values
      if(directed == TRUE){
        edges.pvalmattemp <- matrix(c(edges.pvaltemp),byrow=FALSE, nrow=noc) #Erstelle eine leere Matrix und packe da die p-Werte der einzelnen Edges rein für gerichtete Netzwerke
        diag(edges.pvalmattemp) <- 0 #Die Mittelzeile wird zu NULL (macht NA nicht mehr Sinn, oder weglassen??)
      } else {
        edges.pvalmattemp[upper.tri(edges.pvalmattemp,diag=FALSE)] <- edges.pvaltemp #Erstelle eine Matrix und fülle die obere Hälfte mit den p-Werten der ganzen Edges, weil die Matrixhälften bei ungerichteten Netzwerken ja identisch sind
        edges.pvalmattemp <- edges.pvalmattemp + t(edges.pvalmattemp) #Spiegel die Matrix, dass die p-Werte nun auf beiden Hälften stehen
      }


      for(j in 1:length(edges)) #mache das Folgende für die jeweilige edge (j), die überprüft werden soll
      {
        pairs <- rbind(pairs, c(colnames(nw1)[edges[[j]][1]], colnames(nw1)[edges[[j]][2]]))
        uncorrpvals[j] <- edges.pvalmattemp[edges[[j]][1],edges[[j]][2]]
        einv.real[j] <- diffedges.realoutput[edges[[j]][1],edges[[j]][2]]
        for(l in 1:it){
          einv.perm[l,j] <- einv.perm.all[,,l][edges[[j]][1],edges[[j]][2]]
        }
      } #im letzen Abschnitt werden Objekte mit Ergebnissen erstellt, die sich auf die edge invariance beziehen, z.B. der Permutationen oder realen Daten und nur die spezifischen edges erhalten, die sich angeschaut werden sollen

      corrpvals <- p.adjust(uncorrpvals, method=p.adjust.methods) #die unkorrigierten p-Werte der edge incariance werden korrigiert
      corrpvals<-round(corrpvals,5)
      corrpvals_mat <- matrix(NA,length(edges),3) #es entsteht eine Matrix, mit 3 Spalten und so vielen Zeilen, wie sich edges angeschaut werden
      corrpvals_mat[,3] <- corrpvals #In diese Matrix werden die korrigierten p-Werte gesteckt
      corrpvals_mat[,1:2] <- pairs #?
      einv.pvals <- as.data.frame(corrpvals_mat) #die Matrix wird in einen Datensatz umgewandelt
      einv.pvals <- cbind(einv.pvals, einv.real) #p-values mit differenzmatrix verbinden
      colnames(einv.pvals) <- c('from', 'to', 'p-value', "differences") #die Spalten des Datensatzes werden umbenannt
    }


    res <- list(glstrinv.real = glstrinv.real, #Hier werden alle Ergebnisse in einem Obejkt (res) abgespeichert, z.T. werden noch p-Werte berechnet und dann werden diese abgespeichert
                glstrinv.sep = glstrinv.sep,
                glstrinv.pval = (sum(glstrinv.perm >= glstrinv.real) + 1) / (it + 1),
                glstrinv.perm = glstrinv.perm,
                nwinv.real = nwinv.real,
                nwinv.pval = (sum(nwinv.perm >= nwinv.real) + 1) / (it + 1),
                nwinv.perm = nwinv.perm,
                edges.tested = edges.tested,
                einv.real = einv.real,
                einv.pvals = einv.pvals,
                einv.perm = einv.perm,
                directed = directed
          )

  }

  if (progressbar==TRUE) close(pb) #Die Progressbar schließt sich nun



  #### Ohne Individual Edge Invariance ####
  if(test.edges==FALSE) #wenn die Individual Edge Invariance nicht berechnet werden soll, brauchen wir die Berechnungen darüber nicht
  {
    res <- list( #Hier werden alle Ergebnisse in einem Obejkt (res) abgespeichert, z.T. werden noch p-Werte berechnet und dann werden diese abgespeichert
      glstrinv.real = glstrinv.real,
      glstrinv.sep = glstrinv.sep,
      glstrinv.pval = (sum(glstrinv.perm >= glstrinv.real) + 1) / (it + 1),
      glstrinv.perm = glstrinv.perm,
      nwinv.real = nwinv.real,
      nwinv.pval = (sum(nwinv.perm >= nwinv.real) + 1) / (it + 1),
      nwinv.perm = nwinv.perm,
      directed = directed
    )
  }

  #### Plots ####
  if (plot == TRUE) {
    qgraph::qgraph(nw1, vsize = 4.5,  label.cex = 1.5, title = "network 1", title.cex = 1.3)
    qgraph::qgraph(nw2, vsize = 4.5, label.cex = 1.5, title = "network 2", title.cex = 1.3)
  }

  #### Zentralitäten ####
  if(test.centrality){ #Wenn die Zentralitäten berechnet werden sollen, mache das folgende
    if(cen.nodes[1]=="all"){ #Wenn alle Knoten betrachtet werden
      diffcen.real.vec <- reshape2::melt(diffcen.real[,centrality])$value #Erstelle ein Dataframe mit den Daten der Zentralitätsdifferenzen für alle Knoten
    } else { diffcen.real.vec <- reshape2::melt(diffcen.real[cen.nodes,centrality])$value #oder erstelle ein Dataframe mit den Daten der Zentralitätsdifferenzen für die ausgewählten Knoten #funktioniert noch nicht
    }
    diffcen.realmat <- matrix(diffcen.real.vec, it, c.nodes*length(centrality), byrow = TRUE) #Wandle den eben erstellen Datensatz in eine Matrix um
    diffcen.pvaltemp <- (colSums(abs(diffcen.perm) >= abs(diffcen.realmat)) + 1) / (it + 1) #Nutze diese Matrix und die Matrix der Permutationen um die Signifikanz zu berechnen
    diffcen.HBall <- p.adjust(diffcen.pvaltemp, method = p.adjust.methods) #Korrigiere die p-Werte falls gewünscht
    diffcen.pval <- matrix(diffcen.HBall, c.nodes, length(centrality)) #Verpacke die korrigierten p-Werte als Matrix
    diffcen.real <-  matrix(diffcen.real.vec, nrow=c.nodes,ncol=length(centrality)) #Verpacke die realenDifferenz-Werte der Zentralitäten als Matrix
    colnames(diffcen.pval) <- colnames(diffcen.real) <- centrality #Benenne die Spalten so, wie die Zentralitäten heißen

    res[["diffcen.real"]] <- diffcen.real #Lege das Ergebnis im Objekt "res" ab
    res[["diffcen.perm"]] <- diffcen.perm #Lege das Ergebnis im Objekt "res" ab
    res[["diffcen.pval"]] <- diffcen.pval #Lege das Ergebnis im Objekt "res" ab
                        cen1 <- cen1[c(cen.nodes.all),]
                        cen2 <- cen2[c(cen.nodes.all),]
    res[["cen.nw1"]] <- subset(cen1, select = centrality) #Lege das Ergebnis im Objekt "res" ab
    res[["cen.nw2"]] <- subset(cen2, select = centrality) #Lege das Ergebnis im Objekt "res" ab
    res[["cen.plot1"]] <- plot1 #Lege das Ergebnis im Objekt "res" ab
    res[["cen.plot2"]] <- plot2 #Lege das Ergebnis im Objekt "res" ab

    #Column & row names Zentralitäten
    if(cen.nodes[1]=="all"){ #Hier werden die Spalten & Zeilennamen nochmal umbenannt, je nachdem ob die Zentralitäten aller Knoten überprüft werden sollten...
      rownames(res[["diffcen.real"]]) <- rownames(res[["diffcen.pval"]]) <- colnames(nw1)
      colnames(res[["diffcen.perm"]]) <- apply(expand.grid(colnames(nw1), centrality), 1, paste, collapse=".")
    }
    else {
      rownames(res[["diffcen.real"]]) <- rownames(res[["diffcen.pval"]]) <- cen.nodes #... oder nur bestimmter Knoten #funktioniert noch nicht
    }
  }


  #### Bridge-Zentralitäten ####
  if(test.bridge.centrality){ #Wenn die Bridge-Zentralitäten berechnet werden sollen, mache das folgende, selbes Prinzip wie bei Zentralitäten
    if(brg.nodes[1]=="all"){
      diffbcen.real.vec <- reshape2::melt(diffbcen.real[,bridge.centrality])$value
    } else { diffbcen.real.vec <- reshape2::melt(diffbcen.real[brg.nodes,bridge.centrality])$value
    }
    diffbcen.realmat <- matrix(diffbcen.real.vec, it, b.nodes*length(bridgecen), byrow = TRUE)
    diffbcen.pvaltemp <- (colSums(abs(diffbcen.perm) >= abs(diffbcen.realmat)) + 1) / (it + 1)
    diffbcen.HBall <- p.adjust(diffbcen.pvaltemp, method = p.adjust.methods)
    diffbcen.pval <- matrix(diffbcen.HBall, b.nodes, length(bridgecen))
    diffbcen.real <-  matrix(diffbcen.real.vec, nrow=b.nodes,ncol=length(bridgecen))
    colnames(diffbcen.pval) <- colnames(diffbcen.real) <- bridgecen

    res[["diffbcen.real"]] <- diffbcen.real
    res[["diffbcen.perm"]] <- diffbcen.perm
    res[["diffbcen.pval"]] <- diffbcen.pval
                        b1 <- b1[c(brg.nodes.all),]
                        b2 <- b2[c(brg.nodes.all),]
    res[["bridgecen.nw1"]] <- data.matrix(b1)
    res[["bridgecen.nw2"]] <- data.matrix(b2)
    res[["bridge.plot1"]] <- plot3 #Lege das Ergebnis im Objekt "res" ab
    res[["bridge.plot2"]] <- plot4 #Lege das Ergebnis im Objekt "res" ab

    #Column & row names Bridge-Zentralitäten
    if(brg.nodes[1]=="all"){
      rownames(res[["diffbcen.real"]]) <- rownames(res[["diffbcen.pval"]]) <- colnames(nw1)
      colnames(res[["diffbcen.perm"]]) <- apply(expand.grid(colnames(nw1), bridge.centrality), 1, paste, collapse=".")
    }
    else {
      rownames(res[["diffbcen.real"]]) <- rownames(res[["diffbcen.pval"]]) <- brg.nodes
    }
  }



  res$info$call <- cl #Irgendwas mit der Progressbar vielleicht??

  class(res) <- "CEN" #Gibt dem Objekt "res" die Klasse "CEN"
  return(res) #Gebe nach Ende der Funktion das Objekt "res" aus
}

