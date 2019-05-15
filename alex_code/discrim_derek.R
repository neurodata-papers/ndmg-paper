require(mgc)
require(igraph)
require(stringr)
require(reshape2)
require(ggplot2)
require(BBmisc)
require(readr)

setwd('/tmp')

discr.distance <- function(X, method='euclidian') {
  D <- as.matrix(dist(X, diag=TRUE, method=method))

  return(D)
}

discr.rdf <- function(D, ids) {
  N <- dim(D)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }

  uniqids <- unique(as.character(ids))
  countvec <- vector(mode="numeric",length=length(uniqids))

  for (i in 1:length(uniqids)) {
    countvec[i] <- sum(uniqids[i] == ids) # total number of scans for the particular id
  }

  scans <- max(countvec) # assume that we will worst case have the most
  rdf <- array(NaN, N*(scans-1)) # initialize empty ra

  count <- 1
  for (i in 1:N) {
    ind <- which(ids[i] == ids) # all the indices that are the same subject, but different scan
    for (j in ind) {
      if (!isTRUE(all.equal(j, i))) { # if j != i, then we want j to have a close distance to i, and estimate where it ranks
        di <- D[i,] # get the entire ra for the particular scan
        di[ind] <- Inf # don't want to consider the particular scan itself
        d <- D[i,j] # the distance between the particular scan of a subject and another scan of the subject
        rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind)) # 1 for less than, .5 if equal, then average
        count <-  count + 1
      }
    }
  }
  return(rdf[1:count-1]) # return only the occupied portion
}

discr.mnr <- function(rdf, remove_outliers=FALSE, thresh=0, output=FALSE) {
  if (remove_outliers) {
    discr <- mean(rdf[which(rdf[!is.nan(rdf)] > thresh)]) # mean of the rdf
    ol <- length(which(rdf<thresh))
    if (output) {
      print(paste('Graphs with reliability <',thresh,'(outliers):', ol))
    }
  } else {
    ol <- 0
    discr <- mean(rdf[!is.nan(rdf)])
  }
  nopair <- length(rdf[is.nan(rdf)])
  if (output) {
    print(paste('Graphs with unique ids:',nopair))
    print(paste('Graphs available for reliability analysis:', length(rdf)-ol-nopair))
    print(paste('discr:', discr))
  }
  return(discr)
}

discr.stat <- function(X, ids, remove_outliers=FALSE, thresh=0, verbose=FALSE) {
  X <- as.matrix(X)
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(X) != ncol(X) | sum(diag(X)^2) > 0){
    X <- discr.distance(X)
  }
  return(discr.mnr(discr.rdf(X, ids), remove_outliers=remove_outliers, thresh=thresh, output=(verbose)))
}  
  
fmriu.list2array <- function(list_in, flatten=FALSE) {
  nroi <- max(sapply(list_in, function(graph) dim(graph)[1]))
  nsub <- length(list_in)
  array_out <- array(NaN, dim=c(nsub, nroi, nroi))
  subnames <- names(list_in)
  incl_ar <- logical(nsub)
  for (i in 1:nsub) {
    if (isTRUE(all.equal(dim(list_in[[i]]), c(nroi, nroi)))) {
      array_out[i,,] <-list_in[[i]]
      incl_ar[i] <- TRUE
    }
  }
  array_out <- array_out[incl_ar,,]
  subnames <- subnames[incl_ar]
  if (flatten) {
    dimar <- dim(array_out)
    if (!is.na(dimar[3])){
      dim(array_out) <- c(dimar[1], dimar[2]*dimar[3])
    }
  }
  return(list(array=array_out, incl_ar=incl_ar, names=subnames))
}

fmriu.io.open_graphs <- function(ipath, dataset_id="", atlas_id="",
                                 fmt='elist', verbose=FALSE, rtype='list', flatten=FALSE,
                                 rem.diag=TRUE) {
  if (! (fmt %in% c('adj', 'elist'))) {
    stop('You have passed an invalid format type. Options are: [\'adj\', \'elist\', and \'graphml\'].')
  }

  if (fmt == 'elist') {
    fmt = 'ncol'; ext = "ssv"
  } else if (fmt == "graphml") {
    fmt = "graphml"; ext = "graphml"
  } else if (fmt == "adj") {
    fmt = "adj"; ext="adj"
  }

  fnames <- paste(ipath, list.files(ipath, pattern = atlas_id), sep='/')
  if (! (rtype %in% c('list', 'array'))) {
    stop('You have passed an invalid return type. Options are: [\'list\', \'array\'].')
  }

  print(sprintf("opening graphs for %s dataset and %s parcellation atlas...", dataset_id, atlas_id))
  subjects <- vector("character", length(fnames))
  sessions <- vector("character", length(fnames))
  tasks <- vector("character", length(fnames))
  gr <- list()

  vertices <- c()

  # so that we don't get any annoying errors if particular vertices are empty
  if (fmt != "adj") {
    for (i in 1:length(fnames)) {
      tgr <- igraph::read_graph(fnames[i], format=fmt) # read the graph from the filename
      vertices <- union(vertices, V(tgr)$name)
      #print(fnames[i])
      #print(length(vertices))
    }
  }

  vertices <- vertices[ordered(as.numeric(vertices))]
  counter <- 1
  for (i in 1:length(fnames)) {
    basename <- basename(fnames[i])     # the base name of the file
    if (verbose) {
      print(paste('Loading', basename, '...'))
    }

    igr <- tryCatch({
      igraph::read_graph(fnames[i], format=fmt, predef=vertices)
    }, error = function(e) {
#      igr <- igraph::read_graph(fnames[i], format=fmt)
#      dd <- length(vertices) - length(V(igr))
#      new_verts <- c(as.character(sort(as.numeric(V(igr)$name))))
#      missing <- setdiff(vertices,new_verts)
#      missing <- c(as.numeric(missing))
#      missing <- c(as.character(sort(as.numeric(missing[!is.na(missing)]))))
#      igr <- add_vertices(igr, dd, name=missing)
#      olde = E(igr)
#      igr[V(igr), V(igr)] <- TRUE
#      E(igr)$weight <- 0
#      E(igr)[olde]$weight <- olde$weight
#      igr <- simplify(igr)
#      igr <- permute(igr, match(V(igr)$name, 1:length(vertices)))
#      return(igr)
      return(NA)
    })
    
    if (is.igraph(igr)) {
      tgr <- get.adjacency(igr, type="both", attr="weight", sparse=FALSE) # convert graph to adjacency matrix
      tgr[is.nan(tgr)] <- 0  # missing entries substituted with 0s
      if (rem.diag) {
        diag(tgr) <- 0
      }
      #tgr = log10(tgr)
      tgr <- normalize(tgr, method = "standardize", range = c(0, 1))
      #print(tgr)
      colnames(tgr) <- V(igr)
      rownames(tgr) <- V(igr)
      gr[[basename]] <- t(tgr)
      subjects[counter] <- str_extract(basename, 'sub(.?)+?(?=_)')
      sessions[counter] <- str_extract(basename, 'ses(.?)+?(?=_)')
      tasks[counter] <- str_extract(basename, 'task(.?)+?(?=_)')
      counter <- counter + 1
    }
  }

  dataset <- rep(dataset_id, counter - 1)
  atlas <- rep(atlas_id, counter - 1)
  subjects <- subjects[1:counter - 1]
  sessions <- sessions[1:counter - 1]
  tasks <- tasks[1:counter - 1]

  if (rtype == 'array') {
    aro <- fmriu.list2array(gr, flatten=flatten)
    gr <- aro$array
    dataset <- dataset[aro$incl_ar]
    atlas <- atlas[aro$incl_ar]
    subjects <- subjects[aro$incl_ar]
    sessions <- sessions[aro$incl_ar]
    tasks <- tasks[aro$incl_ar]
  }
  return(list(graphs=gr, dataset=dataset, atlas=atlas, subjects=subjects,
              sessions=sessions, tasks=tasks))
}

studies <- c('HNU1', 'BNU1', 'SWU', 'NKI', 'KKI')

atlases <- c("desikan")
#atlases <- c("CPAC200", "DKT", "HarvardOxfordcort-maxprob-thr25", "HarvardOxfordsub-maxprob-thr25", "Schaefer2018-200-node", "brodmann", "desikan", "glasser", "yeo-17", "yeo-7")

df <- setNames(data.frame(matrix(ncol = length(studies), nrow = length(atlases))), studies)
row.names(df) <- atlases

for (study in studies) {
  ipath <- paste("/Users/derekpisner/Downloads/graphs_", study, sep="")
  print(study)
  for (atlas in atlases) {
    inputs <- fmriu.io.open_graphs(ipath, fmt="elist", rtype="array", dataset_id=study, atlas_id=atlas, flatten=TRUE, rem.diag=TRUE)  
    if (length(inputs$sessions) != length(inputs$subjects)){
      next
    }
    res <- discr.stat(inputs$graphs, inputs$subjects)
    print(res)
    df[atlas, study] <- res
    remove(inputs)
    closeAllConnections()
  }
}
write_csv(df, "/Users/derekpisner/Downloads/discrim_all.csv")