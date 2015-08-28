plot.RVineMatrix <- function(x, tree = "ALL", type = 0, edge.labels = NULL, legend.pos = "bottomleft", interactive = FALSE, ...) {
    
    M <- x$Matrix
    d <- nrow(M)
    
    ## sanity checks
    if (!inherits(x, "RVineMatrix")) 
        stop("'x' has to be an RVineMatrix object.")
    if (tree != "ALL" && tree > d - 1) 
        stop("Selected tree does not exist.")
    if (any(tree == "ALL") )
        tree <- 1:(d - 1)
    if (!all(type %in% c(0, 1, 2)))
        stop("type not implemented")
    stopifnot(is.logical(interactive))
    
    ## set names if empty
    if (is.null(x$names)) 
        x$names <- paste("V", 1:d, sep = "")
    
    #### set up plotting options ----------------------------
    # reduce default margins of plot range 
    usr <- par()$mar
    par(mar = c(1.1,0.1,3.1,0.1)) 
    on.exit(par(mar = usr))
    
    # set plot.network options
    TUMlightblue <- rgb(red = 100, green = 160, blue = 200, maxColorValue = 255)
    dflt <- list(interactive = interactive,
                 displaylabels = TRUE,
                 pad = 1.5e-1,
                 edge.lwd = 0.35,
                 edge.col = "gray43",
                 boxed.labels = TRUE,
                 label.pad = 1.5,
                 label.bg = TUMlightblue,
                 label.pos = 7,
                 label.col = "gray97",
                 label.cex = 1.3,
                 vertex.cex = 0,
                 object.scale = 0.05)
    # Same color for edges, edge labels and label borders 
    dflt <- append(dflt, list(label.border   = dflt$edge.col,
                              edge.label.col = dflt$edge.col,
                              edge.label.cex = dflt$label.cex - 0.2)) 
    
    ## overwrite defaults with ... argument
    lst <- list(...)
    temp.args <- modifyList(dflt, lst) 
    
    #### loop through the trees -----------------------------
    for (i in tree) {
        
        main <- list(main = paste("Tree ", i, sep = ""), 
                     col.main = ifelse("col.main" %in% names(temp.args),
                                       temp.args$col.main, 
                                       temp.args$edge.col)) 
        final.args <- append(temp.args, main)
        
        ## create network object
        g <- makeNetwork(x, i, !(type %in% c(0, 2)))
        final.args$x = g$nw
        
        ## set edge labels
        if (!is.null(edge.labels)) 
            final.args$edge.label <- set_edge_labels(tree = i,
                                                     RVM = x,
                                                     edge.labels = edge.labels, 
                                                     type = type)
        
        do.call(plot, final.args)
        
        ## add legend
        if (type == 2) {
            legend(legend.pos, 
                   legend = paste(1:d, x$name, sep = " \U002194 "),
                   bty = "n", 
                   xjust = 0, 
                   text.col = final.args$edge.col,
                   cex = final.args$label.cex)
        }
        
        ## wait for key stroke 
        if (i != max(tree)) {
            par(ask = TRUE)
        } else {
            par(ask = FALSE)
        }
    }
}


## -----------------------------------------------------------------------------
## contour generic for RVineMatrix objects
contour.RVineMatrix <- function(x, tree = "ALL", xylim = NULL, cex.nums = 1, ...) {
    
    ## check input
    d <- nrow(x$Matrix)
    if (all(tree == "ALL"))
        tree <- seq.int(d-1)
    n.tree <- length(tree)
    if (!is.null(list(...)$type)) 
        stop("Only contour plots allowed. Don't use the type argument!")
    
    ## set up for plotting windows
    mfrow.usr <- par()$mfrow
    mar.usr <- par()$mar
    par(mfrow = c(n.tree, d - min(tree)))
    par(mar = rep(0, 4))
    on.exit(par(mfrow = mfrow.usr, mar = mar.usr))
    
    
    ## default style --------------------------------------------------
    # headings: blue color scale from dichromat pacakge 
    cs <- 1 / 255 * t(col2rgb(c("#E6FFFF",  
                                "#CCFBFF",
                                "#B2F2FF", 
                                "#99E6FF", 
                                "#80D4FF", 
                                "#66BFFF", 
                                "#4CA6FF", 
                                "#3388FF",
                                "#1A66FF",
                                "#0040FF"))) 
    # contours: set limits for plots
    if (!is.null(list(...)$margins)) {
        margins <- list(...)$margins
        if (!(margins %in% c("norm", "unif")))
            stop("margins not supported")
    } else {
        margins <- "norm"
    }
    if (is.null(xylim))
        xylim <- switch(margins,
                        "norm" = c(-3, 3),
                        "unif" = c(1e-1, 1 - 1e-1))
    xlim <- ylim <- xylim
    
    # contours: adjust limits for headings
    offs <- 0.25 
    mult <- 1.5
    ylim[2] <- ylim[2] + offs*diff(ylim)
    
    
    ## run through trees -----------------------------------------------
    # initialize check variables 
    cnt <- 0
    k <- d
    e <- numeric(0)
    class(e) <- "try-error"
    
    while ("try-error" %in% class(e)) {
        e <- try({
            maxnums <- get_num(1, tree = max(tree), RVM = x)
            for (i in tree) {
                for (j in 1:(d - min(tree))) {
                    if (d - i >= j) {
                        # set up list of contour arguments
                        args <- list(x = BiCop(family=x$family[d-i+1,j],
                                               par=x$par[d-i+1,j],
                                               par2=x$par2[d-i+1,j]),
                                     drawlabels = FALSE,
                                     xlab = "",
                                     ylab = "",
                                     xlim = xlim, 
                                     ylim = ylim,
                                     xaxt = "n",
                                     yaxt = "n")
                        
                        # call plot.BiCop with ... arguments
                        do.call(plot, modifyList(args, list(...)))
                        
                        # draw area for headings
                        abline(h = ylim[2] - diff(ylim)/mult*offs)
                        ci <- min(nrow(cs) + 1 - i, 10)
                        polygon(x = c(xlim[1] - diff(xlim),
                                      xlim[1] - diff(xlim), 
                                      xlim[2] + diff(xlim),
                                      xlim[2] + diff(xlim)),
                                y = c(2*ylim[2],
                                      ylim[2] - diff(ylim)/mult*offs, 
                                      ylim[2] - diff(ylim)/mult*offs, 
                                      2*ylim[2]),
                                col = rgb(cs[ci, 1], cs[ci, 2], cs[ci, 3], 0.3))
                        
                        # add pair-copula ID
                        cx1 <- 0.95 * diff(xlim) / strwidth(maxnums)
                        cx1 <- cx1
                        ty <- ylim[2] - diff(ylim)/mult*offs
                        cx2 <- 0.95 * (ylim[2] - ty) / strheight(maxnums)
                        cx2 <- cx2
                        cx <- min(cx1, cx2)
                        text(x = sum(xlim)/2,
                             y = ty + 0.225 / cex.nums * (ylim[2] - ty), 
                             cex    = cex.nums * cx,
                             labels = get_num(j, tree = i, RVM = x), 
                             pos    = 3,
                             offset = 0)
                    } else {
                        plot.new()
                    }
                }
            }
        }
        , silent = TRUE)
        
        ## adjust to figure margins if necessary
        if (length(tree) < 1)
            stop("Error in plot.new() : figure margins too large")
        if ("try-error" %in% class(e)) {
            cnt <- cnt + 1
            tree <- tree[-which(tree == max(tree))]
            par(mfrow = c(n.tree - cnt, d - min(tree)))
        }        
    }
    
    ## message for the user if not all trees could be plotted -----------
    if (length(tree) != n.tree) {
        nmbr.msg <- as.character(tree[1])
        if (length(tree) > 2) {
            for (i in tree[-c(1, length(tree))]) {
                nmbr.msg <- paste(nmbr.msg, i, sep=", ")
            }
        }
        if (length(tree) > 1) {
            s.msg <- "s "
            nmbr.msg <- paste(nmbr.msg,
                              "and",
                              tree[length(tree)],
                              "were plotted. ")
        } else {
            s.msg <- " "
            nmbr.msg <- paste(nmbr.msg, "was plotted. ", sep=" ")
        }    
        msg.space <- "There is not enough space."
        msg.tree <- paste("Only Tree",
                          s.msg,
                          nmbr.msg,
                          "Use the 'tree' argument or enlarge figure margins",
                          " to see the others.",
                          sep = "")
        message(paste(msg.space, msg.tree))
    }
}


## creates a network object for a tree in a given RVineMatrix ------------------
makeNetwork <- function(RVM, tree, use.names = FALSE) {
    M <- RVM$Matrix
    d <- ncol(M)
    
    I <- matrix(0, d - tree + 1, d - tree + 1)
    
    ## extract node and edge labels as numbers
    if (tree > 1) {
        node.lab <- sapply(1:(d - tree + 1),
                           get_num,
                           tree = tree - 1,
                           RVM = RVM)
    } else {
        node.lab <- paste(diag(M))
    }
    edge.lab <- sapply(seq.int(d - tree),
                       get_num,
                       tree = tree, 
                       RVM = RVM)
    
    ## convert to numeric matrices V and E
    V <- t(sapply(strsplit(node.lab,  " *[,;] *"), as.numeric))
    V <- matrix(V, ncol = tree)
    E <- t(sapply(strsplit(edge.lab,  " *[,;] *"), as.numeric))
    
    ## build incident matrix by matching V and E
    for (i in 1:nrow(E)) {
        ind.i <- which(apply(V, 1, function(x) all(x %in% E[i, ])))
        I[ind.i[1], ind.i[2]] <- I[ind.i[1], ind.i[2]] <- 1
    }
    
    ## convert to variable names (if asked for)
    if (use.names) {
        if (tree > 1) {
            node.lab <- sapply(1:(d - tree + 1),
                               get_name,
                               tree = tree - 1,
                               RVM = RVM)
        } else {
            node.lab <- RVM$names
        }
    }
    
    ## create network
    colnames(I) <- rownames(I) <- node.lab
    nw <- network(I, directed = FALSE)
    
    ## return network and labels
    list(nw = nw, vlabs = node.lab)
}


## finds appropriate edge labels for the plot ----------------------------------
set_edge_labels <- function(tree, RVM, edge.labels, type) {
    d <- nrow(RVM$Matrix)
    if (edge.labels[1] == "family") {
        elabel <- sapply(1:(d - tree + 1), 
                         get_family, 
                         tree = tree, 
                         RVM = RVM)
        elabel <- BiCopName(as.numeric(elabel))
    } else if (edge.labels[1] == "par") {
        elabel <- sapply(1:(d - tree + 1),
                         get_par,
                         tree = tree,
                         RVM = RVM)
    } else if (edge.labels[1] == "tau") {
        elabel <- sapply(1:(d - tree + 1),
                         get_tau,
                         tree = tree,
                         RVM = RVM)
    } else if (edge.labels[1] == "family-par") {
        elabel1 <- sapply(1:(d - tree + 1),
                          get_family, 
                          tree = tree, 
                          RVM = RVM)
        elabel1 <- BiCopName(as.numeric(elabel1))
        elabel2 <- sapply(1:(d - tree + 1),
                          get_par,
                          tree = tree, 
                          RVM = RVM)
        elabel <- paste0(elabel1, "(", elabel2, ")")
        elabel <- sapply(elabel, 
                         function(x){
                             tmp <- gsub("((", "(", x, fixed = TRUE)
                             gsub("))", ")", tmp, fixed = TRUE)
                         })
    } else if (edge.labels[1] == "family-tau") {
        elabel1 <- sapply(1:(d - tree + 1),
                          get_family, 
                          tree = tree, 
                          RVM = RVM)
        elabel1 <- BiCopName(as.numeric(elabel1))
        elabel2 <- sapply(1:(d - tree + 1),
                          get_tau,
                          tree = tree, 
                          RVM = RVM)
        elabel <- paste0(elabel1, "(", elabel2, ")")
    } else if (length(edge.labels) > 1) {
        # user may provide own labels
        if (length(edge.labels) == d - tree) {
            elabel <- as.character(edge.labels)
        } else {
            stop("length of edge.labels does not equal the number of edges in the tree")
        }
    } else if (edge.labels[1] == "pair"){
        if (type %in% c(0, 2)) {
            elabel <- sapply(1:(d - tree + 1), 
                             get_num, 
                             tree = tree, 
                             RVM = RVM)
        } else {
            elabel <- sapply(1:(d - tree + 1), 
                             get_name, 
                             tree = tree, 
                             RVM = RVM)
        }
    } else {
        stop("edge.labels not implemented")
    }
    
    elabel
}


## get info for a pair-copula from RVineMatrix object --------------------------
get_num <-  function(j, tree, RVM) {
    M <- RVM$Matrix
    d <- nrow(M)
    # get numbers from structure matrix
    nums <- as.character(M[c(j, (d - tree + 1):d), j])
    # conditioned set
    bef <- paste(nums[2],
                 nums[1],
                 sep = ",",
                 collapse = "")
    # conditioning set
    aft <- if (length(nums) > 2) {
        gsub(" ",
             ",",
             do.call(paste, as.list(as.character(nums[3:length(nums)])))) 
    }  else ""
    # paste together
    sep <- if (length(nums) > 2) " ; " else ""
    paste(bef, aft, sep = sep, collapse = "")
}

get_name <-  function(j, tree, RVM) {
    M <- RVM$Matrix
    d <- nrow(M)
    # variable names
    nams <- RVM$names[M[c(j, (d - tree + 1):d), j]]
    # conditioned set
    bef <- paste(nams[2],
                 nams[1],
                 sep = ",",
                 collapse = "")
    # conditioning set
    aft <- if (length(nams) > 2) {
        gsub(" ",  ",", do.call(paste, as.list(nams[3:length(nams)]))) 
    }  else ""
    # paste together
    sep <- if (length(nams) > 2) " ; " else ""
    paste(bef, aft, sep = sep, collapse = "")
}

get_family <- function(j, tree, RVM) {
    d <- nrow(RVM$family)
    M <- RVM$Matrix
    paste(RVM$family[M[d - tree + 1, j]])
}

get_par <- function(j, tree, RVM) {
    d <- nrow(RVM$family)
    M <- RVM$Matrix
    # get parameters
    par  <- round(RVM$par[M[d - tree + 1, j]], digits = 2)
    par2 <- round(RVM$par2[M[d - tree + 1, j]], digits = 2)
    # add brackets if par2 != 0
    apply(cbind(par, par2), 1, join_par)
}

join_par <- function(x) {
    if (x[2] != 0) 
        return(paste0("(", x[1], ",", x[2], ")"))
    x[1]
}

get_tau <- function(j, tree, RVM) {
    d <- nrow(RVM$family)
    M <- RVM$Matrix
    # get family and parameters
    family <- RVM$family[M[d - tree + 1, j]]
    par  <- RVM$par[M[d - tree + 1, j]]
    par2 <- RVM$par2[M[d - tree + 1, j]]
    # convert to Kendall's tau
    tau <- BiCopPar2Tau(family, par, par2, check.pars = FALSE)
    round(tau, digits = 2)
}

