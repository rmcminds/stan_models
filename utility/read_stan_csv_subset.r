read_stan_csv_subset <- function (csvfiles, col_major = TRUE, params = NULL, keep = TRUE) {
    ##mostly copied from rstan's read_stan_csv, but with additional code to enable the loading of just a subset of parameters
    if (length(csvfiles) < 1)
        stop("csvfiles does not contain any CSV file name")
    g_skip <- 10
    ss_lst <- vector("list", length(csvfiles))
    cs_lst2 <- vector("list", length(csvfiles))
    for (i in seq_along(csvfiles)) {
        header <- rstan:::read_csv_header(csvfiles[i])
        lineno <- attr(header, "lineno")
        vnames <- strsplit(header, ",")[[1]]
        ##
        if(!is.null(params)) {
            oldnames <- vnames
            alwayskeep <- 'lp__|accept_stat__|stepsize__|treedepth__|n_leapfrog__|divergent__|energy__'
            if(keep)
                vnames <- grep(paste(c(params,alwayskeep),collapse='|'), vnames, value=T)
            else
                vnames <- grep(paste(c(params,alwayskeep),collapse='|'), vnames, value=T, invert=T)
            inds <- which(oldnames %in% vnames)
        }
        ##
        iter.count <- attr(header, "iter.count")
        variable.count <- length(vnames)
        df <- structure(replicate(variable.count, list(numeric(iter.count))),
            names = vnames, row.names = c(NA, -iter.count), class = "data.frame")
        comments = character()
        con <- file(csvfiles[[i]], "rb")
        buffer.size <- min(ceiling(1e+06/variable.count), iter.count)
        row.buffer <- matrix(ncol = variable.count, nrow = buffer.size)
        row <- 1
        buffer.pointer <- 1
        while (length(char <- readBin(con, "int", size = 1L)) > 0) {
            seek(con, origin = "current", -1)
            if (char == 35) {
                line <- readLines(con, n = 1)
                comments <- c(comments, line)
                next
            }
            if (char == 108) {
                readLines(con, n = 1)
                next
            }
            if (char == 10) {
                readLines(con, n = 1)
                next
            }
            row.buffer[buffer.pointer, ] <- scan(con, nlines = 1,
                sep = ",", quiet = TRUE)[inds]
            if (buffer.pointer == buffer.size) {
                df[row:(row + buffer.size - 1), ] <- row.buffer
                row <- row + buffer.size
                buffer.pointer <- 0
            }
            buffer.pointer <- buffer.pointer + 1
        }
        if (buffer.pointer > 1) {
            df[row:(row + buffer.pointer - 2), ] <- row.buffer[1:(buffer.pointer -
                1), ]
            df <- df[1:(row + buffer.pointer - 2),]
        } else {
            df <- df[1:(row - 1),]
        }
        close(con)
        cs_lst2[[i]] <- rstan:::parse_stancsv_comments(comments)
        if ("output_samples" %in% names(cs_lst2[[i]]))
            df <- df[-1, ]
        ss_lst[[i]] <- df
    }
    m_name <- sub("(_\\d+)*$", "", rstan:::filename_rm_ext(basename(csvfiles[1])))
    sdate <- do.call(max, lapply(csvfiles, function(csv) file.info(csv)$mtime))
    sdate <- format(sdate, "%a %b %d %X %Y")
    chains <- length(ss_lst)
    fnames <- names(ss_lst[[1]])
    n_save <- nrow(ss_lst[[1]])
    paridx <- rstan:::paridx_fun(fnames)
    lp__idx <- attr(paridx, "meta")["lp__"]
    par_fnames <- c(fnames[paridx], "lp__")
    pars_oi <- rstan:::unique_par(par_fnames)
    dims_oi <- lapply(pars_oi, function(i) {
        pat <- paste("^", i, "(\\.\\d+)*$", sep = "")
        i_fnames <- par_fnames[grepl(pat, par_fnames)]
        rstan:::get_dims_from_fnames(i_fnames, i)
    })
    names(dims_oi) <- pars_oi
    midx <- if (!col_major)
        rstan:::multi_idx_row2colm(dims_oi)
    else 1:length(par_fnames)
    if (chains > 1) {
        if (!all(sapply(ss_lst[-1], function(i) identical(names(i),
            fnames))))
            stop("the CSV files do not have same parameters")
        if (!all(sapply(ss_lst[-1], function(i) identical(length(i[[1]]),
            n_save))))
            stop("the number of iterations are not the same in all CSV files")
    }
    mode <- 0L
    samples <- lapply(ss_lst, function(df) {
        ss <- df[c(paridx, lp__idx)[midx]]
        attr(ss, "sampler_params") <- df[setdiff(attr(paridx,
            "meta"), lp__idx)]
        ss
    })
    par_fnames <- par_fnames[midx]
    for (i in seq_along(samples)) {
        attr(samples[[i]], "adaptation_info") <- cs_lst2[[i]]$adaptation_info
        attr(samples[[i]], "args") <- list(sampler_t = cs_lst2[[i]]$sampler_t,
            chain_id = cs_lst2[[i]]$chain_id)
        if (cs_lst2[[i]]$has_time)
            attr(samples[[i]], "elapsed_time") <- rstan:::get_time_from_csv(cs_lst2[[i]]$time_info)
    }
    save_warmup <- sapply(cs_lst2, function(i) i$save_warmup)
    warmup <- sapply(cs_lst2, function(i) i$warmup)
    thin <- sapply(cs_lst2, function(i) i$thin)
    iter <- sapply(cs_lst2, function(i) i$iter)
    if (!rstan:::all_int_eq(warmup) || !rstan:::all_int_eq(thin) || !rstan:::all_int_eq(iter))
        stop("not all iter/warmups/thin are the same in all CSV files")
    n_kept0 <- 1 + (iter - warmup - 1)%/%thin
    warmup2 <- 0
    if (max(save_warmup) == 0L) {
        n_kept <- n_save
    }
    else if (min(save_warmup) == 1L) {
        warmup2 <- 1 + (warmup[1] - 1)%/%thin[1]
        n_kept <- n_save - warmup2
    }
    if (n_kept0[1] != n_kept) {
        warning("the number of iterations after warmup found (",
            n_kept, ") does not match iter/warmup/thin from CSV comments (",
            paste(n_kept0, collapse = ","), ")")
        if (n_kept < 0) {
            warmup <- warmup + n_kept
            n_kept <- 0
            mode <- 2L
        }
        n_kept0 <- n_save
        iter <- n_save
        for (i in 1:length(cs_lst2)) {
            cs_lst2[[i]]$warmup <- warmup
            cs_lst2[[i]]$iter <- iter
        }
    }
    idx_kept <- if (warmup2 == 0)
        1:n_kept
    else -(1:warmup2)
    for (i in seq_along(samples)) {
        m <- vapply(samples[[i]], function(x) mean(x[idx_kept]),
            numeric(1))
        attr(samples[[i]], "mean_pars") <- m[-length(m)]
        attr(samples[[i]], "mean_lp__") <- m["lp__"]
    }
    perm_lst <- lapply(1:chains, function(id) sample.int(n_kept))
    sim = list(samples = samples, iter = iter[1], thin = thin[1],
        warmup = warmup[1], chains = chains, n_save = rep(n_save,
            chains), warmup2 = rep(warmup2, chains), permutation = perm_lst,
        pars_oi = pars_oi, dims_oi = dims_oi, fnames_oi = rstan:::dotfnames_to_sqrfnames(par_fnames),
        n_flatnames = length(par_fnames))
    null_dso <- new("cxxdso", sig = list(character(0)), dso_saved = FALSE,
        dso_filename = character(0), modulename = character(0),
        system = R.version$system, cxxflags = character(0), .CXXDSOMISC = new.env(parent = emptyenv()))
    null_sm <- new("stanmodel", model_name = m_name, model_code = character(0),
        model_cpp = list(), dso = null_dso)
    nfit <- new("stanfit", model_name = m_name, model_pars = pars_oi,
        par_dims = dims_oi, mode = mode, sim = sim, inits = list(),
        stan_args = cs_lst2, stanmodel = null_sm, date = sdate,
        .MISC = new.env(parent = emptyenv()))
    return(nfit)
}
