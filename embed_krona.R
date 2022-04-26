
#' Creates a Krona chart (https://github.com/marbl/Krona) from a phyloseq object
#' and writes it to HTML, and/or embeds it in RStudio or HTML documents
#' generated from R-Markdown.
#' 
#' Embedding in other document types (PDF, DOCX, etc.) requires taking a snapshot, 
#' which is currently achieved by some Javascript hacks. They are tested with
#' version 2.8.1. and will not work with all older versions.
#' Install 'webshot' to allow embedding snapshots.
#'
#' @param physeq: The phyloseq object
#' @param output: Optional HTML output file. If not set, a temporary file will
#'   be created and embedded without retaining the HTML file.
#' @param display: (boolean) Should the output be displayed? Setting this to `FALSE`
#'   makes only sense if combined with `output` to allow generating a HTML
#'   file from phyloseq without embedding
#' @param tax_cols: (optional) Character vector with names of taxonomy columns.
#'   If not set, assumes all columns to be part of the taxonomic lineage.
#' @param group: (optional) Name of a sample variable for grouping samples
#' @param color_col: (optional) Name of a taxonomy column to color by.
#'   Should be a numeric vector stored as character in the phyloseq taxonomy table.
#'   When aggregating over the taxonomic ranks (`tax_cols`), the average of all
#'   individual taxa is calculated, weighted by their abundance. NAs are currently
#'   not allowed.
#' @param color_name: (optional) Name to display in the Krona chart for the
#'   coloring according to `color_col`.
#' @param install_dir: Installation directory containing the scripts and the
#'   KronaTools directory (https://github.com/marbl/Krona), which has to be manually
#'   downloaded into that directory.
#'   Only required with `color_col`, since we have a custom import script,
#'   otherwise it is enough to have ktImportText installed and accessible in path.
#' @param snapshot: (optional, logical) Directly produce a snapshot instead of
#'   embedding an interactive chart. Will automatically be enabled if not in RStudio.
#'   This will take `snapshot_delay` seconds.
#' @param snapshot_format: (optional) Snapshot format to use (default: pdf). All
#'   formats accepted by whebshot can be used (png, jpeg, pdf). PDF will embed
#'   vector graphics.
#' @param snapshot_dim: Snapshot dimensions (in inches)
#' @param snapshot_res: Snapshot resolution (pixels per inch) used when including
#'   the snapshot with knitr::include_graphics(). Only determines the output size,
#'   if output size is controlled by the chunk, this setting will not have any effect.
#' @param snapshot_fontsize: Default font size for snapshots. By default, this is
#'   larger than normal, in order to allow snapshots at sufficient resolution.
#' @param snapshot_delay: Delay (in seconds) before the snapshot is taken. This is
#'   necessary, because we have to wait for the initial animation to finish.
#'   There seems to exist no simple workaround for this behavior.
#' @param snapshot_chart_size: Factor controlling the chart size for snapshots.
#'   This also relies on a Javascript hack and works the same way as the buttons
#'   increasing / decreasing the chart size.
#' @param iframe_width: (optional, if knitting to HTML) Width of the Iframe (CSS units)
#' @param iframe_height: (optional, if knitting to HTML) Height of the Iframe (CSS units)
#' @param ... More styling of the snapshot. Key-value settings are directly passed
#'   as URL parameters. Currently available: dataset (number), node (number),
#'   collapse (true/false), showMagnitude (false/true), 
#'   color (false/true) [automatically true if using color_col], depth (number)
#'   font (number), key (true/false)
#'
#' @details # Details
#'   Embedding requires the HTML files and snapshots to be stored somewhere. 
#'   If `output` is specified, this file is created and then embedded. Snapshots
#'   will have the same name with the extension specified by `snapshot_format`.
#'   If `output` is not specified, a directory called *_krona_files* will be
#'   created in the current directory, where the HTML files will be deposited,
#'   named with a timestamp. If knitting an R-Markdown document with krona
#'   charts embedded, the knitr working directory is used (knitr::current_input(dir=T)).
#'
plot_krona = function (physeq,
                       output = NULL,
                       display = TRUE,
                       tax_cols = colnames(tax_table(physeq)),
                       group = NULL,
                       color_col = NULL,
                       color_name = color_col,
                       install_dir = NULL,
                       snapshot = NULL,
                       snapshot_format = 'png',
                       snapshot_res = 300,
                       snapshot_dim = c(6, 5),
                       snapshot_delay = 2,
                       snapshot_fontsize = 23,
                       snapshot_chart_size = 0.7,
                       iframe_width='100%',
                       iframe_height='500px',
                       verbose = F,
                       ...)
{
  # check for packages
  required = c('htmltools', 'phyloseq')
  packages = rownames(installed.packages())
  missing = setdiff(required, packages)
  if (length(missing)) {
    stop(sprintf('The following packages are missing: %s.
               embed_krona.R is intended to be used with R-Markdown.',
                 paste(missing, collapse=', ')))
  }
  require(phyloseq)
  
  if (!display && is.null(output)) {
    stop('display = FALSE, but no ouptput file defined.')
  }

  if (!is.null(color_col) && is.null(install_dir)) {
    stop('If using "color_col", it is required to manually place the KronaTools
         directory (downloaded from https://github.com/marbl/Krona/releases/latest)
         next to embed_krona.R')
  }
  
  if (is.null(install_dir)) {
    r = suppressWarnings(try(system2('ktImportText', stdout=F, stderr=F, wait=T), silent=T))
    if (r == 0) {
      bin = 'ktImportText'
    } else { 
      stop('ktImportText not found in path. Is Krona installed?')
    }
  } else if (!is.null(color_col)) {
    bin = file.path(install_dir, 'krona_import.pl')
  } else {
    bin = file.path(install_dir, 'KronaTools', 'scripts', 'ImportText.pl')
    if (!file.exists(bin)) {
      stop('install_dir supplied, but KronaTools not found in the directory')
    }
  }

  # determine output file location
  
  get_file = function(dir) {
    while (T) {
      prefix = paste0('krona_', strftime(Sys.time(), "%H.%M.%OS"))
      fname = paste0(prefix, '.html')
      f = file.path(dir, fname)
      if (!file.exists(f)) {
        return(f)
      }
    }
  }
  
  html_out = isTRUE(suppressWarnings(try(knitr::is_html_output(), silent=T)))
  # FIXME: not sure if this is the best way to detect these settings
  in_rstudio = isTRUE(suppressWarnings(try(is.null(knitr::opts_knit$get('rmarkdown.pandoc.to')), silent=T)))
  in_console = !isTRUE(suppressWarnings(try(rstudioapi::getActiveDocumentContext()$id != "#console", silent=T)))

  if (html_out) {
    outdir = gsub('\\.[^\\.]*$', '_krona_files', knitr::current_input(dir = T), perl = T)
    dir.create(outdir, F, T)
    outfile = get_file(outdir)
  } else {
    if (is.null(output)) {
      outdir = '_krona_files'
      dir.create(outdir, F, T)
      outfile = get_file(outdir)
    } else {
      outdir = NULL # not used
      outfile = output
    }
  }
  taxdir = paste0(gsub('\\.[^\\.]+$', '', outfile, perl=T), '_tax_input')
  dir.create(taxdir, F, T)
  
  # write to text
  
  taxa_with_abund = function(phy) {
    t = as(tax_table(phy), 'matrix')
    t = as.data.frame(t, stringsAsFactors=F)
    n = taxa_sums(phy)
    stopifnot(rownames(t) == names(n))
    cbind(Abundance=n, t,  stringsAsFactors = F)
  }
  
  write_tax = function(phy, val) {
    val = gsub('[ /\\|]', '_', val, perl = T)
    t = taxa_with_abund(phy)
    # combine all taxa with same lineage
    ta = aggregate(t['Abundance'], t[tax_cols], sum)
    cols_before = 'Abundance'
    if (!is.null(color_col)) {
      s = as.numeric(t[, color_col])
      stopifnot(!is.na(s))
      t[, color_col] = s * t['Abundance']
      ta$score = aggregate(t[color_col], t[tax_cols], sum)[[color_col]] / ta$Abundance
      cols_before = c(cols_before, 'score')
    }
    ta = ta[c(cols_before, tax_cols)]
    f = file.path(taxdir, sprintf("%s_taxonomy.txt", val))
    write.table(ta, f, sep='\t', row.names=F, col.names=F, na='', quote=F)
    paste(f, val, sep = ',')
  }
  
  # generate the Krona chart
  
  args = list()
  if (!is.null(color_col)) {
    args = c(args, list('-c', paste0('"', color_name, '"')))
  }
  if (is.null(group)) {
    args = c(args, write_tax(physeq, 'all'))
  } else {
    values = as.factor(get_variable(physeq, group))
    stopifnot(!is.na(values))
    taxfiles = sapply(levels(values), function(val) {
      phy = prune_samples(values == val, physeq)
      write_tax(phy, val)
    })
    args = c(args, taxfiles)
  }
  cmd = paste(c(bin, args, "-o", outfile), collapse = ' ')
  #print(paste(cmd, collapse=' '))
  system(cmd, ignore.stdout = T)
  
  # clean up
  unlink(taxdir, T)

  # Return the output
  
  if (is.null(snapshot)) {
    snapshot = !html_out & !in_rstudio
  }
  
  if (verbose) {
    cat(sprintf('Embedding krona from %s.\nIn RStudio: %s\nIn console: %s\nSnapshot: %s\nHTML output: %s\n',
                outfile, in_rstudio, in_console, snapshot, html_out), file=stderr())
  }

  # take a snapshot if necessary
  
  if (snapshot) {
    packages = rownames(installed.packages())
    if (!('webshot' %in% packages)) {
      stop("Please install webshot (https://github.com/wch/webshot)
           for displaying Krona widgets in PDF documents.")
    }
    html = readChar(outfile, file.info(outfile)$size)
    # modify JavaScript code to automatically take a snapshot after loading
    # tested with version 2.8.1
    html = gsub('win.document.write', 'document.write', html, fixed=T)
    html = gsub('snapshotButton\\.disabled *= *false', 'snapshot()', html, ignore.case=T)
    html = gsub('Download Snapshot', '', html, ignore.case=T)
    html = gsub('var bufferFactor *= *[^;]+', sprintf('var bufferFactor = %.3f', 0.1/snapshot_chart_size), html)
    # make sure that the original HTML is not changed
    prefix = gsub('\\.[^\\.]+$', '', outfile, perl=T)
    snap_html_out = paste0(prefix, '_snapshot.html')
    writeChar(html, snap_html_out)
    snap = paste0(prefix, '.', snapshot_format)
    param = list(...)
    if (!is.null(color_col)) {
      param['color'] = 'true'
    }
    if (is.null(param$font)) {
      param$font = snapshot_fontsize
    }
    args = paste(sapply(names(param), function(p)
      paste(p, param[[p]], sep = '=')), collapse = '&')
    url = sprintf('file:///%s?%s', normalizePath(snap_html_out), args)
    snapshot_dim = round(snapshot_dim * snapshot_res)
    webshot::webshot(url,
                     snap,
                     delay = snapshot_delay,
                     vwidth = snapshot_dim[1],
                     vheight = snapshot_dim[2],
                     debug=verbose)
    unlink(snap_html_out)
    if (display)
      return(knitr::include_graphics(snap, dpi = snapshot_res))
  }
  
  # ... or embed / browse the output
  if (display) {
    if (in_console) {
      browseURL(outfile)
    } else if (html_out & !in_rstudio) {
      # R-Markdown HTML output: use Iframe
      return(htmltools::tags$iframe(
        src = outfile,
        width=iframe_width, 
        height=iframe_height,
        scrolling = 'no',
        seamless = 'seamless',
        frameBorder = '0'
      ))
    } else {
      html = readChar(outfile, file.info(outfile)$size)
      return(htmltools::tags$iframe(
        srcdoc = html,
        width=iframe_width, 
        height=iframe_height,
        scrolling = 'no',
        seamless = 'seamless',
        frameBorder = '0'
      ))
    }    
  }
}

