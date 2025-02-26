#' Creates a Krona chart (https://github.com/marbl/Krona) from a phyloseq object
#' and writes it to HTML, and/or embeds it in RStudio or HTML documents
#' generated from R-Markdown.
#'
#' @param classification: Character matrix or data frame with a hierarchical
#'   classification (e.g. taxonomic lineages).
#' @param magnitude: (optional) Numeric matrix, data frame or a single vector
#'   with magnitudes (abundances) for each of the rows (lineages) in the
#'   `classification` matrix. If not specified, abundances are either taken
#'   from the phyloseq object, or equal magnitudes are assumed for all lineages.
#'   The row order of `classification` and `magnitude` must match, as well as
#'   the row names (if present).
#' @param physeq: Phyloseq object with `tax_table` and `otu_table` components
#'   (and `sample_data` if `group_vars` are specified).
#'   Separate `classification` and `magnitude` components are not needed in this case.
#' @param group: (optional) Character vector of length `ncol(abund_tab)` specifying 
#'   the grouping of columns in `magnitude` into different data sets. 
#'   If not specified, one single global data set is assumed.
#'   If 'all', all columns of the `magnitude` matrix are assumed to represent
#'   separate data sets.
#'   The magnitudes are accumulated by default (but configurable using `summary_fn`). 
#' @param tax_vars: (optional) Character vector with taxonomy variables to be
#'   chosen from the phyloseq object or community matrix.
#'   Default is to use all columns except for `color_var` and all variables
#'   listed in `attributes`.
#' @param group_vars: (optional) Character vector with sample variables that 
#'   should be used used for creating a `group` vector from the phyloseq object.
#'   The columns will be joined together using `group_sep`.
#' @param group_sep: Separator to use for joining `group_vars`
#' @param shorten_group: Function used to shorten overly long strings resulted
#'   from joining different grouping variables (`group_vars`)
#' @param output: Optional output file. If not set, a temporary file will
#'   be created and embedded without retaining the HTML file.
#' @param display: (boolean) Should the output be displayed? Setting this to `FALSE`
#'   makes only sense if combined with `output` to allow generating a HTML
#'   file from phyloseq without embedding
#' @param tax_cols: (optional) Character vector with names of community columns.
#'   If not set, assumes all columns to be part of the taxonomic lineage.
#' @param color_values: (optional) Vector of values used as basis for coloring 
#'   the chart, should be of length nrow(taxonomy) and match the order of the taxa.
#'   At higher taxonomic levels (or in case of duplicate taxa names), the values
#'   will be summarized using `color_summary_fn` (`weighted.mean` by default).
#'   NAs are currently not allowed.
#' @param color_col: (optional) Name of a taxonomy column from the phyloseq
#'   object used to generate `color_values`
#' @param color_label: (optional) Label to display in the Krona chart for the
#'   coloring according to `color_values` (with phyloseq input, the corresponding
#'   column name is used).
#' @param color_value_range: (optional) Limits for the color scale (see
#'   `color_values`/`color_col`). By default, the full range is used, but the
#'   range can be restricted / adjusted here.
#' @param hue_range: Hue at the start and end points of the fill gradient
#'   (mapping to the value minimum/maximum set in `color_value_range`)
#' @param unknown_label: Label to use for NAs in the taxonomy table (if any present).
#'   Default is to leave the NAs in place, resulting in fewer child nodes. Krona
#'   charts handle this situation well (see e.g. https://krona.sourceforge.net/examples/xml.krona.html).
#' @param root_label: Label to use for the root of the taxonomic tree. Only 
#'   displayed if there is more than one taxonomic group in the leftmost
#'   column of the taxonomy table.
#' @param total_label: Label for the total abundance displayed on the top right
#'   of the chart.
#' @param summary_fn: Function used to summarize the abundances from several 
#'   samples (columns of `abund_tab`). `sum` is the default, but `mean` is also
#'   a good choice.
#' @param color_summary_fn: Function used to summarize values from the 
#'  `color_values` vector across taxonomic groups. Takes a vector of color
#'  values as first argument and a vector of taxa abundances as second argument.
#'  The default is `weighted.mean`, which calculates the average of the
#'  color values, weighted by the taxa abundances.
#' @param kronatools_dir: (optional) Path to the 'KronaTools' source directory,
#'   which was extracted from the Krona source code archive 
#'   (https://github.com/marbl/Krona/releases/latest).
#' @param resources_url: (optional) URL where Javascript and image resources should
#'   be obtained from. This means an internet connection is needed when displaying
#'   the charts, but also will reduce the interactive HTML chart files (especially
#'   useful if several charts are embedded).
#'   Example: resources_url = 'http://krona.sourceforge.net'
#'   Specifically, the given address should serve the contents of the 'src'
#'   and 'img' directories from KronaTools.
#' @param snapshot: (optional, logical) Produce a snapshot image instead of
#'   embedding an interactive chart. Will automatically be enabled if not in RStudio.
#'   This will take `snapshot_delay` seconds.
#' @param snapshot_format: (optional) Snapshot format to use (default: png)
#'   formats accepted by webshot can be used (png, pdf, jpeg). PDF generates
#'   vector graphics and can be embedded in PDF or DOCX/ODT documents, while
#'   PNG is best for the remaining formats. JPEG results in a black background
#'   in many cases.
#' @param snapshot_res: Snapshot resolution (pixels per inch)
#' @param snapshot_dim: Snapshot dimensions (in inches)
#' @param snapshot_fontsize: Font size for snapshots
#' @param snapshot_delay: Delay (in seconds) before the snapshot is taken. This is
#'   necessary, because we have to wait for the initial animation to finish.
#'   There seems to exist no simple way to get rid of the animation.
#' @param snapshot_chart_size: Factor controlling the chart size for snapshots.
#'   Works the same way as the buttons increasing / decreasing the chart size.
#' @param iframe_width: (optional, if knitting to HTML) Width of the Iframe (CSS units)
#' @param iframe_height: (optional, if knitting to HTML) Height of the Iframe (CSS units)
#' @param krona_opts: More styling of the snapshot. A list of key-value 
#'   parameters, which are passed to the krona chart as URL parameters. 
#'   Currently available: 
#'     dataset (number), node (number),
#'     collapse (true/false), 
#'     showMagnitude ('false'/'true'),
#'     color ('false'/'true') [automatically true if using color_col], 
#'     depth (number), 
#'     font (number),
#'     key (true/false)
#' @param minify: Compress the included Javascript code using 
#'   UglifyJS (https://lisperator.net/uglifyjs). This makes files somewhat smaller.
#'   This is ignored if 'resources_url' was specified (which has an even larger
#'   effect on file size).
#' @param verbose: Output some more information
#' @param debug: Print debug information
#'
#' @details
#'   Embedding in other document types (PDF, DOCX, etc.) requires taking a snapshot,
#'   which is currently achieved by modifying the Javascript code of the charts.
#'   Ths has been tested with KronaTools version 2.8.1 and will not work with 
#'   all of the older versions.
#'   
#'   Embedding interactive charts requires 'htmltools' (should be present if
#'   rmarkdown is installed) and the 'js' package (if minify=TRUE).
#'   Taking snapshots requires the 'webshot' package (also make sure to run
#'   `webshot::install_phantomjs()` after the installation).
#'   
#'   Temporary output files will be stored in the Knitr plot directory 
#'   (named after the R-Markdown file with a '_files' suffix). This directory
#'   is usually removed after the rendering of the final document is finished.
#'   If using 'plot_krona' from the R console, files will be stored in a directory
#'   named 'krona_files', unless some output file is specified.
plot_krona <- function(community, ...) {
  UseMethod('plot_krona')
}


plot_krona.phyloseq = function (physeq,
                                magnitude = NULL,
                                tax_vars = NULL,
                                group_vars = NULL,
                                color_var = NULL,
                                color_values = NULL,
                                color_label = color_var,
                                group_sep = ' ',
                                shorten_group = function(x) substr(x, 1, 60),
                                ...)
{
  require(phyloseq)
  stopifnot('phyloseq' %in% class(physeq))
  stopifnot(is.character(group_sep) && length(group_sep) == 1)
  stopifnot(is.function(shorten_group) && length(shorten_group) == 1)
  
  taxonomy = as(tax_table(physeq), 'matrix')

  if (is.null(magnitude)) {
    magnitude = as(otu_table(physeq), 'matrix')
    if (!taxa_are_rows(physeq))
      magnitude = t(magnitude)
  } else {
    magnitude = as(magnitude, 'matrix')
    stopifnot(nrow(magnitude) == nrow(taxonomy))
  }
  
  # remove zero-abundance taxa
  taxa_subset = rowSums(magnitude) > 0
  magnitude = magnitude[taxa_subset, , drop=F]
  taxonomy = taxonomy[taxa_subset, , drop=F]
  if (!is.null(color_values)) {
    color_values = color_values[taxa_subset]
  }
  
  group = if (!is.null(group_vars)) {
    stopifnot(is.character(group_vars))
    sdata = as.data.frame(sample_data(physeq))
    group = do.call(paste, c(sdata[, group_vars], list(sep=group_sep)))
    unname(shorten_group(group))
  } else {
    NULL
  }
  
  if (!is.null(color_var)) {
    if (!is.null(color_values)) {
      warning("plot_krona: both 'color_var' and 'color_values' supplied. Using 'color_var'")
    }
    stopifnot(is.character(color_var) && length(color_var) == 1)
    stopifnot(color_var %in% colnames(taxonomy))
    color_values = as.numeric(taxonomy[, color_var])
  }
  
  if (is.null(tax_vars))
    tax_vars = setdiff(colnames(taxonomy), color_var)
  stopifnot(is.character(tax_vars))
  
  taxonomy = taxonomy[, tax_vars, drop=F]
  
  plot_krona.matrix(
    taxonomy,
    magnitude,
    group = group,
    color_values = color_values,
    color_label = color_label,
    ...
  )
}


plot_krona.data.frame = function(classification, ...) {
  plot_krona.matrix(as(classification, 'matrix'), ...)
}

plot_krona.taxonomyTable = function(taxonomy, ...) {
  plot_krona.matrix(as(taxonomy, 'matrix'), ...)
}


plot_krona.matrix = function(classification,
                             magnitude = NULL,
                             group = NULL,
                             output = NULL,
                             display = TRUE,
                             color_values = NULL,
                             color_label = NULL,
                             color_value_range = NULL,
                             hue_range = c(0, 120),
                             unknown_label = NULL,
                             root_label = 'Root',
                             total_label = 'Total',
                             summary_fn = sum,
                             color_summary_fn = weighted.mean,
                             kronatools_dir = NULL,
                             resources_url = NULL,
                             snapshot = NULL,
                             snapshot_format = 'png',
                             snapshot_res = 220,
                             snapshot_dim = c(6, 5),
                             snapshot_delay = 0.8,
                             snapshot_fontsize = 7,
                             snapshot_chart_size = 0.7,
                             iframe_width = '100%',
                             iframe_height = '470px',
                             krona_opts = c(),
                             minify = T,
                             verbose = F,
                             debug = F)
{
  stopifnot(is.logical(display) && length(display) == 1)
  stopifnot(is.character(root_label) && length(root_label) == 1)
  stopifnot(is.character(total_label) && length(total_label) == 1)
  stopifnot(is.function(summary_fn) && length(summary_fn) == 1)
  stopifnot(is.function(color_summary_fn) && length(color_summary_fn) == 1)
  stopifnot(is.numeric(snapshot_res) && length(snapshot_res) == 1)
  stopifnot(is.numeric(snapshot_dim) && length(snapshot_dim) == 2)
  stopifnot(is.numeric(snapshot_delay) && length(snapshot_delay) == 1)
  stopifnot(is.numeric(snapshot_delay) && length(snapshot_delay) == 1)
  stopifnot(is.numeric(snapshot_fontsize) && length(snapshot_fontsize) == 1)
  stopifnot(is.numeric(snapshot_chart_size) && length(snapshot_chart_size) == 1)
  stopifnot(is.null(resources_url) || is.character(resources_url) && length(resources_url) == 1)
  stopifnot(is.character(iframe_width) && length(iframe_width) == 1)
  stopifnot(is.character(iframe_height) && length(iframe_height) == 1)
  stopifnot(is.logical(display) && length(display) == 1)
  stopifnot(is.logical(verbose) && length(verbose) == 1)
  stopifnot(is.logical(minify) && length(minify) == 1)
  stopifnot(length(krona_opts) == 0 || is.character(krona_opts))
  
  if (is.null(kronatools_dir)) {
    r = suppressWarnings(try(system2(
      'ktImportXML',
      stdout=F, stderr=F, wait=T
    ), silent = T))
    if (r == 0) {
      bin = 'ktImportXML'
    } else {
      stop('plot_krona: ktImportXML not found in path. Is Krona installed?')
    }
  } else {
    stopifnot(is.character(kronatools_dir) && length(kronatools_dir) == 1)
    bin = file.path(kronatools_dir, 'scripts', 'ImportXML.pl')
    if (!file.exists(bin)) {
      stop('plot_krona: kronatools_dir supplied, but KronaTools not found in the directory')
    }
  }
  
  # where is rendering happending, and for which format?
  try_run = function(...) suppressWarnings(try(..., silent=T))
  # FIXME: not sure if this is the best way to detect these settings
  outformat = try_run(knitr::opts_knit$get('rmarkdown.pandoc.to'))
  in_rstudio = is.null(outformat)
  html_out = isTRUE(try_run(knitr::is_html_output())) | 
    !is.null(outformat) && grepl('gfm', outformat)  # gfm = github_document
  in_console = !isTRUE(try_run(rstudioapi::getActiveDocumentContext()$id != "#console")) &
      is.null(outformat) |
    isTRUE(try_run(endsWith(rstudioapi::getActiveDocumentContext()$path, '.R')))
  
  if (!is.null(output)) {
    ext = tail(strsplit(output, '.', fixed=T)[[1]], 1)
    if (ext %in% c('png', 'pdf', 'jpeg')) {
      if (!is.null(snapshot) && snapshot == F) {
        stop(paste('plot_krona: The extension of the output file is png/pdf/jepg,',
                   'but no snapshot will be taken. Change the ending to .html or',
                   'set snapshot=TRUE'))
      } else if (verbose && is.null(snapshot) && (html_out || in_rstudio)) {
        message(sprintf('Output file has "%s" extension, therefore taking a snapshot', ext))
      }
      snapshot = T
      snapshot_format = ext
      outfile = gsub('\\.[^\\.]+$', '.html', output, perl=T)
      output = NULL  # ensure HTML file is deleted
    } else {
      outfile = output
    }
  } else if (!display) {
    stop('plot_krona: display = FALSE, but no ouptput file specified')
  } else if (in_rstudio & !in_console) {
    # FIXME: while running chunks and clicking into the console,
    #   in_console will become TRUE and the browser will open
    outfile = tempfile(fileext = '.krona.html')
  } else {
    doc = try_run(knitr::current_input(dir = T))
    outdir = if (!is.null(doc)) {
      # this directory should be automatically cleaned up after rendering
      file.path(gsub('\\.[^\\.]*$', '_files', doc, perl = T), 'krona')
    } else {
      # if in console, we just specify this dir (not cleaned up)
      'krona_files'
    }
    dir.create(outdir, F, T)
    # determine a unique output file location, where neither HTML nor 
    # snapshot images are overwritten to make sure that chunks don't
    # interfere with each other
    i = 1
    while (T) {
      f = file.path(outdir, paste0('krona_', i, c('.html', '.png', '.pdf', '.jpg')))
      if (all(!file.exists(f))) {
        outfile = f[1]
        break
      }
      i = i + 1
    }
  }
  
  if (is.null(snapshot)) {
    snapshot = !html_out & !in_rstudio
  }
  stopifnot(is.logical(snapshot) && length(snapshot) == 1)
  
  if (verbose) {
    cat(sprintf('Document format: %s (HTML: %s), RStudio: %s, Console: %s', 
                if (is.null(outformat)) 'N/A' else outformat,
                html_out, in_rstudio, in_console), 
        sep='\n', file=stderr())
    if (!is.null(output)) {
      cat(paste0('Writing output to ', outfile), sep='\n', file=stderr())
    }
  }
  
  
  # prepare/validate classifications
  
  if (!is.null(unknown_label)) {
    stopifnot(is.character(unknown_label) && length(unknown_label) == 1)
    classification[is.na(classification)] = unknown_label
  }
  
  # prepare/validate magnitudes
  
  if (is.null(magnitude)) {
    # nothing supplied: use equal magnitude for all
    magnitude = cbind(setNames(rep(1, nrow(classification)), rownames(classification)))
  }
  stopifnot(all(apply(magnitude, 2, is.numeric)))
  if (nrow(magnitude) != nrow(classification) || 
      !is.null(rownames(magnitude)) && !is.null(rownames(classification)) &&
        !all.equal(rownames(magnitude), rownames(classification))) {
    stop(paste("The taxa in magnitude matrix and classifications don't match;",
               "their number should be equal and row names are expected to be absent or equal in both."))
  }
  
  magnitude = as(magnitude, 'matrix')

  if (is.null(group)) {
    group = rep(1L, ncol(magnitude))
  }
  if (!(length(group) == 1 && group == 'all')) {
    if (ncol(magnitude) != length(group) ||
        !is.null(colnames(magnitude)) && !is.null(names(group)) &&
          !all.equal(colnames(magnitude), names(group))) {
      stop(paste("The groups and the columns of the magnitude matrix don't match;",
                 "their number should be equal and names are expected to be absent or equal in both."))
    }
    l = lapply(split(1:ncol(magnitude), group),
           function(i) apply(magnitude[, i, drop = F], 1, summary_fn))
    magnitude = simplify2array(l, except=NA)
  }

  if (any(is.na(color_values))) {
    stop("Invalid values (NAs) detected in 'color_values', which cannot be handled.")
  }
  
  # if all names are equal in first column, use first column as root
  if (length(unique(na.omit(classification[, 1]))) == 1 && ncol(classification) > 1) {
    root_label = classification[1, 1]
    classification = classification[, 2:ncol(classification)]
  }
  
  # write XML file in format for KronaTools 2.0
  
  # functions for formatting values
  float_fn = function(example_vec) {
    digits = max(1, round(log10(1 / max(
      example_vec, na.rm = T
    ))) + 2)
    function(x)
      formatC(x, format = 'f', digits = digits)
  }
  
  get_format_fn = function(example) {
    if (isTRUE(all.equal(example, as.integer(example))))
      function(x)
        as.character(x)
    else
      float_fn(example)
  }
  value_fn = get_format_fn(magnitude)
  color_fn = if (!is.null(color_values)) {
    float_fn(as.vector(color_values))
  } else { NULL }
  
  write_xml_node = function(name, tax, abund, color, file) {
    n = colSums(abund)
    cat(
      sprintf('<node name="%s">', name),
      sprintf('<n>%s</n>', paste(
        sprintf('<v>%s</v>', ifelse(n == 0, '', value_fn(n))), collapse = ''
      )),
      sep = '\n',
      file = file
    )
    if (!is.null(color)) {
      col = apply(abund, 2, function(n)
        color_summary_fn(color, n))
      cat(sprintf('<c>%s</c>', paste(
        sprintf('<v>%s</v>', ifelse(is.nan(col), 0, color_fn(col))), collapse =
          ''
      )),
      sep = '\n',
      file = file)
    }
    if (!is.null(tax)) {
      # write child nodes
      for (val in unique(na.omit(tax[, 1]))) {
        sel = which(val == tax[,1] & !is.na(tax[,1]))
        tsub = if (ncol(tax) > 1) {
          tax[sel, 2:ncol(tax), drop = F]
        } else {
          NULL
        }
        write_xml_node(tax[sel[1], 1], tsub, abund[sel, , drop = F], color[sel], file)
      }
    }
    cat('</node>', sep = '\n', file = file)
  }
  
  xml_file = tempfile(fileext = ".krona.xml")
  f = file(xml_file, "w")  # open an output file connection
  cat('<krona>', sep = '\n', file = f)
  cat(
    '<attributes magnitude="n">',
    sprintf('<attribute display="%s">n</attribute>', total_label),
    sep = '\n',
    file = f
  )
  if (!is.null(color_values)) {
    if (is.null(color_label)) {
      stop("plot_krona: 'color_label' must be supplied along with 'color_values'")
    }
    stopifnot(is.character(color_label) && length(color_label) == 1)
    cat(
      sprintf('<attribute display="%s">c</attribute>', color_label),
      sep = '\n',
      file = f
    )
  }
  cat('</attributes>', sep = '\n', file = f)
  if (!is.null(color_values)) {
    if (is.null(color_value_range))
      color_value_range = range(color_values)
    stopifnot(length(color_value_range) == 2)
    stopifnot(!is.null(hue_range) && length(hue_range) == 2)
    cat(
      sprintf(
        '<color attribute="c" valueStart="%s" valueEnd="%s" hueStart="%s" hueEnd="%s"></color>',
        color_value_range[1],
        color_value_range[2],
        hue_range[1],
        hue_range[2]
      ),
      sep = '\n',
      file = f
    )
  }
  if (ncol(magnitude) > 1) {
    d =
      cat(
        '<datasets>',
        paste(sprintf(
          '<dataset>%s</dataset>', colnames(magnitude)
        ), collapse = '\n'),
        '</datasets>',
        sep = '\n',
        file = f
      )
  }
  write_xml_node(root_label, classification, magnitude, color_values, file = f)
  cat('</krona>', sep = '\n', file = f)
  close(f)
  
  # generate Krona chart from XML
  args = c(
    '-o',
    outfile,
    if (!is.null(resources_url)) c( '-u', resources_url) else NULL,
    xml_file
  )
  # print(paste(c(bin, args, '-o', outfile, xml_file), collapse=' '))
  
  ret = system2(bin, args, stdout=TRUE, stderr=TRUE)
  status = attr(ret, 'status')
  if (!is.null(status) && status != 0) {
    stop(ret)
  }
  file.remove(xml_file)
  
  # Return the output
  
  make_url_params = function(opts) {
    paste(
      sapply(names(opts), function(p) paste(p, opts[p], sep = '=')), 
      collapse = '&'
    )
  }
  
  make_file_url = function(path, opts) {
    sprintf('file:///%s?%s', normalizePath(path), make_url_params(opts))
  }
  
  # activate custom coloring if necessary
  if (!is.null(color_values) && !isTRUE(krona_opts['color'] == 'false')) {
    krona_opts['color'] = 'true'
  }
  
  if (snapshot) {
    if (!('webshot' %in% rownames(installed.packages()))) {
      stop(
        paste(
          "Please install webshot (https://github.com/wch/webshot)",
          "for displaying Krona widgets in PDF documents."
        )
      )
    }
    if (!webshot::is_phantomjs_installed()) {
      stop('PhantomJS not installed, please run webshot::install_phantomjs()')
    }
    
    # modify JavaScript code to automatically take a snapshot after loading
    # tested with version 2.8.1
    html = readChar(outfile, file.info(outfile)$size)
    if (is.null(output))
      file.remove(outfile)
    
    # modify snapshot() function to not open in new window
    html = gsub('var +\\w+ += +window.open\\(\\) *; *\r?\n *\\w+.document.write', 
                'document.write', html, perl=T)
    # hide download link from snapshot
    html = gsub('Download Snapshot', '', html, ignore.case = T)
    # make sure the snapshot is directly taken after loading
    inject_code = sprintf(
      'load(); tweenLength=0; bufferFactor=%.3f; setTimeout(snapshot, %.0f);',
      .1 / snapshot_chart_size, snapshot_delay * 1000
    )
    html = gsub('window\\.onload *= *load *',
                sprintf('window.onload = function() { %s };', inject_code),
                html,
                ignore.case = T)
    # write the modified version to file
    prefix = gsub('\\.[^\\.]+$', '', outfile, perl = T)
    snap_html_out = paste0(prefix, '.krona.snapshot.html')
    writeChar(html, snap_html_out)
    snap = paste0(prefix, '.', snapshot_format)
    
    # Krona charts appears to use 11px font size for the 11pt font,
    # so we calculate the font size in px required to obtain the correct
    # size in pt at final resolution
    krona_opts['font'] = round(snapshot_fontsize * snapshot_res / 72.72, 1)
    
    url = make_file_url(snap_html_out, krona_opts)
    snapshot_dim = round(snapshot_dim * snapshot_res)
    
    # due to https://github.com/wch/webshot/issues/115
    Sys.setenv(OPENSSL_CONF="/dev/null")
    webshot::webshot(
      url,
      snap,
      delay = snapshot_delay,
      vwidth = snapshot_dim[1],
      vheight = snapshot_dim[2],
      debug = debug,
    )
    file.remove(snap_html_out)
    
    if (verbose) {
      cat(paste0('Created snapshot: ', snap), sep='\n', file=stderr())
    }
    
    if (display) {
      if (in_console) {
        browseURL(snap)
      } else {
        return(knitr::include_graphics(snap, dpi = snapshot_res, auto_pdf = F, rel_path = F))
      }
    }
  }
  
  # ... or embed / browse the output
  if (display) {
    if (in_console) {
      return(browseURL(outfile))
    }
    if (!('htmltools' %in% rownames(installed.packages()))) {
      stop("plot_krona: the htmltools package is required for HTML output.")
    }
    if (in_rstudio | html_out & minify) {
      if (!('js' %in% rownames(installed.packages()))) {
        stop("plot_krona: the 'js' package is required with 'minify=TRUE'")
      }
      # in RStudio: use 'srcdoc' to embed whole document
      html = readChar(outfile, file.info(outfile)$size)
      # file is not needed anymore
      if (is.null(output))
        file.remove(outfile)
      # minify the included scripts
      if (is.null(resources_url)) {
        js_pos = gregexec('< *script[^>]*>(.*?)(< */ *script *>)', html)[[1]][2:3,,drop=F]
        for (num in ncol(js_pos):1) {
          start = js_pos[1, num]
          end = js_pos[2, num] - 1
          jscript = substr(html, start, end)
          # modify JS code to make sure URL parameters are parsed even if there is no URL
          jscript = gsub('document.location', sprintf('"url?%s"', make_url_params(krona_opts)), jscript, fixed=T)
          # minify the script
          jscript = suppressWarnings(js::uglify_optimize(jscript))
          html = paste0(substr(html, 1, start-1), jscript, substr(html, end+1, nchar(html, type='bytes')))
        }
      }
      
      # remove line breaks from node hierarchy
      html = gsub('(</?(krona|attributes|attribute|color|dataset|node|n|c)[^>]*>)\n', '\\1', html)
      # return the IFrame
      iframe = htmltools::tags$iframe(
        srcdoc = html,
        width = iframe_width,
        height = iframe_height,
        scrolling = 'no',
        seamless = 'seamless',
        frameBorder = '0'
      )
    } else  if (html_out) {
      # Render to HTML without minifying:
      # Use the src attribute instead (file is then encoded as text/html)
      # Minifying removes all the line breaks in the Javascript code, while the
      # presence of many lines makes pandoc slow
      # FIXME: src may have size limitation
      iframe = htmltools::tags$iframe(
        src = outfile,
        width = iframe_width,
        height = iframe_height,
        scrolling = 'no',
        seamless = 'seamless',
        frameBorder = '0'
      )
    } else {
      stop('plot_krona: snapshot = FALSE is not possible in static document formats')
    }
    if (html_out) {
      # We create a new class with an associated 'knit_print' method
      # (see below), which finally attaches a 'knit_asis_htmlwidget'
      # class to the result. This is all done to obtain a proper figure caption.
      # FIXME: is there a better way to allow captions? This strategy exploits
      #  an implementation detail, which may change (see sew.knit_asis in output.R)
      return(structure(list(iframe), class='embedded_krona_html_widget'))
    }
    return(iframe)
  }
}


knit_print.embedded_krona_html_widget = function(x, ...) {
  x = htmltools::knit_print.shiny.tag(x[[1]], ...)
  class(x) = c(class(x), 'knit_asis_htmlwidget')
  x
}
