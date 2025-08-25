make_pipeline_table <- function(df, page_size = 5) {
  reactable(
    df,
    defaultColDef = colDef(
      header = function(value) gsub("_", " ", value, fixed = TRUE),
      cell = function(value) format(value, nsmall = 0),
      align = "center",
      filterable = TRUE,
      sortable = TRUE,
      resizable = TRUE,
      footerStyle = list(fontWeight = "bold")
    ),
    columns = list(
      SampleID = colDef(
        name = "SampleID",
        sticky = "left",
        style = list(borderRight = "1px solid #eee"),
        headerStyle = list(borderRight = "1px solid #eee"),
        align = "left",
        minWidth = 125,
        footer = "Total reads"
      ),
      raw_rc = colDef(
        name = "raw rc",
        footer = function(values) sprintf("%.0f", sum(values))
      ),
      cutadapt_rc = colDef(
        name = "cutadapt rc",
        footer = function(values) sprintf("%.0f", sum(values))
      ),
      final_rc = colDef(
        name = "final rc",
        footer = function(values) sprintf("%.0f", sum(values))
      ),
      per_reads_kept = colDef(name = "per reads retain")
    ),
    searchable = TRUE,
    defaultPageSize = page_size,
    pageSizeOptions = c(5, 10, nrow(df)),
    showPageSizeOptions = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    striped = TRUE,
    compact = FALSE,
    wrap = FALSE,
    showSortable = TRUE,
    fullWidth = TRUE,
    theme = reactableTheme(style = list(fontSize = "0.8em"))
  ) %>%
    reactablefmtr::add_subtitle(
      "Read changes through the LotuS3 pipeline.",
      font_size = 15
    )
}

rm(tmp_react)