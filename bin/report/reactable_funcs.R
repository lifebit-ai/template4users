#' Necessary MaterialUI dependency to call as function
#' @examples
#' muiDependency()
muiDependency <- function() {
    list(
        # Material UI requires React
        reactR::html_dependency_react(),
        htmlDependency(
            name = "mui",
            version = "5.6.3",
            src = ".",
            script = "material-ui.production.min.js"
        )
    )
}

#' Renders the javascript framework react wrapper in R: reactable using functions
#' specific to our tables for report.Rmd
#'
#' Function that generates the reactable table wrapping JS/react and HTML tools
#' @param data the report table we want to render in report.Rmd
#' @return No return as such - simply generates the reactable table
#' @examples
#' generateReactable(data)
generateReactable <- function(data, ...) {
    # create random tag to attach to table
    element_tag <- as.character(rnorm(1))

    # get ellipsis parameters
    user_params <- list(...)

    # Collects all the numerical columns to include range filters on,
    # the rest get generic reactable filters
    num_cols <- select_if(data, is.numeric)
    # Collects all hyperlink columns to add hyperlinks to
    is_hyperlink <- function(column) {
        any(grepl("^https?://", column)) || any(grepl("<a href=", column))
    }
    hyperlink_cols <- select_if(data, is_hyperlink)

    # Apply colDef for hyperlink columns to run JS customCellLink on
    apply_hyperlink_col_def <- lapply(hyperlink_cols, function(col) {
        col <- colDef(
            cell = JS("customCellLink"),
            html = TRUE
        )
        return(col)
    })
    # Apply numerical Material UI range filter functions on numerical columns
    apply_num_col_def <- lapply(num_cols, function(col) {
        col <- colDef(
            filterMethod = JS("filterRange"),
            filterInput = JS("rangeFilter")
        )
        return(col)
    })

    # Set reactable defaults that can be overwritten with user params
    defaults <- list(
        data = data,
        highlight = TRUE,
        bordered = TRUE,
        searchable = TRUE,
        minRows = 10,
        defaultPageSize = 5,
        resizable = TRUE,
        height = 400,
        defaultColDef = colDef(
            minWidth = 80,
            align = "left",
            filterable = TRUE,
        ),
        theme = flatly(),
        columns = c(
            apply_hyperlink_col_def,
            apply_num_col_def
        ),
        elementId = element_tag
    )
    defaults <- defaults[setdiff(names(defaults), names(user_params))]

    # browsable is used to make specific objects such as filter sliders render as HTML
    # tagList provides the printable tags to render in HTML
    # Each filter is provided per selected numerical column as a tag in the tagList
    # rendered to HTML
    browsable(
        tagList(
            # Load MaterialUI dependency for HTML
            muiDependency(),
            # MaterialUI range filter for min and max values to use in place of text
            # filters for numeric columns - all enclosed in HTML() but using JS framework
            # React with React package MaterialUI.
            # RangFilter const changes and returns the state based on the filter actions
            # The React Element returned is the MaterialUI package slider
            # filterRange const filters the table rows
            tags$script(HTML("const rangeFilter = (column, state) => {
                const range = React.useMemo(() => {
                    let min = Infinity
                    let max = -Infinity
                    let stepcheck = 0.1
                    state.data.forEach(row => {
                        const value = row[column.id]
                        stepcheck = value
                        if (value < min) {
                            min = Math.floor(value)
                        } else if (value > max) {
                            max = Math.ceil(value)
                        }
                    })
                    return [min, max, stepcheck]
                }, [state.data])
                const value = column.filterValue ? column.filterValue : range.slice(0, 2);
                const valueLabel = `${value[0]}...${value[1]}`
                const step = Number.isInteger(range[2]) ? 1 : 0.0001
                return React.createElement(
                    'div',
                    { style: { margin: '0 8px' } },
                    [
                        valueLabel,
                        React.createElement(
                            MaterialUI.Slider,
                            {
                                value: value,
                                min: range[0],
                                max: range[1],
                                step: step,
                                valueLabelDisplay: 'auto',
                                onChange: (e, val) => column.setFilter(val),
                                'aria-label': `Filter ${column.name}`
                            }
                        )
                    ]
                )
            }
            const filterRange = (rows, columnId, filterValue) => {
                const [min, max] = filterValue
                console.log(rows)
                return rows.filter(row => {
                    const value = row.values[columnId]
                    return value >= min && value <= max
                })
            }
            function customCellLink(rows, values) {
                const createLink = (url) => {
                    return `<a href='${url}' target='_blank'>${url}</a>`
                }
                const isHyperlink = (potLink) => {
                    const pattern = /^https?:\\/\\//g
                    return pattern.test(potLink) }
                return isHyperlink(rows.value) ? createLink(rows.value) : rows.value
            }")),
            # Download as csv button
            tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('", element_tag, "')")),
            # reactable table object is also provided in the tagList as an object to render
            # to HTML
            do.call(reactable, c(defaults, user_params)),
        )
    )
}
