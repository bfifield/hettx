## Internal helper functions for formula parsing.
## These replace the orphaned formula.tools package with base R equivalents.

## Get variable names from the left-hand side of a formula.
## Returns character(0) for one-sided formulas (e.g., ~ A + B).
lhs_vars <- function(formula) {
    if (length(formula) == 2) {
        return(character(0))
    }
    all.vars(formula[[2]])
}

## Get variable names from the right-hand side of a formula.
rhs_vars <- function(formula) {
    if (length(formula) == 2) {
        all.vars(formula[[2]])
    } else {
        all.vars(formula[[3]])
    }
}

## Get all resolved variable names from a formula, expanding dot notation.
## Equivalent to formula.tools::get.vars(formula, data=data).
get_vars <- function(formula, data = NULL) {
    if (!is.null(data) && "." %in% all.vars(formula)) {
        attr(terms(formula, data = data), "term.labels")
    } else {
        all.vars(formula)
    }
}
