#' Show the details about RNASeqRParam S4 obejct
#'
#'@name show
#'@exportMethod show
#' @docType methods
#' @rdname eexpr
#' @param x a RNASeqRParam object
#' @return None
#'
#' @examples
#' data(yeast)
#' show(yeast)
setMethod('show', 'RNASeqRParam', function(object){
  cat("RNASeqRParam S4 object\n",
      "             os.type :", object@os.type, "\n",
      "     python.variable : (Availability:",
      object@python.variable$check.answer,
      ", Version:", object@python.variable$python.version,")\n",
      "         python.2to3 :", object@python.2to3, "\n",
      "         path.prefix :", object@path.prefix, "\n",
      "   input.path.prefix :", object@input.path.prefix, "\n",
      "         genome.name :", object@genome.name, "\n",
      "      sample.pattern :", object@sample.pattern,"\n",
      "independent.variable :", object@independent.variable,"\n",
      "          case.group :", object@case.group,"\n",
      "       control.group :", object@control.group,"\n",
      "    indices.optional :", object@indices.optional,"\n",
      "independent.variable :", object@independent.variable,"\n"
      )
})
