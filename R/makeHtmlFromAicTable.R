#  makeHtmlFromAicTable.r

makeHtmlFromAicTable <- function(table,fprefix='',htmlfilename='tmp.html') {
  #'  write html file with links based on AIC table
  #'
  #' \code{writeHtmlFromAicTable} Writes html file with links based on AIC table
  #'
  #'@param table An AIC table of model results
  #'@param fprefix optional prefix for model filenames
  #'@param htmlfilename name of html file to write (default='tmp.html')
  #'@export
  #'@return number of lines written to file.
  #'
  #'@author Jim Hines
  #'
  #'@seealso \code{\link{createAicTable}}
  #'
  #'@examples
  #'\dontrun{
  #' i=makeHtmlFromAicTable(aictable$table)
  #'}
  #'
  d=paste0('<html><head></head><body><table><tr><th>',colnames(table)[1],'</th>')
  for (j in 2:ncol(table)) d=paste0(d,'<th>',colnames(table)[j],'</th>')
  d=paste0(d,'</tr>')

  for (i in 1:nrow(table)) {
    b=paste0('<tr><td><a href="',fprefix,table[i,1],'.out">',table[i,1],'</td>')
    for (j in 2:ncol(table)) b=paste0(b,'<td>',table[i,j],'</td>')
    b=paste0(b,'</tr>')
    d=c(d,b)
  }
  d=c(d,'</table></body></html>')
  write(d,file=htmlfilename)
  return(length(d))
}
