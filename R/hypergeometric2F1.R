hypergeometric2F1 = function(a,b,c,z, method="Cephes", log=TRUE) {
  out = 1.0
  if (c < b  | b < 0) {
    warning("Must have c > b > 0 in 2F1 function for integral to converge")
    return(log(0))}
  if (abs(z) > 1) {
    warning("integral in 2F1 diverges")
    return(out)}
  if (z == 1.0  & c - b - a < 0 ) ans = Inf
  else {
      if (method == "Laplace") {
        ans = .C("logHyperGauss2F1", as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(z), out=as.numeric(out), PACKAGE="BAS")$out 
        if (!log) ans= exp(ans)
      }
      else {
        ans = .C("hypergeometric2F1", as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(z), out=as.numeric(out), PACKAGE="BAS")$out
        if (ans < 0) warning("Cephes routine is negative; try Laplace approximation")
        if (log) ans = log(ans)
      }
    }
  return(ans)
}
