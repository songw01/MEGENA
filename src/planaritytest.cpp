// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace Rcpp;
using namespace boost;

// [[Rcpp::export]]
SEXP planaritytest(SEXP N,SEXP rows,SEXP cols){


    NumericVector rr(rows);
    NumericVector cc(cols);
    int nv = INTEGER(N)[0];

    typedef adjacency_list<vecS,vecS,undirectedS> graph;
    graph g(nv);

    int ne = rr.size();

    for (int i = 0;i < ne;i++)
    {
       add_edge(rr[i],cc[i],g);
    }

    bool ptest = boyer_myrvold_planarity_test(g);

    return wrap(ptest);
}
