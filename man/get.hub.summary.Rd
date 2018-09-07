\name{get.hub.summary}
\alias{get.hub.summary}
\title{
summarize hub information.
}
\description{
hubs in different scales are summarized. 
}
\usage{
get.hub.summary(MEGENA.output)
}
\arguments{
  \item{MEGENA.output}{A list object. The output from "do.MEGENA()".}
}
\details{
returns a data.frame object  
}
\value{
A data.frame object with columns:
 \item{node}{hub gene node names}
 \item{S1,...}{binary vector indicating hubs in each scale}
 \item{frequency}{number of scales that respective gene emerges as hub.}
 \item{scale.summary}{list of scales that respective gene as hub.}
}
\author{
Won-Min Song
}