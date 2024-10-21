###We simulated 200× and 400× Illumina sequencing reads of a 10 kb DNA fragment containing an MEI located between 4500 and 5500 bp. 
###The sequencing read length is 2×150bp, with an insert in between the paired reads following a normal distribution (mean=300bp, standard deviation=100bp). 
###The MEI is a heterozygous mutation with a tissue allele frequency of 1%, equivalent to a frequency of 0.005 in bulk sequencing. 
###In each of the 100,000 simulations, we sampled from a Poisson distribution (λ = 0.005) to determine the total number of reads covering the cells with the somatic MEI. 
###We then estimated the likelihood of having more than two supporting reads covering the MEI junctions, assuming the sequencing reads are uniformly distributed across the 10 kb DNA fragment.

simul <- 100000
DNA <- 10000
depth <- 400
support <- c()
left <- 4500
right <- 5500
overhang <- 30
taf <- 0.002
read_length <- 150

for (i in 1:simul)
   {
    read_num= DNA*depth/ (2*150) 
    n= sum(rpois(read_num, taf/2))
    start <- runif(n, min=0, max=DNA)
    size <- rnorm(n, mean=600, sd=100)
    support[i] = 0
    if (n > 0)
      {
       for (j in 1:n)
        {
          if (start[j] < (left - read_length/2) && ((start[j] + size[j]) > (left + overhang)) )
             {
              support[i] = support[i] + 1
             }
          if (((start[j]+read_length) < right) & ((start[j] + size[j]) > (right + read_length/2)) )
             {
              support[i] = support[i] + 1
             }
        }
      }
    }
sum(support>=2)/simul
