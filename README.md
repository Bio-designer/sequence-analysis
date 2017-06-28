# sequence-analysis



## Gibbs-motif-finder
### Motif-finding by Gibbs Sampling

##### Problem:
Given p strings and a length k, find the most “mutually similar”  substring from each string.

##### Assumption:
1.	Each motif appears exactly once in one sequence
2.	The motif has fixed length 

##### Algorithm:
1.	Start with random motif locations and calculate a motif model
2.	Randomly select a sequence, remove its motif and recalculate temporary model
3.	With temporary model, calculate probability of motif at each position on sequence
4.	Select new position based on this distribution
5.	Update model and iterate

##### Reference:
1.	http://veda.cs.uiuc.edu/courses/fa08/cs466/lectures/Lecture16.pdf
2.	https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gibbs.pdf
3.	https://github.com/zshamsi2/motifFinding
