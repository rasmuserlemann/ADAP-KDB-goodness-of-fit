# ADAP-KDB-goodness-of-fit

The p-values and the power study was done for the manuscript:
1) CreateData.py creates a data file with a different sample representation. Instead of "human:3" there's "human,human,human" in the csv cell.
2) chi_square_p_values.R takes the new data file and calculates p-values for each row in the data file.
3) JSdivergence_clustering.R clusters the data file with hierarchical clustering and JS divergence and generates a PCoA plot of the result.
4) We then find the marginal distribution of each cluster in 3) and insert the frequencies into PowerStudy/power_study.R as label probabilities. 

To run the code, the data files should be sent separately since they are 100mb+ and too large for gitlab.
