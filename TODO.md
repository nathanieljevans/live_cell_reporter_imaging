
---
# 11/10/21

Meeting with Samuel ... 

- ~~add `Burn-in` hyper-parameter to sensitivity assignment~~
- ~~write necessary `batch effect checks`~~ and `QC` 
- ~~Switch from `SVM linear kernel classifier` to `gaussian process` to allow for more conservative predictions (sens/res/neither)~~
- need final results by **Late December or early Jan** 

Later analysis steps ... 
- compare cell viability ~ predicted resistance values 
- can we predict cell viability using cluster proportions? 
- EFM192A and SKBR3 behave differently in their resitance mechanisms .. why is this? 
- what is different? (targeted descriptive analysis)
- use CCLE to compare mutation or expression patterns between EFM/SKB... 
- Focus on HER2/ERB2 and family proteins / EGFR pathway 
- Use all cell lines for dimensionality reduction and clustering to identify cell lines with similar resistance patterns? 
- Use the achilles `gene dependency` data to look for differences in protein dependencies? 
- *Samuel action item:* Come up with list of important genes that may be involved in the difference between SKBR3 and EFM192A resistance patterns. 

NOTE: 
Very challenging to be statistically powered to compare expression patterns between EFM192A and SKBR3 ... HOWEVER ... we can use the common expression values of a gene across ALL cell lines to categorize each genes into up-regulated, down-regulated, or neither. Then we could simply compare *differntial expression* between SKBR3 and EFM192A by simple boolean logic... it would still be a descriptive analysis - but results would be much more meaningful... for instance: 

We would still want to use a targeted search (e.g., just interested in genes that are expected to be relevant) but now we could create a list of genes that have statistically significantly different expression (as measured across all cell lines) between EFM192A and SKBR3. 

** ... could use this same approach with `gene dependency` ... **

--- 