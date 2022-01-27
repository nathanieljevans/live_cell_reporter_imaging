
# General 

- write methods (while it's fresh)


# Specifics 

---
# 1/18/22
Meeting with Samuel 

- ~~add a results column for the `cell count` number of each mutant~~
- Discuss with shannon: the data normalization that Samuel does may behave somewhat as a batch effect correction. 


---
# 1/11/22 

Meeting with shannon 

- add batch effect check at live cell reporter imaging raw value point, e.g., see if specific batches have higher mean intesities. 
- what quality metrics do we have for each batch? How do we know if a batch is acceptable or bad? 
- will batch effect correction at this early point induce negative values? e.g., values outside of the 0,1 reporter range? 
- What EDA + QC does samuel perform? 
- The major probability assignment changes (mutants that switch call, res->sens or sens->res, are these commonly in the same batch?)


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