# hettx

Implements Fisherian hypothesis testing methods proposed by [Ding, Feller, and Miratrix (JRSS-B, 2016)](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12124) for testing null hypothesis of no treatment effect variation, and the systematic estimation methods proposed by [Ding, Feller, and Miratrix (JASA, 2019)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1407322?journalCode=uasa20) for measuring the degree of systematic treatment effect variation explained by observed covariates.

## TODO
- Add calls to both functions
- Capture the name of the function being used
- Capture covariates being used, include in output object
- Make summary functions that replace print functions (model after difference_in_means from estimatr for detect.idiosyncratic and lm_robust for estimate.systematic)
- Clean up R2 results - summary instead of print()
