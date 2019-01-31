Date: 2019-01-29
Title: Null test with start at 80 

Data generating process parameterised with 0.3 versus ctl + 25% increase in seroconversion rate. Clinical arm has no difference.

Results:

Reduces type i error in immunological effect/no clinical effect null-case.




Date: 2019-01-30
Title: Increased tte posterior threshold probability to address type i error in immunological effect/no clinical effect null-case.

One way to do this is to change the threshold that the posterior is assessed. 
Another is to change the number of trials that need to `win`.

Approach used:
1. Corrected a possible bug whereby superiority was assessed even if trial deemed futile.
2. Lifted the prob threshold used to assess whether using all resources is likely to result in a win => increase the futility prob thres for tte
3. Transitioned the posterior threshold for assessing superiority in tte (the threshold is the prob used to assess the posterior)
4. Started clinical later (200) and immuno at 80.






Date: 2019-01-30
Title: Ease up on futility test


post_tte_thresh: 0.975




Date:
Title: What to do about strata - at any given interim a statistical trigger may be reached for all children or for one or more strata.




 
