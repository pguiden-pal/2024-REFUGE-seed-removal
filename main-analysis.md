Part 1: Seed removal responses to management
================
Pete Guiden
2024-12-05

This code provides a preliminary inspection of seed removal data
collected at 19 ReFUGE plots at Nachusa Grasslands in July and October
2024. 8 seeds of 7 species (list here) were placed cafeteria-style in
buckets that allowed access to either invertebrates alone, or
invertebrates and rodents. Plots vary in their management history
(restoration age, bison presence, burn status for 2024).

Expect this code to be the main data presented in a manuscript examining
how management influences seed removal of native versus invasive plant
species.

## Inspecting the data

First step is to inspect the data before running any analyses.

There is one tricky things about the data that should be mentioned:
while we only put 8 seeds per species in each bucket, some buckets had
\>8 seeds of a species found. This is a low proportion (n = 14 cases out
of 531 species-bucket combinations, or 2.6%). My current resolution for
this was to assume that these cases had no seed removal of that species,
and to recalculate the total number of starting seeds to be whatever the
number of intact seeds recovered was. Should be a conservative way to
estimate seed removal without losing data, although I suspect running
models with and without these species/bucket combos will yield very
similar outcomes.

We can also check out some quick and dirty figures summarizing the data.
Looks like no issues with zero-inflation, etc., so a binomial glmer()
should work ok.

``` r
summary(seed.data)
```

    ##        ID          Planting            Month            Treatment        
    ##  Min.   :  1.0   Length:531         Length:531         Length:531        
    ##  1st Qu.:133.5   Class :character   Class :character   Class :character  
    ##  Median :266.0   Mode  :character   Mode  :character   Mode  :character  
    ##  Mean   :266.0                                                           
    ##  3rd Qu.:398.5                                                           
    ##  Max.   :531.0                                                           
    ##                                                                          
    ##    Species          Number.of.seeds.intact Number.of.damaged.seeds
    ##  Length:531         Min.   : 0.0           Min.   :1.00           
    ##  Class :character   1st Qu.: 5.0           1st Qu.:1.00           
    ##  Mode  :character   Median : 7.0           Median :1.00           
    ##                     Mean   : 6.2           Mean   :1.86           
    ##                     3rd Qu.: 8.0           3rd Qu.:2.00           
    ##                     Max.   :14.0           Max.   :8.00           
    ##                                            NA's   :424            
    ##  Evidence.of.damaged.seeds.    Notes            total.seeds     intact.seeds 
    ##  Length:531                 Length:531         Min.   : 8.00   Min.   : 0.0  
    ##  Class :character           Class :character   1st Qu.: 8.00   1st Qu.: 5.0  
    ##  Mode  :character           Mode  :character   Median : 8.00   Median : 7.0  
    ##                                                Mean   : 8.04   Mean   : 6.2  
    ##                                                3rd Qu.: 8.00   3rd Qu.: 8.0  
    ##                                                Max.   :14.00   Max.   :14.0  
    ##                                                                              
    ##  removed.seeds  prop.intact.seeds   invasive            Bison          
    ##  Min.   :0.00   Min.   :0.000     Length:531         Length:531        
    ##  1st Qu.:0.00   1st Qu.:0.625     Class :character   Class :character  
    ##  Median :1.00   Median :0.875     Mode  :character   Mode  :character  
    ##  Mean   :1.84   Mean   :0.770                                          
    ##  3rd Qu.:3.00   3rd Qu.:1.000                                          
    ##  Max.   :8.00   Max.   :1.000                                          
    ##                                                                        
    ##   Burn.2024         Year.restored       lat             lon        
    ##  Length:531         Min.   :1986   Min.   :41.87   Min.   :-89.37  
    ##  Class :character   1st Qu.:2002   1st Qu.:41.88   1st Qu.:-89.36  
    ##  Mode  :character   Median :2007   Median :41.89   Median :-89.33  
    ##                     Mean   :2005   Mean   :41.89   Mean   :-89.34  
    ##                     3rd Qu.:2009   3rd Qu.:41.90   3rd Qu.:-89.33  
    ##                     Max.   :2013   Max.   :41.90   Max.   :-89.31  
    ##                     NA's   :56                                     
    ##     rest.age    
    ##  Min.   :11.00  
    ##  1st Qu.:15.00  
    ##  Median :17.00  
    ##  Mean   :19.25  
    ##  3rd Qu.:21.50  
    ##  Max.   :37.50  
    ##  NA's   :56

``` r
# How many buckets had more than 8 seeds recovered per species?
length(filter(seed.data, intact.seeds > 8)$ID)
```

    ## [1] 14

``` r
length(filter(seed.data, intact.seeds > 8)$ID)/length(seed.data$ID)
```

    ## [1] 0.02636535

``` r
ggplot(seed.data, aes(x = Species, fill = Treatment, y = prop.intact.seeds))+
  facet_wrap(~Month)+
  geom_boxplot()
```

![](main-analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(seed.data, aes(x = prop.intact.seeds, fill = Treatment))+
  facet_wrap(Month~Species)+
  geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](main-analysis_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

## Model construction

After inspecting the data, we can run some models.

Before modeling for real, I like to run an lmer() strictly to check the
denominator DF’s to make sure the model is accurately interpreting our
experimental design (following the advice of Arnqvist 2020 TREE). For
example, any management action should be replicated at the planting
level, not the bucket level (so it should have ~18 DDFs instead of ~500
DDFs). After some tinkering, I ended up with a random intercept term
nesting Treatment inside Planting (i.e., (1\|Planting/Treatment)), which
makes sense. I also add a random intercept for the unit of replication
(i.e., species-planting-month combination, which matches the ID column;
(1\|ID)), which helps account for overdispersion in binomial/Poisson
models.

It bears repeating but *we don’t care about the coefficients or P-values
here*; this exercise is simply done to eyeball the denominator degrees
of freedom to increase our confidence about random effects structures.

``` r
# Check with regular lmer: Treatment is nested within Planting!
seed.model.check = lmer(data = seed.data,
                   prop.intact.seeds~Treatment*invasive+Month+Bison*Burn.2024+
                     (1|Planting/Treatment)+
                     (1|Species))
```

    ## boundary (singular) fit: see help('isSingular')

``` r
Anova(seed.model.check, test = 'F')
```

    ## Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
    ## 
    ## Response: prop.intact.seeds
    ##                         F Df Df.res Pr(>F)  
    ## Treatment          2.1068  1  18.00 0.1639  
    ## invasive           0.0025  1   5.00 0.9622  
    ## Month              1.6368  1 485.05 0.2014  
    ## Bison              2.6188  1  15.03 0.1264  
    ## Burn.2024          0.7984  1  15.03 0.3856  
    ## Treatment:invasive 0.0802  1 485.06 0.7772  
    ## Bison:Burn.2024    5.9715  1  15.08 0.0273 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

After the random effects are structure properly, we can get into the
real modeling. Since we’re interested in knowing the proportion of seeds
removed, we’ll use a binomial GLMM withn the cbind() function. Note that
the order of terms in cbind() matters; here we’ll go cbind(removed,
intact) to express results as the **proportion of seeds removed**. We’ll
use this convention throughout.

As far as structuring the model, we have some decisions to make. We have
six key variables we are hypothesizing to be important: age (continuous,
no remnants), bison (present/absent), fire (burned/unburned), treatment
(inverts, rodents+inverts), month (Jul, Oct), and species (7 levels).
**We will keep all of these main effects in the model**. However, we
also have reason to think that interactions might be important (e.g.,
bison might affect rodents and inverts differently). We’ll use model
selection (MuMIn package) to tell us which interactions are the most
important, starting with a global model with….6 main effects and 15(!!!)
interactive effects. Unsurprisingly, the global model has two key
weaknesses: 1) model diagnostics suggest problems (huge VIF, etc.) and
2) model convergence requires a lot of computer power (~90 minutes on my
machine). Code for this model selection process is hidden in a code
chunk below, but is commented out because it takes \>3 hours to perform.
We decided *a priori* in a Dec 5 meeting to retain the simplest top
model. R returned 4 models with \<2 AICc, and no clear winner.
Therefore, we’ll keep the simplest of the “best” 4 models; it has the
added bonus of slightly nicer model diagnostics, too.

    ## boundary (singular) fit: see help('isSingular')

![](main-analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(removed.seeds, intact.seeds) ~ Bison + Burn.2024 + Month +  
    ##     invasive + Treatment + rest.age + Bison * Burn.2024 + Bison *  
    ##     Treatment + Burn.2024 * Treatment + invasive * Month + Month *  
    ##     Treatment + (1 | Planting/Treatment) + (1 | Species) + (1 |      ID)
    ##    Data: seed.data
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1669.9   1736.5   -818.9   1637.9      459 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -0.98148 -0.55899 -0.07241  0.32112  1.46534 
    ## 
    ## Random effects:
    ##  Groups             Name        Variance  Std.Dev. 
    ##  ID                 (Intercept) 1.987e+00 1.409e+00
    ##  Treatment:Planting (Intercept) 6.771e-02 2.602e-01
    ##  Planting           (Intercept) 9.821e-15 9.910e-08
    ##  Species            (Intercept) 2.431e-02 1.559e-01
    ## Number of obs: 475, groups:  
    ## ID, 475; Treatment:Planting, 34; Planting, 17; Species, 7
    ## 
    ## Fixed effects:
    ##                                                Estimate Std. Error z value
    ## (Intercept)                                    -2.96864    0.69093  -4.297
    ## Bisonno bison                                   1.11491    0.38783   2.875
    ## Burn.2024unburned                               1.66189    0.52701   3.153
    ## MonthOctober                                   -1.79337    0.31853  -5.630
    ## invasivenative                                 -0.58760    0.26646  -2.205
    ## TreatmentRodents and Inverts                   -1.63361    0.40748  -4.009
    ## rest.age                                        0.05973    0.02182   2.738
    ## Bisonno bison:Burn.2024unburned                -2.44669    0.52954  -4.620
    ## Bisonno bison:TreatmentRodents and Inverts      1.21169    0.41170   2.943
    ## Burn.2024unburned:TreatmentRodents and Inverts  1.28474    0.38990   3.295
    ## MonthOctober:invasivenative                     1.34812    0.34364   3.923
    ## MonthOctober:TreatmentRodents and Inverts       1.36761    0.33887   4.036
    ##                                                Pr(>|z|)    
    ## (Intercept)                                    1.73e-05 ***
    ## Bisonno bison                                  0.004044 ** 
    ## Burn.2024unburned                              0.001614 ** 
    ## MonthOctober                                   1.80e-08 ***
    ## invasivenative                                 0.027442 *  
    ## TreatmentRodents and Inverts                   6.10e-05 ***
    ## rest.age                                       0.006185 ** 
    ## Bisonno bison:Burn.2024unburned                3.83e-06 ***
    ## Bisonno bison:TreatmentRodents and Inverts     0.003249 ** 
    ## Burn.2024unburned:TreatmentRodents and Inverts 0.000984 ***
    ## MonthOctober:invasivenative                    8.74e-05 ***
    ## MonthOctober:TreatmentRodents and Inverts      5.44e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Bsnnbs B.2024 MnthOc invsvn TrtRaI rest.g Bb:B.2 Bb:TaI
    ## Bisonnbison -0.721                                                        
    ## Brn.2024nbr -0.717  0.526                                                 
    ## MonthOctobr -0.203  0.007  0.000                                          
    ## invasiventv -0.192 -0.011 -0.018  0.373                                   
    ## TrtmntRdnaI -0.253  0.296  0.082  0.227  0.016                            
    ## rest.age    -0.877  0.541  0.653 -0.001 -0.022 -0.006                     
    ## Bbsn:B.2024  0.642 -0.641 -0.823  0.000  0.019  0.057 -0.608              
    ## Bbsn:TrtRaI  0.169 -0.494  0.071 -0.016 -0.017 -0.669  0.004 -0.053       
    ## B.2024:TRaI  0.082  0.070 -0.361 -0.003 -0.003 -0.361  0.009 -0.034 -0.066
    ## MnthOctbr:n  0.116  0.007  0.007 -0.642 -0.624 -0.015  0.017 -0.012  0.014
    ## MnthOc:TRaI  0.123 -0.009 -0.002 -0.573 -0.007 -0.435 -0.006  0.004  0.031
    ##             B.20aI MnthO:
    ## Bisonnbison              
    ## Brn.2024nbr              
    ## MonthOctobr              
    ## invasiventv              
    ## TrtmntRdnaI              
    ## rest.age                 
    ## Bbsn:B.2024              
    ## Bbsn:TrtRaI              
    ## B.2024:TRaI              
    ## MnthOctbr:n  0.009       
    ## MnthOc:TRaI  0.015  0.032
    ## optimizer (bobyqa) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: cbind(removed.seeds, intact.seeds)
    ##                       Chisq Df Pr(>Chisq)    
    ## Bison                5.8198  1  0.0158468 *  
    ## Burn.2024            1.9036  1  0.1676798    
    ## Month                3.1925  1  0.0739750 .  
    ## invasive             0.0953  1  0.7575851    
    ## Treatment            4.1833  1  0.0408235 *  
    ## rest.age             7.4956  1  0.0061851 ** 
    ## Bison:Burn.2024     21.3482  1  3.830e-06 ***
    ## Bison:Treatment      8.6619  1  0.0032493 ** 
    ## Burn.2024:Treatment 10.8575  1  0.0009839 ***
    ## Month:invasive      15.3907  1  8.742e-05 ***
    ## Month:Treatment     16.2880  1  5.441e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Model output

Ok, lots to unpack here! We’ll use the emmeans package to visualize key
interactive effects, including post-hoc tests to determine significance
between categorical groups. There are five significant interactions to
examine:

### 1) Bison x Fire

Averaged across all treatments, plantings that both were burned and
allowed bison access had the lowest amount of seed removal. Some
interesting tidbits: -Looking only at plots with bison access, seed
removal was **7 times lower** in burned vs unburned plots -Looking only
at unburned plots, there’s no difference in seed removal between bison
absent and bison present plots. This suggests that you really might want
both management practices if minimizing seed removal is important.

    ## $emmeans
    ##  Bison    Burn.2024   prob     SE  df asymp.LCL asymp.UCL
    ##  bison    burned    0.0412 0.0116 Inf    0.0235    0.0711
    ##  no bison burned    0.1936 0.0268 Inf    0.1465    0.2516
    ##  bison    unburned  0.3009 0.0686 Inf    0.1851    0.4492
    ##  no bison unburned  0.1724 0.0275 Inf    0.1248    0.2332
    ## 
    ## Results are averaged over the levels of: Month, invasive, Treatment 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the logit scale 
    ## 
    ## $contrasts
    ##  contrast                            odds.ratio     SE  df null z.ratio p.value
    ##  bison burned / no bison burned          0.1789 0.0604 Inf    1  -5.098  <.0001
    ##  bison burned / bison unburned           0.0998 0.0491 Inf    1  -4.688  <.0001
    ##  bison burned / no bison unburned        0.2063 0.0745 Inf    1  -4.368  0.0001
    ##  no bison burned / bison unburned        0.5579 0.1930 Inf    1  -1.687  0.3306
    ##  no bison burned / no bison unburned     1.1531 0.2714 Inf    1   0.605  0.9305
    ##  bison unburned / no bison unburned      2.0667 0.7150 Inf    1   2.098  0.1536
    ## 
    ## Results are averaged over the levels of: Month, invasive, Treatment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log odds ratio scale

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](main-analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### 2) Bison x Guild

Now we can start looking at how management affects different granivore
guilds. First up, a significant bison x guild interaction. Key take-home
here is that in the Inverts only depots, there is no difference in seed
removal between plots with and without bison. However, where rodents and
inverts both are able to access seeds, we see 2.5 times more seed
removal in plots without bison compared to plots with bison. I interpret
this as perhaps evidence that there is a behavioral change in rodents
where bison are present (Guiden et al. 2022 Ecology).

``` r
emmeans(small.model, pairwise~Bison*Treatment, type = 'response', adjust = 'Tukey')
```

    ## $emmeans
    ##  Bison    Treatment            prob     SE  df asymp.LCL asymp.UCL
    ##  bison    Inverts             0.137 0.0298 Inf    0.0882     0.206
    ##  no bison Inverts             0.125 0.0206 Inf    0.0894     0.171
    ##  bison    Rodents and Inverts 0.104 0.0242 Inf    0.0656     0.162
    ##  no bison Rodents and Inverts 0.260 0.0337 Inf    0.1995     0.331
    ## 
    ## Results are averaged over the levels of: Burn.2024, Month, invasive 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the logit scale 
    ## 
    ## $contrasts
    ##  contrast                                                 odds.ratio     SE  df
    ##  bison Inverts / no bison Inverts                              1.115 0.3324 Inf
    ##  bison Inverts / bison Rodents and Inverts                     1.360 0.4636 Inf
    ##  bison Inverts / no bison Rodents and Inverts                  0.451 0.1331 Inf
    ##  no bison Inverts / bison Rodents and Inverts                  1.220 0.3723 Inf
    ##  no bison Inverts / no bison Rodents and Inverts               0.405 0.0951 Inf
    ##  bison Rodents and Inverts / no bison Rodents and Inverts      0.332 0.0991 Inf
    ##  null z.ratio p.value
    ##     1   0.364  0.9836
    ##     1   0.902  0.8039
    ##     1  -2.699  0.0351
    ##     1   0.652  0.9148
    ##     1  -3.850  0.0007
    ##     1  -3.694  0.0013
    ## 
    ## Results are averaged over the levels of: Burn.2024, Month, invasive 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log odds ratio scale

``` r
bison.treat.emm <- summary(emmeans(small.model, ~Bison*Treatment, type = 'response'))

ggplot(bison.treat.emm, aes(x = Bison, y = prob, color = Treatment))+
  geom_point(position = position_dodge(width = 0.25), size = 5)+
  geom_point(data = seed.data, aes(y = 1-prop.intact.seeds), shape = 1,
             position = position_jitterdodge(dodge.width = 0.25,
                                             jitter.height = 0.05,
                                             jitter.width = 0.05))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                size = 1.25, width = 0.25,
                position = position_dodge(width = 0.25))+
  labs(x = '', y = 'Proportion of seeds removed')+
  scale_y_continuous(limits = c(-0.05, 1.05))+
  scale_color_manual(values = c('royalblue', 'goldenrod1'))+
  pal.theme
```

![](main-analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### 3) Fire x Guild

Burning seems to increase seed removal, but only by rodents. There is no
difference between Invert-only depots in burned versus unburned
plantings. On the other hand, depots allowing inverts and rodents
experienced 4x greater seed removal in unburned areas vs burned areas
(thatch = vole city?). There’s also a compelling increase in seed
removal between unburned plantings allowing rodents compared to unburned
plantings with only inverts.

``` r
emmeans(small.model, pairwise~Burn.2024*Treatment, type = 'response', adjust = 'Tukey')
```

    ## $emmeans
    ##  Burn.2024 Treatment             prob     SE  df asymp.LCL asymp.UCL
    ##  burned    Inverts             0.1077 0.0202 Inf    0.0740     0.154
    ##  unburned  Inverts             0.1576 0.0351 Inf    0.1003     0.239
    ##  burned    Rodents and Inverts 0.0788 0.0160 Inf    0.0526     0.116
    ##  unburned  Rodents and Inverts 0.3240 0.0545 Inf    0.2274     0.438
    ## 
    ## Results are averaged over the levels of: Bison, Month, invasive 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the logit scale 
    ## 
    ## $contrasts
    ##  contrast                                                  odds.ratio     SE
    ##  burned Inverts / unburned Inverts                              0.645 0.2216
    ##  burned Inverts / burned Rodents and Inverts                    1.411 0.3651
    ##  burned Inverts / unburned Rodents and Inverts                  0.252 0.0856
    ##  unburned Inverts / burned Rodents and Inverts                  2.187 0.7820
    ##  unburned Inverts / unburned Rodents and Inverts                0.390 0.1207
    ##  burned Rodents and Inverts / unburned Rodents and Inverts      0.178 0.0605
    ##   df null z.ratio p.value
    ##  Inf    1  -1.276  0.5780
    ##  Inf    1   1.329  0.5445
    ##  Inf    1  -4.057  0.0003
    ##  Inf    1   2.188  0.1264
    ##  Inf    1  -3.042  0.0126
    ##  Inf    1  -5.086  <.0001
    ## 
    ## Results are averaged over the levels of: Bison, Month, invasive 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log odds ratio scale

``` r
fire.treat.emm <- summary(emmeans(small.model, ~Burn.2024*Treatment, type = 'response'))

ggplot(fire.treat.emm, aes(x = Burn.2024, y = prob, color = Treatment))+
  geom_point(position = position_dodge(width = 0.25), size = 5)+
  geom_point(data = seed.data, aes(y = 1-prop.intact.seeds), shape = 1,
             position = position_jitterdodge(dodge.width = 0.25,
                                             jitter.height = 0.05,
                                             jitter.width = 0.05))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                size = 1.25, width = 0.25,
                position = position_dodge(width = 0.25))+
  labs(x = '', y = 'Proportion of seeds removed')+
  scale_y_continuous(limits = c(-0.05, 1.05))+
  scale_color_manual(values = c('royalblue', 'goldenrod1'))+
  pal.theme
```

![](main-analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### 4) Month x Provenance

Grouping species together as native vs non-native, we see that the
problematic non-native species are almost 3x more likely to be removed
in summer than fall. There is no difference in seed removal between
native species in summer and fall. This could suggest that seed removal
by inverts is key for non-native species as a whole.

Probably want to keep the old species-by-species figure as a supplement
to unpack these patterns.

``` r
emmeans(small.model, pairwise~Month*invasive, type = 'response', adjust = 'Tukey')
```

    ## $emmeans
    ##  Month   invasive   prob     SE  df asymp.LCL asymp.UCL
    ##  July    invasive 0.2253 0.0366 Inf    0.1616     0.305
    ##  October invasive 0.0875 0.0180 Inf    0.0580     0.130
    ##  July    native   0.1391 0.0232 Inf    0.0995     0.191
    ##  October native   0.1702 0.0263 Inf    0.1247     0.228
    ## 
    ## Results are averaged over the levels of: Bison, Burn.2024, Treatment 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the logit scale 
    ## 
    ## $contrasts
    ##  contrast                          odds.ratio    SE  df null z.ratio p.value
    ##  July invasive / October invasive       3.033 0.793 Inf    1   4.244  0.0001
    ##  July invasive / July native            1.800 0.480 Inf    1   2.205  0.1218
    ##  July invasive / October native         1.418 0.373 Inf    1   1.328  0.5451
    ##  October invasive / July native         0.593 0.164 Inf    1  -1.886  0.2341
    ##  October invasive / October native      0.467 0.128 Inf    1  -2.779  0.0279
    ##  July native / October native           0.788 0.175 Inf    1  -1.072  0.7069
    ## 
    ## Results are averaged over the levels of: Bison, Burn.2024, Treatment 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log odds ratio scale

``` r
month.invasive.emm <- summary(emmeans(small.model, ~Month*invasive, type = 'response'))

ggplot(month.invasive.emm, aes(x = invasive, y = prob, color = Month))+
  geom_point(position = position_dodge(width = 0.25), size = 5)+
  geom_point(data = seed.data, aes(y = 1-prop.intact.seeds), shape = 1,
             position = position_jitterdodge(dodge.width = 0.25,
                                             jitter.height = 0.05,
                                             jitter.width = 0.05))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                size = 1.25, width = 0.25,
                position = position_dodge(width = 0.25))+
  labs(x = '', y = 'Proportion of seeds removed')+
  scale_y_continuous(limits = c(-0.05, 1.05))+
  scale_color_manual(values = c('forestgreen', 'sienna'))+
  pal.theme
```

![](main-analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### 5) Month x Guild

And finally, we see some evidence that the magnitude of seed removal by
the two consumer guilds changed over time, but in different directions.
In July, seed removal was equally strong between the two guilds. In
October, more than twice as much seed removal was observed in depots
allowing both rodents and inverts, compared to depots with inverts
alone. Again, this makes complete sense given the rather cold weather we
experienced in mid-October. This also supports the previous point about
native vs. non-native species.

``` r
emmeans(small.model, pairwise~Month*Treatment, type = 'response', adjust = 'Tukey')
```

    ## $emmeans
    ##  Month   Treatment             prob     SE  df asymp.LCL asymp.UCL
    ##  July    Inverts             0.2082 0.0329 Inf    0.1509     0.280
    ##  October Inverts             0.0791 0.0155 Inf    0.0535     0.115
    ##  July    Rodents and Inverts 0.1517 0.0264 Inf    0.1069     0.211
    ##  October Rodents and Inverts 0.1865 0.0299 Inf    0.1348     0.252
    ## 
    ## Results are averaged over the levels of: Bison, Burn.2024, invasive 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the logit scale 
    ## 
    ## $contrasts
    ##  contrast                                               odds.ratio    SE  df
    ##  July Inverts / October Inverts                              3.063 0.754 Inf
    ##  July Inverts / July Rodents and Inverts                     1.470 0.393 Inf
    ##  July Inverts / October Rodents and Inverts                  1.147 0.302 Inf
    ##  October Inverts / July Rodents and Inverts                  0.480 0.133 Inf
    ##  October Inverts / October Rodents and Inverts               0.374 0.101 Inf
    ##  July Rodents and Inverts / October Rodents and Inverts      0.780 0.184 Inf
    ##  null z.ratio p.value
    ##     1   4.544  <.0001
    ##     1   1.441  0.4736
    ##     1   0.521  0.9540
    ##     1  -2.657  0.0394
    ##     1  -3.647  0.0015
    ##     1  -1.053  0.7182
    ## 
    ## Results are averaged over the levels of: Bison, Burn.2024, invasive 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log odds ratio scale

``` r
month.treatment.emm <- summary(emmeans(small.model, ~Month*Treatment, type = 'response'))

ggplot(month.treatment.emm, aes(x = Month, y = prob, color = Treatment))+
  geom_point(position = position_dodge(width = 0.25), size = 5)+
  geom_point(data = seed.data, aes(y = 1-prop.intact.seeds), shape = 1,
             position = position_jitterdodge(dodge.width = 0.25,
                                             jitter.height = 0.05,
                                             jitter.width = 0.05))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                size = 1.25, width = 0.25,
                position = position_dodge(width = 0.25))+
  labs(x = '', y = 'Proportion of seeds removed')+
  scale_y_continuous(limits = c(-0.05, 1.05))+
  scale_color_manual(values = c('royalblue', 'goldenrod1'))+
  pal.theme
```

![](main-analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
