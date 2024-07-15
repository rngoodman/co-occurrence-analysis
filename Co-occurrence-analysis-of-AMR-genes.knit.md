---
title: "Co-occurrence analysis of AMR genes"
author: "Richard Goodman"
date: "2023-11-03"
output: html_document
---



This analysis takes a binary table of AMR genes present in a defined set of genomes from [abricate](https://github.com/tseemann/abricate) and creates a co-occurrence matrix and subsequent analysis of the co-occurrence interactions between AMR genes across the genomes. 

# Loading packages


```r
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(readr)
library(cooccur)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
```

# Reading in data 


```r
abricate_output = read_csv("data/all_abricate_773_Ec_Kp_CHL_Resfinder_L60_T90.csv")

drug_class = read_csv("data/resfinder_phenotype_table_edit.csv")

ab_tab = abricate_output

# Make sense of dataset 

head(ab_tab)
```

```
## # A tibble: 6 × 15
##   `#FILE`     SEQUENCE     START    END STRAND GENE  COVERAGE COVERAGE_MAP GAPS 
##   <chr>       <chr>        <dbl>  <dbl> <chr>  <chr> <chr>    <chr>        <chr>
## 1 26141-1-134 .26141_1_1… 158811 160043 -      mdf(… 1-1233/… ===========… 0/0  
## 2 26141-1-134 .26141_1_1…    558   1346 -      aadA… 1-789/7… ===========… 0/0  
## 3 26141-1-134 .26141_1_1…   1477   1950 -      dfrA… 1-474/4… ===========… 0/0  
## 4 26141-1-134 .26141_1_1…    489   1121 -      catB… 1-633/6… ===========… 0/0  
## 5 26141-1-134 .26141_1_1…   1259   2089 -      blaO… 1-831/8… ===========… 0/0  
## 6 26141-1-134 .26141_1_1…   2220   2819 -      aac(… 1-600/6… ===========… 0/0  
## # ℹ 6 more variables: `%COVERAGE` <dbl>, `%IDENTITY` <dbl>, DATABASE <chr>,
## #   ACCESSION <chr>, PRODUCT <chr>, RESISTANCE <chr>
```

```r
class(ab_tab)
```

```
## [1] "spec_tbl_df" "tbl_df"      "tbl"         "data.frame"
```

# Wrangling the abricate table

## Pre-processing the dataset

First Remove duplicates. These cause issues with creating a lit rather than a dataframe. When compared against the database there can sometimes be duplicates.


```r
rm_duplicates = function(ab_tab)
{
  
  # First Remove duplicates
  # These cause issues with creating a lit rather than a dataframe
  # When compared against the database there can sometimes be duplciates 
  
  # arrange and group by ab_tab 
  ab_tab_unique <- ab_tab %>% 
    arrange(GENE, -`%COVERAGE`) %>% 
    group_by(GENE)
  
  # Check the colnames
  colnames(ab_tab_unique)
  
  # Check the duplicates 
  duplicated(ab_tab_unique[c("#FILE","GENE")])
  print("Duplicates found:")
  print(c(which(duplicated(ab_tab_unique[c("#FILE","GENE")]) == TRUE))) # see which contain duplicates
  
  # Remove columns that contain duplicates in FILE, GENE and COVERAGE 
  ab_tab_unique <- ab_tab[!duplicated(ab_tab_unique[c("#FILE","GENE")]),] #using colnames
  #ab_tab_unique <- ab_tab_unique[!duplicated(ab_tab_unique[c(1,2)]),] #using col numbers
  
  
  # Check whether duplicates have been delted 
  print("Duplicates remaining:")
  print(which(duplicated(ab_tab_unique[c("#FILE","GENE")]) == TRUE)) # see which contain duplicates
  
  {
    return(ab_tab_unique)
  }
}


# Run function 
ab_tab_unique = rm_duplicates(ab_tab)
```

```
## [1] "Duplicates found:"
##  [1] 1640 1789 2297 2463 4681 4685 5413 5580 6256 7477 7743 7997 8223 8492 8494
## [16] 8496 8498 8501 8503 8505 8507 8510 8513 8519 8531 8534 8537 8539 8541 8543
## [31] 8550 8564 8567 8571 8575 8577 8580 8582 8586 8588 8929 9661
## [1] "Duplicates remaining:"
##  [1]  608 2662 2663 4920 4921 4922 5806 5822 5904 5929 5987 6064 6113 6143 6197
## [16] 6247 6352 6417 6509 6868 6869 6981 7016 7073 7090 7151 7197 7395 7839 7912
## [31] 7982 8120 8166 8307 8322 8358 8389 8538 9709 9758
```





```r
rename_catB4 = function(ab_tab)
{
  
  # Rename any genes called catB3_1 with a % cverage lower than in the GENE column 
  ab_tab_new <- mutate(ab_tab, GENE = ifelse(GENE == "catB3_1" & "%COVERAGE" < 90, "catB4", GENE))
  test = ab_tab_new %>% filter(GENE == "catB4")
  print("amount of catB4 found:")
  print(nrow(test))
  
  {
    return(ab_tab_new)
  }
}


ab_tab_unique = rename_catB4(ab_tab_unique)
```

```
## [1] "amount of catB4 found:"
## [1] 237
```

## Creating a count table

Next we make a count table of the 


```r
mk_count_table = function(ab_tab)
{
  # Funtion to make the Abricate table wide 
  # Create a new value from original ab_tab to wide count table 
  wide = ab_tab %>% 
    # select strain, gene and resistance and rename
    select("#FILE", GENE, RESISTANCE) %>% 
    # rename 
    rename(gene = "GENE",
           file = "#FILE",
           resistance = "RESISTANCE") %>% 
    # count genes vs strains
    count(file, gene) %>% 
    # Convert to wide table with genes as column headers and n as values 
    pivot_wider(names_from = gene,
                values_from = n)
  
  # Shorten strain names (remove the unnecessary string)
  wide$file = gsub(".fa", "", wide$file)
  wide$file = gsub("_assembled.fasta", "", wide$file)
  
  {
    return(wide)
  }
}

# Run function 
wide = mk_count_table(ab_tab_unique)


# Check
head(wide)
```

```
## # A tibble: 6 × 176
##   file    `aac(3)-IId_1` `aac(6')-Ib-cr_1` aadA5_1 `aph(3'')-Ib_5` `aph(6)-Id_1`
##   <chr>            <int>             <int>   <int>           <int>         <int>
## 1 26141-…              1                 1       1               1             1
## 2 26141-…              1                NA       1              NA            NA
## 3 26141-…              1                 1       1               1             1
## 4 26141-…             NA                NA       1               1             1
## 5 26141-…             NA                NA      NA               1             1
## 6 26141-…             NA                NA       1               1             1
## # ℹ 170 more variables: `blaCTX-M-15_1` <int>, `blaOXA-1_1` <int>,
## #   `blaTEM-1B_1` <int>, catB3_2 <int>, dfrA17_1 <int>, `mdf(A)_1` <int>,
## #   `mph(A)_2` <int>, sul2_3 <int>, `tet(A)_6` <int>, `blaCTX-M-27_1` <int>,
## #   sul1_5 <int>, sul2_2 <int>, aadA2_1 <int>, catA1_1 <int>, dfrA12_8 <int>,
## #   `tet(B)_2` <int>, `aac(3)-IIa_1` <int>, catB4 <int>,
## #   `blaCTX-M-14b_1` <int>, `ant(3'')-Ia_1` <int>, dfrA1_8 <int>,
## #   sul2_13 <int>, dfrA14_1 <int>, sul2_6 <int>, `blaSHV-12_1` <int>, …
```

```r
# Check
nrow(wide)
```

```
## [1] 772
```


Next convert this count table to a dataframe and matrix


```r
# Converts to dataframe
count_to_dataframe = function(wide)
{
  wide = wide
  
  # make the genenames a variable 
  rows = wide$file
  
  # Remove the file column
  wide3 = wide %>% select(-file)
  
  # name the rows as strains
  rownames(wide3) = rows
  
  # convert to dataframe
  d = as.data.frame(wide3)
  
  # convert na to 0
  d[is.na(d)] = 0
  
  #name rows
  rownames(d) = rows
  
  # convert to matrix
  d = d
  
  
  {
    return(d)
  }
}

d = count_to_dataframe(wide)

# Converts to matrix
count_to_mat = function(wide)
{
  wide = wide
  
  # make the genenames a variable 
  rows = wide$file
  
  # Remove the file column
  wide3 = wide %>% select(-file)
  
  # name the rows as strains
  rownames(wide3) = rows
  
  # convert to dataframe
  d = as.data.frame(wide3)
  
  # convert na to 0
  d[is.na(d)] <- 0
  
  # name rows
  rownames(d) = rows
  
  # convert to matrix
  mat = as.matrix(d)
  
  
  {
    return(mat)
  }
}

# Run Function 
mat = count_to_mat(wide)


# Check output
head(mat)
nrow(mat)
class(mat)
mat
colnames(mat)
```



# Using the cooccur package

## Creating a cooccur object


```r
# Transpose so gene names are the rows
# We need this to allow for a cooccur object 
d_t = t(d)

cooccur.ESBL.genes <- cooccur(d_t,
                           type = "spp_site",
                           thresh = TRUE,
                           spp_names = TRUE)
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |                                                                      |   1%
  |                                                                            
  |=                                                                     |   1%
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |=====                                                                 |   6%
  |                                                                            
  |=====                                                                 |   7%
  |                                                                            
  |=====                                                                 |   8%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |======                                                                |   9%
  |                                                                            
  |=======                                                               |   9%
  |                                                                            
  |=======                                                               |  10%
  |                                                                            
  |=======                                                               |  11%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |========                                                              |  12%
  |                                                                            
  |=========                                                             |  12%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |=========                                                             |  14%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |==========                                                            |  15%
  |                                                                            
  |===========                                                           |  15%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |============                                                          |  16%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |============                                                          |  18%
  |                                                                            
  |=============                                                         |  18%
  |                                                                            
  |=============                                                         |  19%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |==============                                                        |  20%
  |                                                                            
  |==============                                                        |  21%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |===============                                                       |  22%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |================                                                      |  23%
  |                                                                            
  |================                                                      |  24%
  |                                                                            
  |=================                                                     |  24%
  |                                                                            
  |=================                                                     |  25%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |==================                                                    |  26%
  |                                                                            
  |===================                                                   |  26%
  |                                                                            
  |===================                                                   |  27%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |====================                                                  |  28%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |=====================                                                 |  29%
  |                                                                            
  |=====================                                                 |  30%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |======================                                                |  31%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |=======================                                               |  32%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |=======================                                               |  34%
  |                                                                            
  |========================                                              |  34%
  |                                                                            
  |========================                                              |  35%
  |                                                                            
  |=========================                                             |  35%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |==========================                                            |  36%
  |                                                                            
  |==========================                                            |  37%
  |                                                                            
  |==========================                                            |  38%
  |                                                                            
  |===========================                                           |  38%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |============================                                          |  39%
  |                                                                            
  |============================                                          |  40%
  |                                                                            
  |============================                                          |  41%
  |                                                                            
  |=============================                                         |  41%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |==============================                                        |  42%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |==============================                                        |  44%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |===============================                                       |  45%
  |                                                                            
  |================================                                      |  45%
  |                                                                            
  |================================                                      |  46%
  |                                                                            
  |=================================                                     |  46%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |=================================                                     |  48%
  |                                                                            
  |==================================                                    |  48%
  |                                                                            
  |==================================                                    |  49%
  |                                                                            
  |===================================                                   |  49%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |===================================                                   |  51%
  |                                                                            
  |====================================                                  |  51%
  |                                                                            
  |====================================                                  |  52%
  |                                                                            
  |=====================================                                 |  52%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=====================================                                 |  54%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |======================================                                |  55%
  |                                                                            
  |=======================================                               |  55%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |========================================                              |  56%
  |                                                                            
  |========================================                              |  57%
  |                                                                            
  |========================================                              |  58%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |=========================================                             |  59%
  |                                                                            
  |==========================================                            |  59%
  |                                                                            
  |==========================================                            |  60%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |===========================================                           |  62%
  |                                                                            
  |============================================                          |  62%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |============================================                          |  64%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |=============================================                         |  65%
  |                                                                            
  |==============================================                        |  65%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |===============================================                       |  66%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |===============================================                       |  68%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |================================================                      |  69%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |=================================================                     |  71%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |==================================================                    |  72%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |===================================================                   |  73%
  |                                                                            
  |===================================================                   |  74%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |=====================================================                 |  75%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |======================================================                |  76%
  |                                                                            
  |======================================================                |  77%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |=======================================================               |  78%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |========================================================              |  79%
  |                                                                            
  |========================================================              |  80%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |=========================================================             |  81%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |==========================================================            |  82%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |==========================================================            |  84%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |===========================================================           |  85%
  |                                                                            
  |============================================================          |  85%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |=============================================================         |  86%
  |                                                                            
  |=============================================================         |  87%
  |                                                                            
  |=============================================================         |  88%
  |                                                                            
  |==============================================================        |  88%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |===============================================================       |  89%
  |                                                                            
  |===============================================================       |  90%
  |                                                                            
  |===============================================================       |  91%
  |                                                                            
  |================================================================      |  91%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |=================================================================     |  92%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |=================================================================     |  94%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |==================================================================    |  95%
  |                                                                            
  |===================================================================   |  95%
  |                                                                            
  |===================================================================   |  96%
  |                                                                            
  |====================================================================  |  96%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |====================================================================  |  98%
  |                                                                            
  |===================================================================== |  98%
  |                                                                            
  |===================================================================== |  99%
  |                                                                            
  |======================================================================|  99%
  |                                                                            
  |======================================================================| 100%
```
## Checking a cooccur object 


```r
class(cooccur.ESBL.genes)
```

```
## [1] "cooccur"
```

```r
summary(cooccur.ESBL.genes)
```

```
## Call:
## cooccur(mat = d_t, type = "spp_site", thresh = TRUE, spp_names = TRUE)
## 
## Of 15225 species pair combinations, 12528 pairs (82.29 %) were removed from the analysis because expected co-occurrence was < 1 and 2697 pairs were analyzed
## 
## Cooccurrence Summary:
```

```
##        Species          Sites       Positive       Negative         Random 
##          175.0          772.0          607.0          703.0         1387.0 
## Unclassifiable Non-random (%) 
##            0.0           48.6 
## attr(,"class")
## [1] "summary.cooccur"
```



## Plotting cooccurence using cooccur function


```r
cooccur_object = cooccur.ESBL.genes

    ##
    allargs <- match.call(expand.dots = TRUE)
    plotrand <- allargs$plotrand
    plotrand <- ifelse(test = is.null(plotrand),yes = FALSE,no = plotrand)
    randsummary<- allargs$randsummary
    randsummary <- ifelse(test = is.null(randsummary),yes = FALSE,no = randsummary)
    
    ##
    
    # Change as necessary - this is the number of genes in our case 169 for all 
    dim = cooccur_object$species
    #
    comat_pos = matrix(nrow=dim,ncol=dim)
    comat_neg = comat_pos 
    
    # Change as necessary - this is the full results table 
    co_tab <- cooccur_object$result
    #
    
    # Create co_occurrence 
    for (i in 1:nrow(co_tab)){
      comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
      comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]
      
      row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])
      
    }
    for (i in 1:nrow(co_tab)){
      comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
      comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
    }
    
    # Join positive and negative cooccurence matrix and label as 0.05, 0 or -1
    comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
    colnames(comat) <- 1:dim
    row.names(comat) <- 1:dim
    
    # Name the co-occurence matrix
    
    colnames(comat) = unique(cooccur_object$spp.names)
    rownames(comat) = unique(cooccur_object$spp.names)
    
    #ind <- apply(comat, 1, function(x) all(is.na(x)))
    #comat <- comat[!ind,]
    #ind <- apply(comat, 2, function(x) all(is.na(x)))
    #comat <- comat[,!ind]
    
    comat[is.na(comat)] <- 0
    
    origN <- nrow(comat)
    
    # SECTION TO REMOVE SPECIES INTERACTION WITH NO OTHERS
    
    #rmrandomspp <- function(orimat,plotrand = FALSE,randsummary = FALSE){
    if(plotrand == FALSE){
      ind <- apply(comat, 1, function(x) all(x==0))
      comat <- comat[!ind,]    
      ind <- apply(comat, 2, function(x) all(x==0))
      comat <- comat[,!ind]
      #ind <- apply(orimat, 1, function(x) all(x==0))
      #orimat <- orimat[!ind,]    
      #ind <- apply(orimat, 2, function(x) all(x==0))
      #orimat <- orimat[,!ind]
    }
    #return(orimat)
    #}
    
    #comat <- rmrandomspp(orimat = comat, dots)
    #_____________________________________
    
    postN <- nrow(comat)
    
    
    ##comat <- comat[order(rowSums(comat)),]
    ##comat <- comat[,order(colSums(comat))]
    
    #comat <- rmrandomspp(orimat = comat, ...)
    
    #ind <- apply(comat, 1, function(x) all(x==0))
    #comat <- comat[!ind,]
    #ind <- apply(comat, 2, function(x) all(x==0))
    #comat <- comat[,!ind]
    
    ind <- apply(comat, 1, function(x) all(x==0))
    comat <- comat[names(sort(ind)),]
    ind <- apply(comat, 2, function(x) all(x==0))
    comat <- comat[,names(sort(ind))]
    
    #comat
    data.m = melt(comat)
    colnames(data.m) <- c("X1","X2","value")
    data.m$X1 <- as.character(data.m$X1)
    data.m$X2 <- as.character(data.m$X2)
    
    meas <- as.character(unique(data.m$X2))
    
    dfids <- subset(data.m, X1 == X2)
    
    X1 <- data.m$X1
    X2 <- data.m$X2
    
    df.lower = subset(data.m[lower.tri(comat),],X1 != X2)
    
    ##### testing the rand summary
    if(randsummary == FALSE){  
    }else{
      dim <- nrow(comat)
      ext.dim <- round(dim*0.2,digits = 0)
      if(ext.dim<0){ext.dim<-1}
      placehold <- paste("ext_", rep(c(1:ext.dim),each = dim), sep="")
      
      randcol.df <- data.frame(
        X1 = placehold,
        X2 = rep(meas,times = ext.dim),
        value = rep(x = c(-2), times = dim*ext.dim))
      
      df.lower <- rbind(df.lower,randcol.df)
      meas <- c(meas,unique(placehold))
    }
```

```
## NULL
```

```r
    #_______________________
    
    X1 <- df.lower$X1
    X2 <- df.lower$X2
    value <- df.lower$value
    
    
    
    ####
    if(randsummary == FALSE){  
      p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1))), colour ="white") 
      p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue"), name = "", labels = c("negative","random","positive"),drop=FALSE) + 
        theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) + 
        ggtitle("AMR gene Co-occurrence Matrix") + 
        xlab("") + ylab("") + 
        scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) + 
        scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) 
      p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
      
      
    }else{
      
      p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1,-2))), colour ="white") 
      p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue","light gray"), name = "", labels = c("negative","random","positive","random"),drop=FALSE) + 
        theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) + 
        ggtitle("Species Co-occurrence Matrix") + 
        xlab("") + ylab("") + 
        scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) + 
        scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) 
      p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
      
      dim <- nrow(comat)
      ext_x <- dim + 0.5 #(ext.dim/2)
      ext_y <- dim + 1
      nrem <- origN - postN
      randtext <- paste(nrem, " completely\nrandom species")
      ext_dat <- data.frame(ext_x=ext_x,ext_y=ext_y,randtext=randtext)
      
      p <- p + geom_text(data=ext_dat,aes(x = ext_x,y = ext_y,label=randtext),hjust=0,vjust=0, color="dark gray")
    }
    ####
    
    p
```

<img src="Co-occurrence-analysis-of-AMR-genes_files/figure-html/Plotting a cooccur object-1.png" width="672" />

Save cooccurence matrix as comat 


```r
 write.csv(comat, file = "out_tabs/CoccuR_cooccurence_matrix.csv")
```


## Plotting cooccurence using pheatmap



```r
# Check comat object 

comat_pheatmap = comat

colnames(comat_pheatmap) = gsub("_1", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_2", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_3", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_4", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_5", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_6", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_7", "", colnames(comat_pheatmap))
colnames(comat_pheatmap) = gsub("_8", "", colnames(comat_pheatmap))

rownames(comat_pheatmap) = gsub("_1", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_2", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_3", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_4", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_5", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_6", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_7", "", rownames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("_8", "", rownames(comat_pheatmap))

colnames(comat_pheatmap) = gsub("catB3_2", "catB3", colnames(comat_pheatmap))
rownames(comat_pheatmap) = gsub("catB3_2", "catB3", rownames(comat_pheatmap))

diag(comat_pheatmap) = 2

pheatmap(mat = comat_pheatmap,
         color = c("#FFCC66", "dark gray", "light blue", "white"),
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         fontsize = 10, 
         legend_breaks = c(-1, 0, 1, 2),
         legend_labels = c("negative", "random", "positive", "same gene"))
```

<img src="Co-occurrence-analysis-of-AMR-genes_files/figure-html/Plotting cooccurrence of all genes with pheatmap-1.png" width="672" />





### subsetting from probability table  table 


```r
# sp1 - Numeric label giving the identity of species 1, assigned based on the order in the input matrix
# sp2 - Numeric label for species 2
# sp1_inc - Number of sites (or samples) that have species 1
# sp2_inc - Number of sites that have species 2
# obs_cooccur - Observed number of sites having both species
# prob_cooccur - Probability that both species occur at a site
# exp_cooccur - Expected number of sites having both species
# p_lt - Probability that the two species would co-occur at a frequency less than the observed number of co-occurrence sites if the two species were distributed randomly (independently) of one another
# p_gt - Probability of co-occurrence at a frequency greater than the observed frequency
# sp1_name - If species names were specified in the community data matrix this field will contain the supplied name of sp1
# sp2_name - The supplied name of sp2

prob.table.ESBL.genes = prob.table(cooccur.ESBL.genes)

prob.table.ESBL.genes.filter = prob.table.ESBL.genes

prob.table.ESBL.genes.filter = prob.table.ESBL.genes.filter %>% 
  filter(sp1_name == "blaCTX-M-15_1" |
           sp1_name == "blaOXA-1_1" |
           sp1_name == "blaTEM-1B_1" |
           sp1_name == "catA1_1"|
           sp1_name == "catB3_2" |
           sp1_name == "catB4" |
           sp1_name == "aac(6')-Ib-cr_1")


prob.table.ESBL.genes.filter.full = prob.table.ESBL.genes %>% 
  filter(sp1_name == "blaCTX-M-15_1" |
           sp1_name == "blaOXA-1_1" |
           sp1_name == "blaTEM-1B_1" |
           sp1_name == "catA1_1"|
           sp1_name == "catB3_2" |
           sp1_name == "catB4" |
           sp1_name == "aac(6')-Ib-cr_1") %>% 
  filter(sp2_name == "blaCTX-M-15_1" |
           sp2_name == "blaOXA-1_1" |
           sp2_name == "blaTEM-1B_1" |
           sp2_name == "catA1_1"|
           sp2_name == "catB3_2" |
           sp2_name == "catB4" |
           sp2_name == "aac(6')-Ib-cr_1")


# Rename sp1_name
prob.table.ESBL.genes.filter.full$sp1_name = gsub("blaCTX-M-15_1", "blaCTX-M-15", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("blaOXA-1_1", "blaOXA-1", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("blaTEM-1B_1", "blaTEM-1B", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("catA1_1", "catA1", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("catB3_2", "catB3", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("catB4", "catB4", prob.table.ESBL.genes.filter.full$sp1_name)
prob.table.ESBL.genes.filter.full$sp1_name = gsub("aac(6')-Ib-cr_1", "aac(6')-Ib-cr", prob.table.ESBL.genes.filter.full$sp1_name)

# Rename sp2_name
prob.table.ESBL.genes.filter.full$sp2_name = gsub("blaCTX-M-15_1", "blaCTX-M-15", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("blaOXA-1_1", "blaOXA-1", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("blaTEM-1B_1", "blaTEM-1B", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("catA1_1", "catA1", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("catB3_2", "catB3", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("catB4", "catB4", prob.table.ESBL.genes.filter.full$sp2_name)
prob.table.ESBL.genes.filter.full$sp2_name = gsub("aac(6')-Ib-cr_1", "aac(6')-Ib-cr", prob.table.ESBL.genes.filter.full$sp2_name)

# Plot ggplot heatmap as in cooccur package

# Plotting filtered for catB3, catB4, aac(6')-Ib-cr and bla genes (prob.table.ESBL.genes.filter.full)

ptab <- prob.table.ESBL.genes.filter.full

# p_lt - Probability that the two species would co-occur at a frequency less than the observed number of co-occurrence sites if the two species were distributed randomly (independently) of one another
# p_gt - Probability of co-occurrence at a frequency greater than the observed frequency
# If p_gt is greater than 0.05 (that is not significant) then name 0, if not name 1
# If p_lt is greater than 0.05 (that is not significant) then name 0, if not name -1

ptab$signs <- ifelse(ptab$p_gt >= 0.05, 0, 1) + ifelse(ptab$p_lt >= 0.05, 0, -1)
exp_cooccur <- ptab$exp_cooccur
obs_cooccur <- ptab$obs_cooccur
signs <- ptab$signs



# _______________________
# Plotting with white background 

p <- ggplot(ptab, aes(x = sp1_name, y = sp2_name)) + 
  geom_tile(aes(fill = factor(signs, levels = c(-1, 0, 1))), colour ="white")
p <- p + scale_fill_manual(values = c("#FFCC66", "dark gray", "light blue"), 
                           name = "", 
                           labels = c("negative", "random", "positive"), 
                           drop = FALSE)
p <- p + theme(plot.title = element_text(vjust = 2, 
                                         size = 20, 
                                         face = "bold"), 
               legend.text = element_text(size = 18), 
               axis.title = element_text(size = 20), 
               axis.text = element_text(size = 18), 
               axis.text.x = element_text(angle = 90, hjust = 1), 
               panel.background = element_rect(fill='white', colour='white'),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) + 
  xlab("Gene name") + 
  ylab("Gene name")
p <- p + ggtitle("Co-occurrence across select genes") 

p
```

<img src="Co-occurrence-analysis-of-AMR-genes_files/figure-html/Subsetting from probability table-1.png" width="672" />

