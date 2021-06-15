#' Discriticizing the OTU table
#'
#'
#'
#' @param otu_table OTU table
#' @param reads reads in meta data
#'
#'
#' @return A transposed matrix with normalized OTU
#'
#' @noRd
#'
#'
#'
discriticize <- function(otu_table, reads){
  count_otu <- sapply(1:nrow(otu_table), function(i){
    normalized_otu <- otu_table[i,]/reads
    epsilon <- min(normalized_otu[normalized_otu>0])
    normalized_otu <- round(normalized_otu/(2*epsilon))
    return(normalized_otu)
  })
  count_otu <- t(count_otu)
  rownames(count_otu) <- rownames(otu_table)
  return(count_otu)
}

#' Log Transform the OTU table
#'
#'
#'
#' @param otu_table OTU table
#' @param mode steps we want to conduct
#' @param coef proportional constant
#'
#'
#' @return Count or Continuous number of log transformed OTU table
#'
#' @noRd
#'
#'
#'
normalize_otu<- function(otu_table, mode, coef=1e4){
  otu_table <- otu_table/colSums(otu_table)
  if (mode=='identify')
    otu_table <- ceiling(log(otu_table*coef+1))
  else
    otu_table <- log(otu_table*coef+1)
  return(otu_table)
}


#' Identifying function deciding the false zero part
#'
#'
#'
#' @param meta meta data/ covariates matrix
#' @param value outcome value we used
#' @param ref reference value for imputation
#' @param mode two values \code{dna} or \code{rna}
#' @param reads reads for microbiome data
#' @param is_identify indiator whether we need identify step
#'
#'
#' @return A vector for deciding false zero
#'
#'
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#'
#'
#' @noRd
#'
#'
#'
model_fit <- function(meta, value, ref, mode, reads = NULL, is_identify){
  if(is_identify){
    print("identify")
    df_all <- meta
    if(is.null(reads))
    {
      if(any(grepl("read",colnames(df_all)))==FALSE)
      {
        stop("reads information needed in meta data")
      }
      else
      {
        reads = df_all[,which(grepl("read",colnames(df_all))==TRUE)[1]]
      }
    }
    df_all$value <- value
    df_all$ref <- as.factor(ref!=0)
    df_all$reads = reads
    pos_prob <- rep(0, length(value))
    if(mode=='dna'){
      pos_prob[which(value==0&ref!=0)] <- 1
      df <- subset(df_all, ref==FALSE, select=-ref)
      filter_idx <- which(df_all$ref==FALSE)
    }
    else{
      df <- subset(df_all, ref==TRUE, select=-ref)
      filter_idx <- which(df_all$ref==TRUE)
    }

    if(any(apply(df,2,function(x){length(table(x))})==1))
    {
      df = df[,-which(apply(df,2,function(x){length(table(x))})==1)]
    }
    # if(length(levels(droplevels(df[,"batch"])))==1)
    #   df <- subset(df, select=-batch)
    # if(length(levels(droplevels(df[,"cariesfree"])))==1)
    #   df <- subset(df, select=-cariesfree)
    #
    outlier_bound <- min(mean(df_all$value)+3*sd(df_all$value), mean(df$value)+3*sd(df$value))
    df[which(df$value>outlier_bound), 'value'] <- round(outlier_bound)
    if(length(levels(as.factor(df[, 'value'])))>1 & sum(value[filter_idx]==0)>1){
      model_zinb <- pscl::zeroinfl(value~.-reads|reads, data=df, dist='negbin')
      model_nb <- MASS::glm.nb(value~.-reads, data=df, control=glm.control(maxit=30))
      model_p <- glm(value~.-reads, data=df, family=poisson)
      ##likelihood ratio test
      p1 <- as.numeric(1 - pchisq(2*(logLik(model_zinb)-logLik(model_nb)), df=2))
      p2 <- as.numeric(1 - pchisq(2*(logLik(model_zinb)-logLik(model_p)), df=3))
      p <- max(p1, p2)
      if(p<0.05){
        pi <- predict(model_zinb, type = "zero")
        mu <- predict(model_zinb, type = "count")
        pos_prob[filter_idx] <- ifelse(value[filter_idx]==0, pi/(pi+(1-pi)*exp(-mu)), 0)
      }
    }
  }
  else{
    print("no identify")
    pos_prob <- rep(0, length(value))
    pos_prob[which(value==0)] = 1
  }
  return(pos_prob)
}


#' Pick out the DNA selected to be imputed
#'
#'
#'
#' @param dna_vec DNA we aims to select
#' @param rna_vec corresponding RNA as reference
#' @param meta_dna meta data for DNA
#' @param i genera or species number
#' @param is_identify indicator of whether we need identification step
#' @param thres threshold for determining imputation set
#'
#'
#' @return List of imputed vector and confidence vector
#'
#' @noRd
#'
#'
#'
pick_out_dna <- function(dna_vec, rna_vec, meta_dna, i, is_identify, thres = 0.5){
  prob <- model_fit(meta_dna, dna_vec, rna_vec, 'dna', is_identify = is_identify)
  impute_set <- which(prob>thres)
  confidence_set <- setdiff(1:length(dna_vec), impute_set)
  impute_set = data.frame(sample=impute_set, taxa=rep(i, length(impute_set)))
  confidence_set = data.frame(sample=confidence_set, taxa=rep(i, length(confidence_set)))
  return(list(impute_set, confidence_set))
}

#' Pick out the RNA selected to be imputed
#'
#'
#'
#' @param rna_vec RNA we aims to select
#' @param dna_vec corresponding DNA as reference
#' @param meta_dna meta data for RNA
#' @param i genera or species number
#' @param is_identify indicator of whether we need identification step
#' @param thres threshold for determining imputation set
#'
#'
#' @return List of imputed vector and confidence vector
#'
#' @noRd
#'
#'
#'
pick_out_rna <- function(rna_vec, dna_vec, meta_rna, i, is_identify){
  prob <- model_fit(meta_rna, rna_vec, dna_vec, 'rna', is_identify = is_identify)
  impute_set <- which(prob>0.5)
  confidence_set <- setdiff(1:length(rna_vec), impute_set)
  impute_set = data.frame(sample=impute_set, taxa=rep(i, length(impute_set)))
  confidence_set = data.frame(sample=confidence_set, taxa=rep(i, length(confidence_set)))
  return(list(impute_set, confidence_set))
}

#' Imputation set determination
#'
#'
#'
#' @param otu_impute OTU table we aims to impute
#' @param meta_data meta data
#' @param otu_ref OTU table as reference
#' @param mode DNA or RNA we aims to impute
#' @param is_identify indicator of whether we need identification step
#'
#'
#' @return List of imputed vector and confidence vector
#'
#' @noRd
#'
#'
#'
identify_step <- function(otu_impute, meta_data, otu_ref, mode, is_identify){
  if (mode=='dna')
    results <- lapply(1:nrow(otu_impute), FUN = function(i){
      pick_out_dna(otu_impute[i,], otu_ref[i,], meta_data, i, is_identify)
    })
  else
    results <- lapply(1:nrow(otu_impute), FUN = function(i){
      pick_out_rna(otu_impute[i,], otu_ref[i,], meta_data, i, is_identify)
    })
  impute_set <- c()
  confidence_set <- c()
  for (i in 1:length(results)){
    impute_set <- rbind(impute_set, results[[i]][[1]])
    confidence_set <- rbind(confidence_set, results[[i]][[2]])
  }
  return(list(impute_set, confidence_set))
}


#' Predict meta information
#'
#'
#'
#' @param otu_impute OTU table we aims to impute
#' @param meta_data meta data
#' @param otu_ref OTU table as reference
#' @param mode DNA or RNA we aims to impute
#' @param is_identify indicator of whether we need identification step
#'
#'
#' @return List of imputed vector and confidence vector
#'
#' @noRd
#'
#'
#'
predict_meta <- function(otu_impute, confidence_set, meta){
  otu_smooth <- sapply(1:dim(otu_impute)[1], function(i){
    y <- otu_impute[i, ]
    df <- meta #data.frame(value=y, batch=meta$batch, age=meta$age, reads=meta$reads, cariesfree=meta$cariesfree)
    training_set <- unlist(subset(confidence_set, taxa==i, select=sample))

    ##remove one-level factor
    delete_col <- c()
    for(j in 1:dim(df)[2]){
      if (is.factor(df[,j])){
        level_num <- nlevels(droplevels(df[training_set, j]))
        if (level_num==1)
          delete_col <- c(delete_col, j)
      }
    }
    if(!is.null(delete_col))
      df <- df[, -delete_col]

    model <- lm(get(colnames(df)[1])~., df[training_set,])
    y_hat <- predict(model, df)
    return(y_hat)
  })
}


#' Predict each sample's imputation result
#'
#'
#' @param i ith samples we aims to impute
#' @param otu_impute OTU table we aims to impute
#' @param impute_set index for imputation set (used for prediction)
#' @param confidence_set index for reference et (used for training)
#'
#' @return cell-wise prediction result
#' @importFrom glmnet cv.glmnet
#'
#' @noRd
#'
#'
#'
predict_cellwise <- function(i, otu_impute, impute_set, confidence_set){
  training_set <- unlist(subset(confidence_set, sample==i, select=taxa))
  testing_set <- unlist(subset(impute_set, sample==i, select=taxa))
  if(length(testing_set)<3 | length(training_set)<3){
    return(otu_impute[,i])
  }
  else{
    X <- otu_impute[, -i]
    y <- otu_impute[, i]

    X_train <- X[training_set, ]
    X_test <- X[testing_set, ]
    y_train <- y[training_set]
    cv.result <- glmnet::cv.glmnet(x=X_train, y=y_train, family="gaussian", intercept=TRUE, nfolds=5 )
    mse.min <- cv.result$cvm[cv.result$lambda == cv.result$lambda.min]
    y_test <- predict(cv.result, X_test, s='lambda.min')
    ###### making the result more robust #####
    thred = mean(y_train) + 3*sd(y_train)
    y_test[which(y_test>thred)] <- thred
    ###### making the result more robust #####
    y[testing_set] <- y_test
    return(y)
  }
}

#' Predict each genera/species' imputation result
#'
#'
#' @param i ith genera/specie we aims to impute
#' @param otu_impute OTU table we aims to impute
#' @param impute_set index for imputation set (used for prediction)
#' @param confidence_set index for reference et (used for training)
#'
#' @return gene-wise prediction result
#'
#' @noRd
#'
#'
#'
predict_genewise <- function(i, otu_impute, impute_set, confidence_set){
  training_set <- unlist(subset(confidence_set, taxa==i, select=sample))
  testing_set <- unlist(subset(impute_set, taxa==i, select=sample))
  if(length(testing_set)<3 | length(training_set)<3){
    return(otu_impute[i,])
  }
  else{
    X <- otu_impute[-i,]
    X <- data.matrix(X)
    y <- otu_impute[i,]

    X_train <- t(X[, training_set])
    X_test <- t(X[, testing_set])
    y_train <- y[training_set]
    cv.result <- cv.glmnet(x = X_train, y = y_train, family = "gaussian",
                           intercept = TRUE, nfolds=5)
    y_test <- predict(cv.result, X_test)
    y_test[which(y_test<0)] <- 0
    y[testing_set] <- y_test
    return(y)
  }
}


#' Imputation conduction on each step
#'
#'
#' @param otu_impute m*n OTU table we aims to impute
#' @param meta_data n*p data frame which represents meta data
#' @param reads by default \code{reads=NULL} if \code{meta_data} contains the information of reads; otherwise a numeric vector of reads for each sample should be given
#' @param otu_ref m*n OTU table we use as reference
#' @param mode DNA or RNA we aims to impute
#' @param is_identify indicator of whether we need identification step
#' @param coef proportional constant
#'
#' @return m*n imputed OTU table
#'
#' @noRd
#'
#'
#'
impute_step <- function(otu_impute, meta_data, reads = NULL, otu_ref, mode, is_identify, coef=1e4){
  ##normalized by reads
  if(is.null(reads))
  {
    if(any(grepl("read",colnames(meta_data)))==FALSE)
    {
      stop("reads information needed in meta data")
    }
    else
    {
      reads = meta_data[,which(grepl("read",colnames(meta_data))==TRUE)[1]]
    }
  }

  otu_impute_identify <- discriticize(otu_impute, reads)
  otu_ref_identify <- discriticize(otu_ref, reads)

  identify_list <- identify_step(otu_impute_identify, meta_data, otu_ref_identify, mode, is_identify)
  impute_set <- identify_list[[1]]
  confidence_set <- identify_list[[2]]

  otu_impute_predict <- normalize_otu(otu_impute, "predict", coef)
  otu_ref_predict <- normalize_otu(otu_ref, "predict", coef)
  smooth_value <- t(predict_meta(otu_impute_predict, confidence_set, meta_data))
  rownames(smooth_value) <- rownames(otu_impute)

  confidence_matrix <- matrix(0, dim(otu_impute_predict)[1], dim(otu_impute_predict)[2])
  for(i in 1:dim(confidence_set)[1])
    confidence_matrix[confidence_set[i,]$taxa, confidence_set[i,]$sample] <- 1

  otu_impute_predict[which(confidence_matrix==1)] <- otu_impute_predict[which(confidence_matrix==1)]-smooth_value[which(confidence_matrix==1)]

  otu_impute_predict <- sapply(1:ncol(otu_impute_predict), function(i){return(
    predict_cellwise(i, otu_impute_predict, impute_set, confidence_set))
  })

  # otu_impute_predict <- t(sapply(1:nrow(otu_impute_predict), function(i){
  #   predict_genewise(i, otu_impute_predict, impute_set, confidence_set)
  # }))

  otu_impute_predict <- otu_impute_predict + smooth_value
  otu_impute_predict[which(otu_impute_predict<0)] <- 0
  new_otu_impute <- (exp(otu_impute_predict)-1)*colSums(otu_impute)/coef
  return(new_otu_impute)
}



#' JointImpute: joint imputation method for microbiome data using DNA and RNA iteratively
#'
#' The function proposed an imputation method which output the imputed data as DNA and RNA level at the same time.
#'
#' @param otu_dna A m*n data frame/ matrix that represents OTU count table for DNA
#' @param otu_rna A m*n data frame/ matrix that represents OTU count table for RNA
#' @param covariates_DNA A n*p data frame/ matrix corresponding to covariates matrix/ meta data for DNA
#' @param covariates_RNA A n*p data frame/ matrix corresponding to covariates matrix/ meta data for RNA
#' @param is_identify indicator of whether we need identify step
#' @param epochs Iteration times for imputation
#' @param coef proportional constant for normalizing
#'
#'
#' @return A list of imputed DNA and RNA matrix
#'
#' @examples
#'
#' covariates_DNA  <- data.frame(batch=factor(meta$batch.DNA), age=meta$age_mos, cariesfree=meta$cariesfree, reads=meta$reads.DNA/1e6)
#' covariates_RNA = data.frame(batch=factor(meta$batch.RNA),age=meta$age_mos,cariesfree=meta$cariesfree,reads = meta$reads.RNA/1e6)
#' otu_dna= otu[,,1]
#' otu_rna= otu[,,2]
#' otu.imp = jointimpute(otu_dna, otu_rna, covariates_DNA, covariates_RNA,is_identify = TRUE,epochs = 1)
#'
#'
#'
#' @export
#'
#'
#'
jointimpute <- function(otu_dna, otu_rna, covariates_DNA, covariates_RNA, is_identify, epochs=1, coef=1e4){
  filter_idx <- which(rowSums(otu_dna!=0)>15 & rowSums(otu_dna==0)>5 &
                        rowSums(otu_rna!=0)>5 & rowSums(otu_rna==0)>5)
  otu_dna_filtered <- otu_dna[filter_idx,]
  otu_rna_filtered <- otu_rna[filter_idx,]
  old_impute_dna <- otu_dna_filtered
  old_impute_rna <- otu_rna_filtered
  cnt <- 1
  while(cnt<=epochs){
    #print(cnt)

    print("imputing DNA")
    impute_dna <- impute_step(old_impute_dna, covariates_DNA, reads = NULL, old_impute_rna, "dna", is_identify, coef)

    print("imputing RNA")
    impute_rna <- impute_step(old_impute_rna, covariates_RNA, reads = NULL, impute_dna, "rna", is_identify, coef)

    if(all(old_impute_dna==impute_dna)&all(old_impute_rna==impute_rna)){
      break
    }

    old_impute_dna <- impute_dna
    old_impute_rna <- impute_rna
    cnt <- cnt+1
  }
  otu_dna[filter_idx,] <- impute_dna
  otu_rna[filter_idx,] <- impute_rna
  return(list(otu_dna, otu_rna))
}
