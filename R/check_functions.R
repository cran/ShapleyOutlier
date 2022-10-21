check_data <- function(X, p = NULL, n = NULL, name = NULL){
  if(is.null(name)){
    name <- paste(substitute(X))
  }
  if(!all(apply(X,2,is.numeric))){
    stop(paste("Input --", name, "-- must not contain non-numeric values"))
  }
  if(any(is.na(X))){
    stop(paste("Input --", name, "-- must not contain missing values"))
  }
  if(is.numeric(p)){
    if(ncol(X) != p){
      stop(paste("Number of columns of input --", name, "-- must be equal to", p))
    }
  }
  if(is.numeric(n)){
    if(nrow(X) != n){
      stop(paste("Number of rows of input --", name, "-- must be equal to", n))
    }
  }
}

check_matrix <- function(X, p = NULL, name = NULL){
  if(is.null(name)){
    name <- paste(substitute(X))
  }
  if(!is.numeric(X)){
    stop(paste("Input --", name, "-- must not contain non-numeric values"))
  }
  if(any(is.na(X))){
    stop(paste("Input --", name, "-- must not contain missing values"))
  }
  if(is.numeric(p)){
    if(nrow(X) != ncol(X)){
      stop(paste("Number of rows and columns of input --", name, "-- must be equal"))
    } else {
      if(nrow(X) != p){
        stop(paste("Number of rows of input --", name, "-- must be equal to", p))
      }
      if(ncol(X) != p){
        stop(paste("Number of columns of input --", name, "-- must be equal to", p))
      }
    }
  }
}


check_vector <- function(x, p = NULL, name = NULL){
  if(is.null(name)){
    name <- paste(substitute(x))
  }
  if(!is.numeric(x)){
    stop(paste("Input --", name, "-- must not contain non-numeric values"))
  }
  if(any(is.na(x))){
    stop(paste("Input --", name, "-- must not contain missing values"))
  }
  if(is.numeric(p)){
    if(length(x) != p){
      stop(paste("Length of input --", name, "-- must be equal to", p))
    }
  }
}

check_indicator <- function(x, p = NULL, name = NULL){
  if(is.null(name)){
    name <- paste(substitute(x))
  }
  if(is.numeric(x)){
    if(!all(x %in% c(0,1))){
      stop(paste("If Input --", name, "-- is numeric, it must only contain zeros and ones"))
    }
  } else if(!is.logical(x)){
    stop(paste("Input --", name, "-- must be either numeric (only containing zeros and ones) or logical"))
  }
  if(any(is.na(x))){
    stop(paste("Input --", name, "-- must not contain missing values"))
  }
  if(is.numeric(p)){
    if(length(x) != p){
      stop(paste("Length of input --", name, "-- must be equal to", p))
    }
  }
}

check_indicator_matrix <- function(X, p = NULL, n = NULL, name = NULL){
  if(is.null(name)){
    name <- paste(substitute(X))
  }
  if(is.numeric(X)){
    if(!all(X %in% c(0,1))){
      stop(paste("If Input --", name, "-- is numeric, it must only contain zeros and ones"))
    }
  } else if(!is.logical(X)){
    stop(paste("Input --", name, "-- must be either numeric (only containing zeros and ones) or logical"))
  }
  if(any(is.na(X))){
    stop(paste("Input --", name, "-- must not contain missing values"))
  }
  if(is.numeric(p)){
    if(ncol(X) != p){
      stop(paste("Number of columns of input --", name, "-- must be equal to", p))
    }
  }
  if(is.numeric(n)){
    if(nrow(X) != n){
      stop(paste("Number of rows of input --", name, "-- must be equal to", n))
    }
  }
}

