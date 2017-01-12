





#' Caluclate the probability of each possible number of hydroxylations given a vector of probabilites of each hydroxylation
#' 
#' @param hprobs a vector of probabilites for each hydroxylation site in the peptide
#' @keywords "Mass spectrum", 
#' @export
#' @examples
#' ba_nhyd(c(0.1,0.2,0.9))
ba_nhyd <- function(hprobs, verbose = F){


	N <- length(hprobs)
	

	df <- data.frame(
		prob=double(),
		cbrob=double(),
		nhyd=integer(),
		depth=integer()
	)


	for(i in 1:length(hprobs)){

		if(verbose == T){
			message(sprintf("Prob %d is %0.2f, nrow(df) is %d",i,hprobs[i],nrow(df)))
		}
	
		ht <- data.frame(prob = hprobs[i], cbprob = hprobs[i], nhyd = 1, depth = i)
		ot <- data.frame(prob = 1-hprobs[i], cbprob = 1-hprobs[i], nhyd = 0, depth = i)
	
		newdf <- data.frame(
			prob=double(),
			cbrob=double(),
			nhyd=integer(),
			depth=integer()
		)
	
		if(nrow(df)==0){
	
			df <- rbind(df,ht)
			df <- rbind(df,ot)
		}
		else{
			for(j in 1:nrow(df)){df
		
				#another hydroxylation
				htn <- data.frame(
					prob = hprobs[i], 
					cbprob = hprobs[i] * df$cbprob[j],
					nhyd = df$nhyd[j] + 1,
					depth = i
					)
				newdf <- rbind(newdf,htn)	
	
				#no new hydroxylation
				otn <- data.frame(
					prob = (1-hprobs[i]), 
					cbprob = (1-hprobs[i]) * df$cbprob[j],
					nhyd = df$nhyd[j] + 0,
					depth = i
					)
				newdf <- rbind(newdf,otn)	
			}
		
			df = newdf
			rm(newdf)
		}
	
	
		for(j in 1:nrow(df)){
			#message(sprintf("beep %d",j))
			if(verbose == T){
				message(sprintf("%0.3f\t%0.3f\t%d\t%d",
							df$prob[j],df$cbprob[j],df$nhyd[j],df$depth[j]))
			}
		}
	

		result <- data.frame(
			prob = double(),
			nhyd = integer()
		)
		#finally, let's get the probabilities for each hydroxylation level.
		for( i in 0:length(hprobs)){
			r <- data.frame( prob = sum(df$cbprob[df$nhyd == i]) , nhyd=i) 	
			result <- rbind(result,r)
		}
	}	
	
	return(result)
}

