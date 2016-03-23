slideFunctMean <- function(data, window, step){
	total <- length(data)
	spots <- seq(from=1, to=(total-window), by=step)
	result <- vector(length = length(spots))
	for(i in 1:length(spots)){
		result[i] <- mean(data[spots[i]:(spots[i]+window)])
	}
	return(result)
}

slideFunctMedian <- function(data, window, step){
	total <- length(data)
	spots <- seq(from=1, to=(total-window), by=step)
	result <- vector(length = length(spots))
	for(i in 1:length(spots)){
		result[i] <- median(data[spots[i]:(spots[i]+window)])
	}
	return(result)
}

step = 1;
#window = 100;
window = 100;
max_ks = 2;

data = read.csv("dombeya.csv");
data = data[order(data[,1]),];
data = subset(data, data[,1] < max_ks);

#data2 = read.csv("g_raimondii.csv");
#data2 = data2[order(data2[,1]),];
#data2 = subset(data2, data2[,1] < max_ks);

#head(data);

ks = slideFunctMedian(data[,1], window, step);
support = slideFunctMean(data[,2], window, step);

#ks2 = slideFunctMedian(data2[,1], window, step);
#support2 = slideFunctMean(data2[,2], window, step);

data = data.frame(ks, support);

#data2 = data.frame(ks2, support2);

#ks

plot(data, type="l");
#lines(data);
#lines(data2, col="red");
