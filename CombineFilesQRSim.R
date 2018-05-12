rm(list = ls(all = TRUE))


comb.obj <- function(o1, o2){
	if( is.matrix(o1) | is.data.frame(o1)){
		o <- rbind(o1, o2)
	}
	if( is.vector(o1)){
		o <- c(o1, o2)
	}
	o
}

##### Name of file 1
#name1 <- "SNTransBB1Normal.Rdata"
#name1 <- "SNTransBoot1Laplace.Rdata"
#name1 <- "SNTransBoot1Simple.Rdata"
name1 <- "SNTransBoot3Normal.Rdata"
##### Name of file 2
#file.names <- paste( "SNTransBB", 1:4, "Normal.Rdata", sep = "") 
#file.names <- "SNTransBoot2Laplace.Rdata"
file.names <- paste("SNTransBoot", 4:6, "Normal.Rdata", sep = "")
#file.names <- "SortTransNormal2Step.Rdata"

load(name1)

object.names <- ls()
for(i in 1:length(object.names)){
	objt <- get(object.names[i])
	if( !(is.matrix(objt) | is.vector(objt) | is.data.frame(objt) ) | is.character(objt)){ 
		object.names[i] <- NA
	}
}
object.names <- na.omit(object.names)

objects.list <- vector("list", length(object.names))

for(i in 1:length(object.names)){

assign(paste(object.names[i],1, sep = ""), get(object.names[i]))
objects.list[[i]] <- get(object.names[i])

}

#file.names <- "BaseCodeSimQRSAE/OutputSNTransBoot/SNTransBoot2.Rdata"
#object.names2 <- ls()[which(ls() %in% object.names)]
iter.files <- 0
repeat{

	iter.files <- iter.files + 1
	load(file.names[iter.files])
	for(i in (1:length(object.names))[-30]){
		d.temp2 <- get(object.names[i])
		d.temp1 <- get(paste(object.names[i], 1, sep = ""))
		d.comb <- comb.obj(d.temp1, d.temp2)
		assign(object.names[i], get("d.comb"))
		assign(paste(object.names[i], 1, sep = ""), get("d.comb"))
	}
	if(iter.files == length(file.names)){break}

}

	



